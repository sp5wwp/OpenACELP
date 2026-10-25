/* OpenACELP decoder, TETRA (EN 300 395-2 clause 4.2.3) compatible
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "openacelp.h"
#include "oa_dsp.h"
#include <string.h>

#define HIST (sizeof(((oa_dec_t *)0)->exc) / sizeof(int32_t) - OA_L)   /* 174 */

static const int16_t q_init[10] = { 30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000 };

void oa_dec_init(oa_dec_t *d, const oa_tables_t *tab)
{
    memset(d, 0, sizeof(*d));
    d->tab = tab;
    memcpy(d->q_prev, q_init, sizeof(q_init));
    d->last_t = 40;
}

/* gather and sort (descending cosine) the quantized LSP vector */
/* Gather the quantized LSPs, keep a minimum distance across the two codebook
 * boundaries, and fall back to the previous vector if they are not strictly
 * decreasing (matches the ETSI reference decoder). Returns 0 if the new vector
 * was used, -1 if the previous one was kept. */
static void spread(int16_t *hi, int16_t *lo, int16_t dmin)
{
    int32_t t = dmin - (*hi - *lo);
    if (t > 0) {
        t >>= 1;
        *hi = oa_sat16((int32_t)*hi + t);
        *lo = oa_sat16((int32_t)*lo - t);
    }
}

int oa_lsp_decode(const oa_tables_t *t, const uint16_t idx[3], const int16_t prev[10], int16_t q[10])
{
    memcpy(&q[0], &t->lsp1[3 * idx[0]], 3 * sizeof(int16_t));
    memcpy(&q[3], &t->lsp2[3 * idx[1]], 3 * sizeof(int16_t));
    memcpy(&q[6], &t->lsp3[4 * idx[2]], 4 * sizeof(int16_t));
    spread(&q[2], &q[3], 917);              /* 0.028: about 50 Hz near 1000 Hz */
    spread(&q[5], &q[6], 1245);             /* 0.038: about 50 Hz near 1600 Hz */
    for (int i = 0; i < 9; i++)
        if (q[i] <= q[i + 1]) {
            memcpy(q, prev, 10 * sizeof(int16_t));
            return -1;
        }
    return 0;
}

void oa_lsp_interp(const int16_t prev[10], const int16_t cur[10], int16_t out[4][10])
{
    for (int i = 0; i < 10; i++) {
        out[0][i] = (int16_t)((3 * prev[i] + cur[i] + 2) >> 2);
        out[1][i] = (int16_t)((prev[i] + cur[i] + 1) >> 1);
        out[2][i] = (int16_t)((prev[i] + 3 * cur[i] + 2) >> 2);
        out[3][i] = cur[i];
    }
}

/* log2 of the prediction gain sum(h^2), h = impulse response of 1/A(z); Q16 */
int32_t oa_log2_pg(const int32_t a[11])
{
    int32_t h[OA_L];
oa_impulse(NULL, a, h);
    uint64_t s = 0;
    for (int n = 0; n < OA_L; n++)
        s += (uint64_t)((int64_t)h[n] * h[n]);
    return oa_log2(s, 24);
}

static void decode_sub(oa_dec_t *d, const int32_t a[11], int t, int frac, const uint8_t pos[4],
                       int sign, int gain_idx, int conceal, int16_t *pcm)
{
    int32_t *exc = &d->exc[HIST];
    int32_t f[OA_L], c[OA_L], v[OA_L];

    /* adaptive codebook */
    oa_pred_lt(exc, t, frac);
    memcpy(v, exc, sizeof(v));

    /* innovation through the shaping filter */
    oa_shaping(a, t, f);
    oa_build_code(f, pos, sign, c);

    /* energies, log2, Q16 */
    int32_t lpg = oa_log2_pg(a);
    uint64_t ev = 1u << 8, ec = 0;                   /* v Q4 -> v^2 Q8; epsilon = 1 */
    for (int n = 0; n < OA_L; n++) {
        ev += (uint64_t)((int64_t)v[n] * v[n]);
        ec += (uint64_t)((int64_t)c[n] * c[n]);       /* Q24 */
    }
    int32_t Ea = oa_log2(ev, 8) + lpg;
    int32_t Ef = ec ? oa_log2(ec, 24) + lpg : 0;

    int32_t Ep, Ec;
    if (conceal) {
        Ep = d->e_pitch - (1 << 15);               /* -0.5, floored at 0 */
        Ec = d->e_code  - (1 << 15);
        if (Ep < 0) Ep = 0;
        if (Ec < 0) Ec = 0;
    } else {
        int32_t pp, pc;
        oa_energy_predict(d->e_pitch, d->e_code, &pp, &pc);
        Ep = pp + ((int32_t)d->tab->gain[2 * gain_idx] << 8);
        Ec = pc + ((int32_t)d->tab->gain[2 * gain_idx + 1] << 8);
        if (Ep > OA_EP_MAX) Ep = OA_EP_MAX;
        if (Ec > OA_EC_MAX) Ec = OA_EC_MAX;
    }
    d->e_pitch = Ep;
    d->e_code = Ec;

    /* gains: gp Q14 (<= 1.2), gc Q16 */
    int64_t gp = oa_pow2((Ep + OA_K_REF - Ea) >> 1, 14);
    if (gp > OA_GP_MAX)
        gp = OA_GP_MAX;
    int64_t gc = ec ? oa_pow2((Ec + OA_K_REF - Ef) >> 1, 16) : 0;

    /* excitation (Q4) and synthesis in the output domain: y = 2u - sum a_i y(n-i) */
    for (int n = 0; n < OA_L; n++) {
        int64_t u = oa_rshr(v[n] * gp, 14) + oa_rshr(c[n] * gc, 24);
        exc[n] = oa_sat32(u);
    }

    /* synthesis on a linear buffer: y[0..9] = history (oldest first), y[10..69] = output */
    int32_t y[10 + OA_L];
    for (int i = 0; i < 10; i++)
        y[9 - i] = d->syn[i];                      /* syn[0] = newest */
    for (int n = 0; n < OA_L; n++) {
        int32_t *yn = &y[10 + n];
        int64_t acc = (int64_t)exc[n] * 8192;      /* 2u, Q4 * Q12 */
        acc -= (int64_t)a[1] * yn[-1] + (int64_t)a[2] * yn[-2] + (int64_t)a[3] * yn[-3]
             + (int64_t)a[4] * yn[-4] + (int64_t)a[5] * yn[-5] + (int64_t)a[6] * yn[-6]
             + (int64_t)a[7] * yn[-7] + (int64_t)a[8] * yn[-8] + (int64_t)a[9] * yn[-9]
             + (int64_t)a[10] * yn[-10];
        *yn = oa_sat32(oa_rshr(acc, 12));
        pcm[n] = oa_sat16(oa_rshr(*yn, 4));
    }
    for (int i = 0; i < 10; i++)
        d->syn[i] = y[10 + OA_L - 1 - i];

    /* shift the excitation history by one subframe */
    memmove(d->exc, &d->exc[OA_L], HIST * sizeof(int32_t));
}

void oa_dec_frame(oa_dec_t *d, const uint8_t bits[OA_FRAME_BITS], int bfi, int16_t pcm[OA_FRAME_SAMPLES])
{
    oa_params_t p;
    int16_t q[10], qs[4][10];
    int conceal = 0;

    if (!bfi) {
        oa_unpack_params(bits, &p);
        d->last = p;
        d->have_last = 1;
        oa_lsp_decode(d->tab, p.lsp, d->q_prev, q);
    } else {
        conceal = 1;
        if (!d->have_last)
            memset(&d->last, 0, sizeof(d->last));
        p = d->last;
        memcpy(q, d->q_prev, sizeof(q));
    }
    oa_lsp_interp(d->q_prev, q, qs);
    memcpy(d->q_prev, q, sizeof(q));

    int tmin = 20;
    for (int s = 0; s < 4; s++) {
        int32_t a[11];
        int t, frac;
        uint8_t pos[4];
        oa_lsp_to_a(qs[s], a);
        if (conceal) {
            t = d->last_t;
            frac = 0;
        } else {
            oa_pitch_decode(p.sub[s].pitch, s, tmin, &t, &frac);
            if (s == 0)
                tmin = oa_pitch_tmin(t);
        }
        d->last_t = (int16_t)t;
        oa_pulse_positions(p.sub[s].p1, p.sub[s].p2, p.sub[s].p3, p.sub[s].p4, p.sub[s].shift, pos);
        decode_sub(d, a, t, frac, pos, p.sub[s].sign, p.sub[s].gain, conceal, &pcm[s * OA_L]);
    }
}
