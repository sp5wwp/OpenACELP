/* OpenACELP decoder, TETRA (EN 300 395-2 clause 4.2.3) compatible
 *
 * The decoder reproduces the reference decoder's output sample for sample:
 * every quantity has the word length, scaling, rounding and overflow
 * behaviour the reference decoder exhibits. docs/TETRA_NOTES.md records how
 * each detail was established.
 *
 * All arithmetic is 16/32-bit integer (no 64-bit operations, no libm, no
 * dynamic allocation), so it runs efficiently on Cortex-M4.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "openacelp.h"
#include "oa_dsp.h"
#include "oa_basicop.h"
#include "oa_consts.h"
#include "oa_xdec.h"
#include <string.h>

static const int16_t q_init[10] = { 30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000 };

/* ------------------------------------------------------------------ LSPs */

void oa_lsp_spread(int16_t *hi, int16_t *lo, int16_t gap)
{
    int32_t d = (int32_t)*hi - *lo;
    if (d < gap) {
        int32_t t = (gap - d) >> 1;
        *hi = bo_sat16((int32_t)*hi + t);
        *lo = bo_sat16((int32_t)*lo - t);
    }
}

int oa_lsp_decode(const oa_tables_t *t, const uint16_t idx[3], const int16_t prev[10], int16_t q[10])
{
    memcpy(&q[0], &t->lsp1[3 * idx[0]], 3 * sizeof(int16_t));
    memcpy(&q[3], &t->lsp2[3 * idx[1]], 3 * sizeof(int16_t));
    memcpy(&q[6], &t->lsp3[4 * idx[2]], 4 * sizeof(int16_t));
    oa_lsp_spread(&q[2], &q[3], OA_LSP_GAP_34);
    oa_lsp_spread(&q[5], &q[6], OA_LSP_GAP_67);
    for (int i = 0; i < 9; i++)
        if (q[i] <= q[i + 1]) {
            memcpy(q, prev, 10 * sizeof(int16_t));
            return -1;
        }
    return 0;
}

/* Subframes 1-3 weight the previous vector by 3/4, 1/2, 1/4 (Q15), rounding
 * down; subframe 4 uses the new vector. */
void oa_lsp_interp(const int16_t prev[10], const int16_t cur[10], int16_t out[4][10])
{
    static const int32_t w[3] = { 24576, 16384, 8192 };
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < 10; i++)
            out[k][i] = (int16_t)(((int32_t)prev[i] * w[k] + (int32_t)cur[i] * (32768 - w[k])) >> 15);
    memcpy(out[3], cur, 10 * sizeof(int16_t));
}

/* x * q / 2^15 for x in Q24 held as 2^16 hi + 2 lo (lo: 15 bits): the high
 * part, the top bit of lo and the remaining 14 bits contribute separately */
static int32_t mpy32x16(int32_t x, int16_t q)
{
    int16_t hi = (int16_t)(x >> 16);
    int32_t lo = (int32_t)((uint32_t)x & 0xFFFFu) >> 1;
    int32_t b = (lo >> 14) & 1, l = lo & 0x3FFF;
    int32_t r = bo_L_mult(hi, q);
    r = bo_L_add(r, b * q);
    return bo_L_add(r, 2 * ((l * q) >> 15));
}

/* sum (or difference) polynomial from every other LSP, Q24 */
static void lsp_poly(const int16_t *q, int32_t f[6])
{
    f[0] = 1 << 24;
    f[1] = -(int32_t)q[0] * 1024;
    for (int i = 2; i <= 5; i++) {
        int16_t b = q[2 * (i - 1)];
        f[i] = f[i - 2];
        for (int j = i; j >= 2; j--) {
            f[j] = bo_L_add(f[j], f[j - 2]);
            f[j] = bo_L_sub(f[j], bo_L_shl(mpy32x16(f[j - 1], b), 1));
        }
        f[1] = bo_L_sub(f[1], (int32_t)b * 1024);
    }
}

/* round(x / 2^13), keeping the low 16 bits */
static inline int16_t q24_to_q12(int32_t x)
{
    return bo_wrap16((x >> 13) + ((x >> 12) & 1));
}

void oa_xlsp_to_a(const int16_t q[10], int16_t a[11])
{
    int32_t f1[6], f2[6];
    lsp_poly(&q[0], f1);
    lsp_poly(&q[1], f2);
    for (int i = 5; i > 0; i--) {
        f1[i] = bo_L_add(f1[i], f1[i - 1]);
        f2[i] = bo_L_sub(f2[i], f2[i - 1]);
    }
    a[0] = 4096;
    for (int i = 1; i <= 5; i++) {
        a[i]      = q24_to_q12(bo_L_add(f1[i], f2[i]));
        a[11 - i] = q24_to_q12(bo_L_sub(f1[i], f2[i]));
    }
}

/* ------------------------------------------------------------- filtering */

/* y = x / A(z), A in Q12: 32-bit saturating accumulator, rounded and
 * saturated to 16 bits. mem[0..9] holds y(-10..-1). x and y may alias. */
void oa_xsyn(const int16_t a[11], const int16_t *x, int16_t *y, int n, int16_t mem[10])
{
    int16_t buf[10 + OA_L];
    memcpy(buf, mem, 10 * sizeof(int16_t));
    /* With 4096 + sum |a[i]| <= 65535 no partial sum can leave 32 bits (|x|, |y| <=
     * 32768), so plain multiply-accumulates give the saturating result exactly. */
    int32_t asum = 4096;
    for (int i = 1; i <= 10; i++)
        asum += a[i] < 0 ? -a[i] : a[i];
    if (asum <= 65535) {
        for (int k = 0; k < n; k++) {
            const int16_t *yp = &buf[10 + k];
            int32_t s = (int32_t)x[k] * a[0] + 2048;
            for (int i = 1; i <= 10; i++)
                s -= (int32_t)a[i] * yp[-i];
            buf[10 + k] = bo_sat16(s >> 12);
        }
    } else {
        for (int k = 0; k < n; k++) {
            const int16_t *yp = &buf[10 + k];
            int32_t s = (int32_t)x[k] * a[0];
            for (int i = 1; i <= 10; i++)
                s = bo_L_sub(s, (int32_t)a[i] * yp[-i]);
            buf[10 + k] = bo_sat16(bo_L_add(s, 2048) >> 12);
        }
    }
    memcpy(y, &buf[10], n * sizeof(int16_t));
    memcpy(mem, &buf[n], 10 * sizeof(int16_t));
}

/* ap[i] = a[i] g^i, the powers of g as recursive rounded Q15 products */
static void weight_a(const int16_t a[11], int16_t g, int16_t ap[11])
{
    int16_t gi = g;
    ap[0] = a[0];
    for (int i = 1; i <= 10; i++) {
        ap[i] = bo_mult_r(a[i], gi);
        gi = bo_mult_r(gi, g);
    }
}

/* ------------------------------------------------------- adaptive codebook */

/* exc[0..59] from exc[-T...]: lag t + frac/3, frac in {-1, 0, 1}. A lag
 * shorter than the subframe repeats the samples just written. */
void oa_xpred_lt(int16_t *exc, int t, int frac)
{
    if (frac == 0) {
        for (int n = 0; n < OA_L; n++)
            exc[n] = exc[n - t];
        return;
    }
    const int16_t *tap = &oa_interp_q15[frac > 0 ? OA_INTERP_TAPS : 0];   /* mu = 2/3 or 1/3 */
    const int16_t *x = &exc[-t - (frac > 0) - 15];
    for (int n = 0; n < OA_L; n++, x++) {
        int32_t s = 0;
        for (int k = 0; k < 32; k++)
            s = bo_L_mac(s, x[k], tap[k]);
        exc[n] = bo_round(s);
    }
}

/* ------------------------------------------------------ algebraic codebook */

void oa_xcode(const int16_t a[11], int t, const uint8_t pos[4], int sign_neg, int16_t c[OA_L])
{
    int16_t num[11], den[11], fb[10 + OA_L];
    int16_t *f = &fb[10];
    weight_a(a, OA_G075, num);
    weight_a(a, OA_G085, den);

    /* f: impulse response of F(z) = A(z/0.75) / A(z/0.85), Q12 */
    memset(fb, 0, 10 * sizeof(int16_t));
    for (int n = 0; n < OA_L; n++) {
        int32_t s = n <= 10 ? bo_L_mult(num[n], den[0]) : 0;
        for (int i = 1; i <= 10; i++)
            s = bo_L_msu(s, den[i], f[n - i]);
        f[n] = bo_round(bo_L_shl(s, 3));
    }
    /* pitch sharpening f(n) += 0.8 f(n - T); 0.8 is 6554 in Q13, product floored */
    if (t < OA_L)
        for (int n = t; n < OA_L; n++)
            f[n] = bo_sat16(f[n] + (((int32_t)f[n - t] * 6554) >> 13));

    /* pulses +1.4142 (Q14 23168), -1, +1, -1, sign bit inverting all; the sum
     * keeps its low 32 bits and the code vector the low 16 bits of sum >> 14 */
    static const int16_t amp[4] = { 23168, -16384, 16384, -16384 };
    for (int n = 0; n < OA_L; n++) {
        uint32_t s = 0;
        for (int p = 0; p < 4; p++)
            if (pos[p] <= n)
                s += (uint32_t)((int32_t)amp[p] * f[n - pos[p]]);
        if (sign_neg)
            s = 0u - s;
        c[n] = bo_wrap16((int32_t)s >> 14);
    }
}

/* ------------------------------------------------------------------ gains */

/* log2 of x > 0: integer part in *e, fraction in Q15 in *f */
static void xlog2(int32_t x, int16_t *e, int16_t *f)
{
    if (x <= 0) {
        *e = 0;
        *f = 0;
        return;
    }
    int n = bo_norm_l(x);
    x = bo_L_shl(x, n);
    *e = (int16_t)(30 - n);
    int i = (x >> 25) - 32;
    int16_t a = (int16_t)((x >> 10) & 0x7FFF);
    int32_t y = (int32_t)oa_xlog2_tab[i] << 16;
    y = bo_L_msu(y, (int16_t)(oa_xlog2_tab[i] - oa_xlog2_tab[i + 1]), a);
    *f = bo_extract_h(y);
}

/* 2^(e + f/2^15), f in Q15 */
static int32_t xpow2(int16_t e, int16_t f)
{
    int i = f >> 10;
    int16_t a = (int16_t)((f & 0x3FF) << 5);
    int32_t y = (int32_t)oa_xpow2_tab[i] << 16;
    y = bo_L_msu(y, (int16_t)(oa_xpow2_tab[i] - oa_xpow2_tab[i + 1]), a);
    return bo_L_shr_r(y, 30 - e);
}

/* Q8 log2 of the product of two 16-bit mantissas, minus the normalization shifts */
static int32_t log2_q8(int32_t prod, int shifts)
{
    int16_t e, f;
    xlog2(prod, &e, &f);
    return (int32_t)(e - shifts) * 256 + (f >> 7) + OA_GAIN_OFS;
}

/* 2^(x/2) for x in Q8 (read as Q9), Q12 */
static int32_t gain_q12(int32_t x)
{
    int16_t e = (int16_t)(x >> 9);
    int16_t f = (int16_t)((x & 0x1FF) << 6);
    if (e > 30)
        e = 30;
    if (e < -1)
        return 0;
    return xpow2(e, f);
}

static int32_t sumsq(const int16_t *x, int32_t acc)
{
    for (int n = 0; n < OA_L; n++)
        acc = bo_L_add(acc, (int32_t)x[n] * x[n]);
    return acc;
}

void oa_xgain_logs(const int16_t a[11], const int16_t v[OA_L], const int16_t c[OA_L],
                   int32_t *ea, int32_t *ef)
{
    /* pg: energy of 1/A(z)'s response to an impulse of 1024 */
    int16_t h[OA_L], mem[10] = { 0 };
    memset(h, 0, sizeof(h));
    h[0] = 1024;
    oa_xsyn(a, h, h, OA_L, mem);
    int32_t pg = sumsq(h, 0), ev = sumsq(v, 1), ec = sumsq(c, 0);

    int s1 = bo_norm_l(pg), s2 = bo_norm_l(ev);
    int16_t mpg = bo_extract_h(bo_L_shl(pg, s1));
    int16_t mv = bo_extract_h(bo_L_shl(ev, s2));
    int16_t mc = bo_extract_h(ec);                    /* fixed scale */
    *ea = log2_q8(bo_L_mult(mpg, mv), s1 + s2);
    *ef = log2_q8(bo_L_mult(mpg, mc), s1);
}

void oa_xgains(int32_t ep, int32_t ec, int32_t ea, int32_t ef, int16_t *gp, int16_t *gc)
{
    int32_t g = gain_q12(ep - ea);
    *gp = (int16_t)(g > OA_XGP_MAX ? OA_XGP_MAX : g);
    *gc = bo_wrap16(gain_q12(ec - ef));
}

void oa_xgain_quant(int16_t ep_prev, int16_t ec_prev, int16_t rp, int16_t rc, int32_t ea, int32_t ef,
                    int16_t *ep, int16_t *ec, int16_t *gp, int16_t *gc)
{
    int32_t pp = ((2 * ep_prev + ec_prev) >> 2) - 768;   /* 0.5 E'p + 0.25 E'c - 3 */
    int32_t pc = ((ep_prev + 2 * ec_prev) >> 2) - 768;
    if (pp < 0) pp = 0;
    if (pc < 0) pc = 0;
    pp += rp;
    pc += rc;
    if (pp > OA_EP_MAX_Q8) pp = OA_EP_MAX_Q8;
    if (pc > OA_EC_MAX_Q8) pc = OA_EC_MAX_Q8;
    *ep = (int16_t)pp;
    *ec = (int16_t)pc;
    oa_xgains(pp, pc, ea, ef, gp, gc);
}

void oa_xexcite(const int16_t v[OA_L], const int16_t c[OA_L], int16_t gp, int16_t gc, int16_t u[OA_L])
{
    /* |gp v| < 2^28 and |gc c| <= 2^30: the sum fits 32 bits */
    for (int n = 0; n < OA_L; n++)
        u[n] = bo_wrap16(((int32_t)gp * v[n] + (int32_t)gc * c[n] + 2048) >> 12);
}

/* ------------------------------------------------------------------ frame */

void oa_dec_init(oa_dec_t *d, const oa_tables_t *tab)
{
    memset(d, 0, sizeof(*d));
    d->tab = tab;
    memcpy(d->q_prev, q_init, sizeof(q_init));
    d->last_t = 40;
}

static void decode_sub(oa_dec_t *d, const int16_t a[11], int t, int frac, const uint8_t pos[4],
                       int sign, int gain_idx, int conceal, int16_t *pcm)
{
    int16_t *exc = &d->exc[OA_XHIST];
    int16_t c[OA_L];
    int32_t ea, ef;
    int16_t gp, gc;

    oa_xpred_lt(exc, t, frac);
    oa_xcode(a, t, pos, sign, c);
    oa_xgain_logs(a, exc, c, &ea, &ef);

    if (conceal) {
        /* energies lowered by 0.5 per subframe, not below zero */
        d->e_pitch = d->e_pitch > 128 ? (int16_t)(d->e_pitch - 128) : 0;
        d->e_code = d->e_code > 128 ? (int16_t)(d->e_code - 128) : 0;
        oa_xgains(d->e_pitch, d->e_code, ea, ef, &gp, &gc);
    } else {
        oa_xgain_quant(d->e_pitch, d->e_code, d->tab->gain[2 * gain_idx], d->tab->gain[2 * gain_idx + 1],
                       ea, ef, &d->e_pitch, &d->e_code, &gp, &gc);
    }

    oa_xexcite(exc, c, gp, gc, exc);

    int16_t y[OA_L];
    oa_xsyn(a, exc, y, OA_L, d->syn);
    for (int n = 0; n < OA_L; n++)
        pcm[n] = bo_sat16(2 * (int32_t)y[n]);

    memmove(d->exc, &d->exc[OA_L], OA_XHIST * sizeof(int16_t));
}

void oa_dec_frame(oa_dec_t *d, const uint8_t bits[OA_FRAME_BITS], int bfi, int16_t pcm[OA_FRAME_SAMPLES])
{
    oa_params_t p;
    int16_t q[10], qs[4][10];
    int conceal = 0;

    if (bfi && !d->have_last) {
        /* nothing to conceal from yet: silence, and only the first previous LSP changes */
        memset(pcm, 0, OA_FRAME_SAMPLES * sizeof(int16_t));
        d->q_prev[0] = 0;
        return;
    }
    if (!bfi) {
        oa_unpack_params(bits, &p);
        d->last = p;
        d->have_last = 1;
        oa_lsp_decode(d->tab, p.lsp, d->q_prev, q);
    } else {
        conceal = 1;
        p = d->last;
        memcpy(q, d->q_prev, sizeof(q));
    }
    oa_lsp_interp(d->q_prev, q, qs);
    memcpy(d->q_prev, q, sizeof(q));

    int tmin = 20;
    for (int s = 0; s < 4; s++) {
        int16_t a[11];
        int t, frac;
        uint8_t pos[4];
        oa_xlsp_to_a(qs[s], a);
        if (conceal) {
            t = d->last_t;                         /* last integer lag, no fraction */
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
