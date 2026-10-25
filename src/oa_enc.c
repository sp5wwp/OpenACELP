/* OpenACELP encoder, TETRA (EN 300 395-2 clause 4.2.2) compatible
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "openacelp.h"
#include "oa_dsp.h"
#include "oa_consts.h"
#include <string.h>

#define HIST   174                  /* excitation history kept before the subframe */
#define SP0    10                   /* speech[] offset of the current frame */
#define WSH    143                  /* weighted speech history */

/* Largest pitch gain the encoder will select (Q14). Gains near the decoder's 1.2
 * clamp make the long-term predictor marginally stable, so small numeric
 * differences between decoder implementations grow on sustained tones. */
#ifndef OA_ENC_GP_LIMIT
#define OA_ENC_GP_LIMIT OA_GP_MAX
#endif

static const int16_t q_init[10] = { 30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000 };

void oa_enc_init(oa_enc_t *e, const oa_tables_t *tab)
{
    memset(e, 0, sizeof(*e));
    e->tab = tab;
    memcpy(e->lsp_old, q_init, sizeof(q_init));
    memcpy(e->q_prev, q_init, sizeof(q_init));
    e->a_old[0] = 1 << 24;
}

/* ------------------------------------------------------------- helpers */

/* zero-state convolution y = x * h, h in Q12 */
static void conv(const int32_t *x, const int32_t *h, int32_t *y)
{
    for (int n = 0; n < OA_L; n++) {
        int64_t acc = 0;
        for (int i = 0; i <= n; i++)
            acc += (int64_t)x[i] * h[n - i];
        y[n] = oa_sat32(oa_rshr(acc, 12));
    }
}

static int64_t dot(const int32_t *a, const int32_t *b, int n)
{
    int64_t s = 0;
    for (int i = 0; i < n; i++)
        s += (int64_t)a[i] * b[i];
    return s;
}

/* s'(n) = s(n)/2 - s(n-1)/2 + alpha s'(n-1), clause 4.2.1 */
static void preprocess(oa_enc_t *e, const int16_t *in, int16_t *out)
{
    int32_t x1 = e->hp_x1;
    int64_t y1 = e->hp_y1;
    for (int n = 0; n < OA_FRAME_SAMPLES; n++) {
        int64_t y = (int64_t)(in[n] - x1) * 16384 + ((y1 * 32735) >> 15);   /* Q15 */
        if (y > ((int64_t)32767 << 15)) y = (int64_t)32767 << 15;
        if (y < -((int64_t)32768 << 15)) y = -((int64_t)32768 << 15);
        out[n] = oa_sat16(oa_rshr(y, 15));
        x1 = in[n];
        y1 = y;
    }
    e->hp_x1 = (int16_t)x1;
    e->hp_y1 = (int32_t)y1;
}

/* ------------------------------------------------------------- LP analysis (4.2.2.1) */

/* LP analysis; a[] returned in Q24 (kept at high precision for the LSP conversion) */
static void lpc_analysis(oa_enc_t *e, const int16_t *sp256, int32_t a[11])
{
    int16_t w[256];
    for (int n = 0; n < 256; n++)
        w[n] = (int16_t)oa_rshr((int32_t)sp256[n] * oa_lpc_window_q15[n], 15);

    int64_t r[11];
    for (int k = 0; k <= 10; k++) {
        int64_t s = 0;
        for (int n = k; n < 256; n++)
            s += (int32_t)w[n] * w[n - k];
        r[k] = s;
    }
    if (r[0] <= 0) {
        memcpy(a, e->a_old, 11 * sizeof(int32_t));
        return;
    }
    /* normalise r[0] below 2^31, apply the lag window, then scale to Q30 */
    int sh = 0;
    while ((r[0] >> sh) >= ((int64_t)1 << 31)) sh++;
    for (int k = 0; k <= 10; k++)
        r[k] = ((r[k] >> sh) * oa_lag_window_q30[k]) >> 30;
    int64_t R[11];
    for (int k = 0; k <= 10; k++)
        R[k] = (r[k] * ((int64_t)1 << 30)) / r[0];

    /* Levinson-Durbin, coefficients Q27, reflection Q30 */
    int64_t ac[11] = { (int64_t)1 << 27 }, an[11];
    int64_t alpha = (int64_t)1 << 30;
    for (int i = 1; i <= 10; i++) {
        int64_t acc = R[i];
        for (int j = 1; j < i; j++)
            acc += (ac[j] * R[i - j]) >> 27;
        int64_t k = -(acc * ((int64_t)1 << 30)) / alpha;
        if (k > 1073217536 || k < -1073217536) {        /* |k| > 32750/32767 */
            memcpy(a, e->a_old, 11 * sizeof(int32_t));
            return;
        }
        for (int j = 1; j < i; j++)
            an[j] = ac[j] + ((k * ac[i - j]) >> 30);
        an[i] = k >> 3;
        for (int j = 1; j <= i; j++)
            ac[j] = an[j];
        alpha = (alpha * (((int64_t)1 << 30) - ((k * k) >> 30))) >> 30;
        if (alpha <= 0) {
            memcpy(a, e->a_old, 11 * sizeof(int32_t));
            return;
        }
    }
    a[0] = 1 << 24;
    for (int i = 1; i <= 10; i++)
        a[i] = oa_sat32(oa_rshr(ac[i], 3));
    memcpy(e->a_old, a, 11 * sizeof(int32_t));
}

/* split VQ, MSE in the q domain (4.2.2.3). The splits are searched in order and
 * each is restricted to entries the decoder will accept (strictly decreasing after
 * the boundary spreading in oa_lsp_decode), falling back to the unrestricted best. */
static int ordered(const int16_t *q, int n)
{
    for (int i = 0; i + 1 < n; i++)
        if (q[i] <= q[i + 1])
            return 0;
    return 1;
}

static void lsp_quantize(const oa_tables_t *t, const int16_t lsp[10], uint16_t idx[3])
{
    const int16_t *cb[3] = { t->lsp1, t->lsp2, t->lsp3 };
    static const int size[3] = { 256, 512, 512 }, dim[3] = { 3, 3, 4 }, off[3] = { 0, 3, 6 };
    static const int16_t gap[3] = { 0, 917, 1245 };
    int16_t q[10];
    for (int s = 0; s < 3; s++) {
        int64_t best = INT64_MAX, best_any = INT64_MAX;
        int bi = -1, bany = 0;
        for (int i = 0; i < size[s]; i++) {
            const int16_t *v = &cb[s][i * dim[s]];
            int64_t err = 0;
            for (int j = 0; j < dim[s]; j++) {
                int32_t dlt = lsp[off[s] + j] - v[j];
                err += (int64_t)dlt * dlt;
            }
            if (err < best_any) { best_any = err; bany = i; }
            if (err >= best)
                continue;
            int16_t tq[10];
            memcpy(tq, q, off[s] * sizeof(int16_t));
            memcpy(&tq[off[s]], v, dim[s] * sizeof(int16_t));
            if (s > 0) {                            /* boundary spreading as in the decoder */
                int16_t *hi = &tq[off[s] - 1], *lo = &tq[off[s]];
                int32_t d = gap[s] - (*hi - *lo);
                if (d > 0) {
                    d >>= 1;
                    *hi = oa_sat16((int32_t)*hi + d);
                    *lo = oa_sat16((int32_t)*lo - d);
                }
            }
            if (ordered(tq, off[s] + dim[s])) { best = err; bi = i; }
        }
        idx[s] = (uint16_t)(bi >= 0 ? bi : bany);
        memcpy(&q[off[s]], &cb[s][idx[s] * dim[s]], dim[s] * sizeof(int16_t));
        if (s > 0) {
            int32_t d = gap[s] - (q[off[s] - 1] - q[off[s]]);
            if (d > 0) {
                d >>= 1;
                q[off[s] - 1] = oa_sat16((int32_t)q[off[s] - 1] + d);
                q[off[s]]     = oa_sat16((int32_t)q[off[s]] - d);
            }
        }
    }
}

/* ------------------------------------------------------------- open-loop pitch (4.2.2.4) */

static int open_loop_pitch(const int32_t *sw)       /* sw[-143 .. 239] */
{
    static const int lo[3] = { 20, 40, 80 }, hi[3] = { 39, 79, 142 };
    int best_k[3];
    int32_t score[3];                                 /* log2 of normalised correlation, Q16 */
    for (int r = 0; r < 3; r++) {
        int64_t cmax = -1;
        best_k[r] = lo[r];
        for (int k = lo[r]; k <= hi[r]; k++) {
            int64_t c = 0;
            for (int j = 0; j < 120; j++)
                c += ((int64_t)sw[2 * j] * sw[2 * j - k]) >> 8;
            if (c > cmax) { cmax = c; best_k[r] = k; }
        }
        int64_t en = 1;
        for (int j = 0; j < 120; j++)
            en += ((int64_t)sw[2 * j - best_k[r]] * sw[2 * j - best_k[r]]) >> 8;
        score[r] = cmax > 0 ? oa_log2((uint64_t)cmax, 0) - (oa_log2((uint64_t)en, 0) >> 1) : -(60 << 16);
    }
    /* favour lower ranges: pick r if R_r > 0.85 R_sel (log2 0.85 = -0.2345 -> Q16 -15366) */
    int sel = 2;
    if (score[1] > score[sel] - 15366) sel = 1;
    if (score[0] > score[sel] - 15366) sel = 0;
    return best_k[sel];
}

/* ------------------------------------------------------------- closed-loop pitch */

/* normalised correlation <x, y_k> / sqrt(<y_k, y_k>) of the filtered excitation at
 * integer delay k; ext points at the subframe start with history before it. Q12 */
static int64_t norm_corr(const int32_t *ext, int k, const int32_t *x, const int32_t *h)
{
    int32_t y[OA_L];
    conv(&ext[-k], h, y);
    int64_t c = dot(x, y, OA_L), en = dot(y, y, OA_L);
    if (en <= 0)
        return 0;
    return (c * 4096) / (int64_t)(oa_isqrt64((uint64_t)en) + 1);
}

/* returns the lag in thirds of a sample */
static int closed_loop_pitch(int32_t *ext, const int32_t *x, const int32_t *h,
                             int l3_min, int l3_max, int frac_limit)
{
    int kmin = (l3_min + 2) / 3, kmax = l3_max / 3;
    if (kmin < 20) kmin = 20;
    if (kmax > OA_TMAX) kmax = OA_TMAX;
    if (kmin > kmax) kmin = kmax;
    int64_t R[OA_TMAX + 6];
    for (int k = kmin - 4; k <= kmax + 4; k++)
        R[k] = norm_corr(ext, k, x, h);
    int tb = kmin;
    for (int k = kmin + 1; k <= kmax; k++)
        if (R[k] > R[tb]) tb = k;

    int best = 3 * tb;
    if (tb < frac_limit) {
        int64_t rbest = R[tb];
        for (int l3 = 3 * tb - 2; l3 <= 3 * tb + 2; l3++) {
            if (l3 == 3 * tb || l3 < l3_min || l3 > l3_max)
                continue;
            int x0 = l3 / 3, mu = l3 % 3;
            const int16_t *tap = &oa_cinterp_q15[mu * OA_CINTERP_TAPS];
            int64_t v = 0;
            for (int k = 0; k < OA_CINTERP_TAPS; k++)
                v += R[x0 + OA_CINTERP_K0 + k] * tap[k];
            v >>= 15;
            if (v > rbest) { rbest = v; best = l3; }
        }
    }
    return best;
}

/* ------------------------------------------------------------- algebraic codebook (4.2.2.5) */

typedef struct { int idx[4]; int shift; int sign; } code_t;

static void acelp_search(oa_enc_t *e, const int32_t *x2, const int32_t *hp, code_t *out)
{
    int32_t *phi = e->phi;
    int64_t dd[OA_L];
    int32_t d[64];

    /* backward-filtered target and correlation matrix */
    int64_t dmax = 1;
    for (int n = 0; n < OA_L; n++) {
        int64_t s = 0;
        for (int i = n; i < OA_L; i++)
            s += (int64_t)x2[i] * hp[i - n];
        dd[n] = s;
        if ((s < 0 ? -s : s) > dmax) dmax = s < 0 ? -s : s;
    }
    int dsh = 0;
    while ((dmax >> dsh) > 8191) dsh++;
    for (int n = 0; n < 64; n++)
        d[n] = n < OA_L ? (int32_t)(dd[n] >> dsh) : 0;

    int64_t p00 = 0;
    for (int m = 0; m < OA_L; m++)
        p00 += (int64_t)hp[m] * hp[m];
    int psh = 0;
    while ((p00 >> psh) > (1 << 23)) psh++;
    memset(phi, 0, 64 * 64 * sizeof(int32_t));
    for (int k = 0; k < OA_L; k++) {
        int64_t s = 0;
        for (int m = k; m < OA_L; m++) {          /* S_k(M) = sum_{m=k..M} hp(m) hp(m-k) */
            s += (int64_t)hp[m] * hp[m - k];
            int i = OA_L - 1 - m;                 /* phi(i, i+k) = S_k(59 - i) */
            int32_t v = (int32_t)(s >> psh);
            phi[i * 64 + i + k] = v;
            phi[(i + k) * 64 + i] = v;
        }
    }

    /* amplitudes: +1.4142, -1, +1, -1, Q14 */
    const int64_t A = 23170, A2 = 32767;          /* a, a^2 in Q14 */
    int64_t best_c = 0, best_e = 1;
    int found = 0;

    for (int shift = 0; shift < 2; shift++) {
        int P1[30], P2[8], P3[8], P4[8];
        for (int i = 0; i < 30; i++) P1[i] = 2 * i + shift;
        for (int i = 0; i < 8; i++) {
            P2[i] = 2 + 8 * i + shift;
            P3[i] = 4 + 8 * i + shift;     /* i = 7 -> absent (>= 60) */
            P4[i] = 6 + 8 * i + shift;
        }
        /* thresholds from the maximum two- and three-pulse correlations */
        int64_t max2 = 0, max3 = 0;
        for (int i0 = 0; i0 < 30; i0++)
            for (int i1 = 0; i1 < 8; i1++) {
                int64_t c2 = ((A * d[P1[i0]]) >> 14) - d[P2[i1]];
                if ((c2 < 0 ? -c2 : c2) > max2) max2 = c2 < 0 ? -c2 : c2;
                for (int i2 = 0; i2 < 8; i2++) {
                    int64_t c3 = c2 + d[P3[i2]];
                    if ((c3 < 0 ? -c3 : c3) > max3) max3 = c3 < 0 ? -c3 : c3;
                }
            }
        int64_t thr2 = (max2 * 19202) >> 15, thr3 = (max3 * 19202) >> 15;   /* 0.586 */
        int counter = 350;

        for (int i0 = 0; i0 < 30 && counter >= 0; i0++) {
            int m0 = P1[i0];
            int64_t c1 = (A * d[m0]) >> 14;
            int64_t e1 = (A2 * phi[m0 * 64 + m0]) >> 14;
            for (int i1 = 0; i1 < 8 && counter >= 0; i1++) {
                int m1 = P2[i1];
                int64_t c2 = c1 - d[m1];
                int64_t e2 = e1 + phi[m1 * 64 + m1] - ((2 * A * phi[m0 * 64 + m1]) >> 14);
                if ((c2 < 0 ? -c2 : c2) <= thr2)
                    continue;
                for (int i2 = 0; i2 < 8 && counter >= 0; i2++) {
                    int m2 = P3[i2];
                    int64_t c3 = c2 + d[m2];
                    int64_t e3 = e2 + phi[m2 * 64 + m2] + ((2 * A * phi[m0 * 64 + m2]) >> 14)
                                 - 2 * (int64_t)phi[m1 * 64 + m2];
                    if ((c3 < 0 ? -c3 : c3) > thr3) {
                        for (int i3 = 0; i3 < 8; i3++) {
                            int m3 = P4[i3];
                            int64_t c4 = c3 - d[m3];
                            int64_t e4 = e3 + phi[m3 * 64 + m3] - ((2 * A * phi[m0 * 64 + m3]) >> 14)
                                         + 2 * (int64_t)phi[m1 * 64 + m3] - 2 * (int64_t)phi[m2 * 64 + m3];
                            if (e4 <= 0)
                                continue;
                            /* c4^2 / e4 > best_c^2 / best_e ? */
                            if (!found || c4 * c4 * (best_e >> 8) > best_c * best_c * (e4 >> 8)) {
                                best_c = c4; best_e = e4; found = 1;
                                out->idx[0] = i0; out->idx[1] = i1; out->idx[2] = i2; out->idx[3] = i3;
                                out->shift = shift;
                            }
                        }
                        counter -= 3;
                    }
                }
                counter -= 4;
            }
        }
    }
    if (!found) {                                  /* silence: any valid code */
        out->idx[0] = 0; out->idx[1] = 0; out->idx[2] = 7; out->idx[3] = 7; out->shift = 0;
        best_c = 0;
    }
    out->sign = best_c < 0;
}

/* ------------------------------------------------------------- frame */

void oa_enc_frame(oa_enc_t *e, const int16_t pcm[OA_FRAME_SAMPLES], uint8_t bits[OA_FRAME_BITS])
{
    oa_params_t p;
    int16_t *sp = &e->speech[SP0];

    /* new samples go after the 40 look-ahead samples kept from the last call */
    preprocess(e, pcm, &e->speech[SP0 + 40]);

    /* LP analysis over the last 256 samples (ends with the 40 look-ahead) */
    int32_t a_unq[11];
    int16_t lsp[10], q[10], qs[4][10], ls[4][10];
    lpc_analysis(e, &e->speech[290 - 256], a_unq);
    if (oa_a_to_lsp(a_unq, lsp) != 0)
        memcpy(lsp, e->lsp_old, sizeof(lsp));
    lsp_quantize(e->tab, lsp, p.lsp);
    oa_lsp_decode(e->tab, p.lsp, e->q_prev, q);
    oa_lsp_interp(e->q_prev, q, qs);
    oa_lsp_interp(e->lsp_old, lsp, ls);

    /* weighted speech A(z/0.95)/A(z/0.6), unquantized LP, and open-loop pitch */
    int32_t *sw = &e->wsp[WSH];
    for (int s = 0; s < 4; s++) {
        int32_t a[11], num[11], den[11], x[OA_L];
        oa_lsp_to_a(ls[s], a);
        oa_weight_a(a, oa_g095_q15, num);
        oa_weight_a(a, oa_g060_q15, den);
        for (int n = 0; n < OA_L; n++)
            x[n] = sp[s * OA_L + n] * 16;
        oa_filter(num, den, x, &sw[s * OA_L], OA_L, e->wmem_num, e->wmem_den);
    }
    int t_op = open_loop_pitch(sw);

    int tmin = 20;
    for (int s = 0; s < 4; s++) {
        const int16_t *ss = &sp[s * OA_L];
        int32_t a[11], aw[11], h[OA_L], res[OA_L], x[OA_L], y[OA_L], dummy[10];
        oa_lsp_to_a(qs[s], a);
        oa_weight_a(a, oa_g085_q15, aw);
        oa_impulse(NULL, aw, h);

        /* LP residual (Q4) and target x = 1/A(z/g) residual with weighted-error memory */
        for (int n = 0; n < OA_L; n++) {
            int64_t acc = 0;
            for (int i = 0; i <= 10; i++)
                acc += (int64_t)a[i] * ss[n - i];
            res[n] = (int32_t)oa_rshr(acc, 8);
        }
        int32_t wm[10];
        memcpy(wm, e->werr, sizeof(wm));
        oa_filter(NULL, aw, res, x, OA_L, dummy, wm);

        /* adaptive codebook: search on excitation extended by the residual */
        int32_t *exc = &e->exc[HIST];
        memcpy(exc, res, sizeof(res));
        int l3;
        if (s == 0) {
            int lo = t_op - 2, hi = t_op + 2;
            if (lo < 20) lo = 20;
            if (hi > OA_TMAX) hi = OA_TMAX;
            /* fractions only below 85; the lowest coded lag is 19 1/3 (58 thirds) */
            l3 = closed_loop_pitch(exc, x, h, 3 * lo - 2 < 58 ? 58 : 3 * lo - 2, 3 * hi + 2, 85);
        } else {
            l3 = closed_loop_pitch(exc, x, h, 3 * tmin - 2, 3 * (tmin + 9) + 2, OA_TMAX + 1);
        }
        int t = (l3 + 1) / 3, frac = l3 - 3 * t;
        if (s == 0) {
            if (l3 <= 254) {
                p.sub[0].pitch = (uint8_t)(l3 - 58);
            } else {                                 /* integer lags 85..143 */
                t = (l3 + 1) / 3; frac = 0;
                if (t < 85) t = 85;
                p.sub[0].pitch = (uint8_t)(t + 112);
            }
            tmin = oa_pitch_tmin(t);
        } else {
            p.sub[s].pitch = (uint8_t)(l3 - 3 * tmin + 2);
        }
        /* re-derive the lag exactly as the decoder will */
        oa_pitch_decode(p.sub[s].pitch, s, tmin, &t, &frac);

        int32_t v[OA_L];
        oa_pred_lt(exc, t, frac);
        memcpy(v, exc, sizeof(v));
        conv(v, h, y);
        int64_t xy = dot(x, y, OA_L), yy = dot(y, y, OA_L);
        int64_t gp_u = (yy > 0 && xy > 0) ? (xy * 16384) / yy : 0;     /* Q14 */
        if (gp_u > OA_GP_MAX) gp_u = OA_GP_MAX;
        int32_t x2[OA_L];
        for (int n = 0; n < OA_L; n++)
            x2[n] = x[n] - (int32_t)oa_rshr(gp_u * y[n], 14);

        /* algebraic codebook with shaping F(z) (and pitch sharpening) folded into h' */
        int32_t f[OA_L], hp[OA_L];
        oa_shaping(a, t, f);
        conv(f, h, hp);
        code_t cd = { { 0, 0, 7, 7 }, 0, 0 };
        acelp_search(e, x2, hp, &cd);
        p.sub[s].p1 = (uint8_t)cd.idx[0];
        p.sub[s].p2 = (uint8_t)cd.idx[1];
        p.sub[s].p3 = (uint8_t)cd.idx[2];
        p.sub[s].p4 = (uint8_t)cd.idx[3];
        p.sub[s].shift = (uint8_t)cd.shift;
        p.sub[s].sign = (uint8_t)cd.sign;
        uint8_t pos[4];
        oa_pulse_positions(cd.idx[0], cd.idx[1], cd.idx[2], cd.idx[3], cd.shift, pos);
        int32_t c[OA_L], z[OA_L];
        oa_build_code(f, pos, cd.sign, c);
        conv(c, h, z);                               /* Q12 */

        /* gain VQ in the energy domain: exact decoder computation, weighted-error criterion */
        int32_t lpg = oa_log2_pg(a);
        uint64_t ev = 1u << 8, ec = 0;
        for (int n = 0; n < OA_L; n++) {
            ev += (uint64_t)((int64_t)v[n] * v[n]);
            ec += (uint64_t)((int64_t)c[n] * c[n]);
        }
        int32_t Ea = oa_log2(ev, 8) + lpg, Ef = ec ? oa_log2(ec, 24) + lpg : 0;
        int32_t pp, pc;
        oa_energy_predict(e->e_pitch, e->e_code, &pp, &pc);
        int64_t best_err = INT64_MAX;
        int best = 0;
        int64_t bgp = 0, bgc = 0;
        int32_t bEp = 0, bEc = 0;
        for (int g = 0; g < 64; g++) {
            int32_t Ep = pp + ((int32_t)e->tab->gain[2 * g] << 8);
            int32_t Ec = pc + ((int32_t)e->tab->gain[2 * g + 1] << 8);
            if (Ep > OA_EP_MAX) Ep = OA_EP_MAX;
            if (Ec > OA_EC_MAX) Ec = OA_EC_MAX;
            int64_t gp = oa_pow2((Ep + OA_K_REF - Ea) >> 1, 14);
            if (gp > OA_GP_MAX) gp = OA_GP_MAX;
            if (gp > OA_ENC_GP_LIMIT)
                continue;
            int64_t gc = ec ? oa_pow2((Ec + OA_K_REF - Ef) >> 1, 16) : 0;
            int64_t err = 0;
            for (int n = 0; n < OA_L; n++) {
                int64_t r = x[n] - oa_rshr(gp * y[n], 14) - oa_rshr(gc * z[n], 24);
                err += r * r;
            }
            if (err < best_err) {
                best_err = err; best = g; bgp = gp; bgc = gc; bEp = Ep; bEc = Ec;
            }
        }
        p.sub[s].gain = (uint8_t)best;
        e->e_pitch = bEp;
        e->e_code = bEc;

        /* final excitation, weighted-error memory, history shift */
        for (int n = 0; n < OA_L; n++)
            exc[n] = oa_sat32(oa_rshr(v[n] * bgp, 14) + oa_rshr(c[n] * bgc, 24));
        for (int i = 0; i < 10; i++) {
            int n = OA_L - 1 - i;
            e->werr[i] = oa_sat32(x[n] - oa_rshr(bgp * y[n], 14) - oa_rshr(bgc * z[n], 24));
        }
        memmove(e->exc, &e->exc[OA_L], HIST * sizeof(int32_t));
    }

    memcpy(e->q_prev, q, sizeof(q));
    memcpy(e->lsp_old, lsp, sizeof(lsp));
    memmove(e->wsp, &e->wsp[OA_FRAME_SAMPLES], WSH * sizeof(int32_t));
    memmove(e->speech, &e->speech[OA_FRAME_SAMPLES], (290 - OA_FRAME_SAMPLES) * sizeof(int16_t));

    oa_pack_params(&p, bits);
}
