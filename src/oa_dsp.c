/* OpenACELP fixed-point DSP helpers
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "oa_dsp.h"
#include "oa_consts.h"
#include <string.h>

/* ---------------------------------------------------------------- math */

static inline int clz64(uint64_t x)
{
#if defined(__GNUC__)
    return __builtin_clzll(x);
#else
    int n = 0;
    while (!(x & 0x8000000000000000ull)) { x <<= 1; n++; }
    return n;
#endif
}

int32_t oa_log2(uint64_t x, int q)
{
    if (x == 0)
        return -(64 << 16);
    int n = 63 - clz64(x);
    uint64_t t = (x << (63 - n)) >> 33;            /* Q30 mantissa in [2^30, 2^31) */
    uint32_t m = (uint32_t)t - (1u << 30);
    uint32_t i = m >> 24;                          /* 6 index bits */
    uint32_t r = (m >> 8) & 0xFFFF;                /* 16 interpolation bits */
    int32_t f = (int32_t)(oa_log2_tab[i] + (((uint64_t)(oa_log2_tab[i + 1] - oa_log2_tab[i]) * r) >> 16));
    return (n - q) * 65536 + f;
}

int64_t oa_pow2(int32_t x_q16, int q)
{
    int32_t ip = x_q16 >> 16;                      /* floor */
    uint32_t fr = (uint32_t)x_q16 & 0xFFFF;
    uint32_t i = fr >> 10, r = fr & 0x3FF;
    uint64_t m = oa_pow2_tab[i] + (((uint64_t)(oa_pow2_tab[i + 1] - oa_pow2_tab[i]) * r) >> 10); /* Q30 */
    int s = ip + q - 30;
    if (s >= 0) {
        if (s > 30)
            return INT64_MAX / 4;
        return (int64_t)(m << s);
    }
    if (s < -62)
        return 0;
    return (int64_t)((m + ((uint64_t)1 << (-s - 1))) >> -s);
}

uint64_t oa_isqrt64(uint64_t x)
{
    uint64_t r = 0, b = (uint64_t)1 << 62;
    while (b > x) b >>= 2;
    while (b) {
        if (x >= r + b) { x -= r + b; r = (r >> 1) + b; }
        else r >>= 1;
        b >>= 2;
    }
    return r;
}


/* ---------------------------------------------------------------- LSP / LP */

/* F(z) coefficients from 5 LSPs (every other one), Q24 */
static void lsp_poly(const int16_t *q, int32_t f[6])
{
    f[0] = 1 << 24;
    f[1] = -((int32_t)q[0] * 1024);                /* -2 q, Q15 -> Q24 */
    for (int i = 2; i <= 5; i++) {
        int32_t b = q[2 * (i - 1)];
        f[i] = (int32_t)(-(((int64_t)b * f[i - 1]) >> 14) + 2 * (int64_t)f[i - 2]);
        for (int j = i - 1; j >= 1; j--) {
            int64_t v = f[j] - (((int64_t)b * f[j - 1]) >> 14);
            if (j >= 2)
                v += f[j - 2];
            f[j] = (int32_t)v;
        }
    }
}

void oa_lsp_to_a(const int16_t q[10], int32_t a[11])
{
    int32_t f1[6], f2[6];
    lsp_poly(&q[0], f1);
    lsp_poly(&q[1], f2);
    for (int i = 5; i > 0; i--) {
        f1[i] += f1[i - 1];
        f2[i] -= f2[i - 1];
    }
    a[0] = 4096;
    for (int i = 1; i <= 5; i++) {
        a[i]      = (int32_t)oa_rshr((int64_t)f1[i] + f2[i], 13);
        a[11 - i] = (int32_t)oa_rshr((int64_t)f1[i] - f2[i], 13);
    }
}

/* Chebyshev series evaluation of F(x), x in Q15, f in Q24, result Q24 */
static int64_t cheb(int32_t x, const int32_t f[6])
{
    int64_t b2 = (int64_t)1 << 24;
    int64_t b1 = (int64_t)x * 1024 + f[1];        /* 2x (Q15 -> Q24 is <<9, times 2) */
    for (int i = 2; i < 5; i++) {
        int64_t b0 = ((2 * (int64_t)x * b1) >> 15) - b2 + f[i];
        b2 = b1;
        b1 = b0;
    }
    return (((int64_t)x * b1) >> 15) - b2 + (f[5] >> 1);
}

int oa_a_to_lsp(const int32_t a[11], int16_t q[10])   /* a in Q24 */
{
    int32_t f1[6], f2[6];
    f1[0] = f2[0] = 1 << 24;
    for (int i = 0; i < 5; i++) {
        f1[i + 1] = a[i + 1] + a[10 - i] - f1[i];
        f2[i + 1] = a[i + 1] - a[10 - i] + f2[i];
    }
    const int32_t *coef = f1;
    int found = 0, j = 0;
    int32_t xlow = oa_lsp_grid_q15[0];
    int64_t ylow = cheb(xlow, coef);
    while (found < 10 && j < OA_GRID_POINTS - 1) {
        j++;
        int32_t xhigh = xlow;
        int64_t yhigh = ylow;
        xlow = oa_lsp_grid_q15[j];
        ylow = cheb(xlow, coef);
        if ((ylow <= 0 && yhigh >= 0) || (ylow >= 0 && yhigh <= 0)) {
            for (int k = 0; k < 4; k++) {
                int32_t xm = (xlow + xhigh) >> 1;
                int64_t ym = cheb(xm, coef);
                if ((ylow <= 0 && ym >= 0) || (ylow >= 0 && ym <= 0)) {
                    yhigh = ym; xhigh = xm;
                } else {
                    ylow = ym; xlow = xm;
                }
            }
            /* linear interpolation between (xlow, ylow) and (xhigh, yhigh) */
            int64_t dy = yhigh - ylow;
            int32_t x = xlow;
            if (dy != 0)
                x = (int32_t)(xlow - (ylow * (int64_t)(xhigh - xlow)) / dy);
            q[found++] = (int16_t)x;
            xlow = x;
            coef = (coef == f1) ? f2 : f1;
            ylow = cheb(xlow, coef);
        }
    }
    return found == 10 ? 0 : -1;
}

void oa_weight_a(const int32_t a[11], const int16_t gpow_q15[11], int32_t out[11])
{
    for (int i = 0; i <= 10; i++)
        out[i] = (int32_t)oa_rshr((int64_t)a[i] * gpow_q15[i], 15);
    out[0] = 4096;
}

void oa_impulse(const int32_t *num, const int32_t den[11], int32_t h[OA_L])
{
    for (int n = 0; n < OA_L; n++) {
        int64_t acc = num ? (n <= 10 ? (int64_t)num[n] * 4096 : 0) : (n == 0 ? (int64_t)4096 << 12 : 0);
        for (int i = 1; i <= 10 && i <= n; i++)
            acc -= (int64_t)den[i] * h[n - i];
        h[n] = oa_sat32(oa_rshr(acc, 12));
    }
}

void oa_filter(const int32_t *num, const int32_t den[11], const int32_t *x, int32_t *y, int n,
               int32_t xmem[10], int32_t ymem[10])
{
    for (int k = 0; k < n; k++) {
        int64_t acc;
        if (num) {
            acc = (int64_t)num[0] * x[k];
            for (int i = 1; i <= 10; i++)
                acc += (int64_t)num[i] * (k - i >= 0 ? x[k - i] : xmem[i - k - 1]);
        } else {
            acc = (int64_t)4096 * x[k];
        }
        for (int i = 1; i <= 10; i++)
            acc -= (int64_t)den[i] * (k - i >= 0 ? y[k - i] : ymem[i - k - 1]);
        y[k] = oa_sat32(oa_rshr(acc, 12));
    }
    /* update memories: mem[0] = newest */
    if (num) {
        int32_t t[10];
        for (int i = 0; i < 10; i++)
            t[i] = (n - 1 - i >= 0) ? x[n - 1 - i] : xmem[i - n];
        memcpy(xmem, t, sizeof(t));
    }
    int32_t t[10];
    for (int i = 0; i < 10; i++)
        t[i] = (n - 1 - i >= 0) ? y[n - 1 - i] : ymem[i - n];
    memcpy(ymem, t, sizeof(t));
}

/* ---------------------------------------------------------------- codebooks */

void oa_pred_lt(int32_t *exc, int t, int frac)
{
    if (frac == 0) {
        for (int n = 0; n < OA_L; n++)
            exc[n] = exc[n - t];
        return;
    }
    /* lag t + frac/3: frac=+1 -> x0 = n-t-1, mu = 2/3 (row 1); frac=-1 -> x0 = n-t, mu = 1/3 (row 0) */
    const int16_t *tap = &oa_interp_q15[(frac > 0 ? 1 : 0) * OA_INTERP_TAPS];
    int off = t + (frac > 0 ? 1 : 0);
    for (int n = 0; n < OA_L; n++) {
        const int32_t *x = &exc[n - off + OA_INTERP_K0];
        int64_t acc = 0;
        for (int k = 0; k < OA_INTERP_TAPS; k++)
            acc += (int64_t)x[k] * tap[k];
        exc[n] = oa_sat32(oa_rshr(acc, 15));
    }
}

void oa_pulse_positions(int p1, int p2, int p3, int p4, int shift, uint8_t pos[4])
{
    pos[0] = (uint8_t)(2 * p1 + shift);
    pos[1] = (uint8_t)(2 + 8 * p2 + shift);
    pos[2] = (uint8_t)(4 + 8 * p3 + shift);
    pos[3] = (uint8_t)(6 + 8 * p4 + shift);
}

void oa_build_code(const int32_t f[OA_L], const uint8_t pos[4], int sign_neg, int32_t c[OA_L])
{
    static const int16_t amp_q14[4] = { 23170, -16384, 16384, -16384 };  /* +1.4142, -1, +1, -1 */
    memset(c, 0, OA_L * sizeof(int32_t));
    for (int p = 0; p < 4; p++) {
        int m = pos[p];
        if (m >= OA_L)
            continue;
        int32_t a = sign_neg ? -amp_q14[p] : amp_q14[p];
        for (int n = m; n < OA_L; n++)
            c[n] += (int32_t)oa_rshr((int64_t)a * f[n - m], 14);
    }
}

void oa_shaping(const int32_t a[11], int t_int, int32_t f[OA_L])
{
    int32_t num[11], den[11];
    oa_weight_a(a, oa_g075_q15, num);
    oa_weight_a(a, oa_g085_q15, den);
    oa_impulse(num, den, f);
    if (t_int < OA_L)                              /* f(n) += 0.8 f(n-T), in place (recursive) */
        for (int n = t_int; n < OA_L; n++)
            f[n] += (int32_t)oa_rshr((int64_t)f[n - t_int] * 26214, 15);
}

/* ---------------------------------------------------------------- parameters */

int oa_pitch_tmin(int t1)
{
    int tmin = t1 - 5;
    if (tmin < 20)
        tmin = 20;
    if (tmin + 9 > OA_TMAX)
        tmin = OA_TMAX - 9;
    return tmin;
}

void oa_pitch_decode(int idx, int subframe, int tmin, int *t, int *frac)
{
    if (subframe == 0) {
        if (idx < 197) {
            *t = (idx + 2) / 3 + 19;
            *frac = idx - 3 * *t + 58;
        } else {
            *t = idx - 112;
            *frac = 0;
        }
    } else {
        int i = (idx + 2) / 3 - 1;
        *t = tmin + i;
        *frac = idx - 2 - 3 * i;
    }
}

void oa_energy_predict(int32_t ep_prev, int32_t ec_prev, int32_t *ep_pred, int32_t *ec_pred)
{
    int32_t p = ((2 * ep_prev + ec_prev) >> 2) - (3 << 16);
    int32_t c = ((ep_prev + 2 * ec_prev) >> 2) - (3 << 16);
    *ep_pred = p > 0 ? p : 0;
    *ec_pred = c > 0 ? c : 0;
}
