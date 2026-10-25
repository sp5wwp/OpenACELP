/* OpenACELP internal fixed-point DSP helpers
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OA_DSP_H
#define OA_DSP_H

#include <stdint.h>

#define OA_L      60      /* subframe length */
#define OA_ORDER  10
#define OA_TMAX   143

/* Q formats used throughout:
 *   LSP q      Q15 int16 (cosine domain)
 *   LP a[]     Q12 int32, a[0] = 4096
 *   signals    Q4  int32 (excitation, synthesis memory, weighted speech)
 *   impulse    Q12 int32
 *   energies   Q16 int32, log2 domain
 */

static inline int16_t oa_sat16(int64_t x)
{
    return x > 32767 ? 32767 : (x < -32768 ? -32768 : (int16_t)x);
}

static inline int32_t oa_sat32(int64_t x)
{
    return x > INT32_MAX ? INT32_MAX : (x < INT32_MIN ? INT32_MIN : (int32_t)x);
}

/* round-to-nearest arithmetic right shift */
static inline int64_t oa_rshr(int64_t x, int s)
{
    return s > 0 ? (x + ((int64_t)1 << (s - 1))) >> s : x * ((int64_t)1 << -s);
}

/* log2(x / 2^q) in Q16; x must be > 0 */
int32_t oa_log2(uint64_t x, int q);
/* 2^(x/65536) as a fixed-point value with q fractional bits (saturates to INT64_MAX/4) */
int64_t oa_pow2(int32_t x_q16, int q);

/* integer square root */
uint64_t oa_isqrt64(uint64_t x);

/* LSPs (Q15) -> LP coefficients a[0..10] (Q12) */
void oa_lsp_to_a(const int16_t q[10], int32_t a[11]);
/* LP coefficients (Q24, a[0] = 2^24) -> LSPs (Q15); returns 0 on success, -1 if not all roots found */
int oa_a_to_lsp(const int32_t a[11], int16_t q[10]);
/* a_out[i] = a[i] * g^i, g powers given in Q15 */
void oa_weight_a(const int32_t a[11], const int16_t gpow_q15[11], int32_t out[11]);
/* impulse response (Q12) of num(z)/den(z), both Q12; num may be NULL for 1/den */
void oa_impulse(const int32_t *num, const int32_t den[11], int32_t h[OA_L]);
/* y = num/den filtering of x, with memories (num/den Q12, signals Q4).
 * xmem: last 10 inputs (xmem[0] newest), ymem: last 10 outputs. num may be NULL (=1). */
void oa_filter(const int32_t *num, const int32_t den[11], const int32_t *x, int32_t *y, int n,
               int32_t xmem[10], int32_t ymem[10]);

/* Adaptive codebook vector: writes exc[0..L-1] from exc[-T...] with fractional
 * lag T + frac/3 (frac in {-1,0,1}); exc must have OA_TMAX+17 samples of history. */
void oa_pred_lt(int32_t *exc, int t, int frac);

/* Algebraic code vector from shaped impulse response f (Q12), pulses -> c (Q12) */
void oa_build_code(const int32_t f[OA_L], const uint8_t pos[4], int sign_neg, int32_t c[OA_L]);
/* pulse positions from indices; positions >= 60 mean absent */
void oa_pulse_positions(int p1, int p2, int p3, int p4, int shift, uint8_t pos[4]);
/* shaping filter F(z) = A(z/0.75)/A(z/0.85) impulse response, with pitch sharpening */
void oa_shaping(const int32_t a[11], int t_int, int32_t f[OA_L]);

/* pitch index <-> lag */
void oa_pitch_decode(int idx, int subframe, int tmin, int *t, int *frac);
int  oa_pitch_tmin(int t1);

/* quantized energy prediction / limits (Q16) */
void oa_energy_predict(int32_t ep_prev, int32_t ec_prev, int32_t *ep_pred, int32_t *ec_pred);

/* shared decoder/encoder helpers (oa_dec.c) */
#include "openacelp.h"
int  oa_lsp_decode(const oa_tables_t *t, const uint16_t idx[3], const int16_t prev[10], int16_t q[10]);
void oa_lsp_interp(const int16_t prev[10], const int16_t cur[10], int16_t out[4][10]);
int32_t oa_log2_pg(const int32_t a[11]);

#define OA_EP_MAX   (27 << 16)
#define OA_EC_MAX   (25 << 16)
#define OA_K_REF    348778          /* log2(40) in Q16 */
#define OA_GP_MAX   19661           /* 1.2 in Q14 */

#endif
