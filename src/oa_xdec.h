/* OpenACELP: decoder primitives that reproduce the reference decoder exactly.
 * The encoder uses them too, so that its local decoder state matches the
 * decoder's.
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OA_XDEC_H
#define OA_XDEC_H

#include <stdint.h>

#define OA_XHIST      160        /* excitation history: lag 144 + 16 interpolation taps */
#define OA_LSP_GAP_34 917        /* minimum distance across the split-VQ boundaries, Q15 */
#define OA_LSP_GAP_67 1245
#define OA_G075       24576      /* 0.75, Q15 */
#define OA_G085       27853      /* 0.85, Q15 */
#define OA_XGP_MAX    4915       /* pitch gain limit 1.2, Q12 */
#define OA_EP_MAX_Q8  (27 * 256) /* quantized energy limits, log2 Q8 */
#define OA_EC_MAX_Q8  (25 * 256)
#define OA_GAIN_OFS   (-4690)    /* Q8: the Q formats of the energies, mantissas and gains */

/* push hi and lo apart until they are about gap apart (half the shortfall each, rounded down) */
void oa_lsp_spread(int16_t *hi, int16_t *lo, int16_t gap);
/* interpolated LSPs (Q15) -> A(z) coefficients, Q12 */
void oa_xlsp_to_a(const int16_t q[10], int16_t a[11]);
/* y = x / A(z); mem holds the last 10 outputs, oldest first */
void oa_xsyn(const int16_t a[11], const int16_t *x, int16_t *y, int n, int16_t mem[10]);
/* adaptive codebook vector into exc[0..59] from the history before it */
void oa_xpred_lt(int16_t *exc, int t, int frac);
/* algebraic code vector: pulses through F(z), pitch sharpening for t < 60 */
void oa_xcode(const int16_t a[11], int t, const uint8_t pos[4], int sign_neg, int16_t c[60]);
/* log2 energies of the gain computation (Q8): ea for the pitch gain, ef for the code gain */
void oa_xgain_logs(const int16_t a[11], const int16_t v[60], const int16_t c[60], int32_t *ea, int32_t *ef);
/* gains (Q12) from quantized energies (Q8) */
void oa_xgains(int32_t ep, int32_t ec, int32_t ea, int32_t ef, int16_t *gp, int16_t *gc);
/* energy prediction from the previous energies plus codebook entry (rp, rc), then gains */
void oa_xgain_quant(int16_t ep_prev, int16_t ec_prev, int16_t rp, int16_t rc, int32_t ea, int32_t ef,
                    int16_t *ep, int16_t *ec, int16_t *gp, int16_t *gc);
/* u = gp v + gc c, rounded, low 16 bits (u may alias v) */
void oa_xexcite(const int16_t v[60], const int16_t c[60], int16_t gp, int16_t gc, int16_t u[60]);

#endif
