/*
 * OpenACELP - free fixed-point ACELP speech codec, TETRA (EN 300 395-2) compatible
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#ifndef OPENACELP_H
#define OPENACELP_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define OA_FRAME_SAMPLES 240   /* 30 ms at 8 kHz */
#define OA_FRAME_BITS    137
#define OA_FRAME_BYTES   18    /* 137 bits packed MSB first, last 7 bits zero */

/*
 * Codebook tables.
 *   lsp1: 256 x 3, lsp2: 512 x 3, lsp3: 512 x 4   LSPs, cosine domain, Q15
 *   gain:  64 x 2                                 (R_pitch, R_code), log2 energy, Q8
 * oa_tables_tetra (src/oa_tables_tetra.c) holds the codebooks for TETRA
 * interoperability (EN 300 395-2).
 */
typedef struct {
    const int16_t *lsp1, *lsp2, *lsp3;
    const int16_t *gain;
} oa_tables_t;

/* The TETRA codebooks */
extern const oa_tables_t oa_tables_tetra;

/* Decoded parameters of one frame */
typedef struct {
    uint16_t lsp[3];
    struct {
        uint8_t pitch, p1, p2, p3, p4, sign, shift, gain;
    } sub[4];
} oa_params_t;

/* Decoder state. Treat as opaque; exposed for static allocation. */
typedef struct {
    const oa_tables_t *tab;
    int16_t  q_prev[10];        /* previous quantized LSPs, Q15 */
    int32_t  exc[234];          /* past excitation, Q4 */
    int32_t  syn[10];           /* synthesis memory, output scale, Q4 */
    int32_t  e_pitch, e_code;   /* quantized energies, Q16 log2 */
    int16_t  last_t;            /* last integer pitch lag */
    uint8_t  have_last;
    oa_params_t last;           /* last good frame parameters */
} oa_dec_t;

void oa_dec_init(oa_dec_t *d, const oa_tables_t *tab);

/* Decode one frame. bits: 137 values 0/1 (B1..B137). bfi != 0: frame lost,
 * bits are ignored and error concealment is applied. */
void oa_dec_frame(oa_dec_t *d, const uint8_t bits[OA_FRAME_BITS], int bfi,
                  int16_t pcm[OA_FRAME_SAMPLES]);

/* Encoder state. Treat as opaque; exposed for static allocation. */
typedef struct {
    const oa_tables_t *tab;
    int16_t  hp_x1;             /* pre-processing filter memory */
    int32_t  hp_y1;             /* Q15 */
    int16_t  speech[290];       /* pre-processed speech: 10 history + 240 frame + 40 look-ahead */
    int16_t  lsp_old[10];       /* previous unquantized LSPs, Q15 */
    int16_t  q_prev[10];        /* previous quantized LSPs, Q15 */
    int32_t  a_old[11];         /* previous unquantized LP filter (Levinson fallback), Q24 */
    int32_t  wsp[143 + 240];    /* weighted speech: 143 history + frame, Q4 */
    int32_t  wmem_num[10], wmem_den[10];
    int32_t  exc[234];          /* past excitation: 174 history + subframe, Q4 */
    int32_t  werr[10];          /* weighted error memory, Q4 */
    int32_t  e_pitch, e_code;   /* quantized energies, Q16 */
    int32_t  phi[64 * 64];      /* scratch: algebraic codebook correlation matrix */
} oa_enc_t;

void oa_enc_init(oa_enc_t *e, const oa_tables_t *tab);

/* Encode one frame of 240 samples (16-bit PCM, 8 kHz). Introduces 40 samples
 * of look-ahead delay, like the reference coder. */
void oa_enc_frame(oa_enc_t *e, const int16_t pcm[OA_FRAME_SAMPLES],
                  uint8_t bits[OA_FRAME_BITS]);

/* Parameter (un)packing and bit/byte helpers */
void oa_unpack_params(const uint8_t bits[OA_FRAME_BITS], oa_params_t *p);
void oa_pack_params(const oa_params_t *p, uint8_t bits[OA_FRAME_BITS]);
void oa_bits_to_bytes(const uint8_t bits[OA_FRAME_BITS], uint8_t bytes[OA_FRAME_BYTES]);
void oa_bytes_to_bits(const uint8_t bytes[OA_FRAME_BYTES], uint8_t bits[OA_FRAME_BITS]);

#ifdef __cplusplus
}
#endif

#endif
