/* OpenACELP bitstream helpers (EN 300 395-2 clause 4.2.2.7, table 3)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "openacelp.h"
#include <string.h>

static unsigned take(const uint8_t *b, int *pos, int n)
{
    unsigned v = 0;
    for (int k = 0; k < n; k++)
        v = (v << 1) | (b[(*pos)++] & 1u);
    return v;
}

static void put(uint8_t *b, int *pos, unsigned v, int n)
{
    for (int k = n - 1; k >= 0; k--)
        b[(*pos)++] = (uint8_t)((v >> k) & 1u);
}

void oa_unpack_params(const uint8_t bits[OA_FRAME_BITS], oa_params_t *p)
{
    int pos = 0;
    p->lsp[0] = (uint16_t)take(bits, &pos, 8);
    p->lsp[1] = (uint16_t)take(bits, &pos, 9);
    p->lsp[2] = (uint16_t)take(bits, &pos, 9);
    for (int s = 0; s < 4; s++) {
        p->sub[s].pitch = (uint8_t)take(bits, &pos, s == 0 ? 8 : 5);
        p->sub[s].p4    = (uint8_t)take(bits, &pos, 3);
        p->sub[s].p3    = (uint8_t)take(bits, &pos, 3);
        p->sub[s].p2    = (uint8_t)take(bits, &pos, 3);
        p->sub[s].p1    = (uint8_t)take(bits, &pos, 5);
        p->sub[s].sign  = (uint8_t)take(bits, &pos, 1);
        p->sub[s].shift = (uint8_t)take(bits, &pos, 1);
        p->sub[s].gain  = (uint8_t)take(bits, &pos, 6);
    }
}

void oa_pack_params(const oa_params_t *p, uint8_t bits[OA_FRAME_BITS])
{
    int pos = 0;
    put(bits, &pos, p->lsp[0], 8);
    put(bits, &pos, p->lsp[1], 9);
    put(bits, &pos, p->lsp[2], 9);
    for (int s = 0; s < 4; s++) {
        put(bits, &pos, p->sub[s].pitch, s == 0 ? 8 : 5);
        put(bits, &pos, p->sub[s].p4, 3);
        put(bits, &pos, p->sub[s].p3, 3);
        put(bits, &pos, p->sub[s].p2, 3);
        put(bits, &pos, p->sub[s].p1, 5);
        put(bits, &pos, p->sub[s].sign, 1);
        put(bits, &pos, p->sub[s].shift, 1);
        put(bits, &pos, p->sub[s].gain, 6);
    }
}

void oa_bits_to_bytes(const uint8_t bits[OA_FRAME_BITS], uint8_t bytes[OA_FRAME_BYTES])
{
    memset(bytes, 0, OA_FRAME_BYTES);
    for (int i = 0; i < OA_FRAME_BITS; i++)
        if (bits[i] & 1)
            bytes[i >> 3] |= (uint8_t)(0x80 >> (i & 7));
}

void oa_bytes_to_bits(const uint8_t bytes[OA_FRAME_BYTES], uint8_t bits[OA_FRAME_BITS])
{
    for (int i = 0; i < OA_FRAME_BITS; i++)
        bits[i] = (bytes[i >> 3] >> (7 - (i & 7))) & 1;
}
