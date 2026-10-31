/* oa_decode: decode an ETSI-format serial file (138 16-bit words per frame:
 * BFI + 137 bits) or a packed file (-p, 18 bytes per frame) to raw 16-bit PCM.
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include <stdio.h>
#include <string.h>
#include "openacelp.h"

int main(int argc, char **argv)
{
    int packed = 0, a = 1;
    if (argc > 1 && !strcmp(argv[1], "-p")) { packed = 1; a++; }
    if (argc - a != 2) {
        fprintf(stderr, "usage: %s [-p] <in.ser> <out.raw>\n", argv[0]);
        return 1;
    }
    FILE *fi = fopen(argv[a], "rb"), *fo = fopen(argv[a + 1], "wb");
    if (!fi || !fo) { perror("open"); return 1; }
    static oa_dec_t dec;
    oa_dec_init(&dec, &oa_tables_tetra);
    uint8_t bits[OA_FRAME_BITS];
    int16_t pcm[OA_FRAME_SAMPLES];
    long frames = 0;
    for (;;) {
        int bfi = 0;
        if (packed) {
            uint8_t by[OA_FRAME_BYTES];
            if (fread(by, 1, sizeof(by), fi) != sizeof(by)) break;
            oa_bytes_to_bits(by, bits);
        } else {
            int16_t w[1 + OA_FRAME_BITS];
            if (fread(w, sizeof(int16_t), 1 + OA_FRAME_BITS, fi) != 1 + OA_FRAME_BITS) break;
            bfi = w[0] & 1;
            for (int i = 0; i < OA_FRAME_BITS; i++) bits[i] = w[1 + i] & 1;
        }
        oa_dec_frame(&dec, bits, bfi, pcm);
        fwrite(pcm, sizeof(int16_t), OA_FRAME_SAMPLES, fo);
        frames++;
    }
    fclose(fi); fclose(fo);
    fprintf(stderr, "%ld frames\n", frames);
    return 0;
}
