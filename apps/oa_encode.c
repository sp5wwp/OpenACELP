/* oa_encode: encode raw 16-bit 8 kHz PCM to an ETSI-format serial file
 * (138 16-bit words per frame: BFI=0 + 137 bits) or packed bytes (-p).
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include <stdio.h>
#include <string.h>
#include "openacelp.h"

extern const oa_tables_t oa_tables_tetra;

int main(int argc, char **argv)
{
    int packed = 0, a = 1;
    if (argc > 1 && !strcmp(argv[1], "-p")) { packed = 1; a++; }
    if (argc - a != 2) {
        fprintf(stderr, "usage: %s [-p] <in.raw> <out.ser>\n", argv[0]);
        return 1;
    }
    FILE *fi = fopen(argv[a], "rb"), *fo = fopen(argv[a + 1], "wb");
    if (!fi || !fo) { perror("open"); return 1; }
    static oa_enc_t enc;
    oa_enc_init(&enc, &oa_tables_tetra);
    int16_t pcm[OA_FRAME_SAMPLES];
    uint8_t bits[OA_FRAME_BITS];
    long frames = 0;
    while (fread(pcm, sizeof(int16_t), OA_FRAME_SAMPLES, fi) == OA_FRAME_SAMPLES) {
        oa_enc_frame(&enc, pcm, bits);
        if (packed) {
            uint8_t by[OA_FRAME_BYTES];
            oa_bits_to_bytes(bits, by);
            fwrite(by, 1, sizeof(by), fo);
        } else {
            int16_t w[1 + OA_FRAME_BITS];
            w[0] = 0;
            for (int i = 0; i < OA_FRAME_BITS; i++) w[1 + i] = bits[i];
            fwrite(w, sizeof(int16_t), 1 + OA_FRAME_BITS, fo);
        }
        frames++;
    }
    fclose(fi); fclose(fo);
    fprintf(stderr, "%ld frames\n", frames);
    return 0;
}
