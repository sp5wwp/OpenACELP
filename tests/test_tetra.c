/* OpenACELP TETRA decoder regression tests.
 *
 * Decodes deterministic streams with the TETRA codebooks and checks hashes of
 * the PCM. The expected hashes were checked when they were recorded: every
 * stream below decodes to the same PCM, sample for sample, in an independent
 * implementation that is bit-exact with the ETSI reference decoder (see
 * docs/TETRA_NOTES.md, "Verification"). A changed hash means the decoder no
 * longer reproduces the reference.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include <stdio.h>
#include <string.h>
#include "openacelp.h"

static int failures = 0;

static uint32_t rs;
static uint32_t rnd(void)
{
    rs ^= rs << 13;
    rs ^= rs >> 17;
    rs ^= rs << 5;
    return rs;
}

static uint32_t fnv(uint32_t h, const int16_t *x, int n)
{
    for (int i = 0; i < n; i++) {
        h = (h ^ (uint8_t)x[i]) * 16777619u;
        h = (h ^ (uint8_t)((uint16_t)x[i] >> 8)) * 16777619u;
    }
    return h;
}

/* kind 0: random bits; 1: random bits with random and burst losses;
 * 2: losses before the first good frame; 3: encoder output for a synthetic voiced signal */
uint32_t oa_test_stream(int kind, uint32_t seed, int frames, void (*sink)(const uint8_t *, int, void *), void *ctx);

static int16_t pulse_train(int n)
{
    /* voiced-like test signal: decaying pulses every 57 samples, plus a slow second partial */
    int p = n % 57;
    int32_t v = p < 12 ? (int32_t)(9000 - 700 * p) * ((p & 1) ? -1 : 1) : 0;
    v += (int32_t)((n * 37) % 400) - 200;
    return (int16_t)v;
}

uint32_t oa_test_stream(int kind, uint32_t seed, int frames, void (*sink)(const uint8_t *, int, void *), void *ctx)
{
    static oa_dec_t d;
    static oa_enc_t e;
    oa_dec_init(&d, &oa_tables_tetra);
    oa_enc_init(&e, &oa_tables_tetra);
    rs = seed;
    uint32_t h = 2166136261u;
    for (int f = 0; f < frames; f++) {
        uint8_t bits[OA_FRAME_BITS];
        int bfi = 0;
        if (kind == 3) {
            int16_t in[OA_FRAME_SAMPLES];
            for (int i = 0; i < OA_FRAME_SAMPLES; i++)
                in[i] = pulse_train(f * OA_FRAME_SAMPLES + i);
            oa_enc_frame(&e, in, bits);
        } else {
            for (int i = 0; i < OA_FRAME_BITS; i++)
                bits[i] = rnd() & 1;
            if (kind == 1)
                bfi = (rnd() % 8 == 0) || (f >= 100 && f < 112);
            if (kind == 2)
                bfi = f < 5 || (f >= 40 && f < 43);
        }
        if (sink)
            sink(bits, bfi, ctx);
        int16_t pcm[OA_FRAME_SAMPLES];
        oa_dec_frame(&d, bits, bfi, pcm);
        h = fnv(h, pcm, OA_FRAME_SAMPLES);
        if (kind == 2 && f < 5)
            for (int i = 0; i < OA_FRAME_SAMPLES; i++)
                if (pcm[i]) {
                    printf("FAIL: lost frame %d before the first good frame is not silent\n", f);
                    failures++;
                    break;
                }
    }
    return h;
}

#ifndef OA_TEST_NO_MAIN
static const struct {
    int kind;
    uint32_t seed;
    int frames;
    uint32_t hash;
} cases[] = {
#include "test_tetra_hashes.inc"
};

int main(void)
{
    printf("OpenACELP TETRA decoder regression\n");
    for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        uint32_t h = oa_test_stream(cases[i].kind, cases[i].seed, cases[i].frames, NULL, NULL);
        if (h != cases[i].hash) {
            printf("FAIL: stream %u (kind %d, seed %u): hash %08x, expected %08x\n", i, cases[i].kind,
                   (unsigned)cases[i].seed, (unsigned)h, (unsigned)cases[i].hash);
            failures++;
        }
    }
    if (failures) {
        printf("%d failure(s)\n", failures);
        return 1;
    }
    printf("all %u streams match\n", (unsigned)(sizeof(cases) / sizeof(cases[0])));
    return 0;
}
#endif
