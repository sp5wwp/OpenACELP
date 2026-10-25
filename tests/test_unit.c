/* OpenACELP unit tests. Self-contained: uses synthetic codebooks, no ETSI data.
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "openacelp.h"
#include "oa_dsp.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("FAIL %s:%d: ", __FILE__, __LINE__); \
                                printf(__VA_ARGS__); printf("\n"); } } while (0)

static uint32_t rng = 12345;
static uint32_t rnd(void) { rng = rng * 1664525u + 1013904223u; return rng >> 8; }
static double rndu(void) { return (rnd() & 0xFFFFFF) / 16777216.0; }

/* ---- fixed-point math */
static void test_math(void)
{
    double worst = 0;
    for (int i = 0; i < 20000; i++) {
        uint64_t x = ((uint64_t)rnd() << 24 | rnd()) >> (rnd() % 40);
        if (!x) continue;
        double err = fabs(oa_log2(x, 8) / 65536.0 - (log2((double)x) - 8));
        if (err > worst) worst = err;
    }
    CHECK(worst < 1e-4, "log2 max error %g", worst);
    worst = 0;
    for (int i = 0; i < 20000; i++) {
        int32_t xq = (int32_t)(rnd() % (40 << 16)) - (20 << 16);
        double ref = pow(2.0, xq / 65536.0) * (1 << 20);
        double rel = fabs((double)oa_pow2(xq, 20) - ref) / ref;
        if (ref > 1e6 && rel > worst) worst = rel;
    }
    CHECK(worst < 1e-4, "pow2 max relative error %g", worst);
    for (uint64_t x = 0; x < 100000; x += 7) {
        uint64_t r = oa_isqrt64(x);
        CHECK(r * r <= x && (r + 1) * (r + 1) > x, "isqrt(%llu) = %llu", (unsigned long long)x, (unsigned long long)r);
    }
}

/* ---- LSP <-> LP round trip */
static void test_lsp(void)
{
    double worst = 0;
    int fail = 0;
    for (int t = 0; t < 3000; t++) {
        double w[10];
        for (int i = 0; i < 10; i++) w[i] = rndu();
        for (int i = 0; i < 10; i++) for (int j = i + 1; j < 10; j++) if (w[j] < w[i]) { double s = w[i]; w[i] = w[j]; w[j] = s; }
        int ok = 1;
        int16_t q[10], q2[10];
        for (int i = 0; i < 10; i++) {
            double om = 0.08 + (M_PI - 0.16) * w[i];       /* realistic: >= ~100 Hz spacing */
            if (i && om - acos(q[i - 1] / 32768.0) < 0.08) ok = 0;
            q[i] = (int16_t)lrint(cos(om) * 32767);
        }
        if (!ok) { t--; continue; }
        /* A(z) in Q24 from the exact (double) LSPs */
        double f1[6] = { 1 }, f2[6] = { 1 }, ad[11];
        for (int k = 0; k < 2; k++) {
            double *f = k ? f2 : f1;
            for (int i = 1; i <= 5; i++) {
                double b = q[2 * (i - 1) + k] / 32768.0;
                f[i] = -2 * b * f[i - 1] + 2 * (i >= 2 ? f[i - 2] : 0);
                for (int j = i - 1; j >= 1; j--)
                    f[j] += -2 * b * f[j - 1] + (j >= 2 ? f[j - 2] : 0);
            }
        }
        for (int i = 5; i > 0; i--) { f1[i] += f1[i - 1]; f2[i] -= f2[i - 1]; }
        ad[0] = 1;
        for (int i = 1; i <= 5; i++) { ad[i] = 0.5 * (f1[i] + f2[i]); ad[11 - i] = 0.5 * (f1[i] - f2[i]); }
        int32_t a[11];
        for (int i = 0; i <= 10; i++) a[i] = (int32_t)lrint(ad[i] * (1 << 24));
        if (oa_a_to_lsp(a, q2)) { fail++; continue; }
        for (int i = 0; i < 10; i++) {
            double e = fabs(acos(q[i] / 32768.0) - acos(q2[i] / 32768.0)) * 4000 / M_PI;
            if (e > worst) worst = e;
        }
    }
    CHECK(fail == 0, "a_to_lsp failed %d times", fail);
    CHECK(worst < 2.0, "LSP round trip worst error %.2f Hz", worst);
    printf("  LSP round trip worst error %.2f Hz\n", worst);
}

/* ---- bit packing */
static void test_bits(void)
{
    for (int t = 0; t < 1000; t++) {
        uint8_t b[OA_FRAME_BITS], b2[OA_FRAME_BITS], by[OA_FRAME_BYTES];
        for (int i = 0; i < OA_FRAME_BITS; i++) b[i] = rnd() & 1;
        oa_params_t p;
        oa_unpack_params(b, &p);
        oa_pack_params(&p, b2);
        CHECK(!memcmp(b, b2, sizeof(b)), "param pack/unpack round trip");
        oa_bits_to_bytes(b, by);
        oa_bytes_to_bits(by, b2);
        CHECK(!memcmp(b, b2, sizeof(b)), "byte pack/unpack round trip");
    }
}

/* ---- pitch index mapping is a bijection over the coded range */
static void test_pitch(void)
{
    int seen[3 * 150] = { 0 };
    for (int idx = 0; idx < 256; idx++) {
        int t, f;
        oa_pitch_decode(idx, 0, 0, &t, &f);
        int l3 = 3 * t + f;
        CHECK(l3 >= 58 && l3 <= 3 * 143 && !seen[l3], "subframe-1 index %d -> lag %d/3", idx, l3);
        seen[l3] = 1;
    }
    for (int t1 = 20; t1 <= 143; t1++) {
        int tmin = oa_pitch_tmin(t1);
        for (int idx = 0; idx < 32; idx++) {
            int t, f;
            oa_pitch_decode(idx, 1, tmin, &t, &f);
            CHECK(3 * t + f == 3 * tmin - 2 + idx, "relative index %d", idx);
        }
    }
}

/* ---- synthetic codebooks: LSP vectors around a uniform spectrum, gain grid */
static int16_t cb1[256 * 3], cb2[512 * 3], cb3[512 * 4], gcb[64 * 2];
static void make_tables(oa_tables_t *t)
{
    double base[10];
    for (int i = 0; i < 10; i++) base[i] = (i + 1) * M_PI / 11;
    struct { int16_t *cb; int n, dim, off; } s[3] = { { cb1, 256, 3, 0 }, { cb2, 512, 3, 3 }, { cb3, 512, 4, 6 } };
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < s[k].n; i++) {
            double prev = 0;
            for (int j = 0; j < s[k].dim; j++) {
                double om = base[s[k].off + j] + (rndu() - 0.5) * 0.45;
                if (j && om < prev + 0.05) om = prev + 0.05;
                prev = om;
                s[k].cb[i * s[k].dim + j] = (int16_t)lrint(cos(om) * 32767);
            }
        }
    for (int i = 0; i < 64; i++) {
        gcb[2 * i]     = (int16_t)lrint((i % 8) * 1.7 * 256);
        gcb[2 * i + 1] = (int16_t)lrint((i / 8) * 1.7 * 256);
    }
    t->lsp1 = cb1; t->lsp2 = cb2; t->lsp3 = cb3; t->gain = gcb;
}

/* ---- encoder -> decoder loopback on a synthetic voiced signal; decoder robustness */
static void test_loopback(const oa_tables_t *tab)
{
    static oa_enc_t enc;
    static oa_dec_t dec;
    oa_enc_init(&enc, tab);
    oa_dec_init(&dec, tab);
    enum { NF = 100 };
    static int16_t in[NF * 240], out[NF * 240];
    double ph = 0;
    for (int n = 0; n < NF * 240; n++) {
        double f0 = 110 + 40 * sin(2 * M_PI * n / 16000.0);
        ph += 2 * M_PI * f0 / 8000;
        double s = 0;
        for (int h = 1; h <= 25 && h * f0 < 3600; h++)
            s += sin(h * ph) / h * exp(-pow((h * f0 - 700) / 900, 2));
        in[n] = (int16_t)lrint(6000 * s * (0.6 + 0.4 * sin(2 * M_PI * n / 6000.0)));
    }
    uint8_t bits[OA_FRAME_BITS];
    for (int k = 0; k < NF; k++) {
        oa_enc_frame(&enc, &in[k * 240], bits);
        oa_dec_frame(&dec, bits, 0, &out[k * 240]);
    }
    double se = 0, sn = 0;
    for (int n = 2400; n < NF * 240; n++) {
        double d = in[n - 40] - out[n];
        se += (double)in[n - 40] * in[n - 40];
        sn += d * d;
    }
    double snr = 10 * log10(se / sn);
    printf("  loopback SNR with synthetic tables: %.1f dB\n", snr);
    CHECK(snr > 3.0, "loopback SNR %.1f dB too low", snr);

    /* garbage and erased frames must keep the decoder bounded */
    oa_dec_init(&dec, tab);
    for (int k = 0; k < 2000; k++) {
        for (int i = 0; i < OA_FRAME_BITS; i++) bits[i] = rnd() & 1;
        oa_dec_frame(&dec, bits, (rnd() % 5) == 0, out);
    }
    /* silence after garbage decays */
    static oa_enc_t e2;
    oa_enc_init(&e2, tab);
    memset(in, 0, sizeof(in));
    int peak = 0;
    for (int k = 0; k < 60; k++) {
        oa_enc_frame(&e2, in, bits);
        oa_dec_frame(&dec, bits, 0, out);
        if (k > 40)
            for (int n = 0; n < 240; n++) if (abs(out[n]) > peak) peak = abs(out[n]);
    }
    CHECK(peak < 200, "decoder did not settle on silence after garbage (peak %d)", peak);
}

int main(void)
{
    oa_tables_t tab;
    printf("OpenACELP unit tests\n");
    test_math();
    test_lsp();
    test_bits();
    test_pitch();
    make_tables(&tab);
    test_loopback(&tab);
    printf(failures ? "%d FAILURES\n" : "all tests passed\n", failures);
    return failures != 0;
}
