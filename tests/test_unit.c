/* OpenACELP unit tests. Self-contained: uses synthetic codebooks, no ETSI data.
 * (tests/test_tetra.c checks the decoder with the TETRA codebooks.)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "openacelp.h"
#include "oa_dsp.h"
#include "oa_basicop.h"
#include "oa_xdec.h"

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

/* ---- saturating fixed-point operators against 64-bit references */
static int64_t clamp64(int64_t x, int64_t lo, int64_t hi) { return x < lo ? lo : (x > hi ? hi : x); }
static int32_t rnd32(void)
{
    static const int32_t edge[] = { 0, 1, -1, 32767, -32768, INT32_MAX, INT32_MIN, 0x40000000, -0x40000000, 0x3FFFFFFF };
    uint32_t r = rnd();
    if (r % 8 == 0) return edge[rnd() % 10];
    return (int32_t)((rnd() << 16) ^ rnd()) >> (rnd() % 24);
}
static int16_t rnd16(void)
{
    static const int16_t edge[] = { 0, 1, -1, 32767, -32768, 16384, -16384 };
    return (rnd() % 8 == 0) ? edge[rnd() % 7] : (int16_t)rnd();
}
static void test_basicop(void)
{
    for (int t = 0; t < 200000; t++) {
        int32_t a = rnd32(), b = rnd32();
        int16_t x = rnd16(), y = rnd16();
        int n = (int)(rnd() % 34);
        CHECK(bo_L_add(a, b) == clamp64((int64_t)a + b, INT32_MIN, INT32_MAX), "L_add(%d, %d)", a, b);
        CHECK(bo_L_sub(a, b) == clamp64((int64_t)a - b, INT32_MIN, INT32_MAX), "L_sub(%d, %d)", a, b);
        CHECK(bo_sat16(a) == clamp64(a, -32768, 32767), "sat16(%d)", a);
        CHECK(bo_wrap16(a) == (int16_t)(uint16_t)(uint32_t)a, "wrap16(%d)", a);
        CHECK(bo_L_mult(x, y) == clamp64(2 * (int64_t)x * y, INT32_MIN, INT32_MAX), "L_mult(%d, %d)", x, y);
        int64_t m = clamp64(2 * (int64_t)x * y, INT32_MIN, INT32_MAX);
        CHECK(bo_L_mac(a, x, y) == clamp64(a + m, INT32_MIN, INT32_MAX), "L_mac");
        CHECK(bo_L_msu(a, x, y) == clamp64(a - m, INT32_MIN, INT32_MAX), "L_msu");
        CHECK(bo_mult_r(x, y) == clamp64(((int32_t)x * y + 0x4000) >> 15, -32768, 32767), "mult_r");
        CHECK(bo_L_shl(a, n) == clamp64((int64_t)a * ((int64_t)1 << (n > 31 ? 31 : n)) * (n > 31 ? 2 : 1), INT32_MIN, INT32_MAX)
              || a == 0, "L_shl(%d, %d)", a, n);
        int64_t rr = n == 0 ? a : ((int64_t)a + ((int64_t)1 << (n - 1))) >> n;
        CHECK(n > 31 || bo_L_shr_r(a, n) == rr, "L_shr_r(%d, %d)", a, n);
        CHECK(bo_round(a) == (int16_t)(clamp64((int64_t)a + 0x8000, INT32_MIN, INT32_MAX) >> 16), "round(%d)", a);
        if (a > 0) {
            int k = bo_norm_l(a);
            CHECK(((int64_t)a << k) >= 0x40000000 && ((int64_t)a << k) <= INT32_MAX, "norm_l(%d) = %d", a, k);
        }
    }
}

/* ---- LSP spreading, rejection and interpolation */
static void test_lsp_steps(void)
{
    /* the boundary rule: half the shortfall each way, rounded down; even inputs
     * give the same vector for gaps 916 and 917 */
    for (int t = 0; t < 20000; t++) {
        int16_t hi = (int16_t)(2 * (int)(rnd() % 16000) - 16000), lo = (int16_t)(hi - 2 * (int)(rnd() % 600));
        int16_t h1 = hi, l1 = lo, h2 = hi, l2 = lo;
        oa_lsp_spread(&h1, &l1, 917);
        oa_lsp_spread(&h2, &l2, 916);
        int32_t d = hi - lo, sh = d < 917 ? (917 - d) >> 1 : 0;
        CHECK(h1 == hi + sh && l1 == lo - sh, "spread (%d, %d)", hi, lo);
        CHECK(h1 == h2 && l1 == l2, "gaps 916 and 917 differ on even LSPs (%d, %d)", hi, lo);
    }
    /* a vector that stays crossed after spreading (here inside the second split)
     * is replaced by the previous one */
    static int16_t c1[3] = { 30000, 26000, 20000 }, c2[3] = { 9000, 15000, 1000 }, c3[4] = { -9000, -15000, -21000, -26000 };
    static int16_t gz[2] = { 0, 0 };
    oa_tables_t t = { c1, c2, c3, gz };
    const uint16_t idx[3] = { 0, 0, 0 };
    int16_t prev[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }, q[10];
    CHECK(oa_lsp_decode(&t, idx, prev, q) == -1 && !memcmp(q, prev, sizeof(q)), "crossed vector not rejected");
    c2[0] = 18000;
    c2[1] = 9000;
    CHECK(oa_lsp_decode(&t, idx, prev, q) == 0 && q[2] - q[3] == 2000, "ordered vector rejected");
    c2[0] = 19900;                                    /* 100 apart: pushed to 917 - 1 apart */
    CHECK(oa_lsp_decode(&t, idx, prev, q) == 0 && q[2] == 20408 && q[3] == 19492, "boundary spreading %d %d", q[2], q[3]);
    /* interpolation: 3/4, 1/2, 1/4 of the previous vector, rounded down */
    for (int k = 0; k < 2000; k++) {
        int16_t p[10], c[10], o[4][10];
        for (int i = 0; i < 10; i++) { p[i] = rnd16(); c[i] = rnd16(); }
        oa_lsp_interp(p, c, o);
        for (int i = 0; i < 10; i++) {
            CHECK(o[0][i] == (int16_t)floor((3.0 * p[i] + c[i]) / 4), "interp 1");
            CHECK(o[1][i] == (int16_t)floor((p[i] + (double)c[i]) / 2), "interp 2");
            CHECK(o[2][i] == (int16_t)floor((p[i] + 3.0 * c[i]) / 4), "interp 3");
            CHECK(o[3][i] == c[i], "interp 4");
        }
    }
}

/* ---- LSPs -> A(z) in Q12 against double precision */
static void test_lsp_to_a(void)
{
    int worst = 0;
    for (int t = 0; t < 3000; t++) {
        double w[10];
        for (int i = 0; i < 10; i++) w[i] = rndu();
        for (int i = 0; i < 10; i++) for (int j = i + 1; j < 10; j++) if (w[j] < w[i]) { double x = w[i]; w[i] = w[j]; w[j] = x; }
        int16_t q[10], a[11];
        for (int i = 0; i < 10; i++) q[i] = (int16_t)lrint(cos(0.05 + 3.04 * w[i]) * 32767);
        double f1[6] = { 1 }, f2[6] = { 1 };
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
        oa_xlsp_to_a(q, a);
        CHECK(a[0] == 4096, "a[0]");
        for (int i = 1; i <= 5; i++) {
            double ai = 0.5 * (f1[i] + f2[i]) * 4096, aj = 0.5 * (f1[i] - f2[i]) * 4096;
            if (fabs(ai) < 32000 && fabs(aj) < 32000) {
                int e1 = abs(a[i] - (int)lrint(ai)), e2 = abs(a[11 - i] - (int)lrint(aj));
                if (e1 > worst) worst = e1;
                if (e2 > worst) worst = e2;
            }
        }
    }
    CHECK(worst <= 2, "A(z) worst error %d (Q12)", worst);
}

/* ---- synthesis filter: the plain-arithmetic fast path equals the saturating sum */
static void syn_ref(const int16_t a[11], const int16_t *x, int16_t *y, int n, int16_t mem[10])
{
    int16_t buf[10 + 60];
    memcpy(buf, mem, 20);
    for (int k = 0; k < n; k++) {
        int64_t s = (int64_t)x[k] * a[0];
        for (int i = 1; i <= 10; i++)
            s = clamp64(s - (int64_t)a[i] * buf[10 + k - i], INT32_MIN, INT32_MAX);
        s = clamp64(s + 2048, INT32_MIN, INT32_MAX);
        buf[10 + k] = (int16_t)clamp64(s >> 12, -32768, 32767);
    }
    memcpy(y, &buf[10], n * 2);
    memcpy(mem, &buf[n], 20);
}
static void test_synthesis(void)
{
    for (int t = 0; t < 4000; t++) {
        int16_t a[11], x[60], y1[60], y2[60], m1[10], m2[10];
        int scale = (int)(rnd() % 4);                     /* small coefficients take the fast path */
        a[0] = 4096;
        for (int i = 1; i <= 10; i++) a[i] = (int16_t)((int16_t)rnd() >> (scale + 2 + (int)(rnd() % 4)));
        for (int i = 0; i < 60; i++) x[i] = (int16_t)((int16_t)rnd() >> (rnd() % 12));
        for (int i = 0; i < 10; i++) m1[i] = m2[i] = (int16_t)((int16_t)rnd() >> (rnd() % 12));
        oa_xsyn(a, x, y1, 60, m1);
        syn_ref(a, x, y2, 60, m2);
        CHECK(!memcmp(y1, y2, sizeof(y1)) && !memcmp(m1, m2, sizeof(m1)), "synthesis differs from the saturating reference");
    }
}

/* ---- algebraic code vector: amplitudes, positions, sign, Q13 sharpening */
static void test_code(void)
{
    const int16_t flat[11] = { 4096, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };   /* F(z) = 1: f = impulse of 4096 */
    int16_t c[60];
    uint8_t pos[4];
    oa_pulse_positions(5, 7, 7, 7, 0, pos);                             /* pulse 1 at 10, pulse 2 at 58 */
    CHECK(pos[0] == 10 && pos[1] == 58 && pos[2] >= 60 && pos[3] >= 60, "pulse positions");
    oa_xcode(flat, 100, pos, 0, c);
    CHECK(c[10] == 5792 && c[58] == -4096, "pulse amplitudes %d %d", c[10], c[58]);
    int nz = 0;
    for (int n = 0; n < 60; n++) nz += c[n] != 0;
    CHECK(nz == 2, "flat F(z) gives %d nonzero samples", nz);
    oa_xcode(flat, 100, pos, 1, c);
    CHECK(c[10] == -5792 && c[58] == 4096, "sign bit");
    oa_pulse_positions(0, 7, 7, 7, 0, pos);                             /* pulse 1 at 0 */
    oa_xcode(flat, 20, pos, 0, c);
    /* f(20) = 4096 * 6554 >> 13 = 3277, f(40) = 3277 * 6554 >> 13 = 2621; c = 1.4142 f */
    CHECK(c[0] == 5792 && c[20] == (int16_t)((23168 * 3277) >> 14) && c[40] == (int16_t)((23168 * 2621) >> 14),
          "sharpening %d %d %d", c[0], c[20], c[40]);
}

/* ---- gains and excitation */
static void test_gains(void)
{
    int16_t gp, gc;
    oa_xgains(0, 0, -12 * 512, -12 * 512, &gp, &gc);                  /* 2^(12/2 * 2)... 2^12 = 1.0 in Q12 */
    CHECK(gp == 4096 && gc == 4096, "unit gains %d %d", gp, gc);
    oa_xgains(512, 1024, -12 * 512, -12 * 512, &gp, &gc);
    CHECK(gp == 4915 && gc == 16384, "pitch gain limit / code gain 4.0: %d %d", gp, gc);
    oa_xgains(0, 2048, 0, -12 * 512, &gp, &gc);                        /* code gain 16.0 keeps its low 16 bits */
    CHECK(gp == 1 && gc == 0, "code gain wrap %d %d", gp, gc);
    int16_t prev = 0;
    for (int x = -8000; x < 1000; x += 7) {                            /* monotonic in the energy */
        oa_xgains(x, 0, 0, 0, &gp, &gc);
        CHECK(gp >= prev, "pitch gain not monotonic at %d", x);
        prev = gp;
    }
    /* energy prediction: floored at zero, limits 27 and 25 */
    int16_t ep, ec;
    oa_xgain_quant(0, 0, 100, 200, 0, 0, &ep, &ec, &gp, &gc);
    CHECK(ep == 100 && ec == 200, "prediction from zero %d %d", ep, ec);
    oa_xgain_quant(6912, 6400, 3000, 3000, 0, 0, &ep, &ec, &gp, &gc);
    CHECK(ep == 27 * 256 && ec == 25 * 256, "energy limits %d %d", ep, ec);
    oa_xgain_quant(4000, 2000, 0, 0, 0, 0, &ep, &ec, &gp, &gc);
    CHECK(ep == ((2 * 4000 + 2000) >> 2) - 768 && ec == ((4000 + 2 * 2000) >> 2) - 768, "prediction %d %d", ep, ec);

    int16_t v[60], c[60], u[60];
    for (int n = 0; n < 60; n++) { v[n] = rnd16(); c[n] = rnd16(); }
    oa_xexcite(v, c, 4915, -32768, u);
    for (int n = 0; n < 60; n++)
        CHECK(u[n] == (int16_t)(uint16_t)(((int64_t)4915 * v[n] - 32768LL * c[n] + 2048) >> 12), "excitation keeps low 16 bits");
}

/* ---- lost frames, and the encoder's local decoder against the decoder */
static void test_frames(const oa_tables_t *tab)
{
    static oa_dec_t d;
    uint8_t bits[OA_FRAME_BITS] = { 0 };
    int16_t pcm[OA_FRAME_SAMPLES];
    oa_dec_init(&d, tab);
    for (int k = 0; k < 3; k++) {
        oa_dec_frame(&d, bits, 1, pcm);
        int z = 1;
        for (int i = 0; i < OA_FRAME_SAMPLES; i++) z &= pcm[i] == 0;
        CHECK(z, "lost frame before the first good frame is not silent");
    }
    CHECK(d.q_prev[0] == 0 && d.q_prev[1] == 26000 && !d.have_last && d.e_pitch == 0, "state after a bad start");

    /* concealment lowers both energies by 0.5 (Q8 128) per subframe, not below zero */
    for (int i = 0; i < OA_FRAME_BITS; i++) bits[i] = rnd() & 1;
    oa_dec_frame(&d, bits, 0, pcm);
    d.e_pitch = 600;
    d.e_code = 300;
    oa_dec_frame(&d, bits, 1, pcm);
    CHECK(d.e_pitch == 600 - 4 * 128 && d.e_code == 0, "concealment energies %d %d", d.e_pitch, d.e_code);

    /* the encoder tracks the decoder's excitation and energies exactly */
    static oa_enc_t e;
    oa_enc_init(&e, tab);
    oa_dec_init(&d, tab);
    int16_t in[OA_FRAME_SAMPLES];
    double ph = 0;
    int same = 1;
    for (int k = 0; k < 60; k++) {
        for (int n = 0; n < OA_FRAME_SAMPLES; n++) {
            ph += 2 * M_PI * (120 + 30 * sin(k * 0.2)) / 8000;
            in[n] = (int16_t)lrint(5000 * (sin(ph) + 0.5 * sin(3 * ph) + 0.25 * sin(7 * ph)) + (int)(rnd() % 400) - 200);
        }
        oa_enc_frame(&e, in, bits);
        oa_dec_frame(&d, bits, 0, pcm);
        same &= !memcmp(e.xexc, d.exc, sizeof(d.exc)) && e.e_pitch == d.e_pitch && e.e_code == d.e_code
                && !memcmp(e.q_prev, d.q_prev, sizeof(d.q_prev));
    }
    CHECK(same, "encoder's local decoder state differs from the decoder's");
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
    test_basicop();
    test_lsp();
    test_lsp_steps();
    test_lsp_to_a();
    test_synthesis();
    test_code();
    test_gains();
    test_bits();
    test_pitch();
    make_tables(&tab);
    test_frames(&tab);
    test_loopback(&tab);
    printf(failures ? "%d FAILURES\n" : "all tests passed\n", failures);
    return failures != 0;
}
