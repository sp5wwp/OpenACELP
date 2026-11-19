/* OpenACELP: saturating 16/32-bit fixed-point operators
 *
 * The decoder reproduces the TETRA reference decoder's output exactly, so it
 * computes with the kind of operators fixed-point speech codecs
 * conventionally use: 16-bit data, 32-bit accumulators, saturation on
 * overflow. These are OpenACELP's own definitions of those usual semantics.
 *
 * Every operator works in 32 bits (no 64-bit arithmetic) and does not rely on
 * signed overflow. On Cortex-M4/M7 (__ARM_FEATURE_DSP) the saturating adds and
 * shifts map to QADD/QSUB/SSAT.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OA_BASICOP_H
#define OA_BASICOP_H

#include <stdint.h>

#if defined(__ARM_FEATURE_DSP) && __ARM_FEATURE_DSP
#include <arm_acle.h>
#define OA_HAVE_ARM_DSP 1
#endif

static inline int16_t bo_sat16(int32_t x)
{
#ifdef OA_HAVE_ARM_DSP
    return (int16_t)__ssat(x, 16);
#else
    return x > 32767 ? 32767 : (x < -32768 ? -32768 : (int16_t)x);
#endif
}

/* low 16 bits, two's complement */
static inline int16_t bo_wrap16(int32_t x)
{
    return (int16_t)(uint16_t)((uint32_t)x & 0xFFFFu);
}

static inline int32_t bo_L_add(int32_t a, int32_t b)
{
#ifdef OA_HAVE_ARM_DSP
    return __qadd(a, b);
#else
    int32_t r = (int32_t)((uint32_t)a + (uint32_t)b);
    if (((a ^ r) & (b ^ r)) < 0)
        r = a < 0 ? INT32_MIN : INT32_MAX;
    return r;
#endif
}

static inline int32_t bo_L_sub(int32_t a, int32_t b)
{
#ifdef OA_HAVE_ARM_DSP
    return __qsub(a, b);
#else
    int32_t r = (int32_t)((uint32_t)a - (uint32_t)b);
    if (((a ^ b) & (a ^ r)) < 0)
        r = a < 0 ? INT32_MIN : INT32_MAX;
    return r;
#endif
}

/* (a * b + 2^14) >> 15, saturated */
static inline int16_t bo_mult_r(int16_t a, int16_t b)
{
    return bo_sat16(((int32_t)a * b + 0x4000) >> 15);
}

/* 2 * a * b, saturated to 32 bits (only -32768 * -32768 overflows) */
static inline int32_t bo_L_mult(int16_t a, int16_t b)
{
    int32_t p = (int32_t)a * b;
    return p == 0x40000000 ? INT32_MAX : 2 * p;
}

/* acc + 2ab, acc - 2ab: both steps saturating (QDADD/QDSUB on Cortex-M4) */
static inline int32_t bo_L_mac(int32_t acc, int16_t a, int16_t b)
{
#ifdef OA_HAVE_ARM_DSP
    return __qadd(acc, __qdbl((int32_t)a * b));
#else
    return bo_L_add(acc, bo_L_mult(a, b));
#endif
}

static inline int32_t bo_L_msu(int32_t acc, int16_t a, int16_t b)
{
#ifdef OA_HAVE_ARM_DSP
    return __qsub(acc, __qdbl((int32_t)a * b));
#else
    return bo_L_sub(acc, bo_L_mult(a, b));
#endif
}

/* shift left by n >= 0 with saturation */
static inline int32_t bo_L_shl(int32_t x, int n)
{
    if (n <= 0)
        return x;
    if (n >= 31)
        return x == 0 ? 0 : (x > 0 ? INT32_MAX : INT32_MIN);
    int32_t lim = INT32_MAX >> n;
    if (x > lim)
        return INT32_MAX;
    if (x < -lim - 1)
        return INT32_MIN;
    return (int32_t)((uint32_t)x << n);
}

/* shift right by n with rounding to nearest (ties up); n < 0 shifts left, saturating */
static inline int32_t bo_L_shr_r(int32_t x, int n)
{
    if (n <= 0)
        return bo_L_shl(x, -n);
    if (n > 31)
        return 0;
    if (n == 31)
        return x >= 0x40000000 ? 1 : (x < -0x40000000 ? -1 : 0);
    return (x >> n) + ((x >> (n - 1)) & 1);
}

static inline int16_t bo_extract_h(int32_t x) { return (int16_t)(x >> 16); }

/* upper half of x + 2^15, saturated */
static inline int16_t bo_round(int32_t x)
{
    return bo_extract_h(bo_L_add(x, 0x8000));
}

/* left shift that normalizes x > 0 to [2^30, 2^31); 0 for x <= 0 */
static inline int bo_norm_l(int32_t x)
{
    if (x <= 0)
        return 0;
#if defined(__GNUC__)
    return __builtin_clz((uint32_t)x) - 1;
#else
    int n = 0;
    while (x < 0x40000000) {
        x <<= 1;
        n++;
    }
    return n;
#endif
}

#endif
