# TETRA codec implementation notes

OpenACELP implements the TETRA full-rate speech codec (ETSI EN 300 395-2,
4567 bit/s ACELP) without using any of ETSI's C code. This file records
where each design decision comes from, so the provenance of the
implementation stays traceable.

Sources used:

- **Prose:** [EN 300 395-2 V1.3.3 (2025-02)](https://www.etsi.org/deliver/etsi_en/300300_300399/30039502/01.03.03_60/en_30039502v010303p.pdf), clause 4. Clause numbers below refer to it.
- **Tables:** the numeric values of `CLSP_334.TAB` and `ENER_QUA.TAB` from the
  ETSI ZIP, in OpenACELP's own arrays in `src/oa_tables_tetra.c`. They are
  regenerated with `make regen-tables ETSI_ZIP=...`.
- **Black-box testing:** ETSI's reference `scoder`/`sdecoder` were compiled locally and
  used as test oracles. Crafted bitstreams were fed to them and their output
  compared with ours. Findings that the prose does not state are marked
  **[observed]**.
- **Reading the reference:** two decoder behaviours that testing could not pin
  down were later confirmed by reading the reference source. They are marked
  **[reference]**. They were re-implemented from the description of what the
  reference does. No reference code was copied or transcribed.
- **Bit-exact oracle:** to make the decoder reproduce the reference sample for
  sample, it was compared with an independent decoder that passes the
  EN 300 395-4 decoder test, through its public interface only: frames in,
  PCM out. Its published description of the reference behaviour was used as
  a guide to what to measure; every detail below was then
  settled by measurement. Details found this way are marked **[measured]**.

## Frame and bit layout (clause 4.2.2.7, table 3)

137 bits per 30 ms frame, MSB first per field:

| Field | Bits |
|---|---|
| LSP indices | 8 + 9 + 9 |
| Per subframe (x4) | pitch (8 in subframe 1, 5 otherwise), pulse 4 (3), pulse 3 (3), pulse 2 (3), pulse 1 (5), global sign (1), shift (1), gain (6) |

The ETSI test file format is 138 16-bit words per frame (word 0 = BFI,
then one bit per word, LSB). The tools accept that format.

## LP parameters (4.2.2.2, 4.2.2.3, 4.2.3.1.1)

- **Codebooks:** split VQ (3/3/4) in the cosine domain, Q15.
- **Initial past vector:** {30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000} in Q15.
- **[reference]** After the codebook lookup, a minimum distance is enforced
  across the two split boundaries, in the cosine domain (Q15):
  - q3/q4 by 917 (0.028, about 50 Hz near 1 kHz)
  - q6/q7 by 1245 (0.038, about 50 Hz near 1.6 kHz)

  A closer pair has half the shortfall (rounded down) added to the upper one and
  subtracted from the lower one. If the vector is then not strictly decreasing,
  the previous frame's LSPs are used instead, for the interpolation and as the
  next frame's previous vector.

  Every entry of the TETRA LSP codebooks is even, so the distance across a
  boundary is always even, and 916/1244 would give exactly the same vectors.

  Black-box testing first suggested "sort, then keep 50 Hz everywhere", which
  matches on most frames. It fails on signals that frequently trigger the
  fallback, such as slow sweeps.

  The encoder restricts its split-VQ search to index combinations the decoder
  accepts.
- **[measured] Interpolation:** the previous vector weighted 3/4, 1/2, 1/4
  (Q15) and the new one 1/4, 1/2, 3/4, the sum shifted right by 15 (rounded
  down); subframe 4 uses the new vector.
- **[measured] LSP->A:** the F1/F2 expansion of clause 4.2.2.2 in Q24 with
  saturating 32-bit adds. Each product of a Q24 coefficient with a Q15 LSP is
  formed from the coefficient's high 16 bits, the top bit of its next 15 bits
  and their remaining 14 bits: `2 hi q + b q + 2 floor(l q / 2^15)`, then
  doubled. The coefficients of A(z) are the sums rounded to Q12 (16 bits; a
  value beyond +-8 wraps).

  These were read back from the output of single-pulse frames: with no new
  excitation, the output decays freely through 1/A(z). About 93,000 such samples,
  with random LSP indices (rejected vectors included), match exactly.

## Adaptive codebook (4.2.2.4, 4.2.3.1.2)

**[observed]** Lag index mapping (consistent with all 256 + 32 indices tested):

- **Subframe 1:**
  - idx < 197: `T = (idx+2)/3 + 19`, `frac = idx - 3T + 58` (frac in {-1, 0, 1} thirds)
  - otherwise: `T = idx - 112`, frac 0
- **Subframes 2-4:**
  - `tmin = clamp(T1 - 5, 20, 134)` with T1 the integer lag of subframe 1
  - `i = (idx+2)/3 - 1`, `T = tmin + i`, `frac = idx - 2 - 3i`
  - This spans T1 - (5 2/3) ... T1 + (4 2/3) as clause 4.2.2.4 states.

The lag is `T + frac/3`. Interpolation uses a Hamming-windowed sinc
truncated at +-16 samples (32 taps, clause 4.2.2.4). The excitation is built
in place, so lags shorter than 60 repeat it. A window of
`0.54 + 0.46 cos(pi*t / 16)` reproduces the oracle's interpolated pulses.

**[measured]** The 32 taps are the windowed sinc at -15-mu ... 16-mu samples
(mu = 1/3 for lag T - 1/3, 2/3 for T + 1/3), rounded to Q15. Each output
sample is a saturating 32-bit sum of doubled products, rounded to 16 bits.

## Algebraic codebook (4.2.2.5, 4.2.3.1.3)

- **Pulses:**
  - Pulse 1: amplitude +1.4142 at 2*i.
  - Pulses 2/3/4: amplitudes -1/+1/-1 at 2+8i, 4+8i, 6+8i.
  - Shift adds 1 to all positions. Positions >= 60 are absent.
- **[observed]** Global sign bit 1 means all pulses inverted.
- **Shaping:** F(z) = A(z/0.75)/A(z/0.85), impulse response over 60 samples.
  If the integer lag T < 60, `f(n) += 0.8 f(n-T)`. The code vector is
  the sum of the 4 shifted, scaled copies of f.
- **[measured] Fixed point:**
  - The weighted coefficients are `a[i] g^i` with the powers of g as
    recursive rounded Q15 products, and each product rounded to Q12.
  - f is the response to the numerator's coefficients, computed with doubled
    products in a saturating 32-bit accumulator, shifted left by 3
    (saturating) and rounded to 16 bits: f(0) = 4096.
  - Sharpening multiplies by 0.8 in **Q13** (6554, which is 0.80005), shifts
    right by 13 (rounding down) and adds with 16-bit saturation, in place
    (recursively). Q15 26214 or 26215 fail: f(T) = 3277 after a pulse of 4096,
    and -60 becomes -49. This was read off frames whose code gain is large
    enough (after pumping the energies with high gain indices) that the
    excitation shows the code vector exactly; 4,762 of 4,762 sharpened samples
    agree.
  - Amplitudes in Q14: 23168 (1.4142 to Q12 precision, x4), -16384, 16384,
    -16384. Each code vector sample is the 32-bit sum of `amplitude x f`,
    negated by the sign bit, shifted right by 14 (rounding down) and reduced
    to its low 16 bits.

## Gains (4.2.2.6, 4.2.3.1.4)

- **Prediction gain:** `pg = sum h^2(i)`, i = 0..59, with h the impulse response of 1/A(z).
- **Energies:**
  - `Ea = log2((sum v^2 + 1) * pg)`
  - `Ef = log2(sum c^2 * pg)`
- **Prediction:** `Ep~ = 0.5 Ep' + 0.25 Ec' - 3`, `Ec~ = 0.25 Ep' + 0.5 Ec' - 3`
  (' = previous subframe), with the energies 0 at the start.
- **[observed]** Ep~ and Ec~ are floored at 0.
- **Quantized energies:** `E = R + E~` (R from the 64x2 table, Q8: column 0
  pitch, column 1 code **[observed]**), limited to 27 (pitch) and 25 (code).
- **[observed]** Gains: `gp = min(2^(0.5 (Ep + log2 40 - Ea)), 1.2)`,
  `gc = 2^(0.5 (Ec + log2 40 - Ef))`.
- **[measured] Fixed point** (energies in Q8; `oa_xgain_logs`, `oa_xgains`):
  - h is the response of 1/A(z) (the synthesis filter below) to an impulse of
    1024. pg, `1 + sum v^2` and `sum c^2` are plain (not doubled) saturating
    32-bit sums of squares.
  - pg and the adaptive-vector energy are normalized and their high 16 bits
    taken; the code-vector energy contributes its high 16 bits without
    normalization.
  - Each product of two mantissas (a doubled 16 x 16 product) goes through
    one log2: 33-entry table `round(32767 log2(1 + i/32))`, linear
    interpolation on the next 15 bits. The normalization shifts are subtracted
    (both for Ea, pg's for Ef), the fraction is truncated to Q8, and one offset,
    -4690 (Q8), accounts for log2 40 and the Q formats of both gains.
  - The predictions are `(2 E'p + E'c) >> 2` and `(E'p + 2 E'c) >> 2` minus 768.
  - `E - Ea`, read as Q9, goes through a 33-entry pow2 table
    (`round(16384 * 2^(i/32))`) with linear interpolation to give the gain in
    Q12 (exponents above 30 are limited to 30). The pitch gain is limited to
    4915 (1.2); the code gain keeps the low 16 bits.

  The offset was fitted on 2,397 subframes and then matches all 21,890
  code-gain and 15,645 pitch-gain observations, including the pumped
  high-gain streams.

## Reconstruction (4.2.3.1.5, 4.2.1)

`u = gp v + gc c`, `s = u - sum ai s(n-i)`, output = 2s with saturation.
There is no postfilter.

**[measured]** `u = (gp v + gc c + 2048) >> 12`, keeping the low 16 bits; u
is also the adaptive codebook's memory. The synthesis filter sums `x a0 - sum
a(i) s(n-i)` (plain products, Q12 coefficients) in a saturating 32-bit
accumulator, then adds 2048, shifts right by 12 and saturates to 16 bits; the
output doubles that with saturation. The saturation order was settled on
random-bit streams whose output clips: a doubled-product accumulator
(saturating 2x earlier) or a 32-bit filter memory both fail.

## Error concealment (4.2.3.2)

On BFI:

- Keep the previous LSPs (all four subframes use them).
- Use the last integer lag with frac 0 for all subframes.
- Reuse the previous frame's innovation indices.
- `E = E' - 0.5` per subframe, floored at 0 **[reference]**.
- **[measured]** Bad frames before the first good frame since initialization
  decode to silence and leave the state alone, except that the first LSP of
  the previous vector becomes 0.

Any nonzero bad-frame indicator in the ETSI file format marks a bad frame.

## Encoder

The encoder follows clause 4.2.2:

- **LP analysis:** 256-sample asymmetric window, 60 Hz lag window, Levinson-Durbin.
- **LP->LSP:** Chebyshev root search, then 26-bit split VQ (MSE, full search).
- **Weighting and open-loop pitch:** perceptual weighting A(z/0.95)/A(z/0.6),
  then open-loop pitch with the three-range rule.
- **Closed-loop pitch:** on 1/Â(z/0.85), with 1/3-sample refinement.
- **Algebraic search:** 4-pulse focused search with k2 = k3 = 0.586 and the
  350-count limit, both shifts.
- **Gain VQ:** in the energy domain.

It is an independent implementation. Its own choices, which the bitstream
does not constrain, include:

- The LP->LSP conversion takes Q24 coefficients straight from Levinson-Durbin
  (not Q12), on a 120-point grid (the spec says 60), so the root search
  stays accurate for closely spaced LSPs.
- The gain index is chosen by evaluating every one of the 64 candidates
  with the decoder's own gain functions (`oa_xgain_quant`) and minimizing the
  weighted error. The adaptive and algebraic code vectors, the energies and
  the excitation history are computed with the decoder's primitives too
  (`oa_xdec.h`), so the encoder's local decoder is the decoder.
- Fractional pitch uses an 8-tap Hamming-windowed sinc on the normalized
  correlation, as described.

## Testing

- `make test`:
  - `test_unit`: unit tests with synthetic codebooks (no ETSI data). They
    cover the fixed-point math, the LSP round trip, bit packing, the pitch
    index bijection, encoder/decoder loopback, and decoder robustness on
    garbage and erased frames.
  - `test_tetra`: decodes seven deterministic streams with the TETRA
    codebooks (random bits, random and burst losses, losses before the first
    good frame, encoder output) and checks PCM hashes. Each stream was decoded
    by the bit-exact oracle when the hashes were recorded and gave the same PCM.
- `make conformance ETSI_TEST_DIR=dir`: the EN 300 395-4 decoder test.
- `tests/interop.py <etsi-bin> <speech.raw>...`: black-box comparison
  with the ETSI executables.
- `tools/ref_decoder.py`: float reference model of the decoder, useful
  for checking the structure of fixed-point changes.

### Building the oracle

The ETSI code assumes a 32-bit `long`. On 64-bit hosts, change
`typedef long Word32;` to `typedef int Word32;` in `source.h` and
`channel.h` of a scratch copy, then:

```
gcc -O2 -w -o scoder  scoder.c  scod_tet.c sub_sc_d.c sub_dsp.c fbas_tet.c fexp_tet.c fmat_tet.c tetra_op.c
gcc -O2 -w -o sdecoder sdecoder.c sdec_tet.c sub_sc_d.c sub_dsp.c fbas_tet.c fexp_tet.c fmat_tet.c tetra_op.c
```

These binaries are test tools only and must never be distributed or linked
into OpenACELP.

## Validation status

Decoder, compared sample for sample with the bit-exact oracle:

| Streams | Result |
|---|---|
| codec2 speech (hts1a, hts2a, ve9qrp_10s), 1 kHz tone and sweep, encoded by OpenACELP (previous and current encoder) | identical |
| ve9qrp_10s with 10% of frames lost | identical |
| 16,650 frames of encoded speech | identical |
| single-pulse decay frames, random LSPs (400 frames) | identical |
| random parameters, all pulse/lag/gain combinations (600 + 1,500 frames) | identical |
| pumped high-gain frames (3,000) | identical |
| 2,000 frames of uniform random bits | identical |
| 9,000 random-bit streams of 10 frames from reset (3,000 of them with 10% losses) | 8,995 identical |
| 300 random-bit streams of 100 frames | 299 identical |

The streams that still differ are random bits driving unstable LP filters
whose output is pinned at full scale for whole subframes; the decoder's
state then differs from the reference's without being visible in the
output. No such stream has been found among valid speech.

The decoder is free of undefined behaviour under GCC's ASan and UBSan on
about a million random, all-ones, all-zeros and loss-heavy frames.

Encoder, OpenACELP encoder -> oracle decoder (SNR / segmental SNR, after the
40-sample look-ahead; 1 kHz tone and 100-1000 Hz sweep at amplitude 8000):

| Signal | Previous encoder | Current encoder |
|---|---|---|
| hts1a | 10.03 / 5.96 dB | 10.03 / 5.97 dB |
| hts2a | 11.74 / 8.48 dB | 11.80 / 8.45 dB |
| ve9qrp_10s | 6.68 / 5.17 dB | 6.69 / 5.15 dB |
| 1 kHz tone | 25.75 / 29.07 dB | 27.32 / 32.78 dB |
| sweep | 20.20 / 23.88 dB | 20.46 / 24.60 dB |
