# TETRA codec implementation notes

OpenACELP implements the TETRA full-rate speech codec (ETSI EN 300 395-2,
4567 bit/s ACELP) without using any of ETSI's C code. This file records
where each design decision comes from, so the provenance of the
implementation stays traceable.

Sources used:

- **Prose:** [EN 300 395-2 V1.3.3 (2025-02)](https://www.etsi.org/deliver/etsi_en/300300_300399/30039502/01.03.03_60/en_30039502v010303p.pdf), clause 4. Clause numbers below refer to it.
- **Tables:** the two numeric tables `CLSP_334.TAB` and `ENER_QUA.TAB` from the
  ETSI ZIP, extracted at build time (`make tetra-tables`) and never committed.
- **Black-box testing:** ETSI's reference `scoder`/`sdecoder` were compiled locally and
  used as test oracles. Crafted bitstreams were fed to them and their output
  compared with ours. Findings that the prose does not state are marked
  **[observed]**.
- **Reading the reference:** two decoder behaviours that testing could not pin
  down were later confirmed by reading the reference source. They are marked
  **[reference]**. They were re-implemented from the description of what the
  reference does. No reference code was copied or transcribed.

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
- **Interpolation:** linear in the q domain, 0.75/0.25, 0.5/0.5, 0.25/0.75, then 1.
- **LSP->A:** standard F1/F2 expansion as written in clause 4.2.2.2.
- **[reference]** After the codebook lookup, a minimum distance is enforced
  across the two split boundaries, in the cosine domain (Q15):
  - q3/q4 by 917 (0.028, about 50 Hz near 1 kHz)
  - q6/q7 by 1245 (0.038, about 50 Hz near 1.6 kHz)

  A closer pair has half the shortfall (rounded down) added to the upper one and
  subtracted from the lower one. If the vector is then not strictly decreasing,
  the previous frame's LSPs are used instead.

  Black-box testing first suggested "sort, then keep 50 Hz everywhere", which
  matches on most frames. It fails on signals that frequently trigger the
  fallback, such as slow sweeps (7 dB decoder match, now 38 dB).

  The encoder restricts its split-VQ search to index combinations the decoder
  accepts.

## Adaptive codebook (4.2.2.4, 4.2.3.1.2)

**[observed]** Lag index mapping (consistent with all 256 + 32 indices tested):

- **Subframe 1:**
  - idx < 197: `T = (idx+2)/3 + 19`, `frac = idx - 3T + 58` (frac ∈ {-1, 0, 1} thirds)
  - otherwise: `T = idx - 112`, frac 0
- **Subframes 2-4:**
  - `tmin = clamp(T1 - 5, 20, 134)` with T1 the integer lag of subframe 1
  - `i = (idx+2)/3 - 1`, `T = tmin + i`, `frac = idx - 2 - 3i`
  - This spans T1 - (5 2/3) ... T1 + (4 2/3) as clause 4.2.2.4 states.

The lag is `T + frac/3`. Interpolation uses a Hamming-windowed sinc
truncated at ±16 samples (32 taps, clause 4.2.2.4). The excitation is built
in place, so lags shorter than 60 repeat it. A window of
`0.54 + 0.46 cos(pi*t / 16)` reproduces the oracle's interpolated pulses
(normalized correlation 0.999).

## Algebraic codebook (4.2.2.5, 4.2.3.1.3)

- **Pulses:**
  - Pulse 1: amplitude +1.4142 at 2*i.
  - Pulses 2/3/4: amplitudes -1/+1/-1 at 2+8i, 4+8i, 6+8i.
  - Shift adds 1 to all positions. Positions ≥ 60 are absent.
- **[observed]** Global sign bit 1 means all pulses inverted.
- **Shaping:** F(z) = A(z/0.75)/A(z/0.85), impulse response over 60 samples.
  If the integer lag T < 60, `f(n) += 0.8 f(n-T)`. The code vector is
  the sum of the 4 shifted, scaled copies of f.

## Gains (4.2.2.6, 4.2.3.1.4)

- **Prediction gain:** `pg = Σ ĥ²(i)`, i = 0..59, with ĥ the impulse response of 1/Â(z).
- **Energies:**
  - `Ea = log2((Σ v² + 1) * pg)`
  - `Ef = log2(Σ c² * pg)`
- **Prediction:** `Ẽp = 0.5 Êp' + 0.25 Êc' - 3`, `Ẽc = 0.25 Êp' + 0.5 Êc' - 3`
  (' = previous subframe), with Ê initially 0.
- **[observed]** Ẽp and Ẽc are floored at 0 (steady-state and start-up
  behaviour of the oracle match only with this floor).
- **Quantized energies:** `Ê = R̂ + Ẽ` (R̂ from the 64x2 table, Q8: column 0
  pitch, column 1 code **[observed]**), limited to 27 (pitch) and 25 (code).
- **[observed]** Gains:
  - `gp = min(2^(0.5 (Êp + log2 40 - Ea)), 1.2)`
  - `gc = 2^(0.5 (Êc + log2 40 - Ef))`
  - The log2 40 reference offset was measured as 5.32 ± 0.01 from recovered
    oracle gains, and log2 40 = 5.3219 gives the best match (median per-frame
    SNR vs. oracle 47.7 dB, against ~43 dB for 5.30 and 5.34).
  - The decoder clamps gp at 1.2, the same bound the prose gives for the encoder.

## Reconstruction (4.2.3.1.5, 4.2.1)

`u = gp v + gc c`, `ŝ = u - Σ âi ŝ(n-i)`, output = 2ŝ with saturation.
There is no postfilter.

## Error concealment (4.2.3.2)

On BFI:

- Keep the previous LSPs (all four subframes use them).
- Use the last integer lag with frac 0 for all subframes.
- Reuse the previous frame's innovation indices.
- `Ê = Ê' - 0.5` per subframe, floored at 0 **[reference]**.

This matches the oracle (37-50 dB in concealed frames).

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
  with exactly the decoder's gain computation and minimizing the weighted
  error. The encoder therefore always agrees with the decoder.
- Fractional pitch uses an 8-tap Hamming-windowed sinc on the normalized
  correlation, as described.

## Testing

- `make test`: unit tests with synthetic codebooks (no ETSI data). They cover
  the fixed-point math, the LSP round trip, bit packing, the pitch index
  bijection, encoder/decoder loopback, and decoder robustness on garbage and
  erased frames.
- `tests/interop.py <etsi-bin> <speech.raw>...`: black-box comparison
  with the ETSI executables.
- `tools/ref_decoder.py`: float reference model of the decoder, useful
  for checking fixed-point changes (49 dB match with the C decoder).

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

Measured on codec2's test speech (hts1a, hts2a, ve9qrp_10s), plus a 1 kHz tone
and a 100-1000 Hz sweep:

| Test | Result |
|---|---|
| Decoder match vs ETSI decoder, ETSI bitstreams | speech 38.9-50.0 dB, tone 43.4 dB, sweep 38.4 dB, unchanged with frame erasures |
| ETSI encoder -> ETSI decoder (SNR) | 8.95 / 11.59 / 6.35 dB speech, 15.45 tone, 13.06 sweep |
| OpenACELP encoder -> ETSI decoder (SNR) | 10.03 / 11.74 / 6.68 dB speech, 24.93 tone, 15.06 sweep |
| OpenACELP encoder -> OpenACELP decoder | within 0.2 dB of the above |

The remaining decoder difference is numeric. The reference keeps the excitation
and ŝ as 16-bit integers, the pitch gain in Q12 and the code gain as an integer -
OpenACELP keeps 4 more fractional bits. Rounding like the reference changes the
match by less than 1 dB, so the extra precision is kept.
