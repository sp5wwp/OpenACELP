# OpenACELP
Free ACELP vocoder. It is based on **ETSI EN 300-395-2**<sup>[1]</sup> and **TIA/EIA IS-641**<sup>[2]</sup>, but right now it is *not* compatible with any of them (their codebooks can't be published as a part of this codec). **OpenACELP** is an alternative to (great) [Codec 2](https://github.com/drowe67/codec2). It uses floating point arithmetic. I aim to optimize it for the STM32 Cortex-M7, and they have a hardware floating point unit (FPU).

I'm using TED-LIUM release 1 ([OpenSLR link](http://www.openslr.org/7/), [download](https://projets-lium.univ-lemans.fr/ted-lium/release1/)) as the english speech corpus and [py-lbg](https://github.com/internaut/py-lbg) for codebook generation using Linde-Buzo-Gray (LBG) algorithm.

**Actual phase:** weighted speech signal computation.

Done:
- speech framing and windowing (refer to [2], chapter 2.2.1)
- autocorrelation parameters of the windowed speech, bandwidth expanded (equation 2.4 in [2])
- Levinson-Durbin recursive solver for computing Linear Prediction (LP) filter coefficients ([2] - chapter 2.2.2)
- conversion of LP coefficients to the LSPs (Line Spectral Pair) in the cosine domain ([2] - chapter 2.2.3)
- Chebyshev polynomials generation for LSPs root search (LSP polynomial evaluation)
- prepared test codebooks for LSP split-vector quantization (LSPs in the cosine domain). Refer to [1], chapter 4.2.2.3.
- LSP codebooks full search
- quantized and unquantized LSPs interpolation for each subframe
- LSP->LP conversion
