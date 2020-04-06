# OpenACELP
Free ACELP vocoder. It is based on **ETSI EN 300-395-2** and **TIA/EIA IS-641**, but right now it is *not* compatible with any of them (their codebooks can't be published as a part of this codec). **OpenACELP** is an alternative to (great) [Codec 2](https://github.com/drowe67/codec2).

I'm using TED-LIUM release 1 ([OpenSLR link](http://www.openslr.org/7/), [download](https://projets-lium.univ-lemans.fr/ted-lium/release1/)) as the english speech corpus and [py-lbg](https://github.com/internaut/py-lbg) for codebook generation using Linde-Buzo-Gray (LBG) algorithm.

**Actual phase:** preparing codebooks for LSP split-vector quantization.
