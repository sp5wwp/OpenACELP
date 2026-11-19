# OpenACELP

A free, fixed-point ACELP speech codec that interoperates with the TETRA
full-rate codec ([ETSI EN 300 395-2](https://www.etsi.org/deliver/etsi_en/300300_300399/30039502/01.03.03_60/en_30039502v010303p.pdf), 4567 bit/s, 137 bits per 30 ms frame).
It contains no ETSI code. It was written from the standard's prose and
black-box tests (bitstreams in, PCM out) against the ETSI reference
executables. Two details the standard leaves out were
confirmed by reading the reference and re-implemented. See
[TETRA_NOTES](docs/TETRA_NOTES.md) for where every detail comes from.

- **Decoder output matches the reference decoder sample for sample**: speech,
  frame losses (including losses before the first good frame) and random-bit
  stress streams. `make conformance` runs the EN 300 395-4 decoder test
  against ETSI's test files. See [TETRA_NOTES](docs/TETRA_NOTES.md#validation-status).
- **Encoder and decoder:** both included, with error concealment for lost frames.
  The encoder's local decoder runs the decoder's own arithmetic, so the two
  never drift apart.
- **Arithmetic:** fixed point throughout, C17, no libm, no dynamic allocation.
  The decoder uses 16/32-bit integers only (no 64-bit operations); on
  Cortex-M4/M7 its saturating arithmetic maps to the DSP instructions
  (QADD, QSUB, QDADD, QDSUB, SSAT). State structs can be allocated statically:
  encoder 20.1 KB (16 KB of it search scratch), decoder 536 bytes.
- **Speed (Cortex-M4):** the decoder needs about 80,000 instructions per 30 ms
  frame on average, 98,000 at most (about 3.3 MIPS peak; emulated M4, clang 18
  `-O2 -mcpu=cortex-m4`), with 1.2 KB of stack. Its code is about 10 KB plus
  9 KB of tables.

## Tables

The codec needs two numeric tables, the LSP and gain codebooks. For TETRA
interoperability their values must be those of EN 300 395-2. They are in
[`src/oa_tables_tetra.c`](src/oa_tables_tetra.c), in OpenACELP's own arrays,
names and comments: only the numbers come from the standard, no ETSI code or
text.

Maintainers can regenerate or check that file from the ETSI ZIP:

```
make regen-tables ETSI_ZIP=path/to/en_30039502v010303p0.zip
```

It verifies the ZIP's SHA-256, then runs `tools/extract_etsi_tables.py`. The
ZIP is available from [etsi.org](https://www.etsi.org/deliver/etsi_en/300300_300399/30039502/01.03.03_60/en_30039502v010303p0.zip).

## Build

```
make            # library build/libopenacelp.a + tests
make test       # unit tests and the decoder regression test
make apps       # build/oa_encode, build/oa_decode
make cortex-m4  # compile the library for Cortex-M4 (CROSS=arm-none-eabi-)
make conformance ETSI_TEST_DIR=dir   # EN 300 395-4 decoder test, see below
```

For OpenRTX or another firmware build, compile `src/oa_dec.c`, `src/oa_dsp.c`,
`src/oa_bits.c` and `src/oa_tables_tetra.c` (add `src/oa_enc.c` for the
encoder) with `include/` and `src/` on the include path. With
`-mcpu=cortex-m4` (or any target defining `__ARM_FEATURE_DSP`) the DSP
instructions are used automatically.

`oa_encode in.raw out.ser` and `oa_decode in.ser out.raw` use raw 16-bit 8 kHz
PCM and the ETSI test-file bitstream format (138 16-bit words per frame:
BFI + 137 bits). With `-p` they use 18-byte packed frames.

## API

```c
#include "openacelp.h"     /* declares oa_tables_tetra */

static oa_enc_t enc;  oa_enc_init(&enc, &oa_tables_tetra);
static oa_dec_t dec;  oa_dec_init(&dec, &oa_tables_tetra);

uint8_t bits[OA_FRAME_BITS];            /* 137 bits, one per byte */
oa_enc_frame(&enc, pcm_in, bits);       /* 240 samples in */
oa_dec_frame(&dec, bits, bfi, pcm_out);
```

## Conformance

EN 300 395-4 specifies the decoder test: `TE1_46.CHD` must decode to
`TE2_128.PCM` byte for byte. The test files are not distributed with OpenACELP;
put them in a directory and run `make conformance ETSI_TEST_DIR=dir`.

## Interop testing

`tests/interop.py <etsi-bin-dir> <speech.raw>...` compares OpenACELP with
the ETSI executables, both on decoder output and with all four
encoder/decoder pairings. The ETSI binaries are only run as test tools, see
[TETRA_NOTES](docs/TETRA_NOTES.md) for how to build them.

## Training your own codebooks

`scripts/` holds tools for training non-ETSI LSP codebooks with [LBG](https://en.wikipedia.org/wiki/Linde%E2%80%93Buzo%E2%80%93Gray_algorithm) on a
speech corpus, see its [readme file](scripts/README.md) for what they cover
and what is still missing.

License: GPL-3.0-or-later.
