# OpenACELP

A free, fixed-point ACELP speech codec that interoperates with the TETRA
full-rate codec ([ETSI EN 300 395-2](https://www.etsi.org/deliver/etsi_en/300300_300399/30039502/01.03.03_60/en_30039502v010303p.pdf), 4567 bit/s, 137 bits per 30 ms frame).
It contains no ETSI code. It was written from the standard's prose and
checked against the ETSI reference executables. Two details the standard  leaves
out were confirmed by reading the reference and re-implemented. See
[TETRA_NOTES](docs/TETRA_NOTES.md) for where every detail comes from.

- **Encoder and decoder:** both included, with error concealment for lost frames.
- **Arithmetic:** fixed point throughout (int16/int32 data with 64-bit
  accumulators), C17, no libm, no dynamic allocation. State structs can be
  allocated statically: encoder 19.7 KB (16 KB of it search scratch), decoder 1.1 KB.
- **Speed:** on x86-64, the encoder runs 1.4x and the decoder 2.4x faster than
  the reference. Code size is about 19 KB.

## Tables

The codec needs two numeric tables, the LSP and gain codebooks. They are ETSI
material, so they are not in this repository:

```
make tetra-tables                    # downloads EN 300 395-2 from etsi.org, verifies SHA-256
make tetra-tables ETSI_ZIP=file.zip  # or use a ZIP if you already have it
```

This writes `build/tetra_tables.c`, which defines `oa_tables_tetra`. The file
is git-ignored.

>[!Warning]
> **Do not commit the tables, and do not distribute binaries built with them!**

## Build

```
make          # library build/libopenacelp.a + unit tests (no ETSI data needed)
make test     # run unit tests
make apps     # build/oa_encode, build/oa_decode (needs the tables)
```

`oa_encode in.raw out.ser` and `oa_decode in.ser out.raw` use raw 16-bit 8 kHz
PCM and the ETSI test-file bitstream format (138 16-bit words per frame:
BFI + 137 bits). With `-p` they use 18-byte packed frames.

## API

```c
#include "openacelp.h"
extern const oa_tables_t oa_tables_tetra;

static oa_enc_t enc;  oa_enc_init(&enc, &oa_tables_tetra);
static oa_dec_t dec;  oa_dec_init(&dec, &oa_tables_tetra);

uint8_t bits[OA_FRAME_BITS];            /* 137 bits, one per byte */
oa_enc_frame(&enc, pcm_in, bits);       /* 240 samples in */
oa_dec_frame(&dec, bits, bfi, pcm_out);
```

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
