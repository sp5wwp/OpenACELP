# Codebook training scripts

These scripts are for training your own codebooks instead of using the ETSI
tables. Codebooks you train yourself do **not** interoperate with TETRA: only
the codebooks in `src/oa_tables_tetra.c` do. They are useful for a free,
non-ETSI mode of the codec, or for experiments.

| Script | Purpose |
|---|---|
| `sound_process.sh` | Turns a directory of TED-LIUM `.sph` files into one training corpus (`corpus.raw`): 8 kHz, 16-bit, normalized, low-passed at 3.4 kHz, with a little noise added. Needs `sox` and `ffmpeg`. |
| `q_codebook_generator.py` | Runs Linde-Buzo-Gray training (via [py-lbg](https://github.com/internaut/py-lbg)) for the three LSP split codebooks (3/3/4 dimensions). It expects the training vectors as Python modules `Q1`, `Q2`, `Q3`. |
| `lbg.c` | Standalone C implementation and test of LBG, for training on large data sets. |

## What is missing

The codec loads its tables through `oa_tables_t` (see `include/openacelp.h`):

- three LSP codebooks of 256x3, 512x3 and 512x4 entries, cosine domain, Q15
- a 64x2 gain codebook of (R_pitch, R_code) prediction errors, log2 energy, Q8

Two steps are not covered by these scripts yet:

1. **LSP training vectors.** The 2020 encoder used to dump these and has been
   removed. They can be regenerated from the corpus with the new encoder's LP
   analysis. Store them as Q15 cosine-domain values, already in the 3/3/4 split.
2. **The gain codebook.** Its training data are the energy prediction errors
   (R_p, R_c), which depend on the codebooks in use. Train it iteratively:
   encode the corpus with a provisional gain table, collect the unquantized
   prediction errors, run LBG on them, and repeat.

The LBG output also has to be written as an `oa_tables_t` C source file,
like the one `tools/extract_etsi_tables.py` produces for the ETSI tables.
