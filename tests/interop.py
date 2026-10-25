#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Black-box interoperability test against the ETSI reference executables.

    tests/interop.py <etsi-bin-dir> <speech.raw> [more.raw ...]

<etsi-bin-dir> must contain `scoder` and `sdecoder` built from the ETSI ZIP
(see docs/TETRA_NOTES.md, "Building the oracle"). The ETSI code is only run, never
linked or copied. Inputs are raw 16-bit little-endian 8 kHz PCM.

Reports, per input:
  - decoder match: OpenACELP vs ETSI decoder on the ETSI encoder's bitstream,
    clean and with frame erasures (BFI)
  - codec SNR for all four encoder/decoder combinations
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "..", "build")


def run(exe, *args):
    subprocess.run([exe, *args], check=True, stdin=subprocess.DEVNULL,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def raw(path):
    return np.fromfile(path, "<i2").astype(float)


def snr(ref, test, delay=0):
    n = min(len(ref), len(test) - delay)
    a, b = ref[:n], test[delay:delay + n]
    return 10 * np.log10((np.sum(a * a) + 1) / (np.sum((a - b) ** 2) + 1))


def seg_snr(ref, test, delay=0, thr=1e7):
    n = (min(len(ref), len(test) - delay) // 240) * 240
    a, b = ref[:n].reshape(-1, 240), test[delay:delay + n].reshape(-1, 240)
    s = [10 * np.log10(np.sum(p * p) / (np.sum((p - q) ** 2) + 1)) for p, q in zip(a, b) if np.sum(p * p) > thr]
    return float(np.mean(s)) if s else float("nan")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    etsi = sys.argv[1]
    scoder, sdecoder = os.path.join(etsi, "scoder"), os.path.join(etsi, "sdecoder")
    oa_enc, oa_dec = os.path.join(BUILD, "oa_encode"), os.path.join(BUILD, "oa_decode")
    tmp = tempfile.mkdtemp(prefix="oa_interop_")
    t = lambda name: os.path.join(tmp, name)
    worst = 99.0
    for wav in sys.argv[2:]:
        x = raw(wav)
        print(f"== {os.path.basename(wav)} ({len(x) / 8000:.1f} s)")
        run(scoder, wav, t("e.ser"), t("e.syn"))
        run(oa_enc, wav, t("o.ser"))

        # decoder match on the ETSI bitstream, clean and with erasures
        frames = np.fromfile(t("e.ser"), "<i2").reshape(-1, 138)
        for name, mask in (("clean", []), ("BFI 1/7", list(range(3, len(frames), 7))),
                           ("BFI bursts", [i for i in range(len(frames)) if i % 40 in (20, 21, 22)])):
            f = frames.copy()
            f[mask, 0] = 1
            f.astype("<i2").tofile(t("b.ser"))
            run(sdecoder, t("b.ser"), t("r.raw"))
            run(oa_dec, t("b.ser"), t("d.raw"))
            m = snr(raw(t("r.raw")), raw(t("d.raw")))
            worst = min(worst, m)
            print(f"   decoder match vs ETSI ({name:10s}): {m:5.1f} dB")

        # codec quality, all combinations (40-sample codec delay)
        for enc_name, ser in (("ETSI enc", "e.ser"), ("OA enc  ", "o.ser")):
            row = []
            for dec_name, dec in (("ETSI dec", sdecoder), ("OA dec", oa_dec)):
                run(dec, t(ser), t("y.raw"))
                y = raw(t("y.raw"))
                row.append(f"{dec_name}: SNR {snr(x, y, 40):5.2f} / segSNR {seg_snr(x, y, 40):5.2f} dB")
            print(f"   {enc_name} -> " + "   ".join(row))
    print(f"worst decoder match: {worst:.1f} dB")


if __name__ == "__main__":
    main()
