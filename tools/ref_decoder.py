# SPDX-License-Identifier: GPL-3.0-or-later
"""OpenACELP float reference model of the TETRA (EN 300 395-2) speech decoder.

Usage: tools/ref_decoder.py <src/oa_tables_tetra.c> <in.ser> <out.raw>


Written from the prose of EN 300 395-2 clause 4, plus behaviour established by
black-box testing against the ETSI executables (see docs/TETRA_NOTES.md).
"""
import re
import numpy as np

L = 60
K_REF = np.log2(40.0)          # energy reference offset (empirical)
Q_INIT = np.array([30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000]) / 32768.0


def load_tables(path):
    t = re.sub(r"/\*.*?\*/", " ", open(path).read(), flags=re.S)   # drop comments
    arr = {}
    for name in ["oa_tetra_lsp_cb1", "oa_tetra_lsp_cb2", "oa_tetra_lsp_cb3", "oa_tetra_gain_cb"]:
        body = re.search(name + r"\[[^\]]*\]\s*=\s*\{(.*?)\};", t, re.S).group(1)
        arr[name] = np.array([int(v) for v in body.replace("\n", " ").split(",") if v.strip()], float)
    return (arr["oa_tetra_lsp_cb1"].reshape(-1, 3) / 32768, arr["oa_tetra_lsp_cb2"].reshape(-1, 3) / 32768,
            arr["oa_tetra_lsp_cb3"].reshape(-1, 4) / 32768, arr["oa_tetra_gain_cb"].reshape(-1, 2) / 256)


def unpack(bits137):
    pos = 0
    def take(n):
        nonlocal pos
        v = 0
        for k in range(n):
            v = (v << 1) | (int(bits137[pos + k]) & 1)
        pos += n
        return v
    p = dict(lsp=(take(8), take(9), take(9)), sub=[])
    for sf in range(4):
        s = dict(pitch=take(8 if sf == 0 else 5))
        s["p4"], s["p3"], s["p2"], s["p1"] = take(3), take(3), take(3), take(5)
        s["sign"], s["shift"], s["gain"] = take(1), take(1), take(6)
        p["sub"].append(s)
    return p


def lsp_rule(q, prev):
    """minimum distance across the codebook boundaries; keep prev if not decreasing"""
    q = np.round(q * 32768).astype(int)
    for hi, lo, dmin in ((2, 3, 917), (5, 6, 1245)):
        t = dmin - (q[hi] - q[lo])
        if t > 0:
            t >>= 1
            q[hi] += t; q[lo] -= t
    if np.any(np.diff(q) >= 0):
        return prev.copy()
    return q / 32768.0


def lsp_to_a(q):
    def poly(qs):
        f = np.zeros(6); f[0] = 1.0; f[1] = -2 * qs[0]
        for i in range(2, 6):
            b = qs[i - 1]
            f[i] = -2 * b * f[i - 1] + 2 * f[i - 2]
            for j in range(i - 1, 0, -1):
                f[j] += -2 * b * f[j - 1] + (f[j - 2] if j >= 2 else 0)
        return f
    f1, f2 = poly(q[0::2]), poly(q[1::2])
    for i in range(5, 0, -1):
        f1[i] += f1[i - 1]; f2[i] -= f2[i - 1]
    a = np.zeros(11); a[0] = 1
    for i in range(1, 6):
        a[i] = 0.5 * (f1[i] + f2[i]); a[11 - i] = 0.5 * (f1[i] - f2[i])
    return a


def iir_impulse(num, den, n=L):
    y = np.zeros(n)
    for k in range(n):
        acc = num[k] if k < len(num) else 0.0
        for i in range(1, min(k, 10) + 1):
            acc -= den[i] * y[k - i]
        y[k] = acc
    return y


def interp_taps(mu, half=16):
    k = np.arange(-half + 1, half + 1)
    t = k - mu
    return k, np.sinc(t) * (0.54 + 0.46 * np.cos(np.pi * t / half))


def pitch_first(idx):
    if idx < 197:
        ti = (idx + 2) // 3 + 19
        return ti, idx - 3 * ti + 58
    return idx - 112, 0


def pitch_rel(idx, tmin):
    i = (idx + 2) // 3 - 1
    return i + tmin, idx - 2 - 3 * i


class Decoder:
    def __init__(self, tables):
        self.cb1, self.cb2, self.cb3, self.gcb = tables
        self.qprev = Q_INIT.copy()
        self.exc = np.zeros(154 + 20 + L)
        self.mem = np.zeros(10)
        self.Ep = self.Ec = 0.0
        self.last = None          # last good parameters
        self.last_T = 40

    def frame(self, frame138):
        bfi = int(frame138[0]) & 1
        if not bfi:
            p = unpack(frame138[1:])
            self.last = p
            return self._synth(p, conceal=False)
        if self.last is None:
            self.last = unpack(np.zeros(137, int))
        return self._synth(self.last, conceal=True)

    def _synth(self, p, conceal):
        if conceal:
            q = self.qprev
        else:
            i1, i2, i3 = p["lsp"]
            q = lsp_rule(np.concatenate([self.cb1[i1], self.cb2[i2], self.cb3[i3]]), self.qprev)
        qs = [0.75 * self.qprev + 0.25 * q, 0.5 * self.qprev + 0.5 * q, 0.25 * self.qprev + 0.75 * q, q]
        self.qprev = q
        out = np.zeros(4 * L)
        tmin = 20
        for sf, s in enumerate(p["sub"]):
            a = lsp_to_a(qs[sf])
            if conceal:
                ti, fr = self.last_T, 0
            elif sf == 0:
                ti, fr = pitch_first(s["pitch"])
                tmin = min(max(ti - 5, 20), 134)
            else:
                ti, fr = pitch_rel(s["pitch"], tmin)
            self.last_T = ti
            # adaptive codebook (in place: repeats excitation for short lags)
            buf = np.concatenate([self.exc[L:], np.zeros(L)])
            st = len(buf) - L
            T = ti + fr / 3.0
            for n in range(L):
                x = st + n - T
                x0 = int(np.floor(x)); k, h = interp_taps(x - x0)
                buf[st + n] = np.dot(buf[x0 + k], h)
            v = buf[st:].copy()
            # algebraic code through the shaping filter F(z) = A(z/0.75)/A(z/0.85)
            g = np.arange(11)
            f = iir_impulse(a * 0.75 ** g, a * 0.85 ** g)
            if ti < L:
                for n in range(ti, L):
                    f[n] += 0.8 * f[n - ti]
            sign = -1.0 if s["sign"] else 1.0
            c = np.zeros(L)
            for m, amp in ((2 * s["p1"], 1.4142), (2 + 8 * s["p2"], -1.0), (4 + 8 * s["p3"], 1.0), (6 + 8 * s["p4"], -1.0)):
                m += s["shift"]
                if m < L:
                    c[m:] += sign * amp * f[:L - m]
            # gains in the energy domain
            pg = np.sum(iir_impulse([1.0], a) ** 2)
            Ea = np.log2((np.sum(v * v) + 1.0) * pg)
            Ef = np.log2(np.sum(c * c) * pg + 1e-30)
            if conceal:
                Ep, Ec = max(self.Ep - 0.5, 0.0), max(self.Ec - 0.5, 0.0)
            else:
                Rp, Rc = self.gcb[s["gain"]]
                Ep = min(Rp + max(0.5 * self.Ep + 0.25 * self.Ec - 3.0, 0.0), 27.0)
                Ec = min(Rc + max(0.25 * self.Ep + 0.5 * self.Ec - 3.0, 0.0), 25.0)
            self.Ep, self.Ec = Ep, Ec
            gp = min(2 ** (0.5 * (Ep + K_REF - Ea)), 1.2)
            gc = 2 ** (0.5 * (Ec + K_REF - Ef))
            u = gp * v + gc * c
            buf[st:] = u
            self.exc = buf
            mem = self.mem
            for n in range(L):
                y = u[n] - np.dot(a[1:], mem)
                mem = np.concatenate([[y], mem[:-1]])
                out[sf * L + n] = y
            self.mem = mem
        return np.clip(np.round(2 * out), -32768, 32767)


def decode_file(tables, ser):
    fr = np.fromfile(ser, "<i2").reshape(-1, 138)
    d = Decoder(tables)
    return np.concatenate([d.frame(f) for f in fr])


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    decode_file(load_tables(sys.argv[1]), sys.argv[2]).astype("<i2").tofile(sys.argv[3])
