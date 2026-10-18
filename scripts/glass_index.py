# -*- coding: utf-8 -*-
"""Parse real Schott Sellmeier (CD) data from the installed SCHOTT.AGF.

CD row layout, verified against the nd / vd the catalogue itself carries on
every NM row:   CD B1 C1 B2 C2 B3 C3 B4 C4   (pairs, not blocks).
load_agf() re-computes nd and vd from the coefficients and refuses any glass
that does not match its NM row, so a wrong column guess cannot slip through.
"""
import math
import sys

AGF = r"C:\Users\Lenovo\Documents\Zemax\Glasscat\SCHOTT.AGF"


def _n(BC, lam):
    B, C = BC
    l2 = lam * lam
    s = 1.0
    for b, c in zip(B, C):
        s = s + b * l2 / (l2 - c)
    return math.sqrt(s)


def load_agf(path=AGF, tol_nd=2e-3, tol_vd=0.5):
    """Return {name: (B[3], C[3])} for every record whose CD reproduces its NM."""
    lines = open(path, "r", encoding="latin-1").read().splitlines()
    out, rejected = {}, []
    name = None
    for line in lines:
        if line.startswith("NM "):
            p = line.split()
            name = p[1]
            try:
                nd_cat, vd_cat = float(p[4]), float(p[5])
            except (IndexError, ValueError):
                nd_cat = vd_cat = None
        elif line.startswith("CD ") and name is not None:
            v = [float(x) for x in line.split()[1:7]]
            B, C = v[0::2], v[1::2]
            if nd_cat is not None:
                nd = _n((B, C), 0.5875618)
                vd = (nd - 1.0) / (_n((B, C), 0.4861327) - _n((B, C), 0.6562725))
                if abs(nd - nd_cat) > tol_nd or abs(vd - vd_cat) > tol_vd:
                    rejected.append(name)
                    name = None
                    continue
            out[name] = (B, C)
            name = None
    return out, rejected


_CACHE = {}


def n_index(glass, lam_um, cat=None):
    """Index at wavelength lam (um). Accepts scalars or numpy arrays."""
    if cat is None:
        if "cat" not in _CACHE:
            _CACHE["cat"] = load_agf()[0]
        cat = _CACHE["cat"]
    B, C = cat[glass]
    l2 = lam_um ** 2
    s = 1.0
    for b, c in zip(B, C):
        s = s + b * l2 / (l2 - c)
    return s ** 0.5


if __name__ == "__main__":
    cat, rej = load_agf()
    print("glasses accepted: %d   rejected: %d" % (len(cat), len(rej)))
    if rej:
        print("rejected sample:", rej[:10])
    for g in ("N-BK7", "F2", "SF2", "SF5", "SF6", "SF10", "SF57", "N-SF5",
              "N-SF6", "N-SF11", "N-LAK9", "N-LAK22", "N-SK16", "N-SK2",
              "N-KZFS4", "N-KZFS11", "N-FK5", "SK16", "BALKN3", "N-LASF9",
              "N-BAK4", "N-BAF10", "N-SSK8", "N-F2", "N-SF1"):
        if g not in cat:
            print("%-10s MISSING" % g)
            continue
        nd = float(n_index(g, 0.5875618, cat))
        nF = float(n_index(g, 0.4861327, cat))
        nC = float(n_index(g, 0.6562725, cat))
        print("%-10s nd=%.6f  vd=%.2f" % (g, nd, (nd - 1) / (nF - nC)))
