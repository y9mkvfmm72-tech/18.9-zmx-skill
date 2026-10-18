# -*- coding: utf-8 -*-
"""Naive independent tracer for the Zemax stock doublet 3- doublet.zmx."""

import math
import glass_index as G

CAT = G.load_agf()[0]
N1 = float(G.n_index("N-BK7", 0.633, CAT))
N2 = float(G.n_index("N-SF5", 0.633, CAT))
C = [0.0, 0.016249593760156, -0.022527596305474, -0.007676953784738]
# The last surface of the real 3- doublet.zmx is CONCAVE (negative curvature,
# re-read from the file on 2026-10-18).  An earlier version of this script had
# +0.007676953784738; with that sign the doublet is nearly AFOCAL (EFL ~ -19 m)
# and gate A still PASSED, because gate A only compares two implementations fed
# the same wrong data.  That failure mode is exactly why check_design_gates now
# includes a known-answer check (gate E) against the actual file's numbers.
VERT = [0.0, 5.0, 8.0, 9.5]
BFD = 97.9387


def trace(y0, field_deg=0.0, image_dist=BFD):
    a = math.radians(field_deg)
    X, Y, Z = 0.0, y0, 0.0
    dx, dy, dz = 0.0, math.sin(a), math.cos(a)
    med = [(1.0, 1.0), (1.0, N1), (N1, N2), (N2, 1.0)]
    for i in (1, 2, 3):
        c = C[i]
        zv = VERT[i]
        t = (zv - Z) / dz
        while True:
            xx, yy = X + dx * t, Y + dy * t
            r2 = xx * xx + yy * yy
            qz = math.sqrt(1.0 - c * c * r2)
            s = c * r2 / (1.0 + qz)
            F = (Z + t * dz) - (zv + s)
            dF = dz - c * (xx * dx + yy * dy) / qz
            step = -F / dF
            t += step
            if abs(step) < 1e-14:
                break
        X, Y, Z = X + dx * t, Y + dy * t, Z + dz * t
        nb, na = med[i - 1][1], med[i][1]
        if nb == na:
            continue
        r2 = X * X + Y * Y
        qz = math.sqrt(1.0 - c * c * r2)
        nx, ny, nz = -c * X, -c * Y, qz
        cosi = -(nx * dx + ny * dy + nz * dz)
        if cosi < 0.0:
            nx, ny, nz = -nx, -ny, -nz
            cosi = -cosi
        k = 1.0 - (nb / na) ** 2 * (1.0 - cosi * cosi)
        g = (nb / na) * cosi - math.sqrt(k)
        ux = (nb / na) * dx + g * nx
        uy = (nb / na) * dy + g * ny
        uz = (nb / na) * dz + g * nz
        nn = math.sqrt(ux * ux + uy * uy + uz * uz)
        dx, dy, dz = ux / nn, uy / nn, uz / nn
    return Y + dy * (image_dist / dz)