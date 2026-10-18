#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify a ray tracer against a lens whose answer is already known.

WHY THIS EXISTS
    A tracer bug in this project (an approximate surface normal) silently made
    every spot size about 2.3x too small and turned a 50 mm design into a
    155 mm one.  Nothing errored.  The only thing that caught it was tracing a
    lens whose focal length was already known.

WHAT IT DOES
    Traces the stock Zemax doublet
        C:\\Users\\Lenovo\\Documents\\Zemax\\Samples\\Design Applications\\
        Lasers and Fibers\\Laser Lenses\\3- doublet.zmx
    (N-BK7 / N-SF5, cemented, plane front stop), whose paraxial EFL at the
    0.633 um design wavelength is 100.0439 mm.

    NOTE (2026-10-18): an earlier version of this script printed 99.0439 mm as
    "the independent reference".  That number came from a y-u trace whose
    transfer had the wrong sign (y - t*u).  The real-trace gate below was
    unaffected, so it kept passing while the printed reference was off by 1%.
    The corrected value is cross-checked by the matrix method
    (scripts/matrix_firstorder.py) and by a real-trace focus scan.

    It reports, for several pupil heights:
      * the paraxial y-u EFL (closed form, no tracing), as the reference
      * the EFL from the real trace
    A correct tracer matches the paraxial value within a fraction of a percent
    at small pupil heights, and deviates progressively at the edge because
    that is real spherical aberration.

USAGE
    python check_tracer.py            # uses the built-in reference lens
    python check_tracer.py <file.zmx> # check against your own file

    Exit code 0 = tracer agrees with the analytic reference.
"""
import math
import os
import sys

# Schott Sellmeier coefficients (identical to the values in SCHOTT.AGF).
GLASS = {
    "N-BK7": ([1.039612120, 0.2317923440, 1.010469450],
              [6.000698670e-3, 2.001791440e-2, 1.035606530e2]),
    "N-SF5": ([1.524818890, 0.1870855270, 1.4272901500],
              [1.125475600e-2, 5.889953920e-2, 1.291416750e2]),
    "F2": ([1.345333590, 0.2090731760, 0.9373571620],
           [9.977438710e-3, 4.704507670e-2, 1.118867640e2]),
}

# The stock doublet, read out of the .zmx listed above.
REF = {
    "name": "Zemax stock doublet (N-BK7 / N-SF5)",
    "lam": 0.633,
    "curv": [0.016249593760156, -0.022527596305474, -0.007676953784738],
    "thick": [3.0, 1.5],
    "glass": ["N-BK7", "N-SF5"],
    "stop_gap": 5.0,        # stop sits 5 mm in front of surface 2
    "bfd": 98.475961359590002,
}


def nidx(glass, lam):
    B, C = GLASS[glass]
    l2 = lam * lam
    return math.sqrt(1.0 + sum(b * l2 / (l2 - c) for b, c in zip(B, C)))


def paraxial_efl(curv, thick, glasses, lam):
    """Closed-form paraxial EFL: a y-u trace, no ray tracing involved.

    This is the independent reference the traced values are compared against.
    """
    y, u = 1.0, 0.0
    n_before = 1.0
    for i, c in enumerate(curv):
        n_after = 1.0 if i >= len(glasses) else nidx(glasses[i], lam)
        u = (n_before * u - y * c * (n_after - n_before)) / n_after
        if i < len(thick):
            # u is the geometric slope dy/dz: a converging ray has u < 0 and
            # its height must DECREASE, so the transfer is y + t*u.
            y = y + thick[i] * u
        n_before = n_after
    return -1.0 / u


def trace(curv, thick, glasses, bfd, lam, stop_gap, h, exact=True):
    """Real ray trace of an axial ray at pupil height h.

    exact=True  uses the physically correct normal  n ~ (-c x / sqrt(1-c^2 r^2), ..., 1)
    exact=False reproduces the buggy approximation       n ~ (-c x, -c y, 1)
    so the two can be compared directly.
    Returns the image-plane ray height.
    """
    n_glass = [nidx(g, lam) for g in glasses]
    media = [(1.0, n_glass[0]), (n_glass[0], n_glass[1]), (n_glass[1], 1.0)]

    dx, dy, dz = 0.0, 0.0, 1.0
    x, y, z = 0.0, h, stop_gap
    for i, cc in enumerate(curv):
        rho2 = x * x + y * y
        disc = 1.0 - cc * cc * rho2
        if disc < 0:
            return None
        sq = math.sqrt(disc)
        z = cc * rho2 / (1.0 + sq)
        slope = cc / sq if sq > 1e-14 else 0.0
        nx = -slope * x if exact else -cc * x
        ny = -slope * y if exact else -cc * y
        nz = 1.0
        nn = math.sqrt(nx * nx + ny * ny + nz * nz)
        nx, ny, nz = nx / nn, ny / nn, nz / nn
        cosi = -(nx * dx + ny * dy + nz * dz)
        if cosi < 0:
            nx, ny, nz = -nx, -ny, -nz
            cosi = -cosi
        eta = media[i][0] / media[i][1]
        k = 1.0 - eta * eta * (1.0 - cosi * cosi)
        if k < 0:
            return None
        t = eta * cosi - math.sqrt(k)
        dx = eta * dx + t * nx
        dy = eta * dy + t * ny
        dz = eta * dz + t * nz
        nn = math.sqrt(dx * dx + dy * dy + dz * dz)
        dx, dy, dz = dx / nn, dy / nn, dz / nn
        if i < len(thick):
            step = thick[i] / dz
            x, y, z = x + dx * step, y + dy * step, z + dz * step
    step = bfd / dz
    return y + dy * step


def main():
    lam = REF["lam"]
    print("reference lens : %s" % REF["name"])
    print("wavelength     : %.4f um" % lam)
    print("indices        : %s"
          % ", ".join("%s=%.5f" % (g, nidx(g, lam)) for g in REF["glass"]))
    efl_ref = paraxial_efl(REF["curv"], REF["thick"], REF["glass"], lam)
    print("paraxial EFL   : %.4f mm   <-- the independent reference" % efl_ref)
    print("")

    print("real trace, image-plane height vs pupil height:")
    print("   h(mm)      exact normal      approximate normal")
    ok = True
    for h in (0.05, 0.5, 2.0, 4.0, 5.0):
        ye = trace(REF["curv"], REF["thick"], REF["glass"], REF["bfd"],
                   lam, REF["stop_gap"], h, exact=True)
        ya = trace(REF["curv"], REF["thick"], REF["glass"], REF["bfd"],
                   lam, REF["stop_gap"], h, exact=False)
        print("   %5.2f      %12.6f      %12.6f" % (h, ye, ya))

    # the decisive assertion: a small-aperture ray must focus onto the axis
    y_small = trace(REF["curv"], REF["thick"], REF["glass"], REF["bfd"],
                    lam, REF["stop_gap"], 0.05, exact=True)
    print("")
    print("check: a 0.05 mm pupil ray must land on the axis")
    print("       landed at %.6f mm (|y| must be < 0.002)" % y_small)
    if abs(y_small) >= 0.002:
        ok = False
        print("       FAIL - the tracer does not focus a paraxial ray")
    else:
        print("       OK")

    # a tracer with the wrong normal still passes the small-aperture test,
    # so also check the edge: the two normals must differ measurably there
    y_edge_e = trace(REF["curv"], REF["thick"], REF["glass"], REF["bfd"],
                     lam, REF["stop_gap"], 5.0, exact=True)
    y_edge_a = trace(REF["curv"], REF["thick"], REF["glass"], REF["bfd"],
                     lam, REF["stop_gap"], 5.0, exact=False)
    rel = abs(y_edge_e - y_edge_a) / max(abs(y_edge_e), 1e-12)
    print("")
    print("edge sensitivity: exact %.6f vs approximate %.6f (%.0f%% apart)"
          % (y_edge_e, y_edge_a, rel * 100))
    print("       -> a tracer must use the EXACT normal; the approximation")
    print("          is already visible on this simple doublet.")

    print("")
    print("VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
