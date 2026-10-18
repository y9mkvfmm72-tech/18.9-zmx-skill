# -*- coding: utf-8 -*-
"""Design gate: the full check a tracer/paraxial pair must pass before any
design number is reported.

Run from this directory:  python check_design_gates.py

A) on the Zemax stock doublet (3- doublet.zmx), the tracer and a naive,
   deliberately independent implementation must agree to < 1e-6 mm at six
   pupil heights.  (The naive reference shares ONLY the glass catalogue with
   the tracer - data, not optics.)
B) small apertures must be linear in the pupil height (no hidden sag leak).
C) a known lens - plano-convex N-BK7, R = 20 mm, so f = R/(n-1) = 38.59 mm -
   must come out right, with the axial ray landing on the focus.
D) EXACT-FOCUS: set the image distance to the paraxial l' from the y-u march;
   a small-aperture ray must then land on the axis.  This is the decisive
   first-order check: it catches a transfer-sign bug that a bare EFL printout
   can hide (2026-10-18: a wrong sign made a 100 mm design look like 183 mm).
E) KNOWN ANSWER, THIRD PATH: the 2x2 matrix method (matrix_firstorder.py) must
   reproduce the march's EFL/l' for the stock doublet to < 1e-6, and both must
   equal 100.0439 / 97.9387.  Gates A/B only compare implementations with each
   other - if both are fed the same wrong data they happily agree (this file
   itself shipped with the wrong sign of the doublet's last surface and passed).
F) CONCAVE SURFACE: a meniscus with a strongly concave rear surface must trace
   and focus correctly.  A tracer that picks the non-negative sphere root
   returns "miss" or a wildly wrong focus here (vertex-plane marching meets a
   concave surface at slightly NEGATIVE t).
"""
import math
import sys

import matrix_firstorder as MX
import naive_ref as R
import paraxial as P
import tracer as T

# the ACTUAL curvatures of the stock doublet (re-read from the file 2026-10-18;
# an earlier revision of this gate had the last one with the wrong sign and
# gate A passed on a nearly-afocal system without noticing)
CURV = [0.0, 0.016249593760156, -0.022527596305474, -0.007676953784738]
THICK = [5.0, 3.0, 1.5, R.BFD]
GLASS = ["", "N-BK7", "N-SF5", ""]
H = 5.0
EFL_REF, LP_REF = 100.0439, 97.9387


def main():
    ok = True
    lens = T.Lens(CURV, THICK, GLASS, h=H, stop=0, fields=[0.0, 5.0],
                  lams=[0.633])

    print("A) tracer vs naive independent reference (on axis)")
    worst = 0.0
    for py in (1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0):
        tr = lens.trace_one(0.0, H * py, 0.0, 0.0, 1.0, 0)[1]
        nv = R.trace(H * py)
        worst = max(worst, abs(tr - nv))
        print("   pupil %-6g  tracer %+.6f  naive %+.6f" % (py, tr, nv))
    print("   worst %.2e -> %s" % (worst, "OK" if worst <= 1e-6 else "FAIL"))
    ok = ok and worst <= 1e-6

    print("B) small-aperture linearity")
    rat = [lens.trace_one(0.0, H * p, 0.0, 0.0, 1.0, 0)[1] / (H * p)
           for p in (1e-5, 1e-4, 1e-3)]
    print("   y/pupil", ["%.8f" % v for v in rat])
    lin = max(rat) - min(rat) < 1e-6
    print("   ->", "OK" if lin else "FAIL")
    ok = ok and lin

    print("C) known lens: plano-convex N-BK7 R=20, f = R/(n-1) = 38.59 mm")
    base = T.Lens([0.0, 0.0, -0.05], [0.0, 3.0, 0.0], ["", "N-BK7", ""],
                  h=10.0, stop=0, fields=[0.0, 5.0], lams=[0.555])
    pb = base.paraxial()
    print("   paraxial efl %.3f (expected 38.590)" % pb.efl(0.555))
    fok = base.with_params(thick=[0.0, 3.0, pb.bfd(0.555)])
    y_ax = fok.trace_one(0.0, 1e-3, 0.0, 0.0, 1.0, 0)[1]
    a = math.radians(5.0)
    y_ch = fok.trace_one(0.0, 0.0, 0.0, math.sin(a), math.cos(a), 0)[1]
    print("   axial ray at image %.6f mm (expect ~0), chief %.3f (expect ~3.376)"
          % (y_ax, y_ch))
    c_ok = (abs(pb.efl(0.555) - 38.59) < 0.02 and abs(y_ax) < 1e-4
            and abs(abs(y_ch) - 3.376) < 0.05)
    print("   ->", "OK" if c_ok else "FAIL")
    ok = ok and c_ok

    print("D) exact-focus: image plane at the paraxial l' from the march")
    lp = lens.paraxial().bfd(0.633)
    foc = lens.with_params(thick=[THICK[0], THICK[1], THICK[2], lp])
    y_small = foc.trace_one(0.0, 0.05, 0.0, 0.0, 1.0, 0)[1]
    y_tiny = foc.trace_one(0.0, 1e-3, 0.0, 0.0, 1.0, 0)[1]
    print("   l' = %.4f ;  h=0.05 lands %+.8f ;  h=1e-3 lands %+.8f"
          % (lp, y_small, y_tiny))
    d_ok = abs(y_small) < 1e-4 and abs(y_tiny) < 1e-6
    print("   ->", "OK" if d_ok else "FAIL")
    ok = ok and d_ok

    print("E) known answer, third path: matrix method on the stock doublet")
    efl_m, lp_m = MX.efl_lprime(CURV[1:], [3.0, 1.5], [R.N1, R.N2, 1.0])
    efl_p, lp_p = lens.paraxial().efl(0.633), lens.paraxial().bfd(0.633)
    print("   matrix : EFL %.6f  l' %.6f" % (efl_m, lp_m))
    print("   march  : EFL %.6f  l' %.6f" % (efl_p, lp_p))
    print("   known  : EFL %.6f  l' %.6f" % (EFL_REF, LP_REF))
    e_ok = (abs(efl_m - efl_p) < 1e-6 and abs(lp_m - lp_p) < 1e-6
            and abs(efl_m - EFL_REF) < 0.001 and abs(lp_m - LP_REF) < 0.001)
    print("   ->", "OK" if e_ok else "FAIL")
    ok = ok and e_ok

    print("F) concave surface: strong meniscus must trace and focus")
    # N-BK7 meniscus, c1 = +1/30, c2 = -1/12, t = 5 mm.  The rear surface is
    # strongly concave: vertex-plane marching meets it at NEGATIVE t (~3.7 mm
    # at h = 3).  A tracer insisting on t >= 0 reports a miss here.
    men = T.Lens([0.0, 1.0 / 30.0, -1.0 / 12.0], [0.0, 5.0, 0.0],
                 ["", "N-BK7", ""], h=10.0, stop=0, fields=[0.0], lams=[0.633])
    pm = men.paraxial()
    efl_men = pm.efl(0.633)
    men_foc = men.with_params(thick=[0.0, 5.0, pm.bfd(0.633)])
    x_m, y_m, bad_m = men_foc.trace_one(0.0, 3.0, 0.0, 0.0, 1.0, 0)
    x_m0, y_m0, bad_m0 = men_foc.trace_one(0.0, 1e-3, 0.0, 0.0, 1.0, 0)
    print("   EFL %.3f (expect ~17.34) ; h=3 lands %+.6f ; h=1e-3 %+.8f"
          % (efl_men, y_m, y_m0))
    # the h=3 ray of an f/2.9 meniscus shows real spherical aberration; the
    # gate only requires the ray to TRACE (no miss/TIR) and the small ray to
    # focus on the axis.  A "t >= 0 only" tracer would return miss here.
    f_ok = (abs(efl_men - 17.34) < 0.05 and not bad_m and not bad_m0
            and abs(y_m0) < 1e-6 and abs(y_m) < 0.5)
    print("   ->", "OK" if f_ok else "FAIL")
    ok = ok and f_ok

    print("")
    print("VERDICT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

