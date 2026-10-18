# -*- coding: utf-8 -*-
"""Independent first-order checker: the 2x2 ray-transfer matrix.

State vector is (y, nu) with nu = n * theta (the reduced angle):
    refraction : (y, nu) -> (y, nu - phi*y),       phi = (n' - n) * c
    transfer   : (y, nu) -> (y + d*nu/n, nu)       (n = medium of travel)

For an object at infinity, y0 = 1 at surface 0, nu0 = 0, and after the whole
system (A, C) = M . (1, 0):
    EFL = -n_air / C          (focal length)
    l'  = -A / (C / n_air)    (back focal distance from the last vertex)

This is a THIRD independent way of getting the first-order quantities
(the other two are the y-u march in paraxial.py and the real ray trace).
It exists because 2026-10-18 showed that a sign slip in ONE engine can stay
hidden for weeks when it is the only thing the numbers are compared against.

Run directly to verify the stock doublet:
    python matrix_firstorder.py
"""
import math


def _mul(A, B):
    return [[A[0][0] * B[0][0] + A[0][1] * B[1][0],
             A[0][0] * B[0][1] + A[0][1] * B[1][1]],
            [A[1][0] * B[0][0] + A[1][1] * B[1][0],
             A[1][0] * B[0][1] + A[1][1] * B[1][1]]]


def system_matrix(curv, thick, n_after):
    """curv/thick/n_after are per SURFACE (n_after[i] = medium after surface i;
    1.0 for air).  thick[i] is the vertex-to-vertex distance i -> i+1;
    thick may be one element shorter than curv."""
    M = [[1.0, 0.0], [0.0, 1.0]]
    n = 1.0
    for i, c in enumerate(curv):
        n2 = n_after[i] if i < len(n_after) else 1.0
        phi = (n2 - n) * c
        M = _mul([[1.0, 0.0], [-phi, 1.0]], M)
        if i < len(thick):
            M = _mul([[1.0, thick[i] / n2], [0.0, 1.0]], M)
        n = n2
    return M, n


def efl_lprime(curv, thick, n_after):
    """-> (EFL, l' from the last vertex).  Image medium must be air."""
    M, n_last = system_matrix(curv, thick, n_after)
    A, C = M[0][0], M[1][0]
    if abs(n_last - 1.0) > 1e-9:
        raise ValueError("matrix_firstorder: image medium must be air")
    return -1.0 / C, -A / C


def main():
    # the stock 3- doublet.zmx (stop plane + 3 refracting surfaces)
    import glass_index as G
    cat = G.load_agf()[0]
    n1 = G.n_index("N-BK7", 0.633, cat)
    n2 = G.n_index("N-SF5", 0.633, cat)
    curv = [0.0, 0.016249593760156, -0.022527596305474, -0.007676953784738]
    thick = [5.0, 3.0, 1.5]
    na = [1.0, n1, n2, 1.0]
    efl, lp = efl_lprime(curv, thick, na)
    print("stock doublet (matrix method): EFL = %.6f   l' = %.6f" % (efl, lp))
    print("expected                       EFL = 100.043872   l' = 97.938653")
    ok = abs(efl - 100.0439) < 0.001 and abs(lp - 97.9387) < 0.001
    print("VERDICT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
