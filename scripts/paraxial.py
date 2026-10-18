# -*- coding: utf-8 -*-
"""Paraxial (first-order) quantities: independent cross-check for the tracer.

State at surface i is (y, n*u).
    refraction : n'u' = n*u - y * c * (n' - n)
    transfer   : y' = y + d * (n'u') / n'
"""
import math

import glass_index as G


class Paraxial(object):
    def __init__(self, curv, thick, glass, stop=0, cat=None):
        self.curv = [float(v) for v in curv]
        self.thick = [float(v) for v in thick]
        self.glass = list(glass)
        self.stop = int(stop)
        self.curv0 = list(curv)
        self.glass0 = list(glass)
        self.last = len(self.curv) - 1
        self.cat = G.load_agf()[0] if cat is None else cat
        self.air_before = True

    def n_before(self, i, lam):
        if i == 0:
            return 1.0
        g = self.glass[i - 1]
        return G.n_index(g, lam, self.cat) if g else 1.0

    def n_after(self, i, lam):
        g = self.glass[i] if i < len(self.glass) else ""
        return G.n_index(g, lam, self.cat) if g else 1.0

    def power(self, i, lam):
        return self.curv[i] * (self.n_after(i, lam) - self.n_before(i, lam))

    def total_power(self, lam):
        return sum(self.power(i, lam) for i in range(self.last + 1))

    def march(self, y, u, i0, lam):
        n = self.n_before(i0, lam)
        nu = n * u
        for i in range(i0, self.last + 1):
            nu = nu - y * self.power(i, lam)
            n = self.n_after(i, lam)
            if i < self.last:
                y = y + self.thick[i] * nu / n
        return y, nu / n

    def efl(self, lam):
        # EFL is NOT 1/total_power for separated elements: the plain sum of
        # surface powers is only valid for thin elements in contact.  Use the
        # y-nu march: for a unit-height ray from infinity, u' (in air) = -1/EFL.
        # (measured 2026-10-18: total_power gives 99.54 for the stock doublet,
        # the march and the matrix method give 100.04, and the real-trace focus
        # lands at the march's bfd - so the sum was wrong by 0.5 mm.)
        y, u = self.march(1.0, 0.0, 0, lam)
        return -1.0 / u

    def bfd(self, lam):
        y, u = self.march(1.0, 0.0, 0, lam)
        return -y / u

    def image_height(self, field_deg, lam, h=0.0, pupil_y=0.0, image_z=None):
        if image_z is None:
            image_z = sum(self.thick)
        tanf = math.tan(math.radians(field_deg))
        z_stop = sum(self.thick[:self.stop])
        y0 = h * pupil_y - z_stop * tanf
        y, u = self.march(y0, tanf, 0, lam)
        z_last = sum(self.thick[:self.last])
        return y + u * (image_z - z_last)