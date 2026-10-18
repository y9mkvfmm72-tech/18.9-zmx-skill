# -*- coding: utf-8 -*-
"""Sequential real-ray tracer for spherical surfaces.

Conventions
-----------
* Surface i is a sphere through (0, 0, z_i) on the axis.  Local sag:
      sag(r) = c*r^2 / (1 + sqrt(1 - c^2 r^2))        (bump for c > 0)
  A ray along +z sees the surface at z = z_i + sag(r).
* `thick[i]`   : vertex-to-vertex distance i -> i+1.  `thick[last]` = image gap.
* `glass[i]`   : medium AFTER surface i ('' = air); `glass[last]` unused.
* `stop`       : index of the aperture stop.  Pupil coords (xp, yp) are
  normalised to [-1, 1] over radius `h` measured in the stop plane.

The tracer is deliberately plain scalar Python so that every intermediate
quantity can be inspected; the bundle driver just loops over rays.
"""
import math

import numpy as np

import glass_index as G
import paraxial as PX

RINGS = (0.35, 0.6, 0.8, 0.925, 1.0)
ARMS = 20
F, D, C = 0.4861327, 0.5875618, 0.6562725
LAMS = (F, D, C)

_CAT = None


def catalogue():
    global _CAT
    if _CAT is None:
        _CAT = G.load_agf()[0]
    return _CAT


def pupil_points():
    pts = []
    for k, r in enumerate(RINGS):
        off = (math.pi / ARMS) if k % 2 else 0.0
        for j in range(ARMS):
            th = 2.0 * math.pi * j / ARMS + off
            pts.append((r * math.cos(th), r * math.sin(th)))
    pts.append((0.0, 0.0))
    return np.array(pts)


class Lens(object):
    def __init__(self, curv, thick, glass, h, stop, fields=(0.0, 7.0, 10.0),
                 lams=LAMS, cat=None):
        self.curv = [float(v) for v in curv]
        self.thick = [float(v) for v in thick]
        self.glass = list(glass)
        self.h = float(h)
        self.stop = int(stop)
        self.fields = [float(f) for f in fields]
        self.lams = [float(l) for l in lams]
        self.cat = catalogue() if cat is None else cat
        self.last = len(self.curv) - 1

    # ------------------------------------------------------------- plumbing
    def copy(self):
        return Lens(self.curv, self.thick, self.glass, self.h, self.stop,
                    self.fields, self.lams, self.cat)

    def with_params(self, curv=None, thick=None, glass=None):
        o = self.copy()
        if curv is not None:
            o.curv = [float(v) for v in curv]
        if thick is not None:
            o.thick = [float(v) for v in thick]
        if glass is not None:
            o.glass = list(glass)
        o.last = len(o.curv) - 1
        return o

    def n_of(self, name, lam):
        return G.n_index(name, lam, self.cat) if name else 1.0

    def media(self, wl):
        lam = self.lams[wl]
        return [(1.0 if i == 0 else self.n_of(self.glass[i - 1], lam),
                 1.0 if i >= self.last else self.n_of(self.glass[i], lam))
                for i in range(self.last + 1)]

    def paraxial(self):
        return PX.Paraxial(self.curv, self.thick, self.glass, self.stop, self.cat)

    def vertices(self):
        v = [0.0]
        for t in self.thick[:self.last]:
            v.append(v[-1] + t)
        return v

    # ------------------------------------------------------------- tracing
    @staticmethod
    def _sag(c, x, y):
        r2 = x * x + y * y
        disc = 1.0 - c * c * r2
        if disc <= 0.0:
            raise ArithmeticError
        qz = math.sqrt(disc)
        return c * r2 / (1.0 + qz), qz

    def trace_one(self, x, y, dx, dy, dz, wl, record=None):
        """Trace one ray from (x, y) on surface 0's plane z = 0."""
        verts = self.vertices()
        med = self.media(wl)
        X, Y, Z = x, y, 0.0
        bad = False
        for i in range(self.last + 1):
            c = self.curv[i]
            zv = verts[i]
            t = (zv - Z) / dz
            if c != 0.0:
                for _ in range(100):
                    xx = X + dx * t
                    yy = Y + dy * t
                    try:
                        s, qz = self._sag(c, xx, yy)
                    except ArithmeticError:
                        bad = True
                        break
                    F = (Z + t * dz) - (zv + s)
                    dF = dz - c * (xx * dx + yy * dy) / qz
                    step = -F / dF
                    t += step
                    if abs(step) < 1e-13:
                        break
            X += dx * t
            Y += dy * t
            Z += dz * t
            if record is not None:
                record.append((i, X, Y, Z, dx, dy, dz))
            nb, na = med[i]
            if nb == na:
                # no index change here (air->air): nothing to refract
                if i == self.last:
                    break
                continue
            r2 = X * X + Y * Y
            qz = math.sqrt(max(1.0 - c * c * r2, 0.0))
            nx, ny, nz = -c * X, -c * Y, qz     # points upstream (+z)
            cosi = -(nx * dx + ny * dy + nz * dz)
            if cosi < 0.0:                      # ray overtook the normal: flip
                nx, ny, nz = -nx, -ny, -nz
                cosi = -cosi
            k = 1.0 - (nb / na) ** 2 * (1.0 - cosi * cosi)
            if k < 0.0:
                bad = True
                break
            g = (nb / na) * cosi - math.sqrt(k)
            ux = (nb / na) * dx + g * nx
            uy = (nb / na) * dy + g * ny
            uz = (nb / na) * dz + g * nz
            nn = math.sqrt(ux * ux + uy * uy + uz * uz)
            dx, dy, dz = ux / nn, uy / nn, uz / nn
            if i == self.last:
                break
        t = self.thick[self.last] / dz
        return X + dx * t, Y + dy * t, bad

    def aim(self, field_deg):
        """Start height on surface 0's plane for the chief ray of a field."""
        if self.stop == 0:
            return 0.0
        a = math.sin(math.radians(field_deg))
        b = math.cos(math.radians(field_deg))
        z_stop = sum(self.thick[:self.stop])

        def y_at_stop(v):
            rec = []
            try:
                self.trace_one(0.0, v, 0.0, a, b, 0, record=rec)
            except ArithmeticError:
                return None
            if len(rec) <= self.stop:
                return None            # ray died before reaching the stop
            return rec[self.stop][2]

        v = -z_stop * math.tan(math.radians(field_deg))
        for _ in range(50):
            f0 = y_at_stop(v)
            f1 = y_at_stop(v + 1e-6)
            if f0 is None or f1 is None:
                break
            der = (f1 - f0) / 1e-6
            if abs(der) < 1e-12:
                break
            step = -f0 / der
            v += step
            if abs(step) < 1e-12:
                break
        return v

    def trace(self, xp, yp, field_deg, wl):
        xp = np.asarray(xp, float)
        yp = np.asarray(yp, float)
        v0 = self.aim(field_deg)
        a = math.sin(math.radians(field_deg))
        b = math.cos(math.radians(field_deg))
        X, Y, BAD = [], [], []
        for k in range(len(xp)):
            x, y, bad = self.trace_one(xp[k] * self.h, v0 + yp[k] * self.h,
                                       0.0, a, b, wl)
            X.append(x)
            Y.append(y)
            BAD.append(bad)
        return dict(x=np.array(X), y=np.array(Y), bad=np.array(BAD))

    def spot(self, field_deg, wl, pts=None):
        if pts is None:
            pts = pupil_points()
        r = self.trace(pts[:, 0], pts[:, 1], field_deg, wl)
        if r["bad"].any():
            return dict(ok=False, nbad=int(r["bad"].sum()))
        x, y = r["x"], r["y"]
        cx, cy = x.mean(), y.mean()
        rms = float(np.sqrt(((x - cx) ** 2 + (y - cy) ** 2).mean()))
        return dict(ok=True, rms=rms, cx=cx, cy=cy, nbad=0)
