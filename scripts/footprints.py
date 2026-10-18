# -*- coding: utf-8 -*-
"""Vignetting detector: read a .zmx, trace the design beam, and compare each
surface's BEAM semi-height with the DIAM written in the file.

WHY THIS EXISTS
    On 2026-10-18 a user file claimed f/4.5 but its written clear apertures
    were half the actual beam heights (stop semi 5.27 mm vs the 11.11 mm the
    ENPD needed).  Nothing in Zemax or in the file errors on that -- the lens
    just quietly vignettes to about f/12.8.  The only way to catch it is to
    compute the beam footprint per surface and compare with DIAM.

USAGE
    python footprints.py my_lens.zmx

Reports, per surface: written DIAM/2, beam semi-height (max over all traced
rays), and the ratio.  Exit code 1 if any written semi-diameter is smaller
than the beam.
"""
import math
import sys

import tracer as T

RINGS = (0.25, 0.5, 0.75, 1.0)
ARMS = 8
FIELDS_MAX = 3          # up to this many Y fields are swept


def parse_zmx(path):
    raw = open(path, "rb").read()
    txt = (raw[2:].decode("utf-16-le") if raw[:2] == b"\xff\xfe"
           else raw.decode("utf-8", "replace"))
    surfs, sidx = {}, None
    head = {"enpd": None, "yfln": []}
    for ln in txt.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "SURF":
            sidx = int(t[1])
            surfs[sidx] = {"stop": False, "glass": "", "c": 0.0, "t": 0.0,
                           "diam": 0.0}
        elif t[0] == "ENPD":
            head["enpd"] = float(t[1])
        elif t[0] == "YFLN":
            head["yfln"] = [float(v) for v in t[1:] if float(v) != 0.0]
        elif sidx is not None:
            if t[0] == "STOP":
                surfs[sidx]["stop"] = True
            elif t[0] == "CURV":
                surfs[sidx]["c"] = float(t[1])
            elif t[0] == "DISZ":
                surfs[sidx]["t"] = (0.0 if t[1] == "INFINITY" else float(t[1]))
            elif t[0] == "GLAS":
                surfs[sidx]["glass"] = t[1]
            elif t[0] == "DIAM":
                surfs[sidx]["diam"] = float(t[1])
    return surfs, head


def build_lens(surfs, head):
    """Convert the parsed file into a tracer.Lens in its native convention."""
    order = sorted(surfs)
    last = order[-1]
    curv, thick, glass = [], [], []
    stop = 0
    for k in order[:-1]:                     # skip the image surface
        s = surfs[k]
        curv.append(s["c"])
        thick.append(s["t"])
        glass.append(s["glass"])
        if s["stop"]:
            stop = len(curv) - 1
    return T.Lens(curv, thick, glass, h=head["enpd"] / 2.0, stop=stop,
                  fields=head["yfln"][:FIELDS_MAX] or [0.0], lams=[T.D])


def footprints(lens):
    """Max |ray height| per refracting surface over the pupil/field sweep."""
    verts = lens.vertices()
    semi = [0.0] * (lens.last + 1)
    pts = []
    for r in RINGS:
        for a in range(ARMS):
            th = 2.0 * math.pi * a / ARMS
            pts.append((r * math.cos(th), r * math.sin(th)))
    pts.append((0.0, 0.0))
    for field in lens.fields:
        v0 = lens.aim(field)
        a = math.sin(math.radians(field))
        b = math.cos(math.radians(field))
        for (px, py) in pts:
            rec = []
            lens.trace_one(px * lens.h, v0 + py * lens.h, 0.0, a, b, 0,
                           record=rec)
            for (i, X, Y, Z, dx, dy, dz) in rec:
                h = math.hypot(X, Y)
                if h > semi[i]:
                    semi[i] = h
    return semi


def main(argv):
    path = argv[1]
    surfs, head = parse_zmx(path)
    lens = build_lens(surfs, head)
    semi = footprints(lens)
    order = sorted(surfs)
    # NOTE: the .zmx DIAM field is the clear SEMI-DIAMETER (radius), not the
    # diameter -- the stock doublet proves it: ENPD 10 -> stop DIAM 5.
    # An earlier revision of this script divided DIAM by 2 again and reported
    # vignetting everywhere.
    print("%-4s %12s %12s %8s  %s" % ("surf", "written R", "beam semi",
                                      "ratio", "note"))
    ok = True
    for i, k in enumerate(order[:-1]):
        w = surfs[k]["diam"]
        b = semi[i]
        ratio = (w / b) if b > 1e-9 else float("inf")
        note = ""
        if surfs[k]["stop"]:
            note += "STOP "
        if w + 1e-9 < b and w > 0:
            note += "VIGNETTING (beam larger than written clear semi-diameter)"
            ok = False
        print("%-4d %12.3f %12.3f %8.2f  %s" % (k, w, b, ratio, note))
    print("")
    print("VERDICT:", "PASS (no surface clips the beam)" if ok else
          "FAIL (some written DIAM is smaller than the beam)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
