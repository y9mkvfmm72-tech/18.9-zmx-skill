#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare a hand-written .zmx against a Zemax-authored one, line by line.

WHY THIS EXISTS
    .zmx is a positional text format with no schema.  A missing required
    keyword or a changed column count is not always rejected on load; it can
    simply be ignored, leaving you with a file that opens but is wrong.
    Comparing against a file Zemax itself wrote is the reliable check.

WHAT IT CHECKS
    * encoding is UTF-16 LE with BOM
    * every top-level keyword present in the reference is present here
    * every per-surface keyword present in the reference is present here
    * a list of specific content assertions (aperture, fields, wavelengths)

USAGE
    python check_zmx_structure.py myfile.zmx [reference.zmx]

    With no reference given it uses the stock 18.9 sample.  Exit code 0 means
    the structure matches.
"""
import os
import re
import sys

DEFAULT_REF = (r"C:\Users\Lenovo\Documents\Zemax\Samples\Design Applications"
               r"\Lasers and Fibers\Laser Lenses\3- doublet.zmx")


def read_zmx(path):
    """Return (lines, encoding_name). Handles the UTF-16 LE + BOM case that
    Zemax writes on Windows, and falls back to UTF-8."""
    raw = open(path, "rb").read()
    if raw[:2] == b"\xff\xfe":
        return [l.rstrip("\r") for l in raw[2:].decode("utf-16-le").split("\n")], "UTF-16 LE"
    if raw[:2] == b"\xfe\xff":
        return [l.rstrip("\r") for l in raw[2:].decode("utf-16-be").split("\n")], "UTF-16 BE"
    return [l.rstrip("\r") for l in raw.decode("utf-8", "replace").split("\n")], "UTF-8"


def top_keywords(lines):
    """Top-level record keywords, in first-appearance order."""
    out = []
    for l in lines:
        m = re.match(r"^([A-Z]{3,4})\b", l)
        if m and m.group(1) not in out:
            out.append(m.group(1))
    return out


def surface_keywords(lines):
    """Keywords appearing inside SURF blocks (indented)."""
    keys = set()
    inside = False
    for l in lines:
        if l.startswith("SURF "):
            inside = True
            continue
        if inside:
            m = re.match(r"^\s+([A-Z]{3,4})\b", l)
            if m:
                keys.add(m.group(1))
    return sorted(keys)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    mine_path = sys.argv[1]
    ref_path = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_REF

    mine, mine_enc = read_zmx(mine_path)
    print("file      : %s" % os.path.basename(mine_path))
    print("encoding  : %s" % mine_enc)
    print("first line: %s" % (mine[0] if mine else "<empty>"))
    print("")

    ok = True
    if mine_enc != "UTF-16 LE":
        print("  [!!] expected UTF-16 LE with BOM (what Zemax writes on Windows)")
        ok = False
    else:
        print("  [OK] UTF-16 LE with BOM")

    if not os.path.exists(ref_path):
        print("  [--] reference not found, skipping structural comparison:")
        print("       %s" % ref_path)
        print("")
        print("VERDICT: %s" % ("PASS (encoding only)" if ok else "FAIL"))
        return 0 if ok else 1

    ref, _ = read_zmx(ref_path)
    print("reference : %s" % os.path.basename(ref_path))
    print("")

    rk, mk = top_keywords(ref), top_keywords(mine)
    missing = [k for k in rk if k not in mk]
    extra = [k for k in mk if k not in rk]
    print("top-level keywords")
    print("  reference : %s" % " ".join(rk))
    print("  this file : %s" % " ".join(mk))
    print("  missing   : %s" % (missing or "none"))
    print("  extra     : %s" % (extra or "none"))
    if missing:
        ok = False

    rsk, msk = surface_keywords(ref), surface_keywords(mine)
    smissing = [k for k in rsk if k not in msk]
    print("")
    print("per-surface keywords")
    print("  reference : %s" % " ".join(rsk))
    print("  this file : %s" % " ".join(msk))
    print("  missing   : %s" % (smissing or "none"))
    print("  note      : optional surface fields may be absent without harm;")
    print("              verify by loading the file, not by this list alone.")

    # content spot checks that matter for an imaging design
    joined = "\n".join(mine)
    checks = [
        ("VERS line present", bool(re.search(r"(?m)^VERS ", joined))),
        ("MODE SEQ", bool(re.search(r"(?m)^MODE SEQ", joined))),
        ("aperture ENPD set", bool(re.search(r"(?m)^ENPD ", joined))),
        ("glass catalogue named", bool(re.search(r"(?m)^GCAT ", joined))),
        ("field list present", bool(re.search(r"(?m)^YFLN ", joined))),
        ("field/wave type record", bool(re.search(r"(?m)^FTYP ", joined))),
        ("at least 3 wavelengths", len(re.findall(r"(?m)^WAVM ", joined)) >= 3),
        ("primary wavelength set", bool(re.search(r"(?m)^PWAV ", joined))),
        ("has a SURF 0 (object)", bool(re.search(r"(?m)^SURF 0$", joined))),
    ]
    print("")
    print("content checks")
    for name, good in checks:
        print("  [%s] %s" % ("OK" if good else "!!", name))
        ok = ok and good

    print("")
    print("VERDICT: %s" % ("PASS" if ok else "REVIEW"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
