#!/usr/bin/env python3
"""Clean-room replay of the cubic-seven K12^2 continuation through H169.

This imports no production search package.  It extends the raw arithmetic and
K12 implementation from the independent H33 replay, reconstructs controls and
Gram types, and reruns every shell correction that affects the 7/5 ray.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from functools import lru_cache
from itertools import combinations, product
from math import gcd, isqrt
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "scripts/independent_hostile_audit_cubic7_coxeter_todd_squared_h33.py"
OWNER = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h169.json"
AUDIT = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h169_independent_hostile_audit.json"
REPORT = ROOT / "results/cubic7_coxeter_todd_squared_h169_independent_hostile_audit.md"


spec = importlib.util.spec_from_file_location("cubic7_h169_cleanroom_base", BASE)
h = importlib.util.module_from_spec(spec); sys.modules[spec.name] = h; spec.loader.exec_module(h)
C7, F, T, A = h.C7, h.F, h.T, h.A
CUT_H, LIMIT = h.CUT_H, C7(169, 474, 211)


def controls():
    out = set()
    for q in (C7(2), C7(6)):
        for x in h.bounded((q * LIMIT / A, q / (-A), q / (-A))):
            if x.sign() <= 0 or not (2 * A * x / q).integral:
                continue
            target = q + A * x * x
            if any(target.sign(e) < 0 for e in (1, 2)):
                continue
            height = A * x * x / q
            if height <= LIMIT:
                out.add(h.Level(height, q, x, target))
    return tuple(sorted(out))


def gram(v):
    return ((v[0], v[3], v[4]), (v[3], v[1], v[5]), (v[4], v[5], v[2]))


def types(level):
    if level.norm == C7(2):
        return tuple(gram(v) for v in h.types(level))
    assert level.norm == C7(6) and (level.target / 3).integral
    fake = h.Level(level.height, C7(2), level.control, level.target / 3)
    return tuple(tuple(tuple(3 * x for x in row) for row in gram(v)) for v in h.types(fake))


def dmax(matrix):
    raw = h.Matrix(matrix); rank = raw.rank(); kernel = raw.nullspace()
    inds = next(s for s in combinations(range(3), rank) if raw.extract(s, s).det())
    principal = raw.extract(inds, inds); det = int(principal.det()); adj = principal.adjugate()
    bounds = tuple(isqrt(4 * matrix[i][i] // 3) + 2 for i in range(3)); best = arg = None; count = 0
    for p in product(*(range(-b, b + 1) for b in bounds)):
        if any(sum(k[i] * p[i] for i in range(3)) for k in kernel):
            continue
        z = tuple(p[i] for i in inds)
        q = sum(int(adj[i, j]) * z[i] * z[j] for i in range(rank) for j in range(rank))
        if 3 * q > 4 * det:
            continue
        count += 1; value = C7(*p)
        if best is None or best < value:
            best, arg = value, p
    return best, arg, count


def gp_shell(matrix, bound):
    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required")
    source = h.gp_matrix(matrix)
    program = f"V=qfminim({source},{bound},,2)[3];for(j=1,matsize(V)[2],print(vector({len(matrix)},i,V[i,j])));quit()\n"
    done = subprocess.run([executable, "-q", "-s", "512000000"], input=program, text=True, capture_output=True, timeout=180, check=True)
    reps = np.array([ast.literal_eval(line) for line in done.stdout.splitlines() if line.strip()], dtype=np.int64)
    return np.concatenate((reps, -reps))


def radical_data():
    raw = h.Matrix(h.K12); D, L, R = h.smith_normal_decomp(raw, domain=h.ZZ); assert L * raw * R == D
    factors = tuple(3 // gcd(3, abs(int(D[i, i]))) for i in range(12))
    basis = h.hermite_normal_form(R * h.Matrix.diag(*factors)); q = basis.T * raw * basis
    return np.array(basis.tolist(), dtype=np.int64), np.array(q.tolist(), dtype=np.int64)


@lru_cache(maxsize=1)
def shells():
    ordinary_q = np.array(h.K12, dtype=np.int64); ordinary = gp_shell(h.K12, 8)
    ordinary_norm = np.einsum("bi,ij,bj->b", ordinary, ordinary_q, ordinary)
    rbasis, radical_q = radical_data(); radical = gp_shell(tuple(tuple(map(int, row)) for row in radical_q), 18)
    radical_norm = np.einsum("bi,ij,bj->b", radical, radical_q, radical)
    return {
        "ordinary": (ordinary_q, {n: ordinary[ordinary_norm == n] for n in (4, 6, 8)}, None),
        "radical": (radical_q, {n: radical[radical_norm == n] for n in (12, 18)}, rbasis),
    }


def bitmasks(matrix, value):
    output = []
    for row in matrix:
        mask = 0
        for index in np.flatnonzero(row == value):
            mask |= 1 << int(index)
        output.append(mask)
    return output


def shell_census(lattice, reduced, transport):
    q, by_norm, basis = shells()[lattice]; r = len(reduced); vectors = [by_norm[reduced[i][i]] for i in range(r)]
    actual = [v if basis is None else v @ basis.T for v in vectors]
    transport = np.array(transport, dtype=np.int64)
    if r == 2:
        pairing = vectors[0] @ q @ vectors[1].T; seen = {}; total = 0
        for i, j in np.argwhere(pairing == reduced[0][1]):
            total += 1; new_h = np.array((actual[0][i, 0], actual[1][j, 0])); old_h = tuple(map(int, new_h @ transport)); seen.setdefault(old_h, (int(i), int(j)))
    else:
        ab = vectors[0] @ q @ vectors[1].T; ac = vectors[0] @ q @ vectors[2].T; bc = vectors[1] @ q @ vectors[2].T
        ma = bitmasks(ac, reduced[0][2]); mb = bitmasks(bc, reduced[1][2]); groups = {}
        for hc in set(map(int, actual[2][:, 0])):
            mask = 0
            for index in np.flatnonzero(actual[2][:, 0] == hc): mask |= 1 << int(index)
            groups[hc] = mask
        seen = {}; total = 0
        for i in range(len(vectors[0])):
            for j in np.flatnonzero(ab[i] == reduced[0][1]):
                overlap = ma[i] & mb[int(j)]; total += overlap.bit_count()
                if not overlap: continue
                for hc, group in groups.items():
                    hit = overlap & group
                    if hit:
                        k = (hit & -hit).bit_length() - 1
                        new_h = np.array((actual[0][i, 0], actual[1][j, 0], hc)); old_h = tuple(map(int, new_h @ transport)); seen.setdefault(old_h, (i, int(j), k))
    best = max(seen, key=lambda z: C7(*z)); indices = seen[best]
    chosen = [actual[i][indices[i]] for i in range(r)]; old = np.stack(chosen, axis=1) @ transport
    return C7(*best), best, total, tuple(tuple(map(int, col)) for col in old.T)


CASES = {
    (C7(81, 227, 101), 4): ("ordinary", ((4, 0, 2), (0, 4, 2), (2, 2, 6)), ((0, 2, 1), (2, 5, 2), (1, 2, 1)), C7(5, 12, 5), 9_797_760),
    (C7(126, 353, 157), 2): ("ordinary", ((4, -2, 1), (-2, 4, 1), (1, 1, 8)), ((2, 5, 2), (-1, -2, -1), (1, 4, 2)), C7(6, 16, 7), 55_520_640),
    (C7(126, 353, 157), 3): ("ordinary", ((4, 0, 2), (0, 4, -2), (2, -2, 6)), ((3, 7, 3), (-1, -2, -1), (0, 2, 1)), C7(6, 16, 7), 9_797_760),
    (C7(149, 418, 186), 1): ("ordinary", ((4, 0, 0), (0, 4, 2), (0, 2, 4)), ((3, 9, 4), (0, -2, -1), (1, 4, 2)), C7(6, 18, 8), 6_940_080),
    (C7(156, 438, 195), 9): ("radical", ((12, -3, 0), (-3, 12, 3), (0, 3, 18)), ((2, 5, 2), (3, 7, 3), (1, 4, 2)), C7(11, 31, 14), 29_393_280),
    (C7(156, 438, 195), 18): ("radical", ((18, 3), (3, 18)), ((-2, -5, -2), (3, 7, 3)), C7(13, 31, 13), 3_265_920),
    (LIMIT, 31): ("ordinary", ((4, 0, -1), (0, 6, -2), (-1, -2, 6)), ((1, 4, 2), (-3, -7, -3), (-2, -5, -2)), C7(7, 18, 8), 107_775_360),
    (LIMIT, 35): ("ordinary", ((4, 2, 1), (2, 4, -1), (1, -1, 4)), ((1, 4, 2), (-2, -5, -2), (3, 7, 3)), C7(9, 23, 10), 1_693_440),
}


@lru_cache(maxsize=1)
def build():
    arithmetic = h.local_arithmetic(); levels = []
    for level in controls():
        if not (CUT_H < level.height <= LIMIT): continue
        gram_types = types(level)
        if gram_types: levels.append((level, gram_types))
    assert [len(x[1]) for x in levels] == [4, 1, 8, 13, 16, 19, 37]
    corrections = {}
    for key, (lattice, reduced, transport, expected, expected_count) in CASES.items():
        actual, arg, count, explicit = shell_census(lattice, reduced, transport)
        assert (actual, count) == (expected, expected_count)
        corrections[key] = (actual, arg, count, explicit)
    records = []; all_replacement = []
    for level, grams in levels:
        rows = []
        for index, matrix in enumerate(grams):
            necessary, arg, count = dmax(matrix); actual, actual_arg, census, explicit = corrections.get((level.height, index), (necessary, arg, None, None))
            old = -A * level.control + F(7, 5) * actual; replacement = -A * level.control + F(11, 8) * actual; all_replacement.append(replacement)
            rows.append({"index": index, "necessary": necessary.json(), "necessary_arg": arg, "necessary_count": count, "actual": actual.json(), "actual_arg": actual_arg, "shell_count": census, "explicit_coefficients": explicit, "seven_fifths": old.json(), "eleven_eighths": replacement.json()})
        records.append({"height": level.height.json(), "control": level.control.json(), "norm": level.norm.json(), "type_count": len(grams), "types": rows})
    threat = records[-1]["types"][35]; assert C7(**{k:F(v) for k,v in threat["seven_fifths"].items()}).sign() > 0
    assert all(x.sign() < 0 for x in all_replacement)
    owner = json.loads(OWNER.read_text()); assert owner["replacement"]["scale"] == "11/8"
    payload = {"schema":"coxeter.cubic7_k12_h169.hostile.v1","imports_production_search_code":False,"owner_sha256":hashlib.sha256(OWNER.read_bytes()).hexdigest(),"base_cleanroom_sha256":hashlib.sha256(BASE.read_bytes()).hexdigest(),"local_arithmetic":arithmetic,"control_count":len(controls()),"levels":records,"dangerous_shell_case_count":len(CASES),"first_realized_threat":{"height":LIMIT.json(),"type_index":35,"direction":C7(9,23,10).json(),"status":"not_chronologically_accepted"},"replacement":{"scale":"11/8","norm":(C7(F(121,48))-A).json(),"all_necessary_bounds_strict":True},"terminal_verdict":"PASS_BOUNDED_THREAT"}
    return json.loads(json.dumps(payload))


def report(p):
    return """# Independent hostile audit: cubic-seven K12-squared through H169

Terminal verdict: **PASS_BOUNDED_THREAT**.

The replay imports no production search code.  It independently reconstructs
the seven nonempty levels and 98 Gram types, including the scaled ternary gate
at H156.  Eight initially dangerous relaxed types are replayed against the
complete signed K12 or ternary-radical shells.  The first realized positive
7/5 threat is H169 type 35.  Its chronological Vinberg acceptance is not
asserted.  The 11/8 ray is hyperideal and strictly negative against every
necessary type through H169.  No compactness, finite-volume, or termination
claim is made.
"""


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--verify",action="store_true");args=ap.parse_args();p=build();j=json.dumps(p,indent=2,sort_keys=True)+"\n";r=report(p)
    if args.verify:
        if AUDIT.read_text()!=j or REPORT.read_text()!=r:raise SystemExit("hostile H169 artifact drift")
        print("PASS");return
    AUDIT.write_text(j);REPORT.write_text(r);print(AUDIT);print(REPORT)


if __name__=="__main__":main()
