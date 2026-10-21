#!/usr/bin/env python3
"""Find and exactly verify escapes from selected H20 boundary intersections.

LP only proposes vectors. Every reported escape is an integer vector in
f-perp with positive square and STRICTLY negative pairing with every
selected wall. A strict timelike point is recorded as well. Failed searches
are inconclusive; there is no finite-volume assertion from sampling.
"""

from fractions import Fraction
from functools import reduce
from math import lcm
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog

import h20_gluing_independent_audit_20261013 as audit

ROOT = Path(__file__).resolve().parents[1]
BOUNDARY = ROOT / "data/enumerations/new_dimension_h20_six_patch_boundary_20261013.json"
CANDIDATES = ROOT / "data/enumerations/new_dimension_h20_compatible_subsets_20261013.json"
OUTPUT = ROOT / "data/enumerations/h20_gluing_independent_escape_20261013.json"
CUSPS = ROOT / "data/enumerations/h20_gluing_independent_cusps_20261013.json"


def setup():
    Q, W, G, h = audit.geometry()
    f = W[0]
    center = np.array(audit.primitive(2*h - (h @ Q @ f)*f), dtype=object)
    normals = np.array(json.loads(BOUNDARY.read_text())["six_tile"]["primitive_outward_normals"], dtype=object)
    assert center @ Q @ center < 0 and center @ Q @ f == 0
    assert max(normals @ Q @ center) < 0
    E = np.array([f @ Q, center @ Q], dtype=float)
    spatial = null_space(E)
    assert spatial.shape == (22, 20)
    return Q, f, center, normals, spatial


def corrected_integer_escape(point, Q, f, center, normals):
    raw = np.array([int(round(float(x) * 10**7)) for x in point], dtype=object)
    vector = 2*raw - (raw @ Q @ f)*f
    center_pairings = normals @ Q @ center
    values = normals @ Q @ vector
    shift = max([0] + [int(v) // int(-d) + 1 for v, d in zip(values, center_pairings)])
    vector += shift * center
    vector = np.array(audit.primitive(vector), dtype=object)
    if vector @ Q @ vector <= 0:
        return None
    values = normals @ Q @ vector
    assert vector @ Q @ f == 0 and max(values) < 0
    return {"vector": list(map(int, vector)), "square": int(vector @ Q @ vector),
            "maximum_selected_wall_pairing": int(max(values)),
            "strict_center_correction": shift}


def rational_cusp(point, Q, f, center, normals):
    scale = max(abs(float(x)) for x in point)
    rationals = [Fraction(float(x)/scale).limit_denominator(100000) for x in point]
    denominator = reduce(lcm, (x.denominator for x in rationals), 1)
    vector = np.array([int(x*denominator) for x in rationals], dtype=object)
    if vector @ Q @ center > 0:
        vector = -vector
    if vector @ Q @ f != 0 or vector @ Q @ vector != 0 or max(normals @ Q @ vector) > 0:
        return None
    return np.array(audit.primitive(vector), dtype=object)


def proposals(normals, Q, center, spatial, objectives=80, seed=113):
    A = np.array(normals @ Q, dtype=float) @ spatial
    b = -np.array(normals @ Q @ center, dtype=float)
    scale = np.linalg.norm(A, axis=1)
    A, b = A / scale[:, None], b / scale
    rng = np.random.default_rng(seed)
    directions = list(np.eye(20)) + list(-np.eye(20))
    directions += list(rng.normal(size=(max(0, objectives-40), 20)))
    for i, direction in enumerate(directions[:objectives]):
        result = linprog(-direction, A_ub=A, b_ub=b, bounds=(None, None), method="highs")
        if result.status == 3:
            result = linprog(-direction, A_ub=A, b_ub=b, bounds=(-10000, 10000), method="highs")
        if result.success:
            yield i, np.array(center, dtype=float) + spatial @ result.x


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--discover-cusps", type=int, default=0)
    parser.add_argument("--objectives", type=int, default=80)
    args = parser.parse_args()
    Q, f, center, all_normals, spatial = setup()
    if args.discover_cusps:
        seen = {}
        for objective, point in proposals(all_normals, Q, center, spatial, args.discover_cusps):
            cusp = rational_cusp(point, Q, f, center, all_normals)
            if cusp is not None:
                pairs = all_normals @ Q @ cusp
                seen[tuple(cusp)] = {"vector": list(map(int, cusp)), "square": 0,
                    "incident_boundary_indices": [int(i) for i in np.where(pairs == 0)[0]]}
        result = {"status": "EXACT_IDEAL_POINTS_FOUND_BY_BOUNDED_LP_SAMPLING",
                  "objectives": args.discover_cusps, "cusps": list(seen.values())}
        CUSPS.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps({"cusps": len(seen), "output": str(CUSPS)}), flush=True)
        return
    source = json.loads(CANDIDATES.read_text())
    results = []
    for number, candidate in enumerate(source["candidates"]):
        ids = candidate["indices_zero_based"]
        normals = all_normals[ids]
        witness = None
        for objective, point in proposals(normals, Q, center, spatial, args.objectives, 113+number):
            witness = corrected_integer_escape(point, Q, f, center, normals)
            if witness is not None:
                witness["objective_index"] = objective
                break
        row = {"candidate_id": candidate.get("id", candidate.get("name", str(number))),
               "indices_zero_based": ids, "wall_count": len(ids),
               "status": "EXACT_STRICT_SPACELIKE_ESCAPE" if witness else "INCONCLUSIVE_BOUNDED_SAMPLING",
               "witness": witness}
        results.append(row)
        print(json.dumps({k: v for k, v in row.items() if k not in ("indices_zero_based", "witness")}), flush=True)
    result = {"source_candidates": str(CANDIDATES.relative_to(ROOT)),
        "source_boundary": str(BOUNDARY.relative_to(ROOT)),
        "common_strict_timelike_center": list(map(int, center)),
        "center_square": int(center @ Q @ center), "results": results,
        "scope": "An exact strict spacelike escape proves infinite volume for that literal wall intersection. No-witness outcomes are inconclusive; compatibility selection is not exhaustive."}
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
