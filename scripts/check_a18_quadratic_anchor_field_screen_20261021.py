#!/usr/bin/env python3
"""Join the two reviewed exact quadratic isolated-cap domains, without rerunning them."""

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path

import check_a18_isolated_cap_domain_and_norm_screen_20261021 as exact


ROOT = Path(__file__).resolve().parents[1]
PROOF = "results/arithmetic_whole_a18_quadratic_anchor_field_screen_20261021.md"
FIXED = {
    "results/a18_isolated_cap_domain_norm_screen_20261021.json":
        "671550c23549a9f55af602cd5f26578c5901e16646deecd540b328d8010386fc",
    "results/a18_point50_isolated_domain_anchor_screen_20261021.json":
        "8c3b9ffb07f7a6aef39463b464c3a5c18afe5fbd07929eb2cdf0b0e3b820beae",
    "results/arithmetic_whole_a18_point52_actual_field_degree_reduction_20261021.md":
        exact.DOMAIN_SHA,
    "results/arithmetic_whole_a18_isolated_cap_domain_and_norm_screen_20261021.md":
        "06172cc3d1f1e218a29d7e20916f105366d10344f4e086a4aeb4d4eed45e75ac",
    "results/arithmetic_whole_a18_point50_cubic_inventory_and_anchor_screen_20261021.md":
        "bcefd51d39d814b11a37eb04dd08e17962a1e71caee932ebc0e35af10dd8d125",
    "scripts/check_a18_isolated_cap_domain_and_norm_screen_20261021.py":
        "29da7434634c80fb001784d9d5f0d69c3d6247a7b6de1e55e2a0eb6b7ff22932",
    "scripts/check_a18_point50_isolated_domain_and_anchor_screen_20261021.py":
        "2ef5c3651c1a0883434b3bbc8b45542cb8593629731299fa65db2d501a98f60c",
}


def signature(poly, cell):
    return (cell["degree"], poly["discriminant_squarefree_part"],
            cell["norm_squarefree_part"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    assert os.environ.get("COXETER_CHECK_OWNER"), "Use the bounded parent-lock runner"
    output = (ROOT / args.out).resolve()
    assert output.parent == ROOT / "results" and not output.exists()
    paths = [ROOT / p for p in FIXED] + [ROOT / PROOF, Path(__file__).resolve()]
    pins = [{"path": str(p.relative_to(ROOT)), "sha256": exact.digest(p),
             "bytes": p.stat().st_size} for p in paths]
    assert all(pin["sha256"] == FIXED[pin["path"]] for pin in pins if pin["path"] in FIXED)
    old = json.loads(paths[0].read_text())
    new = json.loads(paths[1].read_text())
    groups = defaultdict(list)
    for cell in new["cells"]:
        if cell["degree"] == 2:
            groups[signature(new["polynomials"][cell["polynomial"]], cell)].append(cell["id"])
    matches = []
    for anchor in old["cells"]:
        if anchor["degree"] != 2:
            continue
        poly = old["polynomials"][anchor["polynomial"]]
        key = signature(poly, anchor)
        ids = groups.get(key, [])
        assert any(new["cells"][j]["index"] == anchor["index"] and
                   new["polynomials"][new["cells"][j]["polynomial"]]["coefficients"] == poly["coefficients"] for j in ids)
        matches.append({"anchor_cell": anchor["id"], "field_squarefree_D": key[1],
            "target_cells": ids,
            "only_same_polynomial_and_k": all(new["cells"][j]["k"] == anchor["k"] and
                new["polynomials"][new["cells"][j]["polynomial"]]["coefficients"] == poly["coefficients"] for j in ids)})
    outside = [m for m in matches if m["field_squarefree_D"] not in (2, 5)]
    result = {"status": "EXACT_SAVED_QUADRATIC_ANCHOR_COMPARISON_COMPLETE",
        "input_pins": pins, "anchor_matches": matches,
        "quadratic_anchor_count": len(matches), "outside_silver_golden_anchor_count": len(outside),
        "outside_silver_golden_only_same_polynomial_and_k":
            all(m["only_same_polynomial_and_k"] for m in outside),
        "scope": "Necessary saved-data match only; field restriction also requires the written complete nonisolated inventory, reflected integrality and actual-foot energy proof."}
    for p, pin in zip(paths, pins):
        assert exact.digest(p) == pin["sha256"] and p.stat().st_size == pin["bytes"]
    with output.open("x") as stream:
        json.dump(result, stream, separators=(",", ":"), sort_keys=True)
        stream.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("input_pins", "anchor_matches")}))
    print(json.dumps({"output": str(output.relative_to(ROOT)), "sha256": exact.digest(output)}))


if __name__ == "__main__":
    main()
