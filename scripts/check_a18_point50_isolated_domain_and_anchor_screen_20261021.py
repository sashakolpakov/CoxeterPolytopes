#!/usr/bin/env python3
"""New half-loss isolated domain and necessary comparison with .52 anchors."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
import json
import os
from pathlib import Path

import check_a18_isolated_cap_domain_and_norm_screen_20261021 as exact


ROOT = Path(__file__).resolve().parents[1]
PROOF = "results/arithmetic_whole_a18_point50_cubic_inventory_and_anchor_screen_20261021.md"
OLD_OUTPUT = "results/a18_isolated_cap_domain_norm_screen_20261021.json"
OLD_OUTPUT_SHA = "671550c23549a9f55af602cd5f26578c5901e16646deecd540b328d8010386fc"
PRIMITIVES = "scripts/check_a18_isolated_cap_domain_and_norm_screen_20261021.py"
PRIMITIVES_SHA = "29da7434634c80fb001784d9d5f0d69c3d6247a7b6de1e55e2a0eb6b7ff22932"


def coefficients(degree):
    for t in range(5, 17 if degree == 2 else 18):
        if degree == 2:
            for m in range(1, 9):
                yield {"T": t, "M": m}, exact.trim((m, -t, 1))
        else:
            for s in range(1, 18):
                for m in range(1, 5):
                    yield {"T": t, "S": s, "M": m}, exact.trim((-m, s, -t, 1))


def signature(polynomial, cell):
    return (cell["degree"], polynomial["discriminant_squarefree_part"],
            cell["norm_squarefree_part"])


def census(old):
    counts, polynomials, cells = Counter(), [], []
    grouped = defaultdict(list)
    for degree in (2, 3):
        tag = str(degree) + ":"
        for c, p in coefficients(degree):
            counts[tag + "visited"] += 1
            if not exact.irreducible(p, c, degree):
                counts[tag + "reducible_or_nonreal_quadratic"] += 1
                continue
            sturm = exact.sturm(p)
            if exact.roots(sturm, Q(0), Q(19, 35)) != degree - 1 or exact.roots(sturm, Q(4), Q(16)) != 1:
                counts[tag + "outside_root_domain"] += 1
                continue
            discriminant = exact.polynomial_discriminant(c, degree)
            assert discriminant > 0
            part, square = exact.squarefree(discriminant)
            polynomial = {"id": len(polynomials), "degree": degree,
                "coefficients": c, "ascending": [int(x) for x in p],
                "discriminant": discriminant, "discriminant_squarefree_part": part,
                "discriminant_square_factor": square}
            polynomials.append(polynomial)
            counts[tag + "retained_polynomials"] += 1
            for i in range(5, 15):
                counts[tag + "index_tests"] += 1
                k = i * (19 - i)
                if exact.roots(sturm, Q(0), Q(38, k)) != degree - 1:
                    counts[tag + "failed_index_schur"] += 1
                    continue
                numerator = exact.norm_numerator(c, degree, k)
                assert numerator and (numerator > 0) == (degree == 3)
                assert exact.multiplication_norm(p, k) == Q(numerator, 19**degree)
                part, square = exact.squarefree(abs(numerator))
                cell = {"id": len(cells), "polynomial": polynomial["id"],
                    "degree": degree, "index": i, "k": k,
                    "norm_numerator": numerator, "norm_denominator": 19**degree,
                    "norm_squarefree_part": part, "norm_square_factor": square}
                cells.append(cell)
                grouped[signature(polynomial, cell)].append(cell["id"])
                counts[tag + "retained_cells_and_norm_checks"] += 1
        assert counts[tag + "visited"] == (96 if degree == 2 else 884)
        assert counts[tag + "visited"] == sum(counts[tag + name] for name in (
            "reducible_or_nonreal_quadratic", "outside_root_domain", "retained_polynomials"))
        assert counts[tag + "index_tests"] == 10 * counts[tag + "retained_polynomials"]
        assert counts[tag + "index_tests"] == counts[tag + "failed_index_schur"] + counts[tag + "retained_cells_and_norm_checks"]
    matches = []
    for anchor in old["cells"]:
        if anchor["degree"] != 3:
            continue
        polynomial = old["polynomials"][anchor["polynomial"]]
        ids = grouped.get(signature(polynomial, anchor), [])
        assert any(cells[j]["index"] == anchor["index"] and
                   polynomials[cells[j]["polynomial"]]["coefficients"] == polynomial["coefficients"] for j in ids)
        matches.append({"anchor_cell": anchor["id"], "target_cells": ids,
            "only_same_polynomial_and_k": all(cells[j]["k"] == anchor["k"] and
                polynomials[cells[j]["polynomial"]]["coefficients"] == polynomial["coefficients"] for j in ids)})
    groups = [{"degree": key[0], "discriminant_squarefree_part": key[1],
               "norm_squarefree_part": key[2], "cells": ids}
              for key, ids in sorted(grouped.items())]
    return {"counts": dict(sorted(counts.items())), "polynomials": polynomials,
            "cells": cells, "necessary_pair_groups": groups,
            "point52_cubic_anchor_matches": matches,
            "all_cubic_anchor_matches_only_same_polynomial_and_k":
                all(m["only_same_polynomial_and_k"] for m in matches)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    assert os.environ.get("COXETER_CHECK_OWNER"), "Use the bounded parent-lock runner"
    output = (ROOT / args.out).resolve()
    assert output.parent == ROOT / "results" and not output.exists()
    paths = [ROOT / PROOF, ROOT / exact.PROOF, ROOT / exact.DOMAIN,
             ROOT / exact.RATIONAL_OBSTRUCTION, ROOT / OLD_OUTPUT,
             ROOT / PRIMITIVES, Path(__file__).resolve()]
    pins = [{"path": str(p.relative_to(ROOT)), "sha256": exact.digest(p),
             "bytes": p.stat().st_size} for p in paths]
    assert pins[2]["sha256"] == exact.DOMAIN_SHA
    assert pins[3]["sha256"] == exact.RATIONAL_OBSTRUCTION_SHA
    assert pins[4]["sha256"] == OLD_OUTPUT_SHA
    assert pins[5]["sha256"] == PRIMITIVES_SHA
    old = json.loads((ROOT / OLD_OUTPUT).read_text())
    result = census(old)
    for path, pin in zip(paths, pins):
        assert exact.digest(path) == pin["sha256"] and path.stat().st_size == pin["bytes"]
    result.update({"status": "EXACT_POINT50_NECESSARY_DOMAIN_AND_ANCHOR_SCREEN_COMPLETE",
        "input_pins": pins,
        "scope": "Necessary isolated caps only; no K-square, realization or energy conclusion is asserted by this worker."})
    with output.open("x") as stream:
        json.dump(result, stream, separators=(",", ":"), sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "counts": result["counts"],
        "cubic_anchor_matches": len(result["point52_cubic_anchor_matches"]),
        "all_cubic_anchor_matches_only_same_polynomial_and_k":
            result["all_cubic_anchor_matches_only_same_polynomial_and_k"],
        "output": str(output.relative_to(ROOT)), "output_sha256": exact.digest(output)}))


if __name__ == "__main__":
    main()
