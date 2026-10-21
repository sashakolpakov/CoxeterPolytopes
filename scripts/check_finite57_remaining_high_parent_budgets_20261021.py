#!/usr/bin/env python3
"""Exact arithmetic for the remaining high-parent Context662 source families.

Replay only: six cap-domain terminal rows, all fifteen (p,q) grid rows,
and one explicitly negative relaxation diagnostic. No geometry or forest
classification is inferred. Run through the shared bounded checker wrapper.
"""

import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from math import comb, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PINS = {
    "manuscripts/finite_volume57.tex":
        "c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0",
    "results/finite57_complete_high_output_finite_classification_20261021.md":
        "162965f7c47aeafb1d0df7bb2181f53ff50e79ee2f3256a9a44b135a68a0261b",
}
CANONICAL = {
    14: F(367, 60), 15: F(411, 64), 16: F(227, 34),
    17: F(119, 18), 18: F(128, 19), 19: F(137, 20),
    21: F(155, 22), 22: F(163, 23), 23: F(85, 12),
    24: F(176, 25), 25: F(181, 26),
}


def require(predicate, message):
    if not predicate:
        raise ValueError(message)


def check_pins():
    for relative, expected in PINS.items():
        digest = hashlib.sha256()
        with (ROOT / relative).open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        require(digest.hexdigest() == expected, "Changed input: " + relative)


def coefficients(ranks):
    h = [1]
    for r in ranks:
        require(isinstance(r, int) and r >= 0, "Invalid simplex rank")
        new = [0] * (len(h) + r)
        for i, value in enumerate(h):
            for j in range(r + 1):
                new[i + j] += value
        h = new
    return h


def inclusion_exclusion(ranks):
    """Coefficient of H/(1-t)^2; retain every numerator subset."""
    answer = 0
    denominator_power = len(ranks) + 2
    for chosen in itertools.product((0, 1), repeat=len(ranks)):
        shift = sum(bit * (r + 1) for bit, r in zip(chosen, ranks))
        degree = 27 - shift
        if degree >= 0:
            term = comb(degree + denominator_power - 1,
                        denominator_power - 1)
            answer += (-1) ** sum(chosen) * term
    return answer


def evaluate(name, ranks, normalized_cost, debt):
    require(len(ranks) == 4 and sum(ranks) == 56, "Incomplete rank profile")
    h = coefficients(ranks)
    require(len(h) == 57 and h == h[::-1], "Wrong link polynomial")
    link_vertex_count = prod(r + 1 for r in ranks)
    require(sum(h) == link_vertex_count, "Wrong link vertex count")
    a = sum((28 - j) * h[j] for j in range(28))
    require(a == inclusion_exclusion(ranks), "Independent coefficient mismatch")
    g = sum(min(j, 56 - j) * value for j, value in enumerate(h))
    require(g == 28 * link_vertex_count - 2 * a, "Wrong index identity")
    kappa = 42 * link_vertex_count - 5 * a
    require(F(5 * g - 56 * link_vertex_count, 2) == kappa,
            "Wrong curvature identity")
    envelope = link_vertex_count * normalized_cost
    surplus = kappa - envelope - debt
    return {
        "name": name, "ranks": ranks, "link_vertices": link_vertex_count,
        "coefficient_a": a, "kappa": kappa,
        "envelope": str(envelope), "debt": str(debt),
        "surplus": str(surplus),
    }, surplus


def replay():
    check_pins()
    safe = []
    for p, k, s, small_cost, expected in [
        (25, 14, 1, F(15, 2), F(6815)),
        (24, 14, 2, F(16, 2), F(13815, 2)),
        (23, 14, 3, F(17, 2), F(5250)),
        (22, 14, 4, F(18, 2), F(4495, 2)),
        (21, 14, 5, CANONICAL[14] + F(5, 2), F(28098)),
        (21, 13, 6, F(19, 2), F(3549)),
    ]:
        row, surplus = evaluate(
            "safe-p%d-cap%d" % (p, k), [16, p, k, s],
            CANONICAL[16] + CANONICAL[p] + small_cost, F(0))
        require(surplus == expected, "Safe terminal value changed")
        safe.append(row)
    require(len(safe) == 6 and len({r["name"] for r in safe}) == 6,
            "Incomplete safe terminal domain")
    require(min(F(r["surplus"]) for r in safe) == F(4495, 2),
            "Wrong safe minimum")

    grid = []
    for p in range(21, 26):
        for q in range(15, min(19, 40 - p) + 1):
            m = 40 - p - q
            require(0 <= m <= 4, "Wrong small mass")
            d = 2 if q == 16 else 1
            high_ports = 4 if p == 21 else 2 * p - 40
            debt = F(32 * d * high_ports * (q + 1) * (m + 1))
            row, surplus = evaluate(
                "grid-p%d-q%d" % (p, q), [16, p, q, m],
                CANONICAL[16] + CANONICAL[p] + CANONICAL[q] + F(m, 2), debt)
            row.update({"varying_d16_upper": d, "large_parent_high_ports": high_ports})
            require(surplus > 0, "Nonpositive grid source: " + row["name"])
            require(p + 1 <= 8 * high_ports,
                    "Affine-A all-port debt not dominated")
            grid.append(row)
    require(len(grid) == 15 and len({r["name"] for r in grid}) == 15,
            "Incomplete grid domain")
    require(min(F(r["surplus"]) for r in grid) == 167,
            "Wrong grid minimum")

    diagnostic, relaxation = evaluate(
        "negative-coarse-p21-cap14", [16, 21, 14, 5],
        CANONICAL[16] + CANONICAL[21] + F(19, 2), F(0))
    require(relaxation == -1635, "Negative relaxation diagnostic changed")
    check_pins()
    return {
        "status": "PASS_EXACT_SOURCE_ARITHMETIC_ONLY",
        "required_positive_rows": 21, "negative_diagnostics": 1,
        "safe_minimum": "4495/2", "grid_minimum": "167",
        "safe_rows": safe, "grid_rows": grid, "diagnostic": diagnostic,
        "source_pins": PINS,
        "scope": "Finite exact arithmetic; geometry, complete hard lists, capped concentration, actual incidence and global allocation remain written obligations.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", required=True)
    parser.parse_args()
    print(json.dumps(replay(), sort_keys=True))


if __name__ == "__main__":
    main()
