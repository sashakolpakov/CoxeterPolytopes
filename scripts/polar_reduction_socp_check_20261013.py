#!/usr/bin/env python3
"""Independent numerical facet-feasibility check of the reduction criterion.

Uses CVXPY/Clarabel, not the production graph/angle search. The SOCP asks
for a point on each proposed polar wall strictly inside the unit Klein
ball and strictly satisfying every other old/polar inequality.
"""

from itertools import product
import json
from pathlib import Path

import cvxpy as cp
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/polar_reduction_socp_check_20261013.json"


def main():
    gram = np.zeros((11, 11))
    gram[:9, :9] = 2 * np.eye(9)
    for i in range(8):
        gram[i, i + 1] = gram[i + 1, i] = -np.sqrt(2) if i == 0 else -1
    gram[9, 9] = gram[10, 10] = 2
    gram[0, 9] = gram[9, 0] = -np.sqrt(3)
    gram[5, 10] = gram[10, 5] = -np.sqrt(3)
    gram[9, 10] = gram[10, 9] = 6 * np.sqrt(2) - np.sqrt(115)
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    null = eigenvectors[:, np.argmin(abs(eigenvalues))]
    inverse = np.linalg.pinv(gram)
    pairings, omitted = [], []
    for i, j in product(np.where(null > 1e-7)[0], np.where(null < -1e-7)[0]):
        row = np.zeros(11)
        row[i], row[j] = -abs(null[j]), -abs(null[i])
        pairings.append(row)
        omitted.append((int(i), int(j)))
    pairings = np.array(pairings)
    squares = np.diag(pairings @ inverse @ pairings.T)
    hyperideal = np.where(squares > 1e-7)[0]
    columns = [0] + list(np.where(eigenvalues > 1e-7)[0])
    old_roots = eigenvectors[:, columns] @ np.diag(np.sqrt(abs(eigenvalues[columns])))
    old_roots[:, 0] *= np.sign(sum(old_roots[:, 0]))
    vertices = pairings @ inverse @ old_roots
    all_roots = np.vstack([old_roots, vertices[hyperideal]])
    matrix, bounds = all_roots[:, 1:], all_roots[:, 0]
    scale = np.linalg.norm(matrix, axis=1)
    matrix, bounds = matrix / scale[:, None], bounds / scale
    point, margin = cp.Variable(9), cp.Variable()
    checks = []
    for polar_position, vertex in enumerate(hyperideal):
        incident = set(range(11)) - set(omitted[vertex])
        seen, todo = {min(incident)}, [min(incident)]
        while todo:
            i = todo.pop()
            for j in incident - seen:
                if abs(gram[i, j]) > 1e-7:
                    seen.add(j)
                    todo.append(j)
        predicted = seen == incident
        target = 11 + polar_position
        other = [j for j in range(len(all_roots)) if j != target]
        problem = cp.Problem(cp.Maximize(margin), [
            matrix[target] @ point == bounds[target],
            matrix[other] @ point + margin <= bounds[other],
            cp.norm(point, 2) + margin <= 1,
        ])
        optimum = problem.solve(solver="CLARABEL", tol_gap_abs=1e-9,
                                tol_feas=1e-9, tol_gap_rel=1e-9)
        assert predicted == (optimum > 1e-6)
        checks.append({"vertex": int(vertex), "omitted_old_walls": omitted[vertex],
                       "predicted_visible": predicted, "socp_margin": float(optimum)})
    result = {
        "status": "PASS_NUMERICAL_INDEPENDENT_SOCP_CROSSCHECK",
        "proof_status": "NUMERICAL_CHECK_OF_THE_EXACT_THEOREM_NOT_A_PROOF_BY_ITSELF",
        "template": {"chain_labels": [4] + [3] * 7,
                     "attachment_positions": [0, 5], "branch_labels": [6, 6],
                     "dotted_weight": "sqrt(115)-6*sqrt(2)"},
        "hyperideal_cuts_checked": len(checks),
        "visible_polar_facets": sum(row["predicted_visible"] for row in checks),
        "invisible_cuts": sum(not row["predicted_visible"] for row in checks),
        "checks": checks,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "checks"}, indent=2))


if __name__ == "__main__":
    main()
