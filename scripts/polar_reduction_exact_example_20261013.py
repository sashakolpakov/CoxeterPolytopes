#!/usr/bin/env python3
"""Check an exact example where a proposed hyperideal polar cut is invisible.

This checks the Gram/cone data underlying the example, independently of the
floating scan. The geometric visibility theorem is in the companion report.
"""

import json
from pathlib import Path

import sympy as s


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/polar_reduction_exact_example_20261013.json"


def main():
    dotted = s.sqrt(115) - 6 * s.sqrt(2)
    # sqrt(115) > 2 + 6 sqrt(2), using positive-side squaring twice.
    assert 115 - 76 == 39 and 39**2 > 24**2 * 2
    gram = 2 * s.eye(11)
    for i in range(8):
        gram[i, i + 1] = gram[i + 1, i] = -s.sqrt(2) if i == 0 else -1
    gram[0, 9] = gram[9, 0] = -s.sqrt(3)
    gram[5, 10] = gram[10, 5] = -s.sqrt(3)
    gram[9, 10] = gram[10, 9] = -dotted
    base = gram[:9, :9]
    leading = [base[:i, :i].det() for i in range(1, 10)]
    assert leading == [2] * 9  # Sylvester: base is positive definite.
    attachment = gram[:9, 9:]
    schur = (gram[9:, 9:] - attachment.T * base.inv() * attachment).applyfunc(s.simplify)
    assert schur == s.Matrix([[-s.Rational(23, 2), -s.sqrt(115)],
                              [-s.sqrt(115), -10]])
    assert schur.det() == 0 and schur.trace() < 0
    # Thus the complete Gram has inertia (9 positive, 1 negative, 1 zero).
    tail = s.Matrix([2 * s.sqrt(115), -23])
    null = (-base.inv() * attachment * tail).col_join(tail).applyfunc(s.simplify)
    assert (gram * null).applyfunc(s.simplify) == s.zeros(11, 1)
    expected_signs = [1, 1, 1, -1, -1, -1, -1, -1, -1, 1, -1]
    assert [int(s.sign(x)) for x in null] == expected_signs

    # Omitting walls 0,3 gives a feasible cone ray: its two nonzero wall
    # pairings can be z_3<0 and -z_0<0, and their scalar product with z is 0.
    pairings = s.zeros(11, 1)
    pairings[0], pairings[3] = null[3], -null[0]
    assert s.simplify((null.T * pairings)[0]) == 0
    incident = [i for i in range(11) if i not in (0, 3)]
    incident_gram = gram.extract(incident, incident)
    incident_determinant = s.simplify(incident_gram.det())
    assert s.simplify(incident_determinant - 18 * (-4 - dotted**2)) == 0
    # The incident diagram splits into elliptic A2={1,2} and a connected
    # Lorentzian component {4,5,6,7,8,9,10}. Its A5 chain is positive,
    # with remaining Schur complement [[2,-dotted],[-dotted,-2]].
    chain = gram.extract([4, 5, 6, 7, 8], [4, 5, 6, 7, 8])
    branch = gram.extract([4, 5, 6, 7, 8], [9, 10])
    rest = (gram.extract([9, 10], [9, 10]) - branch.T * chain.inv() * branch).applyfunc(s.simplify)
    assert rest == s.Matrix([[2, -dotted], [-dotted, -2]])
    assert [chain[:i, :i].det() for i in range(1, 6)] == [2, 3, 4, 5, 6]
    assert s.simplify(rest.det() + 4 + dotted**2) == 0
    # Incident inertia is (8,1): its orthogonal cone vertex is spacelike.
    result = {
        "status": "PASS_EXACT_REDUNDANT_POLAR_EXAMPLE",
        "ambient_gram_inertia": [9, 1, 1],
        "chain_labels": [4] + [3] * 7,
        "attachment_positions": [0, 5],
        "attachment_labels": [6, 6],
        "dotted_weight": str(dotted),
        "kernel": [str(x) for x in null],
        "kernel_signs": expected_signs,
        "omitted_old_walls": [0, 3],
        "vertex_pairing_vector": [str(x) for x in pairings],
        "incident_gram_inertia": [8, 1, 0],
        "incident_gram_determinant": str(incident_determinant),
        "incident_components": [[1, 2], [4, 5, 6, 7, 8, 9, 10]],
        "conclusion": "A feasible simple hyperideal cone vertex has a disconnected incident diagram, so its polar cut is invisible and redundant inside the hyperbolic ball by the companion theorem.",
        "scope": "An exact example of redundancy, not a compact Coxeter-polytope certificate.",
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(result["status"])


if __name__ == "__main__":
    main()
