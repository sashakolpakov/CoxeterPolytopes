#!/usr/bin/env python3
"""Exact finite arithmetic checks for rigidity of a literal Bugaenko facet.

The proof's universal step is monotonicity of sin(pi/m) for integer m>=5;
the remaining comparisons are exact in Q(phi).  The statement concerns a
literal Bugaenko facet, not every compact H9 Coxeter polytope.
"""

from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/known_examples/bugaenko_compact_h8.json"
OUTPUT = ROOT / "data/enumerations/fresh_bugaenko_facet_tilt_rigidity_20261013.json"


def verify():
    data = json.loads(SOURCE.read_text())
    relations = data["coxeter_relations"]
    gram = [[Golden.from_json(x) for x in row] for row in data["wall_gram"]]
    unseen = set(range(11))
    components = []
    while unseen:
        current = {min(unseen)}
        todo = list(current)
        while todo:
            i = todo.pop()
            for j in range(11):
                if relations[i][j] == 3 and j not in current:
                    current.add(j)
                    todo.append(j)
        unseen -= current
        components.append(sorted(current))
    assert sorted(components) == [[0], [1, 2, 3, 4, 5, 6, 7, 9, 10], [8]]
    assert gram[1][3] == Golden(0)
    for endpoint, neighbor in ((0, 1), (8, 7)):
        assert relations[endpoint][neighbor] == 5
        assert all(relations[endpoint][j] == 2 for j in range(11)
                   if j not in (endpoint, neighbor))
    phi = Golden(0, 1)
    m3_square = 3 * phi * phi / 4
    m4_square = phi * phi / 2
    m5_upper_square = phi * phi * (3 - phi) / 4
    assert (m3_square - 1).sign() > 0 and (2 - m3_square).sign() > 0
    assert (m4_square - 1).sign() > 0 and (2 - m4_square).sign() > 0
    assert m5_upper_square.sign() > 0 and (1 - m5_upper_square).sign() > 0
    return {
        "status": "PASS_EXACT_LITERAL_FACET_TILT_RIGIDITY",
        "source_sha256": sha256(SOURCE.read_bytes()).hexdigest(),
        "label_three_components": components,
        "orthogonal_pair_in_large_component": [1, 3],
        "normalization": {
            "root_square": 2,
            "finite_gram_pairing": "-2*cos(pi/k)",
            "endpoint_remaining_pairing_squares": "Squared Gram entries: 4*cos(pi/k)^2, not cos(pi/k)^2.",
        },
        "endpoint_remaining_pairing_squares": {
            "base_angle_pi_over_3": m3_square.to_json(),
            "base_angle_pi_over_4": m4_square.to_json(),
            "base_angle_pi_over_m_for_m_at_least_5_upper_bound": m5_upper_square.to_json(),
        },
        "conclusion": "Every facet adjacent to a literal Bugaenko H8 facet in a hyperbolic Coxeter polytope meets that facet orthogonally.",
        "scope": "This does not exclude compact H9 polytopes, even those with a right Bugaenko facet.",
    }


if __name__ == "__main__":
    result = verify()
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
