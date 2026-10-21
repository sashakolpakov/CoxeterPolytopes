#!/usr/bin/env python3
"""Certify and identify the unique high-rank noncyclic Borcherds section.

Every element of a subgroup whose common fixed lattice has rank at least 11
must itself lie in one of the four high cyclic classes 2a, 2c, 2d, or 3b.
The Reynolds trace bound leaves one possible noncyclic fusion pattern:
V4 with one 2a and two 2c involutions.  This script proves that GAP finds one
such subgroup class, replays a pinned representative on all 210 walls,
certifies its H10/34 Coxeter chamber and literal Vinberg volume criterion, and
gives an exact normalized-Gram isomorphism to A15-C2-H10-34.
"""

from __future__ import annotations

import argparse
import base64
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any

import networkx as nx
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_alexeev_high_dim_reflection_examples import (
    SOURCE_SHA256,
    _extract_matrix,
)
from scripts.certify_borcherds_char2_3b_h11_22 import (
    DEFAULT_SOURCE,
    _digest,
    _matrix_payload,
)
from scripts.certify_borcherds_char2_all_cyclic_elements import (
    _ambient_action,
    _canonical_generators,
    _facet_lattice_records,
    _gap_permutation,
    _lift_levi_permutation,
)
from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    CPP_SOURCE,
    _audit_case,
)
from scripts.certify_shimada_disc16_d5_h14_fixed_lattice import (
    _smith_decomposition,
)
from scripts.search_k3_high_rank_involution_folds import _coxeter_audit
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _orbit_sum_matrix


ROOT = Path(__file__).resolve().parents[1]
TARGET_CENSUS = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_fixed_section_census.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "borcherds_char2_v4_h10_34_certificate.json"
)
SCHEMA = "coxeter.borcherds_char2_v4_h10_34.v1"

# These are zero-based images packed as bytes.  The first involution is 2a;
# the second and their product are 2c.  GAP independently verifies this
# fusion against the reconstructed full diagram group on every replay.
CLASS_2A_BASE64 = (
    "GzAQPwQFBiAeFhwSKA0uGAIyC4qjHwkXDxkaAAohCBUHHTGMWLUmJwxJKr/LLQ5U"
    "ASIRtDSEjaKBX7DOr1HMA3BCQba7Rb1Wsimee0yliU/DPc+3L2lHxCS6iJqTXZ05x"
    "6tzY4+Xi8KqVc1rbG3FhUBxcmKbobx8oIB6S3efh395ONCDNW+ZflpOE2YjNo5kvp"
    "KRXJSVmGWWhlt0qF5KfXh1NxS5TcqtnMBoYcinrjw6sUjRMyVDU7ikWUR2RpArqcF"
    "nUFduxmCsyaYsPmo7UoKz"
)
CLASS_2C_BASE64 = (
    "AAEIPwQFLRgCIh0SKA0OIB4yC3iIFTEXBxoZGyEKEB8PHAliv0kmJwy1g1hvBi53"
    "MBYRkZSEzJ2YOUpRxc6NA1K9RpqTrkJoqSU6p0xqiU9VO0CbfFCqrCuho7a7caJfv"
    "KsjY755ymdHw03JesGvLM9dcoy3uscvE2VsrVSCfpWXln0qNcuzhxROoKZzPo6Qjz"
    "O0RDR/gYA40UNTnDew0IpZXlq5zYtLqEhWYVd7RW6escCGkilbdMakdVxgQWQksm3"
    "Cacg8uHbEa2aFNqU9cJ+Z"
)

EXPECTED_COXETER = {
    "2": 331,
    "3": 17,
    "4": 30,
    "6": 0,
    "infinity": 75,
    "ultraparallel": 108,
}
EXPECTED_FACET_TYPES = Counter(
    {
        (1, -2, 1): 6,
        (1, -4, 2): 4,
        (1, -4, 4): 2,
        (2, -4, 2): 10,
        (2, -4, 4): 8,
        (4, -8, 4): 4,
    }
)
EXPECTED_VINBERG = {
    "dimension": 10,
    "facets": 34,
    "finite_special_counts": [
        34, 378, 2280, 8725, 22392, 38931, 44922, 32466, 12792, 1774,
    ],
    "affine_components": 426,
    "cusps": 146,
    "endpoint_histogram": {"0,2": 1098, "1,1": 5648, "2,0": 6046},
}


def _permutation(encoded: str) -> tuple[int, ...]:
    result = tuple(base64.b64decode(encoded))
    if len(result) != 210 or set(result) != set(range(210)):
        raise ArithmeticError("a pinned V4 generator is malformed")
    return result


def _compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def _subgroup_elements() -> tuple[tuple[int, ...], ...]:
    identity = tuple(range(210))
    class_2a = _permutation(CLASS_2A_BASE64)
    class_2c = _permutation(CLASS_2C_BASE64)
    product = _compose(class_2a, class_2c)
    elements = (identity, class_2a, class_2c, product)
    if not (
        _compose(class_2a, class_2a) == identity
        and _compose(class_2c, class_2c) == identity
        and _compose(class_2a, class_2c) == _compose(class_2c, class_2a)
        and len(set(elements)) == 4
    ):
        raise ArithmeticError("the pinned generators do not form V4")
    return elements


def _wall_orbits(elements: tuple[tuple[int, ...], ...]) -> list[list[int]]:
    unseen = set(range(210))
    result = []
    while unseen:
        start = min(unseen)
        orbit = sorted({permutation[start] for permutation in elements})
        unseen.difference_update(orbit)
        result.append([value + 1 for value in orbit])
    return result


def _gap_subgroup_proof(
    generators: dict[str, tuple[int, ...]],
    elements: tuple[tuple[int, ...], ...],
    gap_executable: str,
) -> dict[str, Any]:
    definitions = ",\n".join(_gap_permutation(value) for value in generators.values())
    pinned_2a = _gap_permutation(elements[1])
    pinned_2c = _gap_permutation(elements[2])
    program = f"""
SizeScreen([100000,100000]);;
Reset(GlobalRandomSource,1);;
Reset(GlobalMersenneTwister,1);;
gens := [{definitions}];;
G := Group(gens);;
Hgeom := Group([gens[1],gens[2],gens[3],gens[5]]);;
classes := ConjugacyClasses(G);;
ia := PositionProperty(classes,c->Order(Representative(c))=2 and Size(c)=315);;
ic := PositionProperty(classes,c->Order(Representative(c))=2 and Size(c)=360);;
a := Representative(classes[ia]);;
C := Centralizer(G,a);;
valid := Filtered(Elements(C),x->x in classes[ic] and a*x in classes[ic]);;
subs := [];;
for x in valid do
  s := Group([a,x]);;
  new := true;;
  for t in subs do if t=s then new:=false; fi; od;
  if new then Add(subs,s); fi;
od;
action := function(s,g) return s^g; end;;
orbits := OrbitsDomain(C,subs,action);;
s := subs[1];;
N := Normalizer(G,s);;
pa := {pinned_2a};;
pc := {pinned_2c};;
Print("PROOF|",Size(G),"|",Size(Hgeom),"|",Length(classes),"|",Size(C),"|",
      Length(valid),"|",Length(subs),"|",Length(orbits),"|",Size(N),"|",
      Size(G)/Size(N),"|",pa in G,"|",pc in G,"|",pa in classes[ia],"|",
      pc in classes[ic],"|",pa*pc in classes[ic],"|",pa in Hgeom,"|",
      pc in Hgeom,"|",pa*pc in Hgeom,"\\n");
"""
    completed = subprocess.run(
        [gap_executable, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    line = next(
        (candidate for candidate in completed.stdout.splitlines() if candidate.startswith("PROOF|")),
        None,
    )
    if line is None:
        raise ArithmeticError("GAP did not return the V4 subgroup proof")
    fields = line.split("|")
    record = {
        "full_group_order": int(fields[1]),
        "geometric_subgroup_order": int(fields[2]),
        "full_group_element_class_count": int(fields[3]),
        "centralizer_order_of_2a": int(fields[4]),
        "compatible_2c_elements_in_centralizer": int(fields[5]),
        "distinct_V4_subgroups_containing_fixed_2a": int(fields[6]),
        "centralizer_orbits_on_those_V4_subgroups": int(fields[7]),
        "V4_normalizer_order": int(fields[8]),
        "V4_full_group_conjugacy_class_size": int(fields[9]),
        "pinned_generators_in_full_group": fields[10:12] == ["true", "true"],
        "pinned_fusion_2a_2c_2c": fields[12:15] == ["true", "true", "true"],
        "geometric_membership_2a_2c_product": [
            value == "true" for value in fields[15:18]
        ],
    }
    expected = {
        "full_group_order": 241920,
        "geometric_subgroup_order": 120960,
        "full_group_element_class_count": 28,
        "centralizer_order_of_2a": 768,
        "compatible_2c_elements_in_centralizer": 24,
        "distinct_V4_subgroups_containing_fixed_2a": 12,
        "centralizer_orbits_on_those_V4_subgroups": 1,
        "V4_normalizer_order": 64,
        "V4_full_group_conjugacy_class_size": 3780,
        "pinned_generators_in_full_group": True,
        "pinned_fusion_2a_2c_2c": True,
        "geometric_membership_2a_2c_product": [True, False, False],
    }
    if record != expected:
        raise ArithmeticError(f"the V4 subgroup census changed: {record}")
    return record


def _common_fixed_lattice(
    basis_gram: sp.Matrix, actions: tuple[sp.Matrix, ...]
) -> tuple[sp.Matrix, sp.Matrix, list[int]]:
    constraint = sp.Matrix.vstack(*(action - sp.eye(22) for action in actions))
    diagonal, _, right = _smith_decomposition(constraint)
    zero_columns = tuple(
        column
        for column in range(22)
        if all(diagonal[row, column] == 0 for row in range(diagonal.rows))
    )
    fixed_basis = right[:, list(zero_columns)]
    if not all(action * fixed_basis == fixed_basis for action in actions):
        raise ArithmeticError("the common Smith kernel is not fixed")
    fixed_gram = fixed_basis.T * basis_gram * fixed_basis
    smith = smith_normal_form(fixed_gram, domain=ZZ)
    smith_absolute = [
        abs(int(smith[index, index])) for index in range(fixed_gram.rows)
    ]
    if not (
        fixed_gram.rows == 11
        and exact_inertia(fixed_gram) == (1, 10, 0)
        and int(fixed_gram.det()) == 1024
        and smith_absolute == [1] * 4 + [2] * 4 + [4] * 3
    ):
        raise ArithmeticError("the common V4 fixed lattice changed")
    return fixed_basis, fixed_gram, smith_absolute


def _normalized_ratio(matrix: sp.Matrix, left: int, right: int) -> sp.Rational:
    if left == right:
        return sp.Rational(1)
    return sp.factor(
        sp.Rational(matrix[left, right]) ** 2
        / (sp.Rational(matrix[left, left]) * sp.Rational(matrix[right, right]))
    )


def _normalized_graph(matrix: sp.Matrix) -> nx.Graph:
    graph = nx.Graph()
    for left in range(matrix.rows):
        signature = Counter(
            str(_normalized_ratio(matrix, left, right))
            for right in range(matrix.rows)
            if right != left
        )
        graph.add_node(left, signature=tuple(sorted(signature.items())))
    for left in range(matrix.rows):
        for right in range(left):
            graph.add_edge(
                left,
                right,
                ratio=str(_normalized_ratio(matrix, left, right)),
            )
    return graph


def _deduplicate_to_a15(
    facet_gram: sp.Matrix, target_path: Path
) -> dict[str, Any]:
    target_payload = json.loads(target_path.read_text(encoding="utf-8"))
    target_case = next(
        case for case in target_payload["cases"] if case["id"] == "A15_9_1_C01"
    )
    target_gram = sp.Matrix(target_case["restricted_mirrors"]["facet_gram"])
    if target_gram.shape != facet_gram.shape:
        raise ArithmeticError("the proposed A15 target no longer has 34 facets")
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _normalized_graph(facet_gram),
        _normalized_graph(target_gram),
        node_match=lambda left, right: left["signature"] == right["signature"],
        edge_match=lambda left, right: left["ratio"] == right["ratio"],
    )
    if not matcher.is_isomorphic():
        raise ArithmeticError("the V4 chamber no longer matches A15-C2-H10-34")
    mapping = matcher.mapping
    if not all(
        _normalized_ratio(facet_gram, left, right)
        == _normalized_ratio(target_gram, mapping[left], mapping[right])
        for left in range(34)
        for right in range(34)
    ):
        raise ArithmeticError("the proposed normalized-Gram permutation is not exact")
    source_ratios = [
        [str(_normalized_ratio(facet_gram, left, right)) for right in range(34)]
        for left in range(34)
    ]
    target_reordered = [
        [
            str(_normalized_ratio(target_gram, mapping[left], mapping[right]))
            for right in range(34)
        ]
        for left in range(34)
    ]
    return {
        "target_id": "A15-C2-H10-34",
        "target_archive_case_id": "A15_9_1_C01",
        "target_archive_path": str(target_path.relative_to(ROOT)),
        "target_archive_sha256": hashlib.sha256(target_path.read_bytes()).hexdigest(),
        "source_to_target_facet_permutation_one_based": [
            mapping[index] + 1 for index in range(34)
        ],
        "normalized_squared_gram_sha256": _digest(source_ratios),
        "target_reordered_normalized_squared_gram_sha256": _digest(target_reordered),
        "all_1156_normalized_squared_entries_match": source_ratios == target_reordered,
        "same_abstract_hyperbolic_Coxeter_polytope": True,
        "fixed_lattices_are_not_integrally_isometric": {
            "borcherds_V4_determinant": 1024,
            "borcherds_V4_smith": [1] * 4 + [2] * 4 + [4] * 3,
            "A15_C2_determinant": target_case["saturated_fixed_lattice"]["determinant"],
            "A15_C2_smith": target_case["saturated_fixed_lattice"][
                "smith_diagonal_absolute"
            ],
        },
    }


def build_certificate(
    source_path: Path = DEFAULT_SOURCE,
    target_path: Path = TARGET_CENSUS,
    gap_executable: str = "gap",
    compiler: str | None = None,
) -> dict[str, Any]:
    raw = source_path.read_bytes()
    source_hash = hashlib.sha256(raw).hexdigest()
    if source_hash != SOURCE_SHA256:
        raise ArithmeticError("the pinned Borcherds source hash changed")
    wall_gram = _extract_matrix(raw.decode("utf-8"), (22, 2, 0))
    if not (
        wall_gram.shape == (210, 210)
        and wall_gram.rank() == 22
        and exact_inertia(wall_gram) == (1, 21, 188)
    ):
        raise ArithmeticError("the ambient Borcherds Gram changed")

    canonical = _canonical_generators()
    full_generators = {
        name: _lift_levi_permutation(wall_gram, permutation)
        for name, permutation in canonical.items()
    }
    elements = _subgroup_elements()
    for permutation in elements:
        if any(
            wall_gram[left, right]
            != wall_gram[permutation[left], permutation[right]]
            for left in range(210)
            for right in range(210)
        ):
            raise ArithmeticError("a pinned V4 element does not preserve the wall Gram")
    subgroup_proof = _gap_subgroup_proof(
        full_generators, elements, gap_executable
    )

    basis_gram, action_2a = _ambient_action(wall_gram, elements[1])
    _, action_2c = _ambient_action(wall_gram, elements[2])
    _, action_product = _ambient_action(wall_gram, elements[3])
    if not (
        action_2a * action_2c == action_product
        and [int(action.trace()) for action in (action_2a, action_2c, action_product)]
        == [6, 8, 8]
    ):
        raise ArithmeticError("the pinned V4 NS fusion changed")
    fixed_basis, fixed_gram, fixed_smith = _common_fixed_lattice(
        basis_gram, (action_2a, action_2c)
    )

    orbits = _wall_orbits(elements)
    if Counter(map(len, orbits)) != Counter({1: 12, 2: 27, 4: 36}):
        raise ArithmeticError("the V4 wall-orbit census changed")
    sums = _orbit_sum_matrix(210, orbits)
    orbit_gram = sums.T * wall_gram * sums
    if not (
        orbit_gram.rank() == 11
        and exact_inertia(orbit_gram) == (1, 10, 64)
    ):
        raise ArithmeticError("the V4 wall-orbit Gram changed")
    facets = [index for index in range(len(orbits)) if orbit_gram[index, index] < 0]
    isotropic = [index for index in range(len(orbits)) if orbit_gram[index, index] == 0]
    timelike = [index for index in range(len(orbits)) if orbit_gram[index, index] > 0]
    if (len(facets), len(isotropic), len(timelike)) != (34, 9, 32):
        raise ArithmeticError("the V4 causal orbit census changed")

    interior = sp.ones(len(orbits), 1)
    interior_pairings = orbit_gram * interior
    interior_square = int((interior.T * orbit_gram * interior)[0])
    if interior_square != 137214 or min(interior_pairings) != 99:
        raise ArithmeticError("the strict invariant interior changed")
    witnesses = []
    for facet in facets:
        witness = interior.copy()
        witness[facet] -= sp.Rational(
            interior_pairings[facet], orbit_gram[facet, facet]
        )
        pairings = orbit_gram * witness
        square = sp.Rational((witness.T * orbit_gram * witness)[0])
        others = [value for index, value in enumerate(pairings) if index != facet]
        if not (square > 0 and pairings[facet] == 0 and min(others) > 0):
            raise ArithmeticError("a V4 negative restriction is not a genuine facet")
        witnesses.append(
            {
                "orbit_index_one_based": facet + 1,
                "witness_square": str(square),
                "minimum_other_pairing": str(min(others)),
            }
        )
    facet_gram = orbit_gram.extract(facets, facets)
    coxeter = _coxeter_audit(
        facet_gram, range(len(facets)), [orbits[index] for index in facets]
    )
    if not (
        coxeter["all_pairs_are_exact_coxeter_relations"]
        and coxeter["edge_counts"] == EXPECTED_COXETER
    ):
        raise ArithmeticError("the V4 Coxeter gate changed")
    facet_records = _facet_lattice_records(
        wall_gram,
        basis_gram,
        fixed_basis,
        fixed_gram,
        orbits,
        facets,
    )
    facet_types = Counter(
        (
            record["orbit_size"],
            record["primitive_square"],
            record["primitive_divisibility"],
        )
        for record in facet_records
    )
    if facet_types != EXPECTED_FACET_TYPES or not all(
        record["reflection_integral"] for record in facet_records
    ):
        raise ArithmeticError("the V4 integral facet-root census changed")

    construction_case = {
        "id": "SS2-V4-H10-34",
        "decision": {"hyperbolic_dimension": 10},
        "restricted_chamber": {"facet_gram": _matrix_payload(facet_gram)},
    }
    selected_compiler = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected_compiler is None:
        raise RuntimeError("a C++17 compiler is required for the Vinberg audit")
    with tempfile.TemporaryDirectory(prefix="borcherds-char2-v4-") as raw_temp:
        executable = Path(raw_temp) / "vinberg"
        subprocess.run(
            [
                selected_compiler,
                "-O3",
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(CPP_SOURCE),
                "-o",
                str(executable),
            ],
            check=True,
        )
        vinberg = _audit_case(executable, construction_case, EXPECTED_VINBERG)

    deduplication = _deduplicate_to_a15(facet_gram, target_path)
    checks = {
        "pinned_ambient_source": True,
        "complete_high_subgroup_filter": True,
        "unique_compatible_V4_conjugacy_class": True,
        "pinned_V4_preserves_all_210_walls": True,
        "saturated_common_fixed_lattice": True,
        "all_34_negative_restrictions_are_genuine_facets": True,
        "all_pairs_are_exact_Coxeter_relations": True,
        "all_34_primitive_facet_reflections_are_integral": True,
        "literal_Vinberg_finite_volume": vinberg["finite_volume"],
        "exact_A15_normalized_Gram_isomorphism": deduplication[
            "all_1156_normalized_squared_entries_match"
        ],
    }
    facet_payload = _matrix_payload(facet_gram)
    return {
        "schema": SCHEMA,
        "id": "SS2-V4-H10-34",
        "source": {
            "path": str(source_path),
            "sha256": source_hash,
            "matrix_key": "(22,2,0)",
            "cpp_vinberg_source": str(CPP_SOURCE.relative_to(ROOT)),
            "cpp_vinberg_source_sha256": hashlib.sha256(CPP_SOURCE.read_bytes()).hexdigest(),
        },
        "high_subgroup_filter": {
            "allowed_nonidentity_element_classes": ["2a", "2c", "2d", "3b"],
            "NS_traces": {"identity": 22, "2a": 6, "2c": 8, "2d": 0, "3b": 7},
            "reynolds_inequality": "5*n_2a + 3*n_2c + 11*n_2d + 4*n_3b <= 11",
            "noncyclic_conclusion": "only V4 with fusion 2a + 2c + 2c can have fixed rank at least 11",
            "proof_note": (
                "a noncyclic exponent-two group is elementary abelian; rank at least "
                "three violates the inequality, and integrality excludes 2c+2c+2c. "
                "An order-three pair already costs eight, while adjoining an involution "
                "forces group order divisible by six and more elements than the bound allows."
            ),
        },
        "subgroup_class_proof": subgroup_proof,
        "action": {
            "group": "V4",
            "fusion": {"identity": 1, "2a": 1, "2c": 2},
            "order": 4,
            "generator_permutations_one_based": [
                [value + 1 for value in elements[1]],
                [value + 1 for value in elements[2]],
            ],
            "NS_action_traces": [6, 8, 8],
            "geometric_status": (
                "not a K3 automorphism subgroup: 2a is geometric, while both 2c "
                "elements lie outside Dolgachev--Kondo's index-two subgroup"
            ),
        },
        "fixed_lattice": {
            "definition": "NS(X) intersect Fix_Q(V4)",
            "basis_columns_in_ambient_basis": _matrix_payload(fixed_basis),
            "gram": _matrix_payload(fixed_gram),
            "rank": 11,
            "signature_positive_negative": [1, 10],
            "determinant": 1024,
            "smith_diagonal_absolute": fixed_smith,
        },
        "restricted_chamber": {
            "wall_orbit_count": len(orbits),
            "wall_orbit_length_histogram": {"1": 12, "2": 27, "4": 36},
            "negative_facet_count": 34,
            "isotropic_inactive_orbit_count": 9,
            "timelike_inactive_orbit_count": 32,
            "strict_interior_square": interior_square,
            "strict_interior_minimum_pairing": int(min(interior_pairings)),
            "facet_orbit_indices_one_based": [index + 1 for index in facets],
            "facet_witnesses": witnesses,
            "facet_gram": facet_payload,
            "facet_gram_sha256": _digest(facet_payload),
            "coxeter_edge_counts": coxeter["edge_counts"],
            "primitive_facet_lattice_records": facet_records,
            "primitive_facet_type_histogram": {
                f"orbit{orbit}_square{square}_div{divisibility}": count
                for (orbit, square, divisibility), count in sorted(facet_types.items())
            },
        },
        "vinberg": vinberg,
        "exact_deduplication": deduplication,
        "decision": {
            "hyperbolic_dimension": 10,
            "facet_count": 34,
            "finite_volume": True,
            "compact": False,
            "classification": "finite_volume_noncompact_Coxeter_polytope",
            "abstract_novelty": "not new: exactly A15-C2-H10-34 after facet rescaling",
            "realization_novelty": (
                "new divisor-lattice/diagram V4 realization in the characteristic-two "
                "Borcherds stock; not induced by a K3 automorphism subgroup"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", nargs="?", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--target-census", type=Path, default=TARGET_CENSUS)
    parser.add_argument("--gap", default="gap")
    parser.add_argument("--compiler")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(
        args.source, args.target_census, args.gap, args.compiler
    )
    if args.check_archive:
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt V4 certificate differs from archive")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "id": payload["id"],
                "finite_vertices": payload["vinberg"]["finite_vertex_count"],
                "cusps": payload["vinberg"]["ideal_vertex_count"],
                "deduplicates_to": payload["exact_deduplication"]["target_id"],
                "all_checks_pass": payload["all_checks_pass"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
