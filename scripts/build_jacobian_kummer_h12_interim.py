"""Freeze the exact 172-facet H12 Jacobian--Kummer fold.

This is deliberately an interim *metric* certificate.  It proves the exact
fixed slice, irredundant facet set, and Coxeter pair gate.  Finite volume is a
separate Vinberg edge/endpoint and full parabolic-neighbourhood calculation.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import json
from pathlib import Path
import re
import subprocess
from typing import Any

import sympy as sp

from scripts.search_apery_fermi_chamber_folds import _exact_nonnegative_dependence
from scripts.search_jacobian_kummer_chamber_folds import (
    ROOT,
    SOURCE,
    _gap_matrix,
    _gap_permutation,
    _matrix_payload,
    _parse_source,
    _rational_text,
    _wall_permutation,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _coxeter_label


OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_interim_facet_gram.json"
)


def _h12_involution(source: dict[str, Any]) -> sp.Matrix:
    """Return a representative of the unique fixed-rank-13 involution class."""

    gram = source["gram"]
    ample = source["ample"]
    wall_vectors = source["wall_vectors"]
    generator_permutations = tuple(
        _wall_permutation(matrix, gram, ample, wall_vectors)
        for matrix in source["stabilizer_generators"]
    )
    script = (
        "pgens:=["
        + ",".join(_gap_permutation(p) for p in generator_permutations)
        + "];;G:=Group(pgens);;mgens:=["
        + ",".join(_gap_matrix(m) for m in source["stabilizer_generators"])
        + "];;MG:=Group(mgens);;"
        + "hom:=GroupHomomorphismByImages(G,MG,pgens,mgens);;found:=0;;"
        + "for c in ConjugacyClasses(G) do r:=Representative(c);;"
        + "if Order(r)=2 then m:=Image(hom,r);;if TraceMat(m)=9 then "
        + "found:=found+1;;Print(\"BEGIN\\n\",m,\"\\nEND\\n\");;fi;;fi;;od;;"
        + "Print(\"COUNT \",found,\"\\n\");;QUIT;\n"
    )
    completed = subprocess.run(
        ["gap", "-q"],
        input=script,
        text=True,
        capture_output=True,
        timeout=120,
        check=False,
    )
    if completed.returncode or completed.stderr.strip():
        raise RuntimeError(completed.stdout + "\n" + completed.stderr)
    matches = re.findall(r"BEGIN\n(.*?)\nEND", completed.stdout, re.DOTALL)
    count = re.search(r"COUNT (\d+)", completed.stdout)
    if count is None or int(count.group(1)) != 1 or len(matches) != 1:
        raise ArithmeticError("the fixed-rank-13 involution class is not unique")
    matrix = sp.Matrix(ast.literal_eval(matches[0]))
    if matrix.trace() != 9 or matrix**2 != sp.eye(17):
        raise ArithmeticError("the selected H12 representative is not the required involution")
    return matrix


def build_interim(path: Path = SOURCE) -> dict[str, Any]:
    source = _parse_source(path)
    gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    walls = tuple(
        sp.Matrix(1, gram.rows, vector) for vector in source["wall_vectors"]
    )
    involution = _h12_involution(source)
    permutation = _wall_permutation(
        involution, gram, ample, source["wall_vectors"]
    )

    unseen = set(range(len(walls)))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        orbit = tuple(sorted({start, permutation[start]}))
        unseen.difference_update(orbit)
        orbits.append(orbit)
    orbit_normals = tuple(
        sum((walls[index] for index in orbit), sp.zeros(1, gram.rows))
        for orbit in orbits
    )
    stacked = sp.Matrix.vstack(*orbit_normals)
    orbit_gram = stacked * gram * stacked.T
    _, pivot_columns = stacked.rref()
    coordinates = stacked[:, list(pivot_columns)]
    if stacked.rank() != 13 or exact_inertia(orbit_gram) != (1, 12, 167):
        raise ArithmeticError("the H12 orbit-normal span has the wrong exact inertia")

    facets: list[int] = []
    facet_witnesses: list[dict[str, Any]] = []
    unresolved: list[int] = []
    for index, normal in enumerate(orbit_normals):
        norm = sp.Rational(orbit_gram[index, index])
        if norm >= 0:
            unresolved.append(index)
            continue
        ample_pairing = sp.Rational((ample * gram * normal.T)[0])
        witness = ample - ample_pairing / norm * normal
        witness_square = sp.Rational((witness * gram * witness.T)[0])
        pairings = tuple(
            sp.Rational(value) for value in list(witness * gram * stacked.T)
        )
        if (
            witness_square > 0
            and pairings[index] == 0
            and all(value > 0 for other, value in enumerate(pairings) if other != index)
        ):
            facets.append(index)
            facet_witnesses.append(
                {
                    "orbit_index": index + 1,
                    "witness_square": _rational_text(witness_square),
                    "minimum_strict_pairing": _rational_text(
                        min(value for other, value in enumerate(pairings) if other != index)
                    ),
                }
            )
        else:
            unresolved.append(index)

    redundancy_certificates: list[dict[str, Any]] = []
    for target in unresolved:
        dependence = _exact_nonnegative_dependence(coordinates, target)
        if dependence is None:
            raise ArithmeticError("an orbit inequality is neither witnessed nor redundant")
        redundancy_certificates.append(
            {
                "orbit_index": target + 1,
                "terms": [
                    {
                        "orbit_index": support + 1,
                        "coefficient": _rational_text(coefficient),
                    }
                    for support, coefficient in dependence
                ],
            }
        )

    facet_gram = orbit_gram.extract(facets, facets)
    labels: Counter[str] = Counter()
    for left_position, left in enumerate(facets):
        for right in facets[left_position + 1 :]:
            labels[
                _coxeter_label(
                    orbit_gram[left, left],
                    orbit_gram[right, right],
                    orbit_gram[left, right],
                )
            ] += 1
    if len(facets) != 172 or len(unresolved) != 8:
        raise ArithmeticError("unexpected irredundant H12 facet count")
    if labels != Counter(
        {"2": 2820, "3": 144, "4": 800, "infinity": 1750, "ultraparallel": 9192}
    ):
        raise ArithmeticError("unexpected exact H12 Coxeter pair histogram")

    return {
        "schema": "jacobian_kummer_h12_interim_facet_gram_v1",
        "status": {
            "exact_fixed_slice": True,
            "exact_irredundant_facets": True,
            "exact_coxeter_pair_gate": True,
            "vinberg_finite_volume_audit_complete": False,
            "finite_volume_claimed": False,
        },
        "source": {
            "surface": "general Jacobian Kummer K3 surface Y16",
            "paper_doi": "10.1007/s12215-019-00464-7",
            "paper_arxiv": "1908.05390",
            "computational_data_assignment": "D16WallRecs",
            "source_file": str(path.relative_to(ROOT)),
        },
        "group_class": {
            "ambient_group": "O(S16,D16) = (C2)^5 : S6",
            "ambient_group_order": 23040,
            "gap_subgroup_class_index": 5,
            "subgroup_structure": "C2",
            "conjugacy_class_size": 30,
            "normalizer_order": 768,
            "representative_selection": "unique involution conjugacy class with lattice trace 9",
            "representative_matrix": _matrix_payload(involution),
            "fixed_rank": 13,
            "hyperbolic_dimension": 12,
        },
        "wall_orbits": [[index + 1 for index in orbit] for orbit in orbits],
        "wall_orbit_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(len, orbits)).items())
        },
        "orbit_normal_rank": stacked.rank(),
        "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
        "facet_count": len(facets),
        "facet_orbit_indices": [index + 1 for index in facets],
        "facet_source_wall_orbits": [
            [source_index + 1 for source_index in orbits[index]] for index in facets
        ],
        "redundant_orbit_count": len(unresolved),
        "redundancy_certificates": redundancy_certificates,
        "facet_witnesses": facet_witnesses,
        "facet_normal_square_histogram": {
            _rational_text(square): count
            for square, count in sorted(
                Counter(sp.Rational(facet_gram[index, index]) for index in range(len(facets))).items()
            )
        },
        "coxeter_edge_counts": dict(sorted(labels.items())),
        "facet_gram_sign_convention": "negative-square orbit-sum wall normals",
        "facet_gram": _matrix_payload(facet_gram),
        "positive_wall_gram": _matrix_payload(-facet_gram),
        "checks": {
            "representative_is_involution": True,
            "representative_preserves_gram_and_ample_class": True,
            "representative_permutes_all_316_published_walls": True,
            "orbit_normals_span_fixed_lorentzian_space": True,
            "every_retained_wall_has_a_strict_timelike_relative_interior_witness": True,
            "every_removed_wall_has_an_exact_nonnegative_dependence": True,
            "all_14706_facet_pairs_pass_the_exact_coxeter_gate": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build_interim(arguments.source)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text() != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "facet_count": payload["facet_count"],
                "redundant_orbit_count": payload["redundant_orbit_count"],
                "coxeter_edge_counts": payload["coxeter_edge_counts"],
                "finite_volume_claimed": payload["status"]["finite_volume_claimed"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
