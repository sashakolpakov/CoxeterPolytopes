"""Exact subgroup-fold census for the 316-wall Jacobian--Kummer chamber.

The input is Ichiro Shimada's GAP-format computational-data file accompanying
Dolgachev--Shimada.  In particular, this program uses the full published
lattice/chamber stabilizer ``O(S16,D16)`` of order 23040.  It does not require
a subgroup to be induced by automorphisms of the K3 surface.

GAP enumerates all conjugacy classes of subgroups of the finite chamber
stabilizer.  The common fixed rank is calculated in the published integral
17-dimensional lattice representation, so classes of fixed rank below 11 can
be discarded rigorously: adding generators can only shrink the fixed space.
The complete high-rank list is frozen here.  The unique H12 class points to
its dedicated exact metric, Vinberg, and cone certificates; the two H10
classes remain explicitly marked for later screening.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
import re
import subprocess
from typing import Any

import sympy as sp

from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    _extract_assignment,
    _matrix_payload,
    _rational_text,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "data"
    / "external"
    / "dolgachev_shimada_15_nodal_2019"
    / "15nodalcompdata.txt"
)
STOCK_OUTPUT = (
    ROOT / "data" / "known_examples" / "jacobian_kummer_d16_chamber_stock.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_d16_subgroup_folds.json"
)


def _parse_source(path: Path = SOURCE) -> dict[str, Any]:
    """Parse exactly the S16/D16 fields from Shimada's GAP data file."""

    text = path.read_text(encoding="utf-8").replace("\\\n", "")
    gram = sp.Matrix(ast.literal_eval(_extract_assignment(text, "GramS16")))
    ample = sp.Matrix(1, 17, ast.literal_eval(_extract_assignment(text, "ample16")))
    stabilizer_generators = tuple(
        sp.Matrix(matrix)
        for matrix in ast.literal_eval(
            _extract_assignment(text, "GeneratorsOS16D16")
        )
    )
    geometric_generators = tuple(
        sp.Matrix(matrix)
        for matrix in ast.literal_eval(
            _extract_assignment(text, "GeneratorsAutY16D16")
        )
    )

    walls_text = _extract_assignment(text, "D16WallRecs")
    vector_texts = re.findall(r"\bvector\s*:=\s*\[([^\]]+)\]", walls_text)
    wall_vectors = tuple(
        tuple(sp.Rational(entry.strip()) for entry in vector.split(","))
        for vector in vector_texts
    )
    wall_numbers = tuple(
        int(value) for value in re.findall(r"\bno\s*:=\s*(\d+)", walls_text)
    )
    wall_orbit_numbers = tuple(
        int(value) for value in re.findall(r"\borbit\s*:=\s*(\d+)", walls_text)
    )
    wall_types = tuple(
        re.findall(r'\binnout\s*:=\s*"([^"]+)"', walls_text)
    )
    wall_squares = tuple(
        sp.Rational(value.strip())
        for value in re.findall(
            r"\bn\s*:=\s*([^,]+),\s*no\s*:=", walls_text
        )
    )
    if not (
        len(wall_vectors)
        == len(wall_numbers)
        == len(wall_orbit_numbers)
        == len(wall_types)
        == len(wall_squares)
        == 316
    ):
        raise ValueError("D16WallRecs does not contain 316 complete wall records")
    if wall_numbers != tuple(range(1, 317)):
        raise ValueError("D16WallRecs wall numbers are not 1,...,316")
    if any(len(vector) != 17 for vector in wall_vectors):
        raise ValueError("a D16 wall vector does not have 17 coordinates")

    return {
        "gram": gram,
        "ample": ample,
        "stabilizer_generators": stabilizer_generators,
        "geometric_generators": geometric_generators,
        "wall_vectors": wall_vectors,
        "wall_numbers": wall_numbers,
        "wall_orbit_numbers": wall_orbit_numbers,
        "wall_types": wall_types,
        "published_wall_squares": wall_squares,
    }


def _wall_permutation(
    matrix: sp.Matrix,
    gram: sp.Matrix,
    ample: sp.Matrix,
    wall_vectors: tuple[tuple[sp.Rational, ...], ...],
) -> Permutation:
    """Derive and verify the wall permutation of one lattice matrix."""

    if matrix * gram * matrix.T != gram:
        raise ArithmeticError("a published chamber generator is not an isometry")
    if ample * matrix != ample:
        raise ArithmeticError("a published chamber generator does not fix ample16")
    index = {vector: position for position, vector in enumerate(wall_vectors)}
    images = tuple(
        tuple(sp.Matrix(1, gram.rows, vector) * matrix)
        for vector in wall_vectors
    )
    if any(image not in index for image in images):
        raise ArithmeticError("a published chamber generator does not permute D16 walls")
    permutation = tuple(index[image] for image in images)
    if len(set(permutation)) != len(wall_vectors):
        raise ArithmeticError("the derived wall action is not a permutation")
    return permutation


def _gap_permutation(permutation: Permutation) -> str:
    return "PermList([" + ",".join(str(image + 1) for image in permutation) + "])"


def _gap_matrix(matrix: sp.Matrix) -> str:
    return str([[int(matrix[row, column]) for column in range(matrix.cols)] for row in range(matrix.rows)]).replace(" ", "")


def _gap_census_script(
    generator_permutations: tuple[Permutation, ...],
    generator_matrices: tuple[sp.Matrix, ...],
    geometric_permutations: tuple[Permutation, ...],
) -> str:
    """Return a deterministic GAP program for the complete subgroup census."""

    permutations = ",".join(_gap_permutation(p) for p in generator_permutations)
    matrices = ",".join(_gap_matrix(m) for m in generator_matrices)
    geometric = ",".join(_gap_permutation(p) for p in geometric_permutations)
    # Keeping each statement on its own line makes GAP diagnostics usable while
    # the machine-readable records themselves are tab-separated.
    return f"""
pgens := [{permutations}];;
G := Group(pgens);;
mgens := [{matrices}];;
MG := Group(mgens);;
hom := GroupHomomorphismByImages(G,MG,pgens,mgens);;
iso := SmallerDegreePermutationRepresentation(G);;
K := Image(iso);;
Geom := Group([{geometric}]);;
C := ConjugacyClassesSubgroups(K);;
Print("HEADER\\t",GAPInfo.Version,"\\t",Size(G),"\\t",StructureDescription(G),
      "\\t",NrMovedPoints(K),"\\t",Size(Image(hom)),"\\t",Size(Geom),"\\n");
for pair in Collected(List(Elements(K),Order)) do
  Print("ELEMENTORDER\\t",pair[1],"\\t",pair[2],"\\n");
od;;
i := 0;;
for c in C do
  i := i+1;;
  H := Representative(c);;
  kg := GeneratorsOfGroup(H);;
  if Length(kg)=0 then
    fixed := 17;; mm := [];;
  else
    og := List(kg,x->PreImagesRepresentative(iso,x));;
    mm := List(og,x->Image(hom,x));;
    fixed := 17-RankMat(Concatenation(List(mm,x->x-IdentityMat(17))));;
  fi;;
  Print("CLASS\\t",i,"\\t",Size(H),"\\t",Size(c),"\\t",fixed,"\\n");
  if fixed>=11 then
    Print("HIGHBEGIN\\t",i,"\\n",mm,"\\nHIGHEND\\t",i,"\\n");
  fi;;
od;;
QUIT;
"""


def _parse_gap_census(output: str) -> dict[str, Any]:
    header_match = re.search(
        r"^HEADER\t([^\t]+)\t(\d+)\t([^\t]+)\t(\d+)\t(\d+)\t(\d+)$",
        output,
        flags=re.MULTILINE,
    )
    if header_match is None:
        raise ArithmeticError("GAP subgroup census did not emit a valid header")
    version, order, structure, reduced_degree, image_order, geometric_order = (
        header_match.groups()
    )
    element_orders = [
        (int(order), int(count))
        for order, count in re.findall(
            r"^ELEMENTORDER\t(\d+)\t(\d+)$", output, re.MULTILINE
        )
    ]
    if not element_orders:
        raise ArithmeticError("GAP subgroup census omitted element orders")

    classes: list[dict[str, Any]] = []
    for match in re.finditer(
        r"^CLASS\t(\d+)\t(\d+)\t(\d+)\t(\d+)$",
        output,
        flags=re.MULTILINE,
    ):
        index, subgroup_order, class_size, fixed_rank = (
            match.groups()
        )
        classes.append(
            {
                "id": f"class_{int(index):04d}",
                "gap_class_index": int(index),
                "subgroup_order": int(subgroup_order),
                "conjugacy_class_size": int(class_size),
                "normalizer_order": 23040 // int(class_size),
                "fixed_rank": int(fixed_rank),
                "hyperbolic_dimension": int(fixed_rank) - 1,
            }
        )

    high_matrices: dict[int, tuple[sp.Matrix, ...]] = {}
    for match in re.finditer(
        r"^HIGHBEGIN\t(\d+)\n(.*?)\nHIGHEND\t\1$",
        output,
        flags=re.MULTILINE | re.DOTALL,
    ):
        index = int(match.group(1))
        matrices = ast.literal_eval(match.group(2))
        high_matrices[index] = tuple(sp.Matrix(matrix) for matrix in matrices)

    return {
        "gap_version": version,
        "group_order": int(order),
        "group_structure": structure,
        "reduced_faithful_permutation_degree": int(reduced_degree),
        "matrix_image_order": int(image_order),
        "geometric_subgroup_order": int(geometric_order),
        "element_order_histogram": {
            str(element_order): count for element_order, count in element_orders
        },
        "classes": classes,
        "high_matrices": high_matrices,
    }


def _run_gap_census(source: dict[str, Any]) -> dict[str, Any]:
    gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    walls = source["wall_vectors"]
    generator_permutations = tuple(
        _wall_permutation(matrix, gram, ample, walls)
        for matrix in source["stabilizer_generators"]
    )
    geometric_permutations = tuple(
        _wall_permutation(matrix, gram, ample, walls)
        for matrix in source["geometric_generators"]
    )
    script = _gap_census_script(
        generator_permutations,
        source["stabilizer_generators"],
        geometric_permutations,
    )
    completed = subprocess.run(
        ["gap", "-q", "-o", "2g"],
        input=script,
        text=True,
        capture_output=True,
        timeout=900,
        check=False,
    )
    if completed.returncode != 0 or completed.stderr.strip():
        raise RuntimeError(
            "GAP subgroup census failed:\n"
            + completed.stdout[-4000:]
            + "\n"
            + completed.stderr[-4000:]
        )
    census = _parse_gap_census(completed.stdout)
    census["generator_wall_permutations"] = generator_permutations
    census["geometric_generator_wall_permutations"] = geometric_permutations
    return census


@lru_cache(maxsize=2)
def build_artifacts(path: Path = SOURCE) -> tuple[dict[str, Any], dict[str, Any]]:
    source = _parse_source(path)
    gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    wall_rows = tuple(
        sp.Matrix(1, gram.rows, vector) for vector in source["wall_vectors"]
    )

    calculated_squares = tuple(
        sp.Rational((wall * gram * wall.T)[0]) for wall in wall_rows
    )
    if calculated_squares != source["published_wall_squares"]:
        raise ArithmeticError("published and calculated D16 wall squares disagree")
    if exact_inertia(gram) != (1, 16, 0):
        raise ArithmeticError("GramS16 is not Lorentzian of signature (1,16)")

    gap = _run_gap_census(source)
    if not (
        gap["group_order"] == gap["matrix_image_order"] == 23040
        and gap["geometric_subgroup_order"] == 32
        and gap["group_structure"] == "(C2 x C2 x C2 x C2 x C2) : S6"
        and len(gap["classes"]) == 2455
    ):
        raise ArithmeticError("unexpected GAP chamber-group census")

    high_classes = [record for record in gap["classes"] if record["fixed_rank"] >= 11]
    if set(gap["high_matrices"]) != {
        record["gap_class_index"] for record in high_classes
    }:
        raise ArithmeticError("GAP omitted a high-fixed-rank subgroup representative")
    high_candidates: list[dict[str, Any]] = []
    for record in high_classes:
        candidate = dict(record)
        candidate["structure"] = "1" if record["subgroup_order"] == 1 else "C2"
        candidate["representative_generator_matrices"] = [
            _matrix_payload(matrix)
            for matrix in gap["high_matrices"][record["gap_class_index"]]
        ]
        if record["subgroup_order"] == 1:
            candidate["status"] = "unfolded_parent_chamber"
        elif record["fixed_rank"] == 13:
            candidate.update(
                {
                    "status": "exact_finite_volume_noncompact_coxeter_polytope",
                    "facet_count": 172,
                    "metric_certificate": (
                        "data/enumerations/"
                        "jacobian_kummer_h12_interim_facet_gram.json"
                    ),
                    "vinberg_certificate": (
                        "data/enumerations/"
                        "jacobian_kummer_h12_vinberg_volume_certificate.json"
                    ),
                    "exact_cone_certificate": (
                        "data/enumerations/"
                        "jacobian_kummer_h12_exact_cone_certificate.json"
                    ),
                }
            )
        else:
            candidate["status"] = "not_screened_lower_priority_h10_candidate"
        high_candidates.append(candidate)

    stock = {
        "schema": "jacobian_kummer_d16_exact_chamber_stock_v1",
        "source": {
            "title": "15-nodal quartic surfaces. Part II: The automorphism group: Computational data",
            "author": "Ichiro Shimada",
            "paper_authors": ["Igor Dolgachev", "Ichiro Shimada"],
            "paper_doi": "10.1007/s12215-019-00464-7",
            "paper_arxiv": "1908.05390",
            "computational_data_url": "https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/15NodalQuart/15nodal/compdata/15nodalcompdata.txt",
            "computational_data_guide_url": "https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/15NodalQuart/15nodal/compdatatex/15nodalCompData.pdf",
            "source_file": path.name,
            "source_assignments": [
                "GramS16",
                "ample16",
                "GeneratorsOS16D16",
                "GeneratorsAutY16D16",
                "D16WallRecs",
            ],
        },
        "algebro_geometric_object": {
            "surface": "minimal resolution Y16 of a general Jacobian Kummer quartic",
            "construction": "Kum(Jac(C)) for a general genus-2 curve C",
            "divisor_lattice": "S16 = NS(Y16)",
            "induced_chamber": "D16 from Kondo's primitive embedding S16 into II_(1,25)",
            "surface_automorphism_required_for_fold_search": False,
        },
        "lattice": {
            "rank": gram.rows,
            "determinant": int(gram.det()),
            "inertia": list(exact_inertia(gram)),
            "gram": _matrix_payload(gram),
            "ample16": [_rational_text(value) for value in ample],
            "ample_square": _rational_text(sp.Rational((ample * gram * ample.T)[0])),
        },
        "walls": [
            {
                "index": index + 1,
                "vector": [_rational_text(value) for value in vector],
                "square": _rational_text(calculated_squares[index]),
                "kind": source["wall_types"][index],
                "published_orbit_number": source["wall_orbit_numbers"][index],
            }
            for index, vector in enumerate(source["wall_vectors"])
        ],
        "wall_summary": {
            "count": len(wall_rows),
            "full_stabilizer_orbit_size_histogram": {
                str(orbit): count
                for orbit, count in sorted(Counter(source["wall_orbit_numbers"]).items())
            },
            "kind_histogram": dict(sorted(Counter(source["wall_types"]).items())),
            "square_histogram": {
                _rational_text(square): count
                for square, count in sorted(Counter(calculated_squares).items())
            },
        },
        "full_lattice_chamber_stabilizer": {
            "published_name": "O(S16,D16)",
            "order": gap["group_order"],
            "structure": gap["group_structure"],
            "element_order_histogram": gap["element_order_histogram"],
            "published_generator_count": len(source["stabilizer_generators"]),
            "generator_matrices": [
                _matrix_payload(matrix) for matrix in source["stabilizer_generators"]
            ],
            "generator_wall_permutations": [
                [image + 1 for image in permutation]
                for permutation in gap["generator_wall_permutations"]
            ],
            "faithful_reduced_permutation_degree": gap[
                "reduced_faithful_permutation_degree"
            ],
            "checks": {
                "every_generator_preserves_gram": True,
                "every_generator_fixes_ample16": True,
                "every_generator_permutes_all_316_walls": True,
                "wall_action_order_equals_matrix_image_order": True,
            },
        },
        "projective_surface_automorphism_subgroup": {
            "published_name": "Aut(Y16,alpha16)",
            "order": gap["geometric_subgroup_order"],
            "structure": "(C2)^5",
            "index_in_full_lattice_chamber_stabilizer": gap["group_order"]
            // gap["geometric_subgroup_order"],
            "used_as_search_gate": False,
        },
    }

    fixed_rank_histogram = Counter(record["fixed_rank"] for record in gap["classes"])
    subgroup_order_histogram = Counter(record["subgroup_order"] for record in gap["classes"])
    nontrivial_high = [
        candidate
        for candidate in high_candidates
        if candidate["subgroup_order"] > 1
    ]
    enumeration = {
        "schema": "jacobian_kummer_d16_subgroup_fold_census_v1",
        "method": {
            "symmetry_scope": "all conjugacy classes of subgroups of the full published O(S16,D16)",
            "surface_automorphism_required": False,
            "high_dimension_cutoff": "fixed rank >= 11, equivalently hyperbolic dimension > 9",
            "subgroup_enumerator": "GAP ConjugacyClassesSubgroups on an exact faithful degree-32 permutation image",
            "fixed_space_gate": "exact common kernel in the published integral rank-17 representation",
            "normal_restriction": "sum wall normals over every subgroup orbit",
            "wall_reduction": "exact positive rational dependences or strict timelike facet witnesses",
            "coxeter_gate": "exact normalized squared products",
            "high_candidate_policy": (
                "the unique H12 class is certified in dedicated exact artifacts; "
                "the two H10 classes are frozen for later screening"
            ),
        },
        "gap_version": gap["gap_version"],
        "group_order": gap["group_order"],
        "subgroup_conjugacy_class_count": len(gap["classes"]),
        "actual_subgroup_count": sum(
            record["conjugacy_class_size"] for record in gap["classes"]
        ),
        "subgroup_class_order_histogram": {
            str(order): count for order, count in sorted(subgroup_order_histogram.items())
        },
        "fixed_rank_class_histogram": {
            str(rank): count for rank, count in sorted(fixed_rank_histogram.items(), reverse=True)
        },
        "subgroup_classes": gap["classes"],
        "high_dimension_candidates": high_candidates,
        "summary": {
            "high_dimension_class_count_including_trivial": len(high_candidates),
            "nontrivial_high_dimension_class_count": len(nontrivial_high),
            "nontrivial_high_dimension_specs": sorted(
                [
                    candidate["subgroup_order"],
                    candidate["fixed_rank"],
                    candidate["hyperbolic_dimension"],
                ]
                for candidate in nontrivial_high
            ),
            "certified_h12_finite_volume_noncompact_count": 1,
            "pending_h10_class_count": 2,
            "nontrivial_high_dimension_outcomes": [
                {
                    "id": candidate["id"],
                    "dimension": candidate["hyperbolic_dimension"],
                    "status": candidate["status"],
                }
                for candidate in nontrivial_high
            ],
        },
    }
    return stock, enumeration


def _write_or_check(path: Path, payload: dict[str, Any], check: bool) -> None:
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if check:
        if not path.is_file() or path.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rendered, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--stock-output", type=Path, default=STOCK_OUTPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    stock, enumeration = build_artifacts(arguments.source)
    _write_or_check(arguments.stock_output, stock, arguments.check_archive)
    _write_or_check(arguments.output, enumeration, arguments.check_archive)
    print(json.dumps(enumeration["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
