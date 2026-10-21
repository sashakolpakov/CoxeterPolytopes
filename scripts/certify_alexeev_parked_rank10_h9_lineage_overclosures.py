#!/usr/bin/env python3
"""Certify arithmetic overclosures and bounded lineage of five parked H9 seeds.

The immutable five-seed construction archive proves the hyperbolic geometry.
This separate certificate does four different jobs:

* exhaust every integral overlattice preserved by all displayed reflections;
* normalize the unique maximal stable lattice to a pinned rank-ten
  Scharlau--Kirschmer class;
* compare the five labelled diagrams with pinned ordinary and index-two
  CoxIter snapshots; and
* record exactly when Allcock's classical redoubling lemma applies.

The certificate deliberately distinguishes an exact lattice/corpus statement
from a literature-scope statement.  It makes no worldwide priority claim.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import tarfile
from typing import Any, Iterable, Sequence
import zipfile

import networkx as nx
import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _canonical_sha256,
    _matrix,
    _overlattice_basis,
    _payload,
    _q_histogram,
    _reflection,
    _smith,
    _stable_subgroups,
)
from scripts.certify_higher_dimensional_rational_lineage import (
    _compare,
    _e8_cartan,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/alexeev_parked_rank10_h9_examples.json"
OUTPUT = (
    ROOT
    / "data/enumerations/alexeev_parked_rank10_h9_lineage_overclosures.json"
)
SCHEMA = "coxeter.alexeev_parked_rank10_h9.lineage_overclosures.v1"
SOURCE_SHA256 = "e181b8c021fe0f37d23b53b33e4e3900559502d0f8c982a4a9e9cbcf0866e430"

COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
COXITER_ROOT = Path("/private/tmp/coxiter-current/CoxIter-master")
DOUBLING_ARCHIVE = COXITER_ROOT / "graphs/testing_index2.tar.gz"
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
DOUBLING_ARCHIVE_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)

SK_ROOT = Path("/private/tmp/Scharlau-Kirschmer-Lorentz")
SK_LIST9 = SK_ROOT / "List9"
SK_GRAM9 = SK_ROOT / "Gram9"
SK_LIST9_SHA256 = "b7f6982a5c03b1ff07a5a9f383517438a926a499d27fb33ea89a8a69d7fb66c3"
SK_GRAM9_SHA256 = "772187d28a067c2d10f2fe00d7cee07a063a6c3d79be5cc8400c2c0ca8c96c36"

ALLCOCK_PDF = Path("/private/tmp/allcock-doubling.pdf")
ALLCOCK_SOURCE = Path("/private/tmp/allcock-source.tar")
ALLCOCK_FIGURE4 = Path("/private/tmp/allcock-source/figs/fig4.pdf")
ALLCOCK_PDF_SHA256 = "df06b0722bd97806b236857e4be65be8ab350434eaee8553415165090a2752d3"
ALLCOCK_SOURCE_SHA256 = "dbf49889ab9a2529aafeff69ee627c7ea4920c9f94d56f093e5db18a52055113"
ALLCOCK_FIGURE4_SHA256 = "02d3952254f07fffbd4a247e700f5461969fefbb87560866e6f037ca7d2dcc5f"

ALEXEEV_ENGEL_SOURCE = Path("/private/tmp/alexeev_engel_2208.10383.tar")
ALEXEEV_ENGEL_SOURCE_SHA256 = (
    "d395cc19a5d32064264b2f498a3bcf8a9f5153df272c37bd67fc28c721570470"
)

IMHOF_MATCH = (
    "CoxIter-master/graphs/doubly truncated orthoscheme/9-imhof85_03.coxiter"
)
IMHOF_HISTOGRAM_COLLISIONS = [
    "CoxIter-master/graphs/doubly truncated orthoscheme/9-imhof85_02.coxiter",
    IMHOF_MATCH,
    "CoxIter-master/graphs/doubly truncated orthoscheme/9-imhof85_04.coxiter",
]
IMHOF_CANDIDATE_TO_REFERENCE_ONE_BASED = [6, 11, 10, 9, 8, 1, 2, 3, 4, 5, 7, 12]

P9_PATH = (
    "CoxIter-master/graphs/group units qf/"
    "9-1,1,1,1,1,1,1,1,1,1.coxiter"
)
TUMARKIN_11_01_PATH = (
    "CoxIter-master/graphs/Tumarkin pyramids n+2/9-tum04_11_01.coxiter"
)

EXPECTED_STABLE = {
    "A14_8_0_C01": ([1, 4], [1, 4], 4, -64, (10, 6, 0)),
    "A14_8_1_C02": ([1, 2, 4, 8], [1, 2, 4, 8], 8, -16, (10, 4, 0)),
    "A14_10_1_C01": ([1, 2, 4], [1, 2, 4], 4, -256, (10, 8, 0)),
    "A14_10_1_C02": ([1, 2, 4, 8], [1, 2, 4, 8], 8, -64, (10, 6, 0)),
    "A16_14_1_C01": ([1, 16, 32], [1, 16, 32], 32, -16, (10, 4, 0)),
}

EXPECTED_DISJOINT_DOUBLING_PAIRS = {
    "A14_8_0_C01": [],
    "A14_8_1_C02": [[15, 16]],
    "A14_10_1_C01": [],
    "A14_10_1_C02": [[6, 15]],
    "A16_14_1_C01": [[1, 14], [8, 15]],
}

SK_TARGETS = {
    7: {
        "genus_symbol": "II_(9,1)(2^2)",
        "determinant": -4,
        "smith_invariants": [1] * 8 + [2] * 2,
        "discriminant_q_mod_2_histogram": {"0": 3, "1": 1},
        "full_chamber_facets": 10,
        "full_chamber_finite_vertices": 8,
        "full_chamber_cusps": 2,
        "classical_chamber": "Vinberg--Kaplinskaja P9 / JKRT S9",
        "coxiter_path": P9_PATH,
    },
    23: {
        "genus_symbol": "II_(9,1)(2^4)",
        "determinant": -16,
        "smith_invariants": [1] * 6 + [2] * 4,
        "discriminant_q_mod_2_histogram": {"0": 10, "1": 6},
        "full_chamber_facets": 11,
        "full_chamber_finite_vertices": 24,
        "full_chamber_cusps": 2,
        "classical_chamber": "Tumarkin 9-tum04_11_01",
        "coxiter_path": TUMARKIN_11_01_PATH,
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Return inertia by exact rational congruence elimination."""

    work = sp.Matrix(matrix)
    positive = negative = 0
    while work.rows:
        pivot = next((index for index in range(work.rows) if work[index, index]), None)
        if pivot is None:
            pair = next(
                (
                    (left, right)
                    for left in range(work.rows)
                    for right in range(left + 1, work.rows)
                    if work[left, right]
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            order = [*pair, *(index for index in range(work.rows) if index not in pair)]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0
        if pivot:
            order = [pivot, *(index for index in range(work.rows) if index != pivot)]
            work = work.extract(order, order)
        value = sp.Rational(work[0, 0])
        positive += int(bool(value > 0))
        negative += int(bool(value < 0))
        if work.rows == 1:
            break
        column = work[1:, 0]
        work = sp.simplify(work[1:, 1:] - column * column.T / value)
    return positive, negative, 0


def _two_elementary_invariants(lattice: sp.Matrix) -> tuple[int, int, int] | None:
    if not all(lattice[index, index] % 2 == 0 for index in range(lattice.rows)):
        return None
    inverse = lattice.inv()
    if not all(2 * value in sp.ZZ for value in inverse):
        return None
    smith = _smith(lattice)
    if any(value not in (1, 2) for value in smith):
        return None
    length = sum(value == 2 for value in smith)
    delta = int(any(value.q == 2 for value in inverse.diagonal()))
    return lattice.rows, length, delta


def _d4_negative() -> sp.Matrix:
    gram = -2 * sp.eye(4)
    for left, right in ((0, 1), (1, 2), (1, 3)):
        gram[left, right] = gram[right, left] = 1
    return gram


def _standard_two_elementary_lattices() -> dict[tuple[int, int, int], dict[str, Any]]:
    u = sp.Matrix([[0, 1], [1, 0]])
    standards = {
        (10, 4, 0): ("U + D4(-1)^2", sp.diag(u, _d4_negative(), _d4_negative())),
        (10, 6, 0): (
            "U(2) + D4(-1)^2",
            sp.diag(2 * u, _d4_negative(), _d4_negative()),
        ),
        (10, 8, 0): ("U + E8(-2)", sp.diag(u, -2 * _e8_cartan())),
    }
    result = {}
    for invariants, (name, gram) in standards.items():
        checks = {
            "signature_1_9": _exact_inertia(gram) == (1, 9, 0),
            "invariants_match": _two_elementary_invariants(gram) == invariants,
            "determinant_absolute_2^a": abs(int(gram.det())) == 2 ** invariants[1],
        }
        if not all(checks.values()):
            raise ArithmeticError(f"standard 2-elementary lattice changed: {name}")
        result[invariants] = {
            "isometry_type": name,
            "gram": _payload(gram),
            "checks": checks,
        }
    return result


def _algebro_geometric_interpretation() -> dict[str, Any]:
    """Record the conditional K3 meaning of the three maximal carriers.

    This is deliberately not a promotion of any of the five fixed sections:
    a reflection-stable abstract overclosure need not be the invariant lattice
    of an involution on the original surface.  The statement here is only the
    classical interpretation of the same abstract 2-elementary lattice when
    it is used as the full invariant lattice of a nonsymplectic K3 involution.
    """

    if _sha256(ALEXEEV_ENGEL_SOURCE) != ALEXEEV_ENGEL_SOURCE_SHA256:
        raise ArithmeticError("pinned Alexeev--Engel source archive changed")
    with tarfile.open(ALEXEEV_ENGEL_SOURCE, "r:gz") as archive:
        handle = archive.extractfile("2el.tex")
        if handle is None:
            raise ArithmeticError("Alexeev--Engel source lacks 2el.tex")
        source = handle.read().decode("utf-8")
    source_markers = {
        "fixed_locus_formula": (
            "g = 11 - \\tfrac12 (r+a) \\text{ and } "
            "k = \\tfrac12 (r-a)"
        ),
        "exceptional_10_8_0": (
            "R=C_0\\sqcup C_1$ is the union of two elliptic curves"
        ),
        "rational_elliptic_model": (
            "\\item $\\I_0^2$: $(10,8,0)$. $R=F_0\\cup F$ is a union "
            "of two smooth"
        ),
    }
    if not all(marker in source for marker in source_markers.values()):
        raise ArithmeticError("Alexeev--Engel K3 source markers changed")

    rows: dict[str, dict[str, Any]] = {}
    lattice_names = {
        (10, 4, 0): "U + D4(-1)^2",
        (10, 6, 0): "U(2) + D4(-1)^2",
        (10, 8, 0): "U + E8(-2)",
    }
    for invariants, lattice_name in lattice_names.items():
        rank, length, delta = invariants
        genus_sum = 11 - (rank + length) // 2
        rational_component_count = (rank - length) // 2
        row: dict[str, Any] = {
            "invariants_r_a_delta": list(invariants),
            "lattice": lattice_name,
            "nonsymplectic_K3_moduli_dimension": 20 - rank,
            "generic_Picard_and_invariant_lattice": lattice_name,
        }
        if invariants == (10, 8, 0):
            row.update(
                {
                    "fixed_locus": "two disjoint smooth elliptic curves",
                    "fixed_locus_exceptional_case": True,
                    "quotient_model": (
                        "rational elliptic surface; the ramification divisor "
                        "is two smooth fibers"
                    ),
                }
            )
        else:
            row.update(
                {
                    "fixed_locus": (
                        f"one smooth genus-{genus_sum} curve and "
                        f"{rational_component_count} disjoint rational curves"
                    ),
                    "fixed_locus_exceptional_case": False,
                    "fixed_curve_genus": genus_sum,
                    "additional_rational_curve_count": rational_component_count,
                }
            )
        rows[str(invariants)] = row

    return {
        "scope_caveat": (
            "This describes the abstract maximal carrier only if it is promoted "
            "to the full invariant lattice S=H^2(X,Z)^+ of a nonsymplectic K3 "
            "involution; the five reflection-stable overclosure computations "
            "alone do not construct such involutions on the original sections."
        ),
        "classification_status": (
            "all three are classical Nikulin 2-elementary K3 involution types"
        ),
        "primary_source": (
            "V. Alexeev and P. Engel, Compactifications of moduli spaces of K3 "
            "surfaces with a nonsymplectic involution, arXiv:2208.10383v2, "
            "Sections 2A and 4B"
        ),
        "pinned_primary_source_archive": str(ALEXEEV_ENGEL_SOURCE),
        "pinned_primary_source_archive_sha256": ALEXEEV_ENGEL_SOURCE_SHA256,
        "source_marker_checks": {key: True for key in source_markers},
        "types": rows,
        "equation_family_status": (
            "known ten-dimensional moduli and quotient constructions are "
            "source-backed; no explicit polynomial normal form is frozen or "
            "claimed in this certificate"
        ),
    }


def _relation_token(ratio: sp.Rational) -> str:
    finite = {
        sp.Rational(0): "2",
        sp.Rational(1, 4): "3",
        sp.Rational(1, 2): "4",
        sp.Rational(3, 4): "6",
        sp.Rational(1): "parallel",
    }
    if ratio in finite:
        return finite[ratio]
    if ratio > 1:
        return f"ultraparallel:{ratio}"
    raise ArithmeticError(f"unsupported Coxeter ratio {ratio}")


def _candidate_relations(case: dict[str, Any]) -> tuple[tuple[str, ...], ...]:
    gram = _matrix(case["section"]["restricted_mirrors"]["positive_wall_gram"])
    return tuple(
        tuple(
            "1"
            if left == right
            else _relation_token(
                sp.factor(gram[left, right] ** 2 / (gram[left, left] * gram[right, right]))
            )
            for right in range(gram.rows)
        )
        for left in range(gram.rows)
    )


def _relation_histogram(relations: Sequence[Sequence[str]]) -> dict[str, int]:
    histogram = Counter(
        relations[left][right]
        for left in range(len(relations))
        for right in range(left)
    )
    return dict(sorted(histogram.items()))


def _relation_graph(relations: Sequence[Sequence[str]]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(relations)))
    for left in range(len(relations)):
        for right in range(left):
            if relations[left][right] != "2":
                graph.add_edge(left, right, relation=relations[left][right])
    return graph


EDGE_MATCH = nx.algorithms.isomorphism.categorical_edge_match("relation", None)


def _isomorphic(
    left: Sequence[Sequence[str]], right: Sequence[Sequence[str]]
) -> bool:
    return nx.is_isomorphic(
        _relation_graph(left), _relation_graph(right), edge_match=EDGE_MATCH
    )


def _parse_coxiter_text(text: str) -> tuple[int, int, tuple[tuple[str, ...], ...]]:
    lines = [
        content
        for line in text.splitlines()
        if (content := line.split("#", 1)[0].strip())
    ]
    facets, dimension = map(int, lines[0].split()[:2])
    matrix = [
        ["1" if left == right else "2" for right in range(facets)]
        for left in range(facets)
    ]
    vertex_line = next(
        (line for line in lines[1:] if line.startswith("vertices labels:")), None
    )
    labels = (
        [str(index + 1) for index in range(facets)]
        if vertex_line is None
        else vertex_line.split(":", 1)[1].split()
    )
    positions = {label: index for index, label in enumerate(labels)}
    if len(positions) != facets:
        raise ArithmeticError("invalid CoxIter vertex labels")
    replacements = {"0": "parallel", "infinity": "parallel", "∞": "parallel"}
    for line in lines[1:]:
        if line.startswith("vertices labels:"):
            continue
        left_label, right_label, raw_relation = line.split()
        left, right = positions[left_label], positions[right_label]
        relation = replacements.get(raw_relation, raw_relation)
        if relation == "1":
            relation = "ultraparallel:coxiter_1"
        matrix[left][right] = matrix[right][left] = relation
    return facets, dimension, tuple(tuple(row) for row in matrix)


def _validate_mapping(
    left: Sequence[Sequence[str]],
    right: Sequence[Sequence[str]],
    mapping_one_based: Sequence[int],
) -> bool:
    if len(left) != len(right) or len(mapping_one_based) != len(left):
        return False
    mapping = [value - 1 for value in mapping_one_based]
    return sorted(mapping) == list(range(len(left))) and all(
        left[i][j] == right[mapping[i]][mapping[j]]
        for i in range(len(left))
        for j in range(len(left))
    )


def _scan_coxiter(
    candidates: dict[str, tuple[tuple[str, ...], ...]],
) -> tuple[dict[str, Any], dict[str, tuple[tuple[str, ...], ...]]]:
    if not (
        _sha256(COXITER_ARCHIVE) == COXITER_ARCHIVE_SHA256
        and _sha256(DOUBLING_ARCHIVE) == DOUBLING_ARCHIVE_SHA256
    ):
        raise ArithmeticError("pinned CoxIter snapshots changed")
    target_facets = {len(relations) for relations in candidates.values()}
    records: list[dict[str, Any]] = []
    named: dict[str, tuple[tuple[str, ...], ...]] = {}
    ordinary_count = 0
    with zipfile.ZipFile(COXITER_ARCHIVE) as archive:
        for path in sorted(name for name in archive.namelist() if name.endswith(".coxiter")):
            ordinary_count += 1
            facets, dimension, relations = _parse_coxiter_text(
                archive.read(path).decode("utf-8", errors="replace")
            )
            if path in (P9_PATH, TUMARKIN_11_01_PATH):
                named[path] = relations
            if dimension == 9 and facets in target_facets:
                records.append(
                    {
                        "lane": "ordinary",
                        "path": path,
                        "facets": facets,
                        "relations": relations,
                    }
                )
    doubling_count = 0
    with tarfile.open(DOUBLING_ARCHIVE, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            doubling_count += 1
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"cannot read {member.name}")
            facets, dimension, relations = _parse_coxiter_text(
                handle.read().decode("utf-8", errors="replace")
            )
            if dimension == 9 and facets in target_facets:
                records.append(
                    {
                        "lane": "doubling",
                        "path": member.name,
                        "facets": facets,
                        "relations": relations,
                    }
                )
    if ordinary_count != 1214 or doubling_count != 2959:
        raise ArithmeticError("pinned CoxIter file inventory changed")

    comparisons = {}
    for identifier, candidate in candidates.items():
        same = [record for record in records if record["facets"] == len(candidate)]
        histogram = _relation_histogram(candidate)
        collisions = [
            record["path"]
            for record in same
            if _relation_histogram(record["relations"]) == histogram
        ]
        matches = [
            record["path"]
            for record in same
            if _isomorphic(candidate, record["relations"])
        ]
        comparisons[identifier] = {
            "hyperbolic_dimension": 9,
            "facet_count": len(candidate),
            "candidate_relation_histogram": histogram,
            "ordinary_same_header_record_count": sum(
                record["lane"] == "ordinary" for record in same
            ),
            "doubling_same_header_record_count": sum(
                record["lane"] == "doubling" for record in same
            ),
            "same_header_record_count": len(same),
            "relation_histogram_collision_paths": collisions,
            "exact_labelled_graph_match_paths": matches,
        }

    expected_matches = {
        identifier: ([IMHOF_MATCH] if identifier == "A14_8_0_C01" else [])
        for identifier in candidates
    }
    if any(
        comparisons[identifier]["exact_labelled_graph_match_paths"] != expected
        for identifier, expected in expected_matches.items()
    ):
        raise ArithmeticError("pinned CoxIter exact-match verdict changed")
    if (
        comparisons["A14_8_0_C01"]["relation_histogram_collision_paths"]
        != IMHOF_HISTOGRAM_COLLISIONS
    ):
        raise ArithmeticError("Im Hof histogram-collision stock changed")
    imhof_relations = next(
        record["relations"] for record in records if record["path"] == IMHOF_MATCH
    )
    if not _validate_mapping(
        candidates["A14_8_0_C01"],
        imhof_relations,
        IMHOF_CANDIDATE_TO_REFERENCE_ONE_BASED,
    ):
        raise ArithmeticError("explicit Im Hof wall permutation changed")
    comparisons["A14_8_0_C01"]["exact_match"] = {
        "path": IMHOF_MATCH,
        "candidate_to_reference_wall_permutation_one_based": (
            IMHOF_CANDIDATE_TO_REFERENCE_ONE_BASED
        ),
        "classical_identification": (
            "Im Hof's doubly truncated orthoscheme, announced in 1985 and "
            "tabulated in the 1990 Napier-cycle paper"
        ),
    }
    return (
        {
            "scope": "bounded corpus comparison, not a worldwide priority theorem",
            "repository": "https://github.com/rgugl/CoxIter",
            "commit": COXITER_COMMIT,
            "ordinary_archive_sha256": COXITER_ARCHIVE_SHA256,
            "doubling_archive_sha256": DOUBLING_ARCHIVE_SHA256,
            "ordinary_coxiter_file_count": ordinary_count,
            "doubling_coxiter_file_count": doubling_count,
            "comparisons": comparisons,
        },
        named,
    )


def _crystallographic_gram(
    relations: Sequence[Sequence[str]], root_squares: Sequence[int]
) -> sp.Matrix:
    gram = sp.diag(*root_squares)
    for left in range(len(relations)):
        for right in range(left):
            relation = relations[left][right]
            if relation == "2":
                value = 0
            elif relation == "3":
                if root_squares[left] != root_squares[right]:
                    raise ArithmeticError("m3 edge joins unequal root squares")
                value = -root_squares[left] // 2
            elif relation == "4":
                if sorted((root_squares[left], root_squares[right])) != [2, 4]:
                    raise ArithmeticError("m4 edge lacks the crystallographic 2:1 ratio")
                value = -2
            else:
                raise ArithmeticError(f"named crystallographic stock has {relation}")
            gram[left, right] = gram[right, left] = value
    return gram


def _sk_snapshot(named: dict[str, tuple[tuple[str, ...], ...]]) -> dict[str, Any]:
    if not (
        _sha256(SK_LIST9) == SK_LIST9_SHA256
        and _sha256(SK_GRAM9) == SK_GRAM9_SHA256
    ):
        raise ArithmeticError("pinned Scharlau--Kirschmer rank-ten snapshot changed")
    rows = ast.literal_eval(SK_GRAM9.read_text(encoding="utf-8"))
    target_grams = {
        entry: sp.Matrix(10, 10, rows[entry - 1]) for entry in SK_TARGETS
    }
    target_rows = {}
    for entry, target in SK_TARGETS.items():
        gram = target_grams[entry]
        checks = {
            "signature_9_1": _exact_inertia(gram) == (9, 1, 0),
            "determinant": int(gram.det()) == target["determinant"],
            "smith": _smith(gram) == target["smith_invariants"],
            "q_histogram": _q_histogram(gram)
            == target["discriminant_q_mod_2_histogram"],
            "two_elementary_invariants": _two_elementary_invariants(gram)
            == (10, len([x for x in target["smith_invariants"] if x == 2]), 0),
        }
        if not all(checks.values()):
            raise ArithmeticError(f"Scharlau--Kirschmer List9 entry {entry} changed")
        target_rows[str(entry)] = {**target, "gram": _payload(gram), "checks": checks}

    p9_relations = named[P9_PATH]
    p9_gram = _crystallographic_gram(p9_relations, [2] * 8 + [4, 2])
    p9_checks = {
        "H9_simplex": len(p9_relations) == 10,
        "root_gram_full_rank": p9_gram.rank() == 10,
        "root_lattice_matches_entry7_invariants": (
            int(p9_gram.det()) == -4
            and _smith(p9_gram) == SK_TARGETS[7]["smith_invariants"]
            and _q_histogram(p9_gram)
            == SK_TARGETS[7]["discriminant_q_mod_2_histogram"]
        ),
    }
    tumarkin_relations = named[TUMARKIN_11_01_PATH]
    tumarkin_gram = _crystallographic_gram(
        tumarkin_relations, [4, 4] + [2] * 7 + [4, 4]
    )
    expected_null = sp.Matrix([-1, -2, -3, -2, -1, 0, 1, 2, 3, 2, 1])
    tumarkin_basis_gram = tumarkin_gram[1:, 1:]
    tumarkin_checks = {
        "H9_with_11_facets": len(tumarkin_relations) == 11,
        "root_gram_rank_10": tumarkin_gram.rank() == 10,
        "primitive_integral_null_relation": tumarkin_gram * expected_null
        == sp.zeros(11, 1),
        "nodes_2_through_11_span_all_roots": abs(int(expected_null[0])) == 1,
        "root_lattice_matches_entry23_invariants": (
            int(tumarkin_basis_gram.det()) == -16
            and _smith(tumarkin_basis_gram) == SK_TARGETS[23]["smith_invariants"]
            and _q_histogram(tumarkin_basis_gram)
            == SK_TARGETS[23]["discriminant_q_mod_2_histogram"]
        ),
    }
    if not all((*p9_checks.values(), *tumarkin_checks.values())):
        raise ArithmeticError("named full H9 chamber stock changed")
    return {
        "official_census_url": (
            "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html"
        ),
        "List9_sha256": SK_LIST9_SHA256,
        "Gram9_sha256": SK_GRAM9_SHA256,
        "classification_scope": (
            "reflective strongly squarefree Lorentzian genera in the published "
            "Scharlau--Kirschmer census; the two exact 2-elementary classes are "
            "isolated by signature and discriminant form"
        ),
        "targets": target_rows,
        "named_full_chamber_replays": {
            "Vinberg_Kaplinskaja_P9": {
                "coxiter_path": P9_PATH,
                "root_squares": [2] * 8 + [4, 2],
                "root_gram": _payload(p9_gram),
                "checks": p9_checks,
            },
            "Tumarkin_9_tum04_11_01": {
                "coxiter_path": TUMARKIN_11_01_PATH,
                "root_squares": [4, 4] + [2] * 7 + [4, 4],
                "primitive_null_relation": [int(value) for value in expected_null],
                "basis_nodes_one_based": list(range(2, 12)),
                "basis_gram": _payload(tumarkin_basis_gram),
                "checks": tumarkin_checks,
            },
        },
        "target_grams_internal": target_grams,
    }


def _case_overclosure(
    case: dict[str, Any],
    relations: tuple[tuple[str, ...], ...],
    standards: dict[tuple[int, int, int], dict[str, Any]],
    sk_targets: dict[int, sp.Matrix],
) -> dict[str, Any]:
    identifier = case["id"]
    expected_integral, expected_even, expected_index, expected_det, expected_invariants = (
        EXPECTED_STABLE[identifier]
    )
    section = case["section"]
    gram = _matrix(section["saturated_fixed_lattice"]["gram"])
    roots = [
        sp.Matrix(record["primitive_integral_same_direction_coordinates"])
        for record in section["restricted_mirrors"]["records"]
        if record["hyperbolic_status"] == "active_facet"
    ]
    reflections = [_reflection(gram, root) for root in roots]
    audit = _stable_subgroups(gram, reflections)
    groups = audit.pop("groups_internal")
    left_inverse = audit.pop("left_smith_inverse")
    if (
        audit["stable_integral_subgroup_orders"] != expected_integral
        or audit["stable_even_subgroup_orders"] != expected_even
    ):
        raise ArithmeticError(f"stable subgroup census changed for {identifier}")
    maximal_positions = audit["maximal_stable_subgroup_positions_zero_based"]
    if len(maximal_positions) != 1:
        raise ArithmeticError(f"maximal stable overclosure not unique for {identifier}")
    position = maximal_positions[0]
    maximal_group = groups[position]
    basis = _overlattice_basis(
        gram,
        left_inverse,
        audit["stable_subgroups"][position]["cosets_in_smith_coordinates"],
    )
    maximal_gram = sp.simplify(basis.T * gram * basis)
    maximal_actions = [sp.simplify(basis.inv() * action * basis) for action in reflections]
    maximal_invariants = _two_elementary_invariants(maximal_gram)
    maximal_checks = {
        "index": len(maximal_group) == expected_index,
        "basis_index": abs(sp.Rational(1, basis.det())) == expected_index,
        "determinant_square_law": maximal_gram.det()
        == gram.det() / expected_index**2,
        "expected_determinant": int(maximal_gram.det()) == expected_det,
        "expected_two_elementary_invariants": maximal_invariants
        == expected_invariants,
        "standard_isometry_type_by_Nikulin_uniqueness": (
            maximal_invariants in standards
            and _exact_inertia(maximal_gram) == (1, 9, 0)
        ),
        "all_reflection_actions_integral": all(
            value.is_Integer for action in maximal_actions for value in action
        ),
        "all_reflection_actions_unimodular": all(
            abs(action.det()) == 1 for action in maximal_actions
        ),
        "all_reflection_actions_preserve_maximal_gram": all(
            action.T * maximal_gram * action == maximal_gram
            for action in maximal_actions
        ),
    }
    if not all(maximal_checks.values()):
        raise ArithmeticError(f"maximal stable overclosure replay failed for {identifier}")

    if abs(expected_det) == 16:
        normalized_gram = -maximal_gram
        normalized_actions = maximal_actions
        operation = "global sign reversal -M"
        target_entry = 23
    else:
        normalized_gram = sp.simplify(-2 * maximal_gram.inv())
        normalized_actions = [sp.simplify(action.inv().T) for action in maximal_actions]
        operation = "scaled dual M^*(2), with Gram -2 M^(-1)"
        target_entry = 7 if abs(expected_det) == 256 else 23
    target = SK_TARGETS[target_entry]
    target_gram = sk_targets[target_entry]
    normalized_checks = {
        "integral": all(value.is_Integer for value in normalized_gram),
        "signature_9_1": _exact_inertia(normalized_gram) == (9, 1, 0),
        "determinant": int(normalized_gram.det()) == target["determinant"],
        "smith": _smith(normalized_gram) == target["smith_invariants"],
        "q_histogram": _q_histogram(normalized_gram)
        == target["discriminant_q_mod_2_histogram"],
        "same_2_elementary_invariants_as_SK_target": (
            _two_elementary_invariants(normalized_gram)
            == _two_elementary_invariants(target_gram)
        ),
        "Nikulin_uniqueness_hypotheses": (
            _two_elementary_invariants(normalized_gram) is not None
            and normalized_gram.rows >= _two_elementary_invariants(normalized_gram)[1] + 2
            and _exact_inertia(normalized_gram)[0] > 0
            and _exact_inertia(normalized_gram)[1] > 0
        ),
        "all_normalized_reflection_actions_integral": all(
            value.is_Integer for action in normalized_actions for value in action
        ),
        "all_normalized_reflection_actions_preserve_gram": all(
            action.T * normalized_gram * action == normalized_gram
            for action in normalized_actions
        ),
    }
    if not all(normalized_checks.values()):
        raise ArithmeticError(f"SK normalization changed for {identifier}")

    computed_doubling = [
        index + 1
        for index, row in enumerate(relations)
        if all(token not in {"3", "5", "7", "9"} for token in row)
    ]
    if computed_doubling != case["doubling_walls_one_based"]:
        raise ArithmeticError(f"doubling-wall census changed for {identifier}")
    disjoint_pairs = [
        [left, right]
        for left, right in itertools.combinations(computed_doubling, 2)
        if relations[left - 1][right - 1] == "parallel"
        or relations[left - 1][right - 1].startswith("ultraparallel:")
    ]
    if disjoint_pairs != EXPECTED_DISJOINT_DOUBLING_PAIRS[identifier]:
        raise ArithmeticError(f"redoubling pair stock changed for {identifier}")
    disjoint_triples = [
        list(triple)
        for triple in itertools.combinations(computed_doubling, 3)
        if all(
            relations[left - 1][right - 1] == "parallel"
            or relations[left - 1][right - 1].startswith("ultraparallel:")
            for left, right in itertools.combinations(triple, 2)
        )
    ]

    odd_reference = sp.diag(1, *([-1] * 9))
    rational = _compare(
        gram,
        odd_reference,
        Fraction(int(gram.det()), int(odd_reference.det())),
    )
    if not rational["rationally_similar"]:
        raise ArithmeticError(f"odd-unimodular rational class changed for {identifier}")
    standard = standards[expected_invariants]
    redoublable = bool(disjoint_pairs)
    return {
        "id": identifier,
        "hyperbolic_dimension": 9,
        "facet_count": len(relations),
        "relation_histogram": _relation_histogram(relations),
        "original_fixed_lattice": {
            "determinant": int(gram.det()),
            "smith_invariants": _smith(gram),
        },
        "exhaustive_discriminant_module_audit": audit,
        "unique_maximal_reflection_stable_overclosure": {
            "unique": True,
            "index_over_original_lattice": expected_index,
            "basis_columns_in_original_lattice": _payload(basis),
            "gram": _payload(maximal_gram),
            "determinant": int(maximal_gram.det()),
            "smith_invariants": _smith(maximal_gram),
            "discriminant_q_mod_2_histogram": _q_histogram(maximal_gram),
            "two_elementary_invariants_r_a_delta": list(expected_invariants),
            "isometry_type": standard["isometry_type"],
            "identification_method": (
                "exact signature/discriminant-form replay plus Nikulin uniqueness "
                "for indefinite even 2-elementary lattices"
            ),
            "facet_reflection_actions_sha256": _canonical_sha256(
                [_payload(action) for action in maximal_actions]
            ),
            "checks": maximal_checks,
        },
        "normalization_to_strongly_squarefree_class": {
            "operation": operation,
            "not_an_overlattice_when_scaled_dual": abs(expected_det) != 16,
            "gram": _payload(normalized_gram),
            "scharlau_kirschmer_list": "List9",
            "scharlau_kirschmer_entry_one_based": target_entry,
            "genus_symbol": target["genus_symbol"],
            "known_full_reflection_chamber": target["classical_chamber"],
            "known_full_chamber_facets": target["full_chamber_facets"],
            "known_full_chamber_finite_vertices": target[
                "full_chamber_finite_vertices"
            ],
            "known_full_chamber_cusps": target["full_chamber_cusps"],
            "section_reflection_group_has_finite_index_in_full_group": True,
            "numerical_index": None,
            "numerical_index_status": (
                "open here: Coxeter Euler characteristics vanish in odd "
                "hyperbolic dimension H9"
            ),
            "normalized_reflection_actions_sha256": _canonical_sha256(
                [_payload(action) for action in normalized_actions]
            ),
            "checks": normalized_checks,
        },
        "rational_similarity_to_odd_unimodular_I_9_1": rational,
        "arithmetic_conclusion": (
            "classical odd-unimodular rational commensurability class; no new "
            "arithmetic commensurability class or maximal reflection group"
        ),
        "redoubling": {
            "doubling_walls_one_based": computed_doubling,
            "disjoint_doubling_wall_pairs_one_based": disjoint_pairs,
            "all_disjoint_pairs_are_parallel": all(
                relations[left - 1][right - 1] == "parallel"
                for left, right in disjoint_pairs
            ),
            "redoublable_by_Allcock_definition": redoublable,
            "Allcock_Lemma_1_3_iterated_doubling_applies": redoublable,
            "infinite_finite_volume_H9_descendants": redoublable,
            "successive_volume_multipliers": "2^k for k >= 0" if redoublable else None,
            "Allcock_Corollary_4_2_reflection_subgroups_of_every_positive_index": (
                redoublable
            ),
            "pairwise_disjoint_doubling_wall_triples_one_based": disjoint_triples,
            "Allcock_Theorem_4_3_exponential_count_applies_from_certified_walls": (
                bool(disjoint_triples)
            ),
            "novelty_scope": (
                "Allcock's doubling/redoubling mechanism is classical; only the "
                "exact Alexeev fixed-section carrier is project-level new"
            ),
        },
    }


def build_certificate() -> dict[str, Any]:
    if _sha256(SOURCE) != SOURCE_SHA256:
        raise ArithmeticError("immutable five-seed source archive changed")
    allcock_hashes = {
        "pdf": _sha256(ALLCOCK_PDF),
        "arxiv_source_tar": _sha256(ALLCOCK_SOURCE),
        "figure4_pdf": _sha256(ALLCOCK_FIGURE4),
    }
    if allcock_hashes != {
        "pdf": ALLCOCK_PDF_SHA256,
        "arxiv_source_tar": ALLCOCK_SOURCE_SHA256,
        "figure4_pdf": ALLCOCK_FIGURE4_SHA256,
    }:
        raise ArithmeticError("pinned Allcock primary-source snapshot changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    candidates = {case["id"]: _candidate_relations(case) for case in source["cases"]}
    coxiter, named = _scan_coxiter(candidates)
    sk = _sk_snapshot(named)
    target_grams = sk.pop("target_grams_internal")
    standards = _standard_two_elementary_lattices()
    ag_interpretation = _algebro_geometric_interpretation()
    cases = [
        _case_overclosure(case, candidates[case["id"]], standards, target_grams)
        for case in source["cases"]
    ]
    by_id = {case["id"]: case for case in cases}

    literature = {
        "Im_Hof": {
            "primary_sources": [
                {
                    "citation": (
                        "H.-C. Im Hof, A class of hyperbolic Coxeter groups, "
                        "Expositiones Mathematicae 3 (1985), 179--186"
                    ),
                    "role": "announcement and original finite-covolume class",
                },
                {
                    "citation": (
                        "H.-C. Im Hof, Napier cycles and hyperbolic Coxeter groups, "
                        "Bull. Soc. Math. Belgique (Serie A) 42 (1990), 523--545"
                    ),
                    "role": "complete tables, including doubly truncated orthoschemes",
                },
            ],
            "exact_project_collision": "A14_8_0_C01",
            "coxiter_path": IMHOF_MATCH,
        },
        "Allcock_2006": {
            "citation": (
                "D. Allcock, Infinitely many hyperbolic Coxeter groups through "
                "dimension 19, Geometry & Topology 10 (2006), 737--758"
            ),
            "doi": "10.2140/gt.2006.10.737",
            "arxiv": "0903.0138",
            "pinned_hashes": allcock_hashes,
            "definition_scope": (
                "walls are disjoint when they do not meet in H^n; parallel walls "
                "therefore qualify even though their closures share an ideal point"
            ),
            "Lemma_1_3": (
                "doubling across a doubling wall preserves the Coxeter property, "
                "and a disjoint second doubling wall makes the double redoublable"
            ),
            "Corollary_4_2": (
                "a redoublable chamber group has reflection-generated subgroups "
                "of every positive index"
            ),
            "Theorem_4_3_requirement": (
                "three or more pairwise disjoint doubling walls for the exponential "
                "subgroup-count conclusion"
            ),
            "Figure_4_1_H9_Q": {
                "construction": "double of the D6D6D4 face F across wall 04|213",
                "facet_count_from_vector_figure": 14,
                "none_of_five_seeds_is_literal_Q_by_facet_count": all(
                    case["facet_count"] != 14 for case in cases
                ),
            },
        },
        "Bredon_ADEG": {
            "primary_source": "N. Bredon, On ADEG-polyhedra in hyperbolic spaces, arXiv:2507.05153v1",
            "classified_hypotheses": (
                "mutually intersecting facets and only dihedral labels 2, 3, 6"
            ),
            "exceptional_P_star": {"hyperbolic_dimension": 9, "facet_count": 14},
            "all_five_outside_ADEG_because_each_uses_label_4": all(
                case["relation_histogram"].get("4", 0) > 0 for case in cases
            ),
            "none_literal_P_star_by_facet_count": all(
                case["facet_count"] != 14 for case in cases
            ),
        },
        "Tumarkin": {
            "n_plus_2_primary_source": (
                "P. Tumarkin, Hyperbolic Coxeter n-polytopes with n+2 facets, "
                "arXiv:math/0301133"
            ),
            "n_plus_2_is_complete_finite_volume_classification": True,
            "literal_H9_n_plus_2_facet_count": 11,
            "entry23_named_full_chamber": "9-tum04_11_01",
            "none_of_five_literal_n_plus_2_by_facet_count": all(
                case["facet_count"] != 11 for case in cases
            ),
            "compact_n_plus_3_primary_source": "P. Tumarkin, arXiv:math/0406226",
            "compact_n_plus_3_does_not_apply_to_noncompact_five": True,
        },
        "Vinberg_Kaplinskaja_P9": {
            "quadratic_form": "I_(9,1), equivalently -x0^2+x1^2+...+x9^2",
            "literal_chamber_facets": 10,
            "coxiter_path": P9_PATH,
            "none_of_five_literal_P9_by_facet_count": all(
                case["facet_count"] != 10 for case in cases
            ),
            "A14_10_1_C01_scaled_dual_class_is_P9_entry7": (
                by_id["A14_10_1_C01"]["normalization_to_strongly_squarefree_class"][
                    "scharlau_kirschmer_entry_one_based"
                ]
                == 7
            ),
        },
        "Alexeev": {
            "primary_source": "V. Alexeev, Reflective hyperbolic 2-elementary lattices, arXiv:2209.09110",
            "published_object": "complete ambient mixed-root reflection matrices",
            "project_operation": (
                "exact fixed-section extraction under finite ambient actions, "
                "followed by irredundant facet and volume certification"
            ),
            "novelty_distinction": (
                "the extraction may be new project-level algebraic data even when "
                "the abstract chamber or arithmetic class is classical"
            ),
        },
    }

    checks = {
        "five_cases": len(cases) == 5,
        "all_unique_maximal_stable_overclosures": all(
            case["unique_maximal_reflection_stable_overclosure"]["unique"]
            for case in cases
        ),
        "all_exactly_classified_2_elementary": all(
            tuple(
                case["unique_maximal_reflection_stable_overclosure"][
                    "two_elementary_invariants_r_a_delta"
                ]
            )
            in standards
            for case in cases
        ),
        "all_normalize_to_SK_entry7_or23": all(
            case["normalization_to_strongly_squarefree_class"][
                "scharlau_kirschmer_entry_one_based"
            ]
            in (7, 23)
            for case in cases
        ),
        "all_in_odd_unimodular_rational_class": all(
            case["rational_similarity_to_odd_unimodular_I_9_1"]["rationally_similar"]
            for case in cases
        ),
        "only_Im_Hof_literal_CoxIter_match": sum(
            bool(row["exact_labelled_graph_match_paths"])
            for row in coxiter["comparisons"].values()
        )
        == 1,
        "three_redoublable_seeds": sum(
            case["redoubling"]["redoublable_by_Allcock_definition"] for case in cases
        )
        == 3,
        "no_certified_three_disjoint_wall_stock": all(
            not case["redoubling"][
                "Allcock_Theorem_4_3_exponential_count_applies_from_certified_walls"
            ]
            for case in cases
        ),
        "no_new_arithmetic_class": all(
            case["arithmetic_conclusion"].startswith(
                "classical odd-unimodular rational commensurability class"
            )
            for case in cases
        ),
        "three_classical_rank10_K3_involution_types": (
            len(ag_interpretation["types"]) == 3
            and all(
                row["nonsymplectic_K3_moduli_dimension"] == 10
                for row in ag_interpretation["types"].values()
            )
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global lineage checks failed: {checks}")
    return {
        "schema": SCHEMA,
        "source": {
            "five_seed_archive": str(SOURCE.relative_to(ROOT)),
            "five_seed_archive_sha256": SOURCE_SHA256,
        },
        "scope": {
            "certified": (
                "exact reflection-stable integral overclosures, exact rational "
                "quadratic class, pinned-corpus labelled-graph comparisons, and "
                "exact applicability of Allcock redoubling to the displayed walls"
            ),
            "literature_derived": (
                "named chamber attributions and theorem scopes are pinned to the "
                "listed primary sources"
            ),
            "not_claimed": (
                "worldwide abstract-chamber priority, a new commensurability class, "
                "a new maximal reflection group, or a new doubling theorem"
            ),
        },
        "method": {
            "overclosures": (
                "enumerate every totally integral discriminant-bilinear subgroup "
                "stable under all facet reflections and close invariant subgroups "
                "under sums to exhaustion"
            ),
            "integral_classification": (
                "signature and the exact (r,a,delta) discriminant form, followed "
                "by Nikulin uniqueness for indefinite even 2-elementary lattices"
            ),
            "rational_classification": (
                "exact determinant square class and every local Hasse invariant"
            ),
            "diagram_comparison": (
                "exact edge-labelled graph isomorphism against every same-header "
                "H9 record in two pinned CoxIter snapshots"
            ),
        },
        "standard_two_elementary_carriers": {
            str(invariants): row for invariants, row in standards.items()
        },
        "conditional_algebro_geometric_interpretation": ag_interpretation,
        "scharlau_kirschmer_snapshot": sk,
        "coxiter_bounded_comparison": coxiter,
        "cases": cases,
        "primary_source_lineage": literature,
        "summary": {
            "case_count": 5,
            "classical_abstract_chamber_identifications": {
                "A14_8_0_C01": "Im Hof 9-imhof85_03"
            },
            "bounded_unmatched_literal_diagram_candidates": [
                "A14_8_1_C02",
                "A14_10_1_C01",
                "A14_10_1_C02",
                "A16_14_1_C01",
            ],
            "redoublable_seed_ids": [
                case["id"]
                for case in cases
                if case["redoubling"]["redoublable_by_Allcock_definition"]
            ],
            "new_arithmetic_commensurability_class_count": 0,
            "new_maximal_reflection_group_count": 0,
            "worldwide_priority_status": "open for four bounded-unmatched diagrams",
            "safe_novelty_wording": (
                "new exact Alexeev fixed-section/divisor-reflection extractions; "
                "one abstract chamber is Im Hof's and all five arithmetic classes "
                "have classical full reflection overgroups"
            ),
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archived parked-H9 lineage/overclosure certificate is stale")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(rendered, encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
