#!/usr/bin/env python3
"""Certify the worldwide lineage of the four triangular Looijenga H4 types.

The upstream compiler produces four integral Coxeter wall systems in
``I_(1,4)``.  This script forgets the wall norms when comparing abstract
Coxeter polytopes, canonically compares their Coxeter matrices with pinned
literature transcriptions, and computes exact Euler characteristics, volumes,
and congruence-subgroup indices.

The distinction between a complete and a bounded search is intentional:

* Tumarkin's ``n+2`` theorem is complete for the two six-facet types;
* Ma--Zheng's 331-row census is complete for the seven-facet type;
* no complete finite-volume, noncompact H4/eight-facet census is known, so
  the CoxIter exclusion for T8 is only a pinned-corpus result.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import tarfile
from typing import Any, Iterable, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
COMPILER = (
    ROOT / "data" / "enumerations" / "triangular_looijenga_h4_compiler.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "triangular_looijenga_h4_worldwide_lineage.json"
)

DEFAULT_COXITER = Path("/private/tmp/coxiter-current/CoxIter-master")
DEFAULT_COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
DEFAULT_TUMARKIN = Path("/private/tmp/cox-lineage.syY6ZA/nplus2")
DEFAULT_TUMARKIN_ARCHIVE = Path("/private/tmp/cox-lineage.syY6ZA/nplus2.tar")
DEFAULT_MA_ZHENG = Path(
    "/private/tmp/cox-lineage.syY6ZA/seven/paperarXiv/hcp47.tex"
)
DEFAULT_MA_ZHENG_ARCHIVE = Path("/private/tmp/cox-lineage.syY6ZA/seven.tar")
DEFAULT_IM_HOF = Path("/private/tmp/cox-lineage.syY6ZA/napier.pdf")
DEFAULT_2026_ARCHIVE = Path("/private/tmp/cox-lineage.syY6ZA/2607.14715.tar")
DEFAULT_2026_TEX = Path("/private/tmp/cox-lineage.syY6ZA/new2026/main.tex")

COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
DOUBLING_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)
TUMARKIN_ARCHIVE_SHA256 = (
    "b58e9938b987c70fe36b495d904a2a393a54e2e05f5972a4bc595e0a6618f456"
)
TUMARKIN_TEX_SHA256 = (
    "2b290ba3004830d8b96a1e58224b4ce2dc1d742710801ca2e60c4c18b886f565"
)
TUMARKIN_PYR4_SHA256 = (
    "5b5e10ff70ca03f0c1501cba3c55aa93d3bf098609283ce13a58538fef5f408c"
)
MA_ZHENG_ARCHIVE_SHA256 = (
    "13f2d88ad740669f35b4d1859e68da4d489a1223b5ff5a99e128eaa71838234f"
)
MA_ZHENG_TEX_SHA256 = (
    "0017931637f017641d4935450cc2a07d070866d9c749443d9e3ede9f3a39f003"
)
IM_HOF_PDF_SHA256 = (
    "2b3e1c0b80296ef0ea561a2ed0c796ff133b888d2ebe719901fc2f32734fd7ea"
)
ARCHIVE_2026_SHA256 = (
    "d3ae2526a758693e1c3123eff39212ff8eb860bf816a9b116202cd418b0d41ed"
)

TUMARKIN_MATCH_FILES = {
    "T6a": Path(
        "graphs/Tumarkin pyramids n+2/4-tum04_07_02-k=2-l=4.coxiter"
    ),
    "T6b": Path(
        "graphs/Tumarkin pyramids n+2/4-tum04_07_01-k=2-l=4.coxiter"
    ),
}
TUMARKIN_MATCH_SHA256 = {
    "T6a": "41139441b6610343feb1435f83bdb0eedebf4c76225325c2f447d3eb443b6736",
    "T6b": "a1f7c0b3fa2e3f35c7eeab9f8bef662baa1b549a7349e55ec41c942759170718",
}
IM_HOF_MATCH_FILE = Path(
    "graphs/doubly truncated orthoscheme/4-imohf85_13.coxiter"
)
IM_HOF_MATCH_SHA256 = (
    "68905931bad045ce4f2ba00c0cdab4822fdf3901b63d36c9fe35c4674223e59d"
)
P4_FILE = Path("graphs/group units qf/4-1,1,1,1,1.coxiter")
P4_FILE_SHA256 = (
    "a9fa900c9ed1af05f53835c765469c98052ebb3feef056157de36ecd79a44c50"
)

# Simple roots for Vinberg's standard chamber of I_(1,4), in the basis
# diag(1,-1,-1,-1,-1).  The four square-minus-two roots form A4; the
# square-minus-one root is invisible modulo two.
P4_ROOTS = [
    [0, 1, -1, 0, 0],
    [0, 0, 1, -1, 0],
    [0, 0, 0, 1, -1],
    [0, 0, 0, 0, 1],
    [1, -1, -1, -1, 0],
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relation_matrix_from_gram(
    gram: Sequence[Sequence[int]],
) -> list[list[str]]:
    """Return Coxeter labels; ``0`` denotes a parallel pair."""

    size = len(gram)
    result = [
        ["1" if left == right else "2" for right in range(size)]
        for left in range(size)
    ]
    for left in range(size):
        for right in range(left + 1, size):
            product = int(gram[left][right])
            norm_product = int(gram[left][left]) * int(gram[right][right])
            if product == 0:
                label = "2"
            elif 4 * product * product == norm_product:
                label = "3"
            elif 2 * product * product == norm_product:
                label = "4"
            elif product * product == norm_product:
                label = "0"
            else:
                raise ArithmeticError(
                    f"unsupported normalized Gram product at {left},{right}"
                )
            result[left][right] = result[right][left] = label
    return result


def _relation_vector(
    matrix: Sequence[Sequence[str]], permutation: Sequence[int] | None = None
) -> tuple[str, ...]:
    size = len(matrix)
    order = tuple(range(size)) if permutation is None else tuple(permutation)
    return tuple(
        matrix[order[left]][order[right]]
        for left in range(size)
        for right in range(left + 1, size)
    )


def _matrix_from_vector(vector: Sequence[str], size: int) -> list[list[str]]:
    if len(vector) != size * (size - 1) // 2:
        raise ArithmeticError("wrong Coxeter-vector length")
    result = [
        ["1" if left == right else "2" for right in range(size)]
        for left in range(size)
    ]
    cursor = 0
    for left in range(size):
        for right in range(left + 1, size):
            result[left][right] = result[right][left] = str(vector[cursor])
            cursor += 1
    return result


def _canonical_relation_vector(
    matrix: Sequence[Sequence[str]],
) -> tuple[str, ...]:
    return min(
        _relation_vector(matrix, permutation)
        for permutation in itertools.permutations(range(len(matrix)))
    )


def _exact_permutations(
    source: Sequence[Sequence[str]], target: Sequence[Sequence[str]]
) -> list[list[int]]:
    if len(source) != len(target):
        return []
    result: list[list[int]] = []
    for permutation in itertools.permutations(range(len(source))):
        if all(
            source[left][right]
            == target[permutation[left]][permutation[right]]
            for left in range(len(source))
            for right in range(len(source))
        ):
            result.append(list(permutation))
    return result


def _parse_coxiter_text(text: str) -> tuple[int, int, list[list[str]]]:
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
        (line for line in lines[1:] if line.startswith("vertices labels:")),
        None,
    )
    labels = (
        vertex_line.split(":", 1)[1].split()
        if vertex_line is not None
        else [str(index + 1) for index in range(facets)]
    )
    if len(labels) != facets or len(set(labels)) != facets:
        raise ArithmeticError("invalid CoxIter vertex labels")
    index = {label: position for position, label in enumerate(labels)}
    for line in lines[1:]:
        if line.startswith("vertices labels:"):
            continue
        tokens = line.split()
        if len(tokens) != 3:
            raise ArithmeticError(f"unsupported CoxIter line: {line}")
        left, right = index[tokens[0]], index[tokens[1]]
        matrix[left][right] = matrix[right][left] = tokens[2]
    return facets, dimension, matrix


def _parse_coxiter(path: Path) -> tuple[int, int, list[list[str]]]:
    return _parse_coxiter_text(path.read_text(encoding="utf-8"))


def _scan_coxiter(
    corpus: Path, candidates: dict[str, list[list[str]]]
) -> dict[str, Any]:
    ordinary_total = 0
    ordinary_h4: Counter[int] = Counter()
    ordinary_matches: dict[str, list[dict[str, Any]]] = {
        type_id: [] for type_id in candidates
    }
    canonical = {
        type_id: _canonical_relation_vector(matrix)
        for type_id, matrix in candidates.items()
    }
    for path in sorted(corpus.rglob("*.coxiter")):
        facets, dimension, matrix = _parse_coxiter(path)
        ordinary_total += 1
        if dimension == 4:
            ordinary_h4[facets] += 1
        for type_id, target in candidates.items():
            if (dimension, facets) != (4, len(target)):
                continue
            if _canonical_relation_vector(matrix) != canonical[type_id]:
                continue
            permutations = _exact_permutations(target, matrix)
            ordinary_matches[type_id].append(
                {
                    "path": path.relative_to(corpus).as_posix(),
                    "source_file_sha256": _sha256(path),
                    "exact_permutation_count": len(permutations),
                    "lexicographically_first_candidate_to_source_permutation": min(
                        permutations
                    ),
                }
            )

    doubling_path = corpus / "graphs" / "testing_index2.tar.gz"
    doubling_total = 0
    doubling_h4: Counter[int] = Counter()
    doubling_matches: dict[str, list[dict[str, Any]]] = {
        type_id: [] for type_id in candidates
    }
    with tarfile.open(doubling_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"cannot read {member.name}")
            facets, dimension, matrix = _parse_coxiter_text(
                handle.read().decode("utf-8")
            )
            doubling_total += 1
            if dimension == 4:
                doubling_h4[facets] += 1
            for type_id, target in candidates.items():
                if (dimension, facets) != (4, len(target)):
                    continue
                if _canonical_relation_vector(matrix) != canonical[type_id]:
                    continue
                permutations = _exact_permutations(target, matrix)
                doubling_matches[type_id].append(
                    {
                        "path": member.name,
                        "exact_permutation_count": len(permutations),
                        "lexicographically_first_candidate_to_source_permutation": min(
                            permutations
                        ),
                    }
                )
    return {
        "repository": "https://github.com/rgugl/CoxIter",
        "commit": COXITER_COMMIT,
        "archive_sha256": COXITER_ARCHIVE_SHA256,
        "ordinary_file_count": ordinary_total,
        "ordinary_H4_facet_histogram": {
            str(facets): count for facets, count in sorted(ordinary_h4.items())
        },
        "ordinary_exact_matches": ordinary_matches,
        "doubling_archive_sha256": _sha256(doubling_path),
        "doubling_file_count": doubling_total,
        "doubling_H4_facet_histogram": {
            str(facets): count for facets, count in sorted(doubling_h4.items())
        },
        "doubling_exact_matches": doubling_matches,
        "scope": "bounded pinned corpus, not a worldwide priority census",
    }


def _parse_ma_zheng_rows(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    expected_field_lengths = [6, 5, 4, 3, 2, 1]
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        name_match = re.search(r"P_\{([^}]+)\}", line)
        if name_match is None or "&" not in line or r"\\" not in line:
            continue
        cells = line.split("&")
        relation_tokens: list[str] = []
        valid = len(cells) >= 7
        for expected, cell in zip(expected_field_lengths, cells[1:7]):
            cleaned = cell.replace(r"\underline{10}", "X")
            cleaned = re.sub(r"\\color\{[^}]+\}", "", cleaned)
            cleaned = cleaned.replace("$", "").replace("{", "").replace("}", "")
            cleaned = re.sub(r"\\[A-Za-z]+", "", cleaned)
            tokens = [
                "10" if character == "X" else character
                for character in cleaned
                if character.isdigit() or character in "abcX"
            ]
            if len(tokens) != expected:
                valid = False
                break
            relation_tokens.extend(tokens)
        if not valid or len(relation_tokens) != 21:
            continue
        rows.append(
            {
                "id": f"P_{{{name_match.group(1)}}}",
                "line_number_in_pinned_tex": line_number,
                "relation_vector": relation_tokens,
                "raw_row": line.strip(),
            }
        )
    if len(rows) != 331 or len({row["id"] for row in rows}) != 331:
        raise ArithmeticError("Ma--Zheng's 331-row result table did not parse")
    return rows


def _ma_zheng_audit(
    path: Path, candidate: Sequence[Sequence[str]]
) -> dict[str, Any]:
    rows = _parse_ma_zheng_rows(path)
    canonical = _canonical_relation_vector(candidate)
    matches: list[dict[str, Any]] = []
    for row in rows:
        matrix = _matrix_from_vector(row["relation_vector"], 7)
        if _canonical_relation_vector(matrix) != canonical:
            continue
        permutations = _exact_permutations(candidate, matrix)
        matches.append(
            {
                **row,
                "exact_permutation_count": len(permutations),
                "lexicographically_first_candidate_to_table_permutation": min(
                    permutations
                ),
            }
        )
    if [match["id"] for match in matches] != ["P_{2,1}"]:
        raise ArithmeticError("T7 did not have the unique expected Ma--Zheng match")
    raw = matches[0]["raw_row"]
    if "&3&$\\frac{1}{48}$&A" not in raw:
        raise ArithmeticError("P_{2,1} cusp/volume/arithmetic fields changed")
    return {
        "source": "Ma--Zheng, Finite-volume hyperbolic Coxeter 4-dimensional polytopes with 7 facets",
        "url": "https://arxiv.org/abs/2401.13698",
        "version": "v3, 23 December 2024",
        "source_tex_sha256": _sha256(path),
        "result_row_count": len(rows),
        "distinct_result_id_count": len({row["id"] for row in rows}),
        "exact_unlabelled_match_count": len(matches),
        "matches": matches,
        "published_invariants": {
            "cusp_count": 3,
            "volume": "pi^2/48",
            "arithmetic": True,
        },
        "classification_complete_for_finite_volume_H4_seven_facets": True,
    }


def _connected_components(
    nodes: Iterable[int], relations: Sequence[Sequence[str]]
) -> list[list[int]]:
    remaining = set(nodes)
    components: list[list[int]] = []
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        queue = [start]
        component: list[int] = []
        while queue:
            node = queue.pop()
            component.append(node)
            neighbors = {
                other
                for other in remaining
                if relations[node][other] != "2"
            }
            remaining -= neighbors
            queue.extend(neighbors)
        components.append(sorted(component))
    return components


def _finite_component_order(
    component: Sequence[int], relations: Sequence[Sequence[str]]
) -> int:
    rank = len(component)
    if rank == 1:
        return 2
    edges = [
        (left, right, int(relations[left][right]))
        for left, right in itertools.combinations(component, 2)
        if relations[left][right] != "2"
    ]
    if rank == 2:
        if len(edges) != 1:
            raise ArithmeticError("connected rank-two finite component malformed")
        return 2 * edges[0][2]
    degrees = {node: 0 for node in component}
    for left, right, _label in edges:
        degrees[left] += 1
        degrees[right] += 1
    labels = sorted(label for _left, _right, label in edges)
    if rank == 3:
        orders = {(3, 3): 24, (3, 4): 48, (3, 5): 120}
        try:
            return orders[tuple(labels)]
        except KeyError as error:
            raise ArithmeticError("unknown rank-three finite Coxeter type") from error
    if rank == 4:
        degree_sequence = sorted(degrees.values())
        if degree_sequence == [1, 1, 1, 3] and labels == [3, 3, 3]:
            return 192  # D4
        if degree_sequence != [1, 1, 2, 2]:
            raise ArithmeticError("unknown rank-four finite Coxeter graph")
        if labels == [3, 3, 3]:
            return 120  # A4
        if labels == [3, 3, 4]:
            four_edge = next(
                (left, right) for left, right, label in edges if label == 4
            )
            if 1 in (degrees[four_edge[0]], degrees[four_edge[1]]):
                return 384  # B4
            return 1152  # F4
    raise ArithmeticError("finite component outside the certified small-rank list")


def _finite_special_order(
    subset: Sequence[int],
    integral_gram: Sequence[Sequence[int]],
    relations: Sequence[Sequence[str]],
) -> int | None:
    if not subset:
        return 1
    matrix = -sp.Matrix(integral_gram).extract(subset, subset)
    if any(matrix[:rank, :rank].det() <= 0 for rank in range(1, len(subset) + 1)):
        return None
    result = 1
    for component in _connected_components(subset, relations):
        result *= _finite_component_order(component, relations)
    return result


def _solomon_volume(
    gram: Sequence[Sequence[int]], relations: Sequence[Sequence[str]]
) -> dict[str, Any]:
    size = len(gram)
    euler = Fraction(0)
    rank_histogram: Counter[int] = Counter()
    order_histogram: Counter[int] = Counter()
    for mask in range(1 << size):
        subset = [node for node in range(size) if mask & (1 << node)]
        order = _finite_special_order(subset, gram, relations)
        if order is None:
            continue
        euler += Fraction((-1) ** len(subset), order)
        rank_histogram[len(subset)] += 1
        order_histogram[order] += 1
    volume_over_pi_squared = Fraction(4, 3) * euler
    return {
        "rational_euler_characteristic": str(euler),
        "volume_over_pi_squared": str(volume_over_pi_squared),
        "finite_special_subset_rank_histogram": {
            str(rank): count for rank, count in sorted(rank_histogram.items())
        },
        "finite_special_group_order_histogram": {
            str(order): count for order, count in sorted(order_histogram.items())
        },
        "formula": "Vol(P)/pi^2 = (4/3) chi(W) in dimension four",
    }


def _inner(left: Sequence[int], right: Sequence[int]) -> int:
    return left[0] * right[0] - sum(
        left[index] * right[index] for index in range(1, 5)
    )


def _reflection(root: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    square = _inner(root, root)
    if square not in {-1, -2}:
        raise ArithmeticError("expected a square-minus-one/two root")
    covector = [root[0], *(-entry for entry in root[1:])]
    coefficient = 2 // square
    return tuple(
        tuple(
            (1 if row == column else 0)
            - coefficient * root[row] * covector[column]
            for column in range(5)
        )
        for row in range(5)
    )


def _multiply_mod_two(
    left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(5))
            % 2
            for column in range(5)
        )
        for row in range(5)
    )


def _mod_two_image_order(roots: Sequence[Sequence[int]]) -> int:
    generators = [
        tuple(tuple(entry % 2 for entry in row) for row in _reflection(root))
        for root in roots
    ]
    identity = tuple(
        tuple(1 if row == column else 0 for column in range(5))
        for row in range(5)
    )
    seen = {identity}
    queue = deque([identity])
    while queue:
        element = queue.popleft()
        for generator in generators:
            product = _multiply_mod_two(element, generator)
            if product not in seen:
                seen.add(product)
                queue.append(product)
    return len(seen)


def _gram(roots: Sequence[Sequence[int]]) -> list[list[int]]:
    return [[_inner(left, right) for right in roots] for left in roots]


def _type_record(
    upstream: dict[str, Any], p4_volume: dict[str, Any]
) -> dict[str, Any]:
    relations = _relation_matrix_from_gram(upstream["canonical_gram"])
    volume = _solomon_volume(upstream["canonical_gram"], relations)
    p4_ratio = Fraction(volume["volume_over_pi_squared"]) / Fraction(
        p4_volume["volume_over_pi_squared"]
    )
    if p4_ratio.denominator != 1:
        raise ArithmeticError("nonintegral P4 covolume ratio")
    image_order = _mod_two_image_order(upstream["representative_walls"])
    expected_names = {
        "T6a": "W(A2+A1)",
        "T6b": "W(A2)",
        "T7": "W(A1+A1)",
        "T8": "W(A1)",
    }
    p4_image_order = 120
    congruence_index = p4_image_order // image_order
    return {
        "type_id": upstream["type_id"],
        "facet_count": upstream["facet_count"],
        "finite_vertex_count": upstream["cone_certificate"][
            "finite_vertex_count"
        ],
        "ideal_vertex_count": upstream["cone_certificate"]["ideal_vertex_count"],
        "finite_volume_noncompact": upstream["cone_certificate"]["decision"]
        == "finite_volume_noncompact",
        "wall_norm_profile": upstream["wall_norm_profile"],
        "coxeter_relation_profile": upstream["coxeter_relation_profile"],
        "relation_vector_in_compiler_order": list(_relation_vector(relations)),
        "canonical_unlabelled_relation_vector": list(
            _canonical_relation_vector(relations)
        ),
        "unlabelled_relation_orbit_size": math.factorial(len(relations))
        // upstream["weighted_gram_automorphism_order"],
        "volume": volume,
        "index_in_vinberg_P4": int(p4_ratio),
        "mod_two_image": {
            "order": image_order,
            "coxeter_type": expected_names[upstream["type_id"]],
            "index_in_W_A4_equals_S5": congruence_index,
        },
        "equals_full_preimage_of_mod_two_image": int(p4_ratio)
        == congruence_index,
        "arithmetic": True,
    }


def build_certificate(args: argparse.Namespace) -> dict[str, Any]:
    upstream = json.loads(COMPILER.read_text(encoding="utf-8"))
    by_id = {record["type_id"]: record for record in upstream["surviving_gram_types"]}
    if set(by_id) != {"T6a", "T6b", "T7", "T8"}:
        raise ArithmeticError("upstream four-type census changed")

    source_checks = {
        "compiler_sha256": _sha256(COMPILER),
        "coxiter_archive_sha256": _sha256(args.coxiter_archive),
        "tumarkin_archive_sha256": _sha256(args.tumarkin_archive),
        "tumarkin_tex_sha256": _sha256(args.tumarkin / "facets.tex"),
        "tumarkin_pyr4_eps_sha256": _sha256(args.tumarkin / "pyr4.eps"),
        "ma_zheng_archive_sha256": _sha256(args.ma_zheng_archive),
        "ma_zheng_tex_sha256": _sha256(args.ma_zheng),
        "im_hof_list_pdf_sha256": _sha256(args.im_hof),
        "july_2026_archive_sha256": _sha256(args.archive_2026),
    }
    expected_hashes = {
        "coxiter_archive_sha256": COXITER_ARCHIVE_SHA256,
        "tumarkin_archive_sha256": TUMARKIN_ARCHIVE_SHA256,
        "tumarkin_tex_sha256": TUMARKIN_TEX_SHA256,
        "tumarkin_pyr4_eps_sha256": TUMARKIN_PYR4_SHA256,
        "ma_zheng_archive_sha256": MA_ZHENG_ARCHIVE_SHA256,
        "ma_zheng_tex_sha256": MA_ZHENG_TEX_SHA256,
        "im_hof_list_pdf_sha256": IM_HOF_PDF_SHA256,
        "july_2026_archive_sha256": ARCHIVE_2026_SHA256,
    }
    if any(source_checks[key] != value for key, value in expected_hashes.items()):
        raise ArithmeticError("a pinned literature/corpus source hash changed")

    tumarkin_text = (args.tumarkin / "facets.tex").read_text(encoding="utf-8")
    normalized_tumarkin_text = re.sub(r"\s+", " ", tumarkin_text)
    required_tumarkin_phrases = [
        "this completes the classification of hyperbolic Coxeter $n$-polytopes",
        "are listed in Tables~\\ref{pyr3}--\\ref{pyr17}",
        "\\epsfig{file=pyr4.eps",
    ]
    if any(
        phrase not in normalized_tumarkin_text
        for phrase in required_tumarkin_phrases
    ):
        raise ArithmeticError("Tumarkin completeness/table statement changed")
    source_2026_text = args.tex_2026.read_text(encoding="utf-8")
    if "five-dimensional polytopes with eight facets" not in source_2026_text:
        raise ArithmeticError("July 2026 eight-facet classification scope changed")

    relations = {
        type_id: _relation_matrix_from_gram(record["canonical_gram"])
        for type_id, record in by_id.items()
    }
    coxiter = _scan_coxiter(args.coxiter, relations)
    if coxiter["ordinary_file_count"] != 1214:
        raise ArithmeticError("CoxIter ordinary corpus count changed")
    if coxiter["doubling_file_count"] != 2959:
        raise ArithmeticError("CoxIter doubling corpus count changed")
    if coxiter["doubling_archive_sha256"] != DOUBLING_SHA256:
        raise ArithmeticError("CoxIter doubling archive hash changed")

    for type_id, relative in TUMARKIN_MATCH_FILES.items():
        path = args.coxiter / relative
        if _sha256(path) != TUMARKIN_MATCH_SHA256[type_id]:
            raise ArithmeticError(f"{type_id} Tumarkin transcription changed")
        matches = coxiter["ordinary_exact_matches"][type_id]
        if [match["path"] for match in matches if "Tumarkin" in match["path"]] != [
            relative.as_posix()
        ]:
            raise ArithmeticError(f"{type_id} Tumarkin match not unique")

    im_hof_path = args.coxiter / IM_HOF_MATCH_FILE
    if _sha256(im_hof_path) != IM_HOF_MATCH_SHA256:
        raise ArithmeticError("Im Hof transcription changed")
    if [match["path"] for match in coxiter["ordinary_exact_matches"]["T7"]] != [
        IM_HOF_MATCH_FILE.as_posix()
    ]:
        raise ArithmeticError("T7 did not have the unique expected CoxIter match")
    if coxiter["ordinary_exact_matches"]["T8"]:
        raise ArithmeticError("T8 unexpectedly entered the pinned ordinary corpus")
    if any(coxiter["doubling_exact_matches"].values()):
        raise ArithmeticError("a candidate unexpectedly matched the doubling corpus")

    p4_path = args.coxiter / P4_FILE
    if _sha256(p4_path) != P4_FILE_SHA256:
        raise ArithmeticError("Vinberg P4 CoxIter record changed")
    p4_relations = _relation_matrix_from_gram(_gram(P4_ROOTS))
    _facets, _dimension, p4_corpus_relations = _parse_coxiter(p4_path)
    if not _exact_permutations(p4_relations, p4_corpus_relations):
        raise ArithmeticError("explicit Vinberg P4 roots do not match CoxIter")
    p4_volume = _solomon_volume(_gram(P4_ROOTS), p4_relations)
    if p4_volume["volume_over_pi_squared"] != "1/1440":
        raise ArithmeticError("unexpected P4 volume")
    if _mod_two_image_order(P4_ROOTS) != 120:
        raise ArithmeticError("unexpected P4 mod-two image")

    type_records = {
        type_id: _type_record(record, p4_volume)
        for type_id, record in sorted(by_id.items())
    }
    expected_indices = {"T6a": 10, "T6b": 20, "T7": 30, "T8": 60}
    if {
        type_id: record["index_in_vinberg_P4"]
        for type_id, record in type_records.items()
    } != expected_indices:
        raise ArithmeticError("unexpected congruence-cover indices")
    if not all(
        record["equals_full_preimage_of_mod_two_image"]
        for record in type_records.values()
    ):
        raise ArithmeticError("a type is not the full expected mod-two preimage")

    ma_zheng = _ma_zheng_audit(args.ma_zheng, relations["T7"])
    if type_records["T7"]["volume"]["volume_over_pi_squared"] != "1/48":
        raise ArithmeticError("Solomon volume does not match Ma--Zheng")

    lineage = {
        "T6a": {
            "abstract_polytope_status": "previously published",
            "earliest_certified_lineage": "Tumarkin n+2 classification (2003 preprint; Math. Notes 75 (2004))",
            "exact_source_family": "H4 pyramid family 07_02 at (k,l)=(2,4)",
            "coxiter_transcription": TUMARKIN_MATCH_FILES["T6a"].as_posix(),
            "complete_classification_applies": True,
            "project_novelty_remaining": "new full-nef triangular Looijenga carrier, not a new abstract Coxeter polytope",
        },
        "T6b": {
            "abstract_polytope_status": "previously published",
            "earliest_certified_lineage": "Tumarkin n+2 classification (2003 preprint; Math. Notes 75 (2004))",
            "exact_source_family": "H4 pyramid family 07_01 at (k,l)=(2,4)",
            "coxiter_transcription": TUMARKIN_MATCH_FILES["T6b"].as_posix(),
            "complete_classification_applies": True,
            "project_novelty_remaining": "new full-nef triangular Looijenga carrier, not a new abstract Coxeter polytope",
        },
        "T7": {
            "abstract_polytope_status": "previously published",
            "earliest_certified_lineage": "Im Hof, A class of hyperbolic Coxeter groups (1985), completed in Napier cycles and hyperbolic Coxeter groups (1990)",
            "exact_source_family": "doubly truncated H4 Napier cycle, n=4 row 3 column 3 (CoxIter item 13)",
            "coxiter_transcription": IM_HOF_MATCH_FILE.as_posix(),
            "independent_complete_census_match": "Ma--Zheng P_{2,1}",
            "complete_classification_applies": True,
            "project_novelty_remaining": "new full-nef triangular Looijenga carrier, not a new abstract Coxeter polytope",
        },
        "T8": {
            "abstract_polytope_status": "worldwide priority unresolved",
            "exact_group_theoretic_lineage": "index-60 full preimage of W(A1) under W(P4)->W(A4)=S5 modulo two",
            "pinned_coxiter_ordinary_match_count": 0,
            "pinned_coxiter_doubling_match_count": 0,
            "complete_noncompact_H4_eight_facet_classification_available": False,
            "worldwide_novelty_claimed": False,
            "project_novelty_remaining": "previously frozen full-nef carrier plus a new canonical congruence-preimage interpretation",
        },
    }

    payload: dict[str, Any] = {
        "schema": "coxeter.triangular_looijenga_h4_worldwide_lineage.v1",
        "kind": "exact_literature_and_congruence_lineage_certificate",
        "audit_date": "2026-08-08",
        "source_checks": source_checks,
        "sources": {
            "tumarkin": {
                "title": "Hyperbolic Coxeter n-polytopes with n+2 facets",
                "url": "https://arxiv.org/abs/math/0301133",
                "journal": "Mathematical Notes 75 (2004), 848--854",
                "scope": "complete finite-volume n+2-facet classification with Kaplinskaja and Esselmann",
            },
            "im_hof_1985": {
                "title": "A class of hyperbolic Coxeter groups",
                "journal": "Expositiones Mathematicae 3 (1985), 179--186",
            },
            "im_hof_1990": {
                "title": "Napier cycles and hyperbolic Coxeter groups",
                "journal": "Bulletin de la Societe Mathematique de Belgique, Serie A 42 (1990), 523--545",
                "list_url": "https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/poly/napier.pdf",
            },
            "ma_zheng": ma_zheng,
            "current_eight_facet_boundary": {
                "compact_H4_eight_facets": "complete (Ma--Zheng and Burcroff, 2022)",
                "finite_volume_H5_eight_facets": "complete (Guo--Ma--Zang--Zheng, arXiv:2607.14715, 2026)",
                "finite_volume_noncompact_H4_eight_facets": "no complete census located",
                "hcpdm_url": "https://github.com/GeoTopChristy/HCPdm",
                "july_2026_url": "https://arxiv.org/abs/2607.14715",
            },
        },
        "vinberg_P4": {
            "lattice": "I_(1,4)",
            "facet_count": 5,
            "explicit_roots": P4_ROOTS,
            "coxiter_record": P4_FILE.as_posix(),
            "relation_vector": list(_relation_vector(p4_relations)),
            "volume": p4_volume,
            "mod_two_image": {
                "coxeter_type": "W(A4)=S5",
                "order": 120,
            },
        },
        "types": type_records,
        "coxiter_bounded_audit": coxiter,
        "lineage_decisions": lineage,
        "main_conclusion": {
            "T6a_T6b_T7_new_abstract_polytopes": False,
            "T6a_T6b_T7_new_AG_carriers_in_this_project": True,
            "T8_worldwide_abstract_priority_resolved": False,
            "all_four_are_canonical_mod_two_congruence_preimages": True,
            "worldwide_claim_is_made_only_where_a_complete_classification_or_prior_exact_item_exists": True,
        },
        "all_checks_pass": True,
    }
    digest_payload = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["certificate_sha256"] = hashlib.sha256(
        digest_payload.encode("utf-8")
    ).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coxiter", type=Path, default=DEFAULT_COXITER)
    parser.add_argument(
        "--coxiter-archive", type=Path, default=DEFAULT_COXITER_ARCHIVE
    )
    parser.add_argument("--tumarkin", type=Path, default=DEFAULT_TUMARKIN)
    parser.add_argument(
        "--tumarkin-archive", type=Path, default=DEFAULT_TUMARKIN_ARCHIVE
    )
    parser.add_argument("--ma-zheng", type=Path, default=DEFAULT_MA_ZHENG)
    parser.add_argument(
        "--ma-zheng-archive", type=Path, default=DEFAULT_MA_ZHENG_ARCHIVE
    )
    parser.add_argument("--im-hof", type=Path, default=DEFAULT_IM_HOF)
    parser.add_argument("--archive-2026", type=Path, default=DEFAULT_2026_ARCHIVE)
    parser.add_argument("--tex-2026", type=Path, default=DEFAULT_2026_TEX)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = build_certificate(args)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload["main_conclusion"], indent=2, sort_keys=True))
    print(f"certificate_sha256={payload['certificate_sha256']}")


if __name__ == "__main__":
    main()
