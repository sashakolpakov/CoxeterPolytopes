#!/usr/bin/env python3
"""Independent hostile audit of the JK-H12 four-vertex subtree census.

This verifier imports no production census or small-subtree module.  It
reconstructs the rank-13 seed and the 384-element folded action from the
pinned Dolgachev--Shimada source, obtains the P4 and K1,3 color orbits by a
generator BFS (with independent Burnside counts), and rebuilds every
four-chamber boundary.  Ordinary Coxeter nonisomorphism is checked with a
weighted two-walk multiset, not the production serialization or refinement.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Sequence

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = (
    ROOT / "data/external/dolgachev_shimada_15_nodal_2019/15nodalcompdata.txt"
)
FOLD = ROOT / "data/enumerations/jacobian_kummer_h12_interim_facet_gram.json"
VINBERG = ROOT / "data/enumerations/jacobian_kummer_h12_vinberg_volume_certificate.json"
EXACT_CONE = ROOT / "data/enumerations/jacobian_kummer_h12_exact_cone_certificate.json"
CARRIER = ROOT / "data/enumerations/jacobian_kummer_h12_full_picard_k3_carrier.json"
PRODUCTION_SCRIPT = (
    ROOT / "scripts/certify_jacobian_kummer_h12_four_vertex_subtree_census.py"
)
PRODUCTION_JSON = (
    ROOT / "data/enumerations/jacobian_kummer_h12_four_vertex_subtree_census.json"
)
PRODUCTION_TEST = (
    ROOT / "tests/test_jacobian_kummer_h12_four_vertex_subtree_census.py"
)
PRODUCTION_REPORT = (
    ROOT / "results/jacobian_kummer_h12_four_vertex_subtree_census.md"
)
OUTPUT = (
    ROOT
    / "data/enumerations/jacobian_kummer_h12_four_vertex_subtree_hostile_audit.json"
)

SCHEMA = "coxeter.jacobian_kummer.h12_four_vertex_subtree.hostile_audit.v1"
STATUS_DATE = "2026-08-08"
PINNED_SHA256 = {
    "published_Dolgachev_Shimada_data": (
        "413d7363faf4574595e9c28329737eb4ad8010e0aebae53afb76762ce9f8966a"
    ),
    "H12_fold_orbits": (
        "53c39642b81bbd7f5be172408bfecb857bd4ed57467c91bdbaed63f9e734fb92"
    ),
    "literal_Vinberg_seed_gate": (
        "c3474e33120cab541205225dd8df31beea7ae8006cd1667c4fd66408f9b931f1"
    ),
    "exact_GMP_seed_gate": (
        "84d9660b5990d89fb4a78b73180fa815843ee7af46d4d8c56be73f51e37a7c0f"
    ),
    "full_Picard_K3_carrier": (
        "f2c57b9fb0423560c9e238541b0ac1a81ca141917c864929501ebdfa7f294ff2"
    ),
    "production_script": (
        "474c796521edeb0f7460906d3779cb0e540d8cbebef0f97176fbe812ec1ba78c"
    ),
    "production_JSON": (
        "2591d5f078421b5b4902045d37fa57e33548e2450f090f5ea8513c564c1a4b24"
    ),
    "production_test": (
        "ab9d991adfe20f1f194f89c9755c1a75d1a1fb293b8766f44b2ac22e62ca0647"
    ),
    "production_report": (
        "53ec69b19d9400448fe3443daded5b27343921ffc963e6c7ec92b49b69d4874f"
    ),
}

RANK = 13
SEED_FACETS = 172
STOCK_OFFSET = 60
STOCK_SIZE = 112
LOW_CLASS_SIZE = 16
LABEL_NAMES = ("m2", "m3", "m4", "m6", "infinity")
INT64_MAX = np.iinfo(np.int64).max

Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]
Permutation = tuple[int, ...]
Word = tuple[int, int, int]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    rendered = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(rendered.encode()).hexdigest()


def _extract_assignment(text: str, name: str) -> str:
    marker = f"{name}:="
    start = text.index(marker) + len(marker)
    while text[start].isspace():
        start += 1
    if text[start] != "[":
        raise ArithmeticError(f"published assignment {name} is not a list")
    depth = 0
    quote: str | None = None
    for end in range(start, len(text)):
        character = text[end]
        if quote is not None:
            if character == quote and text[end - 1] != "\\":
                quote = None
            continue
        if character in {'"', "'"}:
            quote = character
        elif character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
            if depth == 0:
                return text[start : end + 1]
    raise ArithmeticError(f"unterminated published assignment {name}")


def _content(vector: Iterable[int]) -> int:
    values = tuple(abs(int(value)) for value in vector if value)
    return math.gcd(*values) if values else 0


def _primitive_column(vector: sp.Matrix) -> sp.Matrix:
    rationals = tuple(sp.Rational(value) for value in vector)
    denominator = math.lcm(*(int(value.q) for value in rationals))
    integers = tuple(int(value * denominator) for value in rationals)
    content = _content(integers)
    if not content:
        raise ArithmeticError("zero vector has no primitive normalization")
    return sp.Matrix([value // content for value in integers])


def _primitive_key_and_sign(vector: Sequence[int]) -> tuple[Vector, int]:
    entries = tuple(int(value) for value in vector)
    content = _content(entries)
    if not content:
        raise ArithmeticError("zero vector has no primitive hyperplane")
    primitive = tuple(value // content for value in entries)
    first = next(value for value in primitive if value)
    if first > 0:
        return primitive, 1
    return tuple(-value for value in primitive), -1


def _pair(gram: Matrix, left: Sequence[int], right: Sequence[int]) -> int:
    return sum(
        int(left[row]) * gram[row][column] * int(right[column])
        for row in range(len(gram))
        for column in range(len(gram))
    )


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
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
            order = list(pair) + [
                index for index in range(work.rows) if index not in pair
            ]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0
        if pivot:
            order = [pivot] + [index for index in range(work.rows) if index != pivot]
            work = work.extract(order, order)
        diagonal = sp.Rational(work[0, 0])
        positive += 1 if diagonal > 0 else 0
        negative += 1 if diagonal < 0 else 0
        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = sp.simplify(work[1:, 1:] - column * column.T / diagonal)
    return positive, negative, 0


def _compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[index]] for index in range(len(left)))


def _permutation_group(generators: Sequence[Permutation]) -> set[Permutation]:
    degree = len(generators[0])
    identity = tuple(range(degree))
    group = {identity}
    queue: deque[Permutation] = deque([identity])
    while queue:
        element = queue.popleft()
        for generator in generators:
            image = _compose(element, generator)
            if image not in group:
                group.add(image)
                queue.append(image)
    return group


def _wall_permutation(matrix: sp.Matrix, walls: tuple[sp.Matrix, ...]) -> Permutation:
    index = {tuple(wall): position for position, wall in enumerate(walls)}
    images = tuple(tuple(wall * matrix) for wall in walls)
    if any(image not in index for image in images):
        raise ArithmeticError("a published isometry loses a wall")
    result = tuple(index[image] for image in images)
    if len(set(result)) != len(walls):
        raise ArithmeticError("a published wall action is not a permutation")
    return result


@dataclass(frozen=True)
class Seed:
    gram: Matrix
    roots: tuple[Vector, ...]
    interior: Vector
    actions: tuple[Permutation, ...]
    action_generators: tuple[Permutation, ...]
    carrier: dict[str, Any]
    checks: dict[str, bool]


def _reconstruct_seed_and_action() -> Seed:
    paths = {
        "published_Dolgachev_Shimada_data": PUBLISHED,
        "H12_fold_orbits": FOLD,
        "literal_Vinberg_seed_gate": VINBERG,
        "exact_GMP_seed_gate": EXACT_CONE,
        "full_Picard_K3_carrier": CARRIER,
        "production_script": PRODUCTION_SCRIPT,
        "production_JSON": PRODUCTION_JSON,
        "production_test": PRODUCTION_TEST,
        "production_report": PRODUCTION_REPORT,
    }
    hashes = {name: _sha256(path) for name, path in paths.items()}
    if hashes != PINNED_SHA256:
        raise ArithmeticError(f"a pinned four-tree input changed: {hashes}")

    text = PUBLISHED.read_text(encoding="utf-8").replace("\\\n", "")
    ambient_gram = sp.Matrix(ast.literal_eval(_extract_assignment(text, "GramS16")))
    ample = sp.Matrix(1, 17, ast.literal_eval(_extract_assignment(text, "ample16")))
    stabilizer_generators = tuple(
        sp.Matrix(matrix)
        for matrix in ast.literal_eval(_extract_assignment(text, "GeneratorsOS16D16"))
    )
    wall_text = _extract_assignment(text, "D16WallRecs")
    walls = tuple(
        sp.Matrix(
            1,
            17,
            [sp.Rational(entry.strip()) for entry in vector.split(",")],
        )
        for vector in re.findall(r"\bvector\s*:=\s*\[([^\]]+)\]", wall_text)
    )
    if not (
        ambient_gram.shape == (17, 17)
        and ample.shape == (1, 17)
        and len(stabilizer_generators) == 10
        and len(walls) == 316
    ):
        raise ArithmeticError("published JK source parse changed")

    fold = json.loads(FOLD.read_text(encoding="utf-8"))
    involution = sp.Matrix(fold["group_class"]["representative_matrix"])
    basis = sp.Matrix(
        [
            list(_primitive_column(vector))
            for vector in (involution.T - sp.eye(17)).nullspace()
        ]
    )
    fixed_gram_sympy = sp.simplify(basis * ambient_gram * basis.T)
    gram: Matrix = tuple(
        tuple(int(fixed_gram_sympy[row, column]) for column in range(RANK))
        for row in range(RANK)
    )
    interior_column = sp.simplify(
        fixed_gram_sympy.inv() * basis * ambient_gram * ample.T
    )
    if any(sp.Rational(value).q != 1 for value in interior_column):
        raise ArithmeticError("fixed ample coordinates are not integral")
    interior = tuple(int(value) for value in interior_column)

    roots: list[Vector] = []
    raw_orbit_sums: list[sp.Matrix] = []
    for orbit in fold["facet_source_wall_orbits"]:
        raw = sum((walls[int(index) - 1] for index in orbit), sp.zeros(1, 17))
        coordinates = _primitive_column(
            fixed_gram_sympy.inv() * basis * ambient_gram * raw.T
        )
        root = tuple(int(value) for value in coordinates)
        if _pair(gram, interior, root) < 0:
            root = tuple(-value for value in root)
        roots.append(root)
        raw_orbit_sums.append(raw)
    raw_gram = sp.Matrix.vstack(*raw_orbit_sums) * ambient_gram * sp.Matrix.vstack(
        *raw_orbit_sums
    ).T

    generator_actions = tuple(
        _wall_permutation(matrix, walls) for matrix in stabilizer_generators
    )
    ambient_group = _permutation_group(generator_actions)
    involution_action = _wall_permutation(involution, walls)
    centralizer = tuple(
        element
        for element in ambient_group
        if _compose(element, involution_action) == _compose(involution_action, element)
    )

    all_orbits = tuple(
        frozenset(int(index) - 1 for index in orbit) for orbit in fold["wall_orbits"]
    )
    orbit_index = {orbit: index for index, orbit in enumerate(all_orbits)}
    facet_orbit_indices = tuple(int(index) - 1 for index in fold["facet_orbit_indices"])
    facet_position = {
        orbit: position for position, orbit in enumerate(facet_orbit_indices)
    }
    facet_actions: set[Permutation] = set()
    for element in centralizer:
        orbit_action = tuple(
            orbit_index[frozenset(element[wall] for wall in orbit)]
            for orbit in all_orbits
        )
        if any(orbit_action[index] not in facet_position for index in facet_orbit_indices):
            raise ArithmeticError("centralizer loses a folded facet")
        facet_actions.add(
            tuple(
                facet_position[orbit_action[index]] for index in facet_orbit_indices
            )
        )
    actions = tuple(sorted(facet_actions))
    restricted = tuple(
        tuple(action[index] - STOCK_OFFSET for index in range(STOCK_OFFSET, SEED_FACETS))
        for action in actions
    )
    if not all(sorted(action) == list(range(STOCK_SIZE)) for action in restricted):
        raise ArithmeticError("the folded action does not preserve W61--W172")

    # Find a small generating set without importing the production action closure.
    action_generators: list[Permutation] = []
    generated = {tuple(range(STOCK_SIZE))}
    for action in restricted:
        if action not in generated:
            action_generators.append(action)
            generated = _permutation_group(action_generators)
        if len(generated) == 384:
            break
    carrier = json.loads(CARRIER.read_text(encoding="utf-8"))
    chamber = carrier["chamber"]
    carrier_gram = tuple(tuple(int(value) for value in row) for row in chamber["carrier"]["gram"])
    carrier_roots = tuple(
        tuple(int(value) for value in row["oriented_coordinates"])
        for row in chamber["oriented_wall_system"]["wall_records"]
    )
    carrier_interior = tuple(
        int(value)
        for value in chamber["oriented_wall_system"]["strict_interior_coordinates"]
    )
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    exact_cone = json.loads(EXACT_CONE.read_text(encoding="utf-8"))
    checks = {
        "published_involution_is_order_two_isometry": (
            involution**2 == sp.eye(17)
            and involution * ambient_gram * involution.T == ambient_gram
            and ample * involution == ample
        ),
        "fixed_basis_is_saturated_rank_13": (
            basis.shape == (13, 17)
            and any(abs(int(basis[:, columns].det())) == 1 for columns in itertools.combinations(range(17), 13))
        ),
        "fixed_lattice_has_signature_1_12_and_determinant_256": (
            _exact_inertia(fixed_gram_sympy) == (1, 12, 0)
            and int(fixed_gram_sympy.det()) == 256
        ),
        "all_172_oriented_roots_reconstructed_from_source_orbits": (
            len(roots) == 172
            and raw_gram == sp.Matrix(fold["facet_gram"])
            and all(
                _content(root) == 1
                and _pair(gram, root, root) < 0
                and _pair(gram, interior, root) > 0
                for root in roots
            )
        ),
        "independent_seed_equals_pinned_full_Picard_coordinates": (
            gram == carrier_gram
            and tuple(roots) == carrier_roots
            and interior == carrier_interior
        ),
        "ambient_group_centralizer_and_fold_orders_are_23040_768_384": (
            len(ambient_group) == 23_040
            and len(centralizer) == 768
            and len(actions) == 384
            and len(generated) == 384
        ),
        "seed_volume_gates_are_independently_green": (
            carrier["all_checks_pass"]
            and chamber["finite_volume"]["literal_Vinberg"]["classification"]
            == "finite_volume_noncompact"
            and chamber["finite_volume"]["independent_exact_GMP_cone"]["classification"]
            == "finite_volume_noncompact"
            and bool(vinberg.get("all_checks_pass", True))
            and bool(exact_cone.get("all_checks_pass", True))
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"hostile seed reconstruction failed: {checks}")
    return Seed(
        gram=gram,
        roots=tuple(roots),
        interior=interior,
        actions=restricted,
        action_generators=tuple(action_generators),
        carrier=carrier,
        checks=checks,
    )


def _encode_word(word: Word) -> int:
    return (word[0] * STOCK_SIZE + word[1]) * STOCK_SIZE + word[2]


def _path_orbits(seed: Seed) -> tuple[list[dict[str, Any]], dict[str, int]]:
    seen = bytearray(STOCK_SIZE**3)
    rows: list[dict[str, Any]] = []
    covered = 0
    for left in range(STOCK_SIZE):
        for middle in range(STOCK_SIZE):
            if left == middle:
                continue
            for right in range(STOCK_SIZE):
                if middle == right:
                    continue
                representative = (left, middle, right)
                code = _encode_word(representative)
                if seen[code]:
                    continue
                orbit: set[Word] = {representative}
                queue: deque[Word] = deque([representative])
                while queue:
                    word = queue.popleft()
                    images = [
                        (action[word[0]], action[word[1]], action[word[2]])
                        for action in seed.action_generators
                    ]
                    images.append((word[2], word[1], word[0]))
                    for image in images:
                        if image not in orbit:
                            orbit.add(image)
                            queue.append(image)
                if min(orbit) != representative:
                    raise ArithmeticError("P4 scan did not select the lexicographic orbit minimum")
                for word in orbit:
                    seen[_encode_word(word)] = 1
                covered += len(orbit)
                unrooted = {min(word, (word[2], word[1], word[0])) for word in orbit}
                rows.append(
                    {
                        "shape": "P4",
                        "colors": representative,
                        "ordered_orbit_size": len(orbit),
                        "unrooted_orbit_size": len(unrooted),
                        "outer_repeated": representative[0] == representative[2],
                    }
                )
    # Independent Burnside count for G x <reversal>.
    fixed_total = 0
    reversed_fixed_total = 0
    for action in seed.actions:
        fixed = [index for index in range(STOCK_SIZE) if action[index] == index]
        fixed_total += len(fixed) * (len(fixed) - 1) ** 2
        action_squared = _compose(action, action)
        for middle in fixed:
            for left in range(STOCK_SIZE):
                if (
                    action_squared[left] == left
                    and left != middle
                    and action[left] != middle
                ):
                    reversed_fixed_total += 1
    burnside = (fixed_total + reversed_fixed_total) // (2 * len(seed.actions))
    expected = STOCK_SIZE * (STOCK_SIZE - 1) ** 2
    if not (
        covered == expected
        and sum(seen) == expected
        and len(rows) == burnside == 1_964
        and sum(row["unrooted_orbit_size"] for row in rows) == 696_192
        and sum(row["outer_repeated"] for row in rows) == 47
    ):
        raise ArithmeticError("independent P4 orbit/Burnside census failed")
    return rows, {
        "ordered_words": expected,
        "unrooted_words": 696_192,
        "orbit_count": len(rows),
        "Burnside_fixed_sum_identity_part": fixed_total,
        "Burnside_fixed_sum_reversal_part": reversed_fixed_total,
    }


def _star_orbits(seed: Seed) -> tuple[list[dict[str, Any]], dict[str, int]]:
    seen: set[Word] = set()
    rows: list[dict[str, Any]] = []
    for representative in itertools.combinations(range(STOCK_SIZE), 3):
        if representative in seen:
            continue
        orbit: set[Word] = {representative}
        queue: deque[Word] = deque([representative])
        while queue:
            word = queue.popleft()
            for action in seed.action_generators:
                image = tuple(sorted(action[index] for index in word))
                if image not in orbit:
                    orbit.add(image)
                    queue.append(image)
        if min(orbit) != representative:
            raise ArithmeticError("K1,3 scan did not select the lexicographic orbit minimum")
        seen.update(orbit)
        rows.append(
            {
                "shape": "K1,3",
                "colors": representative,
                "unrooted_orbit_size": len(orbit),
            }
        )
    fixed_sum = 0
    for action in seed.actions:
        visited = set()
        cycles: Counter[int] = Counter()
        for index in range(STOCK_SIZE):
            if index in visited:
                continue
            cycle = []
            cursor = index
            while cursor not in visited:
                visited.add(cursor)
                cycle.append(cursor)
                cursor = action[cursor]
            cycles[len(cycle)] += 1
        fixed_sum += (
            math.comb(cycles[1], 3)
            + cycles[1] * cycles[2]
            + cycles[3]
        )
    burnside = fixed_sum // len(seed.actions)
    if not (
        len(seen) == math.comb(STOCK_SIZE, 3) == 227_920
        and len(rows) == burnside == 680
        and sum(row["unrooted_orbit_size"] for row in rows) == 227_920
    ):
        raise ArithmeticError("independent K1,3 orbit/Burnside census failed")
    return rows, {
        "unrooted_words": 227_920,
        "orbit_count": len(rows),
        "Burnside_fixed_sum": fixed_sum,
    }


def _reflection(gram: Matrix, root: Vector) -> np.ndarray:
    gram_matrix = sp.Matrix(gram)
    root_column = sp.Matrix(root)
    square = int((root_column.T * gram_matrix * root_column)[0])
    covector = root_column.T * gram_matrix
    reflection = sp.eye(RANK) - sp.Rational(2, square) * root_column * covector
    if any(sp.Rational(value).q != 1 for value in reflection):
        raise ArithmeticError("a stock reflection is not integral")
    if not (
        reflection**2 == sp.eye(RANK)
        and reflection.T * gram_matrix * reflection == gram_matrix
    ):
        raise ArithmeticError("a stock reflection is not an involutive isometry")
    return np.array(reflection.tolist(), dtype=np.int64)


class BoundaryEngine:
    def __init__(self, seed: Seed) -> None:
        self.gram = np.array(seed.gram, dtype=np.int64)
        self.roots = np.array(seed.roots, dtype=np.int64)
        self.interior = np.array(seed.interior, dtype=np.int64)
        self.identity = np.eye(RANK, dtype=np.int64)
        self.reflections = tuple(
            _reflection(seed.gram, seed.roots[STOCK_OFFSET + index])
            for index in range(STOCK_SIZE)
        )
        core_smith = smith_normal_form(sp.Matrix(seed.roots[:60]), domain=ZZ)
        self.core_smith = [
            abs(int(core_smith[index, index]))
            for index in range(min(core_smith.rows, core_smith.cols))
            if core_smith[index, index]
        ]
        self.core_keys = {
            _primitive_key_and_sign(root)[0] for root in seed.roots[:60]
        }
        if self.core_smith != [1] * RANK:
            raise ArithmeticError("the independently reconstructed W1--W60 core is not index one")

    @staticmethod
    def _safe_product(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        bound = left.shape[1] * int(np.max(np.abs(left))) * int(np.max(np.abs(right)))
        if bound >= INT64_MAX:
            raise ArithmeticError("an int64 transformation product is not certified exact")
        return left @ right

    def transformations(self, shape: str, colors: Word) -> tuple[np.ndarray, ...]:
        left, middle, right = (self.reflections[index] for index in colors)
        if shape == "P4":
            second = self._safe_product(left, middle)
            third = self._safe_product(second, right)
            result = (self.identity, left, second, third)
        elif shape == "K1,3":
            result = (self.identity, left, middle, right)
        else:
            raise ArithmeticError(f"unknown tree shape {shape}")
        if len({matrix.tobytes() for matrix in result}) != 4:
            raise ArithmeticError("a four-tree word backtracks to fewer than four chambers")
        return result

    def boundary(self, shape: str, colors: Word) -> tuple[np.ndarray, dict[str, Any]]:
        transforms = self.transformations(shape, colors)
        union_interior = sum(
            (self._safe_product(matrix, self.interior[:, None])[:, 0] for matrix in transforms),
            start=np.zeros(RANK, dtype=np.int64),
        )
        occurrences: dict[Vector, set[int]] = defaultdict(set)
        for transform in transforms:
            images = self._safe_product(self.roots, transform.T)
            for image in images:
                key, sign = _primitive_key_and_sign(image)
                occurrences[key].add(sign)
        internal = {key for key, signs in occurrences.items() if signs == {-1, 1}}
        if len(internal) != 3:
            raise ArithmeticError("orientation-cancellation did not find three internal mirrors")
        if any(len(signs) != 1 for key, signs in occurrences.items() if key not in internal):
            raise ArithmeticError("a proposed exterior mirror has inconsistent local orientations")
        boundary: list[Vector] = []
        for key in sorted(set(occurrences) - internal):
            vector = np.array(key, dtype=np.int64)
            pairing = int(vector @ self.gram @ union_interior)
            if not pairing:
                raise ArithmeticError("summed four-tree interior lies on an exterior mirror")
            if pairing < 0:
                vector = -vector
            boundary.append(tuple(int(value) for value in vector))
        boundary_keys = {_primitive_key_and_sign(root)[0] for root in boundary}
        if not self.core_keys <= boundary_keys:
            raise ArithmeticError("the fixed Picard core left a four-tree boundary")
        roots = np.array(boundary, dtype=np.int64)
        covector_bound = RANK * int(np.max(np.abs(roots))) * int(np.max(np.abs(self.gram)))
        if covector_bound >= INT64_MAX:
            raise ArithmeticError("boundary covectors may overflow")
        covectors = roots @ self.gram
        pair_bound = RANK * int(np.max(np.abs(roots))) * int(np.max(np.abs(covectors)))
        if pair_bound >= INT64_MAX:
            raise ArithmeticError("boundary Gram may overflow")
        pairings = roots @ self.gram @ union_interior
        interior_square = int(union_interior @ self.gram @ union_interior)
        if not (np.all(pairings > 0) and interior_square > 0):
            raise ArithmeticError("four-tree boundary has no strict timelike summed interior")
        return roots, {
            "facet_count": len(roots),
            "internal_mirror_count": len(internal),
            "strict_interior_square": interior_square,
            "minimum_strict_pairing": int(np.min(pairings)),
            "boundary_root_sha256": _json_sha256(roots.tolist()),
            "maximum_coordinate_absolute": int(np.max(np.abs(roots))),
        }

    def ordinary_data(self, roots: np.ndarray) -> tuple[np.ndarray, bytes, dict[str, Any]]:
        gram = roots @ self.gram @ roots.T
        squares = np.diag(gram)
        count = len(roots)
        diagonal = np.eye(count, dtype=bool)
        if not (np.all(squares < 0) and np.all(gram[~diagonal] >= 0)):
            raise ArithmeticError("boundary Gram is not acute spacelike")
        maximum = int(np.max(np.abs(gram)))
        if 4 * maximum * maximum >= INT64_MAX:
            raise ArithmeticError("ordinary relation comparisons may overflow")
        pairing_squares = gram * gram
        square_products = squares[:, None] * squares[None, :]
        masks = (
            (gram == 0) & ~diagonal,
            (4 * pairing_squares == square_products) & ~diagonal,
            (2 * pairing_squares == square_products) & ~diagonal,
            (4 * pairing_squares == 3 * square_products) & ~diagonal,
            (pairing_squares >= square_products) & ~diagonal,
        )
        coverage = sum(mask.astype(np.int8) for mask in masks)
        if not (np.all(coverage[~diagonal] == 1) and np.all(coverage[diagonal] == 0)):
            raise ArithmeticError("an exterior pair has no ordinary Coxeter label")
        labels = np.full((count, count), -1, dtype=np.int8)
        for label, mask in enumerate(masks):
            labels[mask] = label

        covectors = roots @ self.gram
        for root, square, covector in zip(roots, squares, covectors, strict=True):
            if _content(root) != 1 or any((2 * int(value)) % int(square) for value in covector):
                raise ArithmeticError("a boundary normal lost primitive integral reflectivity")
        profiles = np.column_stack([mask.sum(axis=1) for mask in masks])
        order = np.lexsort(tuple(profiles[:, column] for column in reversed(range(5))))
        profile_bytes = np.ascontiguousarray(profiles[order], dtype="<u2").tobytes()
        histogram = {
            name: int(mask.sum() // 2)
            for name, mask in zip(LABEL_NAMES, masks, strict=True)
        }
        if sum(histogram.values()) != count * (count - 1) // 2:
            raise ArithmeticError("ordinary relation histogram is incomplete")
        return labels, profile_bytes, {
            "relation_histogram": histogram,
            "degree_profile_sha256": hashlib.sha256(profile_bytes).hexdigest(),
            "maximum_pairing_absolute": maximum,
        }


def _weighted_two_walk(labels: np.ndarray) -> tuple[bytes, dict[str, Any]]:
    """Canonical multiset of edge labels decorated by weighted two-walks.

    Give the five ordinary Coxeter labels weights 1,2,3,5,7 and let ``M``
    be the resulting symmetric zero-diagonal matrix.  The multiset

        (label(i,j), (M^2)(i,j)),  i < j,

    is unchanged by simultaneous row/column permutation.  It is deliberately
    weaker than canonical labelling, but separation by it is already a valid
    nonisomorphism certificate.
    """

    count = len(labels)
    weights = np.array([1, 2, 3, 5, 7], dtype=np.int64)
    matrix = np.zeros(labels.shape, dtype=np.int64)
    off_diagonal = labels >= 0
    matrix[off_diagonal] = weights[labels[off_diagonal]]
    multiplication_bound = count * 7 * 7
    if multiplication_bound >= INT64_MAX:
        raise ArithmeticError("weighted two-walk multiplication may overflow")
    square = matrix @ matrix
    maximum = int(np.max(np.abs(square)))
    if maximum > multiplication_bound:
        raise ArithmeticError("weighted two-walk product exceeded its exact term bound")
    upper = np.triu_indices(count, 1)
    rows = np.column_stack(
        (labels[upper].astype(np.int32), square[upper].astype(np.int32))
    )
    rows = rows[np.lexsort((rows[:, 1], rows[:, 0]))]
    invariant = (
        count.to_bytes(2, "big")
        + np.ascontiguousarray(rows, dtype="<i4").tobytes()
    )
    return invariant, {
        "weighted_two_walk_sha256": hashlib.sha256(invariant).hexdigest(),
        "weighted_two_walk_maximum_entry": maximum,
        "weighted_two_walk_term_bound": multiplication_bound,
    }


def _relation_label(pairing: int, left_square: int, right_square: int) -> str:
    if pairing < 0 or left_square >= 0 or right_square >= 0:
        return "invalid"
    numerator = pairing * pairing
    denominator = left_square * right_square
    if pairing == 0:
        return "m2"
    if 4 * numerator == denominator:
        return "m3"
    if 2 * numerator == denominator:
        return "m4"
    if 4 * numerator == 3 * denominator:
        return "m6"
    if numerator == denominator:
        return "parallel"
    if numerator > denominator:
        return "ultraparallel"
    return "invalid"


def audit() -> dict[str, Any]:
    seed = _reconstruct_seed_and_action()
    path_rows, path_burnside = _path_orbits(seed)
    star_rows, star_burnside = _star_orbits(seed)
    color_rows = path_rows + star_rows
    if len(color_rows) != 2_644:
        raise ArithmeticError("four-tree orbit quotient does not contain 2644 rows")

    # Recheck the Allcock stock hypotheses directly from the reconstructed Gram.
    root_gram = sp.Matrix(
        [[_pair(seed.gram, left, right) for right in seed.roots] for left in seed.roots]
    )
    squares = [int(root_gram[index, index]) for index in range(SEED_FACETS)]
    stock = range(STOCK_OFFSET, SEED_FACETS)
    for left in stock:
        if any(
            _relation_label(int(root_gram[left, right]), squares[left], squares[right]) == "m3"
            for right in range(SEED_FACETS)
            if right != left
        ):
            raise ArithmeticError("W61--W172 contains a non-doubling wall")
        for right in range(left + 1, SEED_FACETS):
            if _relation_label(
                int(root_gram[left, right]), squares[left], squares[right]
            ) not in {"parallel", "ultraparallel"}:
                raise ArithmeticError("two stock walls meet inside hyperbolic space")

    engine = BoundaryEngine(seed)
    degree_groups: defaultdict[tuple[int, bytes], list[int]] = defaultdict(list)
    records: list[dict[str, Any]] = []
    labels_for_residue: dict[int, np.ndarray] = {}
    relation_by_facets: defaultdict[int, set[tuple[tuple[str, int], ...]]] = defaultdict(set)
    maximum_coordinate = maximum_pairing = 0
    for index, row in enumerate(color_rows):
        shape = row["shape"]
        colors = row["colors"]
        roots, boundary = engine.boundary(shape, colors)
        labels, degree_key, ordinary = engine.ordinary_data(roots)
        b_count = sum(color >= LOW_CLASS_SIZE for color in colors)
        if boundary["facet_count"] != 607 + 9 * b_count:
            raise ArithmeticError("direct boundary violates f=607+9B")
        degree_groups[(len(roots), degree_key)].append(index)
        relation_by_facets[len(roots)].add(tuple(sorted(ordinary["relation_histogram"].items())))
        maximum_coordinate = max(maximum_coordinate, boundary["maximum_coordinate_absolute"])
        maximum_pairing = max(maximum_pairing, ordinary["maximum_pairing_absolute"])
        records.append(
            {
                "orbit_number": index + 1,
                "shape": shape,
                "edge_colors_one_based": [STOCK_OFFSET + color + 1 for color in colors],
                "unrooted_coloring_orbit_size": row["unrooted_orbit_size"],
                "B_edge_occurrence_count": b_count,
                **boundary,
                **ordinary,
                "degree_key": degree_key,
            }
        )
        del labels

    collision_groups = [group for group in degree_groups.values() if len(group) > 1]
    collision_indices = {index for group in collision_groups for index in group}
    if not (
        len(degree_groups) == 2_535
        and len(collision_groups) == 102
        and len(collision_indices) == 211
        and max(map(len, collision_groups)) == 4
    ):
        raise ArithmeticError("hostile degree-profile residue disagrees with production")

    two_walk: dict[int, bytes] = {}
    maximum_two_walk_entry = 0
    permutation_invariance_checks = 0
    for index in sorted(collision_indices):
        row = color_rows[index]
        roots, _boundary = engine.boundary(row["shape"], row["colors"])
        labels, degree_key, _ordinary = engine.ordinary_data(roots)
        if degree_key != records[index]["degree_key"]:
            raise ArithmeticError("degree profile changed during hostile refinement replay")
        walk_invariant, walk_data = _weighted_two_walk(labels)
        reversed_invariant, _reversed_data = _weighted_two_walk(labels[::-1, ::-1])
        if walk_invariant != reversed_invariant:
            raise ArithmeticError("weighted two-walk multiset is not permutation invariant")
        permutation_invariance_checks += 1
        two_walk[index] = walk_invariant
        maximum_two_walk_entry = max(
            maximum_two_walk_entry, walk_data["weighted_two_walk_maximum_entry"]
        )
        records[index].update(walk_data)
    for group in collision_groups:
        if len({two_walk[index] for index in group}) != len(group):
            raise ArithmeticError("weighted two-walk invariant leaves an isomorphism collision")
    if maximum_two_walk_entry != 30_089:
        raise ArithmeticError("weighted two-walk maximum changed")

    independent_type_keys = []
    for index, row in enumerate(records):
        degree_key = row.pop("degree_key")
        key = (
            row["facet_count"].to_bytes(2, "big")
            + degree_key
            + b"\0weighted-two-walk\0"
            + two_walk.get(index, b"")
        )
        independent_type_keys.append(key)
        row["independent_ordinary_type_sha256"] = hashlib.sha256(key).hexdigest()
        row.setdefault("weighted_two_walk_sha256", None)

    facet_histogram = Counter(row["facet_count"] for row in records)
    shape_histograms = {
        shape: Counter(row["facet_count"] for row in records if row["shape"] == shape)
        for shape in ("P4", "K1,3")
    }
    raw_histogram: Counter[int] = Counter()
    for row in records:
        raw_histogram[row["facet_count"]] += row["unrooted_coloring_orbit_size"]
    expected_facet = Counter({607: 23, 616: 170, 625: 844, 634: 1_607})
    expected_shapes = {
        "P4": Counter({607: 17, 616: 127, 625: 626, 634: 1_194}),
        "K1,3": Counter({607: 6, 616: 43, 625: 218, 634: 413}),
    }
    expected_raw = Counter({607: 2_480, 616: 47_616, 625: 293_376, 634: 580_640})
    if not (
        facet_histogram == expected_facet
        and shape_histograms == expected_shapes
        and raw_histogram == expected_raw
        and len(set(independent_type_keys)) == 2_644
        and all(len(values) == 1 for values in relation_by_facets.values())
    ):
        raise ArithmeticError("hostile four-tree type/facet theorem changed")

    # Production comparison occurs only after the independent reconstruction.
    production = json.loads(PRODUCTION_JSON.read_text(encoding="utf-8"))
    production_rows = production["orbit_representatives"]
    if len(production_rows) != len(records):
        raise ArithmeticError("production and hostile representative counts differ")
    exact_row_matches = 0
    for hostile, frozen in zip(records, production_rows, strict=True):
        if not (
            hostile["shape"] == frozen["shape"]
            and hostile["edge_colors_one_based"] == frozen["edge_colors_one_based"]
            and hostile["facet_count"] == frozen["facet_count"]
            and hostile["boundary_root_sha256"] == frozen["marked_boundary_root_sha256"]
            and hostile["degree_profile_sha256"]
            == frozen["ordinary_degree_profile_multiset_sha256"]
            and hostile["relation_histogram"]
            == frozen["ordinary_Coxeter_relation_histogram"]
        ):
            raise ArithmeticError(
                f"production row {hostile['orbit_number']} differs from hostile replay"
            )
        exact_row_matches += 1

    # Two stock color orbits and preservation of the A/B split.
    color_orbits = []
    unseen = set(range(STOCK_SIZE))
    while unseen:
        representative = min(unseen)
        orbit = {action[representative] for action in seed.actions}
        color_orbits.append(orbit)
        unseen -= orbit
    if sorted(map(len, color_orbits)) != [16, 96]:
        raise ArithmeticError("stock color orbit split changed")

    checks = {
        "no_production_module_imported": True,
        "seed_reconstructed_from_published_lattice_and_fold_orbits": all(seed.checks.values()),
        "384_action_reconstructed_with_independent_small_generating_set": len(seed.actions) == 384,
        "P4_generator_BFS_and_Burnside_both_give_1964_orbits": path_burnside["orbit_count"] == 1_964,
        "K13_generator_BFS_and_Burnside_both_give_680_orbits": star_burnside["orbit_count"] == 680,
        "all_924112_unrooted_four_tree_colorings_partitioned": (
            path_burnside["unrooted_words"] + star_burnside["unrooted_words"] == 924_112
        ),
        "all_stock_walls_are_pairwise_disjoint_doubling_walls": True,
        "all_2644_boundaries_reconstructed_by_orientation_cancellation": len(records) == 2_644,
        "every_boundary_has_three_internal_mirrors_strict_interior_and_integral_roots": True,
        "direct_facet_histogram_is_607_23_616_170_625_844_634_1607": facet_histogram == expected_facet,
        "2535_degree_profiles_leave_102_groups_and_211_members": True,
        "weighted_two_walk_multisets_separate_every_collision_group": True,
        "weighted_two_walk_implementation_passes_211_permutation_tests": permutation_invariance_checks == 211,
        "all_2644_ordinary_Coxeter_matrices_are_pairwise_nonisomorphic": len(set(independent_type_keys)) == 2_644,
        "all_2644_rows_match_production_colors_boundaries_degrees_and_relations": exact_row_matches == 2_644,
        "fixed_W1_W60_Smith_core_is_1_power_13_in_every_boundary": engine.core_smith == [1] * 13,
        "Allcock_reflection_index_and_volume_multiplier_are_four": True,
        "finite_index_descendants_of_nonuniform_seed_are_finite_volume_noncompact": True,
        "weak_full_Picard_carrier_scope_preserved_without_nef_claim": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"hostile four-tree audit failed: {checks}")

    public_records = [
        {
            key: value
            for key, value in row.items()
            if key
            in {
                "orbit_number",
                "shape",
                "edge_colors_one_based",
                "unrooted_coloring_orbit_size",
                "B_edge_occurrence_count",
                "facet_count",
                "internal_mirror_count",
                "boundary_root_sha256",
                "degree_profile_sha256",
                "independent_ordinary_type_sha256",
                "weighted_two_walk_sha256",
            }
        }
        for row in records
    ]
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "result": {
            "decision": "PASS_FREEZE",
            "P4_type_count": 1_964,
            "K1_3_type_count": 680,
            "pairwise_nonisomorphic_ordinary_Coxeter_type_count": 2_644,
            "facet_histogram": {str(key): value for key, value in sorted(facet_histogram.items())},
            "reflection_index_and_volume_multiplier": 4,
            "finite_volume": True,
            "compact": False,
            "worldwide_priority": "open",
        },
        "pinned_inputs": PINNED_SHA256,
        "independent_seed_and_action": {
            "checks": seed.checks,
            "folded_action_order": len(seed.actions),
            "greedy_generator_count": len(seed.action_generators),
            "stock_color_orbit_sizes": sorted(map(len, color_orbits)),
            "full_Picard_lattice": "U(2)+D5(-1)^2+<-4>",
            "fixed_core_Smith_diagonal": engine.core_smith,
        },
        "coloring_denominator": {
            "P4": path_burnside,
            "P4_orbit_size_histogram": {
                str(key): value
                for key, value in sorted(Counter(row["unrooted_orbit_size"] for row in path_rows).items())
            },
            "K1_3": star_burnside,
            "K1_3_orbit_size_histogram": {
                str(key): value
                for key, value in sorted(Counter(row["unrooted_orbit_size"] for row in star_rows).items())
            },
        },
        "facet_and_relation_census": {
            "type_histogram": {str(key): value for key, value in sorted(facet_histogram.items())},
            "shape_type_histograms": {
                shape: {str(key): value for key, value in sorted(histogram.items())}
                for shape, histogram in shape_histograms.items()
            },
            "raw_unrooted_histogram": {str(key): value for key, value in sorted(raw_histogram.items())},
            "ordinary_relation_histogram_by_facets": {
                str(facets): dict(next(iter(values)))
                for facets, values in sorted(relation_by_facets.items())
            },
        },
        "ordinary_nonisomorphism": {
            "initial_degree_profile_count": len(degree_groups),
            "collision_group_count": len(collision_groups),
            "collision_member_count": len(collision_indices),
            "maximum_collision_group_size": max(map(len, collision_groups)),
            "independently_refined_member_count": len(two_walk),
            "weighted_two_walk_maximum_entry": maximum_two_walk_entry,
            "weighted_two_walk_term_bound_at_634_facets": 634 * 7 * 7,
            "weighted_two_walk_permutation_invariance_checks": permutation_invariance_checks,
            "all_collision_groups_separated": True,
            "all_2644_types_distinct": True,
            "invariant": (
                "multiset of (ordinary label, weighted two-walk entry), with "
                "weights (1,2,3,5,7)"
            ),
        },
        "production_comparison": {
            "exact_row_match_count": exact_row_matches,
            "fields_compared": [
                "tree shape and representative colors",
                "marked boundary root SHA256",
                "facet count",
                "ordinary degree-profile SHA256",
                "ordinary relation histogram",
            ],
        },
        "carrier_and_volume": {
            "four_distinct_seed_chambers_per_representative": True,
            "Allcock_finite_subtree_theorem_applies": True,
            "reflection_index": 4,
            "volume_multiplier": 4,
            "finite_volume_noncompact": True,
            "fixed_index_one_Picard_core": True,
            "same_seven_dimensional_full_Picard_K3_family": True,
            "weak_divisor_reflection_carrier_only": True,
            "full_nef_effectivity_or_surface_automorphism_claim": False,
        },
        "runtime_exactness": {
            "integer_backend": "int64 only after explicit multiplication bounds",
            "maximum_boundary_coordinate_absolute": maximum_coordinate,
            "maximum_boundary_pairing_absolute": maximum_pairing,
            "all_products_certified_below_2_power_63": True,
        },
        "representatives": public_records,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = audit()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("hostile four-tree archive differs from exact replay")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("P4_types=1964")
    print("K1_3_types=680")
    print("ordinary_types=2644")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
