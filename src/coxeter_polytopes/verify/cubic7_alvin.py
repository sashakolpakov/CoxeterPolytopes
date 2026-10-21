"""Independent exact replay of conductor-seven AlVin chamber outputs.

The input consists of three logically separate objects:

* an AlVin/CoxIter ``.coxiter`` diagram;
* the root vectors printed by AlVin; and
* the diagonal form over ``O_K = Z[theta]``, where
  ``theta^3 + theta^2 - 2*theta - 1 = 0``.

AlVin writes a conductor-seven integer as ``[u,v,w]`` in the trigonometric
basis ``(lambda_1,lambda_2,lambda_3)`` with
``lambda_j=2*cos(2*j*pi/7)``.  This module immediately converts that
dependent-looking presentation to the power basis ``(1,theta,theta^2)``:

``[u,v,w] = (w-2*v) + (u-w)*theta + (v-w)*theta^2``.

No root enumeration is performed here.  In particular, a successful replay
certifies the supplied Coxeter polyhedron and its reflections, but does *not*
independently certify that AlVin found every reflecting wall of the lattice.
That separation is recorded explicitly in :class:`Cubic7AlVinVerification`.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Mapping, Sequence

from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_ONE,
    CUBIC7_ZERO,
    LABEL14_PRODUCT,
    LABEL7_PRODUCT,
    Cubic7,
    cubic7_inertia,
)
from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.search.cubic7_seed_manifest import alvin_rc7_to_cubic7


Vector = tuple[Cubic7, ...]
Matrix = tuple[tuple[Cubic7, ...], ...]


class Cubic7AlVinError(ValueError):
    """A malformed input or a failed exact verification condition."""


def _element_key(value: Cubic7) -> tuple[tuple[int, int], ...]:
    return tuple(
        (coefficient.numerator, coefficient.denominator)
        for coefficient in (value.a, value.b, value.c)
    )


def _element_json(value: Cubic7) -> dict[str, str]:
    return value.to_json()


def _matrix_json(matrix: Sequence[Sequence[Cubic7]]) -> list[list[dict[str, str]]]:
    return [[_element_json(entry) for entry in row] for row in matrix]


def _vector_json(vector: Sequence[Cubic7]) -> list[dict[str, str]]:
    return [_element_json(entry) for entry in vector]


def _inertia_json(value: Inertia) -> dict[str, int]:
    return {
        "positive": value.positive,
        "negative": value.negative,
        "zero": value.zero,
    }


def alvin_lambda_integer(u: int, v: int, w: int) -> Cubic7:
    """Convert AlVin's ``[u,v,w]`` to the integral power basis.

    Since ``lambda_2=theta^2-2`` and
    ``lambda_3=1-theta-theta^2``, the displayed formula follows.  Notice in
    particular that AlVin's rational integer ``n`` is ``[-n,-n,-n]``.
    """

    return alvin_rc7_to_cubic7((u, v, w))


def _split_top_level(text: str, delimiter: str = ",") -> tuple[str, ...]:
    pieces: list[str] = []
    start = 0
    depths = {"[": 0, "(": 0, "{": 0}
    closing = {
        "]": "[",
        ")": "(",
        "}": "{",
    }
    for index, character in enumerate(text):
        if character in depths:
            depths[character] += 1
        elif character in closing:
            opener = closing[character]
            depths[opener] -= 1
            if depths[opener] < 0:
                raise Cubic7AlVinError(f"unmatched {character!r} in {text!r}")
        elif character == delimiter and all(depth == 0 for depth in depths.values()):
            pieces.append(text[start:index].strip())
            start = index + 1
    if any(depth != 0 for depth in depths.values()):
        raise Cubic7AlVinError(f"unmatched bracket in {text!r}")
    pieces.append(text[start:].strip())
    return tuple(piece for piece in pieces if piece)


_INTEGER_RE = re.compile(r"^[+-]?\d+$")
_FRACTION_RE = re.compile(r"^[+-]?\d+/[1-9]\d*$")
_ALVIN_TRIPLE_RE = re.compile(r"^\[\s*([+-]?\d+)\s*,\s*([+-]?\d+)\s*,\s*([+-]?\d+)\s*\]$")
_PARI_TRIPLE_RE = re.compile(
    r"^RCI7\(\s*([+-]?\d+)\s*,\s*([+-]?\d+)\s*,\s*([+-]?\d+)\s*\)$",
    re.IGNORECASE,
)


def parse_alvin_element(text: str) -> Cubic7:
    """Parse one integer in generic or PARI AlVin conductor-seven syntax."""

    value = text.strip()
    if _INTEGER_RE.fullmatch(value) or _FRACTION_RE.fullmatch(value):
        return Cubic7(Fraction(value))
    match = _ALVIN_TRIPLE_RE.fullmatch(value) or _PARI_TRIPLE_RE.fullmatch(value)
    if match is None:
        raise Cubic7AlVinError(f"unsupported AlVin conductor-seven integer {text!r}")
    return alvin_lambda_integer(*(int(coefficient) for coefficient in match.groups()))


def parse_diagonal_form(text: str) -> Vector:
    """Parse an actual diagonal form such as ``[-1,0,0],1,1,1,1``.

    The signs are the signs of the quadratic form itself.  Thus the standard
    compact four-dimensional control is ``diag(-theta,1,1,1,1)`` and is
    entered exactly as in the example above.
    """

    stripped = text.strip()
    if stripped.startswith("<") and stripped.endswith(">"):
        stripped = stripped[1:-1]
    values = tuple(parse_alvin_element(token) for token in _split_top_level(stripped))
    if len(values) < 3:
        raise Cubic7AlVinError("a hyperbolic diagonal form needs at least three entries")
    if not all(value.is_integral for value in values):
        raise Cubic7AlVinError("the diagonal lattice form must be integral over Z[theta]")
    return values


def _extract_outer_vector(line: str) -> str | None:
    """Extract the coordinate body from one generic or PARI vector line."""

    generic = re.search(r"\be_?\{?(\d+)\}?\s*=\s*([({])", line)
    if generic is not None:
        opener = generic.group(2)
        closer = ")" if opener == "(" else "}"
        start = generic.end()
        end = line.rfind(closer)
        if end < start:
            raise Cubic7AlVinError(f"unterminated root vector line: {line!r}")
        return line[start:end]

    if "vectors" not in line or "=" not in line:
        return None
    # In PARI output every vector is the first square-bracket expression;
    # RCI7 coordinates use parentheses and therefore do not interfere.
    start = line.find("[")
    if start < 0:
        return None
    depth = 0
    for index in range(start, len(line)):
        if line[index] == "[":
            depth += 1
        elif line[index] == "]":
            depth -= 1
            if depth == 0:
                return line[start + 1 : index]
    raise Cubic7AlVinError(f"unterminated PARI root vector line: {line!r}")


def _root_index(line: str, fallback: int) -> int:
    generic = re.search(r"\be_?\{?(\d+)\}?\s*=", line)
    return int(generic.group(1)) if generic is not None else fallback


def _parse_json_element(value: Any, basis: str) -> Cubic7:
    if isinstance(value, Mapping):
        return Cubic7.from_json(value)
    if isinstance(value, int) or isinstance(value, str):
        return parse_alvin_element(str(value))
    if isinstance(value, list) and len(value) == 3:
        if basis == "alvin_lambda":
            if not all(isinstance(entry, int) for entry in value):
                raise Cubic7AlVinError("AlVin lambda coordinates must be integer triples")
            return alvin_lambda_integer(*value)
        if basis == "power":
            return Cubic7(*(Fraction(entry) for entry in value))
    raise Cubic7AlVinError(f"invalid serialized root coordinate {value!r}")


def parse_alvin_roots_text(text: str) -> tuple[Vector, ...]:
    """Parse roots from generic/PARI AlVin output or a small JSON envelope.

    JSON input is useful for archival replay and has the form
    ``{"basis":"power"|"alvin_lambda", "roots":[...]}``.  A bare JSON
    list is treated as power-basis data.  Generic AlVin output is preferred
    for primary capture because it preserves the program transcript.
    """

    stripped = text.lstrip()
    if stripped.startswith("{") or stripped.startswith("["):
        try:
            payload = json.loads(text)
        except json.JSONDecodeError:
            payload = None
        if payload is not None:
            if isinstance(payload, Mapping):
                basis = str(payload.get("basis", "power"))
                raw_roots = payload.get("roots")
            else:
                basis = "power"
                raw_roots = payload
            if basis not in {"power", "alvin_lambda"}:
                raise Cubic7AlVinError("JSON root basis must be 'power' or 'alvin_lambda'")
            if not isinstance(raw_roots, list):
                raise Cubic7AlVinError("JSON roots must be an array")
            roots = tuple(
                tuple(_parse_json_element(entry, basis) for entry in root)
                for root in raw_roots
                if isinstance(root, list)
            )
            if len(roots) != len(raw_roots):
                raise Cubic7AlVinError("every JSON root must be an array")
            return roots

    indexed: dict[int, Vector] = {}
    next_index = 1
    for line in text.splitlines():
        body = _extract_outer_vector(line)
        if body is None:
            continue
        index = _root_index(line, next_index)
        vector = tuple(parse_alvin_element(token) for token in _split_top_level(body))
        if index in indexed:
            # AlVin can print its initial roots once, not twice.  Treat an
            # identical duplicate defensively but reject conflicting captures.
            if indexed[index] != vector:
                raise Cubic7AlVinError(f"conflicting definitions of root {index}")
        else:
            indexed[index] = vector
        next_index = max(next_index, index + 1)
    if not indexed:
        raise Cubic7AlVinError("no AlVin root vectors were found")
    expected = set(range(1, max(indexed) + 1))
    if set(indexed) != expected:
        raise Cubic7AlVinError(
            f"root indices are not contiguous: found {sorted(indexed)}, expected {sorted(expected)}"
        )
    return tuple(indexed[index] for index in sorted(indexed))


def parse_alvin_roots(path: str | Path) -> tuple[Vector, ...]:
    return parse_alvin_roots_text(Path(path).read_text(encoding="utf-8"))


@dataclass(frozen=True, slots=True)
class CoxIterDiagram:
    facet_count: int
    dimension: int
    labels: tuple[tuple[int, ...], ...]
    vertex_names: tuple[str, ...]


def parse_coxiter_text(text: str) -> CoxIterDiagram:
    """Parse the sparse three-column CoxIter graph format emitted by AlVin."""

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise Cubic7AlVinError("empty .coxiter input")
    header = re.fullmatch(r"(\d+)\s+(\d+)", lines[0])
    if header is None:
        raise Cubic7AlVinError("the first .coxiter line must be 'FACETS DIMENSION'")
    facet_count, dimension = (int(value) for value in header.groups())
    if facet_count < dimension + 1:
        raise Cubic7AlVinError("a full-dimensional chamber needs at least dimension+1 facets")

    cursor = 1
    names = tuple(str(index + 1) for index in range(facet_count))
    if cursor < len(lines) and lines[cursor].lower().startswith("vertices labels:"):
        names = tuple(lines[cursor].split(":", 1)[1].split())
        if len(names) != facet_count or len(set(names)) != facet_count:
            raise Cubic7AlVinError("invalid 'vertices labels' line")
        cursor += 1
    name_to_index = {name: index for index, name in enumerate(names)}
    labels = [[2 for _ in range(facet_count)] for _ in range(facet_count)]
    for index in range(facet_count):
        labels[index][index] = 1
    seen: set[tuple[int, int]] = set()
    for line in lines[cursor:]:
        fields = line.split()
        if len(fields) < 3:
            raise Cubic7AlVinError(f"invalid .coxiter edge line {line!r}")
        left_name, right_name = fields[:2]
        if left_name not in name_to_index or right_name not in name_to_index:
            raise Cubic7AlVinError(f"unknown vertex in .coxiter edge line {line!r}")
        try:
            label = int(fields[2])
        except ValueError as exc:
            raise Cubic7AlVinError(f"nonintegral Coxeter label in {line!r}") from exc
        left, right = name_to_index[left_name], name_to_index[right_name]
        if left == right:
            raise Cubic7AlVinError("a .coxiter edge cannot be a loop")
        pair = tuple(sorted((left, right)))
        if pair in seen:
            raise Cubic7AlVinError(f"duplicate .coxiter edge {pair}")
        if label < 0 or label == 2:
            raise Cubic7AlVinError(
                "stored .coxiter edges must have label 0, 1, or at least 3"
            )
        seen.add(pair)
        labels[left][right] = labels[right][left] = label
    return CoxIterDiagram(
        facet_count,
        dimension,
        tuple(tuple(row) for row in labels),
        names,
    )


def parse_coxiter(path: str | Path) -> CoxIterDiagram:
    return parse_coxiter_text(Path(path).read_text(encoding="utf-8"))


def _diagonal_matrix(diagonal: Sequence[Cubic7]) -> Matrix:
    return tuple(
        tuple(diagonal[i] if i == j else CUBIC7_ZERO for j in range(len(diagonal)))
        for i in range(len(diagonal))
    )


def _pair(diagonal: Sequence[Cubic7], left: Sequence[Cubic7], right: Sequence[Cubic7]) -> Cubic7:
    return sum(
        (diagonal[index] * left[index] * right[index] for index in range(len(diagonal))),
        start=CUBIC7_ZERO,
    )


def _wall_gram(diagonal: Sequence[Cubic7], roots: Sequence[Vector]) -> Matrix:
    return tuple(
        tuple(_pair(diagonal, left, right) for right in roots)
        for left in roots
    )


def _rank(rows: Sequence[Sequence[Cubic7]], column_count: int | None = None) -> int:
    if not rows:
        return 0
    width = len(rows[0]) if column_count is None else column_count
    if any(len(row) != width for row in rows):
        raise Cubic7AlVinError("row-rank input has inconsistent widths")
    work = [list(row) for row in rows]
    pivot_row = 0
    for column in range(width):
        pivot = next((row for row in range(pivot_row, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        for row in range(pivot_row + 1, len(work)):
            if not work[row][column]:
                continue
            multiplier = work[row][column] / pivot_value
            work[row] = [
                work[row][index] - multiplier * work[pivot_row][index]
                for index in range(width)
            ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def _one_dimensional_kernel(rows: Sequence[Sequence[Cubic7]]) -> Vector | None:
    """Return a projective kernel generator exactly when the nullity is one."""

    if not rows:
        return None
    width = len(rows[0])
    work = [list(row) for row in rows]
    if any(len(row) != width for row in work):
        raise Cubic7AlVinError("kernel input has inconsistent widths")
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(width):
        pivot = next((row for row in range(pivot_row, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                work[row][index] - multiplier * work[pivot_row][index]
                for index in range(width)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(work):
            break
    if width - pivot_row != 1:
        return None
    free = next(column for column in range(width) if column not in pivot_columns)
    kernel = [CUBIC7_ZERO for _ in range(width)]
    kernel[free] = CUBIC7_ONE
    for row, column in enumerate(pivot_columns):
        kernel[column] = -work[row][free]
    return tuple(kernel)


def _projective_key(vector: Sequence[Cubic7]) -> tuple[tuple[tuple[int, int], ...], ...]:
    pivot = next((entry for entry in vector if entry), None)
    if pivot is None:
        raise Cubic7AlVinError("zero vector does not define a projective ray")
    return tuple(_element_key(entry / pivot) for entry in vector)


def _submatrix(matrix: Matrix, indices: Sequence[int]) -> Matrix:
    return tuple(tuple(matrix[i][j] for j in indices) for i in indices)


def _connected_components(matrix: Matrix) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(len(matrix)))
    components: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        frontier = [start]
        while frontier:
            vertex = frontier.pop()
            neighbours = {other for other in unseen if matrix[vertex][other]}
            unseen.difference_update(neighbours)
            component.update(neighbours)
            frontier.extend(neighbours)
        components.append(tuple(sorted(component)))
    return tuple(components)


def _is_elliptic(matrix: Matrix) -> bool:
    return (
        all(matrix[i][i].sign(0) > 0 for i in range(len(matrix)))
        and all(matrix[i][j].sign(0) <= 0 for i in range(len(matrix)) for j in range(i))
        and cubic7_inertia(matrix, embedding=0) == Inertia(len(matrix), 0, 0)
    )


def _is_affine_component(matrix: Matrix) -> bool:
    size = len(matrix)
    if size < 2 or len(_connected_components(matrix)) != 1:
        return False
    if cubic7_inertia(matrix, embedding=0) != Inertia(size - 1, 0, 1):
        return False
    return all(
        _is_elliptic(
            _submatrix(matrix, tuple(index for index in range(size) if index != omitted))
        )
        for omitted in range(size)
    )


def _is_full_rank_parabolic(matrix: Matrix, dimension: int) -> bool:
    components = _connected_components(matrix)
    return bool(components) and all(
        _is_affine_component(_submatrix(matrix, component)) for component in components
    ) and sum(len(component) - 1 for component in components) == dimension - 1


@dataclass(frozen=True, slots=True)
class ExtremeRay:
    vector: Vector
    active_facets: tuple[int, ...]
    norm: Cubic7
    kind: str
    link_verified: bool

    def to_json(self) -> dict[str, Any]:
        return {
            "vector": _vector_json(self.vector),
            "active_facets": list(self.active_facets),
            "norm": _element_json(self.norm),
            "kind": self.kind,
            "link_verified": self.link_verified,
        }


def _enumerate_extreme_rays(
    diagonal: Sequence[Cubic7], roots: Sequence[Vector], wall_gram: Matrix
) -> tuple[ExtremeRay, ...]:
    ambient_rank = len(diagonal)
    dimension = ambient_rank - 1
    functionals = tuple(
        tuple(diagonal[index] * root[index] for index in range(ambient_rank))
        for root in roots
    )
    rays: dict[tuple[tuple[tuple[int, int], ...], ...], ExtremeRay] = {}
    for support in combinations(range(len(roots)), dimension):
        rows = tuple(functionals[index] for index in support)
        kernel = _one_dimensional_kernel(rows)
        if kernel is None:
            continue
        pairings = tuple(
            sum((row[index] * kernel[index] for index in range(ambient_rank)), start=CUBIC7_ZERO)
            for row in functionals
        )
        if all(value.sign(0) <= 0 for value in pairings):
            oriented = kernel
            oriented_pairings = pairings
        elif all(value.sign(0) >= 0 for value in pairings):
            oriented = tuple(-entry for entry in kernel)
            oriented_pairings = tuple(-entry for entry in pairings)
        else:
            continue
        active = tuple(index for index, value in enumerate(oriented_pairings) if not value)
        if _rank(tuple(functionals[index] for index in active), ambient_rank) != dimension:
            continue
        norm = _pair(diagonal, oriented, oriented)
        sign = norm.sign(0)
        if sign < 0:
            kind = "finite"
            active_gram = _submatrix(wall_gram, active)
            link_verified = cubic7_inertia(active_gram, embedding=0) == Inertia(
                dimension, 0, len(active) - dimension
            )
        elif sign == 0:
            kind = "cusp"
            link_verified = _is_full_rank_parabolic(_submatrix(wall_gram, active), dimension)
        else:
            kind = "hyperideal"
            link_verified = True
        ray = ExtremeRay(oriented, active, norm, kind, link_verified)
        rays.setdefault(_projective_key(oriented), ray)
    return tuple(
        sorted(
            rays.values(),
            key=lambda ray: (_projective_key(ray.vector), ray.active_facets),
        )
    )


_FINITE_LABEL_PRODUCTS = {
    3: Cubic7(1),
    4: Cubic7(2),
    6: Cubic7(3),
    7: LABEL7_PRODUCT,
    14: LABEL14_PRODUCT,
}


def _verify_diagram_labels(diagram: CoxIterDiagram, pair_products: Matrix) -> None:
    for i in range(diagram.facet_count):
        for j in range(i):
            label = diagram.labels[i][j]
            product = pair_products[i][j]
            if label == 2:
                valid = not product
            elif label == 0:
                valid = product == 4
            elif label == 1:
                valid = (product - 4).sign(0) > 0
            elif label in _FINITE_LABEL_PRODUCTS:
                valid = product == _FINITE_LABEL_PRODUCTS[label]
            else:
                raise Cubic7AlVinError(
                    f"finite Coxeter label {label} is not supported over the conductor-seven field"
                )
            if not valid:
                raise Cubic7AlVinError(
                    f".coxiter label mismatch on facets {i + 1},{j + 1}: "
                    f"label {label}, exact Cartan product {product}"
                )


def _fundamental_cycle_products(cartan: Matrix) -> tuple[Cubic7, ...]:
    """A field-generating set of Vinberg cyclic products.

    Length-two products are included for every edge.  One oriented product
    is then included for each non-tree edge of a spanning forest.  Together
    these generate the same field as all closed cyclic products.
    """

    order = len(cartan)
    adjacency = [set() for _ in range(order)]
    edges: list[tuple[int, int]] = []
    products: list[Cubic7] = []
    for i in range(order):
        for j in range(i):
            if cartan[i][j]:
                adjacency[i].add(j)
                adjacency[j].add(i)
                edges.append((j, i))
                products.append(cartan[i][j] * cartan[j][i])

    parent: dict[int, int | None] = {}
    tree_edges: set[tuple[int, int]] = set()
    for start in range(order):
        if start in parent:
            continue
        parent[start] = None
        frontier = [start]
        while frontier:
            vertex = frontier.pop()
            for neighbour in sorted(adjacency[vertex]):
                if neighbour in parent:
                    continue
                parent[neighbour] = vertex
                tree_edges.add(tuple(sorted((vertex, neighbour))))
                frontier.append(neighbour)

    def tree_path(left: int, right: int) -> list[int]:
        left_chain: list[int] = []
        cursor: int | None = left
        while cursor is not None:
            left_chain.append(cursor)
            cursor = parent[cursor]
        right_chain: list[int] = []
        cursor = right
        while cursor is not None:
            right_chain.append(cursor)
            cursor = parent[cursor]
        left_positions = {vertex: index for index, vertex in enumerate(left_chain)}
        common = next(vertex for vertex in right_chain if vertex in left_positions)
        return left_chain[: left_positions[common] + 1] + list(
            reversed(right_chain[: right_chain.index(common)])
        )

    for left, right in edges:
        if (left, right) in tree_edges:
            continue
        path = tree_path(left, right)
        product = CUBIC7_ONE
        for source, target in zip(path, path[1:]):
            product *= cartan[source][target]
        product *= cartan[right][left]
        products.append(product)
    return tuple(products)


# graph6 encoding is small enough to keep this verifier independent of a
# Python graph package.  nauty's labelg is used only for exact canonical
# labelling, not for any arithmetic or geometric assertion.
def _graph6_order(order: int) -> bytes:
    if order <= 62:
        return bytes((order + 63,))
    if order <= 258047:
        return bytes((126, *(((order >> shift) & 63) + 63 for shift in (12, 6, 0))))
    if order <= 68719476735:
        return bytes(
            (
                126,
                126,
                *(((order >> shift) & 63) + 63 for shift in (30, 24, 18, 12, 6, 0)),
            )
        )
    raise Cubic7AlVinError("graph6 cannot encode this many vertices")


def _encode_graph6(adjacency: Sequence[set[int]]) -> bytes:
    order = len(adjacency)
    bits = [1 if i in adjacency[j] else 0 for j in range(1, order) for i in range(j)]
    bits.extend(0 for _ in range((-len(bits)) % 6))
    body = bytes(
        63 + sum(bits[start + offset] << (5 - offset) for offset in range(6))
        for start in range(0, len(bits), 6)
    )
    return _graph6_order(order) + body + b"\n"


def _decode_graph6(data: bytes) -> tuple[set[int], ...]:
    line = data.strip()
    if not line:
        raise Cubic7AlVinError("labelg returned empty output")
    cursor = 0
    if line[0] != 126:
        order = line[0] - 63
        cursor = 1
    elif len(line) > 1 and line[1] != 126:
        order = sum((line[index] - 63) << shift for index, shift in zip(range(1, 4), (12, 6, 0)))
        cursor = 4
    else:
        order = sum(
            (line[index] - 63) << shift
            for index, shift in zip(range(2, 8), (30, 24, 18, 12, 6, 0))
        )
        cursor = 8
    bits: list[int] = []
    for character in line[cursor:]:
        value = character - 63
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))
    adjacency = [set() for _ in range(order)]
    bit = 0
    for right in range(1, order):
        for left in range(right):
            if bits[bit]:
                adjacency[left].add(right)
                adjacency[right].add(left)
            bit += 1
    return tuple(adjacency)


def _partition_symbols(count: int) -> str:
    # labelg gives '^' repetition syntax a special meaning.  Exclude it, NUL,
    # whitespace, '-' and the nonportable high-bit range.
    symbols = "".join(
        chr(value)
        for value in range(33, 127)
        if chr(value) not in {"^", "-"}
    )
    if count > len(symbols):
        raise Cubic7AlVinError("too many exact edge colours for labelg partition syntax")
    return symbols[:count]


def canonical_pair_product_gram(
    pair_products: Matrix, *, labelg: str | Path | None = None
) -> tuple[Matrix, str]:
    """Canonically permute a symmetric Cartan-product Gram with nauty.

    Each nonzero coloured edge is subdivided by a vertex of that edge colour.
    The original facet vertices and every edge-colour class form separate
    initial partition cells.  Canonical labelling of this ordinary coloured
    incidence graph is therefore exactly canonical labelling of the complete
    edge-coloured Gram matrix.
    """

    order = len(pair_products)
    if any(len(row) != order for row in pair_products):
        raise Cubic7AlVinError("pair-product Gram must be square")
    if any(pair_products[i][j] != pair_products[j][i] for i in range(order) for j in range(i)):
        raise Cubic7AlVinError("pair-product Gram must be symmetric")
    colours = sorted(
        {pair_products[i][j] for i in range(order) for j in range(i) if pair_products[i][j]},
        key=_element_key,
    )
    colour_index = {value: index for index, value in enumerate(colours)}
    grouped_edges: list[list[tuple[int, int]]] = [[] for _ in colours]
    for right in range(order):
        for left in range(right):
            value = pair_products[left][right]
            if value:
                grouped_edges[colour_index[value]].append((left, right))
    transformed_order = order + sum(len(edges) for edges in grouped_edges)
    adjacency = [set() for _ in range(transformed_order)]
    cursor = order
    for edges in grouped_edges:
        for left, right in edges:
            adjacency[cursor].update((left, right))
            adjacency[left].add(cursor)
            adjacency[right].add(cursor)
            cursor += 1

    executable = str(labelg) if labelg is not None else shutil.which("labelg")
    if executable is None:
        raise Cubic7AlVinError(
            "exact permutation canonicalization requires nauty's 'labelg' executable"
        )
    symbols = _partition_symbols(len(colours) + 1)
    partition = symbols[0] * order + "".join(
        symbols[index + 1] * len(edges) for index, edges in enumerate(grouped_edges)
    )
    completed = subprocess.run(
        [executable, "-q", f"-f{partition}"],
        input=_encode_graph6(adjacency),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise Cubic7AlVinError(
            "labelg canonicalization failed: "
            + completed.stderr.decode("utf-8", errors="replace").strip()
        )
    canonical_graph = _decode_graph6(completed.stdout)
    if len(canonical_graph) != transformed_order:
        raise Cubic7AlVinError("labelg changed the transformed graph order")

    canonical = [[CUBIC7_ZERO for _ in range(order)] for _ in range(order)]
    for index in range(order):
        canonical[index][index] = Cubic7(4)
    cursor = order
    for colour, edges in zip(colours, grouped_edges, strict=True):
        for _ in edges:
            neighbours = sorted(vertex for vertex in canonical_graph[cursor] if vertex < order)
            if len(neighbours) != 2:
                raise Cubic7AlVinError("canonical incidence graph has a malformed edge vertex")
            left, right = neighbours
            if canonical[left][right]:
                raise Cubic7AlVinError("canonical incidence graph repeats a coloured edge")
            canonical[left][right] = canonical[right][left] = colour
            cursor += 1
    result = tuple(tuple(row) for row in canonical)
    payload = json.dumps(
        {
            "dimensionless_pair_product_gram": _matrix_json(result),
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return result, sha256(payload).hexdigest()


@dataclass(frozen=True, slots=True)
class Cubic7AlVinVerification:
    dimension: int
    diagonal_form: Vector
    roots: tuple[Vector, ...]
    lattice_inertias: tuple[Inertia, Inertia, Inertia]
    wall_gram: Matrix
    wall_inertias: tuple[Inertia, Inertia, Inertia]
    cartan_matrix: Matrix
    pair_product_gram: Matrix
    cyclic_field_generators: tuple[Cubic7, ...]
    ground_field: str
    extreme_rays: tuple[ExtremeRay, ...]
    finite_volume: bool
    compact: bool
    geometry_status: str
    canonical_pair_product_gram: Matrix
    canonical_fingerprint: str
    alvin_declared_termination: bool | None

    @property
    def finite_vertices(self) -> tuple[ExtremeRay, ...]:
        return tuple(ray for ray in self.extreme_rays if ray.kind == "finite")

    @property
    def cusps(self) -> tuple[ExtremeRay, ...]:
        return tuple(ray for ray in self.extreme_rays if ray.kind == "cusp")

    @property
    def hyperideal_rays(self) -> tuple[ExtremeRay, ...]:
        return tuple(ray for ray in self.extreme_rays if ray.kind == "hyperideal")

    def to_json(self) -> dict[str, Any]:
        return {
            "schema": "cubic7-alvin-independent-replay-0.1.0",
            "field": {
                "generator": "theta=2*cos(2*pi/7)",
                "minimal_polynomial": "theta^3+theta^2-2*theta-1",
                "ring_of_integers": "Z[theta]",
                "ground_field_from_cyclic_products": self.ground_field,
            },
            "dimension": self.dimension,
            "facet_count": len(self.roots),
            "diagonal_form": _vector_json(self.diagonal_form),
            "roots_power_basis": [_vector_json(root) for root in self.roots],
            "lattice_inertias": [_inertia_json(value) for value in self.lattice_inertias],
            "wall_gram": _matrix_json(self.wall_gram),
            "wall_inertias": [_inertia_json(value) for value in self.wall_inertias],
            "normalized_cartan_gram": _matrix_json(self.cartan_matrix),
            "scale_invariant_pair_products": _matrix_json(self.pair_product_gram),
            "vinberg_arithmeticity": {
                "crystallographic_reflections_integral": True,
                "cyclic_product_generators": [
                    _element_json(value) for value in self.cyclic_field_generators
                ],
                "all_cyclic_products_integral": True,
                "nonidentity_conjugate_grams_positive_semidefinite": True,
                "passes": True,
            },
            "geometry": {
                "status": self.geometry_status,
                "finite_volume": self.finite_volume,
                "compact": self.compact,
                "finite_vertex_count": len(self.finite_vertices),
                "cusp_count": len(self.cusps),
                "hyperideal_ray_count": len(self.hyperideal_rays),
                "extreme_rays": [ray.to_json() for ray in self.extreme_rays],
            },
            "deduplication": {
                "equivalence": (
                    "simultaneous facet permutation of the scale-invariant "
                    "pair-product Gram"
                ),
                "canonicalizer": "nauty labelg on an exactly edge-coloured incidence graph",
                "canonical_pair_product_gram": _matrix_json(
                    self.canonical_pair_product_gram
                ),
                "sha256": self.canonical_fingerprint,
            },
            "scope": {
                "alvin_transcript_declares_termination": self.alvin_declared_termination,
                "alvin_wall_enumeration_completeness_independently_checked": False,
                "warning": (
                    "This replay verifies the supplied chamber exactly. It does not replay "
                    "AlVin's height-ordered root enumeration or prove that no reflecting wall "
                    "was omitted."
                ),
            },
        }


def verify_cubic7_alvin_candidate(
    diagram: CoxIterDiagram,
    diagonal_form: Sequence[Cubic7],
    roots: Sequence[Sequence[Cubic7]],
    *,
    transcript_text: str | None = None,
    labelg: str | Path | None = None,
) -> Cubic7AlVinVerification:
    """Replay all exact checks available from the supplied chamber data."""

    diagonal = tuple(Cubic7.coerce(value) for value in diagonal_form)
    root_vectors = tuple(tuple(Cubic7.coerce(entry) for entry in root) for root in roots)
    ambient_rank = len(diagonal)
    dimension = ambient_rank - 1
    if diagram.dimension != dimension:
        raise Cubic7AlVinError(
            f".coxiter dimension {diagram.dimension} does not match form dimension {dimension}"
        )
    if diagram.facet_count != len(root_vectors):
        raise Cubic7AlVinError(
            f".coxiter has {diagram.facet_count} facets but {len(root_vectors)} roots were supplied"
        )
    if any(len(root) != ambient_rank for root in root_vectors):
        raise Cubic7AlVinError("every root must have one coordinate per diagonal-form entry")
    if not all(value.is_integral for value in diagonal):
        raise Cubic7AlVinError("the diagonal form is not integral over Z[theta]")
    if not all(entry.is_integral for root in root_vectors for entry in root):
        raise Cubic7AlVinError("a supplied root is not integral over Z[theta]")
    if len(set(root_vectors)) != len(root_vectors):
        raise Cubic7AlVinError("the supplied roots contain an exact duplicate")

    lattice_gram = _diagonal_matrix(diagonal)
    lattice_inertias = tuple(
        cubic7_inertia(lattice_gram, embedding=embedding) for embedding in range(3)
    )
    expected_lattice = (
        Inertia(dimension, 1, 0),
        Inertia(ambient_rank, 0, 0),
        Inertia(ambient_rank, 0, 0),
    )
    if lattice_inertias != expected_lattice:
        raise Cubic7AlVinError(
            "the diagonal form is not admissible: "
            f"inertias {lattice_inertias}, expected {expected_lattice}"
        )
    if _rank(root_vectors, ambient_rank) != ambient_rank:
        raise Cubic7AlVinError("the wall roots do not span the ambient Lorentz space")

    gram = _wall_gram(diagonal, root_vectors)
    norms = tuple(gram[index][index] for index in range(len(root_vectors)))
    if any(norm.sign(0) <= 0 for norm in norms):
        raise Cubic7AlVinError("every wall normal must be spacelike at the hyperbolic embedding")
    if any(
        gram[i][j].sign(0) > 0
        for i in range(len(root_vectors))
        for j in range(i)
    ):
        raise Cubic7AlVinError("the supplied root orientations are not pairwise acute")

    cartan = tuple(
        tuple(2 * gram[i][j] / norms[i] for j in range(len(root_vectors)))
        for i in range(len(root_vectors))
    )
    if any(not entry.is_integral for row in cartan for entry in row):
        raise Cubic7AlVinError("a root-to-root Cartan coefficient is not integral")
    reflection_matrices = tuple(
        tuple(
            tuple(
                Cubic7(int(row == column))
                - 2
                * root[row]
                * diagonal[column]
                * root[column]
                / norms[index]
                for column in range(ambient_rank)
            )
            for row in range(ambient_rank)
        )
        for index, root in enumerate(root_vectors)
    )
    if any(
        not entry.is_integral
        for reflection in reflection_matrices
        for row in reflection
        for entry in row
    ):
        raise Cubic7AlVinError(
            "a reflection fails to preserve the standard O_K lattice in the supplied basis"
        )

    pair_products = tuple(
        tuple(cartan[i][j] * cartan[j][i] for j in range(len(root_vectors)))
        for i in range(len(root_vectors))
    )
    _verify_diagram_labels(diagram, pair_products)

    wall_inertias = tuple(
        cubic7_inertia(gram, embedding=embedding) for embedding in range(3)
    )
    expected_wall = (
        Inertia(dimension, 1, len(root_vectors) - ambient_rank),
        Inertia(ambient_rank, 0, len(root_vectors) - ambient_rank),
        Inertia(ambient_rank, 0, len(root_vectors) - ambient_rank),
    )
    if wall_inertias != expected_wall:
        raise Cubic7AlVinError(
            f"wall Gram signatures are {wall_inertias}, expected {expected_wall}"
        )

    cyclic_generators = _fundamental_cycle_products(cartan)
    if any(not value.is_integral for value in cyclic_generators):
        raise Cubic7AlVinError("a Vinberg cyclic product generator is not integral")
    ground_field = "Q" if all(value.is_rational for value in cyclic_generators) else "Q(theta)"
    # If the cyclic ground field is the cubic K, the two checked conjugate
    # wall Grams are precisely Vinberg's remaining positivity conditions.
    # If it is Q there are no nonidentity real embeddings to check.

    rays = _enumerate_extreme_rays(diagonal, root_vectors, gram)
    if not rays:
        raise Cubic7AlVinError("the supplied half-space cone has no exact extreme rays")
    if any(not ray.link_verified for ray in rays):
        bad = [ray.active_facets for ray in rays if not ray.link_verified]
        raise Cubic7AlVinError(f"a finite/cusp vertex link failed Coxeter verification: {bad}")
    if any(
        _pair(diagonal, left.vector, right.vector).sign(0) > 0
        for index, left in enumerate(rays)
        for right in rays[index + 1 :]
        if left.kind != "hyperideal" and right.kind != "hyperideal"
    ):
        raise Cubic7AlVinError("nonspacelike chamber rays lie on inconsistent time-cone sheets")
    functionals = tuple(
        tuple(diagonal[index] * root[index] for index in range(ambient_rank))
        for root in root_vectors
    )
    irredundant = all(
        any(
            sum(
                (functionals[facet][index] * ray.vector[index] for index in range(ambient_rank)),
                start=CUBIC7_ZERO,
            ).sign(0)
            < 0
            for ray in rays
        )
        for facet in range(len(root_vectors))
    )
    if not irredundant:
        raise Cubic7AlVinError("at least one supplied wall is redundant in the half-space cone")

    has_hyperideal = any(ray.kind == "hyperideal" for ray in rays)
    has_cusp = any(ray.kind == "cusp" for ray in rays)
    finite_volume = not has_hyperideal
    compact = finite_volume and not has_cusp
    geometry_status = (
        "compact"
        if compact
        else "finite_volume_noncompact"
        if finite_volume
        else "infinite_volume_hyperideal"
    )
    canonical, fingerprint = canonical_pair_product_gram(pair_products, labelg=labelg)

    declared: bool | None = None
    if transcript_text is not None:
        if "The algorithm did not terminate" in transcript_text:
            declared = False
        elif "Algorithm ended" in transcript_text:
            declared = True
    return Cubic7AlVinVerification(
        dimension=dimension,
        diagonal_form=diagonal,
        roots=root_vectors,
        lattice_inertias=lattice_inertias,  # type: ignore[arg-type]
        wall_gram=gram,
        wall_inertias=wall_inertias,  # type: ignore[arg-type]
        cartan_matrix=cartan,
        pair_product_gram=pair_products,
        cyclic_field_generators=cyclic_generators,
        ground_field=ground_field,
        extreme_rays=rays,
        finite_volume=finite_volume,
        compact=compact,
        geometry_status=geometry_status,
        canonical_pair_product_gram=canonical,
        canonical_fingerprint=fingerprint,
        alvin_declared_termination=declared,
    )


def verify_cubic7_alvin_files(
    coxiter_path: str | Path,
    roots_path: str | Path,
    diagonal_form: str,
    *,
    labelg: str | Path | None = None,
) -> Cubic7AlVinVerification:
    transcript = Path(roots_path).read_text(encoding="utf-8")
    return verify_cubic7_alvin_candidate(
        parse_coxiter(coxiter_path),
        parse_diagonal_form(diagonal_form),
        parse_alvin_roots_text(transcript),
        transcript_text=transcript,
        labelg=labelg,
    )


__all__ = [
    "CoxIterDiagram",
    "Cubic7AlVinError",
    "Cubic7AlVinVerification",
    "ExtremeRay",
    "alvin_lambda_integer",
    "canonical_pair_product_gram",
    "parse_alvin_element",
    "parse_alvin_roots",
    "parse_alvin_roots_text",
    "parse_coxiter",
    "parse_coxiter_text",
    "parse_diagonal_form",
    "verify_cubic7_alvin_candidate",
    "verify_cubic7_alvin_files",
]
