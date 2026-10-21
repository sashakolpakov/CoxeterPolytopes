"""Independent exact verification of a hyperbolic Coxeter chamber.

This module consumes only serialized certificate data.  It does not call a
root search or trust a supplied vertex list: all extreme rays of the
polyhedral chamber cone are recomputed from the wall inequalities.

Wall completeness for a Vinberg search is deliberately outside this module;
that requires a separate, method-specific transcript verifier.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Any, Mapping, Sequence

from coxeter_polytopes.exact.coxeter import classify_root_pair
from coxeter_polytopes.exact.diagram import induced_subgram, is_parabolic, parabolic_rank
from coxeter_polytopes.exact.golden import Golden, ZERO
from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia


class CertificateError(ValueError):
    """Raised when an exact chamber certificate fails verification."""


@dataclass(frozen=True, slots=True)
class ChamberVerification:
    dimension: int
    facet_count: int
    finite_vertices: tuple[tuple[int, ...], ...]
    cusps: tuple[tuple[int, ...], ...]
    hyperideal_rays: tuple[tuple[int, ...], ...]
    opposite_sheet_rays: tuple[tuple[int, ...], ...]


@dataclass(frozen=True, slots=True)
class HyperidealRayVerification:
    """An exact extreme-ray obstruction to compactness or finite volume."""

    active_facets: tuple[int, ...]
    ray: tuple[Golden, ...]
    norm: Golden
    wall_pairings: tuple[Golden, ...]


def _require_exact_keys(value: Mapping[str, Any], expected: set[str], location: str) -> None:
    actual = set(value)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        details = []
        if missing:
            details.append(f"missing {missing}")
        if extra:
            details.append(f"unexpected {extra}")
        raise CertificateError(f"{location} has invalid properties: {', '.join(details)}")


def _verify_certificate_shape(certificate: Mapping[str, Any]) -> None:
    """Enforce the closed version-0.1.0 certificate schema semantically."""

    _require_exact_keys(
        certificate,
        {
            "schema_version",
            "certificate_kind",
            "field",
            "dimension",
            "lattice_gram",
            "simple_roots",
            "wall_gram",
            "coxeter_relations",
            "completeness",
            "geometry",
        },
        "certificate",
    )
    if certificate["schema_version"] != "0.1.0":
        raise CertificateError("schema_version must be 0.1.0")

    completeness = certificate["completeness"]
    if not isinstance(completeness, Mapping):
        raise CertificateError("completeness must be an object")
    _require_exact_keys(completeness, {"method", "evidence"}, "completeness")
    if not isinstance(completeness["method"], str) or not completeness["method"]:
        raise CertificateError("completeness.method must be a nonempty string")
    if not isinstance(completeness["evidence"], Mapping):
        raise CertificateError("completeness.evidence must be an object")

    geometry = certificate["geometry"]
    if not isinstance(geometry, Mapping):
        raise CertificateError("geometry must be an object")
    _require_exact_keys(
        geometry,
        {"nonempty", "interior_point", "vertices", "cusps"},
        "geometry",
    )
    if geometry["nonempty"] is not True:
        raise CertificateError("geometry.nonempty must be true")


def _golden(value: Any, location: str) -> Golden:
    if not isinstance(value, Mapping):
        raise CertificateError(f"{location} must be a serialized golden element")
    try:
        return Golden.from_json(value)
    except (TypeError, ValueError) as exc:
        raise CertificateError(f"invalid golden element at {location}: {exc}") from exc


def _vector(value: Any, location: str) -> list[Golden]:
    if not isinstance(value, list):
        raise CertificateError(f"{location} must be an array")
    return [_golden(entry, f"{location}[{i}]") for i, entry in enumerate(value)]


def _matrix(value: Any, location: str) -> list[list[Golden]]:
    if not isinstance(value, list):
        raise CertificateError(f"{location} must be an array")
    rows = [_vector(row, f"{location}[{i}]") for i, row in enumerate(value)]
    if rows and any(len(row) != len(rows[0]) for row in rows):
        raise CertificateError(f"{location} has rows of unequal length")
    return rows


def _require_square_symmetric(matrix: Sequence[Sequence[Golden]], name: str) -> None:
    if any(len(row) != len(matrix) for row in matrix):
        raise CertificateError(f"{name} must be square")
    if any(matrix[i][j] != matrix[j][i] for i in range(len(matrix)) for j in range(i)):
        raise CertificateError(f"{name} must be symmetric")


def _mat_vec(matrix: Sequence[Sequence[Golden]], vector: Sequence[Golden]) -> list[Golden]:
    return [sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def _dot(left: Sequence[Golden], right: Sequence[Golden]) -> Golden:
    return sum((left[i] * right[i] for i in range(len(left))), start=ZERO)


def _pair(
    left: Sequence[Golden],
    gram: Sequence[Sequence[Golden]],
    right: Sequence[Golden],
) -> Golden:
    return _dot(left, _mat_vec(gram, right))


def _root_gram(
    roots: Sequence[Sequence[Golden]], gram: Sequence[Sequence[Golden]]
) -> list[list[Golden]]:
    return [[_pair(left, gram, right) for right in roots] for left in roots]


def _cofactor_kernel(rows: Sequence[Sequence[Golden]]) -> list[Golden]:
    """Kernel vector for a full-rank d by (d+1) matrix."""

    dimension = len(rows)
    if any(len(row) != dimension + 1 for row in rows):
        raise ValueError("cofactor kernel requires a d by (d+1) matrix")
    return [
        ((-1) ** column)
        * determinant([list(row[:column]) + list(row[column + 1 :]) for row in rows])
        for column in range(dimension + 1)
    ]


def _full_rank_kernel(rows: Sequence[Sequence[Golden]]) -> list[Golden] | None:
    """Kernel vector for a full-rank ``d`` by ``d+1`` matrix.

    A single exact row reduction is substantially cheaper than separately
    evaluating all ``d+1`` maximal minors.  ``None`` records rank below ``d``.
    Projective ray enumeration does not depend on the kernel vector's scale.
    """

    dimension = len(rows)
    if any(len(row) != dimension + 1 for row in rows):
        raise ValueError("kernel reduction requires a d by (d+1) matrix")
    work = [list(row) for row in rows]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(dimension + 1):
        pivot = next(
            (row for row in range(pivot_row, dimension) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(dimension):
            if row == pivot_row or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                work[row][j] - multiplier * work[pivot_row][j]
                for j in range(dimension + 1)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == dimension:
            break
    if pivot_row != dimension:
        return None
    free_column = next(
        column for column in range(dimension + 1) if column not in pivot_columns
    )
    kernel = [ZERO for _ in range(dimension + 1)]
    kernel[free_column] = Golden(1)
    for row, column in enumerate(pivot_columns):
        kernel[column] = -work[row][free_column]
    return kernel


def _row_rank(rows: Sequence[Sequence[Golden]], column_count: int) -> int:
    """Exact row rank, specialized to the small wall-functional matrices."""

    work = [list(row) for row in rows]
    if any(len(row) != column_count for row in work):
        raise ValueError("row-rank matrix has the wrong column count")
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        for row in range(pivot_row + 1, len(work)):
            if not work[row][column]:
                continue
            multiplier = work[row][column] / pivot_value
            work[row] = [
                work[row][j] - multiplier * work[pivot_row][j]
                for j in range(column_count)
            ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def _canonical_ray(vector: Sequence[Golden]) -> tuple[Golden, ...]:
    pivot = next((entry for entry in vector if entry), None)
    if pivot is None:
        raise ValueError("zero vector does not define a projective ray")
    return tuple(entry / pivot for entry in vector)


def _declared_faces(value: Any, name: str, facet_count: int) -> tuple[tuple[int, ...], ...]:
    if not isinstance(value, list):
        raise CertificateError(f"geometry.{name} must be an array")
    faces = []
    for i, face in enumerate(value):
        if not isinstance(face, list) or not all(isinstance(index, int) for index in face):
            raise CertificateError(f"geometry.{name}[{i}] must be an integer array")
        if len(face) != len(set(face)):
            raise CertificateError(f"geometry.{name}[{i}] repeats a facet")
        if any(index < 0 or index >= facet_count for index in face):
            raise CertificateError(f"geometry.{name}[{i}] has an out-of-range facet")
        faces.append(tuple(sorted(face)))
    if len(faces) != len(set(faces)):
        raise CertificateError(f"geometry.{name} contains duplicate faces")
    return tuple(sorted(faces))


def _expected_relation(relation) -> int | str:
    if relation.kind == "finite":
        assert relation.label is not None
        return relation.label
    if relation.kind == "parallel":
        return "parallel"
    if relation.kind == "ultraparallel":
        return "ultraparallel"
    raise CertificateError(
        f"finite wall angle has unrecognized exact Coxeter invariant q={relation.q}"
    )


def _verify_relations(
    declared: Any, wall_gram: Sequence[Sequence[Golden]]
) -> None:
    facet_count = len(wall_gram)
    if not isinstance(declared, list) or len(declared) != facet_count:
        raise CertificateError("coxeter_relations must be a square facet-count matrix")
    if any(not isinstance(row, list) or len(row) != facet_count for row in declared):
        raise CertificateError("coxeter_relations must be a square facet-count matrix")
    for i in range(facet_count):
        if declared[i][i] != 1:
            raise CertificateError("coxeter_relations must have diagonal 1")
        for j in range(i):
            if declared[i][j] != declared[j][i]:
                raise CertificateError("coxeter_relations must be symmetric")
            try:
                relation = classify_root_pair(
                    wall_gram[i][i], wall_gram[i][j], wall_gram[j][j]
                )
            except ValueError as exc:
                raise CertificateError(f"invalid oriented wall pair ({i},{j}): {exc}") from exc
            expected = _expected_relation(relation)
            if declared[i][j] != expected:
                raise CertificateError(
                    f"wrong Coxeter relation at ({i},{j}): "
                    f"declared {declared[i][j]!r}, exact value is {expected!r}"
                )


def _verify_integral_roots(
    roots: Sequence[Sequence[Golden]],
    gram: Sequence[Sequence[Golden]],
    wall_gram: Sequence[Sequence[Golden]],
) -> None:
    for index, root in enumerate(roots):
        if not all(entry.is_integral for entry in root):
            raise CertificateError(f"simple_roots[{index}] is not in Z[phi]^n")
        norm = wall_gram[index][index]
        if norm.sign(0) <= 0:
            raise CertificateError(f"simple_roots[{index}] is not spacelike")
        coefficients = [2 * entry / norm for entry in _mat_vec(gram, root)]
        if not all(entry.is_integral for entry in coefficients):
            raise CertificateError(f"simple_roots[{index}] does not define an integral reflection")


def _enumerate_extreme_rays(
    roots: Sequence[Sequence[Golden]],
    gram: Sequence[Sequence[Golden]],
    interior: Sequence[Golden],
    dimension: int,
) -> tuple[
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
]:
    facet_count = len(roots)
    functionals = [_mat_vec(gram, root) for root in roots]
    rays: dict[tuple[Golden, ...], tuple[list[Golden], tuple[int, ...]]] = {}

    for subset in combinations(range(facet_count), dimension):
        active_rows = [functionals[index] for index in subset]
        ray = _full_rank_kernel(active_rows)
        if ray is None:
            continue
        values = [_dot(ray, functional) for functional in functionals]
        if all(value.sign(0) <= 0 for value in values):
            oriented = ray
            oriented_values = values
        elif all(value.sign(0) >= 0 for value in values):
            oriented = [-entry for entry in ray]
            oriented_values = [-entry for entry in values]
        else:
            continue
        active = tuple(i for i, value in enumerate(oriented_values) if not value)
        rays[_canonical_ray(oriented)] = (oriented, active)

    if not rays:
        raise CertificateError("the chamber cone has no extreme rays")

    finite: set[tuple[int, ...]] = set()
    cusps: set[tuple[int, ...]] = set()
    hyperideal: set[tuple[int, ...]] = set()
    opposite: set[tuple[int, ...]] = set()
    for ray, active in rays.values():
        norm_sign = _pair(ray, gram, ray).sign(0)
        sheet_sign = _pair(ray, gram, interior).sign(0)
        if norm_sign < 0:
            if sheet_sign < 0:
                finite.add(active)
            else:
                opposite.add(active)
        elif norm_sign == 0:
            if sheet_sign < 0:
                cusps.add(active)
            else:
                opposite.add(active)
        else:
            hyperideal.add(active)
    return tuple(sorted(finite)), tuple(sorted(cusps)), tuple(sorted(hyperideal)), tuple(sorted(opposite))


def _enumerate_extreme_rays_incremental(
    roots: Sequence[Sequence[Golden]],
    gram: Sequence[Sequence[Golden]],
    interior: Sequence[Golden],
    dimension: int,
) -> tuple[
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
]:
    """Enumerate rays by an exact incremental double-description update."""

    ambient_rank = dimension + 1
    facet_count = len(roots)
    functionals = [_mat_vec(gram, root) for root in roots]
    seed: tuple[int, ...] | None = None
    for candidate in combinations(range(facet_count), ambient_rank):
        if _row_rank([functionals[index] for index in candidate], ambient_rank) == ambient_rank:
            seed = candidate
            break
    if seed is None:
        raise CertificateError("wall functionals do not span the ambient dual space")

    seed_set = set(seed)
    ordered_facets = (*seed, *(i for i in range(facet_count) if i not in seed_set))
    processed = list(seed)
    rays: dict[tuple[Golden, ...], tuple[list[Golden], frozenset[int]]] = {}
    for omitted in seed:
        active = tuple(index for index in seed if index != omitted)
        ray = _full_rank_kernel([functionals[index] for index in active])
        if ray is None:
            raise CertificateError("independent seed walls have a dependent proper subset")
        values = [_dot(ray, functionals[index]) for index in seed]
        if all(value.sign(0) <= 0 for value in values):
            oriented = ray
        elif all(value.sign(0) >= 0 for value in values):
            oriented = [-entry for entry in ray]
        else:
            raise CertificateError("seed wall cone has inconsistent orientation")
        key = _canonical_ray(oriented)
        actual = frozenset(
            index for index in seed if not _dot(oriented, functionals[index])
        )
        rays[key] = (oriented, actual)

    rank_cache: dict[frozenset[int], int] = {}
    for new_facet in ordered_facets[ambient_rank:]:
        functional = functionals[new_facet]
        negative: list[tuple[list[Golden], frozenset[int], Golden]] = []
        zero: list[tuple[list[Golden], frozenset[int], Golden]] = []
        positive: list[tuple[list[Golden], frozenset[int], Golden]] = []
        for ray, active in rays.values():
            value = _dot(ray, functional)
            target = negative if value.sign(0) < 0 else positive if value.sign(0) > 0 else zero
            target.append((ray, active, value))

        updated: dict[tuple[Golden, ...], tuple[list[Golden], frozenset[int]]] = {}
        for ray, active, _ in (*negative, *zero):
            new_active = active | ({new_facet} if _dot(ray, functional) == 0 else set())
            updated[_canonical_ray(ray)] = (ray, frozenset(new_active))

        for negative_ray, negative_active, negative_value in negative:
            for positive_ray, positive_active, positive_value in positive:
                common = negative_active & positive_active
                common_rank = rank_cache.get(common)
                if common_rank is None:
                    common_rank = _row_rank(
                        [functionals[index] for index in common], ambient_rank
                    )
                    rank_cache[common] = common_rank
                if common_rank != ambient_rank - 2:
                    continue
                ray = [
                    positive_value * negative_ray[i]
                    - negative_value * positive_ray[i]
                    for i in range(ambient_rank)
                ]
                key = _canonical_ray(ray)
                actual = frozenset(
                    index
                    for index in (*processed, new_facet)
                    if not _dot(ray, functionals[index])
                )
                updated[key] = (ray, actual)
        rays = updated
        processed.append(new_facet)

    finite: set[tuple[int, ...]] = set()
    cusps: set[tuple[int, ...]] = set()
    hyperideal: set[tuple[int, ...]] = set()
    opposite: set[tuple[int, ...]] = set()
    for ray, active_set in rays.values():
        active = tuple(sorted(active_set))
        norm_sign = _pair(ray, gram, ray).sign(0)
        sheet_sign = _pair(ray, gram, interior).sign(0)
        if norm_sign < 0:
            if sheet_sign < 0:
                finite.add(active)
            else:
                opposite.add(active)
        elif norm_sign == 0:
            if sheet_sign < 0:
                cusps.add(active)
            else:
                opposite.add(active)
        else:
            hyperideal.add(active)
    return tuple(sorted(finite)), tuple(sorted(cusps)), tuple(sorted(hyperideal)), tuple(sorted(opposite))


def diagnose_chamber_cone(
    roots: Sequence[Sequence[Golden]],
    gram: Sequence[Sequence[Golden]],
    interior: Sequence[Golden],
    dimension: int,
) -> ChamberVerification:
    """Classify every extreme ray without requiring a finite-volume claim.

    This is the geometry engine used by :func:`verify_chamber_geometry`,
    exposed separately for bounded searches whose current wall cone may still
    have hyperideal directions.
    """

    ambient_rank = dimension + 1
    if len(gram) != ambient_rank or any(len(row) != ambient_rank for row in gram):
        raise CertificateError("diagnostic lattice Gram has the wrong rank")
    if len(interior) != ambient_rank or any(len(root) != ambient_rank for root in roots):
        raise CertificateError("diagnostic vectors have the wrong rank")
    _require_square_symmetric(gram, "diagnostic lattice Gram")
    if _pair(interior, gram, interior).sign(0) >= 0:
        raise CertificateError("diagnostic interior point must be timelike")
    if any(_pair(interior, gram, root).sign(0) >= 0 for root in roots):
        raise CertificateError("diagnostic interior point must satisfy every wall strictly")
    finite, cusps, hyperideal, opposite = _enumerate_extreme_rays_incremental(
        roots, gram, interior, dimension
    )
    return ChamberVerification(
        dimension=dimension,
        facet_count=len(roots),
        finite_vertices=finite,
        cusps=cusps,
        hyperideal_rays=hyperideal,
        opposite_sheet_rays=opposite,
    )


def verify_hyperideal_extreme_ray(
    roots: Sequence[Sequence[Golden]],
    gram: Sequence[Sequence[Golden]],
    active_facets: Sequence[int],
    dimension: int,
) -> HyperidealRayVerification:
    """Recompute one spacelike extreme ray of a wall cone exactly.

    A full-rank set of ``dimension`` active walls has a one-dimensional
    kernel.  If its orientation satisfies every wall inequality and its norm
    is positive, it is a concise certificate that the current chamber cone is
    not yet compact or finite-volume.
    """

    ambient_rank = dimension + 1
    active = tuple(active_facets)
    if len(gram) != ambient_rank or any(len(row) != ambient_rank for row in gram):
        raise CertificateError("hyperideal witness Gram has the wrong rank")
    if any(len(root) != ambient_rank for root in roots):
        raise CertificateError("hyperideal witness roots have the wrong rank")
    if len(active) != dimension or len(set(active)) != dimension:
        raise CertificateError("a hyperideal witness needs distinct dimension-many facets")
    if any(index < 0 or index >= len(roots) for index in active):
        raise CertificateError("hyperideal witness has an out-of-range facet")
    functionals = [_mat_vec(gram, root) for root in roots]
    ray = _full_rank_kernel([functionals[index] for index in active])
    if ray is None:
        raise CertificateError("hyperideal witness facets are not independent")
    pairings = [_dot(ray, functional) for functional in functionals]
    if all(value.sign(0) <= 0 for value in pairings):
        oriented = ray
        oriented_pairings = pairings
    elif all(value.sign(0) >= 0 for value in pairings):
        oriented = [-entry for entry in ray]
        oriented_pairings = [-value for value in pairings]
    else:
        raise CertificateError("hyperideal witness ray violates a wall inequality")
    actual_active = tuple(i for i, value in enumerate(oriented_pairings) if not value)
    if actual_active != tuple(sorted(active)):
        raise CertificateError(
            f"hyperideal witness has active facets {actual_active}, not {tuple(sorted(active))}"
        )
    canonical = _canonical_ray(oriented)
    norm = _pair(canonical, gram, canonical)
    if norm.sign(0) <= 0:
        raise CertificateError("claimed hyperideal extreme ray is not spacelike")
    return HyperidealRayVerification(
        active_facets=actual_active,
        ray=canonical,
        norm=norm,
        wall_pairings=tuple(_pair(canonical, gram, root) for root in roots),
    )


def verify_chamber_geometry(certificate: Mapping[str, Any]) -> ChamberVerification:
    """Verify all algebraic, Coxeter, and polyhedral geometry in a certificate.

    This proves that the listed walls bound the declared compact or cofinite
    Coxeter polyhedron.  It does not prove that a Vinberg search found every
    wall of a larger ambient reflection group.
    """

    _verify_certificate_shape(certificate)

    field = certificate.get("field")
    if field not in {"Q", "Q(sqrt(5))"}:
        raise CertificateError("field must be Q or Q(sqrt(5))")
    dimension = certificate.get("dimension")
    if not isinstance(dimension, int) or dimension < 2:
        raise CertificateError("dimension must be an integer at least 2")
    kind = certificate.get("certificate_kind")
    if kind not in {"compact_reflective", "cofinite_noncompact"}:
        raise CertificateError("unsupported certificate_kind")

    gram = _matrix(certificate.get("lattice_gram"), "lattice_gram")
    _require_square_symmetric(gram, "lattice_gram")
    ambient_rank = dimension + 1
    if len(gram) != ambient_rank:
        raise CertificateError("lattice_gram rank must equal dimension+1")
    expected_identity = Inertia(dimension, 1, 0)
    if inertia(gram, embedding=0) != expected_identity:
        raise CertificateError(f"lattice_gram identity inertia must be {expected_identity}")
    if field == "Q":
        if any(not entry.is_rational for row in gram for entry in row):
            raise CertificateError("a Q certificate has a nonrational lattice entry")
    else:
        expected_conjugate = Inertia(ambient_rank, 0, 0)
        if inertia(gram, embedding=1) != expected_conjugate:
            raise CertificateError(
                "a Q(sqrt(5)) lattice must be positive definite at the conjugate embedding"
            )

    roots_value = certificate.get("simple_roots")
    if not isinstance(roots_value, list) or not roots_value:
        raise CertificateError("simple_roots must be a nonempty array")
    roots = [_vector(root, f"simple_roots[{i}]") for i, root in enumerate(roots_value)]
    if any(len(root) != ambient_rank for root in roots):
        raise CertificateError("every simple root must have dimension+1 coordinates")
    if field == "Q" and any(not entry.is_rational for root in roots for entry in root):
        raise CertificateError("a Q certificate has a nonrational root coordinate")
    facet_count = len(roots)

    computed_wall_gram = _root_gram(roots, gram)
    declared_wall_gram = _matrix(certificate.get("wall_gram"), "wall_gram")
    _require_square_symmetric(declared_wall_gram, "wall_gram")
    if declared_wall_gram != computed_wall_gram:
        raise CertificateError("wall_gram does not equal R^T*lattice_gram*R")
    wall_inertia = inertia(computed_wall_gram)
    expected_wall_inertia = Inertia(dimension, 1, facet_count - ambient_rank)
    if wall_inertia != expected_wall_inertia:
        raise CertificateError(f"wall_gram inertia must be {expected_wall_inertia}")

    _verify_integral_roots(roots, gram, computed_wall_gram)
    _verify_relations(certificate.get("coxeter_relations"), computed_wall_gram)

    geometry = certificate.get("geometry")
    if not isinstance(geometry, Mapping):
        raise CertificateError("geometry must be an object")
    interior = _vector(geometry.get("interior_point"), "geometry.interior_point")
    if len(interior) != ambient_rank:
        raise CertificateError("geometry.interior_point must have dimension+1 coordinates")
    if field == "Q" and any(not entry.is_rational for entry in interior):
        raise CertificateError("a Q certificate has a nonrational interior point")
    if _pair(interior, gram, interior).sign(0) >= 0:
        raise CertificateError("geometry.interior_point must be timelike")
    if any(_pair(interior, gram, root).sign(0) >= 0 for root in roots):
        raise CertificateError("geometry.interior_point must satisfy every wall strictly")

    finite, cusps, hyperideal, opposite = _enumerate_extreme_rays_incremental(
        roots, gram, interior, dimension
    )
    for cusp in cusps:
        cusp_gram = induced_subgram(computed_wall_gram, cusp)
        if not is_parabolic(cusp_gram):
            raise CertificateError(
                f"cusp {cusp} is not a componentwise affine parabolic diagram"
            )
        if parabolic_rank(cusp_gram) != dimension - 1:
            raise CertificateError(
                f"cusp {cusp} does not have parabolic rank {dimension - 1}"
            )
    declared_vertices = _declared_faces(geometry.get("vertices"), "vertices", facet_count)
    declared_cusps = _declared_faces(geometry.get("cusps"), "cusps", facet_count)
    if declared_vertices != finite:
        raise CertificateError(
            f"declared vertices do not match exact extreme rays: {declared_vertices} != {finite}"
        )
    if declared_cusps != cusps:
        raise CertificateError(
            f"declared cusps do not match exact extreme rays: {declared_cusps} != {cusps}"
        )
    if hyperideal:
        raise CertificateError(f"chamber has hyperideal extreme rays: {hyperideal}")
    if opposite:
        raise CertificateError(f"chamber cone reaches the opposite time sheet: {opposite}")
    if kind == "compact_reflective" and cusps:
        raise CertificateError("a compact certificate cannot contain cusps")
    if kind == "cofinite_noncompact" and not cusps:
        raise CertificateError("a cofinite_noncompact certificate must contain a cusp")

    return ChamberVerification(
        dimension=dimension,
        facet_count=facet_count,
        finite_vertices=finite,
        cusps=cusps,
        hyperideal_rays=hyperideal,
        opposite_sheet_rays=opposite,
    )
