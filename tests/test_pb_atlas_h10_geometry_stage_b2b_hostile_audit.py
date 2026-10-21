from __future__ import annotations

import itertools
import json
from pathlib import Path
import shutil
import subprocess
from typing import Any, Mapping, Sequence

import pytest
import sympy as sp

from scripts.pb_atlas_cluster_stage_b1 import compile_endpoint_program


ROOT = Path(__file__).resolve().parents[1]
STAGE_A = ROOT / "data/enumerations/pb_atlas_numerical_contraction_stage_a.json"
ARCHIVE = ROOT / "data/enumerations/pb_atlas_h10_geometry_stage_b2b.json"

X, Y, Z = sp.symbols("X Y Z")


def _forest_component_sizes(packet: Mapping[str, Any]) -> list[int]:
    centres = packet["reverse_blowup_program"]["center_records"]
    adjacency = [set() for _ in centres]
    for child, centre in enumerate(centres):
        for parent in centre["parent_center_positions_zero_based"]:
            adjacency[child].add(parent)
            adjacency[parent].add(child)
    unseen = set(range(len(centres)))
    sizes = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        stack = [seed]
        size = 0
        while stack:
            vertex = stack.pop()
            size += 1
            for neighbour in adjacency[vertex]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    stack.append(neighbour)
        sizes.append(size)
    return sorted(sizes, reverse=True)


@pytest.fixture(scope="module")
def archive() -> dict[str, Any]:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def packet() -> dict[str, Any]:
    stage_a = json.loads(STAGE_A.read_text(encoding="utf-8"))
    matches = []
    for case in stage_a["cases"]:
        for orbit in case["certificate"]["P2_endpoint_orbits"]:
            candidate = compile_endpoint_program(orbit["program"])
            if _forest_component_sizes(candidate) == [2, 2, 2, 1, 1, 1, 1]:
                matches.append(candidate)
    assert len(matches) == 1
    return matches[0]


def _chart_model(
    packet: Mapping[str, Any],
    values: Mapping[sp.Symbol, sp.Rational] | None = None,
) -> dict[str, Any]:
    centres = packet["reverse_blowup_program"]["center_records"]
    proper = [
        index
        for index, centre in enumerate(centres)
        if centre["center_kind"] == "proper"
    ]
    gauge = proper[:4]
    targets = (
        (sp.Integer(1), sp.Integer(0), sp.Integer(0)),
        (sp.Integer(0), sp.Integer(1), sp.Integer(0)),
        (sp.Integer(0), sp.Integer(0), sp.Integer(1)),
        (sp.Integer(1), sp.Integer(1), sp.Integer(1)),
    )
    points: dict[int, sp.Matrix] = {}
    tangent_frames: dict[int, tuple[sp.Matrix, sp.Matrix]] = {}
    patches: dict[int, int] = {}
    for position in proper:
        if position in gauge:
            point = sp.Matrix(targets[gauge.index(position)])
            patch = next(index for index, entry in enumerate(point) if entry)
        else:
            point = sp.Matrix(
                [sp.Symbol(f"q{position}_0"), sp.Symbol(f"q{position}_1"), 1]
            )
            patch = 2
        if values is not None:
            point = point.subs(values, simultaneous=True)
        free_coordinates = [index for index in range(3) if index != patch]
        frame = []
        for coordinate in free_coordinates:
            vector = sp.zeros(3, 1)
            vector[coordinate] = 1
            frame.append(vector)
        points[position] = point
        tangent_frames[position] = (frame[0], frame[1])
        patches[position] = patch
    return {
        "centres": centres,
        "proper": proper,
        "gauge": gauge,
        "points": points,
        "tangent_frames": tangent_frames,
        "patches": patches,
    }


def _line_root_determinant(
    packet: Mapping[str, Any],
    row: Mapping[str, Any],
    modes: Mapping[int, str],
) -> sp.Expr:
    model = _chart_model(packet)
    constraints = []
    multiplicities = row["multiplicities_chronological"]
    for position, centre in enumerate(model["centres"]):
        if multiplicities[position] == 0:
            continue
        assert multiplicities[position] == 1
        if centre["center_kind"] == "proper":
            constraints.append(list(model["points"][position]))
            continue
        assert centre["center_kind"] == "free"
        parent = centre["parent_center_positions_zero_based"][0]
        first, second = model["tangent_frames"][parent]
        slope = sp.Symbol(f"t{position}")
        direction = (
            first + slope * second
            if modes[position] == "A"
            else slope * first + second
        )
        constraints.append(list(direction))
    assert len(constraints) == 3
    return sp.factor(sp.Matrix(constraints).det())


def _root_rows(packet: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    rows = [
        row
        for row in packet["divisor_table"]
        if not row["is_contracted_wall"]
        and row["all_multiplicities_nonnegative"]
        and row["self_intersection"] == -2
        and row["plane_degree"] == 1
    ]
    assert len(rows) == 9
    return rows


def _normalized_projective(point: sp.Matrix) -> tuple[sp.Rational, ...]:
    entries = [sp.Rational(value) for value in point]
    pivot = next(value for value in entries if value)
    return tuple(sp.cancel(value / pivot) for value in entries)


def _split_at_order(
    expression: sp.Expr, u: sp.Symbol, v: sp.Symbol, multiplicity: int
) -> tuple[dict[tuple[int, int], sp.Expr], sp.Expr]:
    low: dict[tuple[int, int], sp.Expr] = {}
    high = sp.Integer(0)
    for exponents, coefficient in sp.Poly(sp.expand(expression), u, v).terms():
        if sum(exponents) < multiplicity:
            low[exponents] = coefficient
        else:
            high += coefficient * u ** exponents[0] * v ** exponents[1]
    return low, sp.expand(high)


def _closed_coordinates(archive: Mapping[str, Any]) -> dict[sp.Symbol, sp.Rational]:
    record = next(
        record
        for record in archive["generated_incidence"]["component_records"]
        if record["collision_free"]
    )
    return {
        sp.Symbol(name): sp.Rational(value)
        for name, value in record["coordinates"].items()
    }


def _surviving_modes(archive: Mapping[str, Any]) -> dict[int, str]:
    record = next(
        record
        for record in archive["generated_incidence"]["all_chart_records"]
        if record["incidence_possible_on_chart"]
    )
    return {int(position): mode for position, mode in record["free_chart_modes"].items()}


def _proper_substitution(
    model: Mapping[str, Any], position: int, u: sp.Symbol, v: sp.Symbol
) -> dict[sp.Symbol, sp.Expr]:
    point = model["points"][position]
    first, second = model["tangent_frames"][position]
    local = point + u * first + v * second
    return {X: local[0], Y: local[1], Z: local[2]}


def _independent_jet_matrix(
    packet: Mapping[str, Any],
    *,
    degree: int,
    multiplicities: Sequence[int],
    values: Mapping[sp.Symbol, sp.Rational],
    modes: Mapping[int, str],
) -> tuple[sp.Matrix, list[sp.Expr]]:
    model = _chart_model(packet, values)
    centres = model["centres"]
    monomials = [
        X**left * Y**middle * Z ** (degree - left - middle)
        for left in range(degree + 1)
        for middle in range(degree - left + 1)
    ]
    local_variables = {
        position: sp.symbols(f"hostile_u{position} hostile_v{position}")
        for position in range(len(centres))
    }
    current: dict[int, list[sp.Expr]] = {}
    for position in model["proper"]:
        u, v = local_variables[position]
        substitution = _proper_substitution(model, position, u, v)
        current[position] = [
            sp.expand(monomial.subs(substitution, simultaneous=True))
            for monomial in monomials
        ]

    children: dict[int, list[int]] = {}
    for position, centre in enumerate(centres):
        if centre["center_kind"] == "free":
            parent = int(centre["parent_center_positions_zero_based"][0])
            children.setdefault(parent, []).append(position)

    rows = []
    for position, centre in enumerate(centres):
        assert position in current
        u, v = local_variables[position]
        split = [
            _split_at_order(expression, u, v, multiplicities[position])
            for expression in current[position]
        ]
        for total in range(multiplicities[position]):
            for first in range(total + 1):
                multiindex = (first, total - first)
                rows.append([low.get(multiindex, 0) for low, _high in split])

        for child in children.get(position, []):
            child_u, child_v = local_variables[child]
            slope = values[sp.Symbol(f"t{child}")]
            if modes[child] == "A":
                substitution = {u: child_u, v: child_u * (slope + child_v)}
                divisor = child_u
            else:
                substitution = {u: child_v * (slope + child_u), v: child_v}
                divisor = child_v
            current[child] = [
                sp.expand(
                    sp.cancel(
                        high.subs(substitution, simultaneous=True)
                        / divisor ** multiplicities[position]
                    )
                )
                for _low, high in split
            ]
    return sp.Matrix(rows), monomials


def _polynomial(text: str) -> sp.Expr:
    return sp.expand(
        sp.sympify(text.replace("^", "**"), locals={"X": X, "Y": Y, "Z": Z})
    )


def _total_order(expression: sp.Expr, u: sp.Symbol, v: sp.Symbol) -> int:
    terms = [
        exponents
        for exponents, coefficient in sp.Poly(sp.expand(expression), u, v).terms()
        if coefficient != 0
    ]
    assert terms
    return min(sum(exponents) for exponents in terms)


def _independent_cluster_orders(
    packet: Mapping[str, Any],
    polynomial: sp.Expr,
    multiplicities: Sequence[int],
    values: Mapping[sp.Symbol, sp.Rational],
    modes: Mapping[int, str],
) -> list[int]:
    model = _chart_model(packet, values)
    centres = model["centres"]
    proper_local: dict[int, tuple[sp.Expr, sp.Symbol, sp.Symbol]] = {}
    orders: list[int | None] = [None] * len(centres)
    for position in model["proper"]:
        u, v = sp.symbols(f"order_u{position} order_v{position}")
        expression = sp.expand(
            polynomial.subs(
                _proper_substitution(model, position, u, v), simultaneous=True
            )
        )
        proper_local[position] = (expression, u, v)
        orders[position] = _total_order(expression, u, v)

    for child, centre in enumerate(centres):
        if centre["center_kind"] != "free":
            continue
        parent = int(centre["parent_center_positions_zero_based"][0])
        expression, u, v = proper_local[parent]
        low, _high = _split_at_order(expression, u, v, multiplicities[parent])
        assert not any(low.values())
        child_u, child_v = sp.symbols(f"order_u{child} order_v{child}")
        slope = values[sp.Symbol(f"t{child}")]
        if modes[child] == "A":
            substitution = {u: child_u, v: child_u * (slope + child_v)}
            divisor = child_u
        else:
            substitution = {u: child_v * (slope + child_u), v: child_v}
            divisor = child_v
        strict_transform = sp.expand(
            sp.cancel(
                expression.subs(substitution, simultaneous=True)
                / divisor ** multiplicities[parent]
            )
        )
        orders[child] = _total_order(strict_transform, child_u, child_v)
    assert all(value is not None for value in orders)
    return [int(value) for value in orders]


def _bcg_absolute_factor_record(polynomial: sp.Expr) -> tuple[int, int, int]:
    executable = shutil.which("Singular")
    assert executable is not None
    source = str(polynomial).replace("**", "^").replace("X", "x").replace("Y", "y").replace("Z", "z")
    program = f'''
ring r=0,(x,y,z),dp;
poly f={source};
LIB "absfact.lib";
def A=absFactorizeBCG(f);
setring A;
print("HOSTILE_FACTOR_COUNT"); print(size(absolute_factors[1]));
print("HOSTILE_COMPONENT_COUNT"); print(absolute_factors[4]);
print("HOSTILE_FIRST_MULTIPLICITY"); print(absolute_factors[2][2]);
exit;
'''
    completed = subprocess.run(
        [executable, "-q", "--no-rc"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    labels = {
        "HOSTILE_FACTOR_COUNT": 0,
        "HOSTILE_COMPONENT_COUNT": 1,
        "HOSTILE_FIRST_MULTIPLICITY": 2,
    }
    parsed = [0, 0, 0]
    for label, position in labels.items():
        index = lines.index(label)
        parsed[position] = int(lines[index + 1])
    return tuple(parsed)


def test_independent_incidence_chart_component_collision_and_open(
    packet: dict[str, Any], archive: dict[str, Any]
) -> None:
    centres = packet["reverse_blowup_program"]["center_records"]
    free = [
        index
        for index, centre in enumerate(centres)
        if centre["center_kind"] == "free"
    ]
    assert free == [3, 6, 9]
    rows = _root_rows(packet)
    viable = []
    for mode_tuple in itertools.product(("A", "B"), repeat=len(free)):
        modes = dict(zip(free, mode_tuple))
        determinants = {
            int(row["wall_index_zero_based"]): _line_root_determinant(
                packet, row, modes
            )
            for row in rows
        }
        unit_walls = sorted(
            wall
            for wall, determinant in determinants.items()
            if not determinant.free_symbols and determinant != 0
        )
        frozen = next(
            record
            for record in archive["generated_incidence"]["all_chart_records"]
            if record["free_chart_modes"]
            == {str(position): mode for position, mode in modes.items()}
        )
        assert unit_walls == frozen["constant_nonzero_obstruction_walls"]
        for wall, determinant in determinants.items():
            frozen_determinant = sp.sympify(
                frozen["root_determinants"][str(wall)].replace("^", "**")
            )
            if determinant == 0:
                assert frozen_determinant == 0
            else:
                ratio = sp.cancel(frozen_determinant / determinant)
                assert not ratio.free_symbols and ratio != 0
        if not unit_walls:
            viable.append((modes, list(determinants.values())))
    assert len(viable) == 1
    assert viable[0][0] == {3: "A", 6: "B", 9: "B"}

    variables = [
        sp.Symbol(name)
        for name in archive["generated_incidence"]["coordinate_variables"]
    ]
    groebner = sp.groebner(
        [equation for equation in viable[0][1] if equation != 0],
        *variables,
        order="lex",
    )
    basis = [sp.factor(value.as_expr()) for value in groebner.polys]
    assert basis == [
        sp.Symbol("q5_0") - sp.Symbol("q8_1"),
        sp.Symbol("q5_1") - 1,
        sp.Symbol("q7_0") - 1,
        sp.Symbol("q7_1") - sp.Symbol("q8_1"),
        sp.Symbol("q8_0") - sp.Symbol("q8_1"),
        (sp.Symbol("q8_1") - 1) * (sp.Symbol("q8_1") + 1),
        sp.Symbol("t3"),
        sp.Symbol("t6"),
        sp.Symbol("t9"),
    ]

    component_by_sign = {
        sp.Rational(record["coordinates"]["q8_1"]): record
        for record in archive["generated_incidence"]["component_records"]
    }
    assert set(component_by_sign) == {-1, 1}
    for sign, record in component_by_sign.items():
        values = {
            sp.Symbol(name): sp.Rational(value)
            for name, value in record["coordinates"].items()
        }
        model = _chart_model(packet, values)
        points = {
            position: _normalized_projective(point)
            for position, point in model["points"].items()
        }
        collisions = sorted(
            [left, right]
            for left, right in itertools.combinations(model["proper"], 2)
            if points[left] == points[right]
        )
        assert collisions == record["proper_center_collision_pairs_zero_based"]
        for position in model["proper"]:
            assert model["points"][position][model["patches"][position]] != 0
        determinants = [
            sp.Matrix.hstack(
                *(model["points"][position] for position in triple)
            ).det()
            for triple in itertools.combinations(model["gauge"], 3)
        ]
        assert all(value != 0 for value in determinants)
        if sign == -1:
            assert collisions == []
        else:
            assert collisions == [
                [4, 5],
                [4, 7],
                [4, 8],
                [5, 7],
                [5, 8],
                [7, 8],
            ]

    assert not any(centre["center_kind"] == "satellite" for centre in centres)
    assert len(
        {
            centre["parent_center_positions_zero_based"][0]
            for centre in centres
            if centre["center_kind"] == "free"
        }
    ) == len(free)


def test_independent_jets_kernels_and_exact_cluster_orders(
    packet: dict[str, Any], archive: dict[str, Any]
) -> None:
    values = _closed_coordinates(archive)
    modes = _surviving_modes(archive)
    geometry = archive["geometric_component"]
    minus_k, _ = _independent_jet_matrix(
        packet,
        degree=3,
        multiplicities=[1] * 10,
        values=values,
        modes=modes,
    )
    minus_two_k, _ = _independent_jet_matrix(
        packet,
        degree=6,
        multiplicities=[2] * 10,
        values=values,
        modes=modes,
    )
    assert (minus_k.shape, minus_k.rank(), minus_k.cols - minus_k.rank()) == (
        (10, 10),
        10,
        0,
    )
    assert (
        minus_two_k.shape,
        minus_two_k.rank(),
        minus_two_k.cols - minus_two_k.rank(),
    ) == ((30, 28), 27, 1)

    curves = [
        geometry["physical_boundary"],
        *geometry["noncontracted_nonboundary_facets"],
    ]
    for curve in curves:
        multiplicities = curve["declared_multiplicities_chronological"]
        matrix, monomials = _independent_jet_matrix(
            packet,
            degree=int(curve["degree"]),
            multiplicities=multiplicities,
            values=values,
            modes=modes,
        )
        polynomial = _polynomial(curve["homogeneous_polynomial"])
        coefficients = sp.Matrix(
            [sp.Poly(polynomial, X, Y, Z).coeff_monomial(monomial) for monomial in monomials]
        )
        assert matrix * coefficients == sp.zeros(matrix.rows, 1)
        assert matrix.cols - matrix.rank() == 1
        assert _independent_cluster_orders(
            packet, polynomial, multiplicities, values, modes
        ) == multiplicities


def test_alternate_absolute_factorization_and_nonreduced_control(
    archive: dict[str, Any]
) -> None:
    geometry = archive["geometric_component"]
    curves = [
        geometry["physical_boundary"],
        *geometry["noncontracted_nonboundary_facets"],
    ]
    for curve in curves:
        if curve["degree"] == 1:
            continue
        assert _bcg_absolute_factor_record(
            _polynomial(curve["homogeneous_polynomial"])
        ) == (2, 1, 1)
    assert _bcg_absolute_factor_record((X + Y) ** 2) == (2, 2, 2)


def test_wall_partition_genus_smoothness_and_strong_coble_inputs(
    packet: dict[str, Any], archive: dict[str, Any]
) -> None:
    geometry = archive["geometric_component"]
    rows = {
        int(row["wall_index_zero_based"]): row for row in packet["divisor_table"]
    }
    assert set(rows) == set(range(30))
    contracted = set(
        geometry["contracted_exceptional_facets"]["wall_indices_zero_based"]
    )
    nonboundary = {
        int(curve["wall_index_zero_based"])
        for curve in geometry["noncontracted_nonboundary_facets"]
    }
    boundary_index = int(geometry["physical_boundary"]["wall_index_zero_based"])
    assert contracted | nonboundary | {boundary_index} == set(rows)
    assert not contracted & nonboundary
    assert boundary_index not in contracted | nonboundary

    centres = packet["reverse_blowup_program"]["center_records"]
    for centre in centres:
        wall = rows[int(centre["source_contracted_wall_index_zero_based"])]
        strict = centre[
            "strict_exceptional_coordinates_in_chronological_total_basis"
        ]
        assert wall["exceptional_coefficients_chronological"] == strict
        assert wall["multiplicities_chronological"] == [-value for value in strict]

    for curve in geometry["noncontracted_nonboundary_facets"]:
        row = rows[int(curve["wall_index_zero_based"])]
        assert curve["degree"] == row["plane_degree"]
        assert (
            curve["declared_multiplicities_chronological"]
            == row["multiplicities_chronological"]
        )
        square = curve["degree"] ** 2 - sum(
            value**2 for value in curve["declared_multiplicities_chronological"]
        )
        canonical = -3 * curve["degree"] + sum(
            curve["declared_multiplicities_chronological"]
        )
        assert (square, canonical, (square + canonical) // 2 + 1) == (
            row["self_intersection"],
            row["canonical_intersection"],
            0,
        )

    root = rows[boundary_index]
    boundary = geometry["physical_boundary"]
    assert root["plane_degree"] == 3
    assert root["multiplicities_chronological"] == [1] * 10
    assert boundary["degree"] == 2 * root["plane_degree"]
    assert boundary["declared_multiplicities_chronological"] == [
        2 * value for value in root["multiplicities_chronological"]
    ]
    boundary_square = boundary["degree"] ** 2 - sum(
        value**2 for value in boundary["declared_multiplicities_chronological"]
    )
    boundary_canonical = -3 * boundary["degree"] + sum(
        boundary["declared_multiplicities_chronological"]
    )
    assert (boundary_square, boundary_canonical) == (-4, 2)
    assert (boundary_square + boundary_canonical) // 2 + 1 == 0

    # Over characteristic zero, a geometrically integral Cartier curve has
    # p_a = g(normalization) + sum(delta).  Hence p_a=0 forces a smooth P1.
    assert geometry["minus_K"]["kernel_dimension"] == 0
    assert geometry["minus_two_K"]["kernel_dimension"] == 1
    assert boundary["absolute_factorization"]["absolutely_irreducible"] is True
    assert boundary["absolute_factorization"][
        "sole_nontrivial_factor_multiplicity"
    ] == 1
    assert geometry["strong_Coble"] is True


def test_k3_rank_pushback_and_scope_are_not_same_p_claims(
    archive: dict[str, Any]
) -> None:
    bridge = archive["dimension_bridge"]
    assert bridge["current_compiler_base_field"] == "Q, with geometric conclusion over C"
    assert bridge["current_calibration_picard_rank"] == 11
    assert bridge["current_calibration_hyperbolic_dimension"] == 10
    assert "n<=19 over C" in bridge["general_statement"]
    assert "does not construct" in bridge["odd_characteristic_analogue"]
    assert bridge["does_not_require_same_nef_cone_for_rank_bound"] is True
    assert bridge["current_calibration_not_claimed_to_be_reflection_maximal_here"] is True
    assert "global raw arithmetic-chamber essential surjectivity" in archive[
        "scope_boundary"
    ]["not_proved"]

    # For the smooth double cover branched in B=-2K, pullback is injective
    # (push-pull is multiplication by two) and doubles the form.  Thus its
    # rank is rho(S)=n+1, while a complex K3 has Picard rank at most 20.
    n = sp.Symbol("n", integer=True, nonnegative=True)
    assert sp.solve_univariate_inequality(n + 1 <= 20, n) == (n <= 19)
