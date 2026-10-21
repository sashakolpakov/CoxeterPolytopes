"""Exact checks for the global exceptional-tail attachment no-go."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ATTACHMENT = ROOT / "results/rational_finite_volume_mixed_attachment_lemma.md"
REPORT = ROOT / "results/rational_global_exceptional_tail_bridge_no_go.md"


def bilinear(left, gram, right):
    # Every ambient Gram used below is U + <-1>^m.
    assert gram[0][1] == gram[1][0] == 1
    return (
        Fraction(left[0]) * Fraction(right[1])
        + Fraction(left[1]) * Fraction(right[0])
        - sum(
            Fraction(left[index]) * Fraction(right[index])
            for index in range(2, len(left))
        )
    )


def mat_vec(gram, vector):
    return [
        sum(Fraction(entry) * Fraction(value) for entry, value in zip(row, vector))
        for row in gram
    ]


def exact_rank(matrix):
    work = [[Fraction(entry) for entry in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0]) if work else 0
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not work[row][column]:
                continue
            coefficient = work[row][column]
            work[row] = [
                entry - coefficient * pivot_entry
                for entry, pivot_entry in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def core_tail_gram(length):
    # Total-transform basis c,d,e_1,...,e_m on a blowup of P1xP1.
    size = 2 + length
    gram = [[0] * size for _ in range(size)]
    gram[0][1] = gram[1][0] = 1
    for index in range(2, size):
        gram[index][index] = -1
    return gram


def affine_roots(length):
    size = 2 + length
    roots = []
    # delta_0 = c-e_1-e_2: the strict transform of a square-zero fibre
    # after the first two centres are taken along it.
    delta_zero = [0] * size
    delta_zero[0] = 1
    delta_zero[2] = delta_zero[3] = -1
    roots.append(delta_zero)
    # delta_i=e_i-e_{i+1}, followed by delta_m=e_m.
    for index in range(length - 1):
        root = [0] * size
        root[2 + index] = 1
        root[3 + index] = -1
        roots.append(root)
    endpoint = [0] * size
    endpoint[-1] = 1
    roots.append(endpoint)
    return roots


def test_arbitrary_affine_b_tail_packet_has_exact_null_cycle_and_rank():
    for length in range(2, 41):
        ambient = core_tail_gram(length)
        roots = affine_roots(length)
        physical = [
            [bilinear(left, ambient, right) for right in roots]
            for left in roots
        ]
        assert [physical[i][i] for i in range(length + 1)] == (
            [-2] * length + [-1]
        )

        # Affine B_m graph: delta_0 and delta_1 attach to delta_2,
        # then a chain runs to the short endpoint.
        assert physical[0][1] == 0
        assert physical[0][2] == 1
        assert all(
            physical[index][index + 1] == 1
            for index in range(1, length)
        )
        allowed = {
            (0, 2),
            *((index, index + 1) for index in range(1, length)),
        }
        for i in range(length + 1):
            for j in range(i + 1, length + 1):
                if (i, j) not in allowed:
                    assert physical[i][j] == 0

        null_coefficients = [1, 1] + [2] * (length - 1)
        assert mat_vec(physical, null_coefficients) == [0] * (length + 1)
        assert exact_rank(physical) == length

        # The null divisor is literally the pullback c.
        null_vector = [
            sum(
                coefficient * root[entry]
                for coefficient, root in zip(null_coefficients, roots)
            )
            for entry in range(2 + length)
        ]
        assert null_vector == [1, 0] + [0] * length


def test_canonical_coloring_and_attachment_energy_do_not_grow_with_tail():
    for length in range(2, 41):
        ambient = core_tail_gram(length)
        roots = affine_roots(length)
        # K_{P1xP1}=-2c-2d and each blowup adds e_i.
        canonical = [-2, -2] + [1] * length
        host_colors = [
            2 * bilinear(canonical, ambient, root) for root in roots
        ]
        host_norms = [
            2 * bilinear(root, ambient, root) for root in roots
        ]
        assert host_colors == [0] * length + [-2]
        assert host_norms == [-4] * length + [-2]

        mixed_projection = roots[0][2:]
        assert mixed_projection[:2] == [-1, -1]
        assert mixed_projection[2:] == [0] * (length - 2)
        assert sum(value * value for value in mixed_projection) == 2

        # One signed-permutation orbit spans every tail coordinate:
        # flip one of the two nonzero coordinates, take the difference,
        # then permute that coordinate.
        vector = mixed_projection
        flipped = vector.copy()
        flipped[0] *= -1
        difference = [left - right for left, right in zip(vector, flipped)]
        assert difference == [-2] + [0] * (length - 1)
        orbit_basis = [
            [int(column == row) * difference[0] for column in range(length)]
            for row in range(length)
        ]
        assert all(
            row[index] == -2 and sum(value != 0 for value in row) == 1
            for index, row in enumerate(orbit_basis)
        )


def test_all_three_standard_host_canonical_colors_have_the_claimed_scale():
    # Adjunction gives K.C=-2-C^2 for every smooth rational curve C.
    physical_squares = {"minus_two": -2, "minus_one": -1, "minus_four": -4}
    roots = {
        "minus_two": Fraction(1),
        "minus_one": Fraction(1),
        "minus_four": Fraction(1, 2),
    }
    host_norms = {
        name: 2 * square * roots[name] ** 2
        for name, square in physical_squares.items()
    }
    host_colors = {
        name: 2 * (-2 - square) * roots[name]
        for name, square in physical_squares.items()
    }
    assert host_norms == {"minus_two": -4, "minus_one": -2, "minus_four": -2}
    assert host_colors == {"minus_two": 0, "minus_one": -2, "minus_four": 2}


def test_endpoint_reflection_does_not_preserve_the_canonical_coloring():
    for length in (2, 3, 8, 20):
        ambient = core_tail_gram(length)
        roots = affine_roots(length)
        canonical = [-2, -2] + [1] * length
        endpoint = roots[-1]
        endpoint_square = bilinear(endpoint, ambient, endpoint)
        coefficient = (
            -2 * bilinear(canonical, ambient, endpoint) / endpoint_square
        )
        reflected_canonical = [
            Fraction(value) + coefficient * root_value
            for value, root_value in zip(canonical, endpoint)
        ]
        assert coefficient == -2
        assert reflected_canonical != canonical

        mixed = roots[0]
        reflected_mixed = [
            Fraction(value)
            - 2
            * bilinear(mixed, ambient, endpoint)
            / endpoint_square
            * root_value
            for value, root_value in zip(mixed, endpoint)
        ]
        # For length>2 the chosen mixed root is endpoint-orthogonal; at
        # length two it is not.  In either case the group action has no
        # theorem preserving the effective/adjunction coloring globally.
        if length == 2:
            assert 2 * bilinear(canonical, ambient, reflected_mixed) != 0


def test_cusp_reflection_integrality_has_only_short_and_long_mirrors():
    # In U(2)+<-2>^m, a primitive root perpendicular to c has
    # r=a*c+sum b_i e_i, r^2=-2d, d=sum b_i^2.  Integrality requires
    # d | 2a and d | 2b_i.  Exhaust a hostile box; the symbolic size
    # argument in the report excludes all larger d.
    surviving_d = set()
    for b1 in range(-6, 7):
        for b2 in range(-6, 7):
            for b3 in range(-6, 7):
                vector = (b1, b2, b3)
                d = sum(value * value for value in vector)
                if d == 0:
                    continue
                if all((2 * value) % d == 0 for value in vector):
                    surviving_d.add(d)
    assert surviving_d == {1, 2, 4}
    # d=4 consists only of a single coefficient +/-2 and gives a
    # nonprimitive multiple of a d=1 mirror once d|2a is imposed.
    d4 = []
    for vector in (
        (a, b, c)
        for a in range(-3, 4)
        for b in range(-3, 4)
        for c in range(-3, 4)
    ):
        if sum(value * value for value in vector) == 4 and all(
            (2 * value) % 4 == 0 for value in vector
        ):
            d4.append(vector)
    assert all(sorted(abs(value) for value in vector) == [0, 0, 2] for vector in d4)


def test_attachment_scope_is_explicit_and_typo_is_repaired():
    attachment = " ".join(ATTACHMENT.read_text(encoding="utf-8").split())
    assert r"t>0,\quad u\in\mathbf H(V_0),\quad" in attachment
    assert r"\boxed{U=E}" in attachment
    assert "not their raw linear span and not their number" in attachment
    assert "one mixed facet have nonzero tail projection" in attachment

    if REPORT.exists():
        report = " ".join(REPORT.read_text(encoding="utf-8").split())
        assert "does not prove an unbounded reflection-maximal family" in report
        assert "individual mixed-facet complexity cannot bound the tail" in report
        assert "total marked-facet degree" in report
        assert "finite numerical contraction orbit" in report
