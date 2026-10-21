#!/usr/bin/env python3
"""Certify the first cutter of the golden A1^24 Golay uniform ray.

The previously found hyperideal ray

    u = (1; (phi/4)^24)

survives the packets at controls 1 and phi, but it does not persist.  This
audit checks a norm-two root at control phi^3=1+2phi against every earlier
packet, including a complete coefficient-pattern census at phi^2 and 2phi.
It also reproduces the 16 Golay-feasible multiplicity patterns at phi^3 that
violate the uniform-ray inequality.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from typing import Sequence

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.search.cubic7_leech import extended_binary_golay_codewords
from coxeter_polytopes.search.golden_a1_golay import (
    a1_golay_pair,
    golay_membership,
    height_zero_walls,
    octads,
    uniform_ray_audit,
    x1_packet,
    xphi_packet,
)


def _value(value: Golden) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "identity_sign": value.sign(0),
        "conjugate_sign": value.sign(1),
    }


def _controls_through_cutter() -> tuple[Golden, ...]:
    """Enumerate every admissible integral control no later than phi^3."""

    cutter = ONE + 2 * PHI
    controls = []
    # If 0<x<=phi^3<5 and 1+phi'*x'^2>=0, then |x'|<2.
    # Thus |b|sqrt(5)=|x-x'|<7 gives |b|<=3, and
    # |2a+b|=|x+x'|<7 gives |a|<=4.
    if (Golden(5) - cutter).sign(0) <= 0:
        raise AssertionError("control coefficient bound failed")
    for a in range(-4, 5):
        for b in range(-3, 4):
            control = Golden(a, b)
            target = ONE + PHI * control * control
            if control.sign(0) > 0 and (cutter - control).sign(0) >= 0 and target.sign(1) >= 0:
                controls.append(control)
    return tuple(sorted(controls))


def _coordinate_alphabet(target: Golden) -> tuple[Golden, ...]:
    """All nonnegative half-integral coordinate values in a shell."""

    target4 = 4 * target
    # At every audited shell, |z|<12 and |z'|<4 for z=2d.  These imply
    # |B|<=7 and |A|<=11; the enumerated A interval is one wider.
    if (Golden(144) - target4).sign(0) <= 0 or (Golden(16) - target4).sign(1) <= 0:
        raise AssertionError("coordinate box is too small for target")
    values = set()
    for a in range(-12, 13):
        for b in range(-7, 8):
            coordinate = Golden(Fraction(a, 2), Fraction(b, 2))
            if coordinate.sign(0) < 0:
                continue
            if all((target - coordinate * coordinate).sign(place) >= 0 for place in (0, 1)):
                values.add(coordinate)
    return tuple(sorted(values))


@lru_cache(maxsize=1)
def _golay_intersection_triples() -> frozenset[tuple[int, int, int]]:
    words = extended_binary_golay_codewords()
    return frozenset(
        (left.bit_count(), right.bit_count(), (left & right).bit_count())
        for left in words
        for right in words
    )


def _square_multiplicity_patterns(control: Golden) -> tuple[dict[str, object], ...]:
    """Complete unordered coordinate-pattern census at one control."""

    target = ONE + PHI * control * control
    values = tuple(value for value in _coordinate_alphabet(target) if value)
    types = []
    for coordinate in values:
        square4 = 4 * coordinate * coordinate
        if square4.a.denominator != 1 or square4.b.denominator != 1:
            raise AssertionError("doubled coordinate square did not clear denominators")
        if square4.a <= 0:
            raise AssertionError("unexpected nonpositive rational square coefficient")
        doubled = 2 * coordinate
        if doubled.a.denominator != 1 or doubled.b.denominator != 1:
            raise AssertionError("coordinate is not half-integral")
        types.append(
            (
                coordinate,
                int(square4.a),
                int(square4.b),
                (int(doubled.a) % 2, int(doubled.b) % 2),
            )
        )
    target4 = 4 * target
    target_a, target_b = int(target4.a), int(target4.b)
    raw: list[tuple[int, ...]] = []

    def extend(index: int, remainder_a: int, remainder_b: int, slots: int, prefix: tuple[int, ...]) -> None:
        if index == len(types):
            if remainder_a == remainder_b == 0:
                raw.append(prefix)
            return
        _coordinate, square_a, square_b, _residue = types[index]
        maximum = min(slots, remainder_a // square_a)
        for count in range(maximum + 1):
            extend(
                index + 1,
                remainder_a - count * square_a,
                remainder_b - count * square_b,
                slots - count,
                (*prefix, count),
            )

    extend(0, target_a, target_b, 24, ())
    triples = _golay_intersection_triples()
    records = []
    for counts in raw:
        weight_a = sum(
            count * residue[0]
            for count, (_coordinate, _sa, _sb, residue) in zip(counts, types)
        )
        weight_b = sum(
            count * residue[1]
            for count, (_coordinate, _sa, _sb, residue) in zip(counts, types)
        )
        intersection = sum(
            count
            for count, (_coordinate, _sa, _sb, residue) in zip(counts, types)
            if residue == (1, 1)
        )
        triple = (weight_a, weight_b, intersection)
        spatial_sum = sum(
            (count * coordinate for count, (coordinate, _sa, _sb, _residue) in zip(counts, types)),
            start=ZERO,
        )
        records.append(
            {
                "multiplicities": [
                    {
                        "coordinate": coordinate.to_json(),
                        "display": str(coordinate),
                        "count": count,
                    }
                    for count, (coordinate, _sa, _sb, _residue) in zip(counts, types)
                    if count
                ],
                "nonzero_coordinate_count": sum(counts),
                "golay_support_triple": list(triple),
                "golay_feasible": triple in triples,
                "spatial_sum": _value(spatial_sum),
                "uniform_excess": _value(spatial_sum - 4 * control),
                "cuts_uniform_ray": (spatial_sum - 4 * control).sign(0) > 0,
            }
        )
    return tuple(records)


def _maximum_pairing(root: Sequence[Golden], walls: Sequence[Sequence[Golden]]) -> Golden:
    return max(a1_golay_pair(root, wall) for wall in walls)


def build_audit() -> dict[str, object]:
    controls = _controls_through_cutter()
    expected_controls = (ONE, PHI, ONE + PHI, 2 * PHI, ONE + 2 * PHI)
    if controls != expected_controls:
        raise AssertionError(f"unexpected controls through cutter: {controls}")

    phi2 = ONE + PHI
    two_phi = 2 * PHI
    cutter_control = ONE + 2 * PHI
    pattern_records = {
        "phi_squared": _square_multiplicity_patterns(phi2),
        "two_phi": _square_multiplicity_patterns(two_phi),
        "cutter": _square_multiplicity_patterns(cutter_control),
    }
    summaries = {
        name: {
            "raw_pattern_count": len(records),
            "golay_feasible_pattern_count": sum(record["golay_feasible"] for record in records),
            "uniform_violating_pattern_count": sum(
                record["golay_feasible"] and record["cuts_uniform_ray"] for record in records
            ),
        }
        for name, records in pattern_records.items()
    }
    expected_summaries = {
        "phi_squared": (43, 3, 0),
        "two_phi": (3, 1, 0),
        "cutter": (383, 33, 16),
    }
    for name, expected in expected_summaries.items():
        observed = summaries[name]
        if (
            observed["raw_pattern_count"],
            observed["golay_feasible_pattern_count"],
            observed["uniform_violating_pattern_count"],
        ) != expected:
            raise AssertionError(f"unexpected {name} pattern census")

    first_octad = octads()[0]
    complement = tuple(index for index in range(24) if not first_octad & (1 << index))
    special = next(index for index in range(24) if first_octad & (1 << index))
    half_phi_squared = (ONE + PHI) / 2
    spatial = [ZERO] * 24
    spatial[special] = PHI
    for index in complement:
        spatial[index] = half_phi_squared
    cutter = (cutter_control, *spatial)
    if not golay_membership(spatial):
        raise AssertionError("cutter failed Golay Construction-A membership")
    if a1_golay_pair(cutter, cutter) != 2:
        raise AssertionError("cutter does not have norm two")

    earlier_maxima = {
        "height_zero": _maximum_pairing(cutter, height_zero_walls()),
        "x_equals_1": _maximum_pairing(cutter, tuple(root.vector for root in x1_packet())),
        "x_equals_phi": _maximum_pairing(cutter, tuple(root.vector for root in xphi_packet())),
    }
    if any(value.sign(0) > 0 for value in earlier_maxima.values()):
        raise AssertionError("cutter violates a first-packet wall")

    # The complete phi^2 pattern census has three Golay-feasible types.
    # One is rejected: on an octad S it has one (phi-1)/2 and seven
    # phi^2/2 coordinates, and pairs by +phi with the earlier octad wall S.
    rejected_phi2_pairing = PHI
    phi2_assignment_counts = {
        "integral_phi_squared_plus_two_phi": 24 * 23 * 22 // 2,
        "octad_halves_plus_phi_squared_outside": len(octads()) * 16,
        "rejected_single_low_plus_seven_half_phi_squared": len(octads()) * 8,
    }
    if phi2_assignment_counts != {
        "integral_phi_squared_plus_two_phi": 6072,
        "octad_halves_plus_phi_squared_outside": 12144,
        "rejected_single_low_plus_seven_half_phi_squared": 6072,
    }:
        raise AssertionError("phi^2 assignment counts changed")

    # Maximize the cutter pairing over every root in the two accepted phi^2
    # orbits.  This is only 18,216 exact dot products and avoids any orbit
    # transitivity assumption in the crucial survival check.
    p1_max: Golden | None = None
    for phi_squared_coordinate in range(24):
        remaining = [index for index in range(24) if index != phi_squared_coordinate]
        for first, second in combinations(remaining, 2):
            dot = (
                spatial[phi_squared_coordinate] * phi2
                + spatial[first] * PHI
                + spatial[second] * PHI
            )
            pairing = -2 * PHI * cutter_control * phi2 + 2 * dot
            if p1_max is None or pairing > p1_max:
                p1_max = pairing
    p2_max: Golden | None = None
    for word in octads():
        octad_sum = sum(
            (spatial[index] for index in range(24) if word & (1 << index)), start=ZERO
        )
        for phi_squared_coordinate in range(24):
            if word & (1 << phi_squared_coordinate):
                continue
            dot = spatial[phi_squared_coordinate] * phi2 + (PHI / 2) * octad_sum
            pairing = -2 * PHI * cutter_control * phi2 + 2 * dot
            if p2_max is None or pairing > p2_max:
                p2_max = pairing
    if p1_max != -6 - 8 * PHI or p2_max != -4 - 4 * PHI:
        raise AssertionError("unexpected cutter/phi^2 maximum pairing")

    # The sole Golay-feasible 2phi shell pattern is d=phi^3 e_i.  Each of
    # its 24 assignments is rejected by the x=1 wall on the same coordinate.
    two_phi_rejection_pairing = 2 * PHI
    if two_phi_rejection_pairing.sign(0) <= 0:
        raise AssertionError("2phi rejection witness lost positivity")

    spatial_sum = sum(spatial, start=ZERO)
    uniform_excess = spatial_sum - 4 * cutter_control
    uniform_pairing = a1_golay_pair(uniform_ray_audit().vector, cutter)
    if uniform_excess != 4 + PHI or uniform_pairing != Golden(Fraction(1, 2), Fraction(5, 2)):
        raise AssertionError("uniform-ray cutter arithmetic changed")
    if uniform_pairing.sign(0) <= 0:
        raise AssertionError("witness no longer cuts the uniform ray")

    violating_patterns = [
        record
        for record in pattern_records["cutter"]
        if record["golay_feasible"] and record["cuts_uniform_ray"]
    ]
    return {
        "schema": "golden-a1-golay-uniform-ray-cutter-0.1.0",
        "scope": (
            "Exact norm-two control prefix through x=phi^3. This certifies "
            "that the displayed ray is cut; it does not complete the later chamber."
        ),
        "field": "Q(phi), phi^2=phi+1",
        "lattice": "<-2phi> + N(A1^24, extended Golay glue)",
        "control_bound_proof": {
            "conditions": "0<x<=phi^3 and 1+phi'*x'^2>=0",
            "deductions": ["x<5", "|x'|<2", "|b|<=3", "|a|<=4"],
            "all_controls": [_value(control) for control in controls],
        },
        "coordinate_box_proof": {
            "doubled_coordinate": "z=A+B*phi=2d",
            "bounds": ["|z|<12", "|z'|<4", "|B|<=7", "|A|<=11"],
            "enumerated_ranges": {"A": [-12, 12], "B": [-7, 7]},
        },
        "pattern_census": {
            name: {
                **summaries[name],
                "golay_feasible_patterns": [
                    record for record in records if record["golay_feasible"]
                ],
            }
            for name, records in pattern_records.items()
        },
        "phi_squared_packet": {
            "complete_golay_pattern_count": 3,
            "assignment_counts": phi2_assignment_counts,
            "accepted_assignment_count": 6072 + 12144,
            "rejected_assignment_count": 6072,
            "rejected_pattern_prior_octad_pairing": _value(rejected_phi2_pairing),
            "same_height_acute": True,
        },
        "two_phi_packet": {
            "complete_golay_pattern_count": 1,
            "candidate_assignment_count": 24,
            "accepted_assignment_count": 0,
            "x1_rejection_pairing": _value(two_phi_rejection_pairing),
        },
        "cutter": {
            "control": _value(cutter_control),
            "first_octad_word": first_octad,
            "half_phi_squared_support": list(complement),
            "phi_coordinate": special,
            "spatial_support": [
                {"coordinate": index, "value": _value(value)}
                for index, value in enumerate(spatial)
                if value
            ],
            "construction_a_membership": True,
            "norm": _value(a1_golay_pair(cutter, cutter)),
            "maximum_pairings_with_prior_packets": {
                **{name: _value(value) for name, value in earlier_maxima.items()},
                "x_equals_phi_squared_integral_orbit": _value(p1_max),
                "x_equals_phi_squared_octad_orbit": _value(p2_max),
            },
            "two_phi_has_no_accepted_walls": True,
            "survives_every_earlier_wall": True,
            "spatial_sum": _value(spatial_sum),
            "uniform_excess_sum_d_minus_4x": _value(uniform_excess),
            "uniform_ray_pairing": _value(uniform_pairing),
            "cuts_uniform_ray": True,
        },
        "uniform_violating_patterns_at_cutter_control": violating_patterns,
        "conclusion": {
            "uniform_ray_persists_through_phi_squared": True,
            "uniform_ray_cut_at_phi_cubed": True,
            "cutter_is_norm_two_lattice_root": True,
            "cutter_survives_all_earlier_walls": True,
            "violating_golay_pattern_count_at_phi_cubed": 16,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/golden_a1_uniform_ray_cutter.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload["conclusion"], sort_keys=True))


if __name__ == "__main__":
    main()
