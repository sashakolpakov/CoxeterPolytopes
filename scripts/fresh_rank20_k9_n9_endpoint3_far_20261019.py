#!/usr/bin/env python3
"""ONE bounded necessary-condition batch: K9, N=9, far endpoint3.

SOURCE PREPARED, NOT EXECUTED. Root must first review and then schedule
under the shared serial lock/resource wrapper (<=59s, <=1GiB RSS).
This script has an additional 55s deadline; it never starts workers.

No other field, N, or contact role is accepted. All arithmetic decisions
are exact Fractions/Sturm signs. A complete nonempty hit set is research
data, NOT an absence proof or an assertion of geometric realizability.
An independent full-box replay remains required after any future run.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
from itertools import product
import json
from pathlib import Path
import time

import exact_totally_real_field_20261019 as E

ROOT = Path(__file__).resolve().parents[1]
SELF = Path(__file__).resolve()
HELPER = Path(E.__file__).resolve()
HEIGHT = ROOT / "results/arithmetic_rank20_isolated_long_dotted_height_bounds_20261019.md"
DESIGN = ROOT / "results/arithmetic_rank20_isolated_long_exact_enumeration_design_20261019.md"
EXPECTED = {
    HELPER: "fbd896a6160f232e0f6b4e8a715c990b05498d74be5c846e91aa699c95d598d2",
    HEIGHT: "9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516",
    DESIGN: "b1dcfe5e74d1372d2250cdbfc3f5dc256bacc7b704aada73760cdad5e3e2dc26",
}
N = 9
RANK = N + 2
F = E.trim((1, -3, 0, 1))
PHYSICAL = 2
MAX_TAILS = 100_000
MAX_HITS = 128
MAX_CERT_BYTES = 512 * 1024
WALL_SECONDS = 55


class Incomplete(RuntimeError):
    pass


def stream_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_map() -> dict[str, str]:
    return {str(path.relative_to(ROOT)): stream_hash(path)
            for path in (SELF, HELPER, HEIGHT, DESIGN)}


def ceiling(q: Q) -> int:
    return -((-q.numerator) // q.denominator)


def floor(q: Q) -> int:
    return q.numerator // q.denominator


def serial_matrix(matrix: list[list[Q]]) -> list[list[str]]:
    return [[str(v) for v in row] for row in matrix]


def field_definition(guard) -> E.RealField:
    # Cubic irreducibility: a rational root of this monic constant-one
    # polynomial would be +/-1; neither is a root. This is separate from
    # the real-root/Sturm validation performed by RealField.
    if F != E.trim((1, -3, 0, 1)) or any(E.peval(F, Q(v)) == 0
                                       for v in (-1, 1)):
        raise ArithmeticError("fixed K9 minimal-polynomial proof failed")
    roots = [E.Interval(Q(-2), Q(-3, 2)),
             E.Interval(Q(0), Q(1, 2)),
             E.Interval(Q(3, 2), Q(7, 4))]
    return E.RealField(F, roots, guard)


def calibrations(field: E.RealField) -> dict:
    """Fixed exact facts, run only inside the future scheduled batch."""
    cases = [
        (E.pmul(E.trim((Q(-1, 3), 1)), E.trim((Q(-2, 3), 1))), 2),
        (E.pmul(E.trim((Q(-1, 2), 1)), E.trim((Q(-1, 2), 1))), 1),
        (E.trim((1, 0, 1)), 0),
    ]
    for polynomial, count in cases:
        for scale in (Q(1), Q(-7)):
            chain = E.sturm_chain(E.pscale(polynomial, scale))
            if E.root_count(chain, Q(0), Q(1)) != count:
                raise ArithmeticError("fixed Sturm calibration failed")
    if field.inverse(field.theta) != field.element((3, 0, -1)):
        raise ArithmeticError("known K9 inverse calibration failed")
    if [field.trace(field.power(field.theta, j)) for j in (0, 1, 2, 3)] != [3, 0, 6, -3]:
        raise ArithmeticError("known K9 trace calibration failed")
    theta_minus_one = field.sub(field.theta, field.one)
    if [field.sign(theta_minus_one, i) for i in range(3)] != [-1, -1, 1]:
        raise ArithmeticError("known K9 embedding-sign calibration failed")
    if field.generated_degree(field.theta) != 3 or field.generated_degree(field.one) != 1:
        raise ArithmeticError("known K9 generated-degree calibration failed")
    return {"Sturm_counts_including_repeated_root_and_negative_scaling": 6,
            "known_K9_inverse_trace_sign_and_degree_facts": "checked exactly"}


def boxes(field: E.RealField, certificate: dict, guard) -> list[dict]:
    certificate["root_intervals_initial"] = field.serial_roots()
    field.refine_width(Q(1, 4096))
    certificate["root_intervals_used_for_boxes"] = field.serial_roots()
    certificate["sturm_chain"] = [E.serial_element(p) for p in field.chain]
    trace_matrix, dual = field.trace_dual()
    certificate["trace_matrix"] = serial_matrix(trace_matrix)
    certificate["trace_dual_coordinates"] = [E.serial_element(v) for v in dual]
    certificate["traces_of_power_basis"] = [str(v) for v in field.trace_powers[:3]]

    # Actual nonphysical values obey |9A_i+1|<sqrt10. Since 10<361/36,
    # this CLOSED rational OUTER interval contains both actual values.
    if not Q(10) < Q(19, 6) ** 2:
        raise ArithmeticError("rational sqrt10 outer bound failed")
    other = E.Interval(Q(-25, 54), Q(13, 54))
    physical = E.Interval(Q(2), Q(24))  # A<24 checked strictly later.
    trace_interval = physical + other + other
    trace_lo, trace_hi = ceiling(trace_interval.lo), floor(trace_interval.hi)
    certificate["physical_A_outer_interval"] = physical.serial()
    certificate["nonphysical_A_outer_intervals"] = [other.serial(), other.serial()]
    certificate["trace_outer_interval"] = trace_interval.serial()
    certificate["integer_trace_range"] = [trace_lo, trace_hi]
    certificate["trace_dual_formula"] = (
        "c_j=T*lambda_j(theta_phys)+sum_i_nonphys "
        "A_i*(lambda_j(theta_i)-lambda_j(theta_phys))")
    output = []
    for trace in range(trace_lo, trace_hi + 1):
        guard()
        coefficient_boxes = []
        interval_boxes = []
        for dual_element in dual:
            physical_dual = field.enclose(dual_element, PHYSICAL)
            value = physical_dual.scaled(trace)
            for i in (0, 1):
                value += other * (field.enclose(dual_element, i) - physical_dual)
            interval_boxes.append(value.serial())
            coefficient_boxes.append([ceiling(value.lo), floor(value.hi)])
        tail_count = 1
        for lo, hi in coefficient_boxes[1:]:
            tail_count *= max(0, hi - lo + 1)
        output.append({"trace": trace,
                       "closed_coefficient_intervals": interval_boxes,
                       "integer_coefficient_boxes": coefficient_boxes,
                       "raw_tail_vectors": tail_count})
    certificate["trace_boxes"] = output
    certificate["raw_tail_vector_count"] = sum(row["raw_tail_vectors"] for row in output)
    return output


def classify(field: E.RealField, coordinates: list[int], hits: list[dict]) -> str:
    A = field.element(coordinates)
    constant = lambda q: field.element((Q(q),))
    if field.sign(field.sub(A, constant(2)), PHYSICAL) < 0:
        return "A_below_2"
    if field.sign(field.sub(constant(24), A), PHYSICAL) <= 0:
        return "A_not_below_24"

    Z = field.add(field.scale(A, N), field.one)
    Z2 = field.mul(Z, Z)
    for i in (0, 1):
        if field.sign(field.sub(constant(N + 1), Z2), i) <= 0:
            return "nonphysical_proper_block_not_PD"
    qc = field.add(constant(1 - Q(N + 1, 2 * N)), field.scale(Z2, Q(1, 2 * N)))
    if field.sign(field.sub(qc, field.one), PHYSICAL) <= 0:
        return "physical_old_plus_c_not_negative"

    C = field.scale(field.sub(field.scale(Z2, 2), constant(RANK * (N + 1))), Q(1, N))
    # C=2N*A^2+4A-N-3 is integral for EVERY enumerated integral A.
    expanded_C = field.sub(field.add(field.scale(field.mul(A, A), 2 * N),
                                     field.scale(A, 4)), constant(N + 3))
    if C != expanded_C or not field.is_integral_power_basis(C):
        raise ArithmeticError("predicted integral C identity failed")
    if C == field.zero:
        return "zero_C_incompatible_with_nonzero_rhs"
    delta = field.scale(field.inverse(C), 4 * (N + 1))
    if not field.is_integral_power_basis(delta):
        return "nonintegral_delta"
    if field.sign(delta, PHYSICAL) <= 0:
        return "nonpositive_physical_delta"
    z = field.add(delta, constant(4))
    for i in (0, 1):
        if field.sign(z, i) <= 0:
            return "nonpositive_nonphysical_z"
        if field.sign(field.sub(constant(Q(8, RANK)), z), i) <= 0:
            return "nonphysical_z_not_below_8_over_r"
    R = field.scale(Z2, Q(1, 2 * N * (N + 1)))
    if field.mul(field.sub(R, constant(Q(RANK, 4 * N))), delta) != field.one:
        raise ArithmeticError("original corank identity failed")
    if field.mul(C, delta) != constant(4 * (N + 1)):
        raise ArithmeticError("C*delta scalar identity failed")
    if field.generated_degree(z) != field.degree:
        return "z_does_not_generate_actual_ground_field"

    # A surviving vector is ONLY a necessary-condition hit. Save all of
    # its exact data; never infer existence of a compact ambient polytope.
    hits.append({"A_coordinates": coordinates,
                 "A_trace": str(field.trace(A)),
                 "Z_coordinates": E.serial_element(Z),
                 "C_coordinates": E.serial_element(C),
                 "delta_coordinates": E.serial_element(delta),
                 "z_coordinates": E.serial_element(z),
                 "qc_coordinates": E.serial_element(qc),
                 "R_coordinates": E.serial_element(R),
                 "status": "UNRESOLVED_NECESSARY_CONDITION_HIT"})
    if len(hits) > MAX_HITS:
        raise Incomplete("unresolved-hit storage cap reached; absence blocked")
    return "unresolved_hit"


def prepare_and_enumerate(certificate: dict, guard) -> None:
    field = field_definition(guard)
    certificate["fixed_calibrations"] = calibrations(field)
    certificate["field"] = {
        "name": "K9", "degree": 3,
        "minimal_polynomial_ascending": E.serial_element(F),
        "irreducibility_input": "monic cubic rational-root theorem; roots +/-1 excluded",
        "theta": "2*cos(2*pi/9)", "physical_embedding_index": PHYSICAL,
        "embedding_order": "increasing real roots; physical is largest",
        "integral_basis": ["1", "theta", "theta^2"],
        "external_integer_ring_input": "O_K9=Z[2*cos(2*pi/9)]; maximal-real-cyclotomic integer-ring theorem",
    }
    all_boxes = boxes(field, certificate, guard)
    if certificate["raw_tail_vector_count"] > MAX_TAILS:
        raise Incomplete("complete trace boxes exceed the fixed visit budget")
    counters = Counter()
    certificate["counters"] = counters
    certificate["unresolved_hits"] = []
    visited_digest = hashlib.sha256()
    certificate["visited_tail_vectors"] = 0
    certificate["completed_tail_vectors"] = 0
    try:
        for row in all_boxes:
            coefficient_boxes = row["integer_coefficient_boxes"]
            tails = [range(lo, hi + 1) for lo, hi in coefficient_boxes[1:]]
            for tail in product(*tails):
                guard()
                certificate["visited_tail_vectors"] += 1
                trace = row["trace"]
                certificate["partial_current_vector"] = {"trace": trace, "tail": list(tail)}
                c0 = (Q(trace) - sum((Q(v) * field.trace_powers[j + 1]
                                      for j, v in enumerate(tail)), Q(0))) / field.degree
                if c0.denominator != 1:
                    outcome = "nonintegral_trace_recovered_c0"
                elif not coefficient_boxes[0][0] <= c0 <= coefficient_boxes[0][1]:
                    outcome = "trace_recovered_c0_outside_its_outer_box"
                else:
                    coordinates = [int(c0), *tail]
                    if field.trace(field.element(coordinates)) != trace:
                        raise ArithmeticError("exact trace recovery failed")
                    outcome = classify(field, coordinates, certificate["unresolved_hits"])
                counters[outcome] += 1
                visited_digest.update(json.dumps([trace, list(tail), outcome],
                                                  separators=(",", ":")).encode("ascii") + b"\n")
                certificate["completed_tail_vectors"] += 1
                certificate.pop("partial_current_vector", None)
        if certificate["completed_tail_vectors"] != certificate["raw_tail_vector_count"]:
            raise ArithmeticError("enumeration did not cover every tail in every trace box")
        certificate["complete"] = True
        certificate["has_unresolved_hits"] = bool(certificate["unresolved_hits"])
        certificate["status"] = ("COMPLETE_WITH_UNRESOLVED_NECESSARY_HITS"
                                 if certificate["has_unresolved_hits"]
                                 else "COMPLETE_WITH_ZERO_NECESSARY_HITS")
        certificate["necessary_hit_set_empty"] = not certificate["has_unresolved_hits"]
    finally:
        certificate["visited_outcome_sha256"] = visited_digest.hexdigest()
        certificate["root_intervals_after_sign_checks"] = field.serial_roots()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="new JSON path under results; existing evidence is never overwritten")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.parent != (ROOT / "results").resolve() or output.suffix != ".json":
        raise SystemExit("output must be a new JSON file directly under results")
    if output.exists():
        raise SystemExit("refusing to overwrite existing evidence")
    started = time.monotonic()

    def guard() -> None:
        if time.monotonic() - started >= WALL_SECONDS:
            raise Incomplete("55-second internal deadline reached")

    certificate = {
        "schema_version": 1, "status": "INCOMPLETE",
        "complete": False, "has_unknown_or_unresolved_work": True,
        "scope": "ONE K9, N=9, far-endpoint3 isolated-d, spectator-orthogonal dotted-x arithmetic row",
        "N": N, "r": RANK, "contact_type": "endpoint3_far", "integral_variable": "A=2x",
        "hypotheses": "actual-quad corank/cyclic-integrality interfaces and pinned written height theorem; external integral-ring input explicit below",
        "independent_full_box_replay": "NOT YET IMPLEMENTED OR RUN",
        "absence_claim_requires_independent_replay": True,
        "resource_contract": {"serial_root_scheduled_wrapper_required": True,
                              "external_max_seconds": 59, "external_max_RSS_MiB": 1024,
                              "internal_deadline_seconds": WALL_SECONDS,
                              "max_raw_tail_vectors": MAX_TAILS},
    }
    before = {}
    try:
        before = source_map()
        certificate["sources"] = before
        for path, expected in EXPECTED.items():
            if before[str(path.relative_to(ROOT))] != expected:
                raise Incomplete("source dependency differs from its reviewed/design pin: " + path.name)
        prepare_and_enumerate(certificate, guard)
    except Incomplete as error:
        certificate["status"] = "INCOMPLETE_OR_DEFERRED"
        certificate["incomplete_reason"] = str(error)
        certificate["complete"] = False
    except Exception as error:
        certificate["status"] = "INCOMPLETE_INTERNAL_ERROR"
        certificate["incomplete_reason"] = type(error).__name__ + ": " + str(error)[:1000]
        certificate["complete"] = False
    try:
        after = source_map()
        certificate["sources_after"] = after
        if before != after:
            certificate["status"] = "INCOMPLETE_SOURCE_CHANGED"
            certificate["complete"] = False
    except Exception as error:
        certificate["status"] = "INCOMPLETE_SOURCE_RECHECK_FAILED"
        certificate["complete"] = False
        certificate["source_recheck_error"] = type(error).__name__ + ": " + str(error)[:1000]
    certificate["has_unknown_or_unresolved_work"] = (not certificate["complete"] or
                                                     bool(certificate.get("unresolved_hits")))
    if not certificate["complete"]:
        certificate["necessary_hit_set_empty"] = False
    certificate["elapsed_seconds"] = round(time.monotonic() - started, 6)
    data = (json.dumps(certificate, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if len(data) > MAX_CERT_BYTES:
        certificate = {"status": "INCOMPLETE_CERTIFICATE_OVERFLOW", "complete": False,
                       "has_unknown_or_unresolved_work": True,
                       "sources": before, "oversized_certificate_bytes": len(data),
                       "oversized_certificate_sha256": hashlib.sha256(data).hexdigest(),
                       "first_unresolved_hits": certificate.get("unresolved_hits", [])[:2],
                       "partial_current_vector": certificate.get("partial_current_vector"),
                       "unknown_hit_record": "certificate output cap exceeded; no absence claim"}
        data = (json.dumps(certificate, indent=2, sort_keys=True) + "\n").encode("utf-8")
    # Exclusive creation prevents accidental mutation of existing/frozen data.
    with output.open("xb") as handle:
        handle.write(data)
    print(json.dumps({"status": certificate["status"],
                      "complete": certificate["complete"],
                      "unresolved_hits": len(certificate.get("unresolved_hits", [])),
                      "certificate": str(output.relative_to(ROOT)),
                      "certificate_sha256": hashlib.sha256(data).hexdigest()}, sort_keys=True))
    return 0 if certificate["complete"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
