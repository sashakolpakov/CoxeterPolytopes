"""Rational PSD separation for the fixed global model, without diagram pricing.

Every added row is v^T M v >= 0 for one original PSD block and a rational
vector v. A certificate for the resulting LP is lifted to the original SDP
by summing its nonnegative multipliers times v v^T in each block.
"""
from fractions import Fraction as F
import time

from . import conic


def quadratic(entries, direction):
    terms = []
    for i, a in direction:
        terms.append(conic.times(a*a, entries[i][i]))
        for j, b in direction:
            if j < i:
                terms.append(conic.times(2*a*b, entries[i][j]))
    return conic.plus(*terms)


def outer_LP(primal):
    blocks, directions = [], []
    for b, block in enumerate(primal["blocks"]):
        M, name = block["entries"], block["name"]
        for i in range(len(M)):
            v = [(i, F(1))]
            blocks.append(dict(name=name+f"_diag{i}", entries=[[quadratic(M, v)]]))
            directions.append(dict(parent_block=b, vector=v))
            for j in range(i):
                for sign in (-1, 1):
                    v = [(j, F(sign)), (i, F(1))]
                    blocks.append(dict(name=name+f"_pair{i}_{j}_{sign}", entries=[[quadratic(M, v)]]))
                    directions.append(dict(parent_block=b, vector=v))
    return conic.model(primal["variables"], blocks, primal["equations"], dict(primal["metadata"],
        PSD_replacement="necessary rational quadratic directions; not an equivalence",
        parent_SDP_sha256=conic.digest_bytes(conic.canonical(primal)),
        cone_directions=directions))


def lift_certificate(primal_dual, outer, certificate):
    result = {n: F(0) for n in primal_dual["variables"]}
    for n in result:
        if n.startswith("potential_"):
            result[n] = certificate[n]
    for k, direction in enumerate(outer["metadata"]["cone_directions"]):
        multiplier = certificate[f"dual_b{k}_0_0"]
        if multiplier < 0:
            raise ValueError("A PSD direction has a negative certificate multiplier")
        b, vector = direction["parent_block"], direction["vector"]
        for i, a in vector:
            for j, c in vector:
                if i <= j:
                    result[f"dual_b{b}_{i}_{j}"] += multiplier*a*c
    checked = conic.check_point(primal_dual, result)
    return result, dict(equalities_checked=checked["equalities_checked"],
        PSD_blocks_checked=len(checked["blocks"]), candidate_sha256=checked["candidate_sha256"],
        method="nonnegative rational sums of outer products; original adjoint and PSD checked")


def explore(primal, dual, output, solve, recover, seconds=360, rounds=250,retain_outer=True,iteration_seconds=40):
    import numpy as np
    out = output
    outer = outer_LP(primal)
    seen = {conic.canonical(b["entries"][0][0]) for b in outer["blocks"]}
    started = time.monotonic()
    history, last_point = [], None
    report = dict(status="NO_VERIFIED_GLOBAL_CERTIFICATE_AT_THIS_LEVEL", global_bound_changed=False)
    for iteration in range(rounds):
        left = seconds-(time.monotonic()-started)
        if left < 3:
            break
        candidate, numerical = solve(outer, seconds=min(iteration_seconds, left))
        item = dict(iteration=iteration, outer_scalar_rows=len(outer["blocks"]), numerical=numerical)
        history.append(item)
        if candidate is None:
            if numerical["status"] == 2:
                certificate = conic.farkas(outer)
                left = seconds-(time.monotonic()-started)
                proposed, status = solve(certificate, seconds=max(1, min(90, left)))
                conic.write_json(out/"PSD_cut_certificate_numerical.json", status)
                if proposed is not None:
                    conic.write_json(out/"PSD_cut_certificate_candidate.json", proposed)
                recovered = recover(certificate, proposed)
                if recovered:
                    point, checked = recovered
                    conic.write_json(out/"PSD_cut_LP_certificate.exact.json", point)
                    full, full_check = lift_certificate(dual, outer, point)
                    conic.write_json(out/"PSD_cut_full_SDP_certificate.exact.json", full)
                    report.update(status="EXACT_GLOBAL_SDP_INFEASIBILITY_CERTIFICATE",
                        LP_certificate_check=checked, original_SDP_certificate_check=full_check,
                        accepted_dimension_bound_update=False)
            break
        last_point = candidate
        violations, bad_blocks = [], 0
        smallest = 0.0
        for b, block in enumerate(primal["blocks"]):
            if len(block["entries"]) == 1:
                continue
            M = np.asarray([[float(conic.value(entry, candidate)) for entry in row]
                for row in block["entries"]])
            eigenvalues, eigenvectors = np.linalg.eigh(M)
            scale = max(1.0, float(np.max(np.abs(M))))
            smallest = min(smallest, float(eigenvalues[0])/scale)
            if eigenvalues[0] < -1e-8*scale:
                bad_blocks += 1
                vector = [(i, F(str(float(v))).limit_denominator(10000))
                    for i, v in enumerate(eigenvectors[:, 0]) if abs(v) > 1e-10]
                row = quadratic(block["entries"], vector)
                if conic.value(row, candidate) < 0 and conic.canonical(row) not in seen:
                    violations.append((float(eigenvalues[0])/scale, b, vector, row))
        violations.sort(key=lambda a: a[0])
        item["PSD_blocks_with_new_separating_direction"] = len(violations)
        item["PSD_blocks_violated_at_numerical_tolerance"] = bad_blocks
        item["most_negative_scaled_eigenvalue"] = smallest
        item["worst_original_PSD_blocks"] = [dict(name=primal["blocks"][b]["name"],
            scaled_minimum_eigenvalue=eigenvalue) for eigenvalue, b, _, _ in violations[:5]]
        if iteration % 10 == 0 or not violations:
            print(f"Global PSD separation round {iteration}: {len(violations)} violated blocks", flush=True)
        if not violations:
            recovered = recover(primal, candidate) if not bad_blocks else None
            if recovered:
                point, checked = recovered
                conic.write_json(out/"PSD_feasible_point.exact.json", point)
                report.update(status="EXACT_GLOBAL_OUTER_SDP_FEASIBLE", original_point_check=checked,
                    no_strict_Farkas_certificate_at_this_level=True)
            else:
                report["status"] = ("PSD_SEPARATION_STALLED_WITH_VIOLATIONS" if bad_blocks else
                    "GLOBAL_SDP_NUMERICALLY_FEASIBLE_NOT_EXACTLY_VERIFIED")
            break
        for _, b, vector, row in violations[:80]:
            seen.add(conic.canonical(row))
            outer["blocks"].append(dict(name=f"PSD_separation_{iteration}_block{b}", entries=[[row]]))
            outer["metadata"]["cone_directions"].append(dict(parent_block=b, vector=vector))
    if last_point is not None:
        conic.write_json(out/"PSD_cut_last_primal_candidate.json", last_point)
    if retain_outer:
        conic.write_json(out/"PSD_cut_outer.exact.json", outer)
    else:
        original_rows = sum(len(b["entries"])**2 for b in primal["blocks"])
        conic.write_json(out/"additional_PSD_directions.exact.json",dict(
            original_SDP_sha256=conic.digest_bytes(conic.canonical(primal)),
            initial_directions="all e_i and e_i plus/minus e_j, in parent block order",
            initial_row_count=original_rows,
            additional_directions=outer["metadata"]["cone_directions"][original_rows:]))
    conic.write_json(out/"PSD_separation_history.json", history)
    report.update(PSD_separation_rounds=len(history), elapsed_seconds=time.monotonic()-started,
        all_cuts_are_rational_tests_of_original_PSD_blocks=True,
        last_saved_candidate_is_unverified_and_may_violate_final_cuts=True,
        new_diagram_types_added=False)
    return report
