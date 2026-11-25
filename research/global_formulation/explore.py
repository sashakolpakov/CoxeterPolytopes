"""Export the global necessary problem and explore its complete conic dual."""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import shutil
import subprocess

from . import conic, data, model


def sampled_psd_LP(problem):
    """A declared outer LP: test e_i and e_i+-e_j in every PSD block."""
    blocks = []
    for b in problem["blocks"]:
        M, name = b["entries"], b["name"]
        if len(M) == 1:
            blocks.append(b); continue
        for i in range(len(M)):
            blocks.append(dict(name=name+f"_diag{i}", entries=[[M[i][i]]]))
            for j in range(i):
                for sign in (-1, 1):
                    row = conic.plus(M[i][i], M[j][j], conic.times(2*sign, M[i][j]))
                    blocks.append(dict(name=name+f"_pair{i}_{j}_{sign}", entries=[[row]]))
    return conic.model(problem["variables"], blocks, problem["equations"], dict(problem["metadata"],
        PSD_replacement="necessary rational test directions, not an equivalence",
        parent_SDP_sha256=conic.digest_bytes(conic.canonical(problem))))


def linprog(problem, seconds=300):
    import numpy as np
    from scipy.optimize import linprog as solve
    from scipy.sparse import coo_matrix
    columns = {n: i for i, n in enumerate(problem["variables"])}
    def array(rows):
        ii, jj, vv, rhs = [], [], [], []
        for i, row in enumerate(rows):
            scale = max((abs(v) for v in row.values()), default=F(1)) or F(1)
            for n, value in row.items():
                if n:
                    ii.append(i); jj.append(columns[n]); vv.append(float(value/scale))
            rhs.append(float(-row.get("", F(0))/scale))
        return coo_matrix((vv, (ii, jj)), shape=(len(rows), len(columns))).tocsr(), np.asarray(rhs)
    if any(len(b["entries"]) != 1 for b in problem["blocks"]):
        raise ValueError("Only the explicitly exported outer LP is sent to HiGHS")
    equality, rhs = array(problem["equations"])
    inequality, bound = array([conic.times(-1, b["entries"][0][0]) for b in problem["blocks"]])
    result = solve(np.zeros(len(columns)), A_ub=inequality, b_ub=bound,
        A_eq=equality, b_eq=rhs, bounds=(None, None), method="highs",
        options=dict(time_limit=seconds, primal_feasibility_tolerance=1e-9,
                     dual_feasibility_tolerance=1e-9, presolve=True, threads=1))
    answer = dict(status=int(result.status), success=bool(result.success), message=str(result.message))
    point = None if result.x is None else {n: F(str(result.x[i])) for n, i in columns.items()}
    return point, answer


def exact_reconstruction(problem, candidate):
    if candidate is None:
        return None
    for denominator in (1000, 1000000, 1000000000):
        point = {n: v.limit_denominator(denominator) for n, v in candidate.items()}
        try:
            check = conic.check_point(problem, point)
            return point, dict(method="rational_reconstruction", denominator=denominator,
                equalities_checked=check["equalities_checked"], PSD_blocks_checked=len(check["blocks"]),
                candidate_sha256=check["candidate_sha256"])
        except ValueError:
            pass
    # Sparse exact projection is a candidate recovery, never a replacement
    # for checking every original equality and cone constraint.
    for tolerance in (F(1, 10**8), F(1, 10**11)):
        supported = {n for n, v in candidate.items() if abs(v) > tolerance}
        if len(supported) > 400:
            continue
        reduced = conic.model([n for n in problem["variables"] if n in supported],
            [dict(b, entries=[[{n: v for n, v in entry.items() if not n or n in supported} for entry in row]
                for row in b["entries"]]) for b in problem["blocks"]],
            [{n: v for n, v in row.items() if not n or n in supported} for row in problem["equations"]])
        try:
            compiled, substitutions = conic.eliminate(reduced)
            raw = [candidate[n] for n in compiled["variables"]]
            point, _ = conic.reconstruct(reduced, compiled, substitutions, raw)
            full = {n: point.get(n, F(0)) for n in problem["variables"]}
            check = conic.check_point(problem, full)
            return full, dict(method="sparse_support_proposal_and_exact_projection", supported_variables=len(supported),
                equalities_checked=check["equalities_checked"], PSD_blocks_checked=len(check["blocks"]),
                candidate_sha256=check["candidate_sha256"])
        except ValueError:
            pass
    return None


def self_checks():
    # These check the conic alternative, including free equality potentials
    # and off-diagonal trace multiplicity. They are not polytope cases.
    P = conic.model(["x"], [dict(name="positive", entries=[[conic.form(x=1)]])], [conic.form(1, x=1)])
    D = conic.farkas(P)
    conic.check_point(D, dict(potential_0=F(-1), dual_b0_0_0=F(1)))
    P = conic.model(["x"], [dict(name="offdiag", entries=[[{}, conic.form(x=1)], [conic.form(x=1), {}]])],
        [conic.form(-1, x=1)])
    D = conic.farkas(P)
    conic.check_point(D, dict(potential_0=F(1), dual_b0_0_0=F(1, 2),
        dual_b0_0_1=-F(1, 2), dual_b0_1_1=F(1, 2)))
    source = data.literal_module()
    finite_checks, affine_checks = 0, 0
    for kind, (minimum, intercept) in data.FINITE_TAIL.items():
        for rank in (minimum, minimum+1, 57, 80):
            observed = sum((source.finite_same_weight((kind, rank), a, b)
                for a, b in combinations(range(rank), 2)), F(0))
            if observed != 5*rank+intercept:
                raise ValueError("The symbolic finite-tail statistic disagrees with the retained intrinsic law")
            finite_checks += 1
    for kind in ("A", "B", "C", "D"):
        for rank in (data.AFFINE_MINIMUM, data.AFFINE_MINIMUM+1, 57, 80):
            order = rank+1
            observed = sum((source.affine_weight(kind, order, a, b)
                for a, b in combinations(range(order), 2)), F(0))
            if observed != 5*order+{"A": 0, "B": 46, "C": 41, "D": 51}[kind]:
                raise ValueError("The symbolic affine-tail statistic disagrees with the retained intrinsic law")
            if kind == "A":
                deletion = data.prokhorov_statistics(data.graph("A", rank))[0]
            else:
                W, J = data.prokhorov_statistics(data.graph(kind, rank, affine=True))
                deletion = W-J/order
            beta, reciprocal = {"A": (-176, 0), "B": (-349, 1687), "C": (-352, 1696), "D": (-346, 1678)}[kind]
            if deletion != F(29*rank+beta, 3)+F(reciprocal, 3*order):
                raise ValueError("The symbolic affine deletion statistic disagrees with its literal graph")
            affine_checks += 1
    return dict(exact_affine_equality_alternative=True, exact_offdiagonal_trace_pairing=True,
        finite_tail_regressions=finite_checks, affine_tail_regressions=affine_checks,
        finite_regressions_do_not_replace_the_uniform_proof=True)


def run(root, output, threshold, parity, stage, upper_dimension=None):
    out = (root/output).resolve()
    if out.exists() or not out.is_relative_to(root/"results"):
        raise ValueError("Use a fresh results directory")
    paths = [*sorted((root/"research/global_formulation").glob("*.py")),
        root/"notes/lp/global_formulation.md", root/"notes/lp/global_sources.md",
        root/"arxiv/linear_programming_dimension_bounds/main.tex",
        root/"research/global_formulation/sources.json", root/data.LITERAL]
    bindings = {str(p.relative_to(root)): dict(bytes=p.stat().st_size,
        sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    out.mkdir()
    conic.write_json(out/"source_bindings.json", bindings)
    checks = self_checks()
    tables = data.intrinsic_tables()
    primal = model.build(threshold, parity, tables, upper_dimension)
    conic.write_json(out/"necessary_primal.exact.json", primal)
    dual = conic.farkas(primal)
    conic.write_json(out/"certificate_dual.exact.json", dual)
    report = dict(status="GLOBAL_NECESSARY_DATA_EXPORTED", threshold=threshold, parity=parity,
        upper_dimension=upper_dimension,
        primal_variables=len(primal["variables"]), primal_equalities=len(primal["equations"]),
        primal_PSD_blocks=sum(len(b["entries"]) > 1 for b in primal["blocks"]),
        primal_scalar_blocks=sum(len(b["entries"]) == 1 for b in primal["blocks"]),
        largest_PSD_block=max(len(b["entries"]) for b in primal["blocks"]),
        dual_variables=len(dual["variables"]), conic_alternative_checks=checks,
        diagram_family_lists_used=False, global_bound_changed=False, separate_author_review=False)
    if stage == "explore-lp":
        outer = sampled_psd_LP(primal)
        certificate = conic.farkas(outer)
        conic.write_json(out/"LP_outer.exact.json", outer)
        conic.write_json(out/"LP_certificate.exact.json", certificate)
        proposal, numerical = linprog(certificate)
        conic.write_json(out/"LP_certificate_numerical.json", numerical)
        recovered = exact_reconstruction(certificate, proposal)
        report["LP_certificate_numerical"] = numerical
        if recovered:
            point, checked = recovered
            conic.write_json(out/"LP_infeasibility_certificate.exact.json", point)
            report.update(status="EXACT_GLOBAL_OUTER_LP_INFEASIBILITY_CERTIFICATE", certificate_check=checked,
                accepted_dimension_bound_update=False)
        else:
            candidate, diagnostic = linprog(outer)
            conic.write_json(out/"LP_primal_numerical.json", diagnostic)
            recovered = exact_reconstruction(outer, candidate)
            report.update(status="NO_VERIFIED_GLOBAL_CERTIFICATE_AT_THIS_LEVEL", LP_primal_numerical=diagnostic)
            if recovered:
                point, checked = recovered
                conic.write_json(out/"LP_outer_feasible_point.exact.json", point)
                report.update(status="EXACT_GLOBAL_OUTER_LP_FEASIBLE", outer_point_check=checked)
                try:
                    full = conic.check_point(primal, point)
                    report.update(status="EXACT_GLOBAL_OUTER_SDP_FEASIBLE",
                        full_PSD_blocks_checked=len(full["blocks"]), no_strict_Farkas_certificate_at_this_level=True)
                except ValueError as error:
                    report["full_SDP_point_not_verified"] = str(error)
    elif stage == "explore-sdpa":
        compiled, substitutions = conic.eliminate(dual)
        conic.write_json(out/"dual_compiled.exact.json", compiled)
        conic.write_json(out/"dual_substitutions.exact.json", substitutions)
        conic.write_sdpa(out/"dual.dat-s", compiled)
        solver = shutil.which("sdpa_gmp") or shutil.which("sdpa")
        if not solver:
            raise RuntimeError("SDPA is not installed")
        command = [solver, "-ds", str(out/"dual.dat-s"), "-o", str(out/"dual.result")]
        with (out/"solver.stdout").open("x") as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=700, check=False)
        report.update(SDPA_command=command, SDPA_returncode=result.returncode,
            status="GLOBAL_SDPA_RESULT_REQUIRES_EXACT_RECONSTRUCTION")
        try:
            raw, phase = conic.native_vector(out/"dual.result")
            point, checked = conic.reconstruct(dual, compiled, substitutions, raw)
            conic.write_json(out/"SDPA_certificate.exact.json", point)
            report.update(status="EXACT_GLOBAL_SDP_INFEASIBILITY_CERTIFICATE", SDPA_phase=phase,
                certificate_check=checked, accepted_dimension_bound_update=False)
        except ValueError as error:
            report["exact_certificate_not_obtained"] = str(error)
    elif stage != "export":
        raise ValueError("Unsupported research stage")
    after = {str(p.relative_to(root)): dict(bytes=p.stat().st_size,
        sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    if after != bindings:
        raise ValueError("A global model source changed during the batch")
    report["source_bindings"] = bindings
    conic.write_json(out/"result.json", report)
    print(json.dumps(conic.plain({k: v for k, v in report.items() if k != "source_bindings"})), flush=True)
    return report
