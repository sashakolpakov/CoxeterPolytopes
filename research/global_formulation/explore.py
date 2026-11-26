"""Export the global necessary problem and explore its complete conic dual."""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from . import conic, data, model, separation, strengthening, white_profiles, supporting_counts


def sampled_psd_LP(problem):
    """A declared outer LP: test e_i and e_i+-e_j in every PSD block."""
    return separation.outer_LP(problem)


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
        if len(supported) > 2000:
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
    outer = separation.outer_LP(P)
    proposed = {n: F(0) for n in conic.farkas(outer)["variables"]}
    proposed.update(potential_0=F(1), dual_b2_0_0=F(1, 2))
    conic.check_point(conic.farkas(outer), proposed)
    separation.lift_certificate(D, outer, proposed)
    # Check SDPA's constant sign and diagonal-cone packing independently
    # of the exact in-memory model used by the certificate checker.
    P = conic.model(["x"], [
        dict(name="scalar1", entries=[[conic.form(1, x=1)]]),
        dict(name="matrix", entries=[[conic.form(1), conic.form(x=1)],
                                    [conic.form(x=1), conic.form(1)]]),
        dict(name="scalar2", entries=[[conic.form(2, x=-1)]])])
    with tempfile.TemporaryDirectory(prefix="coxeter_sdpa_encoding_") as folder:
        path = Path(folder)/"check.dat-s"
        encoding = conic.write_sdpa(path, P)
        lines = path.read_text().splitlines()
        if lines[3] != "2 -2" or encoding["packed_scalar_constraints"] != 2:
            raise ValueError("The SDPA cone packing changed its intended product cone")
        encoded = {}
        for line in lines[5:]:
            matrix, block, i, j, value = line.split()
            key = (int(block), int(i), int(j))
            factor = F(2) if int(matrix) else F(-1)
            encoded[key] = encoded.get(key, F(0))+factor*F(value)
        if encoded != {(1, 1, 1): F(1), (1, 1, 2): F(2), (1, 2, 2): F(1),
                       (2, 1, 1): F(3), (2, 2, 2): F(0)}:
            raise ValueError("The SDPA encoding does not evaluate to its original affine matrices")
    source = data.literal_module()
    h1000 = sum((F(1, k*k) for k in range(1, 1001)), F(0))
    if not (data.HARMONIC2_LOWER <= h1000+F(1, 1001) and h1000+F(1, 1000) <= data.HARMONIC2_UPPER):
        raise ValueError("The exact integral-test harmonic bounds failed")
    for kind, rank in (("A", 2), ("A", 7), ("B", 2), ("B", 7), ("D", 4), ("D", 12)):
        inverse = data.inverse_angle_square_sum(data.rational_root_Gram(kind, rank))
        if data.component_angle_bounds(kind, rank) != (inverse, inverse):
            raise ValueError("A closed inverse-Gram angle formula failed its rational matrix check")
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
        SDPA_constant_sign_and_scalar_packing=True,
        exact_rational_PSD_cut_certificate_lift=True,
        harmonic_envelope_integral_proof_checked=True, closed_angle_square_inverse_checks=6,
        finite_tail_regressions=finite_checks, affine_tail_regressions=affine_checks,
        finite_regressions_do_not_replace_the_uniform_proof=True)


def run(root, output, threshold, parity, stage, upper_dimension=None, level=1):
    if stage!='export':
        raise RuntimeError('This exploratory solver route is retired; use fixed-point preparation and the gated native driver')
    out = (root/output).resolve()
    if out.exists() or not out.is_relative_to(root/"results"):
        raise ValueError("Use a fresh results directory")
    paths = [*sorted((root/"research/global_formulation").glob("*.py")),
        root/"notes/lp/global_formulation.md", root/"notes/lp/global_sources.md",
        root/"arxiv/linear_programming_dimension_bounds/main.tex",
        root/"manuscripts/finite_volume57_consolidated.tex",
        root/"research/global_formulation/sources.json", root/data.LITERAL]
    audit_57 = level > 1 and threshold == 57 and parity == "odd" and upper_dimension == 57
    baseline_input = root/("results/global_shared_cusp_current/necessary_primal.exact.json" if level >= 6 else
        "results/global_level2_current/necessary_primal.exact.json" if level >= 3 else
        "results/global_tail57_current/necessary_primal.exact.json")
    baseline_point = root/("results/global_shared_cusp_scalar_circuits01/original_feasible_point.exact.json" if level >= 6 else
        "results/global_level2_face_recovery/original_feasible_point.exact.json" if level >= 3 else
        "results/global_current_exact_recovery/original_feasible_point.exact.json")
    if audit_57:
        paths += [baseline_input,baseline_point]
    bindings = {str(p.relative_to(root)): dict(bytes=p.stat().st_size,
        sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    out.mkdir()
    conic.write_json(out/"source_bindings.json", bindings)
    checks = self_checks()
    if level > 1:
        checks.update(strengthening.symbolic_gram_checks())
    if level >= 3:
        checks.update(white_profiles.self_checks())
        checks.update(supporting_counts.self_checks())
    if level >= 6:
        from . import edge_hosts, edge_compatibility
        checks.update(edge_hosts.self_checks())
        checks.update(edge_compatibility.self_checks())
    tables = data.intrinsic_tables()
    primal = model.build(threshold, parity, tables, upper_dimension, level=level)
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
    if level > 1:
        if audit_57:
            base = conic.read_model(baseline_input)
            point = {v: F(c) for v,c in json.loads(baseline_point.read_text()).items()}
            audit = strengthening.compare_and_test(base, primal, point)
            conic.write_json(out/"strict_strengthening_audit.json", audit)
            report["strict_strengthening"] = {k: v for k, v in audit.items()
                if k not in ("equality_violations", "cone_violations")}
        report["uniform_strengthening_level"] = level
    after = {str(p.relative_to(root)): dict(bytes=p.stat().st_size,
        sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    if after != bindings:
        raise ValueError("A global model source changed during the batch")
    report["source_bindings"] = bindings
    conic.write_json(out/"result.json", report)
    print(json.dumps(conic.plain({k: v for k, v in report.items() if k != "source_bindings"})), flush=True)
    return report
