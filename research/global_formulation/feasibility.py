"""Solve a frozen global SDP, with exact cone-face reductions and both orientations.

No geometric row is added here. The input is the retained exact SDP rather
than a rebuilt or weakened model. Numerical answers are checked in that input.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

from . import conic


DEFAULT_INPUT = "results/global_tail57_current/necessary_primal.exact.json"


def substitute(expression, substitutions):
    return conic.plus(conic.form(expression.get("", F(0))),
        *(conic.times(c, substitutions[v]) for v, c in expression.items() if v))


def positive_kernel(rows, progress):
    """Find, then verify over Q, a positive identity among cone test rows."""
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix, csc_matrix, eye, hstack
    unique = {}
    for ref, q in rows:
        if not q: continue
        scale = max(abs(c) for c in q.values())
        normalized = conic.times(1/scale, q)
        unique.setdefault(conic.canonical(normalized), (ref, normalized, scale))
    columns = list(unique.values())
    names = sorted({v for _, q, _ in columns for v in q})
    indices = {v: i for i, v in enumerate(names)}
    ii, jj, vv = [], [], []
    for j, (_, q, _) in enumerate(columns):
        for v, c in q.items(): ii.append(indices[v]); jj.append(j); vv.append(float(c))
    m = len(columns)
    if not m: return [], dict(status="empty")
    B = coo_matrix((vv, (ii, jj)), shape=(len(names), m)).tocsc()
    answer = linprog(np.r_[np.zeros(m), -np.ones(m)],
        A_eq=hstack((B, csc_matrix(B.shape))), b_eq=np.zeros(len(names)),
        A_ub=hstack((-eye(m, format="csc"), eye(m, format="csc"))), b_ub=np.zeros(m),
        bounds=[(0, None)]*m+[(0, 1)]*m, method="highs",
        options=dict(time_limit=30, primal_feasibility_tolerance=1e-9,
            dual_feasibility_tolerance=1e-9, threads=1))
    info = dict(status=int(answer.status), message=str(answer.message), test_rows=m)
    if answer.x is None: return [], info
    support = [j for j in range(m) if answer.x[j] > 1e-7]
    info["proposed_support"] = len(support)
    if not support: return [], info
    exact = None
    for denominator in (1000, 10**6, 10**9):
        proposed = {j: F(str(answer.x[j])).limit_denominator(denominator) for j in support}
        if all(c >= 0 for c in proposed.values()) and not conic.plus(*(
            conic.times(c, columns[j][1]) for j, c in proposed.items())):
            exact = proposed; info["rational_denominator_limit"] = denominator; break
    if exact is None and len(support) <= 3500:
        variable = {j: "lambda_"+str(j) for j in support}
        equations = {}
        for j in support:
            for v, c in columns[j][1].items(): equations.setdefault(v, {})[variable[j]] = c
        normalization = {variable[j]: F(1) for j in support}; normalization[""] = -F(1)
        auxiliary = conic.model(list(variable.values()), [dict(name=v, entries=[[conic.form(**{v: 1})]])
            for v in variable.values()], list(equations.values())+[normalization])
        try:
            reduced, substitutions = conic.eliminate(auxiliary)
            total = sum((F(str(answer.x[j])) for j in support), F(0))
            numerical = {variable[j]: F(str(answer.x[j]))/total for j in support}
            checked, _ = conic.reconstruct(auxiliary, reduced, substitutions,
                [numerical[v] for v in reduced["variables"]])
            exact = {j: checked[variable[j]] for j in support if checked[variable[j]] > 0}
            info["rational_reconstruction_with_exact_nullspace_projection"] = True
        except ValueError as error: info["exact_recovery_error"] = str(error)
    if exact is None:
        progress("scalar_face_search_inconclusive", info); return [], info
    if any(c < 0 for c in exact.values()) or conic.plus(*(
        conic.times(c, columns[j][1]) for j, c in exact.items())):
        raise ValueError("The proposed facial LP identity failed its exact check")
    result = [(columns[j][0], c/columns[j][2], columns[j][1]) for j, c in exact.items() if c]
    info["verified_nonnegative_weights"] = len(result)
    progress("verified_scalar_face", info)
    return result, info


def prepare(problem, progress, facial=True, linear_faces=True):
    """Equivalent exact reductions: equalities, PSD zero rows, scalar faces.

    A scalar face is accepted only after checking a nonnegative combination
    of existing scalar cone rows/PSD diagonals is identically zero.
    Every positive-weight summand must then vanish at any feasible point.
    A PSD zero diagonal forces its entire row to vanish.
    """
    current, substitutions = conic.eliminate(problem)
    maps = [dict(original_block=b, indices=list(range(len(block["entries"]))), scale=F(1))
        for b, block in enumerate(problem["blocks"])]
    trace = []
    for step in range(16 if facial else 1):
        derived, proofs = {}, []
        units, scalar_rows, zero_directions = {}, [], {}
        for b, block in enumerate(current["blocks"]):
            M = block["entries"]
            for i in range(len(M)):
                scalar_rows.append(((b, [(i, F(1))]), M[i][i]))
                if facial:
                    for j in range(i):
                        for sign in (-1, 1):
                            scalar_rows.append(((b, [(j, F(sign)), (i, F(1))]),
                                conic.plus(M[i][i], M[j][j], conic.times(2*sign, M[i][j]))))
        for source_index, (_, q) in enumerate(scalar_rows):
            if "" not in q and len(q) == 1:
                v, c = next(iter(q.items()))
                units.setdefault(v, {})[1 if c > 0 else -1] = (source_index, c)
        if facial:
            normalized = {}
            for ref, q in scalar_rows:
                b, direction = ref
                if not q:
                    if len(current["blocks"][b]["entries"]) > 1:
                        zero_directions.setdefault(b, []).append(direction)
                    continue
                scale = max(abs(c) for c in q.values())
                canonical = conic.times(1/scale, q)
                opposite = conic.canonical(conic.times(-1, canonical))
                if opposite in normalized:
                    old_ref, old_q, old_scale = normalized[opposite]
                    combination = [(old_ref, 1/old_scale), (ref, 1/scale)]
                    if conic.plus(conic.times(1/old_scale, old_q), conic.times(1/scale, q)):
                        raise ValueError("Invalid opposite quadratic-cone identity")
                    derived[conic.canonical(q)] = q
                    proofs.append(dict(kind="nonnegative_quadratic_sum_zero", combination=combination))
                else:
                    normalized[conic.canonical(canonical)] = (ref, q, scale)
            for v, signs in units.items():
                if len(signs) == 2:
                    b, c = signs[1]; d, a = signs[-1]
                    combo = [(b, F(1)), (d, -c/a)]
                    if conic.plus(*(conic.times(w, scalar_rows[i][1]) for i, w in combo)):
                        raise ValueError("Invalid scalar-face identity")
                    derived[conic.canonical(conic.form(**{v: 1}))] = conic.form(**{v: 1})
                    proofs.append(dict(kind="nonnegative_scalar_sum_zero",
                        combination=[(scalar_rows[i][0], w) for i, w in combo], forced_variables=[v]))
            for source_index, (_, q) in enumerate(scalar_rows):
                if not q or "" in q:
                    continue
                combo = [(source_index, F(1))]
                for v, c in q.items():
                    sign = -1 if c > 0 else 1
                    if sign not in units.get(v, {}):
                        break
                    source, coefficient = units[v][sign]
                    combo.append((source, -c/coefficient))
                else:
                    if any(w < 0 for _, w in combo) or conic.plus(*(
                        conic.times(w, scalar_rows[i][1]) for i, w in combo)):
                        raise ValueError("Invalid positive scalar-face identity")
                    for v in q:
                        row = conic.form(**{v: 1}); derived[conic.canonical(row)] = row
                    proofs.append(dict(kind="nonnegative_scalar_sum_zero",
                        combination=[(scalar_rows[i][0], w) for i, w in combo], forced_variables=list(q)))
        if facial and linear_faces and not derived:
            exposed, info = positive_kernel(scalar_rows, progress)
            for _, _, q in exposed: derived[conic.canonical(q)] = q
            if exposed:
                proofs.append(dict(kind="LP_verified_nonnegative_quadratic_sum_zero",
                    combination=[(ref, coefficient) for ref, coefficient, _ in exposed],
                    numerical_search=info, exact_identity_checked=True))
        blocks, new_maps, keys, dropped = [], [], set(), 0
        for b, (block, mapping) in enumerate(zip(current["blocks"], maps, strict=True)):
            M = block["entries"]
            keep = list(range(len(M)))
            if facial and len(M) > 1:
                zero = [i for i in keep if not M[i][i]]
                for i in zero:
                    rows = [q for q in M[i] if q]
                    for q in rows: derived[conic.canonical(q)] = q
                    if rows: proofs.append(dict(kind="PSD_zero_diagonal", block=b, index=i, forced_rows=rows))
                for direction in zero_directions.get(b, []):
                    rows = [conic.plus(*(conic.times(c, row[j]) for j, c in direction)) for row in M]
                    for q in rows:
                        if q: derived[conic.canonical(q)] = q
                    # All two-coordinate directions express the larger
                    # column through a smaller one; their union preserves
                    # a principal spanning set, even when relations repeat.
                    zero.append(max(j for j, _ in direction))
                    if any(rows): proofs.append(dict(kind="PSD_zero_quadratic_direction", block=b,
                        direction=direction, forced_rows=rows))
                keep = [i for i in keep if i not in zero]
            if not keep:
                dropped += 1; continue
            M = [[M[i][j] for j in keep] for i in keep]
            if not any(v for row in M for q in row for v in q if v):
                conic.psd_pivots([[q.get("", F(0)) for q in row] for row in M])
                dropped += 1; continue
            scale = max(abs(c) for row in M for q in row for c in q.values())
            M = [[conic.times(1/scale, q) for q in row] for row in M]
            key = conic.canonical(M)
            if key in keys:
                dropped += 1; continue
            keys.add(key)
            blocks.append(dict(name=block["name"], entries=M))
            new_maps.append(dict(original_block=mapping["original_block"],
                indices=[mapping["indices"][i] for i in keep], scale=mapping["scale"]/scale))
        item = dict(step=step, input_variables=len(current["variables"]), input_blocks=len(current["blocks"]),
            remaining_blocks=len(blocks), derived_equalities=len(derived), dropped_blocks=dropped, proofs=proofs)
        trace.append(item)
        if derived:
            augmented = conic.model(current["variables"], blocks, list(derived.values()))
            current, updates = conic.eliminate(augmented)
            substitutions = {v: substitute(q, updates) for v, q in substitutions.items()}
        else:
            current = conic.model(current["variables"], blocks)
        maps = new_maps
        progress("exact_presolve", {k: v for k, v in item.items() if k != "proofs"})
        if not derived:
            break
    return current, substitutions, maps, trace


def reduction_checks():
    toy = conic.model(["x", "y"], [
        dict(name="zero_face", entries=[[conic.form(x=1), conic.form(y=1)],
                                        [conic.form(y=1), conic.form(1)]]),
        dict(name="opposite", entries=[[conic.form(x=-1)]])])
    compiled, substitutions, _, _ = prepare(toy, lambda *_: None, linear_faces=False)
    if compiled["variables"]: raise ValueError("The exact zero-diagonal face was not eliminated")
    conic.check_point(toy, {v: conic.value(q, {}) for v, q in substitutions.items()})
    toy = conic.model(["x", "y"], [dict(name="common_kernel", entries=[
        [conic.form(1), conic.form(1), conic.form(x=1)],
        [conic.form(1), conic.form(1), conic.form(y=1)],
        [conic.form(x=1), conic.form(y=1), conic.form(2)]])])
    compiled, substitutions, _, _ = prepare(toy, lambda *_: None, linear_faces=False)
    if len(compiled["variables"]) != 1 or max(len(b["entries"]) for b in compiled["blocks"]) != 2:
        raise ValueError("The exact noncoordinate PSD kernel was not reduced")
    conic.check_point(toy, {v: conic.value(q, {compiled["variables"][0]: F(0)}) for v, q in substitutions.items()})
    return dict(exact_zero_diagonal_face=True, exact_noncoordinate_kernel_face=True)


def native_cones(path, problem, field="yMat"):
    text = Path(path).read_text()
    if len(text) > 32*1024**2:
        raise ValueError("Oversized native SDP answer")
    match = re.search(r"\b"+field+r"\s*=\s*\{", text)
    if not match: raise ValueError("No native "+field)
    start, depth = match.end()-1, 0
    for end in range(start, len(text)):
        depth += int(text[end] == "{")-int(text[end] == "}")
        if not depth: break
    else: raise ValueError("Unterminated native cone array")
    values = [F(s) for s in conic.NUMBER.findall(text[start:end+1])]
    matrices, position = {}, 0
    for b, block in enumerate(problem["blocks"]):
        r = len(block["entries"])
        if r == 1: continue
        if position+r*r > len(values): raise ValueError("Short native PSD array")
        M = [values[position+i*r:position+(i+1)*r] for i in range(r)]
        matrices[b] = [[(M[i][j]+M[j][i])/2 for j in range(r)] for i in range(r)]
        position += r*r
    for b, block in enumerate(problem["blocks"]):
        if len(block["entries"]) == 1:
            if position >= len(values): raise ValueError("Short native scalar array")
            matrices[b] = [[values[position]]]; position += 1
    if position != len(values): raise ValueError("Unexpected native cone array length")
    return [matrices[b] for b in range(len(problem["blocks"]))]


def numerical_residuals(problem, point):
    import numpy as np
    equality = max((abs(float(conic.value(q, point))) for q in problem["equations"]), default=0.0)
    worst, worst_scaled, name = 0.0, 0.0, None
    for block in problem["blocks"]:
        M = np.array([[float(conic.value(q, point)) for q in row] for row in block["entries"]])
        eigenvalue = float(np.linalg.eigvalsh(M)[0]) if len(M) > 1 else float(M[0, 0])
        scaled = eigenvalue/max(1.0, float(np.max(np.abs(M))))
        if eigenvalue < worst: worst = eigenvalue
        if scaled < worst_scaled: worst_scaled, name = scaled, block["name"]
    return dict(maximum_absolute_original_equality_residual=equality,
        minimum_original_cone_eigenvalue=worst, minimum_scaled_original_cone_eigenvalue=worst_scaled,
        worst_original_block=name, floating_diagnostics_are_not_exact_certificates=True)


def solve_sdpa(problem, out, seconds, progress):
    solver = shutil.which("sdpa_gmp") or shutil.which("sdpa")
    if not solver: raise RuntimeError("SDPA is unavailable")
    encoding = conic.write_sdpa(out/"input.dat-s", problem)
    parameters = "35\n1e-12\n10\n2\n-1e12\n1e12\n0.1\n0.3\n0.9\n1e-12\n128\n"+"%+90.80Fe\n"*4
    (out/"param.sdpa").write_text(parameters)
    command = [solver, "-ds", str(out/"input.dat-s"), "-o", str(out/"native.result"), "-p", str(out/"param.sdpa")]
    progress("native_solver_start", dict(command=command, encoding=encoding))
    code, timeout = None, False
    with (out/"solver.stdout").open("x") as stream:
        try:
            code = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                timeout=seconds, check=False).returncode
        except subprocess.TimeoutExpired: timeout = True
    numerical = dict(solver="SDPA-GMP" if "gmp" in solver else "SDPA", command=command,
        returncode=code, timed_out=timeout, encoding=encoding)
    x, matrices = None, None
    try:
        raw, phase = conic.native_vector(out/"native.result")
        if len(raw) != len(problem["variables"]): raise ValueError("Wrong native vector dimension")
        x = dict(zip(problem["variables"], raw, strict=True)); numerical["phase"] = phase
        matrices = native_cones(out/"native.result", problem)
    except (ValueError, FileNotFoundError) as error: numerical["native_output_error"] = str(error)
    return x, matrices, numerical


def run(root, output, input_path, solver, orientation, seconds=600, facial=True, resume_from=None):
    if not 1 <= seconds <= 700: raise ValueError("Use a bounded numerical stage of at most 700 seconds")
    out, source = (root/output).resolve(), (root/input_path).resolve()
    if out.exists() or not out.is_relative_to(root/"results"): raise ValueError("Use a fresh results directory")
    if not source.is_relative_to(root/"results"): raise ValueError("Use a retained exact results input")
    out.mkdir()
    started = time.monotonic()
    paths = [source, Path(__file__), root/"research/global_formulation/conic.py",
        root/"scripts/run_global_coxeter_formulation.py"]
    resume_path = None
    if resume_from:
        resume_path = (root/resume_from/"original_numerical_candidate.json").resolve()
        if not resume_path.is_relative_to(root/"results"): raise ValueError("Invalid numerical checkpoint")
        paths.append(resume_path)
    bindings = {str(p.relative_to(root)): dict(bytes=p.stat().st_size,
        sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    conic.write_json(out/"source_bindings.json", bindings)
    def progress(stage, data):
        item = dict(data, elapsed_seconds=time.monotonic()-started, stage=stage)
        with (out/"progress.jsonl").open("a") as stream: stream.write(json.dumps(conic.plain(item))+"\n")
        print(json.dumps(conic.plain(item)), flush=True)
    original = conic.read_model(source)
    checks = reduction_checks()
    target = original if orientation == "necessary" else conic.farkas(original)
    progress("frozen_model_loaded", dict(orientation=orientation, variables=len(target["variables"]),
        equalities=len(target["equations"]), cones=len(target["blocks"])))
    compiled, substitutions, maps, trace = prepare(target, progress, facial=facial)
    conic.write_json(out/"compiled.exact.json", compiled)
    conic.write_json(out/"substitutions.exact.json", substitutions)
    conic.write_json(out/"cone_maps.exact.json", maps)
    conic.write_json(out/"exact_face_reduction.json", trace)
    report = dict(status="EXACT_EQUIVALENT_FROZEN_SDP_PREPARED", orientation=orientation,
        frozen_input=str(source.relative_to(root)), frozen_input_sha256=bindings[str(source.relative_to(root))]["sha256"],
        geometric_constraints_added=False, original_variables=len(target["variables"]),
        compiled_variables=len(compiled["variables"]), compiled_cones=len(compiled["blocks"]),
        accepted_bounds_changed=False, source_bindings=bindings, exact_reduction_checks=checks)
    candidate = matrices = None
    if solver == "recover":
        if resume_path is None: raise ValueError("Choose a numerical checkpoint to recover")
        raw = json.loads(resume_path.read_text())
        candidate = {v: F(raw[v]) for v in compiled["variables"]}
        report.update(status="NUMERICAL_CHECKPOINT_RECOVERY", resume_from=str(resume_path.relative_to(root)))
    elif solver == "sdpa":
        candidate, matrices, numerical = solve_sdpa(compiled, out, seconds, progress)
        report.update(status="NO_EXACT_FEASIBILITY_DECISION", numerical=numerical)
    elif solver != "prepare": raise ValueError("Unsupported solver")
    if candidate is not None:
        point = {v: conic.value(q, candidate) for v, q in substitutions.items()}
        conic.write_json(out/"original_numerical_candidate.json", point)
        report["original_numerical_residuals"] = numerical_residuals(target, point)
        try:
            recovered, checked = conic.reconstruct(target, compiled, substitutions,
                [candidate[v] for v in compiled["variables"]])
            conic.write_json(out/"original_feasible_point.exact.json", recovered)
            report.update(status="EXACT_ORIGINAL_SDP_FEASIBLE" if orientation == "necessary" else
                "EXACT_ORIGINAL_SDP_INFEASIBILITY_CERTIFICATE", exact_check=checked)
        except ValueError as error: report["exact_primal_recovery_failure"] = str(error)
    if matrices is not None:
        conic.write_json(out/"numerical_dual_PSD_blocks.json", matrices)
        constant = sum((sum((q.get("", F(0))*Z[i][j]
            for i, row in enumerate(b["entries"]) for j, q in enumerate(row)), F(0))
            for b, Z in zip(compiled["blocks"], matrices, strict=True)), F(0))
        report["native_dual_constant_pairing"] = constant
    after = {str(p.relative_to(root)): dict(bytes=p.stat().st_size,
        sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    if after != bindings: raise ValueError("A fixed solver source or model changed during the run")
    report["elapsed_seconds"] = time.monotonic()-started
    conic.write_json(out/"result.json", report)
    progress("completed", {k: v for k, v in report.items() if k not in ("source_bindings", "exact_check")})
    return report


def bareiss_psd(matrix):
    """Independent integer PSD test: symmetric-pivoted Bareiss elimination.

    A positive common denominator first gives an integer symmetric matrix.
    Fraction-free Schur updates must divide exactly by the preceding pivot.
    Positive pivots preserve PSD equivalence; a zero remaining diagonal
    requires the whole remaining block to be zero.
    """
    n = len(matrix)
    denominator = math.lcm(*(v.denominator for row in matrix for v in row))
    A = [[int(v*denominator) for v in row] for row in matrix]
    if any(A[i][j] != A[j][i] for i in range(n) for j in range(n)):
        raise ValueError("Nonsymmetric exact matrix")
    previous = 1
    for k in range(n):
        if any(A[i][i] < 0 for i in range(k, n)):
            raise ValueError("Negative integer Schur diagonal")
        pivot_index = next((i for i in range(k, n) if A[i][i] > 0), None)
        if pivot_index is None:
            if any(A[i][j] for i in range(k, n) for j in range(k, n)):
                raise ValueError("Nonzero row in a zero-diagonal Schur block")
            return k
        if pivot_index != k:
            A[k], A[pivot_index] = A[pivot_index], A[k]
            for row in A: row[k], row[pivot_index] = row[pivot_index], row[k]
        pivot = A[k][k]
        for i in range(k+1, n):
            for j in range(i, n):
                numerator = pivot*A[i][j]-A[i][k]*A[k][j]
                quotient, remainder = divmod(numerator, previous)
                if remainder: raise ValueError("Nonexact fraction-free division")
                A[i][j] = A[j][i] = quotient
        for i in range(k+1, n): A[i][k] = A[k][i] = 0
        previous = pivot
    return n


def verify(root, output, input_path, point_path):
    """Check a frozen primal witness and a certificate ruling out its Farkas dual."""
    out = (root/output).resolve()
    if out.exists() or not out.is_relative_to(root/"results"): raise ValueError("Use a fresh results directory")
    source, witness = (root/input_path).resolve(), (root/point_path).resolve()
    dual_source = source.with_name("certificate_dual.exact.json")
    if not source.is_relative_to(root/"results") or not witness.is_relative_to(root/"results"):
        raise ValueError("Use retained results inputs")
    out.mkdir()
    paths = [source, witness, dual_source, Path(__file__), root/"research/global_formulation/conic.py",
        root/"scripts/run_global_coxeter_formulation.py"]
    bindings = {str(p.relative_to(root)): dict(bytes=p.stat().st_size,
        sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    conic.write_json(out/"source_bindings.json", bindings)
    for matrix, rank in (([[F(1), F(1)], [F(1), F(1)]], 1),
            ([[F(0), F(0), F(0)], [F(0), F(1), F(1)], [F(0), F(1), F(1)]], 1),
            ([[F(2), F(1, 3)], [F(1, 3), F(2)]], 2)):
        if bareiss_psd(matrix) != rank: raise ValueError("Independent PSD check failed its rank fixture")
    try:
        bareiss_psd([[F(0), F(1)], [F(1), F(0)]])
    except ValueError: pass
    else: raise ValueError("Independent PSD check accepted an indefinite matrix")
    problem = conic.read_model(source)
    point = {v: F(c) for v, c in json.loads(witness.read_text()).items()}
    if set(point) != set(problem["variables"]): raise ValueError("Wrong complete original variable set")
    def evaluate(row, values):
        return sum((c*(F(1) if v == "" else values[v]) for v, c in row.items()), F(0))
    if any(evaluate(q, point) for q in problem["equations"]):
        raise ValueError("An original equality failed independent verification")
    slacks, ranks = [], []
    for block in problem["blocks"]:
        M = [[evaluate(q, point) for q in row] for row in block["entries"]]
        rank = bareiss_psd(M)
        slacks.append(M); ranks.append(dict(name=block["name"], order=len(M), rank=rank))
    dual = conic.read_model(dual_source)
    if conic.canonical(dual) != conic.canonical(conic.farkas(problem)):
        raise ValueError("The retained certificate dual differs from the original primal's exact alternative")
    alternative = conic.farkas(dual)
    certificate = {v: F(0) for v in alternative["variables"]}
    for i, v in enumerate(problem["variables"]): certificate[f"potential_{i}"] = -point[v]
    certificate[f"potential_{len(problem['variables'])}"] = -F(1)
    for b, M in enumerate(slacks):
        for i in range(len(M)):
            for j in range(i, len(M)):
                certificate[f"dual_b{b}_{i}_{j}"] = M[i][j]
    if any(evaluate(q, certificate) for q in alternative["equations"]):
        raise ValueError("The original dual's infeasibility certificate failed its adjoint equations")
    for b, block in enumerate(alternative["blocks"]):
        for i, row in enumerate(block["entries"]):
            for j, q in enumerate(row):
                if evaluate(q, certificate) != slacks[b][i][j]:
                    raise ValueError("Dual-certificate PSD cones do not match the verified primal slacks")
    conic.write_json(out/"normalized_dual_infeasibility_certificate.exact.json", certificate)
    conic.write_json(out/"independent_PSD_ranks.json", ranks)
    report = dict(status="EXACT_FROZEN_PRIMAL_FEASIBLE_AND_NORMALIZED_CERTIFICATE_DUAL_INFEASIBLE",
        original_variables_checked=len(point), original_equalities_checked=len(problem["equations"]),
        original_PSD_blocks_checked=sum(r["order"] > 1 for r in ranks),
        original_scalar_inequalities_checked=sum(r["order"] == 1 for r in ranks),
        independent_PSD_method="integer symmetric-pivoted fraction-free Bareiss; no eigenvalue tolerances",
        original_point_canonical_sha256=conic.digest_bytes(conic.canonical(point)),
        normalized_dual_adjoint_equalities_checked=len(alternative["equations"]),
        dual_certificate_cones_equal_verified_primal_slacks=True,
        normalized_dual_certificate_potentials="minus primal coordinates, followed by -1",
        no_geometric_constraint_changed=True, geometric_realization_asserted=False, accepted_bounds_changed=False,
        source_bindings=bindings)
    after = {str(p.relative_to(root)): dict(bytes=p.stat().st_size,
        sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    if after != bindings: raise ValueError("A verification input changed during the check")
    conic.write_json(out/"result.json", report)
    print(json.dumps({k: v for k, v in report.items() if k != "source_bindings"}), flush=True)
    return report
