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


def original_binomial_faces(problem, progress):
    """Expose exact positive identities before equality elimination changes sparsity.

    A row can cancel a disjoint union of existing one/two-term rows.
    Only positive weights and identically zero rational sums are accepted.
    In particular this catches inactive-count guards and their bound rows
    without a floating LP having to resolve small eliminated coefficients.
    """
    rows = {b: block["entries"][0][0] for b, block in enumerate(problem["blocks"])
        if len(block["entries"]) == 1 and block["entries"][0][0]}
    index = {}
    for b, q in rows.items():
        if len(q) <= 2:
            for v in q: index.setdefault(v, set()).add(b)
    derived, proofs = {}, []
    for b, q in rows.items():
        if not 3 <= len(q) <= 64: continue
        remainder, combination = dict(q), [(b, F(1))]
        while remainder:
            candidates = set().union(*(index.get(v, set()) for v in remainder))
            for j in sorted(candidates, key=lambda j: (-len(rows[j]), j)):
                r = rows[j]
                if not r.keys() <= remainder.keys(): continue
                ratios = {-remainder[v]/c for v, c in r.items()}
                if len(ratios) != 1: continue
                weight = next(iter(ratios))
                if weight <= 0: continue
                remainder = conic.plus(remainder, conic.times(weight, r))
                combination.append((j, weight))
                break
            else: break
        if remainder: continue
        if any(w <= 0 for _, w in combination) or conic.plus(*(
                conic.times(w, rows[j]) for j, w in combination)):
            raise ValueError("Invalid original-coordinate positive binomial identity")
        for j, _ in combination: derived[conic.canonical(rows[j])] = rows[j]
        proofs.append(dict(kind="original_positive_binomial_sum_zero",
            original_scalar_blocks=combination, exact_identity_checked=True))
    if proofs:
        progress("verified_original_scalar_faces", dict(identities=len(proofs),
            forced_equalities=len(derived)))
    return list(derived.values()), proofs


def symbolic_column_kernels(matrix):
    """Exact constant column dependencies, triangular in the original order."""
    pivots, result = {}, []
    for j in range(len(matrix)):
        column = {(i,v):c for i,row in enumerate(matrix) for v,c in row[j].items()}
        representation = {j:F(1)}
        for pivot,(basis,coefficients) in pivots.items():
            if pivot not in column: continue
            factor = column[pivot]
            column = conic.plus(column,conic.times(-factor,basis))
            representation = conic.plus(representation,conic.times(-factor,coefficients))
        if not column:
            direction = sorted(representation.items())
            if any(conic.plus(*(conic.times(c,row[k]) for k,c in direction)) for row in matrix):
                raise ValueError("An exact symbolic column relation failed its identity check")
            result.append(direction)
        else:
            pivot = min(column)
            scale = column[pivot]
            pivots[pivot] = (conic.times(1/scale,column),conic.times(1/scale,representation))
    return result


def prepare(problem, progress, facial=True, linear_faces=True, seconds=None):
    """Equivalent exact reductions: equalities, PSD zero rows, scalar faces.

    A scalar face is accepted only after checking a nonnegative combination
    of existing scalar cone rows/PSD diagonals is identically zero.
    Every positive-weight summand must then vanish at any feasible point.
    A PSD zero diagonal forces its entire row to vanish.
    """
    started = time.monotonic()
    initial_rows, initial_proofs = original_binomial_faces(problem, progress) if facial and linear_faces else ([], [])
    current, substitutions = conic.eliminate(problem, initial_rows)
    progress("exact_equalities_eliminated",dict(free_variables=len(current["variables"])))
    maps = [dict(original_block=b, indices=list(range(len(block["entries"]))), scale=F(1))
        for b, block in enumerate(problem["blocks"])]
    trace = [dict(stage="original_scalar_faces", proofs=initial_proofs)] if initial_proofs else []
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
            for b, block in enumerate(current["blocks"]):
                M = block["entries"]
                if len(M) > 1:
                    directions = symbolic_column_kernels(M)
                    zero_directions.setdefault(b,[]).extend(directions)
                    proofs.extend(dict(kind="PSD_symbolic_column_kernel",block=b,direction=d)
                        for d in directions)
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
                    # All recorded directions express their largest
                    # column through earlier ones; their union preserves
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
        if seconds is not None and time.monotonic()-started >= seconds:
            progress("exact_presolve_budget_checkpoint",dict(remaining_variables=len(current["variables"])))
            break
    return current, substitutions, maps, trace


def reduction_checks():
    toy = conic.model(["x","y"],[dict(name="three_column_kernel",entries=[
        [conic.form(x=1),{},conic.form(x=1)],
        [{},conic.form(y=1),conic.form(y=1)],
        [conic.form(x=1),conic.form(y=1),conic.form(x=1,y=1)]])])
    reduced,_,_,_ = prepare(toy,lambda *_:None,linear_faces=False)
    if len(reduced["blocks"][0]["entries"]) != 2:
        raise ValueError("The full symbolic PSD column kernel was not removed")
    names = ("u", "v", "w", "a", "b", "c")
    toy = conic.model(list(names), [dict(name="unit_"+v, entries=[[conic.form(**{v: 1})]])
        for v in ("u", "v", "w")]+[dict(name="inactive_sum", entries=[[
            conic.form(u=-1, v=-1, w=-1, a=1, b=1, c=1)]])]+
        [dict(name="gap_"+v, entries=[[conic.form(**{v: 1, z: -1})]])
            for v, z in zip(("u", "v", "w"), ("a", "b", "c"), strict=True)])
    rows, proofs = original_binomial_faces(toy, lambda *_: None)
    if len(rows) != 4 or len(proofs) != 1:
        raise ValueError("The exact binomial guard face was not exposed")
    if any(conic.value(q, dict(u=2, v=3, w=5, a=2, b=3, c=5)) for q in rows):
        raise ValueError("The exact binomial face changed valid nonzero populations")
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
    toy = conic.model(["x"], [dict(name="proportional_kernel", entries=[
        [conic.form(x=1), conic.form(x=3)], [conic.form(x=3), conic.form(x=9)]])])
    compiled, substitutions, _, _ = prepare(toy, lambda *_: None, linear_faces=False)
    if len(compiled["blocks"]) != 1 or len(compiled["blocks"][0]["entries"]) != 1:
        raise ValueError("The rational proportional-column PSD kernel was not reduced")
    conic.check_point(toy, {v: conic.value(q, {compiled["variables"][0]: F(1)}) for v, q in substitutions.items()})
    return dict(exact_zero_diagonal_face=True, exact_noncoordinate_kernel_face=True,exact_multi_column_kernel=True,
        exact_rational_proportional_column_face=True, exact_original_binomial_guard_face=True)


def native_array(path, field):
    text = Path(path).read_text()
    if len(text) > 128*1024**2:
        raise ValueError("Oversized native SDP answer")
    match = re.search(r"\b"+field+r"\s*=\s*\{", text)
    if not match: raise ValueError("No native "+field)
    start, depth = match.end()-1, 0
    for end in range(start, len(text)):
        depth += int(text[end] == "{")-int(text[end] == "}")
        if not depth: break
    else: raise ValueError("Unterminated native cone array")
    return text[start:end+1]


def native_cones(path, problem, field="yMat"):
    values = [F(s) for s in conic.NUMBER.findall(native_array(path, field))]
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


def solve_sdpa(problem, out, seconds, progress, resume_dir=None):
    raise ValueError('Retired solver path: use the verified preparation pipeline and run_sdpa_checkpoint_20261125.py')


def run(root, output, input_path, solver, orientation, seconds=600, facial=True, resume_from=None,prepared_from=None):
    raise ValueError('Retired solver path: use the verified preparation pipeline and run_sdpa_checkpoint_20261125.py')


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
