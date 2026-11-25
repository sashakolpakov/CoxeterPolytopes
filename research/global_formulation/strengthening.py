"""Uniform cone and integer-population strengthenings of the global moment model.

These are families for every face/counter in the existing schemas. They
introduce no selected vertex forests or face-diagram catalogue.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

from . import conic
from .polynomial import as_poly, determinant, var


def compound(matrix, degree):
    subsets = list(combinations(range(len(matrix)), degree))
    return [[determinant([[matrix[i][j] for j in J] for i in I]) for J in subsets] for I in subsets]


def symbolic_gram_checks():
    """Universal polynomial checks of the rank and compound-sign formulas.

    The entries of U and W are indeterminates, not sampled configurations.
    Four columns cover both supported face sizes, including mixed minors.
    """
    U = [[var(f"U{i}_{j}") for j in range(4)] for i in range(3)]
    H = [[U[0][i]*U[0][j]+U[1][i]*U[1][j]-U[2][i]*U[2][j]
        for j in range(4)] for i in range(4)]
    triples = list(combinations(range(4), 3))
    v = [determinant([[U[i][j] for j in J] for i in range(3)]) for J in triples]
    C = compound(H, 3)
    if any(C[i][j]+v[i]*v[j] for i in range(4) for j in range(4)):
        raise ValueError("The symbolic Lorentz compound identity failed")
    W = U[:2]
    M = [[W[0][i]*W[0][j]+W[1][i]*W[1][j] for j in range(4)] for i in range(4)]
    if any(q for row in compound(M, 3) for q in row):
        raise ValueError("The symbolic spatial mixed-minor identity failed")
    return dict(universal_Lorentz_third_compound_identity=True,
        universal_spatial_mixed_rank_minors_identity=True,
        symbolic_checks_are_polynomial_identities_not_sampled_diagrams=True)


def enrich(domain, level):
    if level == 1:
        return
    if level != 2:
        raise ValueError("Available global moment levels are 1 and 2")
    t = var("t")
    for name, minimum_rank in sorted(domain.integer_counts.items()):
        u = var(name)
        factorial = as_poly(1)
        for degree in range(1, 5):
            factorial *= u-(degree-1)*t
            if degree >= 3:
                domain.ge(f"level2_factorial{degree}_{name}", factorial)
            if degree >= 2:
                domain.ge(f"level2_rank_support{degree}_{name}",
                    (domain.rank_total-degree*minimum_rank*t)*factorial)
        for k in range(4):
            gap = (u-k*t)*(u-(k+1)*t)
            if k:
                domain.ge(f"level2_integer_gap{k}_{name}", gap)
            domain.matrix(f"level2_integer_localizer{k}_{name}",
                [[gap, u*gap], [u*gap, u*u*gap]])
    spatial = next((M for name, M in domain.matrices if name == "common_spatial_Gram"), None)
    if spatial is None:
        return
    sides = len(spatial)
    h = [var(f"h{i}") for i in range(sides)]
    H = [[spatial[i][j]-h[i]*h[j] for j in range(sides)] for i in range(sides)]
    C3 = compound(spatial, 3)
    for i in range(len(C3)):
        for j in range(i):
            domain.eq(f"level2_all_spatial_rank_minors_{i}_{j}", C3[i][j])
    negative_C3H = [[-q for q in row] for row in compound(H, 3)]
    domain.matrix("level2_negative_Lorentz_third_compound", negative_C3H)
    tensor = [[spatial[i][k]*spatial[j][l] for k in range(sides) for l in range(sides)]
        for i in range(sides) for j in range(sides)]
    domain.matrix("level2_spatial_tensor_square", tensor)
    for i in range(sides):
        c = var(f"c{i}")
        domain.matrix(f"level2_spatial_localizer_cosine{i}", [[c*q for q in row] for row in spatial])
        domain.matrix(f"level2_spatial_localizer_gap{i}", [[(1-c)*q for q in row] for row in spatial])


def compare_and_test(base, strengthened, old_point):
    """Check literal retention and separate an old feasible point where possible.

    Moment names are compiler slots; compare using their measure/monomial
    identities. A violation involving only old moments forbids every extension
    of the old witness, rather than just one guessed higher-moment completion.
    """
    def registry(problem):
        return {name: (measure, tuple(monomial))
            for measure, moments in problem["metadata"]["moment_coordinates"].items()
            for monomial, name in moments}
    before, after = registry(base), registry(strengthened)
    inverse = {key: name for name, key in before.items()}
    if {k: F(v) for k, v in base["metadata"]["measure_mass_scales"].items()} != {
            k: F(v) for k, v in strengthened["metadata"]["measure_mass_scales"].items()}:
        raise ValueError("Moment measure scaling changed")
    abstract = {v: repr(key) for v, key in before.items()}
    abstract_new = {v: repr(key) for v, key in after.items()}
    def rename(q, names): return {names[v] if v else "": c for v, c in q.items()}
    new_eq = {conic.canonical(rename(q, abstract_new)) for q in strengthened["equations"]}
    new_blocks = {conic.canonical([[rename(q, abstract_new) for q in row] for row in b["entries"]])
        for b in strengthened["blocks"]}
    if any(conic.canonical(rename(q, abstract)) not in new_eq for q in base["equations"]):
        raise ValueError("An original equality is absent after strengthening")
    if any(conic.canonical([[rename(q, abstract) for q in row] for row in b["entries"]]) not in new_blocks
            for b in base["blocks"]):
        raise ValueError("An original cone constraint is absent after strengthening")
    available = {v: old_point[inverse[key]] for v, key in after.items() if key in inverse}
    equality_violations, cone_violations, undecided = [], [], 0
    for i, q in enumerate(strengthened["equations"]):
        if any(v and v not in available for v in q):
            undecided += 1; continue
        value = conic.value(q, available)
        if value:
            equality_violations.append(dict(label=strengthened["metadata"]["equality_labels"][i], value=value))
    for b in strengthened["blocks"]:
        if any(v and v not in available for row in b["entries"] for q in row for v in q):
            undecided += 1; continue
        M = [[conic.value(q, available) for q in row] for row in b["entries"]]
        try:
            conic.psd_pivots(M)
        except ValueError as error:
            item = dict(name=b["name"], order=len(M), exact_failure=str(error))
            if len(M) == 1: item["value"] = M[0][0]
            cone_violations.append(item)
    return dict(all_original_equalities_retained=len(base["equations"]),
        all_original_cones_retained=len(base["blocks"]),
        old_witness_has_no_extension=bool(equality_violations or cone_violations),
        equality_violations=equality_violations, cone_violations=cone_violations,
        constraints_needing_new_moments=undecided,
        old_moments=len(before), shared_moments=len(available), strengthened_moments=len(after))


def audit(root, output, input_path, point_path):
    """Bind both compared models and the previously certified baseline witness."""
    out = (root/output).resolve()
    source, witness = (root/input_path).resolve(), (root/point_path).resolve()
    base = root/"results/global_tail57_current/necessary_primal.exact.json"
    prior = root/"results/global_current_independent_feasibility/result.json"
    if out.exists() or not out.is_relative_to(root/"results"):
        raise ValueError("Use a fresh results directory")
    if any(not p.is_relative_to(root/"results") for p in (source, witness)):
        raise ValueError("Use retained exact results inputs")
    paths = [source, base, witness, prior, Path(__file__),
        root/"research/global_formulation/conic.py", root/"research/global_formulation/polynomial.py",
        root/"scripts/run_global_coxeter_formulation.py"]
    def bindings():
        return {str(p.relative_to(root)): dict(bytes=p.stat().st_size,
            sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    before = bindings()
    inherited = json.loads(prior.read_text())
    for p in (base, witness):
        key = str(p.relative_to(root))
        if before[key] != inherited["source_bindings"][key]:
            raise ValueError("A baseline input differs from its completed independent verification")
    point = {v: F(c) for v, c in json.loads(witness.read_text()).items()}
    if conic.digest_bytes(conic.canonical(point)) != inherited["original_point_canonical_sha256"]:
        raise ValueError("The baseline witness's canonical digest changed")
    result = compare_and_test(conic.read_model(base), conic.read_model(source), point)
    result.update(symbolic_gram_checks())
    result.update(status="EXACT_STRICT_STRENGTHENING_OVER_RETAINED_BASELINE" if
        result["old_witness_has_no_extension"] else "NO_SEPARATION_FROM_OLD_WITNESS",
        baseline_feasibility_uses_its_completed_independent_verification=True,
        accepted_bounds_changed=False, source_bindings=before)
    if bindings() != before:
        raise ValueError("A strengthening audit input changed during the run")
    out.mkdir()
    conic.write_json(out/"source_bindings.json", before)
    conic.write_json(out/"result.json", result)
    print(json.dumps(conic.plain({k: v for k, v in result.items() if k not in
        ("source_bindings", "equality_violations", "cone_violations")})), flush=True)
    return result
