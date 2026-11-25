"""Retained-model nesting and a new independent equation-closure fixture."""
from fractions import Fraction as F
import json
from pathlib import Path

from . import conic, equation_closure, strengthening, white_profiles, supporting_counts


def run(root, output):
    out = (root/output).resolve()
    if out.exists() or not out.is_relative_to(root/"results"):
        raise ValueError("Use a fresh results directory")
    models = [root/"results"/name/"necessary_primal.exact.json" for name in (
        "global_level2_current","global_shared_support_current","global_closed_face_current",
        "global_profile_closure_current","global_shared_star_current","global_shared_cusp_current")]
    witness = root/"results/global_level2_face_recovery/original_feasible_point.exact.json"
    paths = [*models,witness,Path(__file__),root/"scripts/run_global_coxeter_formulation.py"]
    paths += [root/"research/global_formulation"/name for name in (
        "white_profiles.py","supporting_counts.py","equation_closure.py","strengthening.py","conic.py","moments.py","polynomial.py","faces.py","data.py")]
    def bindings():
        return {str(p.relative_to(root)):dict(bytes=p.stat().st_size,
            sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    before = bindings()
    out.mkdir()
    conic.write_json(out/"source_bindings.json",before)
    A = [[conic.form(x=1),conic.form(z=1)],[conic.form(z=1),conic.form(y=1)]]
    B = [[conic.form(w=1),{},{}],[{},*A[0]],[{},*A[1]]]
    if strengthening.principal_embedding(A,B) != [1,2]: raise ValueError("Principal PSD retention failed")
    B[1][2] = B[2][1] = conic.form(1,z=1)
    if strengthening.principal_embedding(A,B) is not None: raise ValueError("False PSD retention from diagonals only")
    result = dict(status="SHARED_PROGRAM_CHECKS_PASSED",new_bound=False,
        principal_retention_checks=dict(noninitial_principal_block=True,changed_cross_entry_rejected=True),
        closure_fixture=equation_closure.self_checks(),local_Gram_lift=white_profiles.local_star_fixture(),
        quadrilateral_Gram_lift=white_profiles.quadrilateral_star_fixture(),
        complete_family_deletions=supporting_counts.self_checks(),nesting=[])
    old = conic.read_model(models[0])
    for path in models[1:]:
        new = conic.read_model(path)
        retained = strengthening.compare_and_test(old,new,None)
        result["nesting"].append(dict(input=str(path.relative_to(root)),**retained))
        old = new
    point = {v:F(c) for v,c in json.loads(witness.read_text()).items()}
    result["level2_witness_separation"] = strengthening.compare_and_test(conic.read_model(models[0]),old,point)
    if bindings() != before: raise ValueError("A shared-program check input changed")
    result["source_bindings"] = before
    conic.write_json(out/"result.json",result)
    print(json.dumps({k:v for k,v in result.items() if k not in ("source_bindings","level2_witness_separation")}))
