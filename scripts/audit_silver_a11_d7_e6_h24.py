#!/usr/bin/env python3
"""Freeze the exact silver-unit N(A11 D7 E6) two-packet prefix."""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_a11_d7_e6 import *

ROOT=Path(__file__).resolve().parents[1];OUTPUT=ROOT/"data/enumerations/silver_a11_d7_e6_h24_h15.json"
def silver(x):return Silver.coerce(x).to_json()
def vec(x):return [silver(v) for v in x]
def weight(w):return {"square":silver(w.square),"weights":vec(w.weights),"coordinates":vec(w.coordinates),"residue":list(w.residue)}
def root(r):return {"norm":silver(r.norm),"height":silver(r.height),"control":silver(r.control),"target":silver(r.target),"a11":weight(r.a),"d7":weight(r.d),"e6":weight(r.e),"glue_nontrivial":r.glue_nontrivial,"three_factor_glue":r.three_factor_glue}
def main():
    if not basis_audit():raise ArithmeticError("basis failed")
    prefix=run_prefix(15);nonempty=[p for p in prefix.packets if p.roots]
    if [len(p.roots) for p in nonempty] != [3,35]:raise ArithmeticError("packets changed")
    second=nonempty[1]
    if sum(r.glue_nontrivial for r in second.roots)!=22 or sum(r.three_factor_glue for r in second.roots)!=12:raise ArithmeticError("glue census changed")
    witness=strict_witness()
    if not witness.strict or not witness.hyperideal:raise ArithmeticError("witness failed")
    roots=[root(r) for r in prefix.accepted_roots]
    geometry=first_geometry()
    payload={"scope":{"field":"Q(sqrt(2))","hyperbolic_dimension":24,"spatial_lattice":"Niemeier N(A11 D7 E6), cyclic glue <(1,1,1)>","timelike_coefficient":silver(EPSILON),"height_bound":15,"claim":"Complete first two nonempty packets and strict hyperideal witness; bounded, not global."},"glue":{"discriminant_factors":[12,4,3],"generator":[1,1,1],"generator_weights":["omega_1(A11)","omega_6(D7)","omega_1(E6)"],"generator_square":"4","order":12,"explicit_basis_even_unimodular":True,"height_zero_orbits":len(height_zero()),"height_zero_roots":288},"mirror_norms":{"canonical":[silver(1),silver(2)],"proof":COMPLETE_UNIT_NORM_PROOF},"packets":[{"height":silver(p.height),"controls":len(p.controls),"candidates":p.candidate_count,"rejected":p.rejected_count,"accepted":len(p.roots),"glue_nontrivial":sum(r.glue_nontrivial for r in p.roots),"three_factor_glue":sum(r.three_factor_glue for r in p.roots)} for p in prefix.packets],"accepted_roots":roots,"accepted_roots_sha256":sha256(json.dumps(roots,sort_keys=True,separators=(",",":")).encode()).hexdigest(),"first_geometry":{"rays":geometry.rays,"finite":geometry.finite,"cusps":geometry.cusps,"hyperideal":geometry.hyperideal},"witness":{"a11_node_zero_based":0,"weight":str(witness.weight),"norm":silver(witness.norm),"strict":witness.strict,"hyperideal":witness.hyperideal,"closest_pairing":silver(max(witness.pairings)),"pairing_histogram":dict(sorted(Counter(str(v) for v in witness.pairings).items()))},"outcome":"No compact/FV prefix: 22 glue walls, 12 coupling all factors, but strict hyperideal A11-axis ray survives."}
    OUTPUT.parent.mkdir(parents=True,exist_ok=True);OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n");print(json.dumps({"output":str(OUTPUT),"sha256":sha256(OUTPUT.read_bytes()).hexdigest(),"accepted":len(roots)},indent=2))
if __name__=="__main__":main()
