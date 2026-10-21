#!/usr/bin/env python3
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_d10_e7_squared import *
ROOT=Path(__file__).resolve().parents[1];OUTPUT=ROOT/"data/enumerations/silver_d10_e7_squared_h24_h15.json"
def silver(x):return Silver.coerce(x).to_json()
def vec(x):return [silver(v) for v in x]
def wp(w):return {"square":silver(w.square),"weights":vec(w.weights),"coordinates":vec(w.coordinates),"residue":w.residue}
def rp(r):return {"norm":silver(r.norm),"height":silver(r.height),"control":silver(r.control),"target":silver(r.target),"d10":wp(r.a),"e7a":wp(r.d),"e7b":wp(r.e),"glue_nontrivial":r.glue_nontrivial,"three_factor_glue":r.three_factor_glue}
def main():
 p=run_prefix(15);nonempty=[x for x in p.packets if x.roots];w=strict_witness();g=first_geometry()
 if not basis_audit() or [len(x.roots) for x in nonempty]!=[3,16] or sum(r.glue_nontrivial for r in nonempty[1].roots)!=7 or sum(r.three_factor_glue for r in nonempty[1].roots)!=3 or not(w.strict and w.hyperideal):raise ArithmeticError("frozen D10E7^2 audit failed")
 roots=[rp(r) for r in p.accepted_roots];payload={"scope":{"field":"Q(sqrt(2))","dimension":24,"spatial_lattice":"Niemeier N(D10 E7^2), binary graph glue","height_bound":15,"claim":"Complete first two nonempty packets; strict hyperideal witness; bounded only."},"glue":{"code":["(0,0,0)","(s,1,0)","(c,0,1)","(v,1,1)"],"index":4,"explicit_hnf_basis_even_unimodular":True,"height_zero_orbits":len(height_zero()),"height_zero_roots":432},"norms":{"canonical":[silver(1),silver(2)],"proof":COMPLETE_UNIT_NORM_PROOF},"packets":[{"height":silver(x.height),"candidates":x.candidate_count,"rejected":x.rejected_count,"accepted":len(x.roots),"glue":sum(r.glue_nontrivial for r in x.roots),"three_factor":sum(r.three_factor_glue for r in x.roots)} for x in p.packets],"roots":roots,"roots_sha256":sha256(json.dumps(roots,sort_keys=True,separators=(",",":")).encode()).hexdigest(),"geometry":{"rays":g.rays,"finite":g.finite,"cusps":g.cusps,"hyperideal":g.hyperideal},"witness":{"block":"D10","node_zero_based":0,"weight":str(w.weight),"norm":silver(w.norm),"closest":silver(max(w.pairings)),"strict":w.strict,"hyperideal":w.hyperideal,"histogram":dict(sorted(Counter(str(v) for v in w.pairings).items()))},"outcome":"No compact/FV prefix: exact strict hyperideal ray survives all 19 walls."}
 OUTPUT.parent.mkdir(parents=True,exist_ok=True);OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n");print(json.dumps({"output":str(OUTPUT),"sha256":sha256(OUTPUT.read_bytes()).hexdigest(),"accepted":len(roots)},indent=2))
if __name__=="__main__":main()
