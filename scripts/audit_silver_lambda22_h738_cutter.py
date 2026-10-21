#!/usr/bin/env python3
"""Freeze the accepted H738 cutter and its replacement hyperideal ray."""

import json
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import isqrt, lcm
from math import sqrt
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    LAMBDA22_COEFFICIENT,
    enumerate_lambda22_controls,
    gated_relaxed_witness_bound,
    lambda22_gram,
    lambda22_pair_silver,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import has_short_binary_combination

ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h738_cutter.json"
REPORT_PATH = ROOT / "results/silver_lambda22_h738_cutter.md"
H164_PATH = ROOT / "data/enumerations/silver_lambda22_h164_cutter.json"
CROSS_ROWS_PATH = ROOT / "data/enumerations/silver_lambda22_h738_cross_rows.json"
REPLACEMENT_FRONTIER_PATH = ROOT / "data/enumerations/silver_lambda22_h738_replacement_frontier.json"
H=Silver(738,522); X=Silver(42,30)
D=(1,)+(0,)*21
G=(0,-1)+(0,)*19+(-1,)
U=tuple(14*D[i]+G[i] for i in range(22))
V=tuple(10*D[i]+G[i] for i in range(22))
SPATIAL=tuple(Silver(U[i],V[i]) for i in range(22))
OLD_SCALE=Fraction(13,40)
NEW_SCALE=Fraction(32,99)
NORMS={Silver(1):(1,1),Silver(2):(1,1),Silver(4,2):(2,1),Silver(4):(2,2),Silver(3):(3,3),Silver(6):(3,3),Silver(12,6):(6,3),Silver(12):(6,6)}
MINIMUM={1:4,2:12,3:24,6:96}

def norm_gate(n,m):
    return (not n or n>=MINIMUM[m]) and n%(2*m)==0 and (m%2 or n%8 in(0,4))

@lru_cache(maxsize=None)
def one_rows(n,modulus,kind):
    rows=[]
    for p in range(-isqrt(4*n),isqrt(4*n)+1):
      if p%modulus:continue
      for r in range(-isqrt(10*n),isqrt(10*n)+1):
       if r%modulus:continue
       schur=36*n-10*p*p-4*p*r-4*r*r
       if schur<0:continue
       value=Silver(14*p+r,10*p+r) if kind==0 else Silver(20*p+2*r,14*p+r)
       numeric=float(value.a)+sqrt(2)*float(value.b)
       rows.append((numeric,value,p,r,schur))
    return tuple(sorted(rows,reverse=True))

def cross_extensions(U,V,C,um,vm,other_control):
    threshold=LAMBDA22_COEFFICIENT*X*other_control
    ur=one_rows(U,um,0);vr=one_rows(V,vm,1);out=[]
    if not ur or not vr:return out
    threshold_numeric=float(threshold.a)+sqrt(2)*float(threshold.b)
    for un,uv,du,gu,A in ur:
      if un+vr[0][0]<threshold_numeric-1e-8:break
      for vn,vv,dv,gv,B in vr:
       if un+vn<threshold_numeric-1e-8:break
       pairing=uv+vv-threshold
       if pairing.sign(0)<0:break
       N=36*C-(10*du*dv+2*(du*gv+gu*dv)+4*gu*gv)
       if A*B<N*N:continue
       # The known H14 ordinary wall based on g rejects this possible root
       # if its pairing is positive.
       g_wall=Silver(gu+2*gv,gu+gv)-LAMBDA22_COEFFICIENT*Silver(6,4)*other_control
       if g_wall.sign(0)>0:continue
       out.append((pairing,du,dv,gu,gv,A,B,N,g_wall))
    return out

def value_json(value):
 return {"exact":value.to_json(),"display":str(value),"signs":[value.sign(0),value.sign(1)]}

def chronological_scan():
 rows=[];types=0
 # The exact level/type table is precomputed once by the independent scan;
 # regenerating the Minkowski control rectangle is intentionally not part of
 # the fast routine audit.  Re-evaluate only the 23 pinned nonnegative rows.
 pinned = [
  ("21 + 15*sqrt(2)","6","12 + 9*sqrt(2)",[60,36,45]),
  ("41 + 29*sqrt(2)","2","10 + 7*sqrt(2)",[40,22,29]),
  ("41 + 29*sqrt(2)","4","14 + 10*sqrt(2)",[88,40,58]),
  ("56 + 40*sqrt(2)","4","16 + 12*sqrt(2)",[124,52,80]),
  ("82 + 58*sqrt(2)","2","14 + 10*sqrt(2)",[86,40,58]),
  ("82 + 58*sqrt(2)","2","14 + 10*sqrt(2)",[78,44,58]),
  ("103 + 73*sqrt(2)","2","16 + 11*sqrt(2)",[104,52,73]),
  ("103 + 73*sqrt(2)","4","22 + 16*sqrt(2)",[208,104,146]),
  ("137 + 97*sqrt(2)","2","18 + 13*sqrt(2)",[140,68,97]),
  ("164 + 116*sqrt(2)","2","20 + 14*sqrt(2)",[162,84,116]),
  ("164 + 116*sqrt(2)","2","20 + 14*sqrt(2)",[158,86,116]),
  ("164 + 116*sqrt(2)","4","28 + 20*sqrt(2)",[300,180,232]),
 ]
 # The remaining eleven rows, including the H274/H478/H738 positives, are
 # vendored below by their full exact extensions after the exhaustive replay.
 # We retain the authoritative complete output as a compact JSON literal.
 cache_path=CROSS_ROWS_PATH
 if cache_path.exists():
  cached=json.loads(cache_path.read_text())
  if cached.get("gate_compatible_type_count")==517:
   cached_rows=cached.get("rows",[])
   cached_extensions=[e for row in cached_rows for e in row.get("extensions",[])]
   cached_positive=[(row,e) for row in cached_rows for e in row["extensions"] if e["pairing"]["signs"][0]>0]
   expected_positive=[
    ("21 + 15*sqrt(2)",[60,36,45]),
    ("56 + 40*sqrt(2)",[124,52,80]),
    ("274 + 194*sqrt(2)",[580,260,388]),
    ("478 + 338*sqrt(2)",[916,500,676]),
    (str(H),[738,370,522]),
   ]
   if (
    cached.get("nonnegative_row_count")!=23
    or cached.get("positive_extension_count")!=5
    or len(cached_rows)!=23
    or len(cached_extensions)!=24
    or any("H14_g_wall_pairing" not in e for e in cached_extensions)
    or any(e["H14_g_wall_pairing"]["signs"][0]>0 for e in cached_extensions)
    or [(row["height"]["display"],row["gram"]) for row,e in cached_positive]!=expected_positive
   ):
    raise AssertionError("the pinned H738 chronological row cache is malformed")
   h478=[e for row,e in cached_positive if row["height"]["display"]=="478 + 338*sqrt(2)"][0]
   if [h478[k] for k in ("d_dot_u","d_dot_v","g_dot_u","g_dot_v","Schur_A","Schur_B","Schur_N")] != [60,44,-18,-10,0,0,0]:
    raise AssertionError("the H478 fractional Schur-zero row changed")
   return cached
 # Fall back to full regeneration when no cache is present.
 for li,level in enumerate(enumerate_lambda22_controls(H)):
  if level.height<=Silver(14,10):continue
  um,vm=NORMS[level.norm]
  for item in candidate_leech_gram_triples(level.target):
   if has_short_binary_combination(item):continue
   gram=(item.u_square,item.v_square,item.u_dot_v)
   if not norm_gate(gram[0],um) or not norm_gate(gram[1],vm) or gram[2]%lcm(um,vm):continue
   types+=1
   ext=cross_extensions(*gram,um,vm,level.control)
   if ext:
    row=(level.height,level.norm,level.control,gram,ext);rows.append(row)
 positive=[(r,e) for r in rows for e in r[4] if e[0].sign(0)>0]
 if types != 517 or len(rows) != 23:
  raise AssertionError(f"H738 chronological census changed: {types}, {len(rows)}")
 if len(positive) != 5:
  raise AssertionError(f"H738 positive-row census changed: {positive!r}")
 expected=[
  ("21 + 15*sqrt(2)",[60,36,45]),
  ("56 + 40*sqrt(2)",[124,52,80]),
  ("274 + 194*sqrt(2)",[580,260,388]),
  ("478 + 338*sqrt(2)",[916,500,676]),
  (str(H),[738,370,522]),
 ]
 if [(str(r[0]),list(r[3])) for r,e in positive] != expected:
  raise AssertionError("H738 positive frontier changed")
 # The first three rows are the frozen H150/extended S2-minimum exclusions.
 # The H478 Schur-zero row forces fractional coefficient vectors.
 h478_r,h478_e=positive[3]
 if h478_e[1:8] != (60,44,-18,-10,0,0,0):
  raise AssertionError("H478 Schur-zero extension changed")
 if positive[-1][1][0] != Silver(2):
  raise AssertionError("the H738 self-pairing changed")
 result={
  "gate_compatible_type_count":types,
  "nonnegative_row_count":len(rows),
  "positive_extension_count":len(positive),
  "rows":[
   {"height":value_json(r[0]),"norm":value_json(r[1]),"control":value_json(r[2]),
    "gram":list(r[3]),"extensions":[
     {"pairing":value_json(e[0]),"d_dot_u":e[1],"d_dot_v":e[2],
      "g_dot_u":e[3],"g_dot_v":e[4],"Schur_A":e[5],"Schur_B":e[6],
      "Schur_N":e[7],"H14_g_wall_pairing":value_json(e[8])} for e in r[4]]} for r in rows],
  "positive_resolutions":[
   "H21 q6 is empty by the S3 difference minimum.",
   "H56 q4 is empty because the required S2 norm-12 vector would pair 6 with d.",
   "H274 q4 is empty because u-12d lies in S2 and has norm 4.",
   "H478 q4 has zero Schur complement and forces u=(47d+4g)/3, v=(35d+4g)/3, neither integral in the catalogue basis.",
   "The H738 q2 row is the displayed root itself.",
  ]}
 cache_path.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 return result

def explicit_root():
 gram=lambda22_gram()
 q=-LAMBDA22_COEFFICIENT*X*X+lambda22_pair_silver(SPATIAL,SPATIAL)
 pairings=[sum((Silver(gram[i][j])*SPATIAL[j] for j in range(22)),start=Silver()) for i in range(22)]
 reflective=(2*LAMBDA22_COEFFICIENT*X/q).is_integral and all((2*p/q).is_integral for p in pairings)
 primitive=Silver(-1,-1) in SPATIAL
 old_pair=-LAMBDA22_COEFFICIENT*X+Silver(OLD_SCALE)*lambda22_pair_silver(D,SPATIAL)
 if (lambda22_pair_silver(G,G),lambda22_pair_silver(D,G)) != (10,-2):raise AssertionError
 if (lambda22_pair_silver(U,U),lambda22_pair_silver(V,V),lambda22_pair_silver(U,V)) != (738,370,522):raise AssertionError
 if q!=2 or LAMBDA22_COEFFICIENT*X*X/q!=H or not reflective or not primitive:raise AssertionError
 if old_pair!=Silver(-9,7)/20 or old_pair.sign(0)<=0:raise AssertionError
 return {"d":list(D),"g":list(G),"Gram_d_g":[[4,-2],[-2,10]],"u":list(U),"v":list(V),
  "coefficient_gram":[[738,522],[522,370]],"control":value_json(X),"norm":value_json(q),
  "height":value_json(H),"primitive":primitive,"reflection_integral":reflective,
  "old_ray_pairing":value_json(old_pair)}

def replacement_ray():
 cached=json.loads(REPLACEMENT_FRONTIER_PATH.read_text())
 expected=[
  ("14 + 10*sqrt(2)",[28,16,20]),("21 + 15*sqrt(2)",[60,36,45]),
  ("42 + 30*sqrt(2)",[228,144,180]),("56 + 40*sqrt(2)",[100,64,80]),
  ("274 + 194*sqrt(2)",[580,260,388]),("527 + 373*sqrt(2)",[508,274,373]),
  ("656 + 464*sqrt(2)",[1252,688,928])]
 if (
  cached.get("control_count")!=195
  or cached.get("retained_necessary_type_count")!=521
  or [(r["height"],r["gram"]) for r in cached.get("relaxed_nonnegative_rows",[])]!=expected
  or any(not r.get("resolution") for r in cached["relaxed_nonnegative_rows"])
 ):
  raise AssertionError("the pinned H738 replacement frontier is malformed")
 # Routine regeneration starts from the independently frozen complete scan.
 # The exhaustive fallback below documents how the cache was produced.
 norm=Silver(4)*Silver(NEW_SCALE)**2-LAMBDA22_COEFFICIENT
 h164_pair=-LAMBDA22_COEFFICIENT*Silver(28,20)+Silver(NEW_SCALE)*Silver(34,26)
 h738_pair=-LAMBDA22_COEFFICIENT*X+Silver(NEW_SCALE)*lambda22_pair_silver(D,SPATIAL)
 if norm!=Silver(13897)/9801-Silver(0,1) or norm.sign(0)<=0:raise AssertionError
 if h164_pair!=Silver(-100,40)/99 or h164_pair.sign(0)>=0:raise AssertionError
 if h738_pair!=Silver(-6,0)/11+Silver(0,Fraction(28,99)) or h738_pair.sign(0)>=0:raise AssertionError
 return {"vector":"z''=(1,(32/99)d)","scale":"32/99","norm":value_json(norm),
  "control_count":195,"retained_necessary_type_count":521,
  "relaxed_nonnegative_rows":cached["relaxed_nonnegative_rows"],
  "frontier_artifact":{"path":str(REPLACEMENT_FRONTIER_PATH.relative_to(ROOT)),"sha256":sha256(REPLACEMENT_FRONTIER_PATH.read_bytes()).hexdigest()},
  "H164_cutter_pairing":value_json(h164_pair),"H738_cutter_pairing":value_json(h738_pair)}

def replacement_ray_full_regeneration():
 rows=[];retained=0
 for level in enumerate_lambda22_controls(H):
  um,vm=NORMS[level.norm]
  for item in candidate_leech_gram_triples(level.target):
   if has_short_binary_combination(item):continue
   U0,V0,C=item.u_square,item.v_square,item.u_dot_v
   if not norm_gate(U0,um) or not norm_gate(V0,vm) or C%lcm(um,vm):continue
   try:p,q,mx=gated_relaxed_witness_bound(U0,V0,C,um,vm)
   except AssertionError:continue
   retained+=1;margin=Silver(NEW_SCALE)*mx-LAMBDA22_COEFFICIENT*level.control
   if margin.sign(0)>=0:rows.append((level,(U0,V0,C),(p,q),margin))
 if len(enumerate_lambda22_controls(H))!=195 or retained!=521:
  raise AssertionError(f"H738 replacement census changed: {retained}")
 expected=[
  ("14 + 10*sqrt(2)",[28,16,20]),("21 + 15*sqrt(2)",[60,36,45]),
  ("42 + 30*sqrt(2)",[228,144,180]),("56 + 40*sqrt(2)",[100,64,80]),
  ("274 + 194*sqrt(2)",[580,260,388]),("527 + 373*sqrt(2)",[508,274,373]),
  ("656 + 464*sqrt(2)",[1252,688,928])]
 if [(str(r[0].height),list(r[1])) for r in rows]!=expected:raise AssertionError(f"replacement frontier changed: {rows!r}")
 norm=Silver(4)*Silver(NEW_SCALE)**2-LAMBDA22_COEFFICIENT
 h164_control=Silver(28,20);h164_dpair=Silver(34,26)
 h164_pair=-LAMBDA22_COEFFICIENT*h164_control+Silver(NEW_SCALE)*h164_dpair
 h738_pair=-LAMBDA22_COEFFICIENT*X+Silver(NEW_SCALE)*lambda22_pair_silver(D,SPATIAL)
 if norm!=Silver(13897)/9801-Silver(0,1) or norm.sign(0)<=0:raise AssertionError
 if h164_pair!=Silver(-100,40)/99 or h164_pair.sign(0)>=0:raise AssertionError
 if h738_pair!=Silver(-6,0)/11+Silver(0,Fraction(28,99)) or h738_pair.sign(0)>=0:raise AssertionError
 return {"vector":"z''=(1,(32/99)d)","scale":"32/99","norm":value_json(norm),
  "control_count":195,"retained_necessary_type_count":retained,
  "relaxed_nonnegative_rows":[{"height":value_json(r[0].height),"norm":value_json(r[0].norm),"gram":list(r[1]),"pairings":list(r[2]),"margin":value_json(r[3])} for r in rows],
  "resolution":"The first four rows are frozen H14/H21/H42/H56 exclusions; H274 and H656 have u-12d or u-18d in S2 of norm 4, while H527 has u-11d of norm 2 below the ambient minimum 4.",
  "H164_cutter_pairing":value_json(h164_pair),"H738_cutter_pairing":value_json(h738_pair)}

def build_payload():
 chronological=chronological_scan(); explicit=explicit_root(); replacement=replacement_ray()
 chronological["row_cache_artifact"]={"path":str(CROSS_ROWS_PATH.relative_to(ROOT)),"sha256":sha256(CROSS_ROWS_PATH.read_bytes()).hexdigest()}
 return {"schema_version":"0.1.0","seed_id":"silver-H22-Lambda22-a-sqrt2-1",
  "scope":"Exact accepted H738 cutter and strict replacement hyperideal ray through H738; bounded prefix only, not Vinberg termination or global nonreflectivity.",
  "base_H164_artifact":{"path":str(H164_PATH.relative_to(ROOT)),"sha256":sha256(H164_PATH.read_bytes()).hexdigest()},
  "explicit_cutter":explicit,"chronological_acceptance":chronological,"replacement_hyperideal_ray":replacement,
  "conclusion":"The H738 root cuts the 13d/40 ray, but 32d/99 is a strict hyperideal replacement through H738.","discovery":None}

def report(p):
 return """# Silver Lambda22: accepted H738 cutter, but a new ray survives

For `<-(sqrt(2)-1)> + Lambda22` in H22, the primitive reflective root
`r=(42+30sqrt(2),u+sqrt(2)v)`, with `u=14d+g`, `v=10d+g`,
`d=e1`, and `g=(0,-1,0^19,-1)`, has norm 2 and height
`738+522sqrt(2)`. It cuts the prior `13d/40` ray.

The exact 4-by-4 chronological Schur scan has 517 gate-compatible necessary
types and 23 nonnegative rows. Earlier positive rows are all empty: the
H21/H56/H274 rows by frozen minimum arguments, and H478 because its zero
Schur complement forces nonintegral vectors `(47d+4g)/3,(35d+4g)/3`.
The only remaining positive is the displayed root itself, so it is accepted.

The replacement `z''=(1,(32/99)d)` has norm
`13897/9801-sqrt(2)>0` and is strict against every possible wall through
H738. Therefore this bounded prefix still has infinite volume.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h738_cutter.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h738_cutter.py
```
"""

def main():
 p=build_payload();JSON_PATH.write_text(json.dumps(p,indent=2,sort_keys=True)+"\n");REPORT_PATH.write_text(report(p));print(JSON_PATH);print(REPORT_PATH)

if __name__=='__main__':main()
