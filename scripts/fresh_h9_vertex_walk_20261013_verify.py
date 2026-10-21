#!/usr/bin/env python3
"""Producer-free exact replay of the targeted spherical vertex walk."""
import argparse
import copy
import hashlib
import json
from itertools import combinations
from scripts.fresh_h9_vertex_walk_20261013_base import *
from scripts.fresh_h9_replaced_core_20261013_verify import census_check
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K

def actual_point_pair(row,F,y):
    common,er=embed(F,row["field"])
    return common,dot(common,tuple(er(x)for x in row["c"]),
                      tuple(lift_existing(common,F,x)for x in y))

def cap_check(record,name,G,I,h,old,target,pool):
    row=cap_profile(record,name);F=row["field"];E=F.base
    assert K.sign(E.d)>=0 and E.sign(F.d)>0
    u,v=[tuple(map(E.parse,record[k]))for k in("line_constant","line_direction")]
    IE=tuple(tuple(E.lift(x)for x in r)for r in I)
    A=pairing(E,v,IE,v);B=E.mul(E.rational(2),pairing(E,u,IE,v))
    C=E.sub(pairing(E,u,IE,u),E.rational(2))
    assert E.sign(A)
    assert all(E.sign(E.sub(x,y))==0 for x,y in zip((A,B,C),map(E.parse,record["norm_polynomial"]),strict=True))
    D=E.sub(E.mul(B,B),E.mul(E.rational(4),E.mul(A,C)))
    assert E.sign(E.sub(D,F.d))==0
    t=F.parse(record["parameter"]);branch=record["branch"];assert branch in(-1,1)
    expected=(E.div(E.neg(B),E.mul(E.rational(2),A)),E.div(E.rational(branch),E.mul(E.rational(2),A)))
    assert F.sign(F.sub(t,expected))==0
    p=tuple(F.add(F.lift(x),F.mul(t,F.lift(y)))for x,y in zip(u,v,strict=True))
    assert all(F.sign(F.sub(x,y))==0 for x,y in zip(p,row["p"],strict=True))
    Gf=tuple(tuple(lift_K(F,x)for x in r)for r in G)
    assert all(F.sign(F.sub(x,y))==0 for x,y in zip(matvec(F,Gf,row["c"]),p,strict=True))
    assert F.sign(F.sub(dot(F,row["c"],p),F.rational(2)))==0
    assert all(F.sign(x)<=0 for x in p)
    assert F.sign(dot(F,row["c"],tuple(lift_K(F,x)for x in h)))>0
    target_pair=dot(F,row["c"],tuple(lift_K(F,x)for x in target))
    assert F.sign(target_pair)<0 and F.encode(target_pair)==record["target_pairing"]
    block=record["block_names"];assert len(block)==len(set(block))==8
    assert all(old[n]["field"]is K for n in block)
    values=[Z]*8
    for j,m in zip(record["block_support"],record["block_labels"],strict=True):values[j]=K.neg({3:O,5:P}[m])
    for n,want in zip(block,values,strict=True):
        c=tuple(E.lift(x)for x in old[n]["c"])
        assert E.sign(E.sub(dot(E,c,u),E.lift(want)))==0
        assert E.sign(dot(E,c,v))==0
    anchor=old[record["anchor_name"]]
    common,ea=embed(E,anchor["field"])
    coefficients=tuple(ea(x)for x in anchor["c"])
    want=K.neg({2:Z,3:O,5:P}[record["anchor_label"]])
    assert common.sign(common.sub(dot(common,coefficients,tuple(lift_existing(common,E,x)for x in u)),lift_K(common,want)))==0
    assert common.sign(dot(common,coefficients,tuple(lift_existing(common,E,x)for x in v)))==0
    checks={n:relation(row,r,pool)for n,r in old.items()}
    assert checks==record["pair_checks"]
    assert relation(row,anchor,pool)=={"kind":"finite","label":record["anchor_label"]}
    dropped=[n for n,r in checks.items()if r["kind"]not in("finite","dotted")]
    assert dropped==record["dropped"] and not any(n.startswith("seed")for n in dropped)
    gram=[[dot(K,old[a]["c"],old[b]["p"])for b in block]for a in block]
    full=[r+[x]for r,x in zip(gram,values,strict=True)]+[values+[scale(O,2)]]
    assert full==[[parse(x)for x in r]for r in record["spherical_vertex_gram"]]
    assert inertia_K(full)==(9,0,0)
    return row

def point_check(proof,record,selected,I,h):
    name=proof["point_field_from_cap"];F=selected[name]["field"]
    y=tuple(map(F.parse,proof["y"]));IF=tuple(tuple(lift_K(F,x)for x in r)for r in I)
    q=pairing(F,y,IF,y);cp=pairing(F,y,IF,tuple(lift_K(F,x)for x in h))
    assert F.sign(q)<0 and F.sign(cp)<0
    assert F.encode(q)==proof["norm"]and F.encode(cp)==proof["center_pairing"]
    signs={}
    for n,r in selected.items():
        common,v=actual_point_pair(r,F,y);signs[n]=common.sign(v)
    assert signs==proof["all_wall_signs"]
    order=record["block_names"]+[name]
    assert len(order)==len(set(order))==9 and sorted(order)==proof["zero_walls"]
    assert all(s==(0 if n in order else 1)for n,s in signs.items())
    gram=[[parse(x)for x in r]for r in proof["gram"]]
    assert inertia_K(gram)==(9,0,0)
    for i,n in enumerate(order):
        for j,m in enumerate(order):
            common,v=pair_profiles(selected[n],selected[m])
            assert common.sign(common.sub(v,lift_K(common,gram[i][j])))==0

def run():
    G,I,h,initial,_=chamber(0);pool=angle_pool(ANGLE_SOURCE)
    assert inertia_K(G)==(9,1,0)
    assert all(matvec(K,I,r)==tuple(O if i==j else Z for j in range(10))for i,r in enumerate(G))
    assert K.sign(pairing(K,h,I,h))<0
    counts={"exact_caps":0,"configurations":0,"full_pair_checks":0,"distinct_pair_replays":0,
            "affine_frontier_checks":0,"actual_finite_vertices":0,"strict_spacelike_escapes":0}
    files=[];cache={};last=None
    specs=[(0,"caps","geometry","walk"),(1,"second_caps","second_geometry","bridge"),
           (2,"third_caps","third_geometry","step")]
    for stage,cap_suffix,geom_suffix,nameprefix in specs:
        cap_source=PREFIX+"_"+cap_suffix+".json";geom_source=PREFIX+"_"+geom_suffix+".json"
        if not(ROOT/geom_source).exists():continue
        files.extend((cap_source,geom_source))
        bank=read(cap_source);geometry=read(geom_source)
        _,_,_,old,oldcase=chamber(stage)
        target_proof=bank["target"]["certificate"]
        assert bank["target"]==oldcase["geometry"]["escapes"][0]
        target=tuple(parse(x)for x in target_proof["y"])
        assert strict_escape(target,I,h,list(old.values()))==target_proof
        rows=[]
        for i,record in enumerate(bank["caps"]):
            name=f"{nameprefix}{i}"
            row=cap_check(record,name,G,I,h,old,target,pool)
            selected={n:r for n,r in old.items()if n not in record["dropped"]}|{name:row}
            point_check(record["finite_vertex"],record,selected,I,h)
            rows.append(row);counts["exact_caps"]+=1;counts["actual_finite_vertices"]+=1
        newpairs={f"{i}|{j}":relation(rows[i],rows[j],pool)for i,j in combinations(range(len(rows)),2)}
        assert newpairs==geometry["new_pair_relations"]
        for case in geometry["cases"]:
            ids=case["selected_caps"]
            dropped=sorted(set().union(*(set(bank["caps"][i]["dropped"])for i in ids)))
            assert dropped==case["dropped"]
            selected={n:r for n,r in old.items()if n not in dropped}|{rows[i]["name"]:rows[i]for i in ids}
            assert list(selected)==case["walls"]and len(selected)==case["facet_count"]
            checks={}
            for a,b in combinations(selected.values(),2):
                kk="|".join(sorted((a["name"],b["name"])))
                if kk not in cache:cache[kk]=relation(a,b,pool)
                checks[kk]=cache[kk]
            assert checks==case["pair_checks"]and all(r["kind"]in("finite","dotted")for r in checks.values())
            counts["full_pair_checks"]+=len(checks)
            counts["affine_frontier_checks"]+=census_check(list(selected),checks,case["affine_census"])
            assert len(case["actual_finite_vertices"])==len(ids)
            for i,proof in zip(ids,case["actual_finite_vertices"],strict=True):
                point_check(proof,bank["caps"][i],selected,I,h)
                counts["actual_finite_vertices"]+=1
            if"geometry"in case:
                assert case["status"]=="EXACT_STRICT_ESCAPES"and not case["affine_census"]["affine_circuits"]
                for escape in case["geometry"]["escapes"]:
                    cert=escape["certificate"]
                    assert strict_escape(tuple(map(parse,cert["y"])),I,h,list(selected.values()))==cert
                    counts["strict_spacelike_escapes"]+=1
            else:
                assert case["status"]=="REJECTED_POSITIVE_AFFINE_COMPACTNESS"and case["affine_census"]["affine_circuits"]
            counts["configurations"]+=1
        last=(bank,old,target)
    counts["distinct_pair_replays"]=len(cache)
    bank,old,target=last
    forged=copy.deepcopy(bank["caps"][0]);row=cap_profile(forged,"forged");F=row["field"]
    x=F.add(F.parse(forged["seed_coefficients"][0]),F.one)
    forged["seed_coefficients"][0]=F.encode(x)
    try:cap_check(forged,"forged",G,I,h,old,target,pool)
    except AssertionError:pass
    else:raise AssertionError("Forged cap accepted")
    try:strict_escape(h,I,h,list(initial.values()))
    except AssertionError:pass
    else:raise AssertionError("Timelike escape accepted")
    return {"status":"PASS_EXACT_TARGETED_VERTEX_WALK","counts":counts,
            "source_sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in files},
            "mutations":{"changed_cap_coefficient_rejected":True,"timelike_escape_rejected":True},
            "scope":"Every saved new cap and claimed finite vertex, all selected chamber pair gates and complete affine-frontier replays, plus exact future strict spacelike escapes. No exhaustive continuous cap-family or universal finite-angle exclusion claim."}

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--verify",action="store_true");args=parser.parse_args()
    result=run();out=ROOT/(PREFIX+"_verify.json")
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
