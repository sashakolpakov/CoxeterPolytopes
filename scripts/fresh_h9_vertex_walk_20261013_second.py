#!/usr/bin/env python3
"""Close the next exposed end while retaining the first exact vertex cap."""
import itertools
import json
import mpmath as mp
from scripts.fresh_h9_vertex_walk_20261013_caps import *
from scripts.fresh_h9_replaced_core_20261013_cores import affine_census
from scripts.fresh_h9_replaced_core_20261013_geometry import volume_gate

def current():
    G,I,h,old,_=selected_coupled_case()
    geometry=read(PREFIX+"_geometry.json");case=geometry["cases"][0]
    bank=read(geometry["cap_source"])
    for n in case["dropped"]:old.pop(n)
    for i in case["selected_caps"]:old[f"walk{i}"]=cap_profile(bank["caps"][i],f"walk{i}")
    assert list(old)==case["walls"]
    return G,I,h,old,case

def run():
    mp.mp.dps=90
    G,I,h,old,oldcase=current()
    target=tuple(parse(x)for x in oldcase["geometry"]["escapes"][0]["certificate"]["y"])
    pool=angle_pool(ANGLE_SOURCE);scout=read(PREFIX+"_scout.json")
    caps=[];rows=[]
    for proposal_index in(3,10,6):
        row,rec=make_cap(scout["proposals"][proposal_index],I,h,old,target,pool)
        rec["proposal_index"]=proposal_index;name=f"bridge{len(caps)}";row["name"]=name
        assert "walk0"not in rec["dropped"]and len(rec["dropped"])==1
        selected={n:r for n,r in old.items()if n not in rec["dropped"]}|{name:row}
        rec["finite_vertex"]=finite_vertex(row,rec,I,h,selected)
        caps.append(rec);rows.append(row)
        print({"cap":name,"proposal":proposal_index,"dropped":rec["dropped"],"vertex":True},flush=True)
    bank={"base_geometry":PREFIX+"_geometry.json","base_case":0,"target":oldcase["geometry"]["escapes"][0],
          "caps":caps,"scope":"Three actual finite-vertex replacements cutting the next exact escape and retaining walk0."}
    (ROOT/(PREFIX+"_second_caps.json")).write_text(json.dumps(bank,indent=2)+"\n")
    pairs={(i,j):relation(rows[i],rows[j],pool)for i,j in itertools.combinations(range(len(rows)),2)}
    selections=[(i,)for i in range(len(rows))]
    selections += [ids for r in(2,3)for ids in itertools.combinations(range(len(rows)),r)
                   if all(pairs[i,j]["kind"]in("finite","dotted")for i,j in itertools.combinations(ids,2))]
    selections.sort(key=lambda x:(-len(x),x))
    cases=[]
    for ids in selections:
        dropped=sorted(set().union(*(set(caps[i]["dropped"])for i in ids)))
        selected={n:r for n,r in old.items()if n not in dropped}|{rows[i]["name"]:rows[i]for i in ids}
        checks={k:r for k,r in oldcase["pair_checks"].items()if all(n in selected for n in k.split("|"))}
        for i in ids:checks.update({"|".join(sorted((rows[i]["name"],n))):r for n,r in caps[i]["pair_checks"].items()if n in selected})
        for i,j in itertools.combinations(ids,2):checks["|".join(sorted((rows[i]["name"],rows[j]["name"])))]=pairs[i,j]
        names=list(selected);assert len(checks)==len(names)*(len(names)-1)//2
        assert all(r["kind"]in("finite","dotted")for r in checks.values())
        census=affine_census(names,checks)
        case={"selected_caps":list(ids),"dropped":dropped,"walls":names,"facet_count":len(names),
              "pair_checks":checks,"affine_census":census,
              "actual_finite_vertices":[finite_vertex(rows[i],caps[i],I,h,selected)for i in ids]}
        if census["affine_circuits"]:case["status"]="REJECTED_POSITIVE_AFFINE_COMPACTNESS"
        else:
            case["geometry"]=volume_gate(I,h,list(selected.values()))
            case["status"]=case["geometry"]["status"]
        cases.append(case)
        out={"cap_source":PREFIX+"_second_caps.json","base_geometry":PREFIX+"_geometry.json","base_case":0,
             "new_pair_relations":{f"{i}|{j}":r for (i,j),r in pairs.items()},"cases":cases,
             "scope":"Coherent selected edge-end replacements retaining the first cap. Exact pair, affine and actual-vertex gates; strict escapes certify remaining infinite volume."}
        (ROOT/(PREFIX+"_second_geometry.json")).write_text(json.dumps(out,indent=2)+"\n")
        print({"case":len(cases)-1,"selected":ids,"facets":len(names),"status":case["status"],
               "affine":len(census["affine_circuits"]),"max_norm":case.get("geometry",{}).get("maximum_numeric_vertex_norm")},flush=True)
        if case["status"]=="UNRESOLVED_VOLUME":break
if __name__=="__main__":run()
