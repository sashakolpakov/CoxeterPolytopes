#!/usr/bin/env python3
"""Coherent replacements by exact spherical edge-end caps."""
import itertools
import json
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h9_vertex_walk_20261013_caps import finite_vertex
from scripts.fresh_h9_replaced_core_20261013_cores import affine_census
from scripts.fresh_h9_replaced_core_20261013_geometry import volume_gate

PREFIX="data/enumerations/fresh_h9_vertex_walk_20261013"

def run():
    G,I,h,old,oldcase=selected_coupled_case()
    bank=read(PREFIX+"_caps.json");pool=angle_pool(ANGLE_SOURCE)
    new=[cap_profile(r,f"walk{i}")for i,r in enumerate(bank["caps"])]
    graph=[set()for _ in new];pairs={}
    for i,j in itertools.combinations(range(len(new)),2):
        r=relation(new[i],new[j],pool);pairs[i,j]=r
        if r["kind"]in("finite","dotted"):
            graph[i].add(j);graph[j].add(i)
    cliques=[]
    def bk(R,P,X):
        if not P and not X:
            cliques.append(tuple(sorted(R)));return
        u=max(P|X,key=lambda u:len(P&graph[u]))if P or X else None
        for v in sorted(P-(graph[u]if u is not None else set())):
            bk(R|{v},P&graph[v],X&graph[v]);P.remove(v);X.add(v)
    bk(set(),set(range(len(new))),set())
    cliques.sort(key=lambda x:(-len(x),x))
    print({"maximal_cliques":len(cliques),"sizes":[len(c)for c in cliques]},flush=True)
    cases=[]
    for ids in cliques:
        dropped=sorted(set().union(*(set(bank["caps"][i]["dropped"])for i in ids)))
        selected={n:r for n,r in old.items()if n not in dropped}
        selected.update({new[i]["name"]:new[i]for i in ids})
        names=list(selected)
        checks={k:r for k,r in oldcase["pair_checks"].items()if all(n in selected for n in k.split("|"))}
        for i in ids:
            checks.update({"|".join(sorted((new[i]["name"],n))):r for n,r in bank["caps"][i]["pair_checks"].items()if n in selected})
        for i,j in itertools.combinations(ids,2):
            checks["|".join(sorted((new[i]["name"],new[j]["name"])))]=pairs[i,j]
        assert len(checks)==len(names)*(len(names)-1)//2
        assert all(r["kind"]in("finite","dotted")for r in checks.values())
        census=affine_census(names,checks)
        vertices=[finite_vertex(new[i],bank["caps"][i],I,h,selected)for i in ids]
        case={"selected_caps":list(ids),"dropped":dropped,"walls":names,"facet_count":len(names),
              "pair_checks":checks,"affine_census":census,"actual_finite_vertices":vertices}
        if census["affine_circuits"]:
            case["status"]="REJECTED_POSITIVE_AFFINE_COMPACTNESS"
        else:
            case["geometry"]=volume_gate(I,h,list(selected.values()))
            case["status"]=case["geometry"]["status"]
        cases.append(case)
        result={"base_source":bank["base_source"],"base_case":bank["base_case"],"cap_source":PREFIX+"_caps.json",
                "new_pair_relations":{f"{i}|{j}":r for (i,j),r in pairs.items()},
                "maximal_accepted_cliques":[list(x)for x in cliques],"cases":cases,
                "scope":"Exact coherent finite-vertex cap replacements on the old73 chamber. Full pair/affine gates and saved strict escapes. This is not a complete search over all cap families or finite labels."}
        (ROOT/(PREFIX+"_geometry.json")).write_text(json.dumps(result,indent=2)+"\n")
        print({"case":len(cases)-1,"caps":ids,"facets":len(names),"affine":len(census["affine_circuits"]),
               "status":case["status"],"max_norm":case.get("geometry",{}).get("maximum_numeric_vertex_norm")},flush=True)
        if case["status"]=="UNRESOLVED_VOLUME":break
if __name__=="__main__":run()
