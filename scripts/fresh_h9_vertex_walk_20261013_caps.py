#!/usr/bin/env python3
"""Exact finite-vertex caps targeted at the old73 strict escape."""
import json
import time
import mpmath as mp
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K
from scripts.fresh_h4_blocks_independent_verify import inverse

PREFIX="data/enumerations/fresh_h9_vertex_walk_20261013"

def solve_line(E,rows,rhs):
    a=[list(row)+[z]for row,z in zip(rows,rhs,strict=True)]
    pivots=[]
    for j in range(10):
        pivot=next((i for i in range(len(pivots),len(a))if E.sign(a[i][j])),None)
        if pivot is None:continue
        i=len(pivots);a[i],a[pivot]=a[pivot],a[i]
        p=a[i][j];a[i]=[E.div(x,p)for x in a[i]]
        for r in range(len(a)):
            if r==i:continue
            t=a[r][j]
            a[r]=[E.sub(x,E.mul(t,y))for x,y in zip(a[r],a[i],strict=True)]
        pivots.append(j)
    assert len(pivots)==9
    free=next(i for i in range(10)if i not in pivots)
    c=[E.zero]*10;v=[E.zero]*10;v[free]=E.one
    for i,j in enumerate(pivots):
        c[j]=a[i][-1];v[j]=E.neg(a[i][free])
    return free,tuple(c),tuple(v)

def make_cap(proposal,I,h,ps,target,pool):
    anchor=ps[proposal["anchor"]];E=anchor["field"]
    if E is K:
        E=Quadratic(K,Z)
        anchor=dict(anchor,field=E,c=tuple(E.lift(t)for t in anchor["c"]),
                    p=tuple(E.lift(t)for t in anchor["p"]))
    if E is not K and E.base is not K and E.base.base is K and E.base.d in(Z,O):
        collapse=(lambda x:x[0])if E.base.d==Z else(lambda x:add(x[0],x[1]))
        F0=E;E=Quadratic(K,collapse(F0.d))
        anchor=dict(anchor,field=E,c=tuple(tuple(collapse(t)for t in x)for x in anchor["c"]),
                    p=tuple(tuple(collapse(t)for t in x)for x in anchor["p"]))
    assert E is not K and E.base is K and anchor["q"]==scale(O,2)
    block=proposal["block"]
    rows=[tuple(E.lift(x)for x in ps[n]["c"])for n in block]+[anchor["c"]]
    values=[Z]*8
    for j,m in zip(proposal["block_support"],proposal["block_labels"],strict=True):
        values[j]=K.neg({3:O,5:P}[m])
    rhs=[E.lift(x)for x in values]+[E.lift(K.neg({2:Z,3:O,5:P}[proposal["anchor_label"]]))]
    free,u,v=solve_line(E,rows,rhs)
    IE=tuple(tuple(E.lift(x)for x in row)for row in I)
    A=pairing(E,v,IE,v);B=E.mul(E.rational(2),pairing(E,u,IE,v))
    C=E.sub(pairing(E,u,IE,u),E.rational(2))
    D=E.sub(E.mul(B,B),E.mul(E.rational(4),E.mul(A,C)))
    assert E.sign(A)and E.sign(D)>0
    F=Quadratic(E,D)
    candidates=[]
    for branch in(-1,1):
        t=(E.div(E.neg(B),E.mul(E.rational(2),A)),E.div(E.rational(branch),E.mul(E.rational(2),A)))
        p=tuple(F.add(F.lift(x),F.mul(t,F.lift(y)))for x,y in zip(u,v,strict=True))
        distance=sum((float(value(F,x))-y)**2 for x,y in zip(p,proposal["seed_pairings"],strict=True))
        candidates.append((distance,branch,t,p))
    distance,branch,t,p=min(candidates,key=lambda x:x[0])
    assert distance<1e-9
    cf=matvec(F,tuple(tuple(lift_K(F,x)for x in row)for row in I),p)
    row={"name":"walk","field":F,"c":cf,"p":p,"q":scale(O,2)}
    assert F.sign(F.sub(dot(F,cf,p),F.rational(2)))==0
    assert all(F.sign(x)<=0 for x in p)
    assert F.sign(dot(F,cf,tuple(lift_K(F,x)for x in h)))>0
    cut=dot(F,cf,tuple(lift_K(F,x)for x in target));assert F.sign(cut)<0
    for rr,want in zip(rows,rhs,strict=True):
        assert F.sign(F.sub(dot(F,p,tuple(F.lift(x)for x in rr)),F.lift(want)))==0
    checks={n:relation(row,r,pool)for n,r in ps.items()}
    dropped=[n for n,r in checks.items()if r["kind"]not in("finite","dotted")]
    assert not any(n.startswith("seed")for n in dropped)
    blockgram=[[dot(K,ps[a]["c"],ps[b]["p"])for b in block]for a in block]
    fullgram=[r+[x]for r,x in zip(blockgram,values,strict=True)]+[values+[scale(O,2)]]
    assert inertia_K(fullgram)==(9,0,0)
    rec={"base_radicand":enc(E.d),"discriminant":E.encode(D),
         "block_names":block,"block_support":proposal["block_support"],"block_labels":proposal["block_labels"],
         "anchor_name":proposal["anchor"],"anchor_label":proposal["anchor_label"],
         "free_coordinate":free,"line_constant":list(map(E.encode,u)),"line_direction":list(map(E.encode,v)),
         "norm_polynomial":list(map(E.encode,(A,B,C))),"parameter":F.encode(t),"branch":branch,
         "seed_pairings":list(map(F.encode,p)),"seed_coefficients":list(map(F.encode,cf)),
         "target_pairing":F.encode(cut),"dropped":dropped,"pair_checks":checks,
         "spherical_vertex_gram":[[enc(x)for x in r]for r in fullgram]}
    return row,rec

def finite_vertex(row,record,I,h,selected):
    F=row["field"];block=record["block_names"];incident=[selected[n]for n in block]+[row]
    gram=[[parse(x)for x in r]for r in record["spherical_vertex_gram"]]
    roots=[]
    for r in incident:
        common,er=embed(F,r["field"]);assert key(common)==key(F)
        roots.append((tuple(er(x)for x in r["p"]),tuple(er(x)for x in r["c"])))
    gh=tuple(dot(F,c,tuple(lift_K(F,x)for x in h))for p,c in roots)
    coeff=matvec(F,tuple(tuple(lift_K(F,x)for x in r)for r in inverse(gram)),gh)
    y=tuple(lift_K(F,x)for x in h)
    for (p,c),a in zip(roots,coeff,strict=True):
        y=tuple(F.sub(x,F.mul(a,z))for x,z in zip(y,p,strict=True))
    IF=tuple(tuple(lift_K(F,x)for x in r)for r in I)
    q=pairing(F,y,IF,y);cp=pairing(F,y,IF,tuple(lift_K(F,x)for x in h))
    assert F.sign(q)<0 and F.sign(cp)<0
    signs={}
    for n,r in selected.items():
        common,er=embed(F,r["field"])
        signs[n]=common.sign(dot(common,tuple(er(x)for x in r["c"]),tuple(lift_existing(common,F,x)for x in y)))
    zeros=set(block+[row["name"]])
    assert all(s==(0 if n in zeros else 1)for n,s in signs.items())
    return {"y":list(map(F.encode,y)),"norm":F.encode(q),"center_pairing":F.encode(cp),
            "all_wall_signs":signs,"zero_walls":sorted(zeros),
            "gram":record["spherical_vertex_gram"],"point_field_from_cap":row["name"]}

def run():
    mp.mp.dps=90;start=time.monotonic()
    G,I,h,ps,case=selected_coupled_case()
    target=tuple(parse(x)for x in case["geometry"]["escapes"][0]["certificate"]["y"])
    pool=angle_pool(ANGLE_SOURCE);scout=read(PREFIX+"_scout.json")
    chosen=[15,14,13,12,11,9,8,5,4,2,0]
    caps=[];rejected=[]
    for i in chosen:
        proposal=scout["proposals"][i]
        row,rec=make_cap(proposal,I,h,ps,target,pool)
        rec["proposal_index"]=i
        if len(rec["dropped"])>2:
            rejected.append({"proposal_index":i,"drop":rec["dropped"]});continue
        name=f"walk{len(caps)}";row["name"]=name
        selected={n:r for n,r in ps.items()if n not in rec["dropped"]}|{name:row}
        rec["finite_vertex"]=finite_vertex(row,rec,I,h,selected)
        caps.append(rec)
        result={"status":"EXACT_TARGETED_SPHERICAL_VERTEX_CAPS","base_source":"data/enumerations/fresh_h9_gap_completion_20261013_coupled_geometry.json",
                "base_case":0,"proposal_source":PREFIX+"_scout.json","angle_source":ANGLE_SOURCE,
                "target":case["geometry"]["escapes"][0],"caps":caps,"rejected":rejected,
                "scope":"Exact accepted cap equations, selected finite angles, target cut, all old pair gates after explicit optional-wall removal, and an actual finite vertex. No full-volume or exhaustive family claim."}
        (ROOT/(PREFIX+"_caps.json")).write_text(json.dumps(result,indent=2)+"\n")
        print({"cap":len(caps)-1,"proposal":i,"dropped":rec["dropped"],"actual_vertex":True,
               "seconds":round(time.monotonic()-start,2)},flush=True)
if __name__=="__main__":run()
