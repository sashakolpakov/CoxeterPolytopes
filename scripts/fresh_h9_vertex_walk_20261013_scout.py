#!/usr/bin/env python3
"""Targeted finite spherical extensions of rank-eight edges; numerical proposals."""
import itertools
import json
import time
import argparse
import numpy as np
import mpmath as mp
from scipy.linalg import null_space
from scripts.fresh_h9_gap_completion_20261013_base import *
from scripts.fresh_h9_vertex_walk_20261013_base import chamber

OUT=ROOT/"data/enumerations/fresh_h9_vertex_walk_20261013_scout.json"

def run(stage=0):
    start=time.monotonic();mp.mp.dps=90
    G,I,h,ps,case=chamber(stage)
    names=list(ps);n=len(names)
    C=np.array([[float(value(r["field"],x))for x in r["c"]]for r in ps.values()])
    Pn=np.array([[float(value(r["field"],x))for x in r["p"]]for r in ps.values()])
    q=np.array([float(value(K,r["q"]))for r in ps.values()])
    In=np.array([[float(value(K,x))for x in row]for row in I])
    hn=np.array([float(value(K,x))for x in h])
    yn=np.array([float(value(K,parse(x)))for x in case["geometry"]["escapes"][0]["certificate"]["y"]])
    rawgram=C@Pn.T
    eligible={i for i,r in enumerate(ps.values())if r["field"] is K and abs(q[i]-2)<1e-8}
    neighbors=[set()for _ in names]
    for i,j in itertools.combinations(eligible,2):
        if case["pair_checks"]["|".join(sorted((names[i],names[j])))]["kind"]=="finite":
            neighbors[i].add(j);neighbors[j].add(i)
    blocks=[]
    def walk(ids,choices):
        if len(ids)==8:
            blocks.append(ids);return
        if len(ids)+len(choices)<8:return
        for i in sorted(choices):
            child=ids+(i,)
            if np.linalg.eigvalsh(rawgram[np.ix_(child,child)])[0]<=1e-8:continue
            walk(child,{j for j in choices&neighbors[i]if j>i})
    walk((),eligible)
    print({"K_square_two_rank8_blocks":len(blocks),"walls":n},flush=True)
    target_values={3:1.,4:np.sqrt(2),5:(1+np.sqrt(5))/2}
    anchor_labels=[2,3,4,5,6,8,10,12]
    allowed=np.array([np.cos(np.pi/m)**2 for m in anchor_labels])
    patterns=[]
    for count in range(1,4):
        for support in itertools.combinations(range(8),count):
            for labels in itertools.product(target_values,repeat=count):
                vector=np.zeros(8)
                for j,m in zip(support,labels):vector[j]=-target_values[m]
                patterns.append((support,labels,vector))
    patterns.append(((),(),np.zeros(8)))
    mandatory=np.array([i for i,name in enumerate(names)if name.startswith("seed")])
    seen={};stats={"blocks":len(blocks),"elliptic_extensions":0,"anchor_solutions":0,
                   "cuts_target":0,"seed_compatible":0,"at_most_two_drops":0}
    for block_index,block in enumerate(blocks):
        old=rawgram[np.ix_(block,block)];old_inverse=np.linalg.inv(old)
        N=null_space(C[list(block)])
        assert N.shape==(10,2)
        Q=N.T@In@N
        outside=[j for j in range(n)if j not in block]
        for support,labels,vec in patterns:
            rem=2-vec@old_inverse@vec
            if rem<=1e-7:continue
            stats["elliptic_extensions"]+=1
            origin=Pn[list(block)].T@(old_inverse@vec)
            assert abs(origin@In@origin-(2-rem))<1e-7
            for anchor in outside:
                aa=C[anchor]@N
                if aa@aa<1e-14:continue
                direction=np.array([-aa[1],aa[0]])
                A=direction@Q@direction
                if abs(A)<1e-12:continue
                for am in anchor_labels:
                    target=-np.sqrt(2*q[anchor])*np.cos(np.pi/am)
                    rhs=target-C[anchor]@origin
                    t0=aa*(rhs/(aa@aa))
                    B=2*t0@Q@direction;D=B*B-4*A*(t0@Q@t0-rem)
                    if D<=1e-9:continue
                    for branch in(-1,1):
                        p=origin+N@(t0+direction*((-B+branch*np.sqrt(D))/(2*A)))
                        stats["anchor_solutions"]+=1
                        if p@In@hn<=1e-7 or p@In@yn>=-1e-7:continue
                        stats["cuts_target"]+=1
                        pair=C@p;ratio=pair*pair/(2*q)
                        good=(pair<=1e-7)&((ratio>1+1e-7)|(np.min(abs(ratio[:,None]-allowed),axis=1)<2e-7))
                        if not np.all(good[mandatory]):continue
                        stats["seed_compatible"]+=1
                        bad=np.flatnonzero(~good)
                        if len(bad)>2:continue
                        stats["at_most_two_drops"]+=1
                        key=tuple(np.round(p,7))
                        rec={"block":[names[i]for i in block],"block_support":list(support),
                             "block_labels":list(labels),"anchor":names[anchor],
                             "anchor_label":am,"branch":branch,"drop":[names[i]for i in bad],
                             "seed_pairings":p.tolist(),"target_pairing":float(p@In@yn),
                             "spherical_residual_norm":float(rem),
                             "anchor_field_is_K":ps[names[anchor]]["field"] is K}
                        prior=seen.get(key)
                        score=(len(bad),not rec["anchor_field_is_K"],am not in(2,3,5),4 in labels,abs(rec["target_pairing"]))
                        if prior is None or score<prior[0]:seen[key]=(score,rec)
        print({"block":block_index+1,"stats":stats,"distinct":len(seen),"seconds":round(time.monotonic()-start,1)},flush=True)
    rows=[r for _,r in sorted(seen.values(),key=lambda t:t[0])]
    source="data/enumerations/fresh_h9_gap_completion_20261013_coupled_geometry.json"if stage==0 else "data/enumerations/fresh_h9_vertex_walk_20261013"+("_geometry.json"if stage==1 else"_second_geometry.json")
    out={"status":"NUMERICAL_TARGETED_SPHERICAL_EDGE_EXTENSIONS","source":source,
         "source_case":0,"target":case["geometry"]["escapes"][0],"stats":stats,"proposals":rows,
         "scope":"Numerical proposal scan on K-defined rank-eight spherical blocks. New cap has an elliptic nine-wall extension with one to three finite attachments to that block and a finite angle to another retained wall. Proposals cut the specific old strict escape and require at most two optional-wall removals. No exact positive claim or exhaustive family exclusion."}
    output=OUT if stage==0 else OUT.with_name(OUT.stem+f"_stage{stage}.json")
    output.write_text(json.dumps(out,indent=2)+"\n")
    print({"status":out["status"],"proposals":len(rows),"stats":stats},flush=True)
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--stage",type=int,default=0)
    run(parser.parse_args().stage)
