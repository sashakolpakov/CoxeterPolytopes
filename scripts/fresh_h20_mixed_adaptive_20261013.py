#!/usr/bin/env python3
"""Adaptive Q(sqrt2) mirrors guided by exact escape/SAT obstructions.

Each new mirror rotates one rational blocking pair to an allowed Coxeter
angle, preserves a common strict timelike center, and cuts a witnessed
spacelike escape. It need not support the entire original patch.
This distinction permits actual trimming as well as outward enlargement.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import combinations
import argparse
import json
from math import isqrt,lcm
from pathlib import Path
import numpy as np
import z3
import fresh_h20_mixed_b3_20261013 as mixed
import fresh_h20_mixed_dictionary_20261013 as selection

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'data/enumerations/fresh_h20_mixed_adaptive_20261013.json'

def sqrt_field(x):
    assert x>0
    for factor in (1,2):
        y=x/factor;a,b=isqrt(y.numerator),isqrt(y.denominator)
        if a*a==y.numerator and b*b==y.denominator:return (F(a,b),F(0)) if factor==1 else (F(0),F(a,b))
    return None

def rotations(i,j,R,Q,center,witnesses):
    r,s=R[i],R[j]
    if any(r[1]) or any(s[1]):return []
    rr=F(int(r[0]@Q@r[0]));ss=F(int(s[0]@Q@s[0]));ip=F(int(r[0]@Q@s[0]));delta=rr-ip*ip/ss
    if delta<=0:return []
    rows=[]
    for m,factor in [(2,F(0)),(3,F(1,3)),(4,F(1)),(6,F(3)),(8,F(1))]:
        cotroot=(F(0),F(0)) if m==2 else sqrt_field(factor*delta/ss)
        if cotroot is None:continue
        if m==8:cotroot=(cotroot[0]+2*cotroot[1],cotroot[0]+cotroot[1])
        t=(-ip/ss-cotroot[0],-cotroot[1]);d=lcm(t[0].denominator,t[1].denominator)
        v=mixed.canonical(d*r[0]+int(t[0]*d)*s[0],int(t[1]*d)*s[0])
        if mixed.sign(*mixed.pairing(v,(center,np.zeros(22,dtype=object)),Q))>=0:continue
        if not any(mixed.sign(*mixed.pairing(v,(np.array(w,dtype=object),np.zeros(22,dtype=object)),Q))>0 for w in witnesses):continue
        assert not selection.conflict(mixed.pairing(v,s,Q),mixed.pairing(v,v,Q),mixed.pairing(s,s,Q))
        rows.append((v,{'parent_indices':[i,j],'angle_m':m,'coefficient_t':[str(x) for x in t],
                        'supports_original_intersection_by_positive_combination':mixed.sign(int(t[0]*d),int(t[1]*d))>=0}))
    return rows

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    Q,f,R,sources,edges,neighbors,center=selection.dictionary();base_count=len(R)
    if args.verify:
        out=json.loads(OUTPUT.read_text());archived=np.array(out['normal_pairs'],dtype=object)
        assert np.array_equal(archived[:base_count],R)
        assert out['common_strict_timelike_center']==list(map(int,center))
        for row in out['new_mirrors']:
            index=row['index'];assert index==len(R)
            i,j=row['parent_indices'];assert i<index and j<index
            a,b=map(F,row['coefficient_t']);d=lcm(a.denominator,b.denominator)
            assert not any(R[i,1]) and not any(R[j,1])
            v=mixed.canonical(d*R[i,0]+int(a*d)*R[j,0],int(b*d)*R[j,0])
            assert v==tuple(tuple(z) for z in row['normal_pair'])==tuple(tuple(z) for z in archived[index])
            assert mixed.sign(*mixed.pairing(v,(center,np.zeros(22,dtype=object)),Q))<0
            assert not selection.conflict(mixed.pairing(v,R[j],Q),mixed.pairing(v,v,Q),mixed.pairing(R[j],R[j],Q))
            R=np.concatenate((R,np.array([v],dtype=object)),axis=0)
        assert np.array_equal(R,archived)
        for row in out['hitting_clauses']:assert selection.hitting(row['vector'],Q,f,R)==row['positive_wall_indices']
        last=out['iterations'][-1];assert last['SAT']=='unsat' and not last['added_mirror_indices']
        solver=z3.Solver();x=[z3.Bool(f'x_{i}') for i in range(len(R))]
        for i,j in last['core_conflict_edges']:
            assert selection.conflict(mixed.pairing(R[i],R[j],Q),mixed.pairing(R[i],R[i],Q),mixed.pairing(R[j],R[j],Q))
            solver.add(z3.Or(z3.Not(x[i]),z3.Not(x[j])))
        for k in last['core_clause_indices']:solver.add(z3.Or([x[i] for i in out['hitting_clauses'][k]['positive_wall_indices']]))
        assert solver.check()==z3.unsat
        print('PASS_EXACT_ADAPTIVE_ROTATIONS_AND_FINAL_DICTIONARY_UNSAT',len(R),'planes',len(out['new_mirrors']),'new mirrors',flush=True);return
    seed=json.loads(selection.OUTPUT.read_text());clauses=seed['hitting_clauses']
    out={'status':'IN_PROGRESS','base_dictionary':str(selection.OUTPUT.relative_to(ROOT)),
         'common_strict_timelike_center':list(map(int,center)),
         'new_mirrors':[],'iterations':[],'selected_candidates':[],
         'limits':{'rounds':20,'new_mirrors':120},
         'scope':'Bounded adaptive angle rotations of rational blocking pairs over Q(sqrt2), allowing trimming while preserving the displayed strict timelike center.'}
    lookup={tuple(tuple(x) for x in r):i for i,r in enumerate(R)}
    z3.set_param(proof=True)
    for iteration in range(20):
        solver=z3.Solver();x=[z3.Bool(f'x_{i}') for i in range(len(R))]
        solver.add(*[z3.Or(z3.Not(x[i]),z3.Not(x[j])) for i,j in edges])
        for k,row in enumerate(clauses):solver.assert_and_track(z3.Or([x[i] for i in row['positive_wall_indices']]),z3.Bool(f'h_{k}'))
        answer=solver.check();print('Round',iteration,'planes',len(R),str(answer),flush=True)
        if answer==z3.sat:
            model=solver.model();ids=[i for i in range(len(R)) if z3.is_true(model.eval(x[i]))]
            available=set(range(len(R)))-set(ids)
            for i in ids:available-=neighbors[i]
            while available:
                i=min(available,key=lambda i:(len(neighbors[i]&available),i));ids.append(i);available.remove(i);available-=neighbors[i]
            ids.sort();w=selection.escape(ids,Q,f,R,center,1800+iteration)
            out['selected_candidates'].append({'iteration':iteration,'indices':ids,'escape':w})
            print('Selected',len(ids),'walls', 'exact escape' if w else 'NO ESCAPE FOUND',flush=True)
            if w is None:out['status']='CANDIDATE_REQUIRES_FINITE_VOLUME_GATE';break
            hit=selection.hitting(w['vector'],Q,f,R)
            clauses.append({'vector':w['vector'],'origin':f'adaptive_iteration_{iteration}','positive_wall_indices':hit})
            out['iterations'].append({'round':iteration,'SAT':'sat','new_escape_clause':len(clauses)-1})
        else:
            assert answer==z3.unsat
            core=sorted(int(str(t)[2:]) for t in solver.unsat_core());variables=set(i for k in core for i in clauses[k]['positive_wall_indices'])
            conflicts=[(i,j) for i,j in edges if i in variables and j in variables]
            witnesses=[clauses[k]['vector'] for k in core];added=[]
            for i,j in conflicts:
                for first,second in ((i,j),(j,i)):
                    for v,metadata in rotations(first,second,R,Q,center,witnesses):
                        if v in lookup:continue
                        if len(R)-base_count>=120:break
                        index=len(R);lookup[v]=index
                        old_R=R;R=np.concatenate((R,np.array([v],dtype=object)),axis=0)
                        norms=mixed.pairing(v,v,Q);neighbors.append(set())
                        for k,r in enumerate(old_R):
                            if selection.conflict(mixed.pairing(v,r,Q),norms,mixed.pairing(r,r,Q)):
                                edges.append([index,k]);neighbors[index].add(k);neighbors[k].add(index)
                        for row in clauses:
                            if mixed.sign(*mixed.pairing(v,(np.array(row['vector'],dtype=object),np.zeros(22,dtype=object)),Q))>0:row['positive_wall_indices'].append(index)
                        record={'index':index,'normal_pair':[list(z) for z in v],'creation_round':iteration,**metadata}
                        out['new_mirrors'].append(record);added.append(index)
                    if len(added)>=12:break
                if len(added)>=12:break
            out['iterations'].append({'round':iteration,'SAT':'unsat','core_clause_indices':core,'core_conflict_edges':conflicts,'added_mirror_indices':added})
            print('Core',len(core),'clauses',len(conflicts),'edges; added',len(added),'mirrors',flush=True)
            if not added:
                out['status']='ADAPTIVE_NEW_MIRROR_LIMIT_REACHED' if len(R)-base_count>=120 else 'NO_NEW_FIELD_ROTATION_FROM_CURRENT_CORE';break
        out['normal_pairs']=R.tolist();out['hitting_clauses']=clauses
        OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    else:out['status']='BOUNDED_ADAPTIVE_ROUND_LIMIT'
    out['normal_pairs']=R.tolist();out['hitting_clauses']=clauses
    out['final_plane_count']=len(R);out['new_plane_count']=len(R)-base_count
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],'new planes',len(R)-base_count,flush=True)

if __name__=='__main__':main()
