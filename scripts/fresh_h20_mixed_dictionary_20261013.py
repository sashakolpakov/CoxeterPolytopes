#!/usr/bin/env python3
"""Bounded exact selection across two H20 patches and a new rotation mirror.

The dictionary combines old580, genuinely new mixed767 planes, and the
Q(sqrt2) rotation 3*r0+(2-sqrt2)*r74. Every conflict and escape-hitting
clause uses integer field arithmetic. UNSAT excludes only this dictionary.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import json
from math import gcd,ceil
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
import z3
import fresh_h20_mixed_b3_20261013 as mixed

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'data/enumerations/fresh_h20_mixed_dictionary_20261013.json'
OLD=ROOT/'data/enumerations/new_dimension_h20_six_patch_boundary_20261013.json'
ESCAPES=ROOT/'data/enumerations/h20_gluing_independent_all_deletions_20261013.json'
SELECTED=ROOT/'data/enumerations/h20_gluing_independent_escape_20261013.json'

def conflict(p,q,r):
    if mixed.sign(*p)>0:return True
    n=mixed.multiply(p,p);d=mixed.multiply(q,r)
    if mixed.sign(d[0]-n[0],d[1]-n[1])<=0:return False
    if any((a*n[0],a*n[1])==(b*d[0],b*d[1]) for a,b in [(1,0),(4,1),(2,1),(4,3)]):return False
    return (4*n[0],4*n[1])!=(2*d[0]+2*d[1],d[0]+2*d[1])

def dictionary():
    Q,W,G=mixed.data();zero=np.zeros(22,dtype=object)
    old=np.array(json.loads(OLD.read_text())['six_tile']['primitive_outward_normals'],dtype=object)
    newer=json.loads(mixed.OUTPUT.read_text())['eight_tile']['boundary_normal_pairs']
    roots=[];sources=[];index={}
    def add(r,source):
        r=mixed.canonical(*r)
        if r not in index:index[r]=len(roots);roots.append(r);sources.append([])
        sources[index[r]].append(source)
    for i,r in enumerate(old):add((r,zero),['old580',i])
    for i,r in enumerate(newer):add(r,['mixed767',i])
    add((3*old[0]+2*old[74],-old[74]),['targeted_rotation',0])
    R=np.array(roots,dtype=object);A,B=R[:,0],R[:,1]
    Ga=A@Q@A.T+2*B@Q@B.T;Gb=A@Q@B.T+B@Q@A.T
    norms=[(int(Ga[i,i]),int(Gb[i,i])) for i in range(len(R))]
    edges=[];neighbors=[set() for r in R]
    for i in range(len(R)):
        assert mixed.sign(*norms[i])>0
        for j in range(i):
            if conflict((int(Ga[i,j]),int(Gb[i,j])),norms[i],norms[j]):edges.append([i,j]);neighbors[i].add(j);neighbors[j].add(i)
    center=np.array(json.loads(SELECTED.read_text())['common_strict_timelike_center'],dtype=object)
    ca=A@Q@center;cb=B@Q@center
    assert center@Q@center<0 and center@Q@W[0]==0
    assert all(mixed.sign(a,b)<0 for a,b in zip(ca,cb))
    return Q,W[0],R,sources,edges,neighbors,center

def hitting(v,Q,f,R):
    v=np.array(v,dtype=object);assert v@Q@f==0 and v@Q@v>0
    a=R[:,0]@Q@v;b=R[:,1]@Q@v
    return [i for i in range(len(R)) if mixed.sign(a[i],b[i])>0]

def escape(ids,Q,f,R,center,seed):
    normals=R[ids];Z=null_space(np.array([f@Q,center@Q],dtype=float))
    real=np.array(normals[:,0],dtype=float)+np.sqrt(2)*np.array(normals[:,1],dtype=float)
    A=real@np.array(Q,dtype=float)@Z;b=-real@np.array(Q@center,dtype=float)
    scale=np.linalg.norm(A,axis=1);A=A/scale[:,None];b=b/scale
    rng=np.random.default_rng(seed);directions=list(np.eye(20))+list(-np.eye(20))+list(rng.normal(size=(40,20)))
    c0=normals[:,0]@Q@center;c1=normals[:,1]@Q@center
    for objective,direction in enumerate(directions):
        prop=linprog(-direction,A_ub=A,b_ub=b,bounds=(None,None),method='highs')
        if prop.status==3:prop=linprog(-direction,A_ub=A,b_ub=b,bounds=(-10000,10000),method='highs')
        if not prop.success:continue
        point=np.array(center,dtype=float)+Z@prop.x
        raw=np.array([int(round(float(x)*10**7)) for x in point],dtype=object)
        v=2*raw-(raw@Q@f)*f
        v0=normals[:,0]@Q@v;v1=normals[:,1]@Q@v
        shift=max(0,ceil(max((float(a)+2**.5*float(b))/(-float(c)-2**.5*float(d)) for a,b,c,d in zip(v0,v1,c0,c1))))+2
        while any(mixed.sign(a+shift*c,b+shift*d)>=0 for a,b,c,d in zip(v0,v1,c0,c1)):shift+=1
        v=v+shift*center;v=v//gcd(*map(int,v))
        if v@Q@v<=0:continue
        assert v@Q@f==0
        assert all(mixed.sign(a,b)<0 for a,b in zip(normals[:,0]@Q@v,normals[:,1]@Q@v))
        return {'vector':list(map(int,v)),'square':int(v@Q@v),'objective':objective}
    return None

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--rounds',type=int,default=20);parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    Q,f,R,sources,edges,neighbors,center=dictionary();n=len(R)
    print('Dictionary',n,'planes,',len(edges),'conflicts,',sum(any(r[1]) for r in R),'nonrational',flush=True)
    z3.set_param(proof=True)
    if args.verify:
        out=json.loads(OUTPUT.read_text());assert out['normal_pairs']==R.tolist()
        for row in out['hitting_clauses']:assert hitting(row['vector'],Q,f,R)==row['positive_wall_indices']
        if out['status']=='EXACT_UNSAT_DICTIONARY':
            core=out['unsat_core'];s=z3.Solver();x=[z3.Bool(f'x_{i}') for i in range(n)]
            for i,j in core['conflict_edges']:assert j in neighbors[i];s.add(z3.Or(z3.Not(x[i]),z3.Not(x[j])))
            for k in core['hitting_clause_indices']:s.add(z3.Or([x[i] for i in out['hitting_clauses'][k]['positive_wall_indices']]))
            assert s.check()==z3.unsat
        for row in out['selected_candidates']:
            ids=row['indices'];assert not any(j in neighbors[i] for k,i in enumerate(ids) for j in ids[:k])
            if row['escape'] is not None:assert not set(ids)&set(hitting(row['escape']['vector'],Q,f,R))
        print('PASS_EXACT_DICTIONARY_REPLAY',out['status'],flush=True);return
    out={'status':'BOUNDED_SEARCH_IN_PROGRESS','normal_pairs':R.tolist(),'normal_sources':sources,
         'common_strict_timelike_center':list(map(int,center)),'complete_conflict_count':len(edges),
         'hitting_clauses':[],'selected_candidates':[],
         'scope':'Only subsets of the combined old580 + mixed767 + one targetedrotation dictionary; no global H20 impossibility claim.'}
    solver=z3.Solver();x=[z3.Bool(f'x_{i}') for i in range(n)]
    solver.add(*[z3.Or(z3.Not(x[i]),z3.Not(x[j])) for i,j in edges])
    seen=set()
    def add_escape(v,origin):
        hit=hitting(v,Q,f,R);key=tuple(hit)
        if key in seen:return
        seen.add(key);k=len(out['hitting_clauses'])
        out['hitting_clauses'].append({'vector':list(map(int,v)),'origin':origin,'positive_wall_indices':hit})
        solver.assert_and_track(z3.Or([x[i] for i in hit]),z3.Bool(f'h_{k}'))
    for key,row in json.loads(ESCAPES.read_text())['mandatory'].items():add_escape(row['vector'],f'old580_delete_{key}')
    for row in json.loads(SELECTED.read_text())['results']:add_escape(row['witness']['vector'],row['candidate_id'])
    print('Seed exact escape clauses',len(out['hitting_clauses']),flush=True)
    for iteration in range(args.rounds+1):
        answer=solver.check();print('Round',iteration,str(answer),flush=True)
        if answer==z3.unsat:
            core=sorted(int(str(t)[2:]) for t in solver.unsat_core())
            variables=set(i for k in core for i in out['hitting_clauses'][k]['positive_wall_indices'])
            relevant=[[i,j] for i,j in edges if i in variables and j in variables]
            small=z3.Solver()
            for i,j in relevant:small.add(z3.Or(z3.Not(x[i]),z3.Not(x[j])))
            for k in core:small.add(z3.Or([x[i] for i in out['hitting_clauses'][k]['positive_wall_indices']]))
            assert small.check()==z3.unsat
            out['status']='EXACT_UNSAT_DICTIONARY';out['unsat_core']={'hitting_clause_indices':core,'conflict_edges':relevant,'variable_indices':sorted(variables),'smt2':small.to_smt2(),'z3_proof':str(small.proof())}
            break
        assert answer==z3.sat
        if iteration==args.rounds:out['status']='BOUNDED_ROUND_LIMIT_REACHED';break
        model=solver.model();ids=[i for i in range(n) if z3.is_true(model.eval(x[i]))]
        available=set(range(n))-set(ids)
        for i in ids:available-=neighbors[i]
        while available:
            i=min(available,key=lambda i:(len(neighbors[i]&available),i));ids.append(i);available.remove(i);available-=neighbors[i]
        ids.sort();w=escape(ids,Q,f,R,center,700+iteration)
        out['selected_candidates'].append({'iteration':iteration,'indices':ids,'walls':len(ids),'escape':w})
        print('Candidate',len(ids),'walls','EXACT_ESCAPE' if w else 'NO_BOUNDED_ESCAPE',flush=True)
        if w is None:out['status']='CANDIDATE_REQUIRES_FINITE_VOLUME_GATE';break
        add_escape(w['vector'],f'new_iteration_{iteration}')
        OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'clauses',len(out['hitting_clauses']),'candidates',len(out['selected_candidates']),flush=True)

if __name__=='__main__':main()
