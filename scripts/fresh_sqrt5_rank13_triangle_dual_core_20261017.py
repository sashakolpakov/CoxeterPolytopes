#!/usr/bin/env python3
"""Recover a stronger exact dual on the saved280-triangle support only."""
import argparse,hashlib,json
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_triangle_dual_core_20261017'
INPUT='fresh_sqrt5_rank13_joint_corner_all_nonade_repair_20261017'
def sparse_dual_bound(d,alphas):
    loads=defaultdict(F)
    for i,a in alphas.items():
        assert a>=0
        if a:
            for j,c in d['faces'][i]['coefficients']:
                assert c>0;loads[j]+=a*c
    beta={};maximizers={}
    for j in sorted(loads):
        var=d['variables'][j];assert var['index']==j and var['multiplicity']>0
        vertex=var['vertex'];value=loads[j]/var['multiplicity'];assert value>0
        if value>beta.get(vertex,F(0)):
            beta[vertex]=value;maximizers[vertex]=[j]
        elif value==beta[vertex]:maximizers[vertex].append(j)
    beta=dict(sorted(beta.items()));maximizers={i:maximizers[i]for i in beta}
    numerator=sum((d['faces'][i]['required']*a for i,a in alphas.items()),F(0))
    denominator=sum(beta.values(),F(0));assert denominator>0
    for j,load in loads.items():
        var=d['variables'][j];assert load<=beta[var['vertex']]*var['multiplicity']
    # Every omitted column has zero load, and every beta and multiplicity is nonnegative.
    return numerator/denominator,beta,maximizers,numerator,denominator

def domain():
    original=B.read(INPUT)
    cert=next(c for c in original['certificates']if c['kind']=='EXACT_RESTRICTED_MODEL_DUAL')
    B.verify_certificate(original,cert)
    indices=[i for i,a in cert['face_multipliers']]
    vertex_ids=sorted({original['variables'][j]['vertex']for i in indices for j,a in original['faces'][i]['coefficients']})
    variables=[];vertices=[];lookup={}
    for old_vertex in vertex_ids:
        v=original['vertices'][old_vertex];vid=len(vertices);vv=[]
        for old,m in v['variables']:
            j=len(variables);lookup[old]=j;var=dict(original['variables'][old]);var.update(index=j,vertex=vid,original_variable=old)
            variables.append(var);vv.append([j,m])
        assert sum(m for j,m in vv)==78
        vertices.append(dict(index=vid,forest=v['forest'],variables=vv,original_vertex=old_vertex))
    faces=[]
    for old in indices:
        face=original['faces'][old];assert face['required']==2
        provenance=face['provenance'][0];assert 'triangle' in provenance['kind']
        faces.append(dict(index=len(faces),required=2,coefficients=[[lookup[j],m]for j,m in face['coefficients']],
            provenance=[provenance],original_face=old))
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'results/{INPUT}.json']
    return dict(sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},variables=variables,vertices=vertices,faces=faces,
        counts=dict(vertices=len(vertices),variables=len(variables),distinct_face_constraints=len(faces)))

def recover_unminimized(d,dual):
    candidates=[]
    for denominator in(10**6,10**8,10**10):
        a={i:F(round(-float(x)*denominator),denominator)for i,x in enumerate(dual)if x<0}
        a={i:x for i,x in a.items()if x>0}
        candidates.append((sparse_dual_bound(d,a)[0],a))
    bound,a=max(candidates,key=lambda x:x[0]);bound,beta,tight,num,den=sparse_dual_bound(d,a)
    return dict(kind='EXACT_RESTRICTED_MODEL_DUAL',bound=str(bound),objective_numerator=str(num),normalizer=str(den),
        face_multipliers=[[i,str(x)]for i,x in sorted(a.items())],vertex_multipliers=[[i,str(x)]for i,x in sorted(beta.items())],
        tight_orbits={str(i):v for i,v in tight.items()},excludes_target12=bound>12)

def run(verify):
    d=domain();out=ROOT/f'results/{TAG}.json';B.dual_bound=sparse_dual_bound
    if verify:
        saved=json.loads(out.read_text())
        for k,v in d.items():assert saved[k]==v,k
        for c in saved['certificates']:B.verify_certificate(d,c)
        return saved
    B.recover_dual=recover_unminimized
    solver,certificates=B.optimize(d)
    d.update(status='PASS_EXACT_REOPTIMIZED_TRIANGLE_DUAL_CORE',solver=solver,certificates=certificates,
        scope='This reoptimizes only the280 retained triangle inequalities and their74 complete vertex budgets; it does not repeat the111190-variable full model and makes no dimension claim.')
    out.write_text(json.dumps(d,indent=2)+'\n');return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.verify)
    print(d['status'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound')))for c in d['certificates']],flush=True)
