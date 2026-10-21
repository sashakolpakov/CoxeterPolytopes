#!/usr/bin/env python3
"""Exact whole-orbit A1 padding of a triangle-only rank13 dual into rank14."""
import argparse, hashlib, json
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as C

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank14_triangle_padding_barrier_20261017'
INPUT='fresh_sqrt5_rank13_active_triangle_dual_20261017'

def thresholds(n):
    paired=[F(comb(n-i,2) if n-i>=2 else 0)+F(comb(i,2) if i>=2 else 0) for i in range(n//2+1)]
    ratios=[x/2 for x in paired]
    minimum=min(ratios)
    assert ratios[0]>minimum
    return dict(dimension=n,paired_h_ratios=list(map(str,ratios)),strict_f2_over_f0_lower=str(minimum),
        triangle_quad_payment_target=str(5*minimum-comb(n,2)))

def run():
    old=B.read(INPUT);B.dual_bound=C.sparse_dual_bound
    cert=next(c for c in old['certificates'] if c['kind']=='EXACT_RESTRICTED_MODEL_DUAL')
    B.verify_certificate(old,cert);assert F(cert['bound'])>14
    graphs=dict(B.C.models());variables=[];vertices=[];mapping={};checks=[]
    for v in old['vertices']:
        g={}
        for name in v['forest']:
            h=graphs[name];offset=len(g)
            for i,row in h.items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==13;pieces=B.C.pieces(g)
        old_codes={B.dump(var['marked_forest_code']):var for var in old['variables'] if var['vertex']==v['index']}
        old_pair_codes={pair:B.dump(sorted(B.C.code(p,pair) for p in pieces)) for pair in combinations(range(13),2)}
        g[13]={};pieces=B.C.pieces(g);orbits=defaultdict(list)
        for pair in combinations(range(14),2):orbits[B.dump(sorted(B.C.code(p,pair) for p in pieces))].append(pair)
        vv=[];new_lookup={};forest=sorted(v['forest']+['A1'])
        for code,pairs in sorted(orbits.items()):
            j=len(variables);new_lookup[code]=j;vv.append([j,len(pairs)])
            variables.append(dict(index=j,vertex=v['index'],forest=forest,marked_forest_code=json.loads(code),multiplicity=len(pairs)))
        assert sum(m for j,m in vv)==91
        vertices.append(dict(index=v['index'],forest=forest,variables=vv))
        for pair,old_code in old_pair_codes.items():
            previous=old_codes[old_code]
            code=B.dump(sorted(B.C.code(p,pair) for p in pieces));j=new_lookup[code]
            i=previous['index']
            assert i not in mapping or mapping[i]==j
            mapping[i]=j
            assert variables[j]['multiplicity']>=previous['multiplicity']
        assert len(set(mapping[i] for i,m in v['variables']))==len(v['variables'])
        checks.append(dict(vertex=v['index'],old_forest=v['forest'],padded_forest=forest,
            old_total_pairs=78,new_total_pairs=91,orbit_map=[[i,mapping[i],m,variables[mapping[i]]['multiplicity']] for i,m in v['variables']]))
    assert len(mapping)==len(old['variables']) and len(set(mapping.values()))==len(mapping)
    assert len({tuple(v['forest']) for v in vertices})==len(vertices)
    faces=[dict(index=f['index'],required=2,coefficients=[[mapping[i],m] for i,m in f['coefficients']],
                original_face=f['original_face'],provenance=f['provenance']) for f in old['faces']]
    d=dict(variables=variables,vertices=vertices,faces=faces)
    alphas={i:F(a) for i,a in cert['face_multipliers']}
    bound,beta,tight,num,den=C.sparse_dual_bound(d,alphas)
    assert bound>=F(cert['bound'])>14
    dual=dict(kind='EXACT_RESTRICTED_MODEL_DUAL',bound=str(bound),objective_numerator=str(num),normalizer=str(den),
        face_multipliers=cert['face_multipliers'],vertex_multipliers=[[i,str(x)] for i,x in sorted(beta.items())],
        tight_orbits={str(i):v for i,v in tight.items()},excludes_target12=bound>12)
    B.verify_certificate(d,dual)
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{C.TAG}.py',f'results/{INPUT}.json']
    d.update(status='PASS_EXACT_RANK14_A1_PADDED_TRIANGLE_MODEL_BARRIER',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=dict(vertices=len(vertices),variables=len(variables),triangles=len(faces)),
        thresholds=[thresholds(13),thresholds(14)],original_dual_bound=cert['bound'],padding_checks=checks,certificates=[dual],
        scope='Orthogonal A1 padding preserves each individual compact local triangle and its two-place arithmetic form, but does not establish simultaneous realization in one ambient polytope. The exact dual obstructs only nonnegative whole-forest corner-orbit rules required to pay all these local triangles.')
    assert d['thresholds'][0]['triangle_quad_payment_target']=='12'
    assert d['thresholds'][1]['triangle_quad_payment_target']=='14'
    return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['certificates'][0]['bound'],flush=True)
