#!/usr/bin/env python3
"""Scout integral isometries between finite H3 vertex stars.

Candidate square roots come from PARI; every returned field square and every
isometry is checked in exact project arithmetic. No claim of exhaustive
isometry enumeration is made.
"""
import argparse
from functools import lru_cache
from itertools import combinations,permutations
import json
from pathlib import Path
import subprocess

from fresh_cubic9_dihedral_h9_20261013 import (
    K,T,Z,O,decode_matrix,encode_matrix,enc,pair,mv,dot,inverse)

ROOT=Path(__file__).resolve().parents[1]

def transpose(a):return tuple(zip(*a))
def mm(a,b):return tuple(tuple(dot(r,c) for c in transpose(b)) for r in a)
def eye(n):return tuple(tuple(K(int(i==j)) for j in range(n)) for i in range(n))
def expr(q):return f'({q.a})+({q.b})*y+({q.c})*y^2'

def squares(values):
    values=list(dict.fromkeys(values))
    lines=['nf=nfinit(y^3-3*y-1);']
    for q in values:
        lines.append(f'r=nfroots(nf,x^2-Mod({expr(q)},y^3-3*y-1));if(#r,print(vector(3,i,polcoef(lift(r[1]),i-1))),print("NONE"));')
    proc=subprocess.run(['gp','-q'],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
    output=proc.stdout.strip().splitlines();assert len(output)==len(values),(proc.stderr,proc.stdout)
    result={}
    for q,line in zip(values,output):
        if line=='NONE':result[q]=None;continue
        assert line.startswith('[') and line.endswith(']'),line
        root=K(*map(str.strip,line[1:-1].split(',')));assert root*root==q
        result[q]=root if root.sign()>0 else -root
    return result

def run(source,output):
    d=json.loads(source.read_text());assert d['dimension']==3 and d['completed']
    G=decode_matrix(d['gram']);walls=decode_matrix(d['all_walls']);alpha=K.from_json(d['alpha'])
    H=K.from_json(d['height_bound']);p=tuple(map(K.from_json,d['strict_interior']))
    cov=tuple(mv(G,r) for r in walls)
    vertex=[]
    for row in d['geometry']['vertices']:
        v=tuple(map(K.from_json,row['vector']));q=pair(v,G,v)
        if q.sign()<0:vertex.append((tuple(row['active_facets']),v,q))
    matches=[]
    for si,ti in combinations(range(len(vertex)),2):
        ids,v,q=vertex[si];jds,w,s=vertex[ti]
        for perm in permutations(jds):
            if all(pair(walls[i],G,walls[j])==pair(walls[k],G,walls[l])
                   for i,k in zip(ids,perm) for j,l in zip(ids,perm)):
                matches.append((si,ti,perm,q/s))
    roots=squares(m[-1] for m in matches);isometries=[];seen=set()
    for si,ti,perm,ratio in matches:
        t=roots[ratio]
        if t is None:continue
        ids,v,q=vertex[si];jds,w,s=vertex[ti]
        source_frame=transpose((v,*(walls[i] for i in ids)))
        target_frame=transpose((tuple(t*x for x in w),*(walls[i] for i in perm)))
        U=mm(target_frame,inverse(source_frame))
        if U in seen:continue
        seen.add(U);assert mm(mm(transpose(U),G),U)==G
        if not all(x.is_integral for row in U for x in row):continue
        inv=inverse(U)
        if not all(x.is_integral for row in inv for x in row):continue
        assert mv(U,v)[0].sign()>0
        isometries.append({'source_vertex':si,'target_vertex':ti,'target_wall_order':list(perm),'matrix':U})
    Hp=alpha*p[0]*p[0]/(-pair(p,G,p))-1
    assert Hp<=H and all(dot(p,c).sign()<0 for c in cov)
    def gate(U,point=p):
        u=mv(U,point);q=pair(u,G,u);height=alpha*u[0]*u[0]/(-q)-1
        ps=tuple(dot(u,c) for c in cov);tr=sum((U[i][i] for i in range(4)),Z)
        point_height=alpha*point[0]*point[0]/(-pair(point,G,point))-1
        return {'point':enc(point),'point_height_bound':point_height.to_json(),
                'image':enc(u),'height_bound':height.to_json(),'trace':tr.to_json(),
                'all_prefix_pairings_strict':all(t.sign()<0 for t in ps),
                'height_within_complete_prefix':height<=H and point_height<=H,
                'trace_proves_infinite_order':tr>4,'wall_pairings':enc(ps)}
    reflections=[]
    for r,c in zip(walls,cov):
        q=dot(r,c)
        reflections.append(tuple(tuple(K(int(i==j))-2*r[i]*c[j]/q for j in range(4)) for i in range(4)))
    def reduce_image(U):
        word=[]
        for _ in range(1000):
            u=mv(U,p)
            bad=next((i for i,c in enumerate(cov) if dot(u,c).sign()>0),None)
            if bad is None:return U,word,True
            U=mm(reflections[bad],U);word.append(bad)
        return U,word,False
    checks=[];successful=[];matrices=[]
    for row in isometries:
        U=row['matrix'];matrices.extend((U,inverse(U)))
        reduced,word,completed=reduce_image(U)
        if completed:matrices.extend((reduced,inverse(reduced)))
        result=gate(U);checks.append({**row,'matrix':encode_matrix(U),'gate':result,
            'reduced_matrix':encode_matrix(reduced),'reduction_word':word,
            'reduction_completed':completed,'reduced_gate':gate(reduced)})
    # Products can turn two vertex involutions into an infinite-order map.
    for U in list(dict.fromkeys([*matrices,*(mm(a,b) for a in matrices for b in matrices)])):
        if sum((U[i][i] for i in range(4)),Z)<=4:continue
        probes=[p]
        for _,v,_ in vertex:
            probes.extend(tuple(x+eps*y for x,y in zip(v,p)) for eps in (K('1/16'),K('1/256'),K('1/4096')))
        for point in probes:
            result=gate(U,point)
            if result['all_prefix_pairings_strict'] and result['height_within_complete_prefix']:
                assert all(dot(point,c).sign()<0 for c in cov)
                successful.append({'matrix':encode_matrix(U),'gate':result});break
        if successful:break
    out={'status':'EXACT_CHAMBER_STABILIZER_CANDIDATE' if successful else 'NO_STABILIZER_IN_DECLARED_VERTEX_SCOUT',
         'source':str(source.relative_to(ROOT)),'finite_vertices':len(vertex),'matched_stars':len(matches),
         'integral_isometries':len(isometries),'strict_point':enc(p),'strict_point_height_bound':Hp.to_json(),
         'vertex_isometries':checks,'successful_matrices':successful,
         'scope':'Vertex-star isometry scout and products of length two. Any positive requires independent complete mirror/height verification before a full-host theorem.'}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:out[k] for k in ('source','status','finite_vertices','matched_stars','integral_isometries')}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--alpha',choices=('theta','2theta','theta_beta','2theta_beta'),required=True)
    p.add_argument('--height',type=int,default=300)
    a=p.parse_args();run(ROOT/f'results/fresh_cubic9_dihedral_h3_{a.alpha}_h{a.height}_20261013.json',
                        ROOT/f'results/fresh_cubic9_vertex_symmetry_{a.alpha}_h{a.height}_20261013.json')
