#!/usr/bin/env python3
"""Structural comparison only: E8+1 versus the closed tripled D9 host."""
import argparse
import hashlib
import json
from pathlib import Path
from fresh_h9_glued_host_20261013_arithmetic import *

ROOT=Path(__file__).resolve().parents[1]
PREFIX='fresh_h9_tripled_d9_20261013_e8_triage'
SOURCE=ROOT/'results/fresh_h9_tripled_d9_20261013_symmetry.json'
FIELD=ROOT/'results/fresh_h9_glued_host_20261013_structure.json'
OUT=ROOT/'results'/f'{PREFIX}.json'
BETA=3*ALPHA
G=tuple(tuple(-BETA if i==j==0 else O if i==j else Z for j in range(10))for i in range(10))

def mm(a,b):return tuple(tuple(dot(row,col)for col in zip(*b))for row in a)
def det(a):
    a=[list(row)for row in a];d=O
    for j in range(len(a)):
        k=next(k for k in range(j,len(a))if a[k][j])
        if k!=j:a[j],a[k]=a[k],a[j];d=-d
        p=a[j][j];d*=p
        for k in range(j+1,len(a)):
            t=a[k][j]/p
            a[k]=[x-t*y for x,y in zip(a[k],a[j])]
    return d
def unit(i):return tuple(O if j==i else Z for j in range(10))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def reflect(r,v):return sub(v,tuple(2*pair(r,G,v)/pair(r,G,r)*x for x in r))

def run(verify=False):
    h=(Z,)+(-O/2,)*8+(Z,)
    e8=(h,)+tuple(sub(unit(i),unit(i+1))for i in range(2,8))+(add(unit(7),unit(8)),)
    basis=(unit(0),)+e8+(unit(9),)
    B=tuple(zip(*basis));BI=inverse(B)
    gram=tuple(tuple(pair(r,G,s)for s in basis)for r in basis)
    eg=tuple(tuple(pair(r,G,s)for s in e8)for r in e8)
    assert det(eg)==1 and det(gram)==-BETA and det(B)in(-1,1)
    assert all(eg[i][i]==2 for i in range(8))
    assert all(eg[i][j]in(0,-1)for i in range(8)for j in range(8)if i!=j)
    assert all(x.is_integral for row in gram for x in row)
    gi=inverse(gram)
    assert all((3*x).is_integral for row in gi for x in row)
    assert any(not x.is_integral for row in gi for x in row)
    pool=tuple(v*q for v in(1,3)for q in NORMS)
    assert len(pool)==10
    for q in pool:
        assert all(q.sign(e)>0 for e in range(4))and(6/q).is_integral
    def in_e8(v):return all(x.is_integral for x in mv(BI,v))
    def in_d9(v):return all(x.is_integral for x in v)and(sum(v[1:],Z)/2).is_integral
    def preserves(r):return all(in_e8(reflect(r,b))for b in basis)
    assert preserves(h)and pair(h,G,h)==2
    v=add(unit(1),unit(2));hv=reflect(h,v)
    assert in_d9(v)and not in_d9(hv)
    flip=tuple(2*x for x in unit(1));fh=reflect(flip,h)
    assert in_e8(h)and not in_e8(fh)and pair(flip,G,flip)==4
    # Two distinct maximal integral lattices share time+D8+<1>.
    d8=tuple(sub(unit(i),unit(i+1))for i in range(1,8))+(add(unit(7),unit(8)),)
    common=(unit(0),)+d8+(unit(9),)
    CB=tuple(zip(*common))
    assert det(CB)in(-2,2)and det(mm(BI,CB))in(-2,2)
    assert all(x.is_integral for row in mm(BI,CB)for x in row)
    assert not in_e8(unit(1))and not all(x.is_integral for x in h)
    source=json.loads(SOURCE.read_text())
    F=tuple(tuple(K.from_json(x)for x in row)for row in source['matrix'])
    Fi=inverse(F)
    failures=[]
    for name,A in(('transport',F),('inverse',Fi)):
        C=mm(mm(BI,A),B)
        bad=[{'row':i,'column':j,'value':C[i][j].to_json()}for i in range(10)for j in range(10)if not C[i][j].is_integral]
        assert len(bad)==2
        failures.append({'map':name,'nonintegral_entries':bad})
    assert json.loads(FIELD.read_text())['prime3_irreducibility_table']
    out={'status':'PASS_E8_NEIGHBOR_STRUCTURAL_TRIAGE_ONLY',
         'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,FIELD)},
         'coefficient':BETA.to_json(),'ambient_gram':menc(G),'e8_basis':menc(e8),'e8_gram':menc(eg),
         'full_lattice_basis':menc(basis),'full_lattice_determinant':(-BETA).to_json(),
         'dual_exponent':3,'safe_norm_divisor':6,'safe_norm_pool':enc(pool),
         'zero_chamber':'E8+A1','zero_mirrors':121,
         'maximality':'The determinant ideal is the inert prime (3), with valuation one. A proper integral overclosure would have a nontrivial index ideal whose square divides (3), which is impossible.',
         'neighbor_intersection':'<-3alpha>+D8+<1> is the exact intersection with <-3alpha>+O^9. Its index ideal in each is (2), cardinality16. Both larger lattices are maximal integral; they are unequal in these coordinates.',
         'norm_pool_scope':'Safe primitive pool only, using the established PID and totally-positive-unit-square facts: dual exponent3 gives q|6. No positive-height root census is claimed.',
         'e8_reflection_not_d9':{'root':enc(h),'d9_vector':enc(v),'image':enc(hv)},
         'd9_reflection_not_e8':{'root':enc(flip),'e8_vector':enc(h),'image':enc(fh),'image_module_coordinates':enc(mv(BI,fh))},
         'old_transport_failures':failures,
         'scope':'Available lattice-preserving mirrors differ in both directions in the fixed ambient identification, and the saved infinite chamber transport fails E8 lattice preservation. This does not exclude another indefinite lattice isometry, conjugacy of full reflection groups, or a different nonreflectivity proof. No new polytope or finite-volume chamber is claimed.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'safe_norms',len(pool),'old_transport_nonintegral',[len(x['nonintegral_entries'])for x in failures],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
