#!/usr/bin/env python3
"""All-real escape gate for the source coefficient refinement of cap34."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_e8_interior_caps_20261014_exact_review as K

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_replacement_escape_20261014_continuous'
BASE = 'results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
PARENT = 'results/fresh_h20_e7_intrinsic_caps_20261014.json'
ESCAPE = 'results/fresh_h20_e7_replacement_escape_20261014.json'

def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,(tuple,list)): return [enc(y) for y in x]
    if isinstance(x,dict): return {str(k):enc(v) for k,v in x.items()}
    return x

def run():
    b,p,e = [json.loads((ROOT/x).read_text()) for x in [BASE,PARENT,ESCAPE]]
    for doc in [b,p,e]:
        for path,h in doc['sources'].items(): assert sha(path)==h,path
    Q=[[F(x) for x in row] for row in b['form']]
    sparse=[[(j,a) for j,a in enumerate(row) if a] for row in Q];cache={}
    def dot(x,y):
        if y not in cache: cache[y]=tuple(tuple(sum(a*y[j][k] for j,a in row) for k in range(8)) for row in sparse)
        v=K.Z
        for a,bb in zip(x,cache[y]):
            if any(a) and any(bb):v=K.add(v,K.multiply(a,bb))
        return v
    N=list(map(K.vector,b['normals']));E=K.vector(e['point']);H=K.vector(b['strict_center'])
    ep=[dot(E,n) for n in N]
    assert all(K.sign(v)<0 for v in ep)
    assert K.sign(dot(E,E))>0 and K.sign(dot(E,H))<0 and K.sign(dot(H,H))<0
    a11,a12=p['face_facet_walls'][11:13]
    assert dot(N[a11],N[a11])==dot(N[a12],N[a12])==K.scalar(2)
    assert dot(N[a11],N[a12])==K.scalar(-1)
    L=K.vadd(K.vscale(N[a11],K.scalar(F(1,3))),K.vscale(N[a12],K.scalar(F(2,3))))
    sq2=(F(0),F(1))+(F(0),)*6
    T=K.vscale(K.vector(p['face_basis'][-1]),K.scale(sq2,F(1,2)))
    cap=K.vector(e['selected_cap']['normal'])
    V=K.vadd(K.vadd(cap,K.vscale(L,K.scalar(2))),K.vscale(T,K.scalar(-5)))
    assert [dot(V,V),dot(L,L),dot(T,T)]==list(map(K.scalar,[F(148,3),F(2,3),-2]))
    assert dot(V,L)==dot(V,T)==dot(L,T)==K.Z
    # r(x)=V-xL+sqrt((71+x²)/3)T has norm2, for every real x.
    C,a,B=dot(V,E),K.neg(dot(L,E)),K.neg(dot(T,E))
    assert all(K.sign(x)>0 for x in [C,a,B])
    D=K.sub(K.multiply(B,B),K.scale(K.multiply(a,a),3))
    margin=K.sub(K.scale(D,F(71,3)),K.multiply(C,C))
    assert K.sign(D)>0 and K.sign(margin)>0
    # Lorentz Cauchy inequality gives B sqrt((71+x²)/3)-a x
    # >=sqrt(71D/3)>C, hence E.r(x)<0 on the entire real line.
    paths=[BASE,PARENT,ESCAPE,f'scripts/{TAG}.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
    return enc(dict(status='PASS_ALL_REAL_SOURCE_COEFFICIENT_REFINEMENTS_RETAIN_STRICT_ESCAPE',
        sources={x:sha(x) for x in paths},point=E,norm=dot(E,E),all87_pairings=ep,
        family={'constant_vector':V,'source_weight':L,'time_vector':T,'time_square':'(71+x*x)/3','norm':F(2),'parameter':'all real x, positive time square root'},
        bound={'C':C,'a':a,'B':B,'D=B2-3a2':D,'positive_gap=71D/3-C2':margin,
               'inequality':'E.r(x) <= C-sqrt(71D/3) < 0 for every real x'},
        refined_contact={'x_square':F(86,11),'time_square':F(289,11)},
        scope='The exact escape is strict for all87 original walls and every member of this one-parameter family. Thus no finite coherent completion using any subset of the old walls and only these family walls has finite volume. Other profiles or source coefficients remain outside this statement.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(path.read_text())
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
