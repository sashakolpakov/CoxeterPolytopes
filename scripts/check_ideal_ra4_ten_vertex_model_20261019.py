#!/usr/bin/env python3
"""Exact Figure7 model check; universal classification/box geometry are inputs."""
import argparse, fcntl, json, resource, signal, sys
from fractions import Fraction as F
from hashlib import file_digest
from itertools import combinations
from pathlib import Path
import run_lean_batches_20261018 as guard

ROOT=Path(__file__).resolve().parents[1]
TAG='ideal_ra4_ten_vertex_model_20261019'
REPORT=ROOT/'results/ideal_ra4_ten_vertex_facet_exclusion_research_20261019.md'


class K(tuple):
    """Q[s]/(s^4-s^2-1), in the ordered basis1,s,s^2,s^3."""
    def __new__(cls,x=0):
        v=x if isinstance(x,(tuple,list)) else (x,0,0,0)
        assert len(v)==4
        return tuple.__new__(cls,map(F,v))
    def __add__(a,b): return K(tuple(x+y for x,y in zip(a,K(b))))
    __radd__=__add__
    def __neg__(a): return K(tuple(-x for x in a))
    def __sub__(a,b): return a+-K(b)
    def __rsub__(a,b): return K(b)+-a
    def __mul__(a,b):
        b=K(b); v=[F()]*7
        for i,x in enumerate(a):
            for j,y in enumerate(b): v[i+j]+=x*y
        for i in range(6,3,-1): v[i-2]+=v[i]; v[i-4]+=v[i]
        return K(v[:4])
    __rmul__=__mul__
    def __truediv__(a,b): return a*(1/F(b))
    def __pow__(a,n):
        assert isinstance(n,int) and n>=0
        out=K(1)
        while n:
            if n&1: out=out*a
            a=a*a; n//=2
        return out
    def __bool__(a): return any(x for x in a)


def rs(x): return (x>0)-(x<0)
def qm(a,b): return (a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])
def qs(a):
    u,v=a[0]+a[1]/2,a[1]/2; su,sv=rs(u),rs(v)
    if not su: return sv
    if not sv or su==sv: return su
    return su*rs(u*u-5*v*v)
def sign(x):
    x=K(x); a,b=(x[0],x[2]),(x[1],x[3]); sa,sb=qs(a),qs(b)
    if not sa: return sb
    if not sb or sa==sb: return sa
    aa,bb=qm(a,a),qm((0,1),qm(b,b))
    d=qs((aa[0]-bb[0],aa[1]-bb[1])); assert d
    return sa*d


S=K((0,1,0,0)); SINV=S**3-S
def P(n): return (S if n>=0 else SINV)**abs(n)
def dot(a,b): return sum((a[i]*b[i] for i in range(3)),K())-a[3]*b[3]
def enc(x):
    if isinstance(x,K): return list(map(str,x))
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [enc(v) for v in x]
    return x
def digest(p):
    with p.open('rb') as f: return file_digest(f,'sha256').hexdigest()
def memory():
    v=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin': v*=1024
    assert v<1024**3
    return v


def run():
    one=K(1); zero=K(); a,R,t=P(-3),P(-2),P(-4)
    assert S**4==S*S+1 and S*SINV==one
    assert a**4+4*a*a==one and R==(1+a*a)/2 and t==(1-a*a)/2
    assert R*P(3)==S and a*P(4)==S and sign(a)>0 and sign(1-a)>0
    faces=['ABCD','ABFE','EHIJ','CGHI','AEJ','ADJ','BFG','BCG','CDI','DIJ','EFH','FGH']
    lines={'ABFE':(one,zero,a),'EHIJ':(-one,zero,a),
           'AEJ':(zero,one,one),'EFH':(zero,-one,one)}
    spheres={'ABCD':(a,t,-2),'CGHI':(-a,-t,-2),'BCG':(a*t,-a*a,-5),
             'CDI':(-a*t,a*a,-5),'ADJ':(zero,one,-3),'DIJ':(-a,R,-4),
             'FGH':(zero,-one,-3),'BFG':(a,-R,-4)}
    walls={name:(x,y,c,c) for name,(x,y,c) in lines.items()}
    for name,(x,y,e) in spheres.items():
        q=x*x+y*y-P(2*e)
        walls[name]=tuple(z*P(-e) for z in (x,y,(q-1)/2,(q+1)/2))
    assert set(walls)==set(faces) and all(dot(w,w)==one for w in walls.values())
    points={'A':(a,one),'B':(a,-a*a),'C':(zero,zero),'D':(-a**3,1-2*a*a),
            'F':(a,-one),'G':(a**3,-1+2*a*a),'H':(-a,-one),
            'I':(-a,a*a),'J':(-a,one)}
    vertices={v:(x,y,(x*x+y*y-1)/2,(x*x+y*y+1)/2) for v,(x,y) in points.items()}
    vertices['E']=(zero,zero,one,one)
    incidence={}
    for v,x in vertices.items():
        assert dot(x,x)==zero
        products={f:dot(x,walls[f]) for f in faces}
        assert all(sign(y)<=0 for y in products.values())
        incidence[v]=sorted(f for f,y in products.items() if not y)
        assert incidence[v]==sorted(f for f in faces if v in f) and len(incidence[v])==4
    pairs=[]; edges=set()
    for f,g in combinations(faces,2):
        common=sorted(set(f)&set(g)); product=dot(walls[f],walls[g])
        assert len(common)<=2
        if len(common)==2:
            assert product==zero; edges.add(tuple(common))
        elif len(common)==1: assert product==-one
        else: assert sign(product+1)<0
        pairs.append(dict(faces=[f,g],product=product,ideal_vertices=common))
    expected={tuple(sorted((f[i],f[(i+1)%len(f)]))) for f in faces for i in range(len(f))}
    assert edges==expected and len(edges)==20 and len(vertices)-len(edges)+len(faces)==2
    permutations=[]
    for cycles in ('AI BH CE DJ FG','AB CE DF GJ HI'):
        perm={v:v for v in vertices}
        for xy in cycles.split(): perm[xy[0]],perm[xy[1]]=xy[1],xy[0]
        assert {frozenset(perm[v] for v in f) for f in faces}=={frozenset(f) for f in faces}
        assert {tuple(sorted((perm[u],perm[v]))) for u,v in edges}==edges
        permutations.append(perm)
    half={v:permutations[0][permutations[1][v]] for v in vertices}
    assert half==dict(zip('ABCDEFGHIJ','HICGEJDABF'))
    orbits=[]; remaining=set(vertices)
    while remaining:
        orbit={min(remaining)}
        while True:
            more=orbit|{p[v] for p in permutations for v in orbit}
            if more==orbit: break
            orbit=more
        remaining-=orbit; orbits.append(''.join(sorted(orbit)))
    assert sorted(orbits)==['ABHI','CE','DFGJ']
    cusps={}
    for v in sorted(vertices):
        planes={}
        for f in incidence[v]:
            if f in lines:
                x,y,c=lines[f]; planes[f]=(x,y,c if v=='E' else zero,one)
            else:
                x,y,e=spheres[f]; px,py=points[v]
                x,y=x-px,y-py
                assert x*x+y*y==P(2*e)
                planes[f]=(x,y,K(F(1,2)),P(-2*e))
        opposite=[(f,g) for f,g in combinations(incidence[v],2)
                  if not(planes[f][0]*planes[g][1]-planes[f][1]*planes[g][0])]
        assert len(opposite)==2 and len(set(sum((list(p) for p in opposite),[])))==4
        for f,g in combinations(incidence[v],2):
            if (f,g) not in opposite: assert not(planes[f][0]*planes[g][0]+planes[f][1]*planes[g][1])
        widths=[]
        for f,g in opposite:
            x,y,c,u=planes[f]; X,Y,C,U=planes[g]
            width=c*c*u+C*C*U-2*c*C*(x*X+y*Y)*u*U
            assert sign(width)>0; widths.append(width)
        i=0 if sign(widths[0]-widths[1])>0 else 1
        exponent=3 if v in 'CE' else 1
        assert widths[i]==widths[1-i]*P(2*exponent)
        long=list(opposite[1-i]); quad=[f for f in incidence[v] if len(f)==4]
        if v in 'ABHI':
            assert len(quad)==2 and len(set(long)&set(quad))==1
        elif v in 'DFGJ': assert len(quad)==1 and set(quad)<=set(long)
        else: assert len(quad)==2 and set(quad)<=set(long)
        cusps[v]=dict(widths_squared=widths,opposite_faces=opposite,long_faces=long,
                      quadrilateral_faces=quad,aspect_power_of_s=exponent)
    assert set(cusps['A']['widths_squared'])=={P(6)/4,P(4)/4}
    assert set(cusps['J']['widths_squared'])=={P(6)/4,P(8)/4}
    assert set(cusps['E']['widths_squared'])=={4*P(-6),K(4)}
    assert 'AEJ' in cusps['A']['long_faces']
    assert all('AEJ' not in cusps[v]['long_faces'] for v in 'EJ')
    return enc(dict(status='PASS_EXACT_TEN_VERTEX_B_MODEL',
        sources={str(p.relative_to(ROOT)):digest(p) for p in (Path(__file__),REPORT)},
        field='Q(s), s^4-s^2-1=0, s>0',walls=walls,vertices=vertices,
        incidence=incidence,face_cycles=faces,edges=sorted(edges),wall_pairs=pairs,
        permutations=permutations,half_turn=half,cusp_orbits=sorted(orbits),cusps=cusps,
        a=P(-3),R=P(-2),t=P(-4),AEJ_signature=['long','short','short'],
        counts=dict(walls=12,ideal_vertices=10,edges=20,wall_pairs=66),
        universal_classification_checked=False,universal_box_proof_checked=False,
        scope='Exact finite12-wall model, feasible null points, incidence cycles, right/tangent/separated products, combinatorial symmetries and inversion cusp rectangles. Completeness of the external classification, geometric uniqueness, universal box closure and structural theorem are written inputs, not finite-model conclusions.'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write',action='store_true'); p.add_argument('--verify',action='store_true')
    args=p.parse_args(); assert not(args.write and args.verify)
    target=ROOT/'results'/(TAG+'.json')
    if args.write: assert not target.exists(), 'refuse to overwrite evidence'
    guard.resources(); memory(); signal.alarm(59)
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        data=json.loads(json.dumps(run())); guard.resources(start=False); peak=memory()
        if args.write: target.write_text(json.dumps(data,indent=2)+'\n')
        if args.verify: assert data==json.loads(target.read_text())
    signal.alarm(0)
    print(json.dumps(dict(status=data['status'],counts=data['counts'],peak_rss_bytes=peak)),flush=True)


if __name__=='__main__': main()
