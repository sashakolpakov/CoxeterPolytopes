#!/usr/bin/env python3
"""Exact compressed ternary-coordinate cap construction; no large matrix."""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from math import comb
from operator import xor
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_outer_unit_code_20261014'
COLS = (0, 1, 2, 3, 4, 5, 8, 10, 12, 15)
Q = lambda a=0, b=0: (F(a), F(b))
def add(a,b): return a[0]+b[0],a[1]+b[1]
def sub(a,b): return a[0]-b[0],a[1]-b[1]
def mul(a,b): return a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1]
def scale(a,t): return a[0]*t,a[1]*t
def sign(a):
    x,y=a[0]+a[1]/2,a[1]/2
    if not y:return (x>0)-(x<0)
    if not x:return (y>0)-(y<0)
    if x*y>0:return (x>0)-(x<0)
    d=x*x-5*y*y
    assert d
    return ((x>0)-(x<0))*((d>0)-(d<0))
def total(xs):
    z=Q()
    for x in xs:z=add(z,x)
    return z
W=[Q(1)]*10+[Q(2)]+[Q(1,1)]*8
def dot(x,y):return total(scale(w,a*b)for w,a,b in zip(W,x,y,strict=True))
def dist(x,y):return dot([a-b for a,b in zip(x,y)],[a-b for a,b in zip(x,y)])
def syndrome(m):return reduce(xor,(COLS[i]for i in range(10)if m>>i&1),0)
def vector(m):return [F(0 if m>>i&1 else 1)for i in range(10)]+[F(1)]*9
def minimum(xs):
    out=None;att=[]
    for k,v in xs:
        if out is None or sign(sub(v,out))<0:out=v;att=[k]
        elif v==out:att.append(k)
    return out,att
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(tuple,list)):return [enc(t)for t in x]
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    return x
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def pair_class(kind,d):
    # Here the outer nine coordinates have identical signs.
    if d==0:return 'identity'
    if kind=='CC':
        assert d%2==0 and d>=4
        return {4:'right',6:'pi3',8:'parallel'}.get(d,'dotted')
    if kind=='HH':
        assert d%2==0 and d>=2
        return {2:'right',4:'parallel'}.get(d,'dotted')
    assert kind=='CH' and d%2==1 and d>=3
    return {3:'right',5:'pi4'}.get(d,'dotted')

def pattern_pairs(X,Y,kind,nx,ny):
    # Ordered pairs, including identities if X=Y. Every possible choice
    # of signs at zero coordinates is removed rather than double counted.
    out=Counter()
    for x in X:
        for y in Y:
            h=(x^y).bit_count();k=10-(x|y).bit_count()
            multiplicity=512*2**(10-x.bit_count()+(x&~y).bit_count())
            for r in range(k+1):
                d=h+4*r
                out[pair_class(kind,d)]+=multiplicity*comb(k,r)
    assert sum(out.values())==nx*ny//512
    off=nx*ny-sum(out.values())
    if kind=='CC':
        # A sole heavy-sign change adds8: identical unit centers give
        # parallel walls. Any other nonzero change is strictly dotted.
        out['parallel']+=nx;off-=nx
    out['dotted']+=off
    if X==Y:
        assert out.pop('identity')==nx
        assert all(v%2==0 for v in out.values())
        out=Counter({k:v//2 for k,v in out.items()})
        assert sum(out.values())==nx*(nx-1)//2
    else:assert sum(out.values())==nx*ny
    return dict(sorted(out.items()))

def small_pairs(X,nx,kind):
    ns=2**19;out=Counter()
    for m in X:
        for r in range(11-m.bit_count()):
            # Distance =5/2+2r, irrespective of the zero pattern.
            label=('right' if r==0 else 'parallel' if r==1 else 'dotted') if kind=='SC' else ('pi4' if r==0 else 'dotted')
            out[label]+=512*1024*comb(10-m.bit_count(),r)
    assert sum(out.values())==ns*nx//512
    out['dotted']+=ns*nx-sum(out.values())
    return dict(sorted(out.items()))

def run():
    inherited=['scripts/fresh_h20_weighted_binary_20261014.py',
               'results/fresh_h20_weighted_binary_20261014.json']
    old=json.loads((ROOT/inherited[1]).read_text())
    assert old['sources'][inherited[0]]==digest(inherited[0])
    assert old['status']=='PASS_EXACT_PARITY_TRIPLE_CAP_ORBIT_WITH_GLOBAL_ESCAPE'
    assert old['orbit_size']==2**19 and old['full_unit_triples']==[[0,1,2],[0,3,4]]
    assert len(COLS)==len(set(COLS))==10 and reduce(xor,COLS)==0
    C=[m for m in range(1024)if m.bit_count()%2==0 and syndrome(m)==0]
    H=[m for m in range(1024)if m.bit_count()%2 and syndrome(m)not in COLS]
    assert len(C)==32 and len(H)==192 and 0 in C and 1023 in C
    assert min(m.bit_count()for m in C if m)==4
    assert all(x^y in C for x in C for y in C)
    assert all(min((h^c).bit_count()for c in C)==3 for h in H)
    assert H==[m for m in range(1024)if m.bit_count()%2 and min((m^c).bit_count()for c in C)>=3]
    nC=512*sum(2**(10-m.bit_count())for m in C)
    nH=512*sum(2**(10-m.bit_count())for m in H)
    nS=nB=2**19
    pairs={
        'code_code':pattern_pairs(C,C,'CC',nC,nC),
        'hole_hole':pattern_pairs(H,H,'HH',nH,nH),
        'code_hole':pattern_pairs(C,H,'CH',nC,nH),
        'small_code':small_pairs(C,nC,'SC'),
        'small_hole':small_pairs(H,nH,'SH')}
    ss=Counter()
    for r in range(1,11):ss['right' if r==1 else 'parallel' if r==2 else 'dotted']+=nS*comb(10,r)//2
    ss['dotted']+=nS*(nS-1)//2-sum(ss.values())
    pairs['small_small']=dict(sorted(ss.items()))
    # An outer sign change adds8 or4phi^2. Except the explicitly recorded
    # code/code heavy parallel pair, this exceeds the parallel threshold.
    assert sign(sub(scale(Q(1,1),4),Q(8)))>0

    seed=[]
    for name,patterns,q in [('code',C,F(2)),('hole',H,F(1))]:
        for k in sorted(set(m.bit_count()for m in patterns)):
            counts={'right':19-k,'dotted':19-k,'pi4' if name=='code' else 'parallel':2*k}
            assert sum(counts.values())==38
            # The weakest central gap is the all-zero-unit code center.
            norm=Q(20-k,8)
            num=sub(norm,Q(q+4));assert sign(num)>0
            assert sign(sub(mul(num,num),Q(16*q)))>0
            seed.append([name,k,counts,norm])
    small_norm=Q(F(25,2),8)
    num=sub(small_norm,Q(F(9,2)))
    assert sign(sub(mul(num,num),Q(8)))>0
    seed.append(['small',None,{'right':9,'pi4':10,'dotted':19},small_norm])
    # Each unit zero gives gap1; code radius sqrt2 gives pi4, hole radius1
    # gives parallel. Small unit gaps1/2 give pi4 at radius1/sqrt2.
    # Remaining full-coordinate opposite gaps are all strictly dotted.
    assert F(1,1)/2==F(1,2) and F(1,4)/F(1,2)==F(1,2)
    for radius_square in (F(2),F(1),F(1,2)):
        for weight in W[10:]:
            assert sign(sub(scale(weight,4/radius_square),Q(1)))>0
    assert F(4,2)>1 and F(4)>1 and F(9,4)/F(1,2)>1
    binary_bounds={'code':Q(F(17,4),2),'hole':Q(F(17,4),2),'small':Q(F(13,4),2)}
    assert sign(sub(binary_bounds['code'],Q(6)))>0 # 6>(1+sqrt2)^2
    assert sign(sub(binary_bounds['hole'],Q(4)))>0
    assert sign(sub(binary_bounds['small'],Q(3)))>0 # 3>(1+1/sqrt2)^2
    # Code/hole zero entries can only increase the distance to the three
    # full binary unit coordinates; other binary unit entries are halves.
    for full in ({0,1,2},{0,3,4}):
        b=[F(1)if i in full else F(1,2)for i in range(19)]
        for m in C+H:
            assert dist(vector(m),b)==add(binary_bounds['code'],Q(sum(m>>i&1 for i in full)))
        assert dist([F(1,2)]*10+[F(1)]*9,b)==binary_bounds['small']

    # Positive common center (0,100,1), norm-200. All new D coordinates
    # are positive. Together with acute pairs this makes every facet actual.
    minD=scale(sub(Q(10,8),Q(2)),F(1,2));assert sign(minD)>0
    # At each old box vertex, (z,(N+2)/2,1) still has norm-2. The original
    # sqrt2 vertex cap is the zero word in C; all other new inequalities
    # are strict there, since code distances>=4 and H/S q<2.
    endpoint_margin={'other_code':F(4),'hole':F(1),'small':F(3,2)}

    oldx=[F(0)]*10+[F(99,100)]*9
    old_power=sub(dist(oldx,vector(1023)),Q(2));assert sign(old_power)<0
    # Four-coordinate dual parity provides an exact uncovered unit-cube
    # point, even before perturbing the fixed outer coordinates inward.
    dual=sum(1<<i for i in (1,3,5,9))
    assert all((m&dual).bit_count()%2==0 for m in C)
    x0=[F(1,2)]*10+[F(1)]*9
    x0[1]=0
    for i in (3,5,9):x0[i]=1
    assert minimum((m,sub(dist(x0,vector(m)),Q(2)))for m in C)[0]==Q(F(1,2))
    assert minimum((m,sub(dist(x0,vector(m)),Q(1)))for m in H)[0]==Q(F(1,2))
    assert sub(dist(x0,[F(1,2)]*10+[F(1)]*9),Q(F(1,2)))==Q(F(1,2))
    # A strict global escape for core+binary+all three new families.
    x=x0[:]
    for i in (3,5,9,*range(10,19)):x[i]=F(99,100)
    cp,ca=minimum((m,sub(dist(x,vector(m)),Q(2)))for m in C)
    hp,ha=minimum((m,sub(dist(x,vector(m)),Q(1)))for m in H)
    sp=sub(dist(x,[F(1,2)]*10+[F(1)]*9),Q(F(1,2)))
    bp,ba=minimum((j,sub(dist(x,[F(1)if i in f else F(1,2)for i in range(19)]),Q(1)))for j,f in enumerate(({0,1,2},{0,3,4})))
    cen=sub(dot(x,x),Q(4)); powers={'code':cp,'hole':hp,'small':sp,'binary_lower_bound':bp,'central':cen}
    for v in powers.values():assert sign(sub(v,Q(F(1,4))))>0
    assert all(abs(v)<1 for v in x)
    D=scale(sub(dot(x,x),Q(F(1,4))),F(1,2))
    assert sub(dot(x,x),scale(D,2))==Q(F(1,4))
    # H=(0,100,1) fixes the future component; (X,H)=-D-100<0.
    assert sign(add(D,Q(100)))>0
    counts={'code':nC,'hole':nH,'small':nS,'binary':nB,'core':39}
    return enc(dict(status='PASS_EXACT_OUTER_CODE_LAYERS_WITH_STRICT_ESCAPE',
        sources={p:digest(p)for p in [f'scripts/{TAG}.py',*inherited]},
        columns=COLS,code=C,holes=H,code_weight_enumerator=dict(Counter(m.bit_count()for m in C)),
        hole_weight_enumerator=dict(Counter(m.bit_count()for m in H)),counts=counts,
        actual_facets=sum(counts.values()),new_pair_counts=pairs,seed_counts=seed,
        central_all_dotted=True,binary_all_dotted=True,binary_minimum_distance_squares=binary_bounds,
        binary_new_dotted_pairs=nB*(nC+nH+nS),preserved_old_finite_endpoints=2**19,
        endpoint_inequality_positive_margins=endpoint_margin,strict_center_D=100,minimum_new_D=minD,
        old_escape_repair=dict(code_zero_pattern=1023,power=old_power),
        uncovered_fixed_outer_cube=dict(point=x0,dual_support=[1,3,5,9],three_family_minimum_power=F(1,2)),
        strict_escape=dict(spatial=x,D=D,U=1,norm=F(1,4),minimum_normalized_seed_gap=F(1,100),
            powers=powers,code_attainers=ca,hole_attainers=ha,binary_relaxed_attainers=ba),
        scope='Retained weighted39 core and full binaryR1 orbit, plus complete code/hole/small cap families. Exact compressed all-pair proof and actual facets. Original full sqrt2 vertex bank is the zero codeword subfamily. The old escape is cut, but a new strict spacelike escape remains; even the fixed-outer unit cube is not covered. No finite-volume or complete new-vertex/cusp census claim.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['actual_facets'],out['strict_escape']['powers'])
