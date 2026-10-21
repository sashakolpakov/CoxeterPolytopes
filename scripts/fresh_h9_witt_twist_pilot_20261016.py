#!/usr/bin/env python3
"""Complete bounded Vinberg pilot for <-phi> + I7 + (4-phi) I2."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import cmp_to_key
from hashlib import sha256
from math import ceil,isqrt
from pathlib import Path
import fresh_h9_beta11_vinberg_h300_20261014_review as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_witt_twist_pilot_20261016'
SETUP='results/fresh_h9_witt_twist_20261016.json'
HELPER='scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py'
K=A.K;PHI=K(0,1);BETA=4-PHI;HEIGHT=F(20)
DIAGONAL=[-PHI]+[K(1)]*7+[BETA]*2
NORMS=[K(1),K(2),BETA,2*BETA]
ORDER=[9,8,7,6,5,4,3,2,1]

def sha(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def read(path):
    data=json.loads((ROOT/path).read_text())
    for p,h in data['sources'].items():assert sha(p)==h,p
    return data
def pair(x,y):return sum((d*a*b for d,a,b in zip(DIAGONAL,x,y)),K())
def ceiling_sqrt(q):
    n=isqrt(q.numerator//q.denominator)
    while n*n<q:n+=1
    return n
def box(physical,conjugate):
    return ceil(F(17*conjugate+7*physical,22)),ceil(F(5*(physical+conjugate),11))
def zero_roots():
    roots=[]
    for i in range(1,7):
        r=[K()]*10;r[i]=K(-1);r[i+1]=K(1);roots.append(r)
    r=[K()]*10;r[7]=K(-1);roots.append(r)
    r=[K()]*10;r[8]=K(-1);r[9]=K(1);roots.append(r)
    r=[K()]*10;r[9]=K(-1);roots.append(r)
    return roots

def run(height=HEIGHT):
    height=F(height);assert height>0
    setup=read(SETUP);assert setup['status'].startswith('PASS')
    assert [[A.dec(x) for x in row] for row in setup['form']]==[
        [d if i==j else K() for j in range(10)] for i,d in enumerate(DIAGONAL)]
    assert [A.dec(x) for x in setup['norm_and_zero']['primitive_norm_classes']]==NORMS
    assert [[A.dec(x) for x in row] for row in setup['norm_and_zero']['zero_simple_walls']]==zero_roots()
    assert setup['norm_and_zero']['zero_root_count']==106
    assert BETA*BETA.conjugate()==11
    assert all(q.sign()>0 and q.sign(True)>0 for q in NORMS)
    assert (PHI-F(8,5)).sign()>0 and (-PHI.conjugate()-F(3,5)).sign()>0
    assert (F(17,10)-PHI).sign()>0 and (BETA-2).sign()>0 and (BETA.conjugate()-4).sign()>0
    assert all((5-q).sign()>0 and (10-q.conjugate()).sign()>0 for q in NORMS)
    xb=box(ceiling_sqrt(F(25,8)*height),5)
    yb=box(ceiling_sqrt(5*(1+height)),4)
    zb=box(ceiling_sqrt(F(5,2)*(1+height)),2)
    controls=[]
    for q in NORMS:
        long=q in (BETA,2*BETA)
        for a in range(-xb[0],xb[0]+1):
            for b in range(-xb[1],xb[1]+1):
                x=K(a,b)
                if x.sign()<=0 or (long and not (x/BETA).integral()):continue
                h=PHI*x*x/q;target=q+PHI*x*x
                if (height-h).sign()>=0 and target.sign(True)>=0:controls.append((h,q,x))
    controls.sort(key=lambda row:row[0])
    assert all(a[0]!=b[0] for a,b in zip(controls,controls[1:]))
    simple=zero_roots();previous=[];levels=[];totals=Counter()
    sortkey=cmp_to_key(lambda a,b:A.sign(A.sub(a,b)))
    for level,(h,q,x) in enumerate(controls):
        target=q+PHI*x*x;ti=target.ri();long=q in (BETA,2*BETA)
        alphabets=[]
        for weight,bounds,is_y in ((K(1),yb,True),(BETA,zb,False)):
            values=[]
            for a in range(-bounds[0],bounds[0]+1):
                for b in range(-bounds[1],bounds[1]+1):
                    value=K(a,b)
                    if value.sign()<0 or (is_y and long and not (value/BETA).integral()):continue
                    remainder=target-weight*value*value
                    if remainder.sign()>=0 and remainder.sign(True)>=0:values.append((a,b))
            values.sort(key=sortkey);alphabets.append(values)
        ya,za=alphabets
        coefficients=[[(DIAGONAL[j]*r[j]).ri() for j in ORDER] for r in previous]
        limits=[(PHI*x*r[0]).ri() for r in previous]
        assert all(A.sign(c)>=0 for row in coefficients for c in row)
        visits=[0]*10;energy_prunes=[0]*10;wall_prunes=[0]*10;roots=[]
        def recurse(prefix,energy,partial,last_index=0):
            depth=len(prefix);visits[depth]+=1
            remainder=A.sub(ti,energy)
            if A.sign(remainder)<0 or A.sign(remainder,True)<0:
                energy_prunes[depth]+=1;return
            if depth:
                z=prefix[-1]
                partial=[A.add(v,A.mul(z,row[depth-1])) for v,row in zip(partial,coefficients)]
                if any(A.sign(A.sub(bound,v))<0 for bound,v in zip(limits,partial)):
                    wall_prunes[depth]+=1;return
                remaining=1 if depth==1 else 9-depth if depth>=3 else 0
                weight=BETA if depth==1 else K(1)
                lower=A.scale(A.mul(weight.ri(),A.mul(z,z)),remaining)
                if A.sign(A.sub(remainder,lower))<0:
                    energy_prunes[depth]+=1;return
            if depth==9:
                if energy!=ti:return
                r=[x]+[K()]*9
                for i,z in zip(ORDER,prefix):r[i]=K(*z)
                assert pair(r,r)==q and all(v.integral() for v in r)
                assert all((2*d*v/q).integral() for d,v in zip(DIAGONAL,r))
                assert all(pair(r,n).sign()<=0 for n in simple+previous)
                roots.append(r);return
            alphabet=za if depth<2 else ya
            start=last_index if depth in (1,3,4,5,6,7,8) else 0
            weight=DIAGONAL[ORDER[depth]].ri()
            for index in range(start,len(alphabet)):
                z=alphabet[index]
                next_energy=A.add(energy,A.mul(weight,A.mul(z,z)))
                recurse(prefix+[z],next_energy,partial,index)
        recurse([],(0,0),[(0,0)]*len(previous))
        assert all(pair(a,b).sign()<=0 for i,a in enumerate(roots) for b in roots[i+1:])
        previous+=roots;totals['nodes']+=sum(visits)
        levels.append(dict(height=h,norm=q,control=x,target=target,y_alphabet=[K(*z) for z in ya],
                           z_alphabet=[K(*z) for z in za],visits=visits,energy_prunes=energy_prunes,
                           wall_prunes=wall_prunes,roots=roots))
        print('level',level+1,'of',len(controls),'walls',9+len(previous),'nodes',sum(visits),flush=True)
    normals=simple+previous;gram=[[pair(r,s) for s in normals] for r in normals]
    assert all(gram[i][j].sign()<=0 for i in range(len(normals)) for j in range(i))
    angle_values={K():2,K(1):3,K(2):4,K(1,1):5,K(3):6,K(2,1):10}
    pair_rows=[];pair_counts=Counter()
    for i in range(len(normals)):
        for j in range(i):
            value=4*gram[i][j]*gram[i][j]/(gram[i][i]*gram[j][j])
            assert value!=4
            if (value-4).sign()>0:kind='dotted'
            else:
                assert value in angle_values
                label=angle_values[value];kind='right' if label==2 else f'pi/{label}'
            pair_counts[kind]+=1;pair_rows.append(dict(pair=[j,i],cartan_cycle=value,kind=kind))
    H=[K(1)]+[K()]*9;shift=[K()]+list(map(K,range(7,0,-1)))+[K(2),K(1)]
    epsilon=F(1)
    while True:
        center=[a+epsilon*b for a,b in zip(H,shift)]
        if pair(center,center).sign()<0 and all(pair(center,n).sign()<0 for n in normals):break
        epsilon/=2
    facets=[]
    for i,n in enumerate(normals):
        p=[a-pair(center,n)/gram[i][i]*b for a,b in zip(center,n)]
        assert pair(p,p).sign()<0
        assert all(pair(p,m).sign()==(0 if i==j else -1) for j,m in enumerate(normals))
        facets.append(p)
    paths=[f'scripts/{TAG}.py',HELPER,SETUP]
    return A.enc(dict(status='PASS_COMPLETE_BOUNDED_WITT_TWIST_PREFIX',sources={p:sha(p) for p in paths},
        form=[[d if i==j else K() for j in range(10)] for i,d in enumerate(DIAGONAL)],
        height_bound=height,basepoint=H,norm_classes=NORMS,zero_simple_walls=simple,
        uniform_coefficient_boxes=dict(control=xb,y=yb,z=zb),positive_coordinate_order=ORDER,
        levels=levels,roots=normals,gram=gram,strict_center=center,center_epsilon=epsilon,actual_facet_points=facets,
        pair_relations=pair_rows,pair_counts=dict(pair_counts),
        counts=dict(controls=len(controls),walls=len(normals),positive_walls=len(previous),**totals),
        completeness='Four complete primitive norm classes and complete B7+B2 zero sector from the bound setup. Explicit rational coefficient boxes followed by exact both-embedding norm checks. Dominant coordinates are ordered within their seven- and two-coordinate blocks. All remaining physical wall contributions are nonnegative; no sign restriction is placed on conjugate coordinates.',
        scope='A complete bounded-height Vinberg prefix on the locally distinct diagonal lattice. No volume or reflectivity conclusion from wall counts.'))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text());print('PASS saved Witt-twist pilot replay')
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
