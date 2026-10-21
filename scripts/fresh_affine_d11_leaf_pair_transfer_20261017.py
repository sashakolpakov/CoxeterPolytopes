#!/usr/bin/env python3
"""Universal affine-D11 opposite-leaf transfer; a new finite coupon is explicit."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product,combinations
from pathlib import Path
import sympy as s
import fresh_affine_d12_leaf_pair_transfer_20261017 as B

ROOT=B.ROOT;M=B.M;C=B.C
TAG='fresh_affine_d11_leaf_pair_transfer_20261017'
SPECTATOR='fresh_spectator_five_quarters_transfer_20261017'
MARKS=list(product((8,9),(10,11)))

def affine():
    g=2*s.eye(12)
    for i,j in[(i,i+1)for i in range(7)]+[(0,8),(0,9),(7,10),(7,11)]:g[i,j]=g[j,i]=-1
    nu=s.Matrix([1]*8+[s.Rational(1,2)]*4)
    assert g*nu==s.zeros(12,1) and C.inertia(g.tolist())==[11,0,1]
    return g,nu

def profiles(mark):
    center=0 if mark<10 else 7;mate=17-mark if mark<10 else 21-mark
    t=s.symbols('t'if mark<10 else'u',nonnegative=True)
    return [dict(name='isolated',contacts={mark:t},c0=2,l=0,h=t/2,t=t),
        dict(name='mate3',contacts={mark:t,mate:1},c0=1,l=1,h=(t+1)/2,t=t),
        dict(name='mate4',contacts={mark:t,mate:s.sqrt(2)},c0=0,l=s.sqrt(2),h=(t+s.sqrt(2))/2,t=t),
        dict(name='center3',contacts={mark:t,center:1},c0=0,l=2,h=1+t/2,t=t)]

def formulas():
    g,nu=affine();rows=[];cross=[];R=11
    for a,b in MARKS:
        W=sorted(set(range(12))-{a,b});J=g.extract(W,W).inv()
        va=-nu[a]*g.extract([a],W).T;vb=-nu[b]*g.extract([b],W).T
        old=s.diag(nu[a],nu[b])*g.extract([a,b],[a,b])*s.diag(nu[a],nu[b])
        old-=s.Matrix.hstack(va,vb).T*J*s.Matrix.hstack(va,vb)
        assert old==s.Matrix([[1,-1],[-1,1]])/R
        for mark,far in((a,vb),(b,va)):
            for p in profiles(mark):
                w=s.Matrix([p['contacts'].get(i,0)for i in W])
                assert s.simplify(sum(nu[i]*q for i,q in p['contacts'].items())-p['h'])==0
                assert s.simplify((far.T*J*w)[0]-p['l']/R)==0
                assert s.simplify(2-(w.T*J*w)[0]-p['c0']-p['l']**2/R)==0
                rows.append(dict(marks=[a,b],mark=mark,profile=p['name'],
                    c0=str(p['c0']),moment=str(p['l']),height=str(p['h']),
                    contacts={str(i):str(q)for i,q in p['contacts'].items()}))
        for p,q in product(profiles(a),profiles(b)):
            w=s.Matrix([p['contacts'].get(i,0)for i in W]);z=s.Matrix([q['contacts'].get(i,0)for i in W])
            value=s.simplify((w.T*J*z)[0]);assert s.simplify(value-p['l']*q['l']/R)==0
            cross.append(dict(marks=[a,b],left=p['name'],right=q['name'],value=str(value)))
    h,k,l,m,c0,d0,p,q,r=s.symbols('h k l m c0 d0 p q r',nonzero=True);A=s.Rational(1,R)
    T=s.Matrix([[A,-A,-h+l*A,-m*A],[-A,A,-l*A,-k+m*A],
        [-h+l*A,-l*A,c0+l*l*A-p,-l*m*A-q],[-m*A,-k+m*A,-l*m*A-q,d0+m*m*A-r]])
    identity=1-(2*l/h+c0/h**2+2*m/k+d0/k**2)/R+(p/h**2-2*q/(h*k)+r/k**2)/R
    assert s.factor(T.det()/(h*h*k*k)-identity)==0
    phi=(1+s.sqrt(5))/2
    eligible=[]
    for name,t in(('isolated',s.sqrt(2)),('isolated',phi),('mate3',1),('mate4',0),('center3',0)):
        z=next(p for p in profiles(8)if p['name']==name)
        hh=s.simplify(z['h'].subs(z['t'],t));cost=s.simplify((2*z['l']/z['h']+z['c0']/z['h']**2).subs(z['t'],t))
        energy=s.simplify(cost-3)
        assert 0<hh<=1 and 0<=energy<=1
        eligible.append(dict(profile=name,opposite_magnitude=str(t),height=str(hh),cost=str(cost),weighted_spectator_energy=str(energy)))
    assert [r['cost']for r in eligible]==['4','12 - 4*sqrt(5)','3','4','4']
    derivatives=[]
    for z in profiles(8):
        t=z['t'];cost=2*z['l']/z['h']+z['c0']/z['h']**2
        derivative=s.factor(s.diff(cost,t));assert derivative.subs(t,1)<0
        derivatives.append(dict(profile=z['name'],cost=str(s.factor(cost)),derivative=str(derivative)))
    assert s.Rational(8,3)<3 and 4*s.sqrt(2)/(1+s.sqrt(2))<3
    assert s.Rational(9,8)>1 and s.Rational(16-R,4)==s.Rational(5,4)
    assert len(rows)==32 and len(cross)==64
    return dict(resistance=R,end_rows=rows,cross_rows=cross,determinant=str(identity),
        native_cost_derivatives=derivatives,nonextreme_eligible=eligible,
        both_nonextreme_cost_bound=8,extreme_cost=8,extreme_height='1/2',
        both_extreme_unscaled_spectator_energy='5/4',
        nonextreme_extreme_side_support='Must vanish: a one-sided contact costs at least2; two sides cost strictly more than9/8, while the required weighted energy is at most1.')

def payment(g,W,p,q):
    ids=sorted(set(W)|{p,q});parts=B.simple_finite_components(g,ids)
    for part in parts:part['marks']=[v for v in(p,q)if v in part['nodes']]
    u=next(c for c in parts if p in c['nodes']);v=next(c for c in parts if q in c['nodes'])
    unit=0;half=F();upgrade=False
    if u is v:assert len(u['nodes'])<=17;unit=1
    else:
        for small,large in((u,v),(v,u)):
            if large['type']=='D11' and any(i in large['short_leaves']for i in large['marks']):
                if small['type']=='A1':unit=1;upgrade=True
                elif small['type'].startswith('A')and len(small['nodes'])<=5:half=F(1,2)
    return dict(marks=[p,q],components=parts,unit=unit,old_half=str(half),full_A1_D11_upgrade=upgrade)

def one_sided_corners():
    f,_=affine();rows=[]
    for(a,b),left,right,contacting in product(MARKS,(8,9),(10,11),(12,13)):
        g=2*s.eye(14);g[:12,:12]=f
        for i,j in((12,left),(13,right)):g[i,j]=g[j,i]=-1
        W=sorted(set(range(12))-{a,b});p,q=(13,a)if contacting==12 else(b,12)
        corner=payment(g,W,p,q);assert corner['unit']==1
        assert corner['full_A1_D11_upgrade']==((right==b)if contacting==12 else(left==a))
        rows.append(dict(marks=[a,b],leaf_attachments=[left,right],excluded_black=contacting,
            selected_corner=corner,scope='Arbitrary spectator contacts occur only at the excluded black. The selected finite components are therefore whole. The excluded black profile is irrelevant to this corner.'))
    assert len(rows)==32
    return rows

def two_components():
    f,_=affine();rows=[];tally=Counter()
    for(a,b),left,right,small_side in product(MARKS,(8,9),(10,11),(12,13)):
        g=2*s.eye(18);g[:12,:12]=f
        other=25-small_side
        for i,j in((12,left),(13,right),(small_side,14),(other,15),(15,16),(16,17)):
            g[i,j]=g[j,i]=-1
        W=sorted(set(range(12))-{a,b})+[14,15,16,17]
        assert C.inertia(g.tolist())==[16,1,1]
        assert C.inertia(g.extract(W,W).tolist())==[14,0,0]
        for root in(a,b,12,13):assert C.inertia(g.extract(W+[root],W+[root]).tolist())==[15,0,0]
        assert C.inertia(g.extract(W+[a,b],W+[a,b]).tolist())==[15,0,1]
        corners=[payment(g,W,p,q)for p,q in((b,12),(12,13),(13,a))]
        units=sum(c['unit']for c in corners);halves=sum((F(c['old_half'])for c in corners),F())
        assert units+halves>=1 and not any(c['full_A1_D11_upgrade']for c in corners)
        tally[(units,str(halves))]+=1
        rows.append(dict(marks=[a,b],leaf_attachments=[left,right],A1_contacting_black=small_side,
            gram=[list(map(str,row))for row in g.tolist()],signature=[16,1,1],white=W,
            corners=corners,unit_total=units,half_total=str(halves)))
    assert len(rows)==32
    return dict(rows=rows,payment_classes=[dict(units=k[0],halves=k[1],cases=v)for k,v in sorted(tally.items())])

def run():
    stems=(B.TAG,M.TAG,SPECTATOR,'fresh_finite8_a17_d11_coupon_budget_20261017')
    for stem in stems:
        old=json.loads((ROOT/f'results/{stem}.json').read_text());assert old['status'].startswith('PASS')
        for name,h in old['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
    f,_=affine();adj={i:{j for j in range(12)if i!=j and f[i,j]}for i in range(12)}
    assert all(C.distances(adj,list(range(12)),a)[b]==9 for a,b in MARKS)
    kept=sorted(set(combinations(range(12),2))-set(MARKS));assert len(kept)==62
    paths=[Path(__file__),Path(C.__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths+=[ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    return dict(status='PASS_UNIVERSAL_AFFINE_D11_CROSS_FORK_LEAF_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        formulas=formulas(),arbitrary_one_sided_corners=one_sided_corners(),two_A1_A3_spectators=two_components(),
        removed_pairs=list(map(list,MARKS)),retained_pairs=list(map(list,kept)),old_cost=66,new_cost=62,
        new_finite_coupon='Upgrade existing half to one on whole A1 versus a short leaf of whole D11; the fork-center half is unchanged.',
        theorem='Every quadrilateral using an opposite short-leaf pair of a whole affine D11 factor is paid by a retained path<=5, a whole finite component<=17, the full A1/D11-short-leaf coupon, or two inherited A2/A4 versus D11-short-leaf halves.',
        scope='Universal local transfer, including selected-core, low-core and immediate-factor uses. A new finite coupon is explicit. Global finite budgets, remaining A/B/C order12 costs and dimension comparisons are separate obligations.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
