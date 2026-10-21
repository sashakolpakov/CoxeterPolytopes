#!/usr/bin/env python3
"""Universal order11 B/C zero transfers with one explicit coupon upgrade.

Exact physical profile identities and distinct whole-corner payments.
No geometric census, forest optimization, or ambient theorem.
"""
import argparse,json
from hashlib import sha256
from pathlib import Path
from itertools import product
import sympy as s
import fresh_order11_bcd_half_pair_transfer_20261018 as K
import fresh_finite_rank15_rank16_rank62_budget_20261018 as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order11_bc_zero_transfer_20261018'
STEMS=(K.TAG,'fresh_order11_c_endpoint_zero_transfer_20261018',
 'fresh_order16_d_distance9_half_transfer_20261018',
 'fresh_order16_d105_ambient61_frontier_20261018',
 'fresh_rank15_rank16_dimension61_20261018',B.TAG,
 'fresh_finite7_ideal5_triangle_payments_20261017')


def enc(x):
    if isinstance(x,s.Basic):
        assert not x.has(s.Float)
        return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def positive_polynomial(expr,variables):
    num,den=s.fraction(s.factor(expr))
    assert all(v.is_nonnegative for v in s.Poly(s.expand(num),*variables).coeffs())
    assert all(v.is_nonnegative for v in s.Poly(s.expand(den),*variables).coeffs())
    assert den.subs({v:0 for v in variables}).is_positive
    return dict(numerator=s.expand(num),denominator=s.expand(den))


def c_profiles():
    f,nu,_=K.affine('C');a,b=0,9
    W=[i for i in range(11) if i not in (a,b)]
    J=f.extract(W,W).inv().applyfunc(s.simplify)
    va,vb=[-nu[i]*f.extract([i],W).T for i in (a,b)]
    old=s.diag(nu[a],nu[b])*f.extract([a,b],[a,b])*s.diag(nu[a],nu[b])
    old-=s.Matrix.hstack(va,vb).T*J*s.Matrix.hstack(va,vb)
    assert old.applyfunc(s.simplify)==s.Matrix([[1,-1],[-1,1]])/9
    blocks=[list(range(9)),[10]]
    inv=f.extract(blocks[0],blocks[0]).inv().applyfunc(s.simplify)
    assert inv[8,8]==1 and inv[7,7]==2 and all(inv[i,i]>2 for i in range(7))
    t,q,u=s.symbols('t q u',nonnegative=True)
    heavy=[]
    for p in K.profiles('C',a,nu):
        p=dict(p);p['contacts']={i:s.sympify(v).subs(p['variable'],t) for i,v in p['contacts'].items()}
        p['height']=p['height'].subs(p['variable'],t);p['variable']=t
        p['maximum']={'isolated':4,'end3':3,'end4':2,'next3':4}[p['name']]
        p['lower_t']=1 if p['name']=='isolated' else 0
        heavy.append(p)
    inner=[
      dict(name='isolated',contacts={b:t},height=t,moment=0,c0=2,maximum=2,lower_t=1,tail_contact=0),
      dict(name='tail_only',contacts={b:t,10:q},height=t+q/s.sqrt(2),moment=0,c0=2-q*q/2,maximum=3,lower_t=0,tail_contact=q),
      dict(name='end3',contacts={b:t,8:1,10:q},height=t+1+q/s.sqrt(2),moment=1,c0=1-q*q/2,maximum=3,lower_t=0,tail_contact=q),
      dict(name='end4',contacts={b:t,8:s.sqrt(2)},height=t+s.sqrt(2),moment=s.sqrt(2),c0=0,maximum=2,lower_t=0,tail_contact=0),
      dict(name='next3',contacts={b:t,7:1},height=t+1,moment=2,c0=0,maximum=4,lower_t=0,tail_contact=0)]
    identities=[];bounds=[]
    for side,rows,far in [('heavy',heavy,vb),('inner',inner,va)]:
        for p in rows:
            v=s.Matrix([p['contacts'].get(i,0) for i in W]);h=p['height'];l=s.sympify(p['moment']);c=p['c0']
            assert s.simplify(sum(nu[i]*v for i,v in p['contacts'].items())-h)==0
            assert s.simplify((far.T*J*v)[0]-l/9)==0
            assert s.simplify(2-(v.T*J*v)[0]-c-l*l/9)==0
            cost=s.factor(2*l/h+c/h**2)
            # Multiplication by positive h² avoids an unnecessary rational denominator.
            gap=s.expand(p['maximum']*h*h-2*l*h-c)
            gap=gap.subs(t,u+p['lower_t'])
            if p['name']=='tail_only':gap=gap.subs(q,q+1)
            cert=positive_polynomial(gap,(u,q))
            identities.append(dict(side=side,name=p['name'],height=h,moment=l,c0=c,cost=cost,
                                   contacts=p['contacts'],maximum=p['maximum']))
            bounds.append(dict(side=side,name=p['name'],height_squared_gap=cert,
                domain='t>=1' if p['lower_t'] else 't>=0; tail_only1<=q<2, end3 0<=q<=sqrt2, otherq0'))
    cross=[]
    for p,r in product(heavy,inner):
        x=s.Matrix([p['contacts'].get(i,0) for i in W]);y=s.Matrix([r['contacts'].get(i,0) for i in W])
        value=s.simplify((x.T*J*y)[0]);assert s.simplify(value-s.sympify(p['moment'])*r['moment']/9)==0
        cross.append(dict(heavy=p['name'],inner=r['name'],product=value,
                          strict_gap=1-s.Rational(p['maximum']+r['maximum'],9)))
    assert len(cross)==20 and min(r['strict_gap'] for r in cross)==s.Rational(1,9)
    return dict(gram=f.tolist(),kernel=list(nu),marks=[a,b],white=W,resistance=9,
        deleted_inner_components=blocks,B9_inverse=inv.tolist(),profiles=identities,bounds=bounds,
        cross=cross,uniform_strict_gap=s.Rational(1,9))


def b_equality_corners():
    f,nu,(a,b)=K.affine('B');assert (a,b)==(8,9)
    rows=[]
    for profile,attach in product(('isolated','next3'),(9,10)):
        g=2*s.eye(14);g[:11,:11]=f
        p=next(p for p in K.profiles('B',a,nu) if p['name']==profile)
        t=1 if profile=='isolated' else 0
        for i,v in p['contacts'].items():g[11,i]=g[i,11]=-s.sympify(v).subs(p['variable'],t)
        g[12,attach]=g[attach,12]=-1;g[12,13]=g[13,12]=-1
        white=[i for i in range(14) if i not in (a,b,11,12)]
        K.pd(g.extract(white,white))
        for mark in (a,b,11,12):K.pd(g.extract(sorted(white+[mark]),sorted(white+[mark])))
        first=sorted(white+[a,12]);K.pd(g.extract(first,first));parts=K.components(g,first)
        if attach==b:
            assert sorted(map(len,parts))==[2,10]
            assert next(v for v in parts if 12 in v)==[12,13]
            far=next(v for v in parts if a in v);assert len(far)==10
            adj=K.graph(g.extract(far,far));j=far.index(a)
            assert len(adj[j])==1 and g[a,far[next(iter(adj[j]))]]**2==2
            first_payment='1/2';first_reason='wholeA2/B10 doubled-terminal existinghalf'
        else:
            assert len(parts)==1 and len(parts[0])==12
            first_payment='1';first_reason='whole finite B12 unit'
        second=sorted(white+[b,11]);parts2=K.components(g,second)
        if profile=='isolated':
            K.pd(g.extract(second,second));assert sorted(map(len,parts2))==[1,1,10]
            assert [11] in parts2 and [13] in parts2
            far=next(v for v in parts2 if b in v);adj=K.graph(g.extract(far,far));j=far.index(b)
            assert len(far)==10 and len(adj[j])==1 and len(adj[next(iter(adj[j]))])==3
            second_reason='wholeA1/D10 short-fork-leaf existinghalf'
        else:
            assert sorted(map(len,parts2))==[1,11] and [13] in parts2
            affine=next(v for v in parts2 if b in v);h=g.extract(affine,affine)
            null=h.nullspace();assert len(null)==1
            z=null[0] if null[0][0]>0 else -null[0];assert all(v>0 for v in z)
            assert K.marked_half_type(h,(affine.index(b),affine.index(11)))=='D'
            second_reason='whole affineD10 opposite-fork retainedhalf'
        rows.append(dict(heavy_profile=profile,heavy_t=t,fork_attachment=attach,gram=g.tolist(),white=white,
            first_corner=[a,12],first_components=parts,first_payment=first_payment,first_reason=first_reason,
            second_corner=[b,11],second_components=parts2,second_payment='1/2',second_reason=second_reason,
            total=str(s.Rational(first_payment)+s.Rational(1,2)),original_B_payment='0'))
    assert len(rows)==4
    return rows


def b_profiles():
    f,nu,(a,b)=K.affine('B');W=[i for i in range(11) if i not in (a,b)]
    J=f.extract(W,W).inv().applyfunc(s.simplify)
    va,vb=[-nu[i]*f.extract([i],W).T for i in (a,b)]
    rows=[];u=s.Symbol('u',nonnegative=True)
    for mark,far in ((a,vb),(b,va)):
        for p in K.profiles('B',mark,nu):
            x=s.Matrix([p['contacts'].get(i,0) for i in W])
            assert s.simplify((far.T*J*x)[0]-s.sympify(p['moment'])/10)==0
            assert s.simplify(2-(x.T*J*x)[0]-p['c0']-s.sympify(p['moment'])**2/10)==0
            h,l,c,t=p['height'],p['moment'],p['c0'],p['variable']
            if mark==a:
                maximum=dict(isolated=4,end3=3,end4=2,next3=4)[p['name']]
            else:maximum=8 if p['name'] in ('isolated','end3') else 4
            minimum=1 if p['name']=='isolated' else 0
            gap=s.expand((maximum*h*h-2*l*h-c).subs(t,u+minimum))
            cert=positive_polynomial(gap,(u,))
            extra=None
            if mark==b and p['name'] in ('isolated','end3'):
                nonextreme=s.sqrt(2) if p['name']=='isolated' else 1
                upper=4 if p['name']=='isolated' else 3
                extra=dict(lower_t=nonextreme,maximum=upper,
                    certificate=positive_polynomial(s.expand((upper*h*h-2*l*h-c).subs(t,u+nonextreme)),(u,)))
                assert h.subs(t,minimum)==s.Rational(1,2)
            rows.append(dict(end=p['end'],profile=p['name'],height=h,moment=l,c0=c,
                maximum=maximum,minimum_t=minimum,certificate=cert,nonextreme=extra))
    assert len(rows)==8
    return rows


def determinant(R):
    h,k,l,m,c,d,p,q,r=s.symbols('h k l m c d p q r',nonzero=True)
    a=s.Rational(1,R)
    t=s.Matrix([[a,-a,-h+l*a,-m*a],[-a,a,-l*a,-k+m*a],
      [-h+l*a,-l*a,c+l*l*a-p,-l*m*a-q],[-m*a,-k+m*a,-l*m*a-q,d+m*m*a-r]])
    formula=1-(2*l/h+c/h**2+2*m/k+d/k**2)/R+(p/h**2-2*q/(h*k)+r/k**2)/R
    assert s.factor(t.det()/(h*h*k*k)-formula)==0
    return formula


def run():
    docs={t:json.loads((ROOT/f'results/{t}.json').read_text()) for t in STEMS}
    for name,d in docs.items():
        assert d['status'].startswith('PASS'),name
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(name,p)
    c=c_profiles();b=b_equality_corners();bp=b_profiles()
    # Only the A1/B10 marked coupon increases: doubled total1 ->2.
    assert B.cross(('A',1),('B',10))==1
    assert B.cross(('A',1),('D',10))==2
    assert B.cross(('A',2),('B',10))==2
    old_c=docs['fresh_order11_c_endpoint_zero_transfer_20261018']
    assert old_c['new_order11_costs']==dict(A='99/2',B='54',C='54',D='53')
    old_b=next(r for r in docs[K.TAG]['factor_rules'] if r['kind']=='B')
    rules=[]
    for kind,zero,total in [('B',[(8,9),(8,10)],53),('C',[(0,9),(1,10),(0,10)],52)]:
        f,_,_=K.affine(kind);adj=K.graph(f)
        old={tuple(r['pair']):s.Rational(r['weight']) for r in (old_b if kind=='B' else old_c)['weighted_pairs']}
        assert len(old)==55 and sum(old.values())==54
        new={p:(s.Integer(0) if p in zero else w) for p,w in old.items()}
        assert sum(new.values())==total and all(w==1 for p,w in new.items() if p not in zero)
        assert all(new[p]==old[p] for p in old if K.distances(adj,p[0])[p[1]]<=8)
        aut=list(range(11)) if kind=='B' else list(reversed(range(11)))
        if kind=='B':aut[9],aut[10]=10,9
        assert all(f[i,j]==f[aut[i],aut[j]] for i in range(11) for j in range(11))
        assert {tuple(sorted(aut[i] for i in p)) for p in zero}==set(zero)
        assert all(K.distances(adj,i)[j]>=9 for i,j in zero)
        assert 55-len(zero)==total
        rules.append(dict(kind=kind,zero_pairs=zero,total_cost=total,
                          weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(new.items())],
                          all_other_pairs_weight='1',changed_distance_minimum=9))
    paths=[Path(__file__),Path(K.__file__),Path(B.__file__)]+[ROOT/f'{folder}/{t}.{ext}' for t in STEMS
        for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return enc(dict(status='PASS_ORDER11_BC_ZERO_WITH_A1_B10_UPGRADE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        C_near_terminal=c,B_complete_profiles=bp,B_equality_corners=b,
        determinants={9:determinant(9),10:determinant(10)},
        rules=rules,new_family_totals=dict(A='99/2',B='53',C='52',D='53'),new_maximum='53',
        new_normalized_rank10_envelope='53/11',
        coupon_upgrade=dict(small='whole A1',large='whole B10',large_role='doubled terminal',
            old_marked_weight='1/2',new_marked_weight='1',added_component_pair_cost='1/2'),
        existing_representative_comparison=dict(actual_A1_B10_new_doubled=2,existing_A1_D10_doubled=2,
            other_pairs_unchanged=True,new_optimization_performed=False),
        triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
        scope='Universal complete all-native principal profiles for C near-terminal pairs and B terminal/fork pairs, with arbitrary spectators and actual-corner priorities. One explicit A1/B10 half-to-unit coupon upgrade; no new finite optimization or global dimension conclusion.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_family_totals'],flush=True)
