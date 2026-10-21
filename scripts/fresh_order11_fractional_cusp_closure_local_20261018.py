#!/usr/bin/env python3
"""Exact local identities and all ideal margins for fractional order11 rule."""
import argparse
from fractions import Fraction as F
import fcntl
from itertools import product
import json
from pathlib import Path
import signal
import sympy as s
import fresh_order12_b_center_terminal_transfer_20261018 as B
import fresh_order11_bc_zero_transfer_20261018 as Z
import fresh_a16_rank17_spectator_extension_scout_20261018 as R

ROOT=R.ROOT;K=B.L;U=R.U;TAG='fresh_order11_fractional_cusp_closure_local_20261018'
DELTA=F(11,512)


def enc(x):
    if isinstance(x,s.Basic):assert not x.has(s.Float);return str(x)
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(tuple,list)):return[enc(v)for v in x]
    return x


def increment(fs):
    """Exact ordinary increment, in scale512."""
    from collections import Counter
    co=Counter(fs);a1,a2,a3,a4=(co['A',i]for i in range(1,5))
    return(22*a1+44*(a2+a3+a4))*co['D',10]+22*a2*co['B',10]+(88*a2+44*a3)*co['B',8]


def bprofiles():
    t,u,v=s.symbols('t u v',nonnegative=True)
    center=[dict(name='outer',contacts={0:t,9:u,10:v},c0=2-(u*u+v*v)/2,moment=0,height=t+(u+v)/2),
        dict(name='inner3',contacts={0:t,1:1,9:u,10:v},c0=1-(u*u+v*v)/2,moment=1,height=t+1+(u+v)/2),
        dict(name='inner4',contacts={0:t,1:s.sqrt(2)},c0=0,moment=s.sqrt(2),height=t+s.sqrt(2)),
        dict(name='inner_next3',contacts={0:t,2:1},c0=0,moment=2,height=t+1)]
    terminal=[dict(name='isolated',contacts={8:t},c0=2,moment=0,height=t/s.sqrt(2)),
        dict(name='end3',contacts={8:t,7:1},c0=1,moment=1,height=t/s.sqrt(2)+1),
        dict(name='end4',contacts={8:t,7:s.sqrt(2)},c0=0,moment=s.sqrt(2),height=t/s.sqrt(2)+s.sqrt(2)),
        dict(name='next3',contacts={8:t,6:1},c0=0,moment=2,height=t/s.sqrt(2)+1)]
    return center,terminal,(t,u,v)


def bformulas():
    g,nu,_=K.affine('B');marks=[0,8];white=[i for i in range(11)if i not in marks];J=g.extract(white,white).inv()
    va,vb=[-nu[i]*g.extract([i],white).T for i in marks]
    M=s.diag(nu[0],nu[8])*g.extract(marks,marks)*s.diag(nu[0],nu[8])-s.Matrix.hstack(va,vb).T*J*s.Matrix.hstack(va,vb)
    assert M.applyfunc(s.simplify)==s.Matrix([[1,-1],[-1,1]])/8
    inverse=[]
    for ns,end,nxt in((list(range(1,9)),1,2),(list(range(8))+[9,10],7,6)):
        iv=g.extract(ns,ns).inv();di=dict(zip(ns,iv.diagonal()))
        assert di[end]==1 and di[nxt]==2 and all(x>2 for i,x in di.items()if i not in(end,nxt))
        inverse.append(dict(nodes=ns,endpoint=end,next_node=nxt,diagonal=di))
    center,terminal,variables=bprofiles();rows=[]
    for end,ps,far in(('center',center,vb),('terminal',terminal,va)):
        for p in ps:
            w=s.Matrix([p['contacts'].get(i,0)for i in white]);h=p['height'];l=p['moment'];c=p['c0']
            assert s.simplify(sum(nu[i]*x for i,x in p['contacts'].items())-h)==0
            assert s.simplify(8*(far.T*J*w)[0]-l)==0
            assert s.simplify(2-(w.T*J*w)[0]-c-l*l/8)==0
            rows.append(dict(end=end,profile=p['name'],contacts=p['contacts'],height=h,moment=l,c0=c,cost=s.factor(2*l/h+c/h**2)))
    cross=[]
    for p,q in product(center,terminal):
        x=s.Matrix([p['contacts'].get(i,0)for i in white]);y=s.Matrix([q['contacts'].get(i,0)for i in white])
        assert s.simplify((x.T*J*y)[0]-p['moment']*q['moment']/8)==0
        cross.append(dict(center=p['name'],terminal=q['name'],product=s.factor((x.T*J*y)[0])))
    t,u,v=variables;cc=[s.factor(2*p['moment']/p['height']+p['c0']/p['height']**2)for p in center]
    assert cc[0].subs({t:0,u:1,v:0})==6 and cc[0].subs({t:0,u:s.sqrt(2),v:0})==2
    assert cc[0].subs({t:0,u:1,v:1})==1
    assert cc[1].subs({t:0,u:0,v:0})==3 and cc[2].subs(t,0)==2 and cc[3].subs(t,0)==4
    w,z=s.symbols('w z',nonnegative=True);bounds=[]
    for p,bound,lower in zip(terminal,(4,3,2,4),(1,0,0,0)):
        h,l,c=p['height'],p['moment'],p['c0'];poly=s.expand((bound*h*h-2*l*h-c).subs(t,w+lower))
        bounds.append(dict(end='terminal',profile=p['name'],maximum=bound,substitution_t=w+lower,
            polynomial=Z.positive_polynomial(poly,(w,))))
    for p,bound in zip(center[1:],(3,2,4)):
        h,l,c=p['height'],p['moment'],p['c0']
        bounds.append(dict(end='center',profile=p['name'],maximum=bound,
            polynomial=Z.positive_polynomial(s.expand(bound*h*h-2*l*h-c),(t,u,v))))
    h,c=center[0]['height'],center[0]['c0'];outer=2*h*h-c
    for name,subs,vs in(('opposite_nonzero',{t:1+w},(w,u,v)),
        ('outer_native_gap_left',{u:s.sqrt(2)+w},(t,w,v)),
        ('outer_native_gap_right',{v:s.sqrt(2)+w},(t,u,w)),
        ('two_outer_contacts',{u:1+w,v:1+z},(t,w,z))):
        bounds.append(dict(end='center',profile='outer',branch=name,maximum=2,
            polynomial=Z.positive_polynomial(s.expand(outer.subs(subs)),vs)))
    return dict(gram=g.tolist(),kernel=list(nu),marks=marks,white=white,resistance=8,old_mark_schur=M.tolist(),
        finite_inverse_diagonals=inverse,profile_identities=rows,cross_identities=cross,
        determinant=Z.determinant(8),center_extreme_cost=6,nonextreme_maximum=4,terminal_maximum=4,
        nonextreme_forces_costs=[4,4],nonextreme_spectator_energy=0,
        extreme_energy_interval=[0,2],native_outer_nonextreme_maximum=2,universal_profile_bounds=bounds)


def bcorner_templates():
    f,_,_=K.affine('B');rows=[]
    # Nonextreme center next3 gives whole affine B8/order9, plus two A1.
    g=2*s.eye(12);g[:11,:11]=f;g[11,2]=g[2,11]=-1
    rat,gauge=B.rational_gram(g,{8});ns=list(range(1,9))+[11];a=rat.extract(ns,ns);z=a.nullspace()[0]
    if z[0]<0:z=-z
    assert all(x>0 for x in z)and K.inertia(a)==[8,0,1]
    adj=K.graph(a);assert sorted(len(adj[i])for i in adj).count(3)==1
    assert K.distances(adj,ns.index(11))[ns.index(8)]==7
    rows.append(dict(branch='nonextreme',gram=rat.tolist(),gauge=gauge,corner_nodes=list(range(1,12)),marks=[11,8],
        whole_affine_nodes=ns,positive_kernel=list(z),affine_type='B8_order9',marked_distance=7,payment='1'))
    for leaf,has1 in product((9,10),(False,True)):
        size=12+int(has1);g=2*s.eye(size);g[:11,:11]=f;g[11,leaf]=g[leaf,11]=-1
        if has1:g[11,12]=g[12,11]=-1
        rat,gauge=B.rational_gram(g,{8});ids=list(range(1,size));K.pd(rat.extract(ids,ids))
        bs=B.finite_components(rat,ids,gauge);small=next(b for b in bs if 11 in b['nodes']);large=next(b for b in bs if 8 in b['nodes'])
        assert small['type']==('A3'if has1 else'A2')and large['type']=='B8'
        assert large['roles']['8']=='doubled_terminal'
        if has1:assert len([j for j in small['nodes']if j!=11 and rat[11,j]])==2
        else:assert small['roles']['11']=='endpoint'
        rows.append(dict(branch='extreme_one_A1'if has1 else'extreme_empty',fork_leaf=leaf,gram=rat.tolist(),gauge=gauge,
            corner_nodes=ids,marks=[11,8],components=bs,new_coupon=4*DELTA,selected_ideal=1-4*DELTA,total=F(1)))
    return rows


def dnext_templates():
    f,nu,(a,b)=K.affine('D');assert(a,b)==(7,9);rows=[]
    t=s.symbols('t',nonnegative=True);gnext=4/(1+t/2)
    assert s.simplify(gnext.subs(t,2)-2)==0 and gnext.subs(t,0)==4
    assert s.simplify((gnext-2).subs(t,2))==0 and(gnext-2).subs(t,0)==2
    for has1,attach in product((False,True),(7,8)):
        tn=0 if has1 else 2;size=14 if has1 else 13;g=2*s.eye(size);g[:11,:11]=f
        g[11,attach]=g[attach,11]=-1;g[12,6]=g[6,12]=-1;g[12,b]=g[b,12]=-tn
        if has1:g[11,13]=g[13,11]=-1
        white=[i for i in range(size)if i not in(a,b,11,12)]
        assert K.inertia(g)==[size-2,1,1];K.pd(g.extract(white,white))
        for root in(a,b,11,12):K.pd(g.extract(sorted(white+[root]),sorted(white+[root])))
        corner_signatures=[]
        for marks in((a,b),(b,11),(11,12),(12,a)):
            ns=sorted(white+list(marks));sig=K.inertia(g.extract(ns,ns));assert sig[1]==0
            corner_signatures.append(dict(marks=marks,signature=sig))
        ids=sorted(white+[11,b]);bs=B.finite_components(g,ids,[1]*size)
        small=next(c for c in bs if 11 in c['nodes']);large=next(c for c in bs if b in c['nodes'])
        if attach==a:
            assert small['type']==('A2'if has1 else'A1')and large['type']=='D10'
            assert large['roles'][str(b)]=='short_leaf'and small['roles']['11']=='endpoint'
            pay=F(1,2)+DELTA;reason='whole_A2_D10'if has1 else'whole_A1_D10'
        else:
            assert small is large and large['type']==('D12'if has1 else'D11')
            adj=K.graph(g.extract(large['nodes'],large['nodes']));dist=K.distances(adj,large['nodes'].index(11))[large['nodes'].index(b)]
            assert dist==9 and large['roles'][str(b)]=='short_leaf'
            pay=F(1);reason='retained_whole_D12_leaf9'if has1 else'retained_whole_D11_leaf9'
        rows.append(dict(next_opposite_product=tn,extreme_attachment=attach,whole_A1_spectator=has1,gram=g.tolist(),white=white,
            next_height=1+s.Rational(tn,2),next_cost=s.simplify(gnext.subs(t,tn)),unscaled_extreme_spectator_energy=F(1,2)if has1 else F(),
            corners=corner_signatures,paying_marks=[11,b],components=bs,payment=pay,reason=reason,selected_ideal=F(1,2)-DELTA,total=pay+F(1,2)-DELTA))
    return dict(templates=rows,next_height=1+t/2,next_cost=gnext,weighted_extreme_energy=gnext-2,
        both_next_cost_sum_upper=8,resistance=10,both_next_determinant_gap='1/5',
        determinant=Z.determinant(10),all_opposite_products_covered=True)


def ideal():
    old=R.V.load(R.H.TAG)['ideal'];rows=[]
    for r in old['rows']:
        saving=F(r['ranks'].count(10),128);margin=F(r['normalized_margin'])+saving
        assert margin>0
        rows.append(dict(ranks=r['ranks'],old_margin=r['normalized_margin'],normalized_saving=saving,normalized_margin=margin))
    return dict(complete_profiles=len(rows),rows=rows,old_negative_count=5,new_negative_count=0,zero_count=0,
        normalized_minimum=min(r['normalized_margin']for r in rows),all_margins_strictly_positive=True,
        rank10_envelope_old='53/11',rank10_envelope_new=F(53,11)-F(1,128),normalized_saving_per_rank10='1/128')


def b_zero_repair():
    rows=[]
    for r in Z.b_equality_corners():
        first=F(1)if r['fork_attachment']==10 else F(1,2)+DELTA
        second=F(1,2)+DELTA if r['heavy_profile']=='isolated'else F(1,2)-DELTA
        assert first+second>=1
        rat,gauge=B.rational_gram(s.Matrix(r['gram']),{8,11}if r['heavy_profile']=='isolated'else{8})
        rows.append(dict(heavy_profile=r['heavy_profile'],fork_attachment=r['fork_attachment'],gram=rat.tolist(),gauge=gauge,
            white=r['white'],first_corner=r['first_corner'],second_corner=r['second_corner'],
            first_payment=first,second_payment=second,total=first+second,
            first_reason=r['first_reason'],second_reason=r['second_reason']))
    return rows


def run():
    b=bformulas();bc=bcorner_templates();d=dnext_templates();bz=b_zero_repair();ii=ideal()
    totals=dict(A=F(99,2)+11*DELTA,B=53-4*DELTA,C=F(52),D=53-4*DELTA)
    assert max(totals.values())==F(53)-4*DELTA and DELTA<=F(7,30)
    paths={Path(__file__),Path(B.__file__),Path(K.__file__),Path(Z.__file__),Path(R.__file__),
        ROOT/'results/order11_fractional_cusp_closure_research_20261018.md',ROOT/f'results/{R.H.TAG}.json'}
    for mod in(B,K,Z):paths.update(ROOT/f'results/{mod.TAG}.{ext}'for ext in('md','json'))
    return enc(dict(status='PASS_ORDER11_FRACTIONAL_CUSP_CLOSURE_LOCAL',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(paths)},
        delta=DELTA,order11_totals=totals,B_center=b,B_corner_templates=bc,D_next=d,B_zero_repair=bz,ideal=ii,
        finite_increment_scale=512,finite_increment='(22*a1+44*(a2+a3+a4))*d10+22*a2*b10+(88*a2+44*a3)*b8',
        arbitrary_spectator_proof_in_bound_MD=True,complete_finite_positive_family_budget=False,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);R.N.A.R.M.R.resources();U.memory_guard()
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();R.N.A.R.M.R.resources(start=False)
    if args.verify:assert R.V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['ideal']['normalized_minimum'],len(out['B_corner_templates']),len(out['D_next']['templates']),flush=True)
