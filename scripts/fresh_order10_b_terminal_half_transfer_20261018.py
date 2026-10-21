#!/usr/bin/env python3
"""Order10 B terminal halves with two explicit whole-component coupons.

The report proves the all-native profile and arbitrary-spectator argument.
No selected-core census or forest optimization is run here.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
import fresh_order11_bcd_half_pair_transfer_20261018 as K
import fresh_a15_d15_ambient59_forest_closure_20261018 as C

ROOT=C.ROOT
TAG='fresh_order10_b_terminal_half_transfer_20261018'
CZERO='fresh_order10_c_endpoint_zero_transfer_20261018'
SPECTATOR='fresh_spectator_energy_below_one_20261017'
STEMS=(C.TAG,CZERO,K.TAG,SPECTATOR,'fresh_order12_b_center_terminal_transfer_20261018',
       'fresh_order13_bc_sparse_terminal_transfer_20261018','fresh_a15_d15_dimension58_20261018')
M,R,A,B=10,9,8,7
SQ=s.sqrt(2)


def affine():
    edges=[(i,i+1,SQ if i==6 else 1)for i in range(7)]+[(0,8,1),(0,9,1)]
    g=K.gram(M,edges);nu=s.Matrix([1]*7+[1/SQ,s.Rational(1,2),s.Rational(1,2)])
    assert g*nu==s.zeros(M,1)and g.rank()==R
    return g,nu


def profiles(mark):
    t=s.Symbol('t',nonnegative=True)
    if mark==A:
        positions=(9,0);heights=(t/2,(t+1)/2,(t+SQ)/2,1+t/2);end='fork'
    else:
        assert mark==B
        positions=(6,5);heights=(t/SQ,t/SQ+1,t/SQ+SQ,t/SQ+1);end='heavy'
    return [dict(name=name,contacts={mark:t,**extra},moment=l,c0=c0,height=h,t=t,end=end)
        for name,extra,l,c0,h in zip(('isolated','end3','end4','next3'),
            ({},{positions[0]:1},{positions[0]:SQ},{positions[1]:1}),
            (0,1,SQ,2),(2,1,0,0),heights)]


def type_of_affine(g,marks):
    adj=K.graph(g);heavy=[(i,j)for i,j in combinations(range(g.rows),2)if g[i,j]**2==2]
    branches=[i for i in adj if len(adj[i])==3]
    distance=K.distances(adj,marks[0])[marks[1]]
    assert all(len(adj[i])==1 for i in marks)
    if len(heavy)==2 and not branches:assert distance==9;return 'C'
    if len(heavy)==len(branches)==1:assert distance==8;return 'B'
    assert not heavy and len(branches)==2 and distance==7
    return 'D'


def physical():
    g,nu=affine();W=[i for i in range(M)if i not in(A,B)];J=g.extract(W,W).inv()
    va,vb=[-nu[i]*g.extract([i],W).T for i in(A,B)]
    old=s.diag(nu[A],nu[B])*g.extract([A,B],[A,B])*s.diag(nu[A],nu[B])
    assert (old-s.Matrix.hstack(va,vb).T*J*s.Matrix.hstack(va,vb)).applyfunc(s.simplify)==s.Matrix([[1,-1],[-1,1]])/R
    inverses=[];rows=[];immediate=[];cross=[];u=s.Symbol('u',nonnegative=True)
    for mark,other,far in((A,B,vb),(B,A,va)):
        ids=[i for i in range(M)if i!=mark];deleted=g.extract(ids,ids);pivots=K.pd(deleted)
        inv=deleted.inv();assert all(v.is_nonnegative for v in inv)
        allowed=[(ids[i],inv[i,i])for i in range(R)if inv[i,i]<=2]
        assert sorted(allowed)==([(0,2),(9,1)]if mark==A else[(5,2),(6,1)])
        ii=[ids.index(i)for i,v in allowed];two=s.simplify(sum(inv[i,j]for i in ii for j in ii));assert two>2
        inverses.append(dict(mark=mark,nodes=ids,diagonal=list(map(str,inv.diagonal())),pivots=pivots,two_contact_minimum=str(two)))
        for p in profiles(mark):
            x=s.Matrix([p['contacts'].get(i,0)for i in W]);h,l,c0,t=[p[k]for k in('height','moment','c0','t')]
            assert s.simplify(sum(nu[i]*q for i,q in p['contacts'].items())-h)==0
            assert s.simplify((far.T*J*x)[0]-l/R)==0
            assert s.simplify(2-(x.T*J*x)[0]-c0-l*l/s.Integer(R))==0
            cost=s.factor(2*l/h+c0/h**2)
            minimum=1 if p['name']=='isolated'else 0
            maximum=s.simplify(cost.subs(t,minimum))
            num,den=s.fraction(s.factor((maximum-cost).subs(t,u+minimum)))
            np,dp=s.Poly(s.expand(num),u),s.Poly(s.expand(den),u)
            assert all(q.is_nonnegative for q in np.all_coeffs()+dp.all_coeffs())and dp.eval(0).is_positive
            record=dict(mark=mark,end=p['end'],profile=p['name'],height=str(h),moment=str(l),c0=str(c0),cost=str(cost),minimum_contact=minimum,maximum=str(maximum))
            if p['name']in('isolated','end3'):
                cutoff=SQ if p['name']=='isolated'else 1
                record['nonextreme_maximum']=str(s.simplify(cost.subs(t,cutoff)))
                if mark==A:assert cost.subs(t,cutoff)<=4
                else:assert cost.subs(t,cutoff)<3
            rows.append(record)
            if p['name']in('end4','next3'):
                gg=2*s.eye(M);gg[:R,:R]=deleted
                for i,q in p['contacts'].items():
                    if i!=mark:gg[ids.index(i),R]=gg[R,ids.index(i)]=-q
                kernel=gg.nullspace();assert len(kernel)==1
                v=kernel[0]if kernel[0][0].is_positive else-kernel[0]
                assert all(z.is_positive for z in v)
                kind=type_of_affine(gg,(R,ids.index(other)))
                expected={('fork','end4'):'C',('fork','next3'):'B',('heavy','end4'):'B',('heavy','next3'):'D'}[p['end'],p['name']]
                assert kind==expected
                immediate.append(dict(mark=mark,profile=p['name'],target=kind,kernel=list(map(str,v)),payment='independent_C10_first'if kind=='C'else'1/2'if kind=='B'else'1'))
    for p,q in product(profiles(A),profiles(B)):
        x=s.Matrix([p['contacts'].get(i,0)for i in W]);y=s.Matrix([q['contacts'].get(i,0)for i in W])
        v=s.simplify((x.T*J*y)[0]);assert s.simplify(v-p['moment']*q['moment']/R)==0
        cross.append(dict(left=p['name'],right=q['name'],value=str(v)))
    assert len(rows)==8 and len(cross)==16 and len(immediate)==4
    return dict(gram=[[str(x)for x in row]for row in g.tolist()],positive_kernel=list(map(str,nu)),inverses=inverses,profiles=rows,cross=cross,immediate=immediate,
        resistance=R,extreme_fork_height='1/2',extreme_fork_cost='8',nonextreme_pair_gap='1/9',spectator_energy_maximum='3')


def finite_corners():
    f,nu=affine();rows=[]
    for attachment,r in product((A,9),range(4)):
        c=10;n=11+r;g=2*s.eye(n);g[:M,:M]=f;g[c,attachment]=g[attachment,c]=-1
        if r:
            for i,j in [(c,11)]+[(i,i+1)for i in range(11,10+r)]:g[i,j]=g[j,i]=-1
        ids=[i for i in range(n)if i!=A];corner=g.extract(ids,ids);K.pd(corner)
        parts=K.components(g,ids);adj=K.graph(g)
        if attachment==A:
            cp=next(z for z in parts if c in z);bp=next(z for z in parts if B in z)
            assert len(cp)==r+1 and len(bp)==9 and len(parts)==2
            assert len(adj[B]&set(bp))==1 and g[B,next(iter(adj[B]&set(bp)))]==-SQ
            assert len(adj[c]&set(cp))==(0 if r==0 else 1)
            payment='1/2';role=f'A{r+1}_endpoint/B9_doubled_terminal'
        else:
            assert len(parts)==1 and len(parts[0])==10+r
            assert K.distances(K.graph(corner),ids.index(c))[ids.index(B)]==9
            assert len(adj[B]&set(ids))==1
            payment='1';role=f'whole_B{10+r}_distance9_doubled_terminal'
        rows.append(dict(extreme_attachment=attachment,spectator_path_rank=r,whole_components=parts,marked_nodes=[c,B],role=role,payment=payment,total_with_original_half=str(F(1,2)+F(payment))))
    return rows


def actual_cross(a,b):
    pair=frozenset((a,b))
    return C.actual_cross(a,b)+(8 if pair==frozenset((('A',1),('B',9)))else 16 if pair==frozenset((('A',3),('B',9)))else 0)


def upper_cross(a,b):
    pair=frozenset((a,b))
    return C.upper_cross(a,b)+(8 if pair==frozenset((('A',1),('B',9)))else 16 if pair==frozenset((('A',3),('B',9)))else 0)


def run():
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text())for tag in STEMS}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS')
        for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,(tag,path)
    assert docs[CZERO]['new_cost']=='44'and docs[CZERO]['strict_minimum']=='1/9'
    ph=physical();fc=finite_corners()
    g,nu=affine();adj=K.graph(g);old={pair:F(1)for pair in combinations(range(M),2)};new=dict(old)
    for pair in((7,8),(7,9)):assert K.distances(adj,pair[0])[pair[1]]==8;new[pair]=F(1,2)
    assert sum(new.values(),F())==44
    for i,j in new:
        if K.distances(adj,i)[j]<=7:assert new[i,j]==1
    swap={i:(9 if i==8 else 8 if i==9 else i)for i in range(M)}
    assert all(new[tuple(sorted((swap[i],swap[j])))]==w for(i,j),w in new.items())
    # One-sided energy <=3/4 permits exactly whole A1,A2,A3 endpoints.
    energies=[F(r,r+1)for r in range(1,4)];assert energies==[F(1,2),F(2,3),F(3,4)]
    assert 2*F(1,2)>F(3,4)  # two active components or one nonsimple contact
    assert 4*F(1,2)==2 and F(4+2,2)==3
    # Every pair in these small finite B factors remains a unit.
    joined_costs={str(n):str(F(C.actual_base(('B',n)),16))for n in range(10,14)}
    assert all(F(joined_costs[str(n)])==F(n*(n-1),2)for n in range(10,14))
    engine=C.C.C.P.Z.B.B;types=engine.types();rep=engine.representative
    for a in types:
        for b in types:assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    costs={str(r):dict(previous=str(F(C.actual_cross(('A',r),('B',9)),16)),new=str(F(actual_cross(('A',r),('B',9)),16)))for r in range(1,5)}
    assert [costs[str(r)]['previous']for r in range(1,5)]==['0','1','0','1']
    assert [costs[str(r)]['new']for r in range(1,5)]==['1/2','1','1','1']
    paths={Path(__file__),Path(K.__file__),Path(C.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}'for tag in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths.update(ROOT/f'results/{tag}.md'for tag in('fresh_cusp_factor_payment_framework_20261017','fresh_finite7_ideal5_triangle_payments_20261017'))
    return dict(status='PASS_ORDER10_B_TERMINAL_HALF_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        physical=ph,finite_corners=fc,changed_pairs=[[7,8],[7,9]],old_weight='1',new_weight='1/2',old_cost='45',new_cost='44',
        old_weighted_pairs=[dict(pair=list(p),weight=str(w))for p,w in old.items()],weighted_pairs=[dict(pair=list(p),weight=str(w))for p,w in new.items()],
        new_coupons=[dict(small='A1',small_role='unique vertex',large='B9',large_role='doubled terminal',weight='1/2',component_pair_increment='1/2'),
                     dict(small='A3',small_role='endpoint',large='B9',large_role='doubled terminal',weight='1/2',component_pair_increment='1')],
        component_pair_costs=costs,retained_joined_B_factor_costs=joined_costs,simultaneous_type_pairs=len(types)**2,new_order10_costs=dict(A='45',B='44',C='44',D='45'),order10_envelope_unchanged='9/2',
        priority='Independent actual C10 endpoint corner first; then actual B10 changed pair with its retained half and a distinct half/unit; otherwise pointwise completed58 rule.',
        triangle_minimum=2,quadrilateral_minimum=1,new_forest_count=False,new_dimension_bound=None,
        scope='Universal local half transfer with all native profiles and arbitrary supporting spectators; two ordinary finite coupons explicitly charged. No new forest bound or dimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    result=json.loads(json.dumps(run()))
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['new_cost'],result['component_pair_costs'],flush=True)
