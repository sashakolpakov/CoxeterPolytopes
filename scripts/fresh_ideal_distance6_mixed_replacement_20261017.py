#!/usr/bin/env python3
"""Mixed finite9/ideal6 replacements and a sealed cycle14 equality case."""
import argparse,json
from collections import Counter,deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
import fresh_distance8_ideal_corner_replacement_20261017 as B
import fresh_cycle16_half_coupon_transfer_20261017 as C

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ideal_distance6_mixed_replacement_20261017'
A=B.A


def large_cycles():
    profiles=[p for p in B.cycle_profiles() if p['name'] in
              ('isolated','both_ends3','end3_1','end4_1','next3_1')]
    rows=[]
    for c,d in product(profiles,repeat=2):
        Q=11 if any(p['name']=='next3_1'for p in(c,d))else 10
        for ub,vb in product((0,1),repeat=2):
            z=A.value(A.Q(F(1,7)),A.Q(F(1,Q)),c,d,
                      c['umax']if ub else A.Q(),d['umax']if vb else A.Q())
            assert(z-A.Q(F(1,35))).sign()>=0
            rows.append(dict(c=c['name'],d=d['name'],P=7,Q_lower=Q,
                             upper_reciprocals=[ub,vb],value=z.wire()))
    assert len(rows)==100 and min(F(r['value'][0])for r in rows if r['value'][1]=='0')==F(1,35)
    return dict(profiles=[B.wire(p)for p in profiles],rows=rows,minimum_gap=['1/35','0'])


def small_cycles():
    rows=[];single=[]
    ps=[p for p in B.cycle_profiles()if p['name']in('isolated','both_ends3')]
    for m in(14,15,16):
        P,Q=7,m-7;alpha=F(1,P)+F(1,Q)
        # Both neighboring corners for a single profile are finite by the
        # fixed simple leaf/degree-two-neighbor argument in the report.
        for side in(0,1):
            for name,pos,w2,kind in [('end3',1,1,'A'),('end4',1,2,'B'),('next3',2,1,'D')]:
                contact=pos if side==0 else m-pos
                ids=list(range(1,m))+[m]
                edges=[(i,i+1,1)for i in range(1,m-1)]+[(m,contact,w2)]
                adj={i:set()for i in ids}
                for i,j,_ in edges:adj[i].add(j);adj[j].add(i)
                fixed=m-1 if side==0 else 1;near=next(iter(adj[fixed]))
                assert len(adj[fixed])==1 and len(adj[near])==2
                dist=C.distances(adj,ids,m)[7];assert dist<=9
                single.append(dict(cycle_order=m,profile=f'{name}_{side}',
                    finite_type=f'{kind}{m}',positive_pivots=C.tree_pd(ids,edges),
                    marked_distance=dist,fixed_leaf=fixed,degree_two_neighbor=near))
        for c,d in product(ps,repeat=2):
            vals=[]
            for ub,vb in product((0,1),repeat=2):
                z=A.value(A.Q(F(1,P)),A.Q(F(1,Q)),c,d,
                          c['umax']if ub else A.Q(),d['umax']if vb else A.Q())
                assert z.b==0
                vals.append(dict(upper_reciprocals=[ub,vb],value=str(z.a)))
            lo=min(F(v['value'])for v in vals)
            if c['name']!='isolated' or d['name']!='isolated':assert lo>0
            else:assert lo==1-4*alpha
            rows.append(dict(cycle_order=m,c=c['name'],d=d['name'],values=vals,minimum=str(lo)))
        assert 1-3*alpha>0  # one isolated height >=sqrt(2) already gives contradiction
    required={m:4-1/(F(1,7)+F(1,m-7))for m in(14,15,16)}
    assert required=={14:F(1,2),15:F(4,15),16:F(1,16)}
    t=s.symbols('t',real=True)
    assert s.simplify(4/(4+2*t)-s.Rational(1,2)-(2-t)/(2*(t+2)))==0
    return dict(single_profile_checks=single,profile_rows=rows,
        isolated_height_gate=[dict(cycle_order=m,one_height_at_least_sqrt2_gap=str(1-3*(F(1,7)+F(1,m-7))))for m in(14,15,16)],
        required_unscaled_spectator_difference={str(m):str(v)for m,v in required.items()},
        spectator_gap=dict(two_disjoint_nonempty_supports='D>=4/(4+2t)>1/2 for0<=t<2',
            exactly_one_nonempty_support='D>=1/2; equality is a single unit contact to an isolated A1',
            application_order='Both native isolated heights are proved1 before applying the unscaled contact gap.'))


def tree_checks():
    ps=B.tree_profiles();rows=[]
    for P in(7,8,9):
        for c,d in product(ps,repeat=2):
            if P in(7,8)and any(p['name']=='inner_next3'for p in(c,d)):
                rows.append(dict(P=P,c=c['name'],d=d['name'],kind='immediate_payment',
                    path_length=P-1,ideal_reason='distance<=6'if P==7 else'new marked fork leaf is within1 of a feature'))
                continue
            costs=[2*p['moment']/p['minimum_height']+p['c0']/p['minimum_height']**2 for p in(c,d)]
            gap=A.Q(1)-sum(costs,A.Q())/P
            lower={7:F(1,7),8:F(1,4),9:F(1,9)}[P]
            assert(gap-A.Q(lower)).sign()>=0
            rows.append(dict(P=P,c=c['name'],d=d['name'],kind='positive_gap',gap=gap.wire()))
    assert Counter(r['kind']for r in rows)==dict(immediate_payment=30,positive_gap=162)
    return dict(profiles=[B.wire(p)for p in ps],rows=rows,feature_radius=6,
                minimum_outer_B_rank=8,minimum_outer_D_rank=9)


def cycle14_core():
    n=17;g=2*s.eye(n)
    for i in range(14):g[i,(i+1)%14]=g[(i+1)%14,i]=-1
    for i,j in((0,14),(7,15),(14,16)):g[i,j]=g[j,i]=-1
    assert C.inertia(g.tolist())==[15,1,1]
    z=s.Matrix([3,2,1,0,-1,-2,-3,-4,-3,-2,-1,0,1,2,2,-2,1])
    assert g*z==s.zeros(n,1)
    inv=g[1:,1:].inv()
    def lift(y):
        x=s.Matrix([0,*(-inv*y[1:,:])]);assert g*x==-y;return x
    y=s.ones(n,1);y[0]=3;h=lift(y);assert(h.T*g*h)[0]==-204
    for i in range(n):
        e=s.zeros(n,1);e[i]=y[i]/2;v=h+e;p=g*v
        assert p[i]==0 and all(p[j]<0 for j in range(n)if j!=i)
        assert(v.T*g*v)[0]==-204-y[i]**2/2
    W=[i for i in range(14)if i not in(0,7)]+[16]
    assert C.inertia(g.extract(W,W).tolist())==[13,0,0]
    adj={i:{j for j in range(n)if i!=j and g[i,j]}for i in range(n)}
    corners=[]
    for a,b in((0,7),(7,14),(14,15),(15,0)):
        ids=sorted(set(W)|{a,b});sig=C.inertia(g.extract(ids,ids).tolist())
        if(a,b)==(0,7):
            assert sig==[14,0,1]
            corners.append(dict(marks=[a,b],kind='ideal_cycle14_with_A1_spectator',coupon='0'));continue
        assert sig==[15,0,0]
        parts=[]
        for mark in(a,b):
            ds=C.distances(adj,ids,mark);leaves=[i for i in ds if len(adj[i]&ds.keys())<=1]
            assert max(len(adj[i]&ds.keys())for i in ds)<=2
            center=len(ds)==13 and len(leaves)==2 and all(ds[j]==6 for j in leaves)
            parts.append(dict(nodes=sorted(ds),type=f'A{len(ds)}',mark=mark,unique_middle=center))
        assert not(set(parts[0]['nodes'])&set(parts[1]['nodes']))
        sizes=sorted(len(p['nodes'])for p in parts)
        coupon=F(1,2)if sizes in([1,13],[2,13])else F()
        if coupon:assert any(p['unique_middle']for p in parts)
        corners.append(dict(marks=[a,b],kind='finite',signature=sig,marked_components=parts,coupon=str(coupon)))
    assert sum(F(r['coupon'])for r in corners)==1
    rays=[]
    for i in range(n):
        if not z[i]:y=s.zeros(n,1);y[i]=1;rays.append((['axis',i],y))
    for i,j in product([i for i in range(n)if z[i]>0],[i for i in range(n)if z[i]<0]):
        y=s.zeros(n,1);y[i]=-z[j];y[j]=z[i];rays.append((['pair',i,j],y))
    rows=[];counts=Counter();face=[]
    for label,y in rays:
        x=lift(y);q=(x.T*g*x)[0];kind='finite'if q<0 else'ideal'if q==0 else'exterior'
        counts[kind]+=1;row=dict(label=label,y=list(map(str,y)),norm=str(q),kind=kind)
        rows.append(row)
        if all(y[i]==0 for i in W):face.append(row)
        if q<=0:assert(h.T*g*x)[0]<0
    assert len(rows)==58 and counts==dict(finite=27,ideal=5,exterior=26)
    assert len(face)==4 and Counter(r['kind']for r in face)==dict(finite=3,ideal=1)
    r=next(r for r in rows if r['kind']=='exterior');v=lift(s.Matrix(list(map(s.Rational,r['y']))))
    q=(v.T*g*v)[0];hv=(h.T*g*v)[0];t=(hv+s.sqrt(hv**2+204*q))/204;u=v+t*h
    assert t>0 and s.simplify((u.T*g*u)[0])==0 and s.simplify((h.T*g*u)[0])<0
    assert all(s.simplify(x)<0 for x in g*u)
    return dict(gram=[[str(x)for x in row]for row in g.tolist()],signature=[15,1,1],radical=list(map(str,z)),
        strict_center=list(map(str,h)),center_norm='-204',white=W,corners=corners,
        rays=rows,ray_counts=dict(counts),face_rays=face,
        strict_escape=dict(y=r['y'],v_norm=str(q),h_dot_v=str(hv),t=str(t)),finite_volume=False,
        scope='An actual local equality quadrilateral in an explicitly infinite-volume H15 Coxeter chamber; no new finite-volume example.')


def graph(kind,m):
    if kind=='A':edges=[(i,(i+1)%m,3)for i in range(m)]
    elif kind=='C':edges=[(i,i+1,4 if i in(0,m-2)else 3)for i in range(m-1)]
    elif kind=='B':edges=[(i,i+1,4 if i==m-4 else 3)for i in range(m-3)]+[(0,m-2,3),(0,m-1,3)]
    else:
        assert kind=='D'
        edges=[(i,i+1,3)for i in range(m-5)]+[(0,m-4,3),(0,m-3,3),(m-5,m-2,3),(m-5,m-1,3)]
    adj={i:set()for i in range(m)}
    for i,j,_ in edges:adj[i].add(j);adj[j].add(i)
    feature={i for i in range(m)if len(adj[i])>=3}|{i for a,b,k in edges if k>=4 for i in(a,b)}
    ds=[C.distances(adj,range(m),i)for i in range(m)]
    near={i for i in range(m)if any(ds[i][j]<=6 for j in feature)}
    paid=[(i,j)for i,j in combinations(range(m),2)if ds[i][j]<=6 or(ds[i][j]in(7,8,9)and({i,j}&near))]
    return len(paid)


def affine_costs():
    rows=[]
    for m in range(18,101):
        for kind in'ABCD':
            actual=graph(kind,m)
            formula=6*m if kind=='A'else 6*m-{'B':14,'C':21,'D':7}[kind]+sum(min(m-k,{'B':15,'C':16,'D':14}[kind])for k in{'B':(9,10,11),'C':(7,8,9),'D':(11,12,13)}[kind])
            assert actual==formula
            rows.append(dict(kind=kind,order=m,cost=actual))
    maxima={kind:max((F(r['cost'],r['order']),r['order'])for r in rows if r['kind']==kind)for kind in'BCD'}
    assert maxima=={'B':(F(36,5),25),'C':(F(85,12),24),'D':(F(95,13),26)}
    branches=[]
    for kind,k,b0,b1,b2,b3 in [('C',23,-45,-22,2,27),('B',24,-44,-20,5,31),('D',25,-43,-18,8,35)]:
        bound=maxima[kind][0]
        for lo,hi,a,b in [(18,k,9,b0),(k,k+1,8,b1),(k+1,k+2,7,b2),(k+2,None,6,b3)]:
            slope=bound-a
            assert hi is not None or slope>0
            at=hi if slope<0 else lo;slack=slope*at-b
            assert slack>=0
            branches.append(dict(kind=kind,lower_order=lo,upper_order=hi,slope=a,intercept=b,
                normalized_bound=str(bound),minimum_slack_at=at,minimum_slack=str(slack)))
    return dict(rows=rows,formulas={'A':'6m','B':'6m-14+min(m-9,15)+min(m-10,15)+min(m-11,15)',
        'C':'6m-21+min(m-7,16)+min(m-8,16)+min(m-9,16)',
        'D':'6m-7+min(m-11,14)+min(m-12,14)+min(m-13,14)'},
        normalized_maxima={k:dict(value=str(v),order=m)for k,(v,m)in maxima.items()},all_rank_branches=branches,
        scope='Exact finite graph checks corroborate the all-rank spine proof in the report; this basic rule alone does not replace inherited selective preferred-component payments.')


def run():
    out=dict(large_cycle_checks=large_cycles(),small_cycle_checks=small_cycles(),tree_checks=tree_checks(),
             cycle14_equality=cycle14_core(),basic_affine_cost=affine_costs())
    stems=('fresh_ideal_distance7_mixed_replacement_20261017','fresh_cycle16_half_coupon_transfer_20261017',
           'fresh_distance8_ideal_corner_replacement_20261017')
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths+=[ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    out.update(status='PASS_MIXED_FINITE9_IDEAL6_QUADRILATERAL_LOCAL_REPLACEMENTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        scope='Universal quadrilateral replacements: ideal cycle distance7 at order>=15; cycle14 sealed A13-center half-coupons; classical tree distances7,8,9 outside feature radius6. Basic all-rank graph costs included. No triangle transfer, complete selective small-factor rule, global finite budget or dimension theorem.')
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'100large-cycle vertices,48small-cycle vertices,192tree cases,58rays,332graph counts',flush=True)
