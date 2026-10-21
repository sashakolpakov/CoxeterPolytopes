#!/usr/bin/env python3
"""Mixed finite9/ideal7 local lemmas and exact order16 obstruction.

No ambient dimension bound or finite-volume construction is asserted.
"""
import argparse,json
from collections import Counter,deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
import fresh_distance8_ideal_corner_replacement_20261017 as B
A=B.A
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ideal_distance7_mixed_replacement_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'

def signature(g):
    g=g.copy();out=[0,0,0]
    while g.rows:
        i=next((i for i in range(g.rows)if g[i,i]),None)
        if i is not None:
            d=g[i,i];assert d.is_positive or d.is_negative
            out[0 if d>0 else 1]+=1;rest=[j for j in range(g.rows)if j!=i]
            g=g.extract(rest,rest)-g.extract(rest,[i])*g.extract([i],rest)/d
        else:
            edge=next(((i,j)for i in range(g.rows)for j in range(i+1,g.rows)if g[i,j]),None)
            if edge is None:out[2]+=g.rows;break
            i,j=edge;out[0]+=1;out[1]+=1;rest=[k for k in range(g.rows)if k not in edge]
            g=g.extract(rest,rest)-g.extract(rest,list(edge))*g.extract(list(edge),list(edge)).inv()*g.extract(list(edge),rest)
    return out

def pd(g):return signature(g)==[g.rows,0,0]
def wire_matrix(g):return [[str(g[i,j])for j in range(g.cols)]for i in range(g.rows)]

def cycle_checks():
    ps=B.cycle_profiles();rows=[];minimum=None
    for c,d in product(ps,repeat=2):
        lower=[8,8]
        for p in(c,d):
            if p['name']not in('isolated','both_ends3'):
                # The corresponding adjacent corner is finite; its marked
                # distance is P-j+1, so unpaidness forces P>=9+j.
                for j,pos in enumerate(p['position']):
                    if pos:lower[j]=max(lower[j],9+pos)
        boxes=[lower]if sum(lower)>=17 else [[8,9],[9,8]]
        for lo in boxes:
            for ub,vb in product((0,1),repeat=2):
                z=A.value(A.Q(F(1,lo[0])),A.Q(F(1,lo[1])),c,d,
                    c['umax']if ub else A.Q(),d['umax']if vb else A.Q())
                assert(z-A.Q(F(1,18))).sign()>=0
                if minimum is None or(z-minimum).sign()<0:minimum=z
                rows.append(dict(c=c['name'],d=d['name'],arc_lower_bounds=lo,
                    upper_reciprocals=[ub,vb],normalized_determinant=z.wire()))
    assert len(rows)==272 and minimum==A.Q(F(1,18))
    return dict(profiles=[B.wire(p)for p in ps],rectangle_checks=rows,minimum_gap=minimum.wire())

def tree_checks():
    ps=B.tree_profiles();rows=[]
    for P in(8,9):
        for c,d in product(ps,repeat=2):
            # At P8 an inner-next3 root has a white-interior path of length7
            # to the opposite old mark, immediately paid at any corner.
            if P==8 and any(p['name']=='inner_next3'for p in(c,d)):
                rows.append(dict(P=P,c=c['name'],d=d['name'],kind='immediate_path7'))
                continue
            costs=[2*p['moment']/p['minimum_height']+p['c0']/p['minimum_height']**2 for p in(c,d)]
            gap=A.Q(1)-sum(costs,A.Q())/P
            expected=F(1,4)if P==8 else F(1,9)
            assert(gap-A.Q(expected)).sign()>=0
            rows.append(dict(P=P,c=c['name'],d=d['name'],kind='positive_gap',gap=gap.wire()))
    counts=Counter(r['kind']for r in rows);assert counts=={'positive_gap':113,'immediate_path7':15}
    return dict(profiles=[B.wire(p)for p in ps],rows=rows,counts=dict(counts),
        feature_radius=7,minimum_outer_B_rank=9,minimum_outer_D_rank=10)

def triangles():
    diagrams=json.loads((ROOT/f'results/{CENSUS}.json').read_text())['diagrams'];rows=[]
    for d in diagrams:
        if d['order']!=10:continue
        adj=[set()for _ in range(10)]
        for i,j,m in d['edges']:adj[i].add(j);adj[j].add(i)
        for marks in combinations(range(10),3):
            paths=[]
            for a,b in combinations(marks,2):
                forbidden=set(marks)-{a,b};q=deque([a]);prev={a:None}
                while q:
                    i=q.popleft()
                    for j in sorted(adj[i]-forbidden):
                        if j not in prev:prev[j]=i;q.append(j)
                path=None
                if b in prev:
                    path=[];v=b
                    while v is not None:path.append(v);v=prev[v]
                    path=path[::-1]
                paths.append(dict(pair=[a,b],path=path,length=len(path)-1 if path else None))
            assert sum(p['length']is not None and p['length']<=7 for p in paths)>=2
            rows.append(dict(diagram=d['index'],marks=list(marks),paths=paths))
    assert len(rows)==360
    return dict(order10_diagrams=[d['index']for d in diagrams if d['order']==10],
        rows=rows,max_second_path=7,lower_order_proof='A Steiner tree with at most8 edges connecting three marked vertices has two marked-to-marked paths of length at most7 whose interiors contain no other marked vertex.')

def order16_star():
    g=2*s.eye(18)
    for i in range(16):j=(i+1)%16;g[i,j]=g[j,i]=-1
    g[0,16]=g[16,0]=g[8,17]=g[17,8]=-1
    assert signature(g)==[16,1,1]
    z=g.nullspace()[0];assert sum(z)==0
    inv=g[1:,1:].inv()
    def lift(y):
        x=s.Matrix([0,*(-inv*y[1:,:])]);assert g*x==-y;return x
    h=lift(s.ones(18,1));assert(h.T*g*h)[0]==-204
    # Each normal has an actual strict facet witness.
    for i in range(18):
        e=s.zeros(18,1);e[i]=s.Rational(1,2);v=h+e;p=g*v
        assert p[i]==0 and all(p[j]<0 for j in range(18)if j!=i)
        assert(v.T*g*v)[0]<0
    positive=[i for i in range(18)if z[i]>0];negative=[i for i in range(18)if z[i]<0]
    rays=[]
    for i in range(18):
        if not z[i]:y=s.zeros(18,1);y[i]=1;rays.append((['zero',i],y))
    for i,j in product(positive,negative):
        y=s.zeros(18,1);y[i]=-z[j];y[j]=z[i];rays.append((['pair',i,j],y))
    data=[];counts=Counter()
    for label,y in rays:
        x=lift(y);norm=(x.T*g*x)[0]
        kind='finite'if norm<0 else'ideal'if norm==0 else'exterior';counts[kind]+=1
        if norm<=0:assert(h.T*g*x)[0]<0
        data.append(dict(label=label,y=list(map(str,y)),norm=str(norm),kind=kind))
    assert len(data)==66 and counts==dict(finite=35,ideal=7,exterior=24)
    W=[i for i in range(16)if i not in(0,8)];black=[0,8,16,17]
    assert pd(g.extract(W,W))
    S=g.extract(black,black)-g.extract(black,W)*g.extract(W,W).inv()*g.extract(W,black)
    expected=s.Matrix([[s.Rational(1,4),-s.Rational(1,4),-1,0],[-s.Rational(1,4),s.Rational(1,4),0,-1],[-1,0,2,0],[0,-1,0,2]])
    assert S==expected
    corners=[]
    for pair in((0,8),(8,16),(16,17),(17,0)):
        sig=signature(g.extract(W+list(pair),W+list(pair)))
        assert sig in([16,0,0],[15,0,1])
        corners.append(dict(pair=list(pair),signature=sig,old_payment=1 if pair==(0,8)else 0,
            ideal7_whole15_payment=0))
    face=[r for r in data if all(s.Rational(r['y'][i])==0 for i in W)]
    assert len(face)==4 and Counter(r['kind']for r in face)==dict(finite=3,ideal=1)
    y=s.zeros(18,1);y[4]=1;v=lift(y)
    assert(v.T*g*v)[0]==1 and(h.T*g*v)[0]==-9
    t=(s.sqrt(285)-9)/204;u=v+t*h
    assert s.simplify((u.T*g*u)[0])==0
    assert s.simplify((h.T*g*u)[0])==-s.sqrt(285)
    assert all(s.simplify(x).is_negative for x in g*u)
    return dict(gram=wire_matrix(g),signature=[16,1,1],ambient_dimension=16,
        radical=list(map(str,z)),strict_center=list(map(str,h)),center_norm='-204',
        facet_witnesses='h+e_i/2, for all18 roots',extreme_rays=data,ray_counts=dict(counts),
        supporting_whites=W,black_order=black,projected_gram=wire_matrix(S),corners=corners,face_rays=face,
        strict_escape=dict(exterior_ray_y=[str(x)for x in y],v=list(map(str,v)),v_norm='1',h_dot_v='-9',
            t='(sqrt(285)-9)/204',future_null='v+t*h',all_native_contacts_strict=True),
        finite_volume=False,scope='An actual acute Coxeter chamber and actual3-finite/1-ideal quadrilateral in H16, with complete cone census and strict escape. It obstructs a local payment rule, not a theorem restricted to finite-volume chambers.')

def run():
    cyc=cycle_checks();tree=tree_checks();tri=triangles();star=order16_star()
    stems=['fresh_distance8_ideal_corner_replacement_20261017','fresh_selective_distance9_dimension108_20261017',CENSUS]
    paths=[Path(__file__)]
    for stem in stems:paths += [ROOT/f'scripts/{stem}.py',ROOT/f'results/{stem}.json',ROOT/f'results/{stem}.md']
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    return dict(status='PASS_MIXED_FINITE9_IDEAL7_LOCAL_LEMMAS_AND_ORDER16_OBSTRUCTION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        cycle_mixed_replacement=cyc,tree_mixed_replacement=tree,triangle_two_path7=tri,
        order16_local_obstruction=star,mixed_replacement_complete=True,
        mixed_replacement_hypotheses=dict(actual_quadrilateral=True,factor_order_at_least=17,
            marked_distance_in=[8,9],feature_radius=7,finite_short_distance=9,ideal_short_distance=7),
        scope='Universal local replacement outside a feature neighborhood, two path7 payments per triangle, and an infinite-volume local countermodel at factor order16. No ambient dimension bound and no finite-volume counterexample.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],'272cycle,128tree,360triples,66rays',flush=True)
