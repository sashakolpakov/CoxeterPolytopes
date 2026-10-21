#!/usr/bin/env python3
"""Cycle12 antipodal transfer with sealed A1/A11 and A2/A11 coupons."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import sympy as s
import fresh_finite8_cycle_replacement_20261017 as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_cycle12_ideal5_transfer_20261017'
A,B,C=R.A,R.B,R.C

def inertia(matrix):
    g=[[F(x)for x in row]for row in matrix];positive=negative=zero=0;pivots=[]
    while g:
        n=len(g);i=next((i for i in range(n)if g[i][i]),None)
        if i is not None:
            ids=[i]+[j for j in range(n)if j!=i];g=[[g[a][b]for b in ids]for a in ids]
            p=g[0][0];positive+=p>0;negative+=p<0;pivots.append(str(p))
            g=[[g[a][b]-g[a][0]*g[0][b]/p for b in range(1,n)]for a in range(1,n)]
        else:
            pair=next(((i,j)for i in range(n)for j in range(i+1,n)if g[i][j]),None)
            if pair is None:zero+=n;break
            i,j=pair;ids=[i,j]+[k for k in range(n)if k not in pair];g=[[g[a][b]for b in ids]for a in ids]
            p=g[0][1];positive+=1;negative+=1;pivots.append(['zero_diagonal_pair',str(p)])
            g=[[g[a][b]-(g[a][0]*g[1][b]+g[a][1]*g[0][b])/p for b in range(2,n)]for a in range(2,n)]
    return dict(positive=positive,negative=negative,zero=zero,pivots=pivots)

def components(matrix,ids):
    left=set(ids);out=[]
    while left:
        seen={min(left)};todo=list(seen)
        while todo:
            i=todo.pop()
            for j in left-seen:
                if matrix[i][j]:seen.add(j);todo.append(j)
        left-=seen;out.append(sorted(seen))
    return out

def literal(name,zorder=0,zedges=(),xc=(),yd=(),hc=1,hd=1,both=None):
    # Square-root2 edges are made rational by positive root congruences.
    n=14+zorder;g=[[F(2)if i==j else F(0)for j in range(n)]for i in range(n)]
    def edge(i,j,x=1):g[i][j]=g[j][i]=-F(x)
    for i in range(12):edge(i,(i+1)%12)
    scales={}
    for v,old,h in [(12,0,hc),(13,6,hd)]:
        if h==2:scales[v]='sqrt2';g[v][v]=4;edge(v,old,2)
        else:edge(v,old)
    if both is not None:
        v,old=(12,0)if both==12 else(13,6)
        g[v][old]=g[old][v]=0
        edge(v,(old-1)%12);edge(v,(old+1)%12)
    for i,j,x in zedges:edge(14+i,14+j,x)
    for root,contacts in [(12,xc),(13,yd)]:
        for j,x in contacts:edge(root,14+j,x)
    # zedges/contacts supplied above may use an explicitly rationalized root.
    if name=='one_sqrt2_A1':g[14][14]=4;scales[14]='sqrt2'
    if name=='one_B4':g[14][14]=4;scales[14]='sqrt2'
    whites=[i for i in range(n)if i not in[0,6,12,13]]
    corners=[]
    for pair in [(0,6),(12,6),(0,13),(12,13)]:
        ids=whites+list(pair);sub=[[g[i][j]for j in ids]for i in ids];sig=inertia(sub)
        comps=components(g,ids);coupon=[]
        for mark,center,old in [(12,6,0),(13,0,6)]:
            if mark not in pair or center not in pair:continue
            chain=[(old+i)%12 for i in range(1,12)]
            large=next(c for c in comps if center in c)
            small=next(c for c in comps if mark in c)
            if set(large)==set(chain):
                assert center==chain[5]
                assert all(g[a][b]==(-1 if abs(i-j)==1 else 0)for i,a in enumerate(chain)for j,b in enumerate(chain)if i!=j)
                if len(small)==1:coupon.append(dict(kind='whole_A1_A11_middle',small=small,center=center,weight='1'))
                if len(small)==2 and g[small[0]][small[1]]==-1 and all(g[i][i]==2 for i in small):
                    coupon.append(dict(kind='whole_A2_A11_middle',small=small,center=center,weight='1/2'))
        corners.append(dict(black_pair=pair,inertia=sig,components=comps,coupons=coupon))
    full=inertia(g)
    edge_checks=[]
    for b in [0,6,12,13]:
        ids=whites+[b];sig=inertia([[g[i][j]for j in ids]for i in ids])
        assert(sig['positive'],sig['negative'],sig['zero'])==(len(ids),0,0)
        edge_checks.append(dict(black=b,inertia=sig))
    rejected=any(c['inertia']['negative']for c in corners)
    if name in ['twins_D6','twins_D7']:assert rejected
    else:
        assert not rejected and(full['positive'],full['negative'],full['zero'])==(n-2,1,1)
        assert corners[0]['inertia']['zero']==1
        if not name.startswith('twins_'):
            total=sum(F(c['weight'])for r in corners for c in r['coupons'])
            assert total>=1
        else:
            comp=next(c for c in corners[-1]['components']if 12 in c);assert 13 in comp
            adj={i:{j for j in comp if i!=j and g[i][j]}for i in comp}
            assert C.distances(adj,comp,12)[13]==4
    return dict(name=name,gram=[[str(x)for x in row]for row in g],rational_congruence=scales,
        white=whites,full_inertia=full,edge_checks=edge_checks,corners=corners,
        outcome='excluded by a negative direction at the proposed actual corner'if rejected else'retained short path'if name.startswith('twins_')else'sealed finite coupon')

def run():
    stems=['fresh_cycle13_ideal5_transfer_20261017','fresh_spectator_energy_one_boundary_20261017',
        'fresh_spectator_energy_below_one_20261017','fresh_finite8_cycle_replacement_20261017',
        'fresh_finite8_ideal5_local_replacement_20261017']
    for stem in stems:
        d=json.loads((ROOT/'results'/f'{stem}.json').read_text());assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    m=12;single=[]
    for p in B.cycle_profiles():
        if p['name']in['isolated','both_ends3']:continue
        side=next(i for i,v in enumerate(p['position'])if v);pos=p['position'][side]
        contact=pos if side==0 else m-pos;w2=2 if p['name'].startswith('end4')else 1
        ids=list(range(1,m+1));edges=[(i,i+1,1)for i in range(1,m-1)]+[(m,contact,w2)]
        adj={i:set()for i in ids}
        for i,j,_ in edges:adj[i].add(j);adj[j].add(i)
        leaf=m-1 if side==0 else 1;near=next(iter(adj[leaf]));assert len(adj[leaf])==1 and len(adj[near])==2
        dist=C.distances(adj,ids,m)[6];assert dist==7-pos and dist<=6
        single.append(dict(profile=p['name'],marked_distance=dist,fixed_leaf=leaf,
            fixed_degree_two_neighbor=near,positive_pivots=C.tree_pd(ids,edges)))
    keep=[p for p in B.cycle_profiles()if p['name']in['isolated','both_ends3']];rect=[];minima=[]
    for c,d in product(keep,repeat=2):
        vv=[]
        for u,v in product((0,1),repeat=2):
            x=A.value(A.Q(F(1,6)),A.Q(F(1,6)),c,d,c['umax']if u else A.Q(),d['umax']if v else A.Q())
            vv.append(x);rect.append(dict(c=c['name'],d=d['name'],upper_reciprocals=[u,v],value=x.wire()))
        expected=F(-1,3)if c['name']==d['name']=='isolated'else F(1,3)if c['name']==d['name']=='both_ends3'else F(0)
        assert R.minimum(vv)==A.Q(expected)
        minima.append(dict(c=c['name'],d=d['name'],minimum=str(expected)))
    h,k,p,q,r=s.symbols('h k p q r',nonzero=True);alpha=s.Rational(1,3)
    t=s.Matrix([[alpha,-alpha,-h,0],[-alpha,alpha,0,-k],[-h,0,2-p,-q],[0,-k,-q,2-r]])
    expected=1-2*alpha/h**2-2*alpha/k**2+alpha*(p/h**2-2*q/(h*k)+r/k**2)
    assert s.factor(t.det()/(h*h*k*k)-expected)==0
    records=[literal('isolated_1_sqrt2',hd=2),literal('isolated_sqrt2_1',hc=2),
        literal('mixed_isolated_both',both=13),literal('mixed_both_isolated',both=12),
        literal('two_A1',2,xc=[(0,1)],yd=[(1,1)]),
        literal('one_sqrt2_A1',1,xc=[(0,2)]),
        literal('one_two_A1',2,xc=[(0,1),(1,1)]),
        literal('one_A3_middle',3,[(0,1,1),(1,2,1)],xc=[(1,1)]),
        literal('one_B4',4,[(0,1,2),(1,2,1),(2,3,1)],xc=[(3,1)]),
        literal('one_D4',4,[(0,1,1),(1,2,1),(1,3,1)],xc=[(0,1)]),
        literal('twins_A3',3,[(0,1,1),(1,2,1)],xc=[(0,1)],yd=[(2,1)])]
    for rank in range(4,8):
        edges=[(i,i+1,1)for i in range(rank-3)]+[(rank-3,rank-2,1),(rank-3,rank-1,1)]
        records.append(literal(f'twins_D{rank}',rank,edges,xc=[(rank-2,1)],yd=[(rank-1,1)]))
    retained=[[i,j]for i in range(12)for j in range(i+1,12)if min(j-i,12-j+i)<=5]
    removed=[[i,(i+6)%12]for i in range(6)];assert len(retained)==60 and len(records)==15
    paths=[Path(__file__),Path(R.__file__),Path(A.__file__),Path(B.__file__),Path(C.__file__)]
    paths +=[ROOT/folder/f'{stem}.{ext}'for stem in stems for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    paths +=[ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    return dict(status='PASS_UNIVERSAL_CYCLE12_IDEAL5_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        single_profiles=single,rectangles=rect,rectangle_minima=minima,literal_cases=records,
        isolated_identity=dict(alpha='1/3',unit_height_energy='1',one_sqrt2_boundary='0',two_sqrt2_gap='1/3'),
        finite_coupons=[dict(small='whole_A1',large='whole_A11',large_position=5,weight='1'),dict(small='whole_A2',large='whole_A11',large_position=5,weight='1/2')],
        pair_rule=dict(retained=retained,transferred=removed,cost=60),
        counts=dict(single_profiles=6,rectangle_values=16,literal_cases=15,transferred_pairs=6),
        theorem='Every actual cycle12 antipodal quadrilateral transfers to a retained finite8/ideal5 payment or the specified sealed whole-component finite coupons. Arbitrary spectators and all native labels are covered by the written classification and energy-one theorem.',
        scope='Local transfer only. No finite vertex budget, triangle-retention/global composition or ambient dimension bound is claimed. Literal checks corroborate the all-rank proof; they do not enumerate all spectator ranks.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
