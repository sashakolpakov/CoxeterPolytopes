#!/usr/bin/env python3
"""Radius-four affine-tree transfer, including the D7/E8 extension."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import fresh_distance8_ideal_corner_replacement_20261017 as B
import fresh_spectator_energy_below_one_20261017 as E

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_radius4_affine_tree_transfer_20261017'
A=B.A


def graph(kind,n):
    if kind=='A':edges=[(i,i+1,1)for i in range(n-1)]
    elif kind=='B':edges=[(i,i+1,2 if i==n-2 else 1)for i in range(n-1)]
    else:
        assert kind=='D'and n>=3
        edges=[(i,i+1,1)for i in range(n-3)]+[(n-3,n-2,1),(n-3,n-1,1)]
    return E.matrix(n,edges,kind=='B')


def extend_path(g,L,x,y):
    r=len(g);n=r+L+1
    out=[[F(2 if i==j else 0)for j in range(n)]for i in range(n)]
    for i in range(r):
        for j in range(r):out[i][j]=g[i][j]
    for i in range(r,r+L):out[i][i+1]=out[i+1][i]=F(-1)
    for i,v in x.items():out[i][r]=out[r][i]=-F(v)
    for i,v in y.items():out[i][-1]=out[-1][i]=-F(v)
    return out


def distance(g,i,j):
    ds={i:0};todo=[i]
    for k in todo:
        for l,v in enumerate(g[k]):
            if l!=k and v and l not in ds:ds[l]=ds[k]+1;todo.append(l)
    return ds.get(j)


def energy_one_corridors():
    """Corroborate every classified family; arbitrary ranks use graph identities."""
    rows=[]
    cases=[('sqrt2_A1',[[F(4)]],{0:2},{},'B',2),
           ('two_A1_one_side',[[F(2),F(0)],[F(0),F(2)]],{0:1,1:1},{},'D',3),
           ('A3_middle',graph('D',3),{0:1},{},'D',4),
           ('two_A1_opposite',[[F(2),F(0)],[F(0),F(2)]],{0:1},{1:1},'A',3)]
    for kind,ns in [('B',(2,6,12,64)),('D',(4,7,12,64))]:
        for r in ns:cases.append((f'{kind}{r}_long_leaf',graph(kind,r),{0:1},{},kind,r+1))
    for L in (1,5,9):
        for name,g,x,y,kind,offset in cases:
            corner=extend_path(g,L,x,y);n=len(corner)
            assert n==L+offset and distance(corner,len(g),n-1)==L
            pivots=E.tree_pivots(corner)
            rows.append(dict(corridor_length=L,spectator=name,finite_type=f'{kind}{n}',
                marked_distance=L,positive_pivots=pivots))
        # Whole A3 endpoints close a simple affine cycle, with alternate path4.
        g=graph('A',3);corner=extend_path(g,L,{0:1},{2:1});n=len(corner)
        assert n==L+4 and all(sum(row)==0 for row in corner)
        assert all(sum(x!=0 for j,x in enumerate(row)if i!=j)==2 for i,row in enumerate(corner))
        principal=[row[:-1]for row in corner[:-1]];E.tree_pivots(principal)
        assert distance(corner,3,n-1)==min(L,4)
        rows.append(dict(corridor_length=L,spectator='A3_endpoints',affine_cycle_order=n,
            positive_kernel=[1]*n,shorter_marked_distance=min(L,4)))
        # The same cycle plus a D tail has a zero quadratic vector not in the
        # matrix kernel, which is impossible for any PSD matrix.
        for r in range(4,8):
            g=graph('D',r);corner=extend_path(g,L,{r-2:1},{r-1:1})
            ids={r-3,r-2,r-1}|set(range(r,r+L+1))
            z=[F(int(i in ids))for i in range(len(corner))]
            gz=[sum((a*b for a,b in zip(row,z)),F())for row in corner]
            assert sum((a*b for a,b in zip(z,gz)),F())==0 and any(gz)
            assert gz[r-4]==-1
            rows.append(dict(corridor_length=L,spectator=f'D{r}_twins',excluded_PSD=True,
                zero_quadratic_vector=list(map(str,z)),gram_product=list(map(str,gz))))
    assert len(rows)==51
    return dict(rows=rows,arbitrary_length_and_rank=True,
        classification='Finite B/D/A extensions of the path, or a whole affine cycle of order L+4 with alternate marked path4. D4..D7 twin tails contradict PSD.',
        scope='Exact representative identities corroborate the complete all-rank energy-one classification and explicit diagram concatenations in the written proof.')


def e8_profile():
    g=graph('D',7)
    # D7 short-leaf inverse column and long-port pairing, checked literally.
    v=[F(i+1,2)for i in range(5)]+[F(7,4),F(5,4)]
    E.check_column(g,v,5)
    assert v[5]==F(7,4)and v[0]==F(1,2)
    ext=extend_path(g,0,{5:1},{})
    pivots=E.tree_pivots(ext)
    deg=[sum(x!=0 for j,x in enumerate(row)if j!=i)for i,row in enumerate(ext)]
    assert sorted(deg)==[1,1,1,2,2,2,2,3]
    center=deg.index(3);leaves=[i for i,d in enumerate(deg)if d==1]
    assert sorted(distance(ext,center,i)for i in leaves)==[1,2,4]
    c0=2-v[5];h=v[0];assert c0/h**2==1
    return dict(outer_type='D7',new_finite_type='E8',contact=5,
        inverse_contact_column=list(map(str,v)),outer_projection_norm='7/4',
        long_port_pairing='1/2',minimum_height='1/2',moment='0',c0='1/4',maximum_cost='1',
        positive_pivots=pivots,arm_lengths=[1,2,4])


def tree_cases():
    ps=B.tree_profiles()+[dict(name='outer_D7_short_E8',moment=A.Q(),c0=A.Q(F(1,4)),minimum_height=A.Q(F(1,2)))]
    maxima={p['name']:2*p['moment']/p['minimum_height']+p['c0']/p['minimum_height']**2 for p in ps}
    assert all((A.Q(4)-v).sign()>=0 for v in maxima.values())
    rows=[]
    for P in range(5,10):
        for c,d in product(ps,repeat=2):
            names=(c['name'],d['name'])
            if P<=8 and 'inner_next3'in names:
                rows.append(dict(P=P,c=names[0],d=names[1],kind='immediate_payment',
                    marked_path=P-1,ideal_reason='distance4'if P==5 else'persistent feature at distance1'))
                continue
            gap=A.Q(1)-(maxima[names[0]]+maxima[names[1]])/P
            if gap.sign()>0:
                rows.append(dict(P=P,c=names[0],d=names[1],kind='positive_gap',minimum=gap.wire()))
            elif gap.sign()==0:
                if P==6:assert names==('inner_end3','inner_end3')
                else:assert P==5 and 'inner_end3'in names and any(n in names for n in('isolated','inner_end4'))
                rows.append(dict(P=P,c=names[0],d=names[1],kind='sealed_finite_equality',
                    finite_marked_distance=P,sealing='No shared spectator contact and zero positive-definite correction force all contacts to vanish.'))
            else:
                assert P==5 and names==('inner_end3','inner_end3')and gap==A.Q(F(-1,5))
                high=F(1)-(F(5,4)+3)/5;assert high==F(3,20)>0
                rows.append(dict(P=P,c=names[0],d=names[1],kind='energy_one_corridor',
                    high_height_gap=str(high),forced_heights=['1','1'],required_energy='1',
                    middle_path_length=5,finite_payment_distance=5,ideal_payment_distance=4))
    counts=Counter(r['kind']for r in rows)
    assert len(rows)==405 and counts==dict(immediate_payment=68,positive_gap=331,sealed_finite_equality=5,energy_one_corridor=1)
    return dict(profiles=[dict(**B.wire(p),maximum_cost=maxima[p['name']].wire())for p in ps],
        rows=rows,counts=dict(counts),minimum_outer_B_rank=6,minimum_outer_D_rank=7,
        finite_payment_radius=8,ideal_direct_radius=4,ideal_feature_radius=4)


def run():
    stems=('fresh_spectator_energy_one_boundary_20261017',E.TAG,B.TAG,
        'fresh_finite8_ideal5_local_replacement_20261017')
    paths=[Path(__file__)]
    for stem in stems:
        doc=json.loads((ROOT/f'results/{stem}.json').read_text())
        assert doc['status'].startswith('PASS')
        for name,h in doc['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,name
        paths +=[ROOT/f'{folder}/{stem}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths +=[ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    return dict(status='PASS_RADIUS4_AFFINE_TREE_LOCAL_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        energy_one_corridor=energy_one_corridors(),outer_exception=e8_profile(),tree=tree_cases(),
        theorem='Actual classical affine tree corner with marked distance5..9, both marks farther than4 from every feature, transfers to another finite corner of marked distance<=8 or ideal corner of marked distance<=4 or distance5..9 with a mark within4 of a feature.',
        scope='Complete local tree transfer with arbitrary supporting spectators and no simple-cusp hypothesis. No corresponding cycle4 theorem, complete small-factor payment rule, triangle theorem or new global dimension bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    data=json.loads(json.dumps(run()));target=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(target.read_text())
    else:target.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['tree']['counts'],flush=True)
