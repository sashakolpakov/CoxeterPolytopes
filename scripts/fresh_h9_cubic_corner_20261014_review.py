#!/usr/bin/env python3
"""Cubic7 replay independent of the producer's Fraction-polynomial engine."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from coxeter_polytopes.exact.cubic7 import Cubic7 as C

TAG='fresh_h9_cubic_corner_20261014_review'
ZERO=C(0);ONE=C(1)


def matrix(data):return [[C.from_json(x)for x in row]for row in data]
def vector(data):return [C.from_json(x)for x in data]
def mv(A,x):return [sum((a*b for a,b in zip(row,x)),ZERO)for row in A]
def mm(A,B):return [[sum((a*b for a,b in zip(row,col)),ZERO)for col in zip(*B)]for row in A]
def transpose(A):return list(map(list,zip(*A)))


def run():
    source='results/fresh_h9_cubic_corner_20261014_extension.json'
    base='results/fresh_cubic7_h9_seed1_h100_20261013.json'
    original=json.loads((ROOT/base).read_text());d=json.loads((ROOT/source).read_text())
    for doc in (d,json.loads((ROOT/'results/fresh_h9_cubic_corner_20261014.json').read_text())):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert original['seed_index']==1 and original['search_complete'] is True
    assert original['height_bound']==100
    assert original['completed_control_levels']==original['total_control_levels']==len(original['levels'])
    a= C.from_json(original['timelike_coefficient'])
    assert a==2*C(0,1) and C.from_json(d['timelike_coefficient'])==a
    assert [a.sign(i)for i in range(3)]==[1,-1,-1]
    Q=[[ZERO for _ in range(10)]for _ in range(10)]
    Q[0][0]=-a;Q[9][9]=ONE
    for i in range(1,9):Q[i][i]=C(2)
    for i,j in ((1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(3,8)):
        Q[i][j]=Q[j][i]=-ONE
    pair=lambda x,y:sum((a*b for a,b in zip(x,mv(Q,y))),ZERO)
    T=matrix(d['isometry']);Id=[[C(i==j)for j in range(10)]for i in range(10)]
    assert mm(T,T)==Id
    assert all(x.is_integral for row in T for x in row)
    assert mm(mm(transpose(T),Q),T)==Q
    assert sum((T[i][i]for i in range(10)),ZERO)==2
    # Coordinates are in the actual O_K basis of <-2theta>+E8+<1>.
    # T is its own inverse, so the same ten integral columns prove both directions.
    u=[ONE]+[ZERO]*9;Tu=mv(T,u)
    assert Tu[0]>0 and pair(Tu,Tu)==-a
    H=[x+y for x,y in zip(u,Tu)]
    assert H==vector(d['fixed_point']) and mv(T,H)==H and H[0]>0
    hnorm=pair(H,H);assert hnorm==C.from_json(d['fixed_point_norm'])<0
    height=a*H[0]*H[0]/(-hnorm)-ONE
    assert height==C.from_json(d['fixed_point_height'])==C(2,7,3)
    assert 0<height<16<100
    N=matrix(original['all_walls']);assert len(N)==17
    assert N[:9]==[[-ONE if j==i else ZERO for j in range(10)]for i in range(1,10)]
    products=[pair(H,r)for r in N];assert all(x<=0 for x in products)
    incident=[i for i,x in enumerate(products)if not x]
    assert incident==d['incident_original_indices']==[0,1,2,3,4,7,8]
    mapped=[N.index(mv(T,N[i]))for i in incident]
    assert mapped==d['incident_permutation']==[4,3,2,1,0,7,8]
    assert sorted(mapped)==incident
    roots=N[:]
    for r in N:
        image=mv(T,r)
        if image not in roots:roots.append(image)
    assert roots==matrix(d['normals']) and len(roots)==d['wall_count']==22
    assert all(mv(T,r)in roots for r in roots)
    norms=[pair(r,r)for r in roots];assert all(x in (ONE,C(2))for x in norms)
    finite={ZERO,ONE,C(2),C(3),C(2,1),C(1,1,1)}
    counts=dict(right=0,finite_nonright=0,dotted=0)
    for i,r in enumerate(roots):
        for j,s in enumerate(roots[:i]):
            g=pair(r,s);assert g<=0
            v=4*g*g/(norms[i]*norms[j])
            if not g:counts['right']+=1
            elif v>4:counts['dotted']+=1
            else:assert v in finite;counts['finite_nonright']+=1
    assert counts==d['pair_counts']==dict(right=128,finite_nonright=19,dotted=84)
    center=vector(d['strict_center']);assert pair(center,center)<0 and center[0]>0
    assert all(pair(center,r)<0 for r in roots)
    escape=vector(original['geometry']['vector']);ep=[pair(escape,r)for r in roots]
    assert ep==vector(d['old_escape_pairings'])
    assert any(x>0 for x in ep) and d['old_escape_preserved'] is False
    # The far corner is now certified by the chamber symmetry, without a
    # distance-bound assertion at the far point itself.
    tip=[pair(Tu,r)for r in roots]
    assert all(x<=0 for x in tip)
    tipids=[i for i,x in enumerate(tip)if not x]
    assert len(tipids)==9 and {tuple(roots[i])for i in tipids}=={tuple(mv(T,r))for r in N[:9]}
    paths=[source,base,'results/fresh_h9_cubic_corner_20261014.json',
           'scripts/fresh_h9_cubic_corner_20261014_extension.py',
           'scripts/fresh_h9_cubic_corner_20261014.py',
           'scripts/audit_fresh_cubic7_h9_no_import_20261013.py',
           'scripts/search_fresh_cubic7_h9_20261013.py',
           'src/coxeter_polytopes/search/cubic7_e8_cubes.py',
           'src/coxeter_polytopes/search/cubic7_root_blocks.py',
           'src/coxeter_polytopes/exact/cubic7.py',f'scripts/{TAG}.py']
    return dict(status='PASS_SEPARATE_CUBIC7_FIXED_POINT_SECTOR_AND_TRANSPORT_REPLAY',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                lattice_basis_columns=10,involution=True,trace='2',form_preserved=True,
                fixed_point_height=height.to_json(),height_strictly_below16=True,
                incident_indices=incident,positive_incident_permutation=mapped,
                old_walls=17,transported_union=22,pair_counts=counts,
                far_corner_incident_indices=tipids,old_escape_cutters=[i for i,x in enumerate(ep)if x>0],
                scope='Different Cubic7 field arithmetic from producer Fraction-polynomial core. '
                      'Independent exact matrix/module/sector/transport and all231 pair checks. '
                      'Full-chamber identification inherits the complete height100 Vinberg prefix and norm alphabet; '
                      'no positive-height root enumeration is repeated. Finite-order symmetry proves neither '
                      'reflectivity nor nonreflectivity, and cutting the old escape is not a volume certificate.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    dest=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['transported_union'],d['pair_counts'])
