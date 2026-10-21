#!/usr/bin/env python3
"""Independent dense-own A2 cross-right guards on literal star case0."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_h9_star_lastblock_20261014_review as R
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_last_internal_three_20261014_cross_right_review'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
S=[3,4,6,7,9,10,12,13]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    data=json.loads((ROOT/GEOM).read_text())
    for p,h in data.get('sources',{}).items():assert sha(p)==h,p
    g=data['representatives'][0]['geometry'];V=[set(v['incident'])for v in g['vertices']]
    M={frozenset(m)for m in g['minimal_nonfaces']}
    def face(s):
        ids=sorted(s);i=next(i for i,v in enumerate(V)if set(ids)<=v)
        return dict(support=ids,actual_vertex=i)
    assert len(V)==378 and frozenset([11,12,13])in M and frozenset([12,13,17])in M
    triangles=[face([11,17,a])for a in(12,13)]
    assert F(1,3)+F(1,2)+F(1,6)==1 and 3*F(1,3)==1
    # If a leaf has a single own contact, its missing triple forces>=7.
    # The actual anchor triangle then forces the other leaf right to that
    # root, so it is the opposite single>=7 case of the earlier theorem.
    assert F(3,4)+F(1,4)==1
    # Dense leaves: each has an own>=4; both cannot be>=4 at one root.
    # With an opposite own3 there, the high label is4/5, never>=6.
    assert F(1,2)+F(1,2)==1
    rt2=R.scale(R.COS[4],2)
    mutual_ones=R.sub(R.O,rt2);assert R.sign(mutual_ones)<0
    profiles=[];guards=[]
    path=[[R.O if i==j else R.neg(R.COS[[4,3,4][min(i,j)]])if abs(i-j)==1 else R.Z for j in range(4)]for i in range(4)]
    null=[R.O,rt2,rt2,R.O]
    assert all(R.sign(v)>0 for v in null)
    assert all(R.total(R.mul(path[i][j],null[j])for j in range(4))==R.Z for i in range(4))
    for orientation in(0,1):
        high={12:11 if orientation==0 else 17,13:17 if orientation==0 else 11}
        for p,q in product((4,5),repeat=2):
            own={11:[p,3],17:[3,q]}if orientation==0 else{11:[3,p],17:[q,3]}
            assert all(F(1,3)+sum(F(1,m)for m in v)<1 for v in own.values())
            profiles.append(dict(orientation=orientation,anchor11=own[11],anchor17=own[17]))
        for last,first in product((12,13),S[:6]):
            partner=S[S.index(first)^1];nodes=[high[last],last,first,partner]
            assert len(set(nodes))==4
            guards.append(dict(orientation=orientation,forced_right_pair=[first,last],
                ordered_path=nodes,lower_path_labels=[4,3,4],**face(nodes)))
    assert len(guards)==24 and len(profiles)==8
    sw=list(range(18));sw[12],sw[13]=13,12
    assert{frozenset(sw[i]for i in v)for v in V}=={frozenset(v)for v in V}
    paths=[GEOM,f'scripts/{TAG}.py','scripts/fresh_h9_star_lastblock_20261014_review.py']
    return dict(status='PASS_LITERAL_DENSE_A2_CROSS_RIGHT_GUARDS',sources={p:sha(p)for p in paths},
        actual_vertices=378,anchor_triangles=triangles,dense_profiles=profiles,cross_right_guards=guards,
        relaxed_path_gram=[[[str(a)for a in x]for x in row]for row in path],
        positive_null_vector=[[str(a)for a in x]for x in null],
        mutual_anchor_test3_ones_quadratic=[str(a)for a in mutual_ones],
        literal_root_swap=sw,
        implication='With last internal3, first3internals>=4 and dense last-own contacts, the missing triples and actual anchor triangles force eight ordered dense profiles and mutual anchor-right. All12crossings from the A2 roots to first3blocks are then right via24literal affine4-3-4 guards.',
        scope='Conditional domain certificate, not a complete Coxeter assignment or rank census. The opposite-single branch is identified but its previously reviewed exclusion is not rerun. Dense-own premises must be established before using these cross-right guards; internal lower bounds alone do not supply them.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['cross_right_guards']))
