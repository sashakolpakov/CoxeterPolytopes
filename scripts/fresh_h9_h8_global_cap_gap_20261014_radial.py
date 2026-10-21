#!/usr/bin/env python3
"""Uniform radial-cover reduction for every all-dotted cap direction."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
from fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate import I,ei,SCALE,BITS
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_global_cap_gap_20261014_radial'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
SOURCE='data/known_examples/bugaenko_compact_h8.json'
GAP='results/fresh_h9_h8_facet_cover_20261014_gap_gate.json'
PAIR='results/fresh_h9_h8_corner_anchors_20261014_pair_gate.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for q,h in d.get('sources',{}).items():assert sha(q)==h,q
    return d
def run():
    base,src,gap,pair=[read(p)for p in(BASE,SOURCE,GAP,PAIR)]
    assert pair['status']=='PASS_ALL_LABEL_TWO_EXCEPTIONAL_CENTER_ANCHORS_EXCLUDED'
    assert gap['source_dimension']==8 and gap['source_facets']==11
    a=A.E(A.Q(F(-1,4),F(1,4)));B=4*a;DA=A.E(A.Q(9,4));DB=A.E(A.Q(17,8))
    assert gap['exceptional_minimum_norm']==a.a.enc()and gap['ordinary_minimum_norm']==B.a.enc()
    def sf(x):return A.E(A.Q(F(x['a'])+F(x['b'])/2,F(x['b'])/2))
    Qs=[[sf(x)for x in r]for r in src['lattice_gram']]
    Ns=[[sf(x)for x in r]for r in src['simple_roots']]
    Hs=[sf(x)for x in src['geometry']['interior_point']]
    verts={}
    for vi in(0,32,33):
        ids=src['geometry']['vertices'][vi];v=A.kernel_one([A.mv(Qs,Ns[i])for i in ids]);assert v is not None
        if A.pair(v,Hs,Qs).sign()>0:v=[-x for x in v]
        h=[-A.pair(v,n,Qs)for n in Ns];assert all(x.sign()>=0 for x in h)
        positive=[x for x in h if x.sign()>0];m=positive[0]
        for x in positive:
            if(x-m).sign()<0:m=x
        v=[x/m for x in v];h=[x/m for x in h]
        assert A.pair(v,v,Qs)==-a and any(x==1 for x in h)
        verts[vi]=(v,h)
    corner_pairs=[]
    for i,j in combinations((0,32,33),2):
        ch=-A.pair(verts[i][0],verts[j][0],Qs)/a
        assert ch==(DA if i==0 else DB)
        corner_pairs.append(dict(pair=[i,j],cosh=A.enc(ch)))
    Q,N=A.dm(base['form']),A.dm(base['final_normals']);assert len(N)==12
    norms=[A.pair(n,n,Q)for n in N]
    strong=[];products=[];adj=[[]for _ in range(12)]
    for i,j in combinations(range(12),2):
        g=A.pair(N[i],N[j],Q);assert g.sign()<=0
        C=-2*ei(g)/(ei(norms[i])*ei(norms[j])).sqrt()
        # Comparison to1 is exact even if normalization adjoins a radical.
        sg=(4*g*g-norms[i]*norms[j]).sign()
        if sg>0:
            assert C.lo>SCALE;strong.append([i,j]);adj[i].append(j);adj[j].append(i)
        elif sg<0:assert C.hi<SCALE
        else:assert C.overlap(I(1))
        products.append(dict(pair=[i,j],normalized_negative_product=C.enc(),comparison_to_one=sg))
    assert strong==[[0,1],[4,11],[5,6],[8,11]]
    assert all(len(adj[i])<=1 for i in range(11))and adj[11]==[4,8]
    aa,BB=ei(a),ei(B)
    def reach(low_norm,other_norm,C):
        dc=((low_norm+2)*(other_norm+2)).sqrt()+C
        dc=dc/(low_norm*other_norm).sqrt()
        return dc*((other_norm+2)/other_norm).sqrt()+(dc.square()-1).sqrt()*(2/other_norm).sqrt()
    def union_same(other_norm,C):
        dc=(other_norm+2+C)/other_norm
        cr=((other_norm+2)/other_norm).sqrt();sr=(2/other_norm).sqrt()
        return dc*(cr.square()+sr.square())+(dc.square()-1).sqrt()*2*cr*sr
    ordinary_reach=reach(aa,BB,I(1));assert ordinary_reach.hi<ei(DA).lo
    own_radius=((aa+2)/aa).sqrt();assert own_radius.hi<ei(DA).lo
    high_diameter=1+4/BB;assert high_diameter.hi<ei(DA).lo
    # Only4 and8 can cover the two remote corners from low11.
    assert not A.pair(N[4],N[8],Q)
    union48=union_same(BB,I(0));assert union48.hi<ei(DB).lo
    C4=next(I.raw(*map(int,r['normalized_negative_product']))for r in products if r['pair']==[4,11])
    C8=next(I.raw(*map(int,r['normalized_negative_product']))for r in products if r['pair']==[8,11])
    assert C4.lo>C8.hi
    phi=ei(A.PHI)
    phi_reach=reach(aa*phi.square(),BB,C4);assert phi_reach.hi<ei(DB).lo
    sqrt2_ball8=reach(2*aa,BB,C8);assert sqrt2_ball8.hi<ei(DB).lo
    # The solelow scale is a source contact: 1,sqrt2, or>=phi.
    assert (A.PHI*A.PHI-2).sign()>0
    paths=[BASE,SOURCE,GAP,PAIR,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
        'scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py',f'scripts/{TAG}.py']
    return dict(status='PASS_UNIFORM_SINGLE_LOW_CENTER_RADIAL_REDUCTION',sources={p:sha(p)for p in paths},
        interval_bits=BITS,interval_scale=str(SCALE),small_norm=A.enc(a),high_norm_threshold=A.enc(B),
        source_exceptional_points={str(i):dict(point=A.enc(v),contacts=A.enc(h))for i,(v,h)in verts.items()},
        source_exceptional_pair_cosh=corner_pairs,old_pair_products=products,
        greater_than_one_graph=strong,neighbor_lists=adj,
        ordinary_neighbor_max_reach_cosh=ordinary_reach.enc(),own_ball_max_radius_cosh=own_radius.enc(),
        high_ball_max_diameter_cosh=high_diameter.enc(),high_balls_4_8_union_diameter_cosh=union48.enc(),
        scale_at_least_phi_max_reach_cosh=phi_reach.enc(),sqrt2_scale_ball8_max_reach_cosh=sqrt2_ball8.enc(),
        eliminated_low_old_walls=list(range(11)),eliminated_low11_source_corner=0,
        surviving_low_branches=[dict(old_wall=11,source_corner=k,scale=s,
            required_assignments=([[4,'opposite_E8'],[8,0]]if s=='sqrt2'else'Both assignments to the other two exceptional corners remain possible'))for k in(32,33)for s in('1','sqrt2')],
        conclusion='Every compatible perpendicular placement of the fixedH8source COVERED BY THE TWELVE OLD BALLS has either no projected norm belowB, or its unique low center is oldwall11 at sourcecorner32/33 with scale1 orsqrt2. At sqrt2, ball4 covers the oppositeE8corner and ball8 coverscorner0.',
        scope='Uniform over all real all-dotted cap directions retaining the frozen12 old walls, conditional on source coverage by those twelve balls. This is necessary for the raw attachment with the source perpendicular sidewalls; additional nonperpendicular covering walls are outside the coverage claim. Inherits the exact at-most-one-low-center theorem and norm-gap source lemma. No placement existence, no-low branch exclusion, remaining four low-branch exclusion, or generalH9 nonexistence is claimed.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['greater_than_one_graph'])
