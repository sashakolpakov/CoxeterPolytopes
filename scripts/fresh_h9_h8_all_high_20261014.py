#!/usr/bin/env python3
"""Three exact norm obstructions to an all-high fixed-H8 raw attachment.

No optimization, cap-direction discretization, or finite angle cutoff is used.
The covering balls are precisely the twelve frozen old-wall balls.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
from fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate import I,ei,SCALE,BITS

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_h8_all_high_20261014'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
SOURCE='data/known_examples/bugaenko_compact_h8.json'
GAP='results/fresh_h9_h8_facet_cover_20261014_gap_gate.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for q,h in d.get('sources',{}).items():assert sha(q)==h,q
    return d
def union_cosh(x,y,C):
    """cosh(distance between centers + the two ball radii)."""
    d=(((x+2)*(y+2)).sqrt()+C)/(x*y).sqrt()
    sx=(2/x).sqrt();cx=((x+2)/x).sqrt()
    sy=(2/y).sqrt();cy=((y+2)/y).sqrt()
    assert d.lo>SCALE
    return d*(cx*cy+sx*sy)+(d.square()-1).sqrt()*(sx*cy+cx*sy)

def run():
    base,src,gap=[read(p)for p in(BASE,SOURCE,GAP)]
    B=A.E(A.Q(-1,1));DB=A.E(A.Q(17,8));L2=2*(B+2)
    assert gap['ordinary_minimum_norm']==B.a.enc()
    Q,N=A.dm(base['form']),A.dm(base['final_normals'])
    assert len(N)==12 and A.inertia(Q)==[9,1]
    norms=[A.pair(n,n,Q)for n in N]
    assert all(q.sign()>0 for q in norms)
    # Independently recover the two actual source vertices used by the proof.
    def sf(x):return A.E(A.Q(F(x['a'])+F(x['b'])/2,F(x['b'])/2))
    SQ=[[sf(x)for x in r]for r in src['lattice_gram']]
    SN=[[sf(x)for x in r]for r in src['simple_roots']]
    SH=[sf(x)for x in src['geometry']['interior_point']]
    vertices=[]
    for vi in(32,33):
        ids=src['geometry']['vertices'][vi]
        v=A.kernel_one([A.mv(SQ,SN[i])for i in ids]);assert v is not None
        if A.pair(v,SH,SQ).sign()>0:v=[-x for x in v]
        nv=A.pair(v,v,SQ);assert nv.sign()<0
        assert all(A.pair(v,n,SQ).sign()<=0 for n in SN)
        vertices.append((v,nv))
    prod=A.pair(vertices[0][0],vertices[1][0],SQ)
    assert prod.sign()<0 and prod*prod==DB*DB*vertices[0][1]*vertices[1][1]
    BB,DD,LL=ei(B),ei(DB),ei(L2)
    # A single high ball cannot contain both selected source vertices.
    single_diameter=1+4/BB;assert single_diameter.hi<DD.lo
    old_pairs=[];dotted=[]
    for i,j in combinations(range(12),2):
        g=A.pair(N[i],N[j],Q);assert g.sign()<=0
        rel=A.relation(g,norms[i],norms[j])
        C=-2*ei(g)/(ei(norms[i])*ei(norms[j])).sqrt()
        if rel[0]=='dotted':dotted.append([i,j]);assert C.lo>2*SCALE
        else:
            assert rel[0]in('right','finite')and rel[1]<=5
            assert (4*g*g-A.PHI*A.PHI*norms[i]*norms[j]).sign()<=0
        old_pairs.append(dict(pair=[i,j],physical_pairing=A.enc(g),normalized_C=C.enc(),relation=list(rel)))
    assert dotted==[[4,11],[5,6],[8,11]]
    finite_union=union_cosh(BB,BB,ei(A.PHI));assert finite_union.hi<DD.lo
    records=[]
    for i,j,support,cutoff in[(4,11,[0,2,9],4),(5,6,[0,1,2,3],3),(8,11,[0,2,5],3)]:
        assert all(norms[k]==2 for k in support)
        assert all(not A.pair(N[k],N[z],Q)for k in support for z in(i,j))
        G=[[A.pair(N[k],N[l],Q)for l in support]for k in support]
        assert A.inertia(G)==[len(support),0]
        J=A.inverse(G);assert all(x.sign()>=0 for row in J for x in row)
        inv_sum=sum((sum(row,A.Z)for row in J),A.Z)
        C=I.raw(*map(int,next(r['normalized_C']for r in old_pairs if r['pair']==[i,j])))
        kk=A.E(F(cutoff*cutoff,2)-2);assert (kk-B).sign()>0
        diameter_at_cutoff=union_cosh(ei(kk),BB,C)
        assert diameter_at_cutoff.hi<DD.lo
        # If the two balls cover32/33, both cap-anchor products are<cutoff.
        # The anchor Lorentz block contributes at worst-2*cutoff^2/(C-2).
        positive_cost=LL*ei(inv_sum)
        penalty=2*cutoff*cutoff/(C-2)
        lower=positive_cost-penalty;assert lower.lo>2*SCALE
        records.append(dict(pair=[i,j],normalized_C=C.enc(),positive_support=support,
            support_gram=A.enc(G),support_inverse=A.enc(J),inverse_entry_sum=A.enc(inv_sum),
            anchor_contact_upper=cutoff,projected_norm_at_cutoff=A.enc(kk),
            ball_union_diameter_at_cutoff_cosh=diameter_at_cutoff.enc(),
            positive_projection_cost_lower=positive_cost.enc(),
            lorentz_projection_penalty_upper=penalty.enc(),cap_norm_lower=lower.enc()))
    paths=[BASE,SOURCE,GAP,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
        'scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py',f'scripts/{TAG}.py']
    return dict(status='PASS_ALL_HIGH_FIXED_H8_RAW_ATTACHMENT_EXCLUDED',
        sources={p:sha(p)for p in paths},interval_bits=BITS,interval_scale=str(SCALE),
        cap_normalization=2,old_normalization=2,high_threshold=A.enc(B),minimum_contact_square=A.enc(L2),
        source_vertex_ids=[32,33],source_vertices=[dict(point=A.enc(v),norm=A.enc(q))for v,q in vertices],
        source_vertex_distance_cosh=A.enc(DB),single_high_ball_diameter_cosh=single_diameter.enc(),
        old_pair_records=old_pairs,finite_pair_union_diameter_cosh=finite_union.enc(),
        required_dotted_pairs=dotted,three_pair_obstructions=records,
        theorem='No real all-dotted cap with every normalized projected squared norm at least sqrt5-1 admits a placement of the fixed compact H8 source covered by the twelve frozen old-wall balls.',
        proof='The two actual source vertices32/33 require distinct balls. All nondotted old pairs fail their distance. For each remaining pair, monotonicity of ball-union diameter bounds both normalized cap contacts. Its Lorentz projection plus an exactly orthogonal positive support then has squared norm strictly greater than2, impossible for a norm2cap.',
        scope='All real cap directions and all Coxeter labels, conditional on coverage by precisely the twelve frozen old-wall balls. This excludes the all-high branch of raw attachment with the source perpendicular sidewalls. Additional nonperpendicular covering walls, changed old walls, and other H8 sources are outside the theorem. No general H9 nonexistence claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(result['three_pair_obstructions']))
