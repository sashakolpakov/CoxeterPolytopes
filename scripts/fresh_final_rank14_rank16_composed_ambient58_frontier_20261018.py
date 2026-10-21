#!/usr/bin/env python3
"""Final monotone A14 composition with the frozen sharp rank16 frontier.

Physical chosen-output compatibility, exact cost identities and unchanged
attainment; no new dynamic program or obstruction-class census.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
import fresh_rank16_a1_trade_composed_ambient58_frontier_20261018 as R
import fresh_a14_endpoint_a1_context_half_transfer_20261018 as U

ROOT=R.ROOT;C=R.C;A=R.A;D=R.D
TAG='fresh_final_rank14_rank16_composed_ambient58_frontier_20261018'
STEMS=(R.TAG,U.TAG,R.N.TAG,R.K.TAG,R.G.TAG,R.L.TAG,U.S.TAG)
PINS={R.TAG:('b8a988f7c37bb1f09f908d576d06f94b4c64421295c1e26d3058d9246f50abc6','e634e3a45f5cd65481bdaabb715e7926e6b952f6381bf3ecd45da1b94f479e14','d2690c302c0e69bc2bf02390f816aa011fda52d982bf088bcd6d9005ccdc81ff'),
 U.TAG:('e283dd353991a13108df4d81d34c0b89f02b68cbefafe39b2f8d1849bffdbdb8','14af340c03fe54cdd6a39c4110e45b701cb2839f941539c3fac116097777e6da','bdfacca37a915bd77a234ed2dd065971d0d1a4ebf1048f0fd520e37618af49dc')}

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return R.actual_base(t)-(64 if t==('A',14)else 0)
actual_cross=R.actual_cross
def context(fs):return R.context(fs)+(32*fs.count(('A',14))if('A',1)in fs else 0)
def value(fs):
    z=sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)
    assert z==R.value(fs)-(32 if('A',1)in fs else 64)*fs.count(('A',14))
    return z

def compatibility(docs):
    old=docs[U.TAG];p=old['local_proof'];g=A.matrix(p['gram']);selected=[]
    assert p['affected_contexts']==[1639]and p['complete_patterns']==600 and p['complete_maps']==552
    for c in p['corners']:
        blocks=D.finite_blocks(g,sorted(set(p['white'])|set(c['marks'])))
        assert blocks==c['components'];ch=R.changes(g,p['white'],c['marks']);assert ch['delta']==0
        if c['reason']=='retained_A16_zero':assert c['images']==[[1,14]]and F(c['payment'])==0
        selected.append(dict(marks=c['marks'],components=blocks,payment=c['payment'],rank16_delta=ch['delta']))
    assert sorted(F(c['payment'])for c in selected)==[0,F(1,4),F(1,4),F(1,2)]
    seals=[];positive=0
    for p in docs[R.TAG]['local_proof']['whole_proofs']:
        g=A.matrix(p['gram']);rows=[]
        for c in p['corners']:
            blocks=D.finite_blocks(g,sorted(set(p['white'])|set(c['marks'])))
            assert blocks==c['components']
            hit=U.endpoint(g,blocks,c['marks'])
            if F(c['payment'])>0:assert not hit;positive+=1
            rows.append(dict(marks=c['marks'],payment=c['payment'],A14_endpoint=hit))
        assert sum(F(c['payment'])for c in rows)==1
        seals.append(dict(index=p['index'],corners=rows))
    assert len(seals)==positive==6
    variants=[];current={(p['index'],p['extra_contact']):p for p in old['interfaces']['whole_D14_variants']}
    r16={(p['index'],p['extra_contact']):p for p in docs[R.TAG]['current_interfaces']['retained_quarter_variants']}
    for p in docs[U.S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            key=p['index'],v['extra_contact'];now=current[key];rr=r16[key];g=A.matrix(v['gram']);rows=[]
            for c,paid,prev in zip(v['corners'],now['corners'],rr['corners'],strict=True):
                assert c['marks']==paid['marks']==prev['marks']
                blocks=D.finite_blocks(g,sorted(set(v['white'])|set(c['marks'])))
                assert blocks==c['components']
                assert R.changes(g,v['white'],c['marks'])['delta']==F(prev['delta'])==0
                assert U.endpoint(g,blocks,c['marks'])==paid['A14_endpoint']
                rows.append(dict(marks=c['marks'],payment=paid['retained_lower'],rank16_delta='0'))
            total=sum(F(c['payment'])for c in rows);assert total==F(now['total'])>=1
            variants.append(dict(index=p['index'],extra_contact=v['extra_contact'],corners=rows,total=total))
    assert len(variants)==len(current)==len(r16)==24
    d17=[]
    for p in docs[R.N.TAG]['local_proof']['proofs']:
        g=A.matrix(p['gram']);rows=[]
        for c in p['corners']:
            blocks=D.finite_blocks(g,sorted(set(p['white'])|set(c['marks'])))
            assert blocks==c['components']
            assert not U.endpoint(g,blocks,c['marks'])and R.changes(g,p['white'],c['marks'])['delta']==0
            rows.append(dict(marks=c['marks'],payment=c['payment']))
        assert sum(F(c['payment'])for c in rows)==1
        d17.append(dict(index=p['index'],corners=rows))
    assert len(d17)==2
    return dict(A14_seal=selected,rank16_seals=seals,whole_D14_variants=variants,D17_seals=d17,
        triangle_minimum='2',quadrilateral_minimum='1',complete_chosen_output_composition=True)

def run():
    for s,hs in PINS.items():
        for (folder,ext),h in zip((('scripts','py'),('results','json'),('results','md')),hs):assert digest(f'{folder}/{s}.{ext}')==h
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    assert U.C is C and U.G is R.G and U.L is R.L and U.H is R.K
    inter=compatibility(docs)
    for t in C.TYPES:
        assert actual_base(t)>=0 and actual_base(t)==R.actual_base(t)-(64 if t==('A',14)else 0)
        for u in C.TYPES:assert actual_cross(t,u)==R.actual_cross(t,u)
    assert F(actual_base(('A',14)),64)==88
    old=docs[R.TAG];assert old['finite']['exact_actual_maximum']=='757/2'
    assert old['increasing_family_upper']['maximum']=='757/2'and old['finite']['complete_actual_fixed_assignments']==748
    witness=[('D',17),('D',17),('D',24)]
    assert value(witness)==R.value(witness)==64*F(757,2)and ('A',14)not in witness
    assert old['inherited_positive_ideal_profiles']==2552 and old['actual_cusp_surplus']=='36'
    assert docs[U.TAG]['ideal_weights_changed']is False and not docs[U.TAG]['new_cross_coupons']
    assert old['finite_curvature']=='377'and old['strict_global_constant']=='4205'
    excess=F(757,2)-377;ratio=F(4205)/excess;minimum=-(-ratio.numerator//ratio.denominator);assert minimum==2804
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return R.enc(dict(status='PASS_FINAL_RANK14_RANK16_COMPOSED_AMBIENT58_FRONTIER',sources={p:digest(p)for p in sorted(paths)},
        baseline=R.TAG,additional_local_rule=U.TAG,scale=64,physical_composition=inter,
        ordinary_A14_cost='88',A14_cost_with_A1='177/2',ordinary_base_identities=len(C.TYPES),cross_identities=len(C.TYPES)**2,
        cost_identity='new=R16-(1/2 if A1 present else1)*count(A14)',witness_multiplicity_not_counted=True,
        pointwise_monotone=True,exact_actual_maximum='757/2',actual_attaining_forest=witness,
        complete_upper_inherited=R.TAG,new_DP_or_obstruction_census=False,
        ideal_weights_changed=False,positive_ideal_profiles=2552,actual_cusp_surplus='36',
        finite_curvature='377',finite_excess='3/2',strict_global_constant='4205',
        necessary_ambient58_inequality='(3/2)*v_f >= 4205 + 36*v_i',minimum_finite_vertices=minimum,
        completed_dimension_bound=58,new_dimension_bound=None,new_geometric_example=False,
        scope='Physical composition of the new A14 seal, all24 D14 variants, six rank16 seals and two D17 seals; exact monotone inheritance of sharp378.5 with unchanged attainment. No dynamic program, new exact obstruction census or dimension57 theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['exact_actual_maximum'],out['minimum_finite_vertices'],flush=True)
