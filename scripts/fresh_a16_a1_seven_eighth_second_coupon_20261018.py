#!/usr/bin/env python3
"""Guarded A16 seven-eighths and an A1-guarded A4/A12 second coupon.

Two existing full radical seals prove the exchange. Existing complete finite
upper certificates suffice: no dynamic program or new obstruction census.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import fresh_d16_a1_long_leaf_d17_zero_ambient58_frontier_20261018 as H

ROOT=H.ROOT;R=H.R;N=H.N;C=H.C;L=H.L;B=H.B;W=H.W;T=H.T;A=H.A;Q=H.Q;D=H.D
TAG='fresh_a16_a1_seven_eighth_second_coupon_20261018'
STEMS=tuple(dict.fromkeys((H.TAG,*H.STEMS)))
enc=H.enc

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
actual_base=H.actual_base
actual_cross=H.actual_cross
def context(fs):
    co=Counter(fs)
    return H.context(fs)+(-16*co['A',16]+32*co['A',4]*co['A',12]if co['A',1]else 0)
def value(fs):return sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)

def change(g,white,marks):
    ids=sorted(set(white)|set(marks))
    if Q.inertia(A.sub(g,ids))!=[len(ids),0,0]:return dict(delta=F(),roles=[],components=[])
    blocks=D.finite_blocks(g,ids)
    witnesses=[b for b in blocks if b['type']=='A1'and not set(marks)&set(b['nodes'])]
    aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
    roles=[];delta=F()
    if aa is bb and aa['type']=='A16':
        im=T.mark_pairs(g,aa['nodes'],marks,16,'A')
        if set(im)&L.ACHANGE:
            assert set(im)<=L.ACHANGE
            old=F(1)if witnesses else F();new=F(7,8)if witnesses else F()
            delta=new-old
            roles.append(dict(role='A16_guarded_seven_eighths',images=im,old=old,new=new,whole_A1=witnesses))
    elif aa is not bb:
        for sm,la,m,lm in((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
            if(sm['type'],la['type'])!=('A4','A12'):continue
            si=T.mark_pairs(g,sm['nodes'],[m],4,'A');li=T.mark_pairs(g,la['nodes'],[lm],12,'A')
            if si==[(0,),(3,)]and li==[(1,),(10,)]:
                delta=F(1,8)if witnesses else F()
                roles.append(dict(role='A4end_A12second_guarded_eighth',small_images=si,large_images=li,old=F(),new=delta,whole_A1=witnesses))
    assert len(roles)<=1
    return dict(delta=delta,roles=roles,components=blocks)

def local(docs):
    old=docs[L.TAG]['local_proofs'];sup=old['domains']['A']['complete_selected_support']
    assert len(sup['rows'])==394 and sum(len(r['maps'])for r in sup['rows'])==238
    assert old['domains']['A']['affected_contexts']==L.ACASES
    over=docs[A.OVER];rows=[]
    for index in L.ACASES:
        c=over['isolation']['cases'][index];p=L.a_seal(c,over['isolation']['cores'][c['core_key']])
        assert enc(p)==next(x for x in old['whole_proofs']if x['index']==index)
        current=next(x for x in docs[R.TAG]['local_proof']['whole_proofs']if x['index']==index)
        corners=[]
        for q in current['corners']:
            ch=change(A.matrix(p['gram']),p['white'],q['marks']);new=F(q['payment'])+ch['delta']
            if ch['roles']:
                assert len(ch['roles'][0]['whole_A1'])==1
            corners.append(dict(marks=q['marks'],previous=q['payment'],**ch,payment=new))
        assert sorted(x['payment']for x in corners)==[0,0,F(1,8),F(7,8)]
        chosen=[x for x in corners if x['payment']]
        assert chosen[0]['marks']!=chosen[1]['marks']
        assert chosen[0]['roles'][0]['whole_A1']==chosen[1]['roles'][0]['whole_A1']
        rows.append(dict(index=index,gram=p['gram'],white=p['white'],full_seal=p,corners=corners,total='1'))
    return dict(complete_A16_patterns=394,complete_A16_maps=238,affected_contexts=L.ACASES,
        affected_maps=4,whole_proofs=rows,changed_A16_pairs=sorted(L.ACHANGE),
        A16_ordinary_cost=F(actual_base(('A',16)),64),A16_cost_with_A1='110',
        A4_A12_ordinary_cost='8',A4_A12_cost_with_A1='17/2',
        coupon_weight='1/8',coupon_component_pair_increment='1/2',guard_counts_presence_only=True)

def interfaces(docs):
    literals=[]
    for old,(stem,case,g,white,cs)in zip(docs[H.TAG]['current_interfaces']['literals'],T.literal_inputs(docs),strict=True):
        assert(old['source'],old['case'])==(stem,enc(case));corners=[]
        for q,(marks,*_)in zip(old['corners'],cs,strict=True):
            assert q['marks']==marks;ch=change(g,white,marks)
            # All eight earlier A16 obligations already use the retained
            # whole D13 unit; their A16 components have no A1 and remain0.
            assert ch['delta']>=0
            corners.append(dict(marks=marks,previous=q['lower'],**ch,lower=F(q['lower'])+ch['delta']))
        assert sum(q['lower']for q in corners)>=1
        literals.append(dict(source=stem,case=case,corners=corners,total=sum(q['lower']for q in corners)))
    assert len(literals)==216
    later=[];count=0
    for stem,keys in B.FINITE.items():
        rows=[]
        for section in keys:
            for path,p,g,white,blocks in B.corners(docs[stem][section],section):
                ch=change(A.matrix(g),white,p['marks']);assert ch['delta']>=0
                rows.append(dict(path=path,marks=p['marks'],delta=ch['delta'],roles=ch['roles']));count+=1
        later.append(dict(source=stem,corners=rows))
    assert count==564
    newer=[]
    collections=[(R.TAG,docs[R.TAG]['local_proof']['whole_proofs']),
        (H.TAG,docs[H.TAG]['local_proof']['whole_proofs']),
        (N.TAG,docs[N.TAG]['local_proof']['proofs']),
        (R.R.TAG,docs[R.R.TAG]['local_seals']),
        (H.Z.U.TAG,[docs[H.Z.U.TAG]['local_proof']])]
    for stem,proofs in collections:
        for p in proofs:
            rows=[];physical=p.get('old_complete_seal',p)
            for q in p['corners']:
                ch=change(A.matrix(physical['gram']),physical['white'],q['marks'])
                if ch['delta']<0:assert stem==R.TAG and p['index']in L.ACASES
                rows.append(dict(marks=q['marks'],delta=ch['delta'],roles=ch['roles']))
            newer.append(dict(source=stem,index=p.get('index',p.get('name','A14')),corners=rows))
    variants=[]
    for p in docs[L.C.R.S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            rows=[]
            for q in v['corners']:
                ch=change(A.matrix(v['gram']),v['white'],q['marks']);assert ch['delta']>=0
                rows.append(dict(marks=q['marks'],delta=ch['delta']))
            variants.append(dict(index=p['index'],extra_contact=v['extra_contact'],corners=rows))
    assert len(variants)==24
    d12=[]
    for stem,expected in((B.PART,280),(B.REST,52)):
        ps=[(path,p)for path,p in W.walk(docs[stem])if F(p['payment'])>0];assert len(ps)==expected
        for _,p in ps:assert not any(b['type']=='A16'for b in p['roles'])
        d12.append(dict(source=stem,positive_outputs=expected,no_A16_output=True))
    return dict(literals=literals,later_packets=later,later_corners=count,newer_seals=newer,
        retained_D14_variants=variants,chosen_D12_outputs=d12,
        sole_decreased_chosen_outputs='The two own A16 selected units, replaced by distinct7/8+1/8 corners.')

def finite(docs):
    old=docs[R.TAG]['increasing_family_upper'];assert len(old['rows'])==539
    rows=[]
    for r in old['rows']:
        upper=F(r['upper'])-F(1,4)*r['rank16_count']
        rows.append(dict(counts=r['counts'],fixed_A5_count=r['fixed_A5_count'],
            rank16_count=r['rank16_count'],upper=upper))
    high=[r for r in rows if r['upper']>377]
    assert len(high)==1 and high[0]['rank16_count']==2 and high[0]['upper']==378
    gaps=[]
    for kind in'ABD':
        t=(kind,16)
        guarded=actual_base(t)+(112 if kind=='A'else 64 if kind=='D'else 0)
        gap=C.upper_base(('D',16))+32-guarded+C.upper_cross(('A',1),('D',16))-actual_cross(('A',1),t)
        assert gap>=32;gaps.append(dict(type=t,guarded_base=F(guarded,64),one_A1_gap=F(gap,64)))
    assert high[0]['upper']-2*min(x['one_A1_gap']for x in gaps)==377
    # Both C's representative allowance and N's old actual A4/A12 cost
    # already include8.5, so the guarded new8.5 needs no altered table.
    assert F(actual_cross(('A',4),('A',12)),64)==8
    assert F(C.upper_cross(('A',4),('A',12)),64)==F(17,2)
    assert F(N.actual_cross(('A',4),('A',12)),64)==F(17,2)
    for a in C.TYPES:
        for b in C.TYPES:
            guarded=actual_cross(a,b)+(32 if L.pair(a,b)else 0)
            assert guarded<=N.actual_cross(a,b)
            assert guarded<=C.upper_cross(C.representative(a),C.representative(b))
    overlap=F(docs[H.TAG]['finite']['inherited_D17_increasing_family_upper']);assert overlap<377
    # The same N-complement base drops and all-residual table remain valid:
    # only a nonnegative A16 guard shrinks, and all new coupons are<=N.
    complement=docs[H.TAG]['finite']['complement_classes'];assert len(complement)==39
    assert sum(len(r['actual_fixed_fibers'])for r in complement)==748
    tab=docs[H.TAG]['finite']['residual_table'];types=list(map(tuple,tab['types']))
    assert ('A',16)not in types and not {('A',4),('A',12)}<=set(types)
    assert all(actual_base(t)==w for t,w in zip(types,tab['weights'],strict=True))
    assert all(context([t])==H.context([t])for t in types)
    assert F(tab['values'][58],64)==378
    witnesses=[[('D',17),('D',17),('D',24)],[('A',13)]+[('A',15)]*3]
    assert all(value(fs)==64*378 for fs in witnesses)
    plateau=[('A',16),('A',16),('A',1),('D',25)]
    assert F(H.value(plateau),64)==F(755,2)and F(value(plateau),64)==377
    return dict(complete_conditional_classes=539,conditional_rows=rows,one_A1_type_gaps=gaps,
        rank16_with_A1_upper='377',D17_increasing_family_upper=overlap,
        inherited_complete_complement_classes=39,inherited_actual_fixed_assignments=748,
        inherited_all_residual_table=H.TAG,exact_actual_maximum='378',actual_attainers=witnesses,
        removed_plateau=plateau,removed_plateau_old_cost='755/2',removed_plateau_new_cost='377',
        new_DP=False,new_obstruction_census=False)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    proof=local(docs);inter=interfaces(docs);count=finite(docs)
    assert F(actual_base(('A',16)),64)==F(433,4)
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_A16_A1_SEVEN_EIGHTH_SECOND_COUPON_AMBIENT58',sources={p:digest(p)for p in sorted(paths)},
        baseline=H.TAG,scale=64,local_proof=proof,current_interfaces=inter,finite=count,
        ideal_weights_changed=False,positive_ideal_profiles=2552,actual_cusp_surplus='36',
        finite_curvature='377',strict_global_constant='4205',triangle_minimum='2',quadrilateral_minimum='1',
        necessary_ambient58_inequality='v_f >= 4205 + 36*v_i',minimum_finite_vertices=4205,
        completed_dimension_bound=58,new_dimension_bound=None,
        scope='Two full radical seals and complete chosen-output retention. Inherited complete conditional/fiber bounds and unchanged residual table; no new dynamic program, convolution or obstruction census.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['exact_actual_maximum'],flush=True)
