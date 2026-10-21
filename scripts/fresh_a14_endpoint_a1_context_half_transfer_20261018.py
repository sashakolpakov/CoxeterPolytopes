#!/usr/bin/env python3
"""An A14 endpoint half only with a whole unmarked A1.

Complete selected support and whole-core projection seal, retained current
outputs, and monotone residual-zero type refinement. No global DP.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
import fresh_d14_small_component_context_retention_20261018 as H
import fresh_finite_a14_near_endpoint_zero_transfer_20261018 as V

ROOT=H.ROOT;G=H.G;C=H.C;L=H.L;S=H.S;B=S.G;A=S.A;T=S.T;Q=S.Q;D=S.D
TAG='fresh_a14_endpoint_a1_context_half_transfer_20261018'
PAIR=(0,13)
ROUTES=('fresh_finite8_cycle_replacement_20261017','fresh_finite8_ideal5_local_replacement_20261017',
 'fresh_cycle12_ideal5_transfer_20261017','fresh_cycle13_ideal5_transfer_20261017','fresh_cycle14_ideal5_transfer_20261017',
 'fresh_affine_d14_long_pair_transfer_20261017','fresh_depth6_finite6_zero_transfer_20261018',
 'fresh_depth6_finite7_adaptive_transfer_20261018','fresh_depth6_finite8_zero_transfer_20261018',
 'fresh_order10_d_opposite_fork_fractional_transfer_20261018',
 'fresh_order10_cycle_half_transfer_20261018')
STEMS=tuple(dict.fromkeys((H.TAG,G.TAG,C.TAG,L.TAG,G.Q.TAG,G.R.TAG,S.TAG,V.TAG,V.L.R.TAG,A.OVER,
 *B.FINITE,*ROUTES,H.P.TAG,H.N.TAG)))
FINAL_REPORTS={H.TAG:'10fb43eb7646aeb36b1f14c776d3defc76a6f63923d8d96c2f29bbbe01170f02',
 G.TAG:'f0d7201ff57e4c23969b20db72acfce6707457124f45babedbddbd2025092bd9'}
enc=S.enc

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return H.actual_base(t)-(64 if t==('A',14)else 0)
actual_cross=H.actual_cross
def context(fs):return H.context(fs)+(32*fs.count(('A',14))if('A',1)in fs else 0)
def value(fs):
    z=sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)
    assert z==H.value(fs)-(32 if('A',1)in fs else 64)*fs.count(('A',14))
    return z

def endpoint(g,blocks,marks):
    shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
    if shared is None or shared['type']!='A14':return False
    images=T.mark_pairs(g,shared['nodes'],marks,14,'A')
    if PAIR in images:assert images==[PAIR];return True
    return False

def local(docs):
    over=docs[A.OVER];rec=next(r for r in docs[V.L.R.TAG]['records']if r['kind']=='A')
    patterns=[c for c in over['isolation']['cases']if(c.get('component_order_bound')or 0)>=14 and len(c['marked_component'])<=14]
    assert len(patterns)==600 and [c['index']for c in patterns]==[r['index']for r in rec['selected_patterns']]
    rows=[];hits=[];maps_count=0
    for c,old in zip(patterns,rec['selected_patterns'],strict=True):
        core=over['isolation']['cores'][c['core_key']];g=A.matrix(core['gram'])
        maps=T.T.mapped(g,c['marked_component'],c['preferred_corner'],rec['target'])
        assert maps==old['maps'];maps_count+=len(maps)
        affected=[m for m in maps if tuple(m['marked_images'])==PAIR]
        rows.append(dict(index=c['index'],maps=maps,affected_maps=affected))
        if affected:hits.append(c)
    assert maps_count==552 and [c['index']for c in hits]==[1639]
    assert sum(len(r['affected_maps'])for r in rows)==2
    assert PAIR not in list(map(tuple,rec['retained_basic_pairs']))
    c=hits[0];core=over['isolation']['cores'][c['core_key']]
    g,white,left,right=V.L.basic(c,core);gg,J,du,dv,end=V.X.constraints(c,core)
    assert gg==g and c['method']=='lorentz_residual'and len(g)==18 and len(c['marked_component'])==14
    assert Q.inertia(g)==[17,1,0]
    remaining=sorted(set(c['possible_contacts'])-set(c['marked_component']))
    assert remaining==[2,3,12,13]and[J[i][i]for i in remaining]==[14,4,14,4]
    assert all(J[i][j]>0 for e in c['end_sets']for i in e for j in e)
    assert end['both_ends_shortcut_bound']<=4
    corners=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        blocks=D.finite_blocks(g,ids);shared=next(b for b in blocks if set(marks)<=set(b['nodes']))
        images=T.mark_pairs(g,shared['nodes'],marks,len(shared['nodes']),shared['type'][0])
        a1=[b for b in blocks if b['type']=='A1'and not set(marks)&set(b['nodes'])]
        if marks==c['preferred_corner']:
            assert shared['type']=='A14'and images==[PAIR]and len(a1)==2
            payment=F(1,2);reason='guarded_A14_endpoint_half'
        elif shared['type']=='A15':
            assert set(images)=={(0,13),(1,14)}and len(a1)==1
            payment=F(1,4);reason='retained_guarded_A15_quarter'
        else:
            assert shared['type']=='A16'and images==[(1,14)]
            payment=F();reason='retained_A16_zero'
        corners.append(dict(marks=marks,components=blocks,images=images,whole_A1_witnesses=a1,payment=payment,reason=reason))
    assert sorted(x['payment']for x in corners)==[F(),F(1,4),F(1,4),F(1,2)]
    return dict(complete_patterns=600,complete_maps=maps_count,selected_patterns=rows,
        affected_contexts=[1639],affected_maps=2,gram=g,white=white,corner_marks=[r['marks']for r in corners],
        end_projection=end,remaining_contacts=remaining,remaining_inverse_diagonals=[[i,J[i][i]]for i in remaining],
        minimum_nonzero_same_end_projection='4',all_further_whites_core_orthogonal_unless_short4=True,
        corners=corners,total='1')

def interfaces(docs):
    # L reconstructed every finite component in the complete216 interface.
    literal=[]
    for row in docs[L.TAG]['current_interfaces']['literal_rows']:
        assert all(not any(b['type']=='A14'and set(c['marks'])<=set(b['nodes'])for b in c['components'])for c in row['corners'])
        literal.append(dict(source=row['source'],case=row['case'],changed_A14_outputs=0))
    assert len(literal)==216
    later=[]
    for stem,keys in B.FINITE.items():
        rows=[]
        for section in keys:
            for path,p,gg,white,blocks in B.corners(docs[stem][section],section):
                g=A.matrix(gg);assert sorted(A.components(g,sorted(set(white)|set(p['marks']))))==sorted(b['nodes']for b in blocks)
                assert not endpoint(g,blocks,p['marks']),(stem,path)
                rows.append(dict(path=path,marks=p['marks'],changed_A14_output=False))
        later.append(dict(source=stem,corners=rows))
    assert len(later)==22 and sum(len(r['corners'])for r in later)==564
    variants=[];changed=0
    current={(p['index'],p['extra_contact']):p for p in docs[G.TAG]['local_proof']['whole_variants']}
    retained={(p['index'],p['extra_contact']):p for p in docs[H.TAG]['selected_whole_variants']}
    for p in docs[S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            g=A.matrix(v['gram']);now=current[p['index'],v['extra_contact']]
            assert retained[p['index'],v['extra_contact']]['current_total']==now['total']
            entries=[]
            for old,before in zip(v['corners'],now['corners'],strict=True):
                assert old['marks']==before['marks']
                hit=endpoint(g,old['components'],old['marks'])
                a1=any(b['type']=='A1'for b in old['components'])
                debit=(F(1,2)if a1 else F(1))if hit else F()
                val=F(before['joint_lower'])-debit;assert val>=0
                if hit:assert old['reason']=='retained_A14_endpoint_unit';changed+=1
                entries.append(dict(marks=old['marks'],A14_endpoint=hit,whole_unmarked_A1=a1,
                    before=before['joint_lower'],debit=debit,retained_lower=val))
            total=sum(c['retained_lower']for c in entries);assert total>=1
            variants.append(dict(index=p['index'],extra_contact=v['extra_contact'],corners=entries,total=total))
    assert len(variants)==24 and changed==8
    # The original near-endpoint A14 theorem spends A15(1,13), not PAIR.
    olda=docs[V.TAG]
    # Explicit current positive-output role checks, beyond the older564 list.
    newer=[]
    for tag,key in((L.TAG,'local_proofs'),(G.R.TAG,'local_seals')):
        rows=[]
        for path,p,gg,white,blocks in B.corners(docs[tag][key],key):
            assert not endpoint(A.matrix(gg),blocks,p['marks'])
            rows.append(dict(path=path,changed_A14_output=False))
        assert rows;newer.append(dict(source=tag,corners=rows))
    d12=[]
    for tag in(H.P.TAG,H.N.TAG):
        rows=[(path,p)for path,p in G.W.walk(docs[tag])if F(p['payment'])>0]
        assert all(not any(b['type']=='A14'for b in p['roles'])for _,p in rows)
        d12.append(dict(source=tag,positive_outputs=len(rows),A14_outputs=0))
    assert sum(r['positive_outputs']for r in d12)==332
    return dict(historical_literals=literal,later_finite=later,whole_D14_variants=variants,
        changed_D14_alternative_corners=8,newer_seals=newer,chosen_D12_outputs=d12,
        near_endpoint_A14_transfer_retains_A15_unit=True)

def finite(docs):
    old=docs[H.TAG]['finite'];affected=[r for r in old['changed_complete_classes']if F(r['new_maximum'])>377]
    assert len(affected)==1
    rows=[]
    for r in affected:
        costs=[]
        for z in r['complete_actual_assignments']:
            fs=list(map(tuple,z['forest']));assert sum(n for _,n in fs)==58 and ('A',1)not in fs
            assert F(H.value(fs),64)==F(z['new_cost'])
            costs.append(dict(forest=fs,previous_cost=z['new_cost'],new_cost=F(value(fs),64)))
        maximum=max(z['new_cost']for z in costs);assert maximum==377
        rows.append(dict(counts=r['counts'],previous_maximum=r['new_maximum'],new_maximum=maximum,
            complete_actual_assignments=costs,maximizers=[z['forest']for z in costs if z['new_cost']==maximum]))
    unchanged=old['unchanged_above_curvature_classes'];assert len(unchanged)==6
    for r in unchanged:
        fs=list(map(tuple,r['witness']));assert ('A',14)not in fs
        assert F(value(fs),64)==F(r['maximum'])
    witness=[('D',17),('D',17),('D',24)]
    assert F(value(witness),64)==F(old['exact_actual_maximum'])==F(1519,4)
    return dict(previous_above_curvature_classes=7,new_above_curvature_classes=6,
        exact_actual_maximum='1519/4',global_attaining_forest=witness,
        changed_complete_classes=rows,unchanged_above_curvature_classes=unchanged,
        all_other_classes_bounded_by_pointwise_nonincrease=True,global_DP_repeated=False)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    for s,h in FINAL_REPORTS.items():assert digest(f'results/{s}.md')==h
    geometry=local(docs);inter=interfaces(docs);count=finite(docs)
    assert F(actual_base(('A',14)),64)==88
    assert all(actual_base(t)==H.actual_base(t)-(64 if t==('A',14)else 0)for t in C.TYPES)
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_A14_ENDPOINT_A1_CONTEXT_HALF_TRANSFER_AND_MONOTONE_FRONTIER',
        sources={p:digest(p)for p in sorted(paths)},baseline=H.TAG,
        changed_pair=PAIR,guarded_weight='1/2',unguarded_weight='0',guard='Whole unmarked A1 in the same finite corner',
        ordinary_A14_cost='88',scale64_context='32*count(A14) if count(A1)>0, otherwise0',
        witness_multiplicity_not_counted=True,local_proof=geometry,interfaces=inter,finite=count,
        all_cross_costs_unchanged=True,ideal_weights_changed=False,new_cross_coupons=[],
        triangle_minimum='2',quadrilateral_minimum='1',
        positive_ideal_profiles=2552,actual_cusp_surplus='36',finite_curvature='377',strict_global_constant='4205',
        necessary_ambient58_inequality='(11/4)*v_f >= 4205 + 36*v_i',minimum_finite_vertices=1530,
        completed_dimension_bound=58,new_dimension_bound=None,
        scope='Complete600/552 selected A14 search, sole full projection seal and two guarded A15 quarter alternatives; full current chosen-output retention, eight D14 alternative adjustments, and monotone complete residual-zero type refinement. No global DP or cusp replay.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['exact_actual_maximum'],out['finite']['new_above_curvature_classes'],flush=True)
