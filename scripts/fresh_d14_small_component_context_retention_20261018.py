#!/usr/bin/env python3
"""Retain D14 long-leaf halves only with a whole A1, A2 or A3 witness.

The report proves retention through the complete current chosen-output
routes. The finite calculation refines only the two affected residual-zero
classes of the frozen complete count, using pointwise monotonicity elsewhere.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement,product
import fresh_rank15_a1_context_composed_ambient58_frontier_20261018 as G
import fresh_d14_long_leaf_a16_half_selected_transfer_20261018 as S
import fresh_order11_bcd_half_pair_transfer_20261018 as K
import fresh_d12_radical_saturated_half_replacements_20261018 as P
import fresh_d12_remaining_selected_half_replacements_20261018 as N

ROOT=G.ROOT;C=G.C;L=G.L;T=S.T;A=S.A;D=S.D
TAG='fresh_d14_small_component_context_retention_20261018'
JOINT='fresh_d14_a16_d15_joint_face_composition_20261018'
KD='fresh_affine_d13_long_pair_transfer_20261017'
SPEC='fresh_spectator_energy_below_one_20261017'
LITERAL=(T.A.OVER,T.A.OLD16,T.A.OLD17,T.A.D13,T.A.D14,T.T.PARTIAL,T.D.TAG,T.T.NEW_D)
STEMS=tuple(dict.fromkeys((G.TAG,C.TAG,L.TAG,G.W.TAG,G.Q.TAG,G.R.TAG,S.TAG,JOINT,
 K.TAG if hasattr(K,'TAG')else 'fresh_order11_bcd_half_pair_transfer_20261018',
 KD,SPEC,P.TAG,N.TAG,*LITERAL)))
FINAL_REPORTS={G.TAG:'f0d7201ff57e4c23969b20db72acfce6707457124f45babedbddbd2025092bd9',
 C.TAG:'37ead4131f4676475d44bb92fd621eba68a4859c9715383d37cbba2cfaac2123',
 L.TAG:'2701281f8d5d3660a7ff792c5829c59c8e30b5f168ddea1aa33d509def66271f'}
SMALL={('A',r)for r in(1,2,3)}
SMALL_NAMES={'A1','A2','A3'}
enc=S.enc

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return G.actual_base(t)-(64 if t==('D',14)else 0)
actual_cross=G.actual_cross
def context(fs):
    co=Counter(fs)
    return G.context(fs)+(64*co['D',14]if any(co[t]for t in SMALL)else 0)
def value(fs):
    out=sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)
    assert out==G.value(fs)-(0 if any(t in SMALL for t in fs)else 64*fs.count(('D',14)))
    return out

def witness(blocks,marks):
    out=[b for b in blocks if b['type']in SMALL_NAMES and not set(marks)&set(b['nodes'])]
    return [dict(type=b['type'],nodes=b['nodes'])for b in out]

def selected(docs):
    current={(p['index'],p['extra_contact']):p for p in docs[G.TAG]['local_proof']['whole_variants']}
    rows=[]
    for p in docs[S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            g=A.matrix(v['gram']);now=current[p['index'],v['extra_contact']];hits=[]
            for c,paid in zip(v['corners'],now['corners'],strict=True):
                assert c['marks']==paid['marks']
                ids=sorted(set(v['white'])|set(c['marks']))
                blocks=D.finite_blocks(g,ids);assert blocks==c['components']
                aa,bb=[next(b for b in blocks if m in b['nodes'])for m in c['marks']]
                if aa is bb and aa['type']=='D14':
                    images=T.mark_pairs(g,aa['nodes'],c['marks'],14,'D')
                    if set(images)&S.CHANGES:
                        assert set(images)<=S.CHANGES
                        ws=witness(blocks,c['marks']);assert ws and all(w['type']in('A1','A2')for w in ws)
                        assert F(paid['joint_lower'])==F(1,2)
                        hits.append(dict(marks=c['marks'],images=images,whole_witnesses=ws,retained_half='1/2'))
            assert len(hits)==1
            assert F(now['total'])>=1
            rows.append(dict(index=p['index'],extra_contact=v['extra_contact'],guarded_outputs=hits,current_total=now['total']))
    assert len(rows)==24 and len(current)==24
    assert docs[S.TAG]['complete_patterns']==600 and docs[S.TAG]['complete_maps']==1104
    assert docs[S.TAG]['affected_contexts']==17 and docs[S.TAG]['affected_maps']==36
    assert docs[S.TAG]['proper_projection_exclusions']==4
    return rows

def historical(docs):
    old=docs[S.TAG]['literal_interfaces']['records'];cur=docs[L.TAG]['current_interfaces']['literal_rows'];rows=[];hits=[]
    for s,l,(stem,case,g,white,cs)in zip(old,cur,T.literal_inputs(docs),strict=True):
        assert(s['source'],s['case'])==(l['source'],l['case'])==(stem,enc(case))
        assert len(s['corners'])==len(l['corners'])==len(cs)
        guarded=[]
        for a,b,(marks,*_)in zip(s['corners'],l['corners'],cs,strict=True):
            assert a['marks']==b['marks']==marks
            if a['debit']=='0':continue
            assert stem==KD and set(map(tuple,a['D14_images']))<=S.CHANGES
            blocks=D.finite_blocks(g,sorted(set(white)|set(marks)))
            assert blocks==b['components']
            ws=witness(blocks,marks);assert len(ws)==1 and ws[0]['type']=='A3'
            shared=next(z for z in blocks if set(marks)<=set(z['nodes']))
            images=T.mark_pairs(g,shared['nodes'],marks,14,'D');assert images==list(map(tuple,a['D14_images']))
            guarded.append(dict(marks=marks,images=images,whole_witnesses=ws,retained_lower=b['lower']))
        assert sum(F(c['lower'])for c in l['corners'])==F(l['total'])>=1
        rows.append(dict(source=stem,case=case,guarded_outputs=guarded,retained_total=l['total']))
        hits.extend((stem,case,r)for r in guarded)
    assert len(rows)==216 and len(hits)==16
    return dict(records=rows,guarded_outputs=16,all_guarded_outputs_from=KD,
        arbitrary_support_input='Complete affine-D13 energy3/4 whole-A3 sealing, not saved absence')

def order11_boundary():
    # The only newly constrained output has the chosen side A3 energy3/4.
    # The total is3/2, and a second active component costs at least1/2.
    # With disjoint unmixed supports the opposite side is therefore exactly
    # one whole A3 at its endpoint (the inherited strict-below-one theorem).
    assert F(3,2)-F(3,4)==F(3,4)and 2*F(1,2)>F(3,4)
    f,_,(a,b)=K.affine('D');adj=K.graph(f);rows=[]
    for reverse in(False,True):
        old,far=(b,a)if reverse else(a,b)
        near_center=next(iter(adj[old]));far_center=next(iter(adj[far]))
        mate=next(i for i in adj[near_center]if i!=old and len(adj[i])==1)
        far_mate=next(i for i in adj[far_center]if i!=far and len(adj[i])==1)
        for far_attach in(far,far_mate):
            g=2*K.s.eye(19);g[:11,:11]=f
            # c=11,d=12, whole spectator A3 chains13..15 and16..18.
            for i,j in((11,mate),(12,far_attach),(11,13),(13,14),(14,15),(12,16),(16,17),(17,18)):
                g[i,j]=g[j,i]=-1
            white=[i for i in range(19)if i not in(old,far,11,12)]
            assert K.inertia(g)==[17,1,1]
            K.pd(g.extract(white,white))
            for black in(old,far,11,12):K.pd(g.extract(white+[black],white+[black]))
            ids=sorted(white+[11,far]);K.pd(g.extract(ids,ids))
            gg=[[F(g[i,j])for j in range(19)]for i in range(19)]
            blocks=D.finite_blocks(gg,ids);shared=next(z for z in blocks if {11,far}<=set(z['nodes']))
            assert shared['type']=='D14'
            images=T.mark_pairs(gg,shared['nodes'],[11,far],14,'D')
            assert images==[(3,12),(3,13)]and not set(images)&S.CHANGES
            ws=witness(blocks,[11,far]);assert len(ws)==1 and ws[0]['type']=='A3'and ws[0]['nodes']==[16,17,18]
            rows.append(dict(reverse=reverse,far_attachment=far_attach,gram=gg,white=white,
                marks=[11,far],components=blocks,images=images,whole_witnesses=ws,
                selected_ideal_half='1/2',distinct_unchanged_finite_unit='1',used_as_half='1/2',total='1'))
    return dict(chosen_side_energy='3/4',opposite_side_energy='3/4',
        opposite_whole_factor='A3',canonical_cases=rows,marked_distance=9,
        no_changed_D14_role=True,all_small_A3_automorphisms_inherited=True)

def d12_outputs(docs):
    out=[]
    def walk(x,ix=None,path=''):
        if isinstance(x,dict):
            ix=x.get('index',ix)
            if x.get('reason')=='retained_whole_D14_weight':
                images=list(map(tuple,x['images']));assert not set(images)&S.CHANGES
                out.append(dict(index=ix,path=path,marks=x['marks'],images=images,payment=x['payment']))
            for k,v in x.items():walk(v,ix,path+'/'+k)
        elif isinstance(x,list):
            for i,v in enumerate(x):walk(v,ix,path+'/'+str(i))
    for tag,key in((P.TAG,'proofs'),(N.TAG,'new_proofs')):walk(docs[tag][key],path=tag+'/'+key)
    assert len(out)==16
    assert {p for c in out for p in c['images']}=={(0,11),(2,12),(2,13)}
    return out

def current_outputs(docs):
    later=docs[S.TAG]['later_finite_interfaces'];assert len(later)==22
    assert sum(len(p['corners'])for p in later)==564
    assert all(not set(map(tuple,c['D14_images']))&S.CHANGES for p in later for c in p['corners'])
    positive=[]
    for p in docs[L.TAG]['local_proofs']['whole_proofs']:
        for c in p['corners']:
            if F(c['payment'])<=0:continue
            aa,bb=[next(b for b in c['components']if m in b['nodes'])for m in c['marks']]
            assert not(aa is bb and aa['type']=='D14')
            positive.append(dict(index=p['index'],marks=c['marks'],marked_types=[aa['type'],bb['type']]))
    assert len(positive)==12
    return dict(inherited_later_corners=564,inherited_later_target_hits=0,new_rank16_positive_outputs=positive,
        rank15_guard_already_checked_on_all24_variants=True,
        unchanged_rank13_cross_and_A1_D5_chosen_outputs=True)

def finite(docs):
    ref=docs[C.TAG]['actual_typed_refinement'];ff=docs[C.TAG]['finite'];fibers={}
    for t in C.TYPES:fibers.setdefault(C.representative(t),[]).append(t)
    old_above=[r for r in ref['classes']if r['actual_obstruction_retained']]
    assert len(old_above)==9
    affected=[r for r in old_above if r['counts'][5]>0 and r['counts'][6]!=1]
    assert len(affected)==2 and all(r['residual_rank']==0 for r in affected)
    rows=[]
    for r in affected:
        key=(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count'])
        u=next(z for z in ff['above_curvature']if(tuple(z['counts']),z['fixed_A5_count'],z['fixed_A12_count'])==key)
        reps=C.fixed_components(tuple(r['counts']),u['key'])+[('A',5)]*r['fixed_A5_count']
        expected=sorted(tuple(sorted(t for group in choice for t in group))for choice in product(*(
            combinations_with_replacement(fibers[t],n)for t,n in sorted(Counter(reps).items()))))
        assert expected==sorted(tuple(map(tuple,z['fixed_actual']))for z in r['rows'])
        costs=[]
        for z in r['rows']:
            fs=list(map(tuple,z['fixed_actual']));assert sum(n for _,n in fs)==58
            assert z['actual_A1_count']==z['actual_A2_count']==0 and z['values']==[0]
            assert not any(t in SMALL for t in fs)
            costs.append(dict(forest=fs,previous_cost=F(G.value(fs),64),new_cost=F(value(fs),64)))
        oldmax=max(z['previous_cost']for z in costs);newmax=max(z['new_cost']for z in costs)
        rows.append(dict(counts=r['counts'],previous_maximum=oldmax,new_maximum=newmax,
            complete_actual_assignments=costs,maximizers=[z['forest']for z in costs if z['new_cost']==newmax]))
    assert sorted((r['previous_maximum'],r['new_maximum'])for r in rows)==[(F(378),F(753,2)),(F(379),F(378))]
    # G already made the third rank14 class at most377. Every other above
    # class contains no rank14 component, so its old attaining witness stays.
    unchanged=[]
    for r in old_above:
        if r['counts'][5]>0:continue
        assert all(('D',14)not in list(map(tuple,fs))for fs in r['maximizers'])
        vals=[F(value(list(map(tuple,fs))),64)for fs in r['maximizers']]
        mx=F(r['actual_maximum'])-F(1,2)*r['counts'][6];assert max(vals)==mx
        unchanged.append(dict(counts=r['counts'],maximum=mx,witness=r['maximizers'][vals.index(mx)]))
    # One of these six rows was changed by G: A13+3rank15, now378.
    for row in unchanged:
        if row['counts'][6]:
            assert row['maximum']==F(378)
            row['witness']=[['A',13],['A',15],['A',15],['A',15]]
            assert F(value(list(map(tuple,row['witness']))),64)==378
    assert len(unchanged)==6
    maximum=max([r['new_maximum']for r in rows]+[r['maximum']for r in unchanged])
    witness0=[('D',17),('D',17),('D',24)]
    assert maximum==F(1519,4)and F(value(witness0),64)==maximum
    assert sum(r['new_maximum']>377 for r in rows)+len(unchanged)==7
    return dict(previous_exact_maximum='1519/4',exact_actual_maximum=maximum,
        global_attaining_forest=witness0,previous_above_curvature_classes=8,new_above_curvature_classes=7,
        changed_complete_classes=rows,unchanged_above_curvature_classes=unchanged,
        complete_predecessor_upper_classes=39,global_dynamic_program_repeated=False,
        pointwise_nonincrease_covers_all_other_classes=True)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    for s,h in FINAL_REPORTS.items():assert digest(f'results/{s}.md')==h
    selected_rows=selected(docs);hist=historical(docs);boundary=order11_boundary()
    d12=d12_outputs(docs);current=current_outputs(docs);count=finite(docs)
    for t in C.TYPES:
        assert actual_base(t)>=0
        assert G.actual_base(t)-actual_base(t)==(64 if t==('D',14)else 0)
    assert F(actual_base(('D',14)),64)==F(177,2)
    assert docs[G.TAG]['actual_cusp_surplus']=='36'and docs[G.TAG]['positive_ideal_profiles']==2552
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_D14_SMALL_COMPONENT_CONTEXT_RETENTION_AND_MONOTONE_FRONTIER',
        sources={p:digest(p)for p in sorted(paths)},baseline=G.TAG,
        guarded_pairs=sorted(S.CHANGES),retained_weight='1/2',unguarded_weight='0',
        guard='At least one whole unmarked A1, A2 or A3 component in the same actual finite corner',
        ordinary_D14_cost='177/2',context_scale64='64*count(D14) if A1/A2/A3 is present, otherwise0',
        witness_multiplicity_not_counted=True,selected_whole_variants=selected_rows,
        historical_interfaces=hist,coupled_order11_boundary=boundary,
        unchanged_D12_D14_outputs=d12,current_interfaces=current,finite=count,
        ideal_weights_changed=False,ideal_verification='Frozen G/C envelope and positive low8 extension inherited; no polynomial or convolution replay',
        positive_ideal_profiles=2552,actual_cusp_surplus='36',finite_curvature='377',strict_global_constant='4205',
        necessary_ambient58_inequality='(11/4)*v_f >= 4205 + 36*v_i',minimum_finite_vertices=1530,
        completed_dimension_bound=58,new_dimension_bound=None,
        scope='Complete current chosen-output retention with inherited arbitrary-support routing/sealing; whole-A3 equality witness, explicit unchanged coupled-order11 D14 marks, and exact residual-zero type refinements. No global DP or new dimension theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['exact_actual_maximum'],out['finite']['new_above_curvature_classes'],flush=True)
