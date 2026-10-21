#!/usr/bin/env python3
"""Zero two affine D16 center halves using a witnessed A2/D5 half.

Reconstruct both full-core seals; check the current literal interface and
the complete finite contextual upper using unchanged frozen residual tables.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_order17_d5_d12_center_quarter_transfer_20261018 as Z
import fresh_order13_bc_inner_terminal_transfer_20261018 as H
import fresh_small_d16_coupon_and_ideal17_rebalance_20261018 as S

ROOT=H.ROOT;C=H.C;P=Z.P;Q=Z.Q;D=Z.D;L=Z.L;T=C.T;A=C.A
TAG='fresh_order17_center_contextual_zero_transfer_20261018'
CHANGE=Z.CHANGE;EXPECTED=Z.EXPECTED;N=59
STEMS=tuple(dict.fromkeys((*H.STEMS,H.TAG,*Z.STEMS,Z.TAG,S.TAG)))
FROZEN_H_MD='4354857b1e4ce2b5b2a67e58d31ccd9cbe1cecf944e1f94502f00007483b6ea0'
enc=H.enc


def contextual_roles(g,blocks,marks):
    selected=[next(b for b in blocks if m in b['nodes']) for m in marks]
    if selected[0] is selected[1]:return []
    witnesses={kind:[b for b in blocks if b['type']==kind and not set(marks)&set(b['nodes'])] for kind in ('A1','A9')}
    if not all(witnesses.values()):return []
    out=[]
    for sm,la,m,n in ((selected[0],selected[1],*marks),(selected[1],selected[0],marks[1],marks[0])):
        if sm['type']=='A2' and la['type']=='D5':
            a=T.mark_pairs(g,sm['nodes'],[m],2,'A');b=T.mark_pairs(g,la['nodes'],[n],5,'D')
            if a==[(0,),(1,)] and b==[(0,)]:
                out.append(dict(small=sm,large=la,marks=[m,n],canonical_images=[a,b],
                    whole_unmarked_witnesses=witnesses,payment='1/2',witness_multiplicity_not_counted=True))
    assert len(out)<=1
    return out


def local(docs):
    prior=docs[Z.TAG];oldcases={r['index']:r for r in docs[Z.O.TAG]['cases']};over=docs[Z.O.OVER]
    weights={tuple(r['pair']):F(r['weight']) for r in docs[S.TAG]['local_proofs']['weighted_pairs']}
    assert len(weights)==136 and sum(weights.values())==F(229,2)
    assert F(docs[C.TAG]['ideal']['current_envelopes']['16'])==F(229,34)
    for p in CHANGE:assert weights[p]==F(1,2);weights[p]=F()
    assert sum(weights.values())==F(227,2)
    target=prior['target'];autos=P.embeddings(target,target);_,ds,_=P.distance_data(target)
    assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGE}==CHANGE for mp in autos)
    assert all(ds[i][j]==10 for i,j in CHANGE)
    for field in ('retained_basic_pairs','retained_short8','retained_triangle_pairs'):
        assert all(weights[tuple(p)]==1 for p in prior[field])
    support=prior['reconstructed_support'];small=docs[Z.O.SMALL]
    patterns=[p for p in small['patterns'] if p['order_bound']>=17 and len(p['nodes'])<=17]
    rec=next(r for r in small['records'] if r['order']==17 and r['kind']=='D')
    saved={r['case']:r for r in rec['case_embeddings']};assert len(patterns)==140 and len(saved)==len(support)==71
    assert sum(len(r['embeddings']) for r in support)==488
    assert [r['index'] for r in support]==[r['case'] for r in rec['case_embeddings']]
    hits=[]
    for r in support:
        assert r['embeddings']==saved[r['index']]['embeddings'] and r['marked_images']==saved[r['index']]['marked_images']
        p=next(p for p in patterns if p['index']==r['index'])
        changed=[mp for mp in r['embeddings'] if tuple(sorted(mp[i] for i in p['marks'])) in CHANGE]
        if changed:
            assert len(changed)==8 and len(p['nodes'])==17
            hits.append(dict(index=r['index'],maps=changed))
    assert [r['index'] for r in hits]==EXPECTED and sum(len(r['maps']) for r in hits)==16
    proofs=[]
    for hit in hits:
        index=hit['index'];c=over['isolation']['cases'][index];core=over['isolation']['cores'][c['core_key']]
        seal=Z.inspect(c,core,oldcases[index])
        frozen=next(r for r in prior['selected_proofs'] if r['index']==index)
        assert seal=={k:v for k,v in frozen.items() if k!='maps'} and hit['maps']==frozen['maps']
        g=[[F(x) for x in row] for row in core['gram']];corners=[]
        for r in seal['corners']:
            if r['kind']=='affineD16_order17':
                corners.append(dict(marks=r['marks'],kind=r['kind'],payment=F(),changed_selected_pair=True));continue
            blocks=D.finite_blocks(g,sorted(set(seal['white'])|set(r['marks'])))
            assert blocks==r['components']
            new=contextual_roles(g,blocks,r['marks']);old=bool(r['coupon_roles'])
            assert not (old and new)
            if old:assert r['coupon_roles'][0]['small']['role']=='unique_long_arm_endpoint' and r['coupon_roles'][0]['large']['role']=='fork_center'
            corners.append(dict(marks=r['marks'],kind='finite',whole_components=blocks,
                old_center_coupon='1/2' if old else '0',new_contextual_roles=new,payment=F(bool(old or new),2),
                uses_finite_same_component_payment=False))
        paying=[r for r in corners if r['payment']]
        assert len(paying)==2 and paying[0]['marks']!=paying[1]['marks'] and sum(r['payment'] for r in corners)==1
        assert sum(bool(r.get('new_contextual_roles')) for r in corners)==1
        proofs.append(dict(index=index,inherited_full_seal=seal,maps=hit['maps'],corners=corners,total='1'))
    return dict(target=target,automorphisms=autos,old_cost='229/2',new_cost='227/2',saving='1',
        changed_pairs=sorted(CHANGE),old_weight='1/2',new_weight='0',
        weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
        complete_inherited_patterns=140,complete_inherited_maps=488,affected_contexts=EXPECTED,affected_maps=16,
        inherited_routing=prior['inherited_complete_high_overlap_routing'],historical_embedding_search_repeated=False,
        proofs=proofs,retained_basic_pairs=prior['retained_basic_pairs'],retained_short8=prior['retained_short8'],
        retained_triangle_pairs=prior['retained_triangle_pairs'])


def context_cost(fs,upper=False):
    co=Counter(fs)
    guard=co[('A',1)] and (any(k==9 and n for (typ,k),n in co.items()) if upper else co[('A',9)])
    return 16*co[('A',2)]*co[('A' if upper else 'D',5)] if guard else 0


def actual_value(fs):
    co=Counter(fs)
    return C.value(fs)+8*co[('A',1)]*co[('B',12)]+context_cost(fs)


def finite(docs):
    saved=docs[C.TAG]['finite'];old=C.finite(docs,saved)
    assert enc(old)==saved
    # The preceding order13 ordinary half is actual-only and fits exactly
    # the same simultaneous representative inequality as in its frozen proof.
    types=[t for t in C.P.Z.B.B.types() if t[1]<=N];rep=C.P.Z.B.B.representative
    for a in types:
        assert C.actual_base(a)<=C.upper_base(rep(a))
        for b in types:
            inc=8 if {a,b}=={('A',1),('B',12)} else 0
            assert C.actual_cross(a,b)+inc<=C.upper_cross(rep(a),rep(b))
    tabs={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in saved['tables']}
    previous={(tuple(r['counts']),r['fixed_A5_count']):r for r in saved['expanded_alternatives']}
    expanded=[];alternatives=[];checks=0;trigger_rows=0
    for cs,used,key in C.P.fixed_domain():
        large=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        candidates=[]
        for k5 in range((N-used)//5+1):
            oldrow=previous[cs,k5];best=oldrow['sixteenths_upper'];fixed=large+[('A',5)]*k5;rem=N-used-5*k5
            fixed_cost=C.value(fixed,True);options=[]
            if cs[0] and k5 and rem>=1:
                trigger_rows+=1
                # With exactly one A1 the old A2/A2 context is absent.
                # Allowing additional A1 in this table is a relaxation of
                # that class, not an assertion of equality for its witness.
                # With two or more A1 and rank10, reserve two and add the
                # complete old interaction. These classes cover every forest.
                for reserve in ([1,2] if cs[1] else [1]):
                    if rem<reserve:continue
                    for j in range((rem-reserve)//2+1):
                        rest=rem-reserve-2*j
                        value=(fixed_cost+reserve*C.item(1,key,k5)+j*C.item(2,key,k5)
                            +tabs[key,k5,True]['values'][rest]+16*j*k5
                            +(4*j*(j-1) if reserve==2 else 0))
                        checks+=1;best=max(best,value)
                        options.append(dict(reserved_A1=reserve,fixed_A2_count=j,residual_rank=rest,
                            old_context_sixteenths=4*j*(j-1) if reserve==2 else 0,
                            new_context_sixteenths=16*j*k5,sixteenths_upper=value))
            row=dict(counts=list(cs),key=list(key),fixed_A5_count=k5,previous_sixteenths_upper=oldrow['sixteenths_upper'],
                sixteenths_upper=best,new_context_options=options)
            expanded.append(row);candidates.append(row)
        alternatives.append(max(candidates,key=lambda r:(r['sixteenths_upper'],r['fixed_A5_count'])))
    maximum=max(r['sixteenths_upper'] for r in alternatives)
    witness=[('A',16),('A',16),('D',27)]
    assert maximum==6344==actual_value(witness) and sum(n for _,n in witness)==N
    assert len(alternatives)==961 and len(expanded)==2574 and checks==3405
    above=[r for r in alternatives if r['sixteenths_upper']>6264]
    return dict(exact_maximum=F(maximum,16),curvature='783/2',actual_attainer=['A16','A16','D27'],
        complete_fixed_large_vectors=961,complete_fixed_large_and_A5_vectors=2574,
        new_context_options=checks,triggered_fixed_vectors=trigger_rows,alternatives=alternatives,
        expanded_alternatives=expanded,above_curvature=above,
        inherited_tables_rechecked=len(saved['tables']),Bellman_inequalities_rechecked=old['Bellman_upper_inequalities'],
        attaining_predecessors_rechecked=old['attaining_predecessors'],residual_item_identities_rechecked=old['residual_item_checks'],
        simultaneous_type_pairs=len(types)**2,new_residual_weights=False,DP_optimization_repeated=False,
        ordinary_order13_half_included=True,formal_attainment_not_geometric_realization=True,
        context_upper_relaxation='whole A9 replaced by any rank9; whole D5 replaced by the rank5 A5 representative; witnesses not multiplied',
        scope='Exact current ambient59 maximum with this contextual coupon; no pending finite6 rule is assumed.')


def literal(docs,target):
    prev={C.R.key(r['source'],r['case']):r for r in docs[H.TAG]['current_literal_interface']['rows']};rows=[];hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        before=prev[C.R.key(stem,index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids));debit=F();roles=[];images=[]
            if sig==[len(ids),0,0]:roles=contextual_roles(g,D.finite_blocks(g,ids),marks)
            elif sig[1]==0:
                b=next((b for b in A.components(g,ids) if set(marks)<=set(b)),None)
                if b and len(b)==17:
                    labels=T.T.labels(g,b)
                    if labels is not None:images=sorted({tuple(sorted(mp[b.index(m)] for m in marks)) for mp in P.embeddings(labels,target)})
                    if set(images)&CHANGE:assert set(images)<=CHANGE;debit=F(1,2)
            inc=F(len(roles),2);pay=F(before['corners'][j]['retained_payment_lower_bound'])-debit+inc
            assert pay>=0
            if debit:hits['ideal_debits']+=1
            if inc:hits['context_credits']+=1
            corners.append(dict(marks=marks,images=images,debit=debit,increment=inc,context_roles=roles,retained_payment_lower_bound=pay))
        total=sum(r['retained_payment_lower_bound'] for r in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(rows)==216
    return dict(configurations=216,changes=dict(hits),rows=rows,minimum_total=min(r['total_lower_bound'] for r in rows))


def ideal(docs):
    current=C.ideal(docs);assert enc(current)==docs[C.TAG]['ideal']
    hrows=docs[H.TAG]['ideal_frontier_update']['rows'];by={tuple(r['ranks']):r for r in hrows};rows=[]
    for r in current['rows']:
        rs=r['ranks'];old=by[tuple(rs)];margin=r['normalized_margin']+F(rs.count(12),13)
        assert margin==F(old['normalized_margin']) and r['link_vertices']==old['link_vertices']
        saving=F(rs.count(16),17);m=margin+saving
        rows.append(dict(ranks=rs,link_vertices=r['link_vertices'],previous_normalized_margin=margin,
            normalized_saving=saving,normalized_margin=m,absolute_margin=m*r['link_vertices']))
    summary=C.U.summary(rows)
    assert summary['negative_count']==1 and summary['zero_count']==0
    assert summary['negative_rows'][0]['ranks']==(9,9,9,15,16) and summary['normalized_minimum']==-F(69,5440)
    threshold=C.U.threshold(rows,15);assert threshold['strict_total_cost_saving_threshold']==F(69,340)
    return dict(complete_profiles=1890,current_envelopes={**current['current_envelopes'],12:F(143,26),16:F(227,34)},
        current_order17_family_costs=dict(A='85',B='111',C='108',D='227/2'),rows=rows,summary=summary,
        next_order16_threshold=threshold,conditional_order16_quarter=C.U.conditional(rows,{15:F(1,4)}),
        no_new_convolution=True,saved_complete_domain_and_moments_checked=True,negative_absolute_margins_not_extended=True)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert sha256((ROOT/f'results/{H.TAG}.md').read_bytes()).hexdigest()==FROZEN_H_MD
    assert docs[H.TAG]['new_order13_maximum']=='143/2'
    proof=local(docs)
    paths={Path(__file__),Path(Z.__file__),Path(H.__file__),Path(S.__file__),Path(C.__file__),Path(T.__file__),Path(D.__file__),Path(Q.__file__),Path(P.__file__),Path(L.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_ORDER17_CENTER_CONTEXTUAL_ZERO_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        local_proof=proof,new_coupon=dict(small='whole_A2',small_marks=[0,1],large='whole_D5',large_marks=[0],large_role='unique_long_arm_endpoint',
            required_whole_unmarked_components=['A1','A9'],old_weight='0',new_weight='1/2',component_pair_increment='1',
            actual_forest_cost='number(A2)*number(D5) if number(A1)>0 and number(A9)>0, otherwise0',witnesses_not_multiplied=True),
        unchanged_D5_D12_center_weight='1/2',finite=finite(docs),ideal=ideal(docs),current_literal_interface=literal(docs,proof['target']),
        triangle_minimum=2,quadrilateral_minimum=1,completed_unrestricted_bound=59,new_dimension_bound=None,
        pending_finite6_rule_assumed=False,scope='Two arbitrary-support whole-core replacements, a conditional whole-factor half, complete unchanged-table finite bound and saved exact ideal update. No bound58.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'negative',out['ideal']['summary']['negative_count'],'literal',out['current_literal_interface']['changes'],flush=True)
