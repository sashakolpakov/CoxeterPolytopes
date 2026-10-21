#!/usr/bin/env python3
"""Delete the C13/order14 basic distance9 singleton with no new coupon.

Complete600 selected native searches, all74 unbounded protected bases,
and current216 literal interfaces. No forest/global improvement claimed.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_order14_c_distance10_terminal_projection_20261018 as C
import fresh_small_d16_coupon_and_ideal17_rebalance_20261018 as S

X=C.X;P=C.P;ROOT=C.ROOT;T=S.T;A=S.A;Q=S.Q
TAG='fresh_order14_c_basic_distance9_transfer_20261018'
CURRENT='fresh_depth6_dimension59_20261018'
DROP=(2,11)
EXPECTED=[361,378,390,525,1095,1096,1097,1108,1145,1160,1620]
PROTECTED=[367,1504]
STEMS=(CURRENT,C.TAG,S.TAG,X.SMALL,X.OVER,X.COMBINED,
       'fresh_depth6_finite7_adaptive_transfer_20261018',
       'fresh_order13_bc_sparse_terminal_transfer_20261018',
       'fresh_finite8_cycle_replacement_20261017',
       'fresh_finite8_ideal5_local_replacement_20261017',
       'fresh_affine_d13_long_pair_transfer_20261017')


def inspect(c,core,maps,target,protected=False):
    g,J,du,dv,ends=X.constraints(c,core);base=c['marked_component'];out=[]
    left=c['left'][1];right=[c['right_copy_map'][i] for i in c['right'][1]]
    white=sorted(set(range(len(g)))-set(left+right))
    white_signature=Q.inertia(A.sub(g,white))
    edge_signatures={str(b):Q.inertia(A.sub(g,sorted(set(white)|{b}))) for b in left+right}
    assert white_signature==[len(white),0,0]
    assert all(sig==[len(white)+1,0,0] for sig in edge_signatures.values())
    assert c['method'] in (('finite_leaf9_or_ideal_feature5',) if protected else ('lorentz_residual','relevant_component'))
    for mp in maps:
        image={mp[k]:node for k,node in enumerate(base)};entries=[]
        for new in sorted(set(range(14))-set(image)):
            for old,node in sorted(image.items()):
                bond=target[new][old]
                if not bond:continue
                assert bond in (3,4) and node in c['possible_contacts']
                q2=F(2 if bond==4 else 1)
                entries.append(dict(new_target_node=new,old_target_node=old,core_node=node,
                    native_label=bond,physical_contact_square=str(q2),inverse_diagonal=str(J[node][node]),
                    projection_lower_bound=str(q2*J[node][node])))
        bad=[e for e in entries if F(e['projection_lower_bound'])>2]
        assert bad and len(entries)==2
        if protected:
            assert c['forced_finite'] and c['protected_marked_leaves']
            chosen=next(e for e in bad if e['core_node'] in c['protected_marked_leaves'])
            assert chosen['native_label']==3 and chosen['projection_lower_bound']=='7/2'
            assert chosen['old_target_node'] in DROP
        else:
            chosen=next(e for e in bad if e['native_label']==4)
            assert chosen['new_target_node'] in (0,13) and F(chosen['projection_lower_bound'])>=F(5,2)
        out.append(dict(map=mp,displayed_target_nodes=sorted(image),forced_entries=entries,witness=chosen,
            outcome='genuinely_new_norm2_first_entry_projection_above2'))
    return dict(index=c['index'],method=c['method'],core_key=c['core_key'],core_gram=core['gram'],
        core_inverse=core['inverse'],base_nodes=base,preferred_corner=c['preferred_corner'],
        white_nodes=white,white_signature=white_signature,actual_edge_signatures=edge_signatures,
        end_projection_proof=ends,protected_marked_leaves=c.get('protected_marked_leaves',[]),maps=out)


def literal(previous,target):
    tags=S.B.C.R.B.C.L.STEMS
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in tags}
    norm=S.B.C.R.B.C.L.N.normalize
    before={(r['source'],norm(r['case'])):r for r in previous['combined_literal_interface']['rows']}
    rows=[];hits=[]
    for stem,index,g,white,cs in T.literal_inputs(docs):
        old=before[stem,norm(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==old['corners'][j]['marks'];ids=sorted(set(white)|set(marks));images=[]
            sig=Q.inertia(A.sub(g,ids))
            if sig[1]==0 and sig[2]>0:
                shared=next((b for b in A.components(g,ids) if set(marks)<=set(b)),None)
                if shared and len(shared)==14:
                    labels=T.T.labels(g,shared)
                    if labels is not None:
                        images=sorted({tuple(sorted(mp[shared.index(m)] for m in marks)) for mp in P.embeddings(labels,target)})
            debit=F(int(DROP in images));value=F(old['corners'][j]['retained_payment_lower_bound'])-debit
            assert value>=0
            row=dict(marks=marks,ideal_images=images,new_debit=str(debit),retained_payment_lower_bound=str(value));corners.append(row)
            if debit:hits.append(dict(source=stem,case=index,corner=j))
        total=sum(F(r['retained_payment_lower_bound']) for r in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=str(total)))
    assert len(rows)==216 and not hits
    return dict(configurations=216,changed_corners=hits,rows=rows),tags


def run():
    data={t:json.loads((ROOT/f'results/{t}.json').read_text()) for t in STEMS}
    for tag,d in data.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    assert data[CURRENT]['dimension_bound']==59
    prior=data[C.TAG];assert prior['new_cost']=='81'
    target=prior['target'];paid=set(map(tuple,prior['retained_unit_pairs']));assert len(paid)==81 and DROP in paid
    retained=paid-{DROP};autos=P.embeddings(target,target);_,dist,_=P.distance_data(target)
    assert len(retained)==80 and dist[2][11]==9
    assert all(tuple(sorted(mp[i] for i in DROP))==DROP for mp in autos)
    short8={p for p in combinations(range(14),2) if dist[p[0]][p[1]]<=8};assert short8<=retained
    tri=next(r for r in data[X.COMBINED]['triangle_pair_retention'] if(r['order'],r['kind'])==(14,'C'))
    assert set(map(tuple,tri['required_pairs']))<=retained
    small=data[X.SMALL];over=data[X.OVER]['isolation']
    rec=next(r for r in small['records'] if(r['order'],r['kind'])==(14,'C'));assert rec['target']==target
    patterns=[p for p in small['patterns'] if p['order_bound']>=14 and len(p['nodes'])<=14]
    saved={r['case']:r for r in rec['case_embeddings']};support=[];proofs=[]
    for p in patterns:
        maps=P.embeddings(p['labels'],target);pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if not maps:assert p['index'] not in saved;continue
        assert maps==saved[p['index']]['embeddings'] and pairs==sorted(map(tuple,saved[p['index']]['marked_images']))
        support.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        chosen=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks']))==DROP]
        if chosen:
            c=over['cases'][p['index']];proofs.append(inspect(c,over['cores'][c['core_key']],chosen,target))
    assert len(patterns)==600 and len(support)==43 and sum(len(r['embeddings']) for r in support)==92
    assert json.loads(json.dumps(support))==prior['reconstructed_support']
    assert [r['index'] for r in proofs]==EXPECTED and sum(len(r['maps']) for r in proofs)==22
    protected=[];protected_proofs=[]
    for c in over['cases']:
        if c['method']!='finite_leaf9_or_ideal_feature5':continue
        core=over['cores'][c['core_key']];g=A.matrix(core['gram']);base=c['marked_component'];labels=T.T.labels(g,base)
        maps=P.embeddings(labels,target) if labels is not None and len(base)<=14 else []
        chosen=[mp for mp in maps if tuple(sorted(mp[base.index(i)] for i in c['preferred_corner']))==DROP]
        protected.append(dict(index=c['index'],base_nodes=base,native_labels=labels,all_embeddings=maps,changed_maps=chosen))
        if chosen:protected_proofs.append(inspect(c,core,chosen,target,True))
    assert len(protected)==74 and [r['index'] for r in protected_proofs]==PROTECTED
    assert sum(len(r['maps']) for r in protected_proofs)==4
    routing=Counter()
    for c in over['cases']:
        if c['method']=='finite_leaf9_or_ideal_feature5':routing['protected']+=1
        elif c['method']=='nine_context_adaptive_finite_payment':routing['adaptive']+=1
        elif c['component_order_bound']<=13:routing['bounded_at_most13']+=1
        elif len(c['marked_component'])>14:routing['displayed_above14']+=1
        else:routing['eligible']+=1
    assert routing==dict(protected=74,adaptive=9,bounded_at_most13=1378,displayed_above14=553,eligible=600)
    lit,tags=literal(data[S.TAG],target)
    paths={Path(__file__),Path(C.__file__),Path(X.__file__),Path(X.L.__file__),Path(X.Q.__file__),Path(P.__file__),
           Path(S.__file__),Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(Q.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    paths.update(ROOT/f'results/{tag}.json' for tag in tags)
    return dict(status='PASS_ORDER14_C_BASIC_DISTANCE9_ZERO_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        target=target,automorphisms=autos,changed_pair=list(DROP),old_weight='1',new_weight='0',old_cost='81',new_cost='80',
        current_order14_family_costs=dict(A=70,B=82,C=80,D=82),current_order14_maximum='82',
        weighted_pairs=[dict(pair=p,weight=str(int(p in retained))) for p in combinations(range(14),2)],
        retained_unit_pairs=sorted(retained),unchanged_short8_pairs=sorted(short8),triangle_required_pairs=tri['required_pairs'],
        complete_eligible_patterns=600,complete_successful_patterns=43,complete_successful_maps=92,
        reconstructed_support=support,selected_proofs=proofs,affected_selected_maps=22,
        all_protected_bases=protected,protected_proofs=protected_proofs,affected_protected_maps=4,
        complete_routing=dict(routing),current_literal_interface=lit,
        finite_weights_unchanged=True,new_cross_component_coupons=[],triangle_minimum=2,quadrilateral_minimum=1,
        current_dimension_bound_unchanged=59,new_dimension_bound=None,
        scope='Coupon-free local C14 singleton deletion. Complete600 selected searches and all74 protected native bases;26 affected maps have arbitrary-contact first-entry obstructions. Current216 literals unchanged. B/D82 still control the all-family envelope; no forest/global theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['new_cost'],'selected',result['affected_selected_maps'],'protected',result['affected_protected_maps'],flush=True)
