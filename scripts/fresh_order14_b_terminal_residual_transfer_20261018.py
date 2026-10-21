#!/usr/bin/env python3
"""Delete the last two nonbasic B13 pairs by first-entry/two-root obstructions."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

ROOT=X.ROOT;P=X.P
TAG='fresh_order14_b_terminal_residual_transfer_20261018'
CURRENT='fresh_sparse_coupon_order11_13_dimension62_20261018'
PRIOR='fresh_order14_affine_b_fork_projection_transfer_20261017'
DROP={(9,12),(9,13)}
TAIL=[630,638,1648,1668,1849,1935,2242,2328]
STEMS=(CURRENT,PRIOR,X.SMALL,X.OVER,X.COMBINED,
       'fresh_order11_order17_dimension63_20261018',
       'fresh_order13_d_center_pair_zero_transfer_20261018',
       'fresh_order16_affine_b_singleton_projection_transfer_20261018',
       'fresh_order17_e6_unit_coupon_zero_transfer_20261018',
       'fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018')


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for tag,d in data.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    small=data[X.SMALL];rec=next(r for r in small['records'] if (r['order'],r['kind'])==(14,'B'))
    target=rec['target'];prior=data[PRIOR];assert target==prior['target']
    paid=set(map(tuple,prior['retained_pairs']));assert len(paid)==84 and DROP<=paid
    retained=paid-DROP;basic=X.L.S.basic_pairs(target);assert retained==basic and len(basic)==82
    autos=P.embeddings(target,target);_,ds,_=P.distance_data(target)
    assert all({tuple(sorted(a[i] for i in p)) for p in DROP}==DROP for a in autos)
    assert all(ds[i][j]==10 for i,j in DROP)
    short8={p for p in combinations(range(14),2) if ds[p[0]][p[1]]<=8};assert short8<=retained
    triangle=next(r for r in data[X.COMBINED]['triangle_pair_retention'] if (r['order'],r['kind'])==(14,'B'))
    assert set(map(tuple,triangle['required_pairs']))<=retained
    patterns=[p for p in small['patterns'] if p['order_bound']>=14 and len(p['nodes'])<=14]
    saved={r['case']:r for r in rec['case_embeddings']};support=[];proofs=[];outcomes=Counter();bounds=Counter()
    over=data[X.OVER]
    for p in patterns:
        maps=P.embeddings(p['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if not maps:assert p['index'] not in saved;continue
        assert maps==saved[p['index']]['embeddings'] and pairs==sorted(map(tuple,saved[p['index']]['marked_images']))
        support.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        selected=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks'])) in DROP]
        if not selected:continue
        c=over['isolation']['cases'][p['index']];core=over['isolation']['cores'][c['core_key']]
        assert c['method'] in ('lorentz_residual','relevant_component')
        g,J,du,dv,ends=X.constraints(c,core);base=c['marked_component'];rows=[]
        for mp in selected:
            image={mp[k]:i for k,i in enumerate(base)};missing=sorted(set(range(14))-set(image));entries=[]
            for new in missing:
                for old,i in sorted(image.items()):
                    bond=target[new][old]
                    if bond:
                        assert bond in (3,4) and i in c['possible_contacts']
                        q2=F(2 if bond==4 else 1)
                        entries.append(dict(new_target_node=new,old_target_node=old,core_node=i,
                            native_label=bond,physical_contact_square=str(q2),inverse_diagonal=str(J[i][i]),
                            projection_lower_bound=str(q2*J[i][i])))
            bad=next((r for r in entries if F(r['projection_lower_bound'])>2),None)
            if bad:
                outcome='first_entry_projection_above_norm';bounds[bad['projection_lower_bound']]+=1
                row=dict(map=mp,forced_entries=entries,witness=bad,outcome=outcome)
            else:
                assert c['index'] in TAIL and missing==[10,11]
                assert len(entries)==1 and entries[0]['new_target_node']==10
                assert entries[0]['old_target_node']==9 and entries[0]['native_label']==3
                i=image[9];assert c['relevant_nodes']==[i] and 1<J[i][i]<2
                rho=2-J[i][i];assert rho in (F(1,2),F(3,4))
                two=[]
                for end in c['end_sets']:
                    for a,b in combinations(end,2):
                        lower=J[a][a]+J[b][b]+2*J[a][b];assert lower>2
                        two.append(dict(nodes=[a,b],unit_pair_projection_lower_bound=str(lower)))
                assert 2*J[i][i]>2 # e's native contact is necessarily label3.
                contacts=[]
                for j in c['possible_contacts']:
                    value=J[i][j]
                    shortcut=min(du[i]+dv[j],du[j]+dv[i])+3
                    if value<0:assert shortcut<=8
                    contacts.append(dict(core_node=j,projection_coefficient=str(value),
                        disposition='opposite_black_short_path' if value<0 else 'nonnegative_projection',
                        path_length=shortcut if value<0 else None))
                # Once every negative coefficient gives a paid path, arbitrary
                # remaining native magnitudes give projection(e,f)>=0. The
                # norm2 root f satisfies ||f_perp||^2<=2 in the positive complement.
                assert 2*rho<2
                outcome='doubled_terminal_residual_cauchy_contradiction'
                row=dict(map=mp,forced_entries=entries,outcome=outcome,
                    first_new_target_node=10,doubled_terminal_target_node=11,old_entry_node=i,
                    first_entry_residual_square=str(rho),same_end_two_contact_exclusions=two,
                    terminal_complete_contact_signs=contacts,
                    terminal_physical_norm='2',native_terminal_bond=4,
                    required_residual_inner_product_square_at_least='2',
                    Cauchy_upper_bound=str(2*rho),strict_gap=str(2-2*rho),
                    terminal_is_new_reason='The first root has exactly one old-core contact; an old terminal in a different corner block would be a second old-core contact. Omitted blacks are absent from the actual corner.')
            outcomes[outcome]+=1;rows.append(row)
        proofs.append(dict(index=c['index'],method=c['method'],core_key=c['core_key'],
            core_gram=core['gram'],core_inverse=core['inverse'],base_nodes=base,
            preferred_corner=c['preferred_corner'],end_projection_proof=ends,maps=rows))
    assert len(patterns)==600 and len(support)==len(saved)==207
    assert sum(len(r['embeddings']) for r in support)==506 and len(proofs)==53
    assert outcomes==dict(first_entry_projection_above_norm=92,doubled_terminal_residual_cauchy_contradiction=16)
    assert sorted(r['index'] for r in proofs if any(m['outcome']=='doubled_terminal_residual_cauchy_contradiction' for m in r['maps']))==TAIL
    costs=dict(prior['new_order14_costs']);assert costs==dict(A=70,B=84,C=84,D=82);costs['B']=82
    paths=[Path(__file__),Path(X.__file__),Path(X.L.__file__),Path(X.Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER14_B_TWO_TERMINAL_RESIDUAL_DELETIONS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        affine_type='B13',order=14,target=target,automorphisms=autos,deleted_pairs=sorted(DROP),
        old_cost=84,new_cost=82,retained_unit_pairs=sorted(retained),basic_pairs_retained=sorted(basic),
        short8_pairs_retained=sorted(short8),triangle_pairs_retained=triangle['required_pairs'],
        complete_eligible_patterns=600,successful_patterns=207,complete_successful_maps=506,
        reconstructed_support=support,affected_contexts=53,affected_maps=108,cases=proofs,
        outcomes=dict(outcomes),first_entry_projection_histogram=dict(bounds),
        new_order14_family_costs=costs,new_order14_maximum_without_separate_C_change=84,
        finite_weights_unchanged=True,added_finite_coupons=[],triangle_minimum=2,quadrilateral_minimum=1,
        new_dimension_bound=None,
        scope='Complete600-pattern/506-map support reconstruction and108 affected maps.92 forced first entries exceed norm2;16 require a doubled terminal whose arbitrary-contact residual contradicts positive-complement Cauchy. No completed zeros, bounded spectator census, finite coupon, or global theorem is used.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['affected_maps'],out['outcomes'],flush=True)
