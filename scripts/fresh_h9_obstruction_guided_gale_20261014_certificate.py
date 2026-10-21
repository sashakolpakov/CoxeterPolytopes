#!/usr/bin/env python3
"""Final bounded ledger plus separate two-cap Gale realization replay."""
import argparse,json
from collections import Counter
from fresh_h9_obstruction_guided_gale_20261014_independent import one,load,sha,ROOT

TAG='fresh_h9_obstruction_guided_gale_20261014'


def run():
    first=load(f'results/{TAG}.json');second=load(f'results/{TAG}_circuits.json')
    audit=load(f'results/{TAG}_independent.json')
    flip=load(f'results/{TAG}_flip_diagnostic.json')
    old=load('results/fresh_h9_rank5_gale_20261014.json')
    assert first['status']=='FINISHED_BOUNDED_OBSTRUCTION_GUIDED_GALE_WALK'
    assert second['status']=='FINISHED_BOUNDED_BALANCED_POSITIVE_CIRCUIT_PROPOSALS'
    assert audit['status']=='PASS_ALL_SAVED_GUIDED_GALE_CASES'
    assert first['new_saved_nerves']==audit['cases_checked']==72
    assert second['new_saved_nerves']==2 and first['counts']['proposals']==12000
    assert second['completed_boundary_attempts']==2000 and second['counts']['perturbed_proposals']==1296
    assert first['first_full_category_survivor'] is second['first_full_category_survivor'] is None
    for suffix,doc in [('',first),('_circuits',second)]:
        state=load(f'results/{TAG}{suffix}_state.json')
        assert state['status']==doc['status'] and state['counts']==doc['counts'] and state['survivor'] is None
    known=set()
    for ref in old['cases']:
        assert sha(ref['path'])==ref['sha256']
        known.add(tuple(load(ref['path'])['geometry']['canonical_minimal_nonfaces']))
    assert len(known)==1063
    records={r['case']:r for r in audit['records']};counts=Counter();new=[];vertices=0
    for ref in first['cases']:
        assert sha(ref['path'])==ref['sha256']
        d=load(ref['path']);key=tuple(d['geometry']['canonical_minimal_nonfaces'])
        assert key not in known;known.add(key)
        r=records[ref['case']]
        assert r['input']==ref['path'] and r['input_sha256']==ref['sha256']
        assert load(f'results/{TAG}_audit_cases/case_{ref["case"]:05d}.json')['result']==r
        counts[d['gate']['gate']]+=1;vertices+=r['vertices'];new.append(key)
    extra=[]
    for ref in second['cases']:
        d=load(ref['path']);key=tuple(d['geometry']['canonical_minimal_nonfaces'])
        assert key not in known;known.add(key);new.append(key)
        r=one(ref['path'],ref['sha256']);extra.append(r);vertices+=r['vertices']
        counts[d['gate']['gate']]+=1
        origin=d['origin'];parent=load(origin['seed']);V=parent['geometry']['gale']
        B=origin['positive_circuit'];i,j=origin['changed_rows'];n=origin['hyperplane']
        assert len(B)==len(set(B))==5 and i in B and j not in B
        shift=[sum(V[k][c] for k in B) for c in range(5)]
        W=[row[:] for row in V]
        W[i]=[a-b for a,b in zip(W[i],shift)];W[j]=[a+b for a,b in zip(W[j],shift)]
        assert all(sum(W[k][c] for k in B)==0 for c in range(5))
        assert all(sum(a*b for a,b in zip(W[k],n))==0 for k in B)
        positive=[k for k in range(15) if k not in B and sum(a*b for a,b in zip(W[k],n))>0]
        negative=[k for k in range(15) if k not in B and sum(a*b for a,b in zip(W[k],n))<0]
        assert len(positive)==len(negative)==5
        assert positive==origin['positive_complement'] and negative==origin['negative_complement']
        target={k for k in range(15) if origin['target_missing_face']>>k&1}
        assert target<=set(positive) or target<=set(negative)
        assert origin['target_missing_face'] in parent['geometry']['minimal_nonfaces']
        col,side=origin['perturbation_column'],origin['perturbation_side'];assert side in (-1,1)
        W[i][col]+=side;W[j][col]-=side
        assert W==d['geometry']['gale']
    assert len(new)==74 and counts=={'literal_induced_pattern':70,'two_lanner_blocks':4}
    assert min(ref['recognized_support_count'] for ref in first['cases'])==4
    assert flip['counts']==dict(exact_bistellar_recipes=625,large_minimal_nonface=543,
                               passes_incidence_filters=82,exact_literal_obstruction=82)
    paths=[f'scripts/{TAG}_certificate.py',f'scripts/{TAG}_independent.py',
           'scripts/fresh_h9_rank5_gale_20261014_independent.py',f'results/{TAG}.json',
           f'results/{TAG}_circuits.json',f'results/{TAG}_independent.json',
           f'results/{TAG}_state.json',f'results/{TAG}_circuits_state.json',
           f'results/{TAG}_flip_diagnostic.json']
    return dict(status='PASS_BOUNDED_GUIDED_PACKAGE_AND_SEPARATE_TWO_CIRCUIT_REPLAYS',
                sources={p:sha(p) for p in paths},new_realized_nerves=74,prior_keys_avoided=1063,
                exact_exclusions=dict(counts),unresolved_saved_nerves=0,full_label_searches=0,metric_solves=0,
                separate_incidence_counts=dict(five_minors=74*3003,six_circuits=74*5005,vertex_links=74*15,
                                               actual_vertices=vertices,guidance_witnesses=audit['guidance_witnesses_rechecked']+
                                               sum(r['all_guidance_witnesses_rechecked'] for r in extra)),
                circuit_case_independent_replays=extra,best_recognized_support_count=4,
                exact_combinatorial_diagnostic=flip['counts'],
                scope='The72 guided-walk cases inherit their completed separate Bareiss/face audit; this command freshly replays both circuit cases and their literal boundary recipes. Inherited local obstruction arithmetic is not repeated. The625 combinatorial flip recipes are a separate exact producer/replay diagnostic, with no realization or independent full-census claim. All searches are bounded.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    path=ROOT/f'results/{TAG}_certificate.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_realized_nerves'],out['exact_exclusions'])
