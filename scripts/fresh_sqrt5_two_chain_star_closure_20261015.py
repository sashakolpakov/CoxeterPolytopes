#!/usr/bin/env python3
"""Containment and complete first-white domains of the112 two-chain stars.

The all-length producer is immutable. This file inherits its finite output,
literal marked canonicalization, and already reviewed contact/closure helpers.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_singular_attachment_closure_20261015 as ATT
import fresh_sqrt5_sealed_weights_20261015_geometry as GEO

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_two_chain_star_closure_20261015'
F=ATT.F;Q=ATT.Q;C=ATT.C


def catalogue(first,base,old,attachment):
    known=ATT.known_catalogue(first,base,old)
    for rec in known.values():rec['sealed']=True
    for r in attachment['records']:
        if r['status']!='NEW_SEALED_PAID':continue
        key,order=C.canonical(r['row']);assert key not in known
        known[key]=dict(family='attachment',index=r['index'],row=r['row'],order=order,sealed=True)
    assert len(known)==169
    for r in first['records']:
        if not r['accepted']:continue
        key,order=C.canonical(r['row']);assert key not in known
        known[key]=dict(family='closed_start',index=r['index'],row=r['row'],order=order,sealed=False)
    assert len(known)==174
    return known


def contained(row,known):
    g=Q.edges(row);checked=0
    sizes=sorted({k[0]-4 for k in known})
    for size in sizes:
        for white in combinations(row['white'],size):
            ids=list(range(4))+list(white)
            sub=dict(nodes=len(ids),white=list(range(4,len(ids))),
                     labels=[g[Q.pair(i,j)] for i,j in combinations(ids,2)])
            key,order=C.canonical(sub);checked+=1
            if key not in known:continue
            target=known[key];mp=[None]*len(ids)
            for a,b in zip(order,target['order']):mp[a]=b
            C.check_map(sub,target['row'],mp)
            outside=sorted(set(row['white'])-set(white))
            if target['sealed']:
                assert all(g[Q.pair(i,j)]==2 for i in outside for j in ids)
            return dict(family=target['family'],index=target['index'],sealed=target['sealed'],
                        nodes=ids,to_target=mp,outside_white=outside),checked
    return None,checked


def run():
    load=lambda tag:json.loads((ROOT/f'results/{tag}.json').read_text())
    source=load('fresh_sqrt5_two_chain_completion_20261015')
    first=load(F.TAG);base=load('fresh_sqrt5_quad_exhaustion_20261015')['labelled_union_supports']
    old=load('fresh_sqrt5_quad_contacts_20261015');attachment=load(ATT.TAG)
    original=load('fresh_sqrt5_forest_sensitive_weights_20261015')
    sealed=load('fresh_sqrt5_sealed_weights_20261015');extended=load('fresh_sqrt5_extended_forest_weights_20261015')
    known=catalogue(first,base,old,attachment)
    classes=[];keys={};literal=0
    for ri,r in enumerate(source['records']):
        for si,star in enumerate(r.get('completed_stars',[])):
            assert star['status']=='UNPAID_PRINCIPAL_STAR'
            row=star['row'];key,order=C.canonical(row)
            if key not in keys:
                keys[key]=len(classes);classes.append(dict(index=len(classes),row=row,order=order,origins=[]))
            rec=classes[keys[key]];mp=[None]*row['nodes']
            for a,b in zip(order,rec['order']):mp[a]=b
            C.check_map(row,rec['row'],mp)
            rec['origins'].append(dict(source_record=ri,completed_star_index=si,
                source_class_index=r['source_class_index'],trial_indices=r['trial_indices'],
                lengths=star['lengths'],to_representative=mp))
            literal+=1
    assert literal==112 and len(classes)==70
    counts=Counter();trialcounts=Counter();subsets=0;records=[]
    for r in classes:
        row=r['row'];match,checked=contained(row,known);subsets+=checked
        r['known_subcore']=match;r['white_subsets_checked']=checked
        r['payment']=ATT.payment(row,original,sealed,extended)
        if match:
            r['status']='KNOWN_SEALED_CONTAINMENT' if match['sealed'] else 'KNOWN_CLOSED_START_CONTAINMENT'
            assert r['payment']['paid']
        else:
            domain=ATT.next_domain(row);r.update(domain)
            trialcounts.update(t['status'] for t in r['trials'])
            r['status']=('OPEN_WITH_CONTACTS' if r['accepted'] else
                         'NEW_SEALED_PAID' if r['payment']['paid'] else 'NEW_SEALED_UNPAID')
        counts[r['status']]+=1;records.append(r)
        print('class',r['index'],r['status'],'profiles',len(r.get('trials',[])),
              'accepted',len(r.get('accepted',[])),'payment',r['payment']['sums_without_with_external_a1'],flush=True)
    # The two nonsealed classes have only one permitted unpaid contact each.
    # Close those literal children, retaining the original first-stage counts.
    child_closures=[]
    for r in records:
        if r['status']!='OPEN_WITH_CONTACTS':continue
        assert len(r['accepted'])==1
        contact=r['accepted'][0];row=F.extension(r['row'],contact)
        domain=ATT.next_domain(row);assert not domain['accepted']
        pay=ATT.payment(row,original,sealed,extended);assert pay['paid']
        child_closures.append(dict(parent=r['index'],contacts=contact,row=row,
                                   **domain,payment=pay,status='SEALED_PAID_CHILD'))
    geometry=[];obligations=[]
    for r in records:
        if r['status']=='NEW_SEALED_UNPAID':
            geo=GEO.geometry(r['row'])
            geometry.append(dict(index=r['index'],**geo))
            connected=[c for c in r['payment']['corners'] if c['marked_component'] is not None]
            assert len(connected)==1
            c=connected[0];a,b=c['marked'];g=Q.edges(r['row'])
            H=ATT.T.graph(r['row']['white']+[a,b],g);ids=ATT.T.component(H,a)
            J={i:{j:m for j,m in H[i].items() if j in ids} for i in ids}
            branch=next(i for i in J if len(J[i])==3)
            obligations.append(dict(index=r['index'],component=c['marked_component'],
                marked=[a,b],pair_orbit=c['pair_orbit'],
                distances_from_trivalent_node=sorted([Q.distance(J,branch,a),Q.distance(J,branch,b)]),
                core_has_a1=c['core_has_a1'],minimum_required_weight='1',
                required_A1_contexts=[True] if c['core_has_a1'] else [False,True],
                scope='Component-local weights with zero weight on disconnected marked pairs; the other three corners are disconnected.'))
    assert counts==Counter(KNOWN_SEALED_CONTAINMENT=32,KNOWN_CLOSED_START_CONTAINMENT=7,
                           NEW_SEALED_PAID=26,NEW_SEALED_UNPAID=3,OPEN_WITH_CONTACTS=2)
    assert sum(trialcounts.values())==303 and trialcounts['ACCEPT_UNPAID']==2
    assert len(child_closures)==2 and sum(len(r['trials']) for r in child_closures)==10
    assert [r['index'] for r in geometry]==[50,55,62]
    paths=[f'scripts/{TAG}.py',f'scripts/{ATT.TAG}.py',f'scripts/{F.TAG}.py',
           'scripts/fresh_sqrt5_all_quad_unions_20261015_classes.py',
           'scripts/fresh_sqrt5_sealed_weights_20261015_geometry.py',
           'results/fresh_sqrt5_sealed_weights_20261015_geometry.json',
           'scripts/fresh_sqrt5_quad_filters_20261015_review.py',
           'results/fresh_sqrt5_quad_filters_20261015_review.json',
           'results/fresh_sqrt5_sealed_weights_20261015_review.json']
    for tag in('fresh_sqrt5_two_chain_completion_20261015',F.TAG,ATT.TAG,ATT.TAG+'_review',
               'fresh_sqrt5_quad_exhaustion_20261015','fresh_sqrt5_quad_contacts_20261015',
               'fresh_sqrt5_forest_sensitive_weights_20261015','fresh_sqrt5_extended_forest_weights_20261015',
               'fresh_sqrt5_sealed_weights_20261015'):
        paths.append(f'results/{tag}.json')
    Q.manifest_check(paths)
    return dict(status='PASS_TWO_CHAIN_STAR_CONTAINMENT_AND_ONE_WHITE_DOMAINS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        input_literal_stars=literal,marked_classes=len(classes),known_sealed_cores=169,known_closed_starts=5,
        white_subsets_checked=subsets,counts=dict(counts),trial_counts=dict(trialcounts),records=records,
        child_closures=child_closures,unpaid_geometry=geometry,new_weight_obligations=obligations,
        open_attachment_states=0,
        new_sealed_unpaid_marked_classes=[50,55,62],
        scope='All112 saved finite two-chain principal stars are quotiented by literal marked diagram bijections. Every induced marked white subset is searched until a known sealed or fully closed start is found; an unmatched result checks all sizes. Unmatched classes receive complete finite-alphabet next-white domains, kernel/inertia tests and both externalA1 contexts of the frozen forest-sensitive rule. The only two surviving contact profiles lead to sealed paid children, so the explicit attachment queue is exhausted. Sealing is conditional on no persistent short-path/H4 payment; it fixes marked components for all unpaid external attachments. Known containment inherits the separately reviewed169 sealed cores or5 closed-start proof. Three new unpaid sealed classes are realized as exact compact quadrilateral faces of local arithmetic-compatible chambers. They remain new weight obligations, not exclusions. Other initial attachment topologies and a new global dimension bound remain outside scope.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],out['trial_counts'])
