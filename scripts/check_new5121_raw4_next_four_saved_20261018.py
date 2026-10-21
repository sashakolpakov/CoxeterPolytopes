#!/usr/bin/env python3
"""Saved-domain accounting and singleton-cover bound; no project imports.

Native projections, complete embeddings and historical graph reconstruction
remain producer inputs. No physical pair/clique/LP work is repeated.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_new5121_raw4_next_four_domain_20261018'
BASE='fresh_new5121_raw4_two_face_local_covers_20261018'
OLD='fresh_cross_vertex_historical32_local_covers_20261018'
MIG='fresh_cross_vertex_support_migration_20261018'

def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def okey(o):return o['face'],o['provenance'],tuple(o['omitted']),o['vertex'],o['raw_family']

def run():
    d,b,u,m=map(read,(TAG,BASE,OLD,MIG));iv=d['cover_inventory']
    assert d['status']=='PASS_NEW5121_RAW4_NEXT_FOUR_EXACT_DOMAIN_AND_BENEFIT_CEILING'
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    rep=set(m['current290_classes'])|{r['selected_face']['face'] for r in u['added_faces']}|{r['face'] for r in b['selected_faces']}
    assert len(rep)==d['baseline_represented_classes']==324
    missing=[r for r in m['quad_rows'] if r['face'] not in rep]
    rank=sorted((r for r in missing if r.get('raw_family')==4),key=lambda r:(-F(r['demand_mass']),r['face']))
    assert len(missing)==d['missing_positive_Q_classes']==1954 and len(rank)==d['missing_raw4_classes']==89
    assert Counter(r['provenance_kind'] for r in missing)==d['missing_Q_provenance_counts']
    assert F(d['missing_raw4_demand_mass'])==sum((F(r['demand_mass']) for r in rank),F())
    assert [r['face'] for r in rank[:4]]==[r['face'] for r in d['selected_faces']]==[14935,14368,14684,14670]
    alpha={i:F(x) for i,x in b['face_multipliers']};mu={v:F(x) for v,x in b['vertex_multipliers']}
    sourceorigins={}
    for face,rankrow in zip(d['selected_faces'],rank[:4]):
        assert face['source_context']==rankrow['source_context'] and face['canonical_prism_core']==rankrow['canonical_core']==4
        assert face['required']==3 and len(face['corners'])==6
        assert alpha[face['face']]==F(face['mass']) and 3*F(face['mass'])==F(rankrow['demand_mass'])==F(face['demand_mass'])
        assert len({tuple(c['omitted']) for c in face['corners']})==6
        for o in face['corners']:
            key=okey(o);assert key not in sourceorigins;sourceorigins[key]=o
    rows={r['vertex']:r for r in d['forest_records']}
    inventory={r['vertex']:r for r in iv['forests']}
    assert len(rows)==len(inventory)==12 and set(rows)==set(inventory)
    origins=set();maps=orient=stars=newstars=neworbits=pairs=0;total=inc=F();oldcount=newcount=0
    for v,raw in rows.items():
        r=inventory[v];assert raw['forest']==r['forest'];oldn=r['old_literal_stars']
        native=raw['literal_stars'];mapping=r['raw4_to_merged_star_indices']
        assert len(native)==len(mapping)==len(set(mapping))
        assert [x['index'] for x in native]==list(range(len(native)))
        assert len({json.dumps(x['native']['contacts'],sort_keys=True) for x in native})==len(native)
        # In this saved domain every native star is new; there is no old-orbit overlap.
        assert set(mapping)==set(range(oldn,oldn+r['new_literal_stars']))
        assert len(r['actual_occurrences'])==len(raw['occurrences'])
        groups={};extra=Counter()
        for a,o in zip(r['actual_occurrences'],raw['occurrences']):
            origin=o['origin'];key=okey(origin)
            assert origin==sourceorigins[key]==a['origin'] and key not in origins;origins.add(key)
            ids=sorted(mapping[i] for i in o['literal_star_indices'])
            assert ids==a['literal_star_indices'] and len(ids)==len(set(ids))==origin['literal_orientations']
            gi=a['orbit']
            if gi in groups:assert groups[gi]==ids
            else:assert not set(ids)&{i for old in groups.values() for i in old};groups[gi]=ids
            extra[gi]+=alpha[origin['face']];maps+=origin['active_embeddings'];orient+=origin['literal_orientations']
        assert set(groups)==set(range(r['old_corner_orbits'],r['old_corner_orbits']+r['new_corner_orbits']))
        assert {i for ids in groups.values() for i in ids}==set(mapping)
        assert [[i,str(x)] for i,x in sorted(extra.items())]==r['extra_orbit_demands']
        value=F(r['baseline_local_cover']);prior=r['baseline']
        if prior is None:assert oldn==r['old_corner_orbits']==value==0;newcount+=1
        else:
            oldcount+=1;record=next(x for x in read(prior['source'])['forest_records'] if x['vertex']==v)
            assert jhash(record)==prior['record_sha256'] and value==F(record['exact_local_optimum'])==F(prior['exact_local_optimum'])
        assert value<=mu.get(v,F())==F(r['baseline_vertex_mass'])
        add=sum(extra.values(),F());upper=max(mu.get(v,F()),value+add)
        assert add==F(r['singleton_extension_mass']) and upper==F(r['possible_vertex_mass_upper'])
        delta=upper-mu.get(v,F());assert delta==F(r['maximum_possible_mu_increase'])
        n=r['new_literal_stars'];np=oldn*n+n*(n-1)//2
        assert np==r['candidate_new_pair_positions']
        total+=add;inc+=delta;stars+=len(native);newstars+=n;neworbits+=r['new_corner_orbits'];pairs+=np
    assert origins==set(sourceorigins) and len(origins)==24
    assert d['counts']==dict(faces=4,actual_corner_occurrences=24,vertex_types=12,active_embeddings=maps,literal_orientation_occurrences=orient,distinct_native_stars=stars)
    assert (oldcount,newcount,newstars,neworbits,pairs)==(iv['affected_existing_forests'],iv['potential_new_forests'],iv['new_literal_stars'],iv['new_corner_orbits'],iv['candidate_new_pair_positions'])==(11,1,102,16,18147)
    assert total==6*sum((alpha[r['face']] for r in d['selected_faces']),F())==F(iv['singleton_extension_mass'])
    assert inc==F(iv['maximum_possible_mu_increase'])
    num,den=F(b['objective_numerator']),F(b['normalizer']);assert den==sum(mu.values(),F())
    assert num==F(iv['current_objective_numerator']) and den==F(iv['current_normalizer'])
    assert num/den==F(iv['current_ratio']) and num/(den+total)==F(iv['coarse_ratio_lower_bound'])
    assert num/(den+inc)==F(iv['slack_aware_ratio_lower_bound'])>12
    assert not d['current_compatibility_graphs_changed'] and d['pair_compatibility_tests']==0 and not d['new_global_optimization']
    print('PASS_RAW4_NEXT_FOUR_SAVED_ACCOUNTING',dict(origins=len(origins),stars=stars,orbits=neworbits,forests=len(rows),potential_pairs=pairs,ratio_lower_bound=str(num/(den+inc))))

if __name__=='__main__':run()
