#!/usr/bin/env python3
"""Independent saved-support/count/slack checks; no project imports.

Native projections, complete embeddings and old graph completeness remain
the bound producer inputs. This checker does not test any physical pair.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_new5073_raw76_support_frontier_20261018'


def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())
def key(x):return json.dumps(x,sort_keys=True,separators=(',',':'))


def run():
    d=read(TAG);assert d['status']=='PASS_NEW5073_SUPPORT_AND_RAW76_TWO_CLASS_PREFLIGHT'
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    local=read('fresh_current324_new5073_local_covers_20261018')
    current=read('fresh_current324_full_model_20261018')
    previous=read('fresh_cross_vertex_current_full_model_20261018')
    cert=next(c for c in current['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    oldcert=next(c for c in previous['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x) for f,x in cert['face_multipliers']};oldalpha={f:F(x) for f,x in oldcert['face_multipliers']}
    represented=set(local['represented_classes']);fresh=set(alpha)-set(oldalpha);outside=fresh-represented
    rows={r['face']:r for r in d['support']['rows']}
    assert len(rows)==len(d['support']['rows'])==982 and set(rows)==outside
    assert len(fresh)==1002 and len(alpha)==5073 and len(oldalpha)==5121
    counts=Counter(r['kind'] for r in rows.values())
    assert dict(counts)==d['support']['counts']['outside_kinds']
    assert sum(bool(r['saved_historical_sources']) for r in rows.values())==11
    assert len(set(alpha)&represented)==d['support']['counts']['positive_represented']==240
    forests=read('fresh_sqrt5_rank13_combined_disjoint_ports_20261017')
    faces=list(forests['faces'])
    for r in read('fresh_sqrt5_forced_base_simplex_domain_20261017')['new_faces']:
        assert r['index']==len(faces);faces.append(r)
    quads=read('fresh_sqrt5_rank13_packing_quad_separation_20261017')
    cores=quads['canonical_cores'];pr=read('fresh_sqrt5_actual_provenance_dual_20261017')['provenance_masses']
    pmap={f:i for i,f in current['prepared']['provenance_to_face']}
    for f,r in rows.items():
        assert F(r['mass'])==alpha[f] and F(r['demand_mass'])==alpha[f]*r['required']
        if f in pmap:
            p=pr[pmap[f]];ci=p['canonical_prism_core']
            assert r['kind']=='unambiguous_Q' and r['provenance']==pmap[f]
            assert (r['source_context'],r['canonical_core'],r['raw_family'],r['required'])==(
                p['source_context'],ci,cores[ci]['source_core'],cores[ci]['required'])
        else:
            ff=faces[f];assert r['required']==ff['required']
            if ff['kind']=='canonical_dotted_quad_aggregate':
                if len(ff['source_contexts'])==1:
                    ctx,ci=ff['source_contexts'][0]
                    assert r['kind']=='unambiguous_Q' and r['provenance'] is None
                    assert (r['source_context'],r['canonical_core'],r['raw_family'])==(ctx,ci,cores[ci]['source_core'])
                else:assert r['kind']=='merged_Q' and r['source_contexts']==ff['source_contexts']
            else:assert r['kind']==ff['kind']
    actual_history={f:[] for f in rows}
    for s in d['support']['historical_scan_order']:
        doc=read(s)
        for r in doc['selected_faces']:
            if r['face'] not in rows:continue
            q=rows[r['face']]
            assert (r['source_context'],r['canonical_prism_core'],r['raw_family'],r['provenance'])==(
                q['source_context'],q['canonical_core'],q['raw_family'],q['provenance'])
            actual_history[r['face']].append(s)
        del doc
    assert all(r['saved_historical_sources']==actual_history[f] for f,r in rows.items())
    ranked=sorted((r for r in rows.values() if r['kind']=='unambiguous_Q'),key=lambda r:(-F(r['demand_mass']),r['face']))
    assert [r['face'] for r in ranked]==d['support']['unambiguous_ranking']
    assert [r['face'] for r in ranked[:2]]==[43028,43066]
    domain=d['domain'];selected={r['face']:r for r in domain['selected_faces']}
    assert set(selected)=={43028,43066} and len(domain['selected_faces'])==2
    expected=set()
    for f,r in selected.items():
        assert r['required']==1 and len(r['corners'])==4 and F(r['mass'])==alpha[f]
        assert len({tuple(o['omitted']) for o in r['corners']})==4
        for o in r['corners']:expected.add((f,tuple(o['omitted']),o['vertex']))
    native=orientation=maps=occ=0;actual=set()
    for r in domain['forest_records']:
        n=len(r['literal_stars']);native+=n
        assert [s['index'] for s in r['literal_stars']]==list(range(n))
        assert len({key(s['native']['contacts']) for s in r['literal_stars']})==n
        target={int(i):{int(j):m for j,m in a.items()} for i,a in r['target'].items()}
        assert set(target)==set(range(13))
        assert all(i!=j and target[j][i]==m and m in (3,5) for i,a in target.items() for j,m in a.items())
        used=set()
        for o in r['occurrences']:
            origin=o['origin'];ids=o['literal_star_indices'];v=r['vertex'];f=origin['face']
            assert origin['vertex']==v and origin['raw_family']==76 and origin['provenance'] is None
            assert len(ids)==len(set(ids))==origin['literal_orientations'] and all(0<=i<n for i in ids)
            assert not used.intersection(ids);used.update(ids)
            k=(f,tuple(origin['omitted']),v);assert k not in actual;actual.add(k)
            assert origin in selected[f]['corners']
            occ+=1;orientation+=len(ids);maps+=origin['active_embeddings']
        assert used==set(range(n))
    assert actual==expected and domain['counts']==dict(faces=2,actual_corner_occurrences=occ,
        vertex_types=len(domain['forest_records']),active_embeddings=maps,
        literal_orientation_occurrences=orientation,distinct_native_stars=native)
    mu={v:F(x) for v,x in local['vertex_multipliers']};saved={r['vertex']:r for r in local['forest_records']}
    inv=d['inventory'];irs={r['vertex']:r for r in inv['forests']};assert len(irs)==len(inv['forests'])==6
    total=increase=F();pairs=0
    for r in domain['forest_records']:
        v=r['vertex'];z=irs[v];old=saved.get(v);cover=F(old['exact_local_optimum']) if old else F()
        extra=sum((alpha[o['origin']['face']] for o in r['occurrences']),F());oldmu=mu.get(v,F())
        upper=max(oldmu,cover+extra);gain=upper-oldmu
        assert F(z['old_local_cover'])==cover<=oldmu==F(z['old_vertex_mass'])
        assert F(z['extra_singleton_mass'])==extra and F(z['maximum_mu_increase'])==gain
        assert F(z['hypothetical_vertex_mass_upper'])==upper and z['candidate_native_stars']==len(r['literal_stars'])
        n=z['candidate_native_stars'];assert z['new_pair_positions_upper']==z['old_literal_stars']*n+n*(n-1)//2
        total+=extra;increase+=gain;pairs+=z['new_pair_positions_upper']
    assert total==4*sum((alpha[f] for f in selected),F())==F(inv['total_extra_singleton_mass'])
    assert increase==F(inv['maximum_mu_increase']) and pairs==inv['candidate_new_pair_positions_upper']==7304
    num=F(local['objective_numerator']);den=sum(mu.values(),F())
    assert num==F(inv['numerator']) and den==F(inv['baseline_normalizer'])
    assert num/den==F(inv['baseline_ratio']) and num/(den+increase)==F(inv['singleton_cover_ratio_lower_bound'])>12
    assert num/12-den==F(inv['required_increment_to_reach12'])>increase
    assert all(d[k]==0 for k in ('pair_tests','clique_enumerations','local_LPs','global_optimizer_invocations'))
    print('PASS_SAVED_NEW5073_RAW76_SUPPORT_ACCOUNTING',len(rows),'support rows',occ,'corners',native,
          'saved stars',len(irs),'slack records','ratio_lower',float(num/(den+increase)),flush=True)


if __name__=='__main__':run()
