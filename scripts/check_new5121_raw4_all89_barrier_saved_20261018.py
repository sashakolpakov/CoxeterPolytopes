#!/usr/bin/env python3
"""Independent saved actual-corner/cardinality accounting; standard library only.

Native template projections and frozen old cover feasibility are theorem inputs.
No wall tests, clique enumeration, covering LP or global optimization is run.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_new5121_raw4_all89_cover_barrier_20261018'
BASE='fresh_new5121_raw4_two_face_local_covers_20261018'
U='fresh_cross_vertex_historical32_local_covers_20261018'
L='fresh_cross_vertex_current290_local_covers_20261018'
M='fresh_cross_vertex_support_migration_20261018'
MODEL='fresh_cross_vertex_current_full_model_20261018'
ACTUAL='fresh_sqrt5_actual_provenance_dual_20261017'
FORESTS='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'

@lru_cache(None)
def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def h(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def run():
    d,b,u,l,m=map(read,(TAG,BASE,U,L,M))
    assert d['status']=='PASS_ALL89_RAW4_CORNER_MASS_AND_PAIR_ONLY_CARDINALITY_BARRIER'
    for p,x in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==x,p
    alpha={f:F(x) for f,x in b['face_multipliers']};mu={v:F(x) for v,x in b['vertex_multipliers']}
    assert b['face_multipliers']==d['face_multipliers'] and b['vertex_multipliers']==d['vertex_multipliers']
    latest={r['vertex']:(L,r) for r in l['forest_records']}
    latest.update({r['vertex']:(U,r) for r in u['forest_records']});latest.update({r['vertex']:(BASE,r) for r in b['forest_records']})
    assert len(latest)==600
    groups={}
    for r in l['forest_records']:
        source=read(r['group_source']);records=source.get('forest_records',[source])
        saved=next(x for x in records if x['vertex']==r['vertex'])
        groups[r['vertex']]=saved.get('corner_orbits') or saved['union']['corner_orbits']
    groups.update({r['vertex']:r['domain']['corner_orbits'] for r in u['forest_records']})
    groups.update({r['vertex']:r['domain']['corner_orbits'] for r in b['forest_records']})
    oldn={v:sum(len(g['literal_stars']) for g in gs) for v,gs in groups.items()}
    assert set(oldn)==set(latest)
    for v,gs in groups.items():assert oldn[v]==len({s for g in gs for s in g['literal_stars']})
    t={v:F(r['exact_local_optimum']) for v,(tag,r) in latest.items()}
    assert all(x<=mu.get(v,F()) for v,x in t.items())
    rep=set(m['current290_classes'])|{r['selected_face']['face'] for r in u['added_faces']}|{r['face'] for r in b['selected_faces']}
    expected=sorted((r for r in m['quad_rows'] if r['face'] not in rep and r.get('raw_family')==4),key=lambda r:(-F(r['demand_mass']),r['face']))
    assert len(expected)==89 and [r['face'] for r in d['faces']]==[r['face'] for r in expected]
    pmap={f:j for j,f in read(MODEL)['prepared']['provenance_to_face']};pr=read(ACTUAL)['provenance_masses']
    fdata=read(FORESTS);vertex={tuple(r['forest']):r['index'] for r in fdata['vertices']}
    templates={tuple(r['omitted']):r for r in d['templates']};assert len(templates)==6
    rows={};embedding_rows={};byv=defaultdict(set)
    for row,q in zip(d['faces'],expected):
        f=row['face'];assert row['provenance_kind']==q['provenance_kind'] and row['provenance']==q['provenance']
        assert row['source_context']==q['source_context'] and row['required']==q['required']==3
        assert fdata['faces'][f]['source_contexts']==[[row['source_context'],4]]
        if q['provenance_kind']=='explicit':
            j=q['provenance'];assert pmap[f]==j and pr[j]['source_context']==q['source_context']
            assert pr[j]['merged_face']==f and pr[j]['canonical_prism_core']==4 and pr[j]['prism_dimension']==3
        else:assert q['provenance_kind']=='unique_saved_context' and f not in pmap and q['provenance'] is None
        assert F(row['mass'])==alpha[f]==F(q['mass'])>0 and F(row['demand_mass'])==3*alpha[f]
        co=Counter();emb=Counter();assert len(row['actual_corners'])==6
        assert {tuple(c['omitted']) for c in row['actual_corners']}==set(templates)
        for corner in row['actual_corners']:
            temp=templates[tuple(corner['omitted'])];vf=tuple(sorted(row['spectator_forest']+temp['finite_components']))
            v=vertex[vf];assert v==corner['vertex'];co[v]+=1;byv[v].add(f)
            need=Counter(temp['finite_components']);have=Counter(vf);count=1
            for name,n in need.items():
                for k in range(n):count*=have[name]-k
                count*=temp['component_automorphisms'][name]**n
            assert count==corner['whole_component_embedding_count']>0;emb[v]+=count
        assert sorted(co.items())==list(map(tuple,row['corner_coefficients']))
        assert sorted(emb.items())==list(map(tuple,row['embedding_count_upper_by_vertex']))
        rows[f]=co;embedding_rows[f]=emb
    assert Counter(r['provenance_kind'] for r in d['faces'])==d['counts']['provenance_kinds']=={'explicit':18,'unique_saved_context':71}
    assert (len(byv),len(set(byv)&set(latest)),len(set(byv)-set(latest)))==(238,183,55)
    def gain(fs):
        e=Counter()
        for f in fs:
            for v,n in rows[f].items():e[v]+=n*alpha[f]
        return sum((max(F(),t.get(v,F())+x-mu.get(v,F())) for v,x in e.items()),F())
    # Independently derive connected components by union-find on shared forests.
    parent={f:f for f in rows}
    def find(f):
        while parent[f]!=f:f=parent[f]
        return f
    for fs in byv.values():
        fs=sorted(fs)
        for f in fs[1:]:parent[find(f)]=find(fs[0])
    comp=defaultdict(list)
    for f in rows:comp[find(f)].append(f)
    comps=sorted((sorted(fs) for fs in comp.values()),key=lambda fs:fs[0])
    assert comps==[c['faces'] for c in d['component_certificates']] and len(comps)==64 and max(map(len,comps))==7
    previous={0:(F(),())};processed=set();subsets=inequalities=0
    assert len(d['cardinality_stages'])==len(comps)
    for cc,certificate,stage in zip(comps,d['component_certificates'],d['cardinality_stages']):
        expectedsets=[fs for k in range(len(cc)+1) for fs in combinations(cc,k)]
        assert [tuple(x['faces']) for x in certificate['complete_subsets']]==expectedsets
        local={}
        for fs,entry in zip(expectedsets,certificate['complete_subsets']):
            z=gain(fs);assert z==F(entry['gain']);subsets+=1;k=len(fs)
            if k not in local or z>local[k][0]:local[k]=(z,fs)
        declared={r['count']:r for r in certificate['cardinality_maxima']};assert set(declared)==set(local)
        for k,r in declared.items():assert F(r['gain'])==local[k][0]==gain(r['faces']) and len(r['faces'])==k
        current={r['count']:(F(r['gain']),tuple(r['faces'])) for r in stage};processed.update(cc)
        assert len(current)==len(stage) and set(current)==set(range(len(processed)+1))
        # Every split supplies an upper inequality; each recorded actual subset attains.
        for a,(x,fs) in previous.items():
            for k,(y,gs) in local.items():assert current[a+k][0]>=x+y;inequalities+=1
        for k,(z,fs) in current.items():
            assert len(fs)==len(set(fs))==k and set(fs)<=processed and z==gain(fs)
        previous=current
    assert subsets==d['counts']['enumerated_component_subsets']==318
    assert set(previous)==set(range(90))
    num,den=F(b['objective_numerator']),F(b['normalizer']);need=num/12-den
    assert str(need)==d['needed_mu_increase_for_ratio12']
    assert [r['count'] for r in d['cardinality_frontier']]==list(range(90))
    for r in d['cardinality_frontier']:
        k=r['count'];assert F(r['gain'])==previous[k][0] and tuple(r['faces'])==previous[k][1]
        assert F(r['ratio_lower_bound'])==num/(den+previous[k][0])
    minimum=min(k for k,(z,fs) in previous.items() if z>need);assert minimum==d['minimum_cardinality_permitting_ratio_below12']==42
    assert num/(den+previous[41][0])>12>num/(den+previous[42][0])
    assert [r['count'] for r in d['exact_mass_prefix']]==list(range(1,90))
    ordered=[r['face'] for r in d['faces']]
    for r in d['exact_mass_prefix']:
        z=gain(ordered[:r['count']])
        assert F(r['gain'])==z and F(r['ratio_lower_bound'])==num/(den+z)
    for name,fs in [('all89',sorted(rows)),('best_minimum_batch',d['best_minimum_batch']['faces'])]:
        r=d[name];z=gain(fs);assert F(r.get('gain',r.get('maximum_possible_mu_increase')))==z
        e=Counter();emb=Counter()
        for f in fs:
            for v,n in rows[f].items():e[v]+=n*alpha[f]
            emb.update(embedding_rows[f])
        assert len(r['forests'])==len(e) and {x['vertex'] for x in r['forests']}==set(e)
        pairupper=0
        for x in r['forests']:
            v=x['vertex'];assert F(x['actual_corner_demand'])==e[v]
            assert F(x['baseline_local_cover'])==t.get(v,F()) and F(x['baseline_vertex_mass'])==mu.get(v,F())
            upper=max(mu.get(v,F()),t.get(v,F())+e[v]);assert F(x['possible_vertex_mass_upper'])==upper
            assert F(x['maximum_possible_mu_increase'])==upper-mu.get(v,F())
            if v in latest:
                tag,old=latest[v];assert x['baseline_record']==dict(source=tag,record_sha256=h(old))
            else:assert x['baseline_record'] is None and x['old_literal_stars']==0
            n=x['candidate_new_stars_upper'];assert n==emb[v] and x['old_literal_stars']==oldn.get(v,0)
            np=x['old_literal_stars']*n+n*(n-1)//2
            assert np==x['candidate_new_pair_tests_upper'];pairupper+=np
        assert pairupper==r['candidate_pair_tests_upper'] and F(r['ratio_lower_bound'])==num/(den+z)
    assert F(d['all89']['remaining_saving_to_restore_barrier'])==previous[89][0]-need==F(8639262097,665280000000000)
    assert all(d[k]==0 for k in ('new_physical_pair_tests','new_clique_enumerations','covering_LPs','global_model_optimizations'))
    assert d['actual_represented_classes']==324 and d['actual_represented_forests']==600 and d['new_dimension_bound'] is None
    print('PASS_ALL89_SAVED_CARDINALITY_BARRIER',dict(classes=89,corners=534,components=64,subsets=subsets,upper_inequalities=inequalities,minimum_cardinality=minimum,all89_ratio=d['all89']['ratio_lower_bound']))

if __name__=='__main__':run()
