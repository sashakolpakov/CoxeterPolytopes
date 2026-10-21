#!/usr/bin/env python3
"""Independent Fraction checks of saved raw4 capacity and repair accounting.

No project imports, field arithmetic, graph enumeration or optimizer. Physical
exclusions and complete count-vector domains remain producer theorem inputs.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_new5121_raw4_two_face_local_covers_20261018'
BASE='fresh_cross_vertex_historical32_local_covers_20261018'
OLD='fresh_cross_vertex_current290_local_covers_20261018'
RAW='fresh_new5121_raw4_two_face_domain_20261018'
MODEL='fresh_cross_vertex_current_full_model_20261018'


def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def okey(x):return x['face'],x['provenance'],tuple(x['omitted']),x['vertex'],x['raw_family']


def run():
    d=read(TAG);b=read(BASE);q=read(RAW);old=read(OLD);model=read(MODEL)
    assert d['status']=='PASS_NEW5121_RAW4_TWO_FACE_EXACT_LOCAL_COVERS'
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert d['face_multipliers']==b['face_multipliers'] and d['initial_vertex_multipliers']==b['vertex_multipliers']
    alpha={i:F(x) for i,x in d['face_multipliers']};mu={i:F(x) for i,x in d['initial_vertex_multipliers']}
    after=dict(mu);rows32={r['vertex']:r for r in b['forest_records']};rows290={r['vertex']:r for r in old['forest_records']}
    qr={r['vertex']:r for r in q['forest_records']};rows=d['forest_records']
    assert len(rows)==len(qr)==6 and {r['vertex'] for r in rows}==set(qr)
    counts=Counter();columns=[];inequalities=terms=0
    for r in rows:
        v=r['vertex'];u=r['domain'];groups=u['corner_orbits'];nodes=u['literal_stars'];oldn=u['old_literal_stars']
        base=rows32.get(v,rows290.get(v));assert jhash(base)==u['baseline']['record_sha256']
        oldcoeff=[z['positive_face_coefficients'] for z in base['domain']['corner_orbits']] if v in rows32 else base['orbit_coefficients']
        assert len(oldcoeff)==u['baseline']['corner_orbits']
        assert [z['index'] for z in groups]==list(range(len(groups)))
        owner={i:z['index'] for z in groups for i in z['literal_stars']}
        assert set(owner)==set(range(len(nodes))) and sum(len(z['literal_stars']) for z in groups)==len(nodes)
        assert len({json.dumps(s['contacts'],sort_keys=True) for s in nodes})==len(nodes)
        mapping=u['raw4_to_merged_star_indices'];raw=qr[v]
        assert len(mapping)==len(set(mapping))==len(raw['literal_stars'])
        for s,i in zip(raw['literal_stars'],mapping):
            for k,x in s['native'].items():assert nodes[i][k]==x
        assert len(u['added_occurrences'])==len(raw['occurrences'])==2
        extras=Counter();origin_keys=set()
        for a,o in zip(u['added_occurrences'],raw['occurrences']):
            assert a['origin']==o['origin'];ids=sorted(mapping[i] for i in o['literal_star_indices'])
            assert a['literal_star_indices']==ids==groups[a['orbit']]['literal_stars']
            assert len(ids)==a['origin']['literal_orientations']
            k=okey(a['origin']);assert k not in origin_keys;origin_keys.add(k)
            extras[a['orbit'],a['origin']['face']]+=1
        allorigins=set()
        for j,z in enumerate(groups):
            cc=Counter(dict(oldcoeff[j])) if j<len(oldcoeff) else Counter()
            cc.update({f:n for (i,f),n in extras.items() if i==j})
            assert cc==Counter(dict(z['positive_face_coefficients']))==Counter(o['face'] for o in z['origins'])
            for o in z['origins']:
                k=okey(o);assert k not in allorigins;allorigins.add(k)
        assert origin_keys<=allorigins
        expected=[(i,j) for i,j in combinations(range(len(nodes)),2) if j>=oldn]
        assert [tuple(z[:2]) for z in r['newly_owned_pair_domain']]==expected
        bad={tuple(z[:2]):z for z in r['incompatibility_certificates']};assert len(bad)==r['new_exclusions']
        reasons=Counter()
        for i,j,code in r['newly_owned_pair_domain']:
            if code is None:assert (i,j) not in bad
            else:
                assert bad[i,j][2]==code;reasons[str(code)]+=1
                if code==0:assert set(map(tuple,nodes[i]['occupied_pairs']))&set(map(tuple,nodes[j]['occupied_pairs']))
                else:assert 0<=bad[i,j][5]<len(r['wall_expressions'])
        assert dict(reasons)==r['reason_counts'] and len(expected)-len(bad)==r['new_unknown_pairs']
        vectors=r['maximal_count_vectors'];reps=r['literal_representatives'];c=r['certificate']
        assert len(vectors)==len(reps)>0 and all(len(a)==len(groups) for a in vectors)
        for vec,cl in zip(vectors,reps):
            assert len(cl)==len(set(cl)) and all(i in owner for i in cl)
            assert vec==[Counter(owner[i] for i in cl).get(j,0) for j in range(len(groups))]
            assert not any(tuple(sorted(p)) in bad for p in combinations(cl,2))
        demands=[sum((alpha.get(f,F())*F(m) for f,m in z['positive_face_coefficients']),F()) for z in groups]
        assert list(map(str,demands))==u['demands']
        w=list(map(F,c['weights']));assert len(w)==len(groups) and all(x>=0 for x in w)
        assert all(sum((F(a)*x for a,x in zip(vec,w)),F())<=1 for vec in vectors);inequalities+=len(vectors)
        mixture=[(j,F(x)) for j,x in c['mixture']]
        assert len({j for j,x in mixture})==len(mixture) and all(0<=j<len(vectors) and x>0 for j,x in mixture)
        assert all(sum((F(vectors[j][i])*x for j,x in mixture),F())>=demands[i] for i in range(len(groups)))
        value=sum((x for j,x in mixture),F());terms+=len(mixture)
        assert value==sum((x*y for x,y in zip(w,demands)),F())==F(c['optimal_mass'])==F(r['exact_local_optimum'])
        assert value>=F(u['baseline']['exact_local_optimum']) and F(r['baseline_vertex_mass'])==mu.get(v,F())
        after[v]=max(mu.get(v,F()),value);assert F(r['required_increase'])==after[v]-mu.get(v,F())
        co=Counter()
        for x,z in zip(w,groups):
            for f,m in z['positive_face_coefficients']:co[f]+=x*F(m)
        coeffs=[[f,str(x)] for f,x in sorted(co.items()) if x]
        assert sum((alpha.get(f,F())*F(x) for f,x in coeffs),F())==value
        if value>mu.get(v,F()):columns.append(dict(vertex=v,vertex_cost='1',positive_face_coefficients=coeffs,orbit_weights=c['weights'],exact_load=str(value)))
        counts.update(new_literal_stars=len(nodes)-oldn,new_actual_corners=2,new_pair_tests=len(expected),new_exclusions=len(bad),new_unknown_pairs=len(expected)-len(bad))
    assert dict(counts)==d['counts'] and columns==d['columns']
    inherited=d['inherited_forests'];assert len(inherited)==594
    assert len({r['vertex'] for r in inherited}|set(qr))==600 and not {r['vertex'] for r in inherited}&set(qr)
    for r in inherited:
        prev=rows32[r['vertex']] if r['source']==BASE else rows290[r['vertex']]
        assert jhash(prev)==r['record_sha256']
    assert [[v,str(x)] for v,x in sorted(after.items()) if x]==d['vertex_multipliers']
    assert all(after.get(v,F())>=x for v,x in mu.items())
    multi=[]
    for col in model['prepared']['new_columns']:
        load=sum((alpha.get(i,F())*m for i,m in col['integer_face_coefficients']),F())
        before=sum((mu.get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        final=sum((after.get(v,F())*m for v,m in col['integer_vertex_capacities']),F())
        assert load<=before<=final
        if len(col['integer_vertex_capacities'])>1:multi.append(dict(variable=col['variable'],all_vertex_costs=col['integer_vertex_capacities'],load=str(load),old_cost=str(before),new_cost=str(final)))
    assert len(multi)==14 and multi==d['preserved_multi_owner_columns']
    num=F(b['objective_numerator']);den=sum(after.values(),F())
    assert F(d['objective_numerator'])==num and F(d['normalizer'])==den and F(d['bound'])==num/den
    assert F(d['baseline_bound'])==F(b['bound'])
    assert (d['represented_classes'],d['positive_represented_classes'],d['represented_forests'])==(324,261,600)
    return dict(status='PASS_RAW4_SAVED_FRACTION_ACCOUNTING',forests=6,counts=dict(counts),
        vector_inequalities=inequalities,mixture_terms=terms,exported_columns=len(columns),
        inherited_forests=594,preserved_two_owner_columns=14,bound=str(num/den),
        scope='Independent saved multiplicity, complete new-pair index accounting, count-vector charge/cover, export and monotone-repair arithmetic only. No historical ownership/graph reconstruction, physical field gates, clique enumeration or optimizer.')


if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
