#!/usr/bin/env python3
"""Exact diagnostic audit of nine frozen rank21 restricted weight models.

Rechecks literal rational weak duality and selected finite-graph incidence
facts. Does not duplicate the producers' face-coefficient reconstruction.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_dimension21_h3_cross_20261016 as H

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_dimension21_diagnostic_synthesis_20261016'
STEMS=['h3_cross','e7_context','e8_context','e7e8_context','e7e8_exact','a4_cross','a4_exceptional','a4_pairs','a4_pair_exceptional']
G=H.G

def floor_decimal(v,den=10000):return str(F(v.numerator*den//v.denominator,den))
def weak_dual(d,drop=()):
    excluded=set(drop);total=Counter();value=F(0);z=d['variables']
    for i,t in enumerate(d['certificate']['terms']):
        if i in excluded:continue
        y=F(t['multiplier']);assert y<=0
        value+=y*F(t['rhs'])
        for j,c in t['coefficients']:total[j]+=y*F(c)
    assert all(F(j==z)-v>=0 for j,v in total.items())
    return value,len(total)

def corner_pairs(f):
    r=f['row'];es=dict(zip(combinations(range(r['nodes']),2),r['labels']));out=[]
    for a,b in product((0,2),(1,3)):
        parts=G.pieces(G.graph(r['white']+[a,b],es))
        left=next(c for c in parts if a in c);right=next(c for c in parts if b in c)
        if b not in left:out.append(tuple(sorted((G.family(left),G.family(right)))))
    return out

def delete_root(g,v):
    rest={u:{w:l for w,l in row.items() if w!=v} for u,row in g.items() if u!=v}
    return tuple(sorted(G.family(c) for c in G.pieces(rest)))

def finite_extension_diagnostic(forests):
    aliases={name:G.family(g) for name,g in G.models() if len(g)<=21}
    graphs={G.family(g):g for name,g in G.models() if len(g)<=21}
    assert len(graphs)==len(aliases)
    forests=[tuple(sorted(aliases[k] for k in forest)) for forest in forests]
    ranks={name:len(g) for name,g in graphs.items()}
    patterns=defaultdict(set)
    for name,g in graphs.items():
        assert G.family(g)==name
        for v in g:patterns[delete_root(g,v)].add(name)
    # No irreducible finite diagram in this complete rank range can absorb
    # an E8 component when one new root is added. E7 can only become E8.
    assert all('E8' not in p for p in patterns)
    e7_patterns=[(p,sorted(names)) for p,names in patterns.items() if 'E7' in p]
    assert e7_patterns==[(('E7',),['E8'])]
    def extensions(s):
        answer=set();cs=Counter(s)
        for pat,names in patterns.items():
            cp=Counter(pat)
            if any(cp[k]>cs[k] for k in cp):continue
            rest=list((cs-cp).elements())
            for name in names:
                t=tuple(sorted(rest+[name]));assert sum(ranks[k] for k in t)==21
                answer.add(t)
        return sorted(answer)
    rows=[]
    for forest in forests:
        assert sum(ranks[k] for k in forest)==21
        deletions=Counter()
        for i,name in enumerate(forest):
            for v in graphs[name]:
                s=tuple(sorted(list(forest[:i]+forest[i+1:])+list(delete_root(graphs[name],v))))
                deletions[s]+=1
        cases=[]
        for s,m in sorted(deletions.items()):
            ts=extensions(s);assert tuple(sorted(forest)) in ts
            assert all(t.count('E8')>=s.count('E8') for t in ts)
            cases.append(dict(edge_diagram=s,incident_root_count=m,finite_endpoint_types=ts,
                original_endpoint_type_permitted=True))
        rows.append(dict(vertex_type=forest,root_deletion_patterns=cases))
    dp=[0]*22;dp[0]=1
    for name in graphs:
        for n in range(ranks[name],22):dp[n]+=dp[n-ranks[name]]
    assert dp[21]==311574
    return dict(finite_component_types=len(graphs),all_rank21_finite_forests=dp[21],
        E8_cannot_join_any_connected_finite_one_root_extension=True,
        E7_connected_one_root_extension_patterns=e7_patterns,vertex_deletion_cases=rows,
        scope='Complete finite-diagram classification through rank21 only. An allowed finite endpoint type is not a realized Lorentzian/arithmetic edge. The original vertex type is always a finite reextension, so these finite checks alone exclude none of the listed vertices.')

def hub_source_diagnostic(faces):
    graphs=dict(G.models());seen=defaultdict(set);counts={}
    for hub in ('H3','A4'):
        counts[hub]=Counter(json.dumps(G.code(graphs[hub],(r,)),separators=(',',':')) for r in graphs[hub])
    for f in faces:
        row=f['row'];es=dict(zip(combinations(range(row['nodes']),2),row['labels']))
        for a,b in product((0,2),(1,3)):
            parts=G.pieces(G.graph(row['white']+[a,b],es))
            left=next(c for c in parts if a in c);right=next(c for c in parts if b in c)
            if b in left:continue
            lf,rf=G.family(left),G.family(right)
            for hub in ('H3','A4'):
                excluded=('A1','I2(5)','H3') if hub=='H3' else ('A1','I2(5)','H3','A4')
                if lf==hub and rf not in excluded:source,sroot,dest,droot=left,a,right,b
                elif rf==hub and lf not in excluded:source,sroot,dest,droot=right,b,left,a
                else:continue
                key=(hub,G.family(dest),json.dumps(G.code(dest,(droot,)),separators=(',',':')))
                seen[key].add(json.dumps(G.code(source,(sroot,)),separators=(',',':')))
    rows=[]
    for hub in ('H3','A4'):
        destinations=[dict(destination_type=k[1],destination_root_orbit=k[2],
            occurring_source_orbits=sorted(v),occurring_source_root_multiplicity=sum(counts[hub][z] for z in v),
            all_source_roots=len(graphs[hub])) for k,v in sorted(seen.items()) if k[0]==hub]
        histogram=Counter(len(v['occurring_source_orbits']) for v in destinations)
        assert histogram==(Counter({3:6}) if hub=='H3' else Counter({1:44,2:31}))
        rows.append(dict(hub=hub,source_root_orbit_multiplicities=sorted(counts[hub].values()),
            destination_orbits=len(destinations),source_orbits_per_destination_histogram=sorted(histogram.items()),
            destinations=destinations))
    return dict(hubs=rows,
        proposed_rule='Replace v_destination/(k*c) by v_source_orbit,destination/(k*c). Summing all c hub components contributes m_source_orbit/k times the new bonus per destination root, where k is hub rank and m_source_orbit its root multiplicity. Thus existing spectator counts and additive forest optimization suffice.',
        scope='Exact census of the 1979 complete-model face rows, not an optimization or a complete global face-coverage theorem. For 44 A4 destination orbits only one of two size-two source orbits occurs in these rows; setting the absent source bonus to zero halves the corresponding vertex contribution while preserving these input payments.')

def run():
    base=H.W.build();faces=base['faces'];assert len(faces)==1979
    paircache={i:corner_pairs(f) for i,f in enumerate(faces)}
    models=[];sources={};forest_recurrence=Counter();face_recurrence=Counter();last_forests=[]
    for stem in STEMS:
        tag=f'fresh_sqrt5_dimension21_{stem}_20261016';path=ROOT/f'results/{tag}_cap2.json'
        d=json.loads(path.read_text());assert d['dimension']==21 and d['target']=='40'
        assert d['certificate']['kind']=='exact_model_obstruction'
        lower,cols=weak_dual(d);assert lower==F(d['certificate']['lower_bound'])>40
        terms=d['certificate']['terms'];fts=[];face_mass=defaultdict(F);drops=[];unpaid=Counter()
        has_a4=stem.startswith('a4_');hubs=('A1','I2(5)','H3','A4') if has_a4 else ('A1','I2(5)','H3')
        for i,t in enumerate(terms):
            m=t['metadata'];y=F(t['multiplier'])
            if m['kind']=='forest':
                fts.append((y,tuple(sorted(m['types'])),i));continue
            assert all(F(c)<=0 for j,c in t['coefficients'])
            first=m['witnesses'][0];fi=first['face'];f=faces[fi]
            assert F(t['rhs'])==sum((F(c['baseline']) for c in f['corners']),F(0))-1
            ids={(faces[w['face']]['family'],faces[w['face']]['index']) for w in m['witnesses']}
            for k in ids:face_mass[k]-=y
            un=set()
            for w in m['witnesses']:
                for pair in paircache[w['face']]:
                    if not any(k in hubs for k in pair):un.add(pair)
            if un:
                drops.append(i)
                for pair in un:unpaid[pair]+=1
        topforests=sorted(fts)[:8]
        for _,f,_ in topforests:forest_recurrence[f]+=1
        topfaces=sorted(face_mass.items(),key=lambda z:(-z[1],z[0]))[:8]
        for f,_ in topfaces:face_recurrence[f]+=1
        rec=dict(prefix=tag,dimension=21,target='40',exact_lower_bound=str(lower),
            lower_bound_floor_2dp=floor_decimal(lower,100),literal_dual_terms=len(terms),
            literal_columns_checked=cols,face_rows=d['faces'],common_spectator_cases=d['face_cases'],
            variables=d['variables'],unique_face_constraints=d['unique_face_constraints'],
            largest_forest_multipliers=[dict(term=i,types=f,multiplier_floor=floor_decimal(-y)) for y,f,i in topforests],
            largest_face_multipliers=[dict(family=k[0],index=k[1],total_multiplier_floor=floor_decimal(v)) for k,v in topfaces])
        if has_a4:
            trimmed,tcols=weak_dual(d,drops)
            rec['unpaid_pair_trim']=dict(dropped_term_indices=drops,dropped_terms=len(drops),
                unpaid_pair_types=[dict(types=k,face_terms=v) for k,v in sorted(unpaid.items())],
                exact_trimmed_lower_bound=str(trimmed),trimmed_lower_floor_2dp=floor_decimal(trimmed,100),
                remains_above_target=trimmed>40,literal_columns_checked=tcols,
                extension_scope='This trimmed dual remains valid upon adding nonnegative variables whose coefficients vanish on every retained face and are nonnegative on every forest. In particular it obstructs adding bonuses only for previously unpaid disconnected marked component pairs while leaving the old contexts and coefficients unchanged. It does not cover context refinements or new root dependence on already-paid pair families.')
        models.append(rec)
        if stem=='a4_pair_exceptional':last_forests=[f for _,f,_ in topforests[:4]]
        for ext,folder in [('.py','scripts'),('_cap2.json','results'),('.md','results')]:
            q=f'{folder}/{tag}{ext}';sources[q]=sha256((ROOT/q).read_bytes()).hexdigest()
        print(stem,'lower',floor_decimal(lower,100),'trim',rec.get('unpaid_pair_trim',{}).get('trimmed_lower_floor_2dp'),flush=True)
    # Bind and check complete recursive manifests without re-running old work.
    seen=set();refs=0;queue=[p for p in sources if p.endswith('.json')]
    while queue:
        p=queue.pop()
        if p in seen:continue
        seen.add(p);d=json.loads((ROOT/p).read_text())
        for q,h in d.get('sources',{}).items():
            assert sha256((ROOT/q).read_bytes()).hexdigest()==h,(p,q);refs+=1
            if q.endswith('.json'):queue.append(q)
    sp=f'scripts/{TAG}.py';sources[sp]=sha256((ROOT/sp).read_bytes()).hexdigest()
    return dict(status='PASS_EXACT_NINE_RESTRICTED_MODELS_DIAGNOSTIC_SYNTHESIS',sources=sources,
        models=models,recurring_top_forest_types=[dict(types=k,models=v) for k,v in forest_recurrence.most_common(12)],
        recurring_top_face_rows=[dict(family=k[0],index=k[1],models=v) for k,v in face_recurrence.most_common(12)],
        finite_endpoint_diagnostic=finite_extension_diagnostic(list(dict.fromkeys(last_forests+[k for k,v in forest_recurrence.most_common(2)]))),
        hub_source_orbit_diagnostic=hub_source_diagnostic(faces),
        recursive_manifests_checked=len(seen),recursive_source_references_checked=refs,
        scope='Independent literal rational weak-duality checks and finite-diagram extension diagnostics; relies on the frozen producers/replays for exact face and forest coefficient reconstruction. No improved dimension bound, nonexistence theorem, global incidence realization or finite-volume example is claimed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
