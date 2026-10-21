#!/usr/bin/env python3
"""Exact order13 selected/immediate payment diagnostic; no new theorem."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_ideal5_small_factor_payments_20261017 as S

P=S.P;ROOT=P.ROOT
TAG='fresh_order13_payment_support_20261017'
OVERLAP='fresh_distance8_high_overlap_20261017'
D14='fresh_affine_d14_long_pair_transfer_20261017'

def immediate_domains(targets,source):
    old_edges=[(i,i+1,3)for i in range(10)]+[(0,11,3),(0,12,3),(10,13,3),(10,14,3)]
    rows=[]
    for name,node,label,expected in(('inner4',1,4,'B'),('inner_next3',2,3,'D')):
        saved=next(r for r in source['profile_formulas']['small_affine_rows']if r['profile']==name)
        ids=saved['nodes'];assert len(ids)==saved['order']==13
        full=P.labelled(16,old_edges+[(15,node,label)])
        pattern=[[full[i][j]for j in ids]for i in ids]
        for kind,target in targets.items():
            maps=P.embeddings(pattern,target)
            if kind!=expected:assert not maps;continue
            assert maps
            ds=P.distance_data(target)[1]
            for b in(10,13,14):
                marks=[ids.index(b),ids.index(15)]
                pairs={tuple(sorted(m[i]for i in marks))for m in maps}
                rows.append(dict(source_profile=name,target_kind=kind,nodes=ids,
                    literal_marks=[b,15],labels=pattern,embeddings=maps,
                    marked_images=[list(p)for p in sorted(pairs)],
                    distances=sorted({ds[a][b]for a,b in pairs}),
                    inherited_positive_kernel=saved['kernel']))
    assert len(rows)==6
    return rows

def run():
    stems=[OVERLAP,P.CENSUS,S.TAG,D14,'fresh_finite8_ideal5_combined_payments_20261017']
    inputs={s:S.read(s)for s in stems}
    for stem,data in inputs.items():
        assert data['status'].startswith('PASS'),stem
        for name,h in data['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
    over=inputs[OVERLAP];diags={r['index']:r for r in inputs[P.CENSUS]['diagrams']}
    assert len(over['isolation']['cases'])==2614
    patterns=[]
    for r in over['isolation']['cases']:
        if r.get('component_order_bound',0)<13 or len(r['marked_component'])>13:continue
        g,_=P.R.unite(P.R.bond_graph(diags[r['left'][0]]),P.R.bond_graph(diags[r['right'][0]]),tuple(map(tuple,r['match'])))
        ids=r['marked_component'];labels=[[int(2-g[i][j])if i!=j and g[i][j]else 0 for j in ids]for i in ids]
        patterns.append(dict(index=r['index'],left=r['left'],right=r['right'],match=r['match'],
            nodes=ids,labels=labels,marks=[ids.index(i)for i in r['preferred_corner']],
            bound=r['component_order_bound'],method=r['method']))
    assert len(patterns)==801
    assert Counter(r['bound']for r in patterns)=={13:556,14:108,15:87,16:50}
    targets={k:P.labelled(13,P.A.P.graph(k,13)[0])for k in 'ABCD'}
    immediate=immediate_domains(targets,inputs[D14]);cache={};records=[]
    for kind,g in targets.items():
        basic=S.basic_pairs(g);autos=P.embeddings(g,g);ds=P.distance_data(g)[1];rows=[];union=set()
        for p in patterns:
            key=(kind,str(p['labels']))
            if key not in cache:cache[key]=P.embeddings(p['labels'],g)
            maps=cache[key];pairs={tuple(sorted(a[i]for i in p['marks']))for a in maps}
            if maps:rows.append(dict(case=p['index'],embeddings=maps,
                marked_images=[list(p)for p in sorted(pairs)],left=p['left'],right=p['right'],match=p['match'],method=p['method']))
            union|=pairs
        immediate_pairs={tuple(p)for r in immediate if r['target_kind']==kind for p in r['marked_images']}
        assert immediate_pairs<=basic|union
        todo=union-basic;orbits=[]
        while todo:
            p=min(todo);orb={tuple(sorted(a[i]for i in p))for a in autos};assert orb<=todo;todo-=orb
            cases=[r['case']for r in rows if set(map(tuple,r['marked_images']))&orb]
            orbits.append(dict(representative=list(p),pairs=[list(p)for p in sorted(orb)],
                distance=ds[p[0]][p[1]],affected_cases=cases))
        required=basic|union|immediate_pairs
        missing=set(combinations(range(13),2))-required
        records.append(dict(kind=kind,order=13,target=g,basic_pairs=[list(p)for p in sorted(basic)],
            selected_pairs=[list(p)for p in sorted(union)],immediate_pairs=[list(p)for p in sorted(immediate_pairs)],
            required_union=[list(p)for p in sorted(required)],not_required_by_these_interfaces=[list(p)for p in sorted(missing)],
            orbits=orbits,case_embeddings=rows,embedding_count=sum(len(r['embeddings'])for r in rows)))
        print('ORDER13_SUPPORT',kind,len(basic),len(union),len(required),len(rows),flush=True)
    assert [len(r['basic_pairs'])for r in records]==[65,73,72,74]
    assert [len(r['selected_pairs'])for r in records]==[39,9,5,8]
    assert [len(r['required_union'])for r in records]==[65,78,73,78]
    assert [len(r['case_embeddings'])for r in records]==[67,107,13,305]
    assert [r['embedding_count']for r in records]==[1742,240,30,2344]
    assert [len(r['affected_cases'])for r in records[1]['orbits']]==[2,23,4]
    assert [len(r['affected_cases'])for r in records[3]['orbits']]==[127]
    high={'center':6,'leaf':8};low={'center':3,'leaf':4};d12=[]
    for a,b in product(('center','leaf'),repeat=2):
        resistance=8+2*(a=='leaf')+2*(b=='leaf')
        extreme=high[a]+high[b];one=max(low[a]+high[b],high[a]+low[b])
        assert extreme-resistance==4
        d12.append(dict(ends=[a,b],resistance=resistance,extreme_cost_sum=extreme,
            conditional_extreme_heights=['1/2','1/2'],required_extreme_unscaled_energy='1',
            inherited_nonextreme_cost_bound=one,inherited_nonextreme_gap=str(F(1)-F(one,resistance))))
    paths=[Path(__file__),ROOT/f'scripts/{S.TAG}.py',ROOT/f'scripts/{P.TAG}.py']
    paths += [ROOT/f'{folder}/{s}.{ext}'for s in stems for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER13_SELECTED_AND_IMMEDIATE_SUPPORT_DIAGNOSTIC',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        patterns=patterns,eligible_patterns=801,pattern_bound_counts=dict(Counter(r['bound']for r in patterns)),
        records=records,immediate_D14_alternatives=immediate,
        conditional_D12_profile_comparison=d12,
        unchanged_whole_order13_cost=78,
        scope='Exact complete selected-pattern support at order13 plus literal immediate affine alternatives from the frozen D14 transfer. No new pair set is promoted: short-path replacement, whole-factor alternatives and any changed local transfer remain separate proof obligations. B and D still require all78 pairs under the saved interfaces. No new dimension bound.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
