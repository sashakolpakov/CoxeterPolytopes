#!/usr/bin/env python3
"""Order12 support diagnostic including all66 low-core chosen corners.

No proposed pair removal is promoted. Universal short8 provenance is kept
separate from selected high/low patterns and immediate D13 alternatives.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
import fresh_ideal5_small_factor_payments_20261017 as S
import fresh_distance8_low_core_20261017 as L

P=S.P;ROOT=P.ROOT
TAG='fresh_order12_payment_support_20261017'
OVERLAP='fresh_distance8_high_overlap_20261017'
D13='fresh_affine_d13_long_pair_transfer_20261017'
LOW=L.TAG

def labels(g):
    vals={s.Integer(-1):3,-s.sqrt(2):4,-(1+s.sqrt(5))/2:5,-s.sqrt(3):6}
    return [[0 if i==j or not g[i,j]else vals[s.simplify(g[i,j])]for j in range(g.rows)]for i in range(g.rows)]

def immediate_domains(targets,source):
    edges=[(i,i+1,3)for i in range(9)]+[(0,10,3),(0,11,3),(9,12,3),(9,13,3)]
    rows=[]
    for name,node,label,expected in(('inner4',1,4,'B'),('inner_next3',2,3,'D')):
        saved=next(r for r in source['formulas']['whole_affine12_payments']if r['profile']==name)
        ids=saved['nodes'];assert len(ids)==saved['order']==12
        full=P.labelled(15,edges+[(14,node,label)])
        pattern=[[full[i][j]for j in ids]for i in ids]
        for kind,target in targets.items():
            maps=P.embeddings(pattern,target)
            if kind!=expected:assert not maps;continue
            assert maps
            ds=P.distance_data(target)[1]
            for b in(9,12,13):
                marks=[ids.index(b),ids.index(14)]
                pairs={tuple(sorted(m[i]for i in marks))for m in maps}
                rows.append(dict(source_profile=name,target_kind=kind,nodes=ids,literal_marks=[b,14],
                    labels=pattern,embeddings=maps,marked_images=[list(p)for p in sorted(pairs)],
                    distances=sorted({ds[a][b]for a,b in pairs}),inherited_positive_kernel=saved['positive_kernel']))
    assert len(rows)==6
    return rows

def low_corners(doc,targets):
    assert len(doc['extension_closure']['cases'])==len(doc['extension_geometry']['cases'])==66
    cases={json.dumps([r['low'],r['high'],r['match']]):r for r in doc['extension_geometry']['cases']}
    out=[];cache={};counts=Counter()
    for index,row in enumerate(doc['extension_closure']['cases']):
        old=cases[json.dumps([row['low'],row['high'],row['match']])]
        assert old['kind']=='fixed_core'and all(x['kind']=='actual_edge_PD_rejection'for x in old['extra_bonds'])
        corner=next(c for c in old['corners']if c['marks']==row['corner'])
        ids=corner['component'];assert ids and set(row['corner'])<=set(ids)
        block=next(b for b in corner['blocks']if b['nodes']==ids)
        core=doc['extension_geometry']['cores'][row['core_key']]
        full=s.Matrix([[s.sympify(x)for x in z]for z in core['gram']]);g=full.extract(ids,ids)
        variants=[('existing',g)]
        if row['method']=='unique_sealed_label3_root':
            node=ids.index(row['contact_node']);gg=2*s.eye(len(ids)+1);gg[:len(ids),:len(ids)]=g
            gg[node,len(ids)]=gg[len(ids),node]=-1;variants.append(('one_sealed_root',gg))
        saved=[]
        for name,m in variants:
            key=str(m.tolist())
            if key not in cache:cache[key]=L.inertia_small(m)
            sig=cache[key]['signature'];assert sig[1]==0 and sig[2]in(0,1)
            assert (sig[2]==1)==(block['kind']=='affine')if name=='existing'else True
            assert m.rows<=row['component_order_bound']<=12
            gg=labels(m);marks=[ids.index(i)for i in row['corner']]
            affine_kinds=[];target_maps={}
            if sig[2]:
                if name=='existing':kernel=s.Matrix(list(map(s.sympify,block['positive_kernel'])))
                else:
                    null=m.nullspace();assert len(null)==1;kernel=null[0]
                    if kernel[0]<0:kernel=-kernel
                assert m*kernel==s.zeros(m.rows,1)and all(x>0 for x in kernel)
                for kind in'ABCD':
                    same=P.labelled(m.rows,P.A.P.graph(kind,m.rows)[0])
                    if P.embeddings(gg,same):affine_kinds.append(kind)
                assert len(affine_kinds)==1
                if m.rows==12:
                    for kind,target in targets.items():
                        maps=P.embeddings(gg,target)
                        if maps:target_maps[kind]=dict(embeddings=maps,marked_images=[list(p)for p in sorted({tuple(sorted(a[i]for i in marks))for a in maps})])
                    assert set(target_maps)==set(affine_kinds)
            saved.append(dict(variant=name,order=m.rows,labels=gg,marks=marks,signature=sig,
                affine_kind=affine_kinds[0]if affine_kinds else None,order12_target_maps=target_maps))
            counts[('affine'if sig[2]else'finite')+'_'+str(m.rows)]+=1
        out.append(dict(index=index,low=row['low'],high=row['high'],match=row['match'],
            method=row['method'],chosen_marks=row['corner'],component_nodes=ids,
            component_order_bound=row['component_order_bound'],variants=saved))
    assert len(out)==66 and sum(len(r['variants'])for r in out)==69
    print('LOW_ACTUAL_VARIANTS',dict(counts),flush=True)
    return dict(cases=out,counts=dict(counts),unique_exact_grams=len(cache),
        scope='All66 frozen chosen-corner closures and all3 optional sealed simple-root extensions. Existing affine factors are whole by their positive kernel; frozen closure forbids any further relevant supporting root. Actual optional variants are checked exactly, rather than embedding finite proper patterns into arbitrary affine completions.')

def run():
    stems=[OVERLAP,P.CENSUS,S.TAG,D13,LOW,'fresh_finite8_ideal5_combined_payments_20261017',
        'fresh_finite8_ideal5_local_replacement_20261017','fresh_cycle13_ideal5_transfer_20261017']
    inputs={tag:S.read(tag)for tag in stems}
    for tag,data in inputs.items():
        assert data['status'].startswith('PASS')
        for path,h in data['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,(tag,path)
    over=inputs[OVERLAP];diags={r['index']:r for r in inputs[P.CENSUS]['diagrams']}
    assert len(over['isolation']['cases'])==2614
    patterns=[]
    for r in over['isolation']['cases']:
        if r.get('component_order_bound',0)<12 or len(r['marked_component'])>12:continue
        g,_=P.R.unite(P.R.bond_graph(diags[r['left'][0]]),P.R.bond_graph(diags[r['right'][0]]),tuple(map(tuple,r['match'])))
        ids=r['marked_component'];gg=[[int(2-g[i][j])if i!=j and g[i][j]else 0 for j in ids]for i in ids]
        patterns.append(dict(index=r['index'],left=r['left'],right=r['right'],match=r['match'],nodes=ids,labels=gg,
            marks=[ids.index(i)for i in r['preferred_corner']],bound=r['component_order_bound'],method=r['method']))
    targets={k:P.labelled(12,P.A.P.graph(k,12)[0])for k in'ABCD'}
    low=low_corners(inputs[LOW],targets);immediate=immediate_domains(targets,inputs[D13]);cache={};records=[]
    for kind,g in targets.items():
        basic=S.basic_pairs(g);autos=P.embeddings(g,g);ds=P.distance_data(g)[1];rows=[];union=set()
        for p in patterns:
            key=(kind,str(p['labels']))
            if key not in cache:cache[key]=P.embeddings(p['labels'],g)
            maps=cache[key];pairs={tuple(sorted(a[i]for i in p['marks']))for a in maps}
            if maps:rows.append(dict(case=p['index'],embeddings=maps,marked_images=[list(p)for p in sorted(pairs)],
                left=p['left'],right=p['right'],match=p['match'],method=p['method']))
            union|=pairs
        immediate_pairs={tuple(p)for r in immediate if r['target_kind']==kind for p in r['marked_images']}
        low_pairs=set();low_refs=[]
        for r in low['cases']:
            for z in r['variants']:
                if kind in z['order12_target_maps']:
                    maps=z['order12_target_maps'][kind];low_pairs.update(map(tuple,maps['marked_images']))
                    low_refs.append(dict(case=r['index'],variant=z['variant'],**maps))
        selected=union|low_pairs|immediate_pairs;todo=selected-basic;orbits=[]
        while todo:
            pair=min(todo);orb={tuple(sorted(a[i]for i in pair))for a in autos};assert orb<=todo;todo-=orb
            orbits.append(dict(representative=list(pair),pairs=[list(p)for p in sorted(orb)],distance=ds[pair[0]][pair[1]],
                high_cases=[r['case']for r in rows if set(map(tuple,r['marked_images']))&orb],
                low_cases=[r['case']for r in low_refs if set(map(tuple,r['marked_images']))&orb],
                immediate_profiles=sorted({r['source_profile']for r in immediate if r['target_kind']==kind and set(map(tuple,r['marked_images']))&orb})))
        short8={p for p in combinations(range(12),2)if ds[p[0]][p[1]]<=8}
        missing_short=short8-basic
        required=basic|selected
        records.append(dict(kind=kind,order=12,target=g,basic_pairs=[list(p)for p in sorted(basic)],
            selected_high_pairs=[list(p)for p in sorted(union)],selected_low_pairs=[list(p)for p in sorted(low_pairs)],
            immediate_pairs=[list(p)for p in sorted(immediate_pairs)],required_by_selected_interfaces=[list(p)for p in sorted(required)],
            short8_pairs=[list(p)for p in sorted(short8)],short8_outside_basic=[list(p)for p in sorted(missing_short)],
            selected_plus_unreplaced_short8=[list(p)for p in sorted(required|short8)],nonbasic_orbits=orbits,
            high_case_embeddings=rows,low_case_embeddings=low_refs,embedding_count=sum(len(r['embeddings'])for r in rows),
            universal_short8_scope='Cycle12 antipodal distance6 is outside the frozen cycle13 and cycle>=15 replacements; absence from selected patterns does not discharge it.'if kind=='A'else'Every short8 pair is retained in the basic order12 set.'))
        print('ORDER12',kind,'basic/high/low/immediate/required/short-unpaid',len(basic),len(union),len(low_pairs),len(immediate_pairs),len(required),len(missing_short),
            'highcases',len(rows),'nonbasic',[(r['pairs'],r['high_cases'],r['low_cases'],r['immediate_profiles'])for r in orbits],flush=True)
    assert [len(r['basic_pairs'])for r in records]==[60,64,63,66]
    assert [len(r['short8_outside_basic'])for r in records]==[6,0,0,0]
    paths=[Path(__file__),Path(S.__file__),Path(P.__file__),Path(L.__file__)]
    paths.extend(ROOT/f'{folder}/{tag}.{ext}'for tag in stems for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return dict(status='PASS_ORDER12_COMPLETE_SELECTED_LOW_AND_IMMEDIATE_SUPPORT_DIAGNOSTIC',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        high_patterns=patterns,eligible_high_patterns=len(patterns),pattern_bound_counts=dict(Counter(r['bound']for r in patterns)),
        low_corner_domain=low,immediate_D13_alternatives=immediate,records=records,unchanged_whole_order12_cost=66,
        scope='Diagnostic only. Complete selected high support with cutoff12, all66 actual low-core chosen corners and their optional sealed extensions, literal immediate D13 affine alternatives, and separate universal short8 provenance. No pair removal, new payment theorem, ideal envelope or dimension bound.')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');x=a.parse_args();out=json.loads(json.dumps(run()))
    path=ROOT/'results'/f'{TAG}.json'
    if x.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],flush=True)
