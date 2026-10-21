#!/usr/bin/env python3
"""Saved-role and exact-core audit for a proposed physical D17 face capacity.

Source preparation only until the coordinating agent schedules a guarded run.
No producer, optimizer or historical census is imported. The reused reader
module imports only the standard library and the resource guard. Universal
chosen-output completeness and actual core seals remain written premises.
"""
import argparse
from collections import Counter
import fcntl
from fractions import Fraction as F
import json
from pathlib import Path
import signal

import check_a15_inner_unit_a1_absence_retention_20261018 as W

ROOT = W.ROOT
TAG = 'd17_node2_face_capacity_retention_20261019'
SELECTED = 'fresh_d17_existing_three_corner_transfer_20261018'
MAPS = {
    1062: ([3,2,0,4,5,6,7,8,9,18,17,16,15,14,10,11,12], 1, 13),
    2525: ([13,12,10,14,15,16,17,18,9,8,7,6,5,4,0,1,2], 11, 3),
}


def role(raw, white, marks, ids=None):
    ids = sorted(set(white or ()) | set(marks)) if ids is None else sorted(ids)
    if len(ids) < 17:
        return None
    g = W.mat(raw)
    block = next((b for b in W.parts(g, ids) if set(marks) <= set(b)), None)
    if block is None or len(block) != 17:
        return None
    adj = {i:[j for j in block if i != j and g[i][j]] for i in block}
    centers = [i for i in block if len(adj[i]) == 3]
    if len(centers) != 1 or any(len(adj[i]) not in (1,2,3) for i in block):
        return None
    center = centers[0]
    leaves = [i for i in block if adj[i] == [center]]
    if len(leaves) != 2:
        return None
    end = [i for i in block if len(adj[i]) == 1 and i not in leaves]
    if len(end) != 1:
        return None
    order = [end[0]]
    while order[-1] != center:
        nxt = [j for j in adj[order[-1]] if j not in order]
        if len(nxt) != 1:
            return None
        order.append(nxt[0])
    if len(order) != 15 or sorted(order+leaves) != sorted(block):
        return None
    if not set(marks) <= set(order) or sorted(order.index(m) for m in marks) != [2,12]:
        return None
    for i in block:
        if g[i][i] <= 0:
            return None
        for j in adj[i]:
            if g[i][j] >= 0 or 4*g[i][j]**2 != g[i][i]*g[j][j]:
                return None
    if not W.pd(g, ids):
        return dict(nonfinite=True)
    return dict(nonfinite=False, nodes=block, canonical_spine=order,
                fork_leaves=sorted(leaves), corner_nodes=ids)


def dot(a, b):
    return sum(x*y for x,y in zip(a[:-1], b[:-1]))-a[-1]*b[-1]


def geometry():
    roots = []
    for i in range(16):
        v = [F()]*18; v[i]=1; v[i+1]=-1; roots.append(v)
    v = [F()]*18; v[15]=v[16]=1; roots.append(v)
    omega2 = [F(int(i < 3)) for i in range(17)]+[F()]
    omega_leaf = [F(1,2)]*17+[F()]
    assert [dot(omega2,a) for a in roots] == [int(i == 2) for i in range(17)]
    assert [dot(omega_leaf,a) for a in roots] == [int(i == 16) for i in range(17)]
    assert (dot(omega2,omega2),dot(omega_leaf,omega_leaf),dot(omega2,omega_leaf)) == (3,F(17,4),F(3,2))
    c=[-x for x in omega2]; c[-1]=1
    d=[-x for x in omega_leaf]; d[-1]=F(3,2)
    assert dot(c,c) == dot(d,d) == 2 and dot(c,d) == 0
    models = {}; proofs = []
    saved = {p['index']:p for p in W.load(SELECTED)['selected_proofs']}
    for index,(mp,ci,di) in MAPS.items():
        case,core = W.native(index)
        assert core['gauge'] == [1]*19 and set(mp+[ci,di]) == set(range(19))
        native = [None]*19
        for j,k in enumerate(mp):
            native[k]=roots[j]
        native[ci]=c; native[di]=d
        g = W.mat(core['gram'])
        assert all(dot(native[i],native[j]) == g[i][j] for i in range(19) for j in range(19))
        white = saved[index]['white']; marks = case['preferred_corner']
        out = role(core['gram'],white,marks)
        assert out and not out['nonfinite'] and set(out['nodes']) == set(mp)
        assert marks == [mp[2],mp[12]] or marks == [mp[12],mp[2]]
        z = list(map(F,core['radical'][0])); allowed=case['possible_contacts']
        assert all(sum(g[i][j]*z[j] for j in range(19)) == 0 for i in range(19))
        assert all(z[i] for i in allowed)
        left=case['left'][1]; right=[case['right_copy_map'][i] for i in case['right'][1]]
        assert sorted(white) == sorted(set(range(19))-set(left+right))
        assert W.pd(g,white) and all(W.pd(g,sorted(white+[i])) for i in left+right)
        for corner in case['all_corners']:
            assert W.pd(g,sorted(set(white)|set(corner['preferred_corner'])))
        du,dv=W.distances(g,left),W.distances(g,right)
        lengths=[min(du[i]+dv[j],du[j]+dv[i])+2 for i in allowed for j in allowed if z[i]*z[j] < 0]
        assert lengths and max(lengths) <= 8
        models[index] = dict(gram=core['gram'],white=white,marks=marks)
        proofs.append(dict(index=index,canonical_nodes=mp,omitted=[ci,di],
            positive_squared_norms=['3','17/4'],timelike_coefficients=['1','3/2'],
            allowed_contacts=allowed,short_bound=max(lengths)))
    support=W.load(SELECTED)['complete_selected_support']
    hits=[r for r in support['rows'] if [2,12] in r['marked_images']]
    assert support['complete_eligible_patterns'] == 140 and support['complete_maps'] == 168
    assert [r['index'] for r in hits] == [1062,2525]
    assert all(len(r['base_nodes']) == 17 for r in hits)
    assert sum(len(r['maps']) for r in hits) == 4
    # Cross-factor positive projections vanish; both omitted walls face away
    # from the same finite vertex, fixing the signs of their time coordinates.
    mixed_product=-F(3,2)
    assert mixed_product == -c[-1]*d[-1]
    # cos^2(pi/4)=1/2 <9/16; cos^2(pi/5)=(3+sqrt5)/8 >9/16
    # since sqrt5>3/2. Thus no integer Coxeter angle has cosine3/4.
    assert F(9,16)>F(1,2) and F(3,2)**2 < 5
    return models,proofs


def run():
    models,proofs=geometry(); counts=Counter(); hits=[]; unknown=[]; nonfinite=0
    def one(stem,where,g,w,row,index=None):
        nonlocal nonfinite
        counts[stem]+=1
        out=role(g,w,row['marks'],row.get('corner_nodes'))
        if not out:
            return
        if out['nonfinite']:
            nonfinite+=1; return
        matched=[i for i,p in models.items() if W.mat(g)==W.mat(p['gram'])
                 and w is not None and sorted(w)==sorted(p['white'])
                 and sorted(row['marks'])==sorted(p['marks'])
                 and out['corner_nodes']==sorted(set(p['white'])|set(p['marks']))]
        record=dict(source=stem,path=where,index=index,marks=row['marks'],**out)
        if len(matched)==1:
            record['selected_core']=matched[0]
        else:
            unknown.append(record)
        hits.append(record)
    legacy=0
    for stem,keys in W.finite_inventory().items():
        doc=W.load(stem)
        for key in keys:
            for where,row,g,w,index in W.physical_rows(doc[key],key,require_components=True):
                one(stem,where,g,w,row,index); legacy+=1
        del doc
    assert legacy==564
    common=W.native(1142)[1]['gram']; literal=literal_hits=0
    for stem,index,g,w,cs in W.literal_rows(common):
        before=len(hits)
        for j,row in enumerate(cs):
            one(stem,'literal/'+str(index)+'/'+str(j),g,w,row,index)
        literal+=1; literal_hits+=len(hits)-before
    assert literal==216 and literal_hits==0
    sections=[]
    for stem,keys in W.NEWER:
        doc=W.load(stem); section=doc
        for key in keys:
            section=section[key]
        before=counts[stem]
        for where,row,g,w,index in W.physical_rows(section,'/'.join(keys)):
            one(stem,where,g,w,row,index)
        assert counts[stem]>before, ('missing physical collection',stem,keys)
        sections.append(dict(source=stem,section=list(keys),corners=counts[stem]-before))
        del doc,section
    chosen=[]
    def outputs(x):
        if isinstance(x,dict):
            if 'payment' in x and isinstance(x.get('roles'),list):
                yield x
            for v in x.values():
                yield from outputs(v)
        elif isinstance(x,list):
            for v in x:
                yield from outputs(v)
    for stem,expected in zip(W.CHOSEN,(280,52)):
        rs=[r for r in outputs(W.load(stem)) if F(r['payment'])>0]
        assert len(rs)==expected
        assert all(not any(b.get('type')=='D17' for b in r['roles']) for r in rs)
        chosen.append(dict(source=stem,positive_outputs=expected,D17_outputs=0))
    for stem in W.EXTENSIONS:
        doc=W.load(stem)
        for where,row,g,w,index in W.physical_rows(doc):
            one(stem,where,g,w,row,index)
        del doc
    if unknown:
        raise AssertionError(json.dumps(dict(status='UNRESOLVED_CURRENT_D17_NODE2_ROLE',hits=unknown[:12])))
    assert hits and {h['selected_core'] for h in hits}=={1062,2525}
    files=sorted(set(W.input_paths()+[W.path(SELECTED)]))
    sources={str(p.relative_to(ROOT)):W.digest(p) for p in files}
    for p in (Path(__file__),Path(W.__file__),ROOT/'scripts/fresh_d12_half_face_composition_20261018.py',
              ROOT/'results/fresh_d17_node2_reverse_wall_capacity_research_20261019.md'):
        sources[str(p.relative_to(ROOT))]=W.digest(p)
    return dict(status='PASS_D17_NODE2_SAVED_ROLE_CAPACITY_INTERFACE',sources=sources,
        input_files=len(files),input_bytes=sum(p.stat().st_size for p in files),
        core_models=proofs,later_packets=22,later_corners=legacy,literal_configurations=literal,
        literal_target_hits=literal_hits,newer_sections=sections,corner_counts=dict(counts),
        target_hits=hits,chosen_D12=chosen,excluded_nonfinite_displays=nonfinite,
        activated_payment_upper_per_vertex='(5/8)*[d17>0]',
        guaranteed_refund_lower='(5/8)*max(d17-1,0)',new_positive_families=0,
        new_dimension_bound=None,universal_geometry_formalized=False,
        scope='Exact saved physical-role audit and two complete native core identities. Universal chosen-output completeness, arbitrary-white contact/short-path theorem, finite-vertex simplicity, actual M-matrix facet projection and per-face activation are written premises. A saved PASS is not their proof.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args(); assert not(args.write and args.verify)
    target=ROOT/'results'/(TAG+'.json')
    if args.write:
        assert not target.exists(), 'refuse to overwrite evidence'
    W.guard.resources(); W.memory(); signal.alarm(59)
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        data=json.loads(json.dumps(run()))
        W.guard.resources(start=False); peak=W.memory()
        if args.write:
            target.write_text(json.dumps(data,indent=2)+'\n')
        if args.verify:
            assert data==json.loads(target.read_text())
    signal.alarm(0)
    print(json.dumps(dict(status=data['status'],input_files=data['input_files'],
        later_corners=data['later_corners'],literal=data['literal_configurations'],
        target_hits=len(data['target_hits']),peak_rss_bytes=peak)),flush=True)


if __name__=='__main__':
    main()
