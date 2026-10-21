#!/usr/bin/env python3
"""Intersect all actual corner/edge contact necessities on the590 open cases.

Uses the frozen supporting methods after strengthening their necessary
contact domain. This is not a finite7 payment theorem.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_distance7_high_supporting_frontier_20261017 as S

ROOT,H=S.ROOT,S.H
TAG='fresh_distance7_joint_contact_domains_20261017'


def run():
    parent=json.loads((ROOT/f'results/{S.TAG}.json').read_text())
    for p,h in parent['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    pending=[(i,r)for i,r in enumerate(parent['cases'])if r['status']=='unresolved_supporting_extensions']
    assert len(pending)==590
    matrices={}

    @lru_cache(None)
    def frame(key,black):
        core=parent['cores'][key]
        native=[[F(x)for x in row]for row in core['native_graph']]
        g,sc=H.native_to_rational(native)
        assert H.encode(g)==core['gram']and sc==core['gauge']and H.inertia(g)==core['signature']
        zs=H.nullspace(g);iv=None if zs else H.inverse(g)
        assert [[str(x)for x in z]for z in zs]==core['radical']
        assert (H.encode(iv)if iv else None)==core['inverse']
        matrices[key]=(native,g,sc,zs,iv)
        white=sorted(set(range(len(g)))-set(black));edges=[];forbidden=set()
        for b in black:
            ids=sorted(white+[b]);a=H.extract(g,ids)
            assert H.inertia(a)==[len(ids),0,0]
            j=H.inverse(a);assert all(x>=0 for row in j for x in row)
            vals=[sc[i]*j[k][k]for k,i in enumerate(ids)]
            rejected=[i for i,v in zip(ids,vals)if v>=2]
            forbidden.update(rejected)
            edges.append(dict(black=b,nodes=ids,normalized_inverse_diagonals=list(map(str,vals)),forbidden_contact_nodes=rejected))
        return edges,forbidden

    def options(g,sc,zs,iv,corners,allowed,du,dv):
        out=[]
        for c in corners:
            if c['marked_component']is None:continue
            c=dict(c,possible_contacts=sorted(allowed))
            r=S.ordinary(sc,iv,zs,c,du,dv)
            if r:out.append(r)
            if iv and set(sc)=={1}:
                r=S.relevant(iv,c,du,dv)
                if r:out.append(r)
            if iv:
                r=S.isolated_or_leaf(g,sc,iv,c,du,dv)
                if r:out.append(r)
        return sorted(out,key=lambda c:(c.get('component_order_bound',100),len(c['possible_contacts']),c['method'],c['preferred_corner']))

    rows=[];counts=Counter();bounds=Counter()
    for i,old in pending:
        key=old['core_key'];u,v=old['left_marks'],old['right_marks'];black=tuple(u+v)
        edges,forbidden=frame(key,black);native,g,sc,zs,iv=matrices[key]
        corners=[H.inspect_corner(g,sc,list(black),p,q)for p,q in product(u,v)]
        assert corners==old['all_corners']
        allowed=set.intersection(*(set(c['possible_contacts'])for c in corners))
        joint=allowed-forbidden
        du,dv=H.graph_dist(native,u),H.graph_dist(native,v)
        assert min(du[i]for i in v)>S.THRESHOLD
        corner_opts=options(g,sc,zs,iv,corners,allowed,du,dv)
        opts=options(g,sc,zs,iv,corners,joint,du,dv)
        assert not corner_opts or opts
        rec=dict(parent_index=i,base_context=old['base_context'],extra_edges=old['extra_edges'],core_key=key,
            left_marks=u,right_marks=v,corner_intersection=sorted(allowed),actual_edge_certificates=edges,
            edge_removed_from_intersection=sorted(allowed&forbidden),joint_contact_domain=sorted(joint),
            corner_intersection_alone_certified=bool(corner_opts))
        if corner_opts:counts['corner_intersection_alone']+=1
        if allowed!=joint:counts['edge_shrinks_domain']+=1
        if opts:
            chosen=opts[0];rec.update(status='structural_certificate',chosen=chosen)
            counts['method_'+chosen['method']]+=1
            if 'component_order_bound'in chosen:bounds[chosen['component_order_bound']]+=1
            if not corner_opts:counts['additional_after_strict_edges']+=1
        else:rec['status']='unresolved_supporting_extensions'
        counts[rec['status']]+=1;rows.append(rec)
        if len(rows)%100==0:print('JOINT_CONTACTS',len(rows),dict(counts),flush=True)
    assert len(rows)==590 and counts['structural_certificate']+counts['unresolved_supporting_extensions']==590
    paths=[Path(__file__)]+[ROOT/f'{d}/{S.TAG}.{e}'for d,e in [('scripts','py'),('results','json'),('results','md')]]
    paths += [Path(H.__file__),ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    return dict(status='PASS_DISTANCE7_JOINT_ACTUAL_CONTACT_DOMAINS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        parent_contexts=6703,previous_certified_contexts=6113,previous_unresolved_contexts=590,
        shortcut_threshold=S.THRESHOLD,counts=dict(counts),component_bounds={str(k):v for k,v in sorted(bounds.items())},
        total_certified_contexts=6113+counts['structural_certificate'],remaining_unresolved_contexts=counts['unresolved_supporting_extensions'],
        distinct_core_black_frames=frame.cache_info().currsize,cases=rows,
        scope='Every590 previously unresolved case. A further supporting white must satisfy all four actual-corner PSD necessities and all four actual-edge strict-PD necessities simultaneously. Exact normalized edge inverse diagonal>=2 excludes every nonzero contact there. Frozen arbitrary-support residual/radical/protected-feature methods then run at threshold7. No new weights, selected maps, low/disjoint completion, finite7 theorem, dimension bound or example.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],out['component_bounds'],flush=True)
