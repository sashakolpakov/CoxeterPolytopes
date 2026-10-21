#!/usr/bin/env python3
"""Numerical old finite-edge frontier of the actual target H3 face."""
import argparse,hashlib,itertools,json
from pathlib import Path
import numpy as np
from fresh_h20_intrinsic_facet_incidence_20261016 import bank
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158
from fresh_h20_release144_cap_20261016_residual_scout import close

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_open_edge_frontier_20261016'
CORE='results/fresh_h20_intrinsic_facet_incidence_20261016.json'

def run():
    G=geometry158();M=G['metric'];X=G['X'].astype(np.longdouble)
    u=X[151];qu=u@M@u;hu=X@M@u
    B=X@M@X.T-np.outer(hu,hu)/qu;oq=np.diag(B)
    old=np.flatnonzero(oq>1e-7);K=B[np.ix_(old,old)]/np.sqrt(np.outer(oq[old],oq[old]))
    rows,_=bank();R=np.array([r['root']for r in rows],dtype=np.longdouble)
    labs=np.array([r['label151']for r in rows]);h=-2*np.cos(np.pi/labs);h[labs==2]=0;rq=2-h*h/2
    C=(R@M@X[old].T-h[:,None]*hu[old]/qu)/np.sqrt(rq[:,None]*oq[old])
    core=json.loads((ROOT/CORE).read_text())['surviving_degree_core'];core_set=set(core)
    edges=[];vertices=set()
    for a,b in itertools.combinations(range(len(old)),2):
        kab=K[a,b]
        if 1-kab*kab<=1e-6:continue
        existing=[]
        for c in range(len(old)):
            if c in (a,b):continue
            kac,kbc=K[a,c],K[b,c]
            det=1+2*kab*kac*kbc-kab*kab-kac*kac-kbc*kbc
            if min(1-kac*kac,1-kbc*kbc,det)>1e-6:
                existing.append(int(old[c]));vertices.add(tuple(sorted(map(int,[old[a],old[b],old[c]]))))
        ra,rb=C[:,a],C[:,b]
        det=1+2*kab*ra*rb-kab*kab-ra*ra-rb*rb
        candidates=np.flatnonzero((1-ra*ra>1e-6)&(1-rb*rb>1e-6)&(det>1e-6)).tolist()
        retained=[i for i in candidates if i in core_set]
        assert len(existing)<=2
        # Adding an acute plane cannot subdivide an already closed edge.
        if len(existing)==2:assert not candidates
        edges.append(dict(edge=[int(old[a]),int(old[b])],old_third_facets=existing,
                          finite_endpoints_needed=2-len(existing),all_bank_candidate_count=len(candidates),
                          degree_core_candidate_ids=retained,
                          degree_core_candidates=[dict(bank_id=i,**{k:v for k,v in rows[i].items()if k!='root'})for i in retained]))
    assert len(edges)==13 and len(vertices)==7
    assert sum(e['finite_endpoints_needed']for e in edges)==5
    src=[f'scripts/{TAG}.py',CORE,'scripts/fresh_h20_intrinsic_facet_incidence_20261016.py',
         'scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py',
         'results/fresh_h20_compact_quad_cap_20261016.json']
    return dict(status='PASS_NUMERICAL_OLD_H3_FINITE_EDGE_FRONTIER',
                sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},
                actual_old_facets=old.tolist(),old_finite_vertices=[list(x)for x in sorted(vertices)],
                edges=edges,open_edges=[e['edge']for e in edges if e['finite_endpoints_needed']],
                scope='Numerical positive-Gram finite-incidence census on seven old facets of H3 S16+151. Existing seven finite vertices and thirteen finite edges; five edges need a second finite endpoint for compact closure. Candidate inventories are finite banks only, and no limitation is asserted for the new continuous edge-driven domain.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['edges']),d['open_edges'],flush=True)
