#!/usr/bin/env python3
"""Numerical necessary facet-incidence screen for a compact target H3 face.

All candidate walls are right to S16. Incidence is computed after projecting
off 151, and uses projected old norms, not ambient native norms. The screen
is a finite floating-point diagnostic, never a continuous-domain exclusion.
"""
import argparse
import hashlib
import json
from collections import Counter, deque
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158
from fresh_h20_release144_cap_20261016_residual_scout import close

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_intrinsic_facet_incidence_20261016'
BASE = 'results/fresh_h20_two_anchor_families_20261016.json'
PAIRS = 'results/fresh_h20_two_anchor_families_20261016_pairs.json'
DISCRETE = 'results/fresh_h20_nonright151_h3_closure_20261016.json'
COORDINATED = 'results/fresh_h20_two_anchor_families_20261016_coordinated_angles.json'
COVER = 'results/fresh_h20_continuous_sample_joint_cover_20261016.json'


def bank():
    d = json.loads((ROOT / BASE).read_text())
    p = json.loads((ROOT / PAIRS).read_text())
    e = json.loads((ROOT / DISCRETE).read_text())
    c = json.loads((ROOT / COORDINATED).read_text())
    rows = []
    for r in d['candidates']:
        f = d['families'][r['family']]
        rows.append(dict(origin='continuous_sample', origin_id=r['id'],
                         family=r['family'], parameter=r['parameter'],
                         label151=f['labels'][0], root=r['orthonormal_chart']))
    for i, r in enumerate(p['strong_pairs']):
        rows.append(dict(origin='solved_partner', origin_id=i, family=r['family'],
                         parameter=r['parameter'], label151=d['families'][r['family']]['labels'][0],
                         root=r['partner_orthonormal_chart']))
    for r in e['candidates']:
        rows.append(dict(origin='three_contact_discrete', origin_id=r['id'],
                         label151=r['labels'][0], root=r['orthonormal_chart']))
    groups = []
    for i, r in enumerate(c['successes']):
        group = []
        for j, (f, t, x) in enumerate(zip(r['family_ids'], r['parameters'], r['roots'])):
            group.append(len(rows))
            rows.append(dict(origin='coordinated_cluster', origin_id=[i,j], family=f,
                             parameter=t, label151=d['families'][f]['labels'][0], root=x))
        groups.append(group)
    return rows, groups


def histogram(values):
    return {str(k):v for k,v in sorted(Counter(map(int,values)).items())}


def run():
    G = geometry158()
    M = G['metric']
    X = G['X'].astype(np.longdouble)
    metric = M.astype(np.longdouble)
    u = X[151]
    qu = u @ metric @ u
    hu = X @ metric @ u
    oq = np.einsum('ij,jk,ik->i',X,metric,X)-hu*hu/qu
    old = np.where(oq > 1e-7)[0]
    assert old.tolist() == [95,142,145,150,152,156,157]
    rows, groups = bank()
    R = np.array([r['root'] for r in rows])
    labs = np.array([r['label151'] for r in rows])
    # These contacts and norms were prescribed analytically in each family.
    h = -2*np.cos(np.pi/labs); h[labs==2] = 0
    qr = 2-h*h/2
    raw_old = R.astype(np.longdouble) @ metric @ X[old].T
    proj_old = raw_old-h[:,None]*hu[old]/qu
    det = qr[:,None]*oq[old]-proj_old*proj_old
    old_edges = det > 1e-6
    old_degree = old_edges.sum(axis=1)
    n = len(rows)
    edge_i, edge_j = [], []
    # A finite intrinsic edge requires both the ambient Coxeter gate and a
    # positive definite 2-by-2 projected Gram. No dotted pair can qualify.
    for start in range(0,n,256):
        stop = min(start+256,n)
        co = -(R[start:stop] @ M @ R.T)/2
        intrinsic = qr[start:stop,None]*qr[None,:]-(2*co+h[start:stop,None]*h[None,:]/2)**2
        eligible = (co > -2e-6)&(co < 1-2e-6)&(intrinsic > 1e-6)
        ii,jj = np.where(eligible)
        greater = jj > ii+start
        ii,jj = ii[greater],jj[greater]
        vals = co[ii,jj]
        angle = np.pi/np.arccos(np.clip(vals,-1,1))
        good = (abs(vals)<2e-6)|((angle>=2.9999)&(abs(angle-np.rint(angle))<2e-4))
        ii,jj = ii[good]+start,jj[good]
        edge_i.extend(ii.tolist());edge_j.extend(jj.tolist())
        if start%4096==0:
            print('INCIDENCE',start,n,'edges',len(edge_i),flush=True)
    ii=np.array(edge_i,dtype=np.int32);jj=np.array(edge_j,dtype=np.int32)
    adjacency=coo_matrix((np.ones(2*len(ii),dtype=np.int8),
                        (np.r_[ii,jj],np.r_[jj,ii])),shape=(n,n)).tocsr()
    degree=np.diff(adjacency.indptr).astype(int)
    alive=np.ones(n,dtype=bool)
    queue=deque(np.where(old_degree+degree<3)[0].tolist())
    removed=[]
    while queue:
        v=queue.popleft()
        if not alive[v]:continue
        alive[v]=False;removed.append(v)
        for w in adjacency.indices[adjacency.indptr[v]:adjacency.indptr[v+1]]:
            if alive[w]:
                degree[w]-=1
                if old_degree[w]+degree[w]<3:queue.append(int(w))
    cover=json.loads((ROOT/COVER).read_text())
    cluster_rows=[]
    def summarize(selected,origin):
        for v in selected:
            neighbors=adjacency.indices[adjacency.indptr[v]:adjacency.indptr[v+1]]
            new=sorted(set(selected)&set(map(int,neighbors)))
            cluster_rows.append(dict(cluster=origin,root=v,old_neighbors=old[old_edges[v]].tolist(),
                                     new_neighbors=new,total_degree=int(old_degree[v])+len(new)))
    for step in cover['steps']:
        if 'selected' in step:summarize(step['selected'],['sample_cover',step['iteration']])
    for i,group in enumerate(groups):summarize(group,['coordinated',i])
    src=[f'scripts/{TAG}.py',BASE,PAIRS,DISCRETE,COORDINATED,COVER,
         'scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py',
         'scripts/fresh_h20_release144_cap_20261016_residual_scout.py']
    return dict(status='PASS_NUMERICAL_INTRINSIC_FACET_INCIDENCE_SCREEN',
                sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},
                candidate_count=n,origins=dict(Counter(r['origin']for r in rows)),
                positive_old_projected_norms={str(i):float(oq[i])for i in old},
                old_degree_histogram=histogram(old_degree),new_edge_count=len(ii),
                initial_total_degree_histogram=histogram(old_degree+np.diff(adjacency.indptr)),
                surviving_degree_core=np.where(alive)[0].tolist(),removed_count=len(removed),
                survivor_origins=dict(Counter(rows[i]['origin']for i in np.where(alive)[0])),
                cluster_incidence=cluster_rows,
                cluster_low_degree_counts=dict(Counter(str(r['cluster'])for r in cluster_rows if r['total_degree']<3)),
                scope='Necessary minimum facet-degree3 for a compact H3 completion, finite numerical union of two-anchor samples, solved partners, three-contact roots and coordinated clusters. Uses projected old norms and all pairwise finite Coxeter incidences, followed by degree-core peeling. No exact exclusion or continuous-domain coverage; no global H20 finite-volume conclusion.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['candidate_count'],'edges',d['new_edge_count'],
          'core',len(d['surviving_degree_core']),flush=True)
