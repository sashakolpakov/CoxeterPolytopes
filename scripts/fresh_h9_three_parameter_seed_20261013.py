#!/usr/bin/env python3
"""Discovery: enforce two cap norms and their angle while changing seed metric.

All conclusions here are numerical and restricted to the saved proposal set.
An exact verifier is required for any positive Coxeter construction claim.
"""
import argparse
import itertools
import json
from pathlib import Path
import time

import networkx as nx
import numpy as np
from scipy.optimize import root

import fresh_h9_multiblock_seed_20261013 as base
from fresh_h9_asymmetric_blocks_20261013 import reduced_geometry

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/enumerations/fresh_h9_asymmetric_four_labels_20261013_scout.json'
OUT = ROOT / 'data/enumerations/fresh_h9_three_parameter_seed_20261013_scout.json'


def symprod(x, y):
    return np.array([x[0]*y[0], x[0]*y[1]+x[1]*y[0], x[1]*y[1]])


def matrices(r, ns):
    u, s, v = r + np.array([ns[0]/4, 0, ns[1]/4])
    det = u*v-s*s
    if abs(det) < 1e-15:
        raise ValueError('singular positive block')
    b = np.array([v, -s, u])/det
    p, q, t = b
    db = -np.array([[p*p, 2*p*q, q*q],
                    [p*q, p*t+q*q, q*t],
                    [q*q, 2*q*t, t*t]])
    return b, db


def metric(r, ns):
    p, q, t = r
    det = p*t-q*q
    if abs(det) < 1e-12:
        return None
    M = np.array([[t, -q], [-q, p]])/det
    if np.min(M) <= 2+1e-6 or np.max(M) > 60:
        return None
    D = np.diag(np.sqrt(ns))
    vals = np.linalg.eigvalsh(D@M@D)
    if not (-4+1e-6 < vals[0] < -1e-6 and vals[1] > 1e-6):
        return None
    return M


def seed(ns, M):
    Z = np.zeros((8, 2)); start = 0
    for g, n in enumerate(ns):
        Z[start:start+n, g] = 1; start += n
    K = Z@M@Z.T
    return np.block([[2*np.eye(8), -2*np.eye(8)-K],
                     [-2*np.eye(8)-K, 2*np.eye(8)]]), Z


def complete_bank(ns, M, P):
    _, _, tv, dv, sq, _ = P
    inv = np.linalg.inv(M)
    R = np.array([inv[0, 0], inv[0, 1], inv[1, 1]])
    B, _ = matrices(R, ns)
    TT = np.array([symprod(t, t) for t in tv])
    DD = np.array([symprod(d, d) for d in dv])
    norms = -.5*TT@R + sq/8 - DD@B/32
    ii = np.flatnonzero(abs(norms-2) < 2e-7)
    G, _ = seed(ns, M); inverse = np.linalg.pinv(G, rcond=1e-10)
    types, inds, *_ = P; covs = []
    for i in ii:
        ts = [types[g][j] for g, j in enumerate(inds[i])]
        for parts in itertools.product(*[list(base.variants(n, t)) for n, t in zip(ns, ts)]):
            covs.append(-np.r_[np.concatenate([p[0] for p in parts]),
                              np.concatenate([p[1] for p in parts])])
    covs = np.unique(np.round(covs, 13), axis=0)
    Rcap = covs@inverse
    assert np.max(abs(Rcap@G-covs)) < 3e-6
    assert np.max(abs(np.sum(Rcap*covs, axis=1)-2)) < 3e-6
    return G, Rcap, covs, ii


def proposals(rows):
    """All relative placements, fixing first orbit member by seed symmetry."""
    cases = {}
    for source_index, row in rows:
        ns = tuple(row['groups']); P = base.profiles(ns)
        G, caps, refs = base.bank(ns, row['a'], row['b'], P, row['profile_indices'])
        covs = np.unique(np.round(-np.c_[np.array(refs)[:, 0], np.array(refs)[:, 1]], 13), axis=0)
        features = []
        for cov in covs:
            x, y = -cov[:8], -cov[8:]; d = x-y
            t = np.array([x[0]+y[0], x[ns[0]]+y[ns[0]]])
            dv = np.array([sum(d[:ns[0]]), sum(d[ns[0]:])])
            features.append((t, dv, d, float(d@d)))
        representatives = {}
        for i, (t, dv, d, sq) in enumerate(features):
            key = tuple(np.round(np.r_[t, dv, sq], 9))
            representatives.setdefault(key, i)
        M = np.full((2, 2), row['a']); M -= row['b']*np.eye(2)
        inv = np.linalg.inv(M); start = np.array([inv[0, 0], inv[0, 1], inv[1, 1]])
        for i in representatives.values():
            t1, d1, v1, q1 = features[i]
            for j, (t2, d2, v2, q2) in enumerate(features):
                U = -.5*np.array([symprod(t1, t1), symprod(t2, t2), symprod(t1, t2)])
                V = -np.array([symprod(d1, d1), symprod(d2, d2), symprod(d1, d2)])/32
                c = np.array([q1/8-2, q2/8-2, v1@v2/8])
                if np.linalg.norm(np.r_[U[0]-U[1], V[0]-V[1], c[0]-c[1]]) < 1e-8:
                    continue
                key = (ns, tuple(np.round(np.r_[U.ravel(), V.ravel(), c], 9)))
                if key not in cases:
                    cases[key] = (ns, U, V, c, start, source_index, [covs[i].tolist(), covs[j].tolist()])
    return list(cases.values())


def run(output, max_systems, seed_value, labels):
    started = time.monotonic()
    base.VALS = np.array([0., 1., np.sqrt(2), base.PHI])
    data = json.loads(SOURCE.read_text())
    rows = [(i, r) for i, r in enumerate(data['records']) if r['geometry_all_caps']['outside'] == 0]
    cases = proposals(rows)
    rng = np.random.default_rng(seed_value); rng.shuffle(cases)
    print({'proposal_systems':len(cases), 'running':min(len(cases), max_systems)}, flush=True)
    points = {}; attempts = 0; converged = 0; signatures = set()
    for index, (ns, U, V, c, start, source_index, covs) in enumerate(cases[:max_systems]):
        for m in labels:
            target = 0. if m == 2 else 2*np.cos(np.pi/m)
            C = c + np.array([0., 0., target])
            def fun(r):
                try:
                    b, db = matrices(r, ns)
                    return U@r + V@b + C, U + V@db
                except ValueError:
                    return np.full(3, 1e12), np.eye(3)*1e12
            sol = root(fun, start, jac=True, method='hybr', options={'xtol':1e-10, 'maxfev':120})
            attempts += 1
            if np.max(abs(fun(sol.x)[0])) > 1e-8:
                continue
            M = metric(sol.x, ns)
            if M is None:
                continue
            converged += 1
            key = (ns, tuple(np.round(M.ravel(), 7)))
            if key in signatures:
                continue
            signatures.add(key)
            points[key] = {'groups':list(ns), 'metric':M.tolist(), 'inverse_metric':sol.x.tolist(),
                           'source_record':source_index, 'anchor_label':m,
                           'anchor_covectors':covs, 'max_equation_residual':float(max(abs(fun(sol.x)[0])))}
        if (index+1) % 200 == 0:
            print({'systems':index+1, 'attempts':attempts, 'valid_solutions':converged,
                   'distinct_metrics':len(points), 'seconds':round(time.monotonic()-started, 1)}, flush=True)
    result = {'status':'NUMERICAL_THREE_PARAMETER_SEED_PROPOSALS', 'source':str(SOURCE.relative_to(ROOT)),
              'proposal_system_count':len(cases), 'processed_systems':min(max_systems,len(cases)),
              'anchor_labels':labels, 'rng_seed':seed_value, 'attempts':attempts,
              'valid_solutions':converged, 'records':list(points.values()),
              'scope':'Numerical roots from one inherited start per distinct prescribed pair system, drawn from the 54 previously all-bank-bounded four-label configurations. Two norm constraints and one pair angle are imposed simultaneously in a symmetric 2x2 core with three independent entries. Entries lie in (2,60], with negative weighted-core eigenvalue in (-4,0). No exhaustion of roots, parameter space, or finite-volume claim.'}
    output.write_text(json.dumps(result, indent=2)+'\n')
    print({'done':len(points), 'seconds':round(time.monotonic()-started, 1)}, flush=True)


def screen(source):
    data = json.loads(source.read_text())
    base.VALS = np.array(data.get('pairing_magnitudes',[0., 1., np.sqrt(2), base.PHI]))
    records = data['records']; start = time.monotonic()
    for i, row in enumerate(records):
        if 'geometry_all_caps' in row:
            continue
        ns = row['groups']; M = np.array(row['metric']); P = base.profiles(ns)
        G, R, covs, ii = complete_bank(ns, M, P)
        geo = reduced_geometry(ns, G, R); geo.pop('N')
        row.update({'cap_count':len(R), 'profile_indices':ii.tolist(), 'geometry_all_caps':geo})
        if geo['outside'] == 0:
            C = base.compat(R@G@R.T)
            cliques = list(nx.find_cliques(nx.from_numpy_array(C)))
            best = None
            for clique in cliques:
                cgeo = base.geometry(G, R, clique)
                if best is None or cgeo['max_klein_square'] < best['geometry']['max_klein_square']:
                    best = {'cap_indices':clique, 'covectors':covs[clique].tolist(), 'geometry':cgeo}
            row.update({'compatible_pair_count':int(C.sum()//2), 'maximal_cliques':len(cliques),
                        'maximum_clique_size':max(map(len,cliques)), 'best_selection':best})
            print({'record':i, 'all_bank_bounded':True, 'caps':len(R),
                   'cliques':len(cliques), 'best':best['geometry']['max_klein_square']}, flush=True)
        if (i+1) % 100 == 0:
            source.write_text(json.dumps(data,indent=2)+'\n')
            print({'screened':i+1, 'total':len(records), 'seconds':round(time.monotonic()-start,1)},flush=True)
    data['status'] = 'NUMERICAL_THREE_PARAMETER_SEED_SCREEN_COMPLETE'
    source.write_text(json.dumps(data,indent=2)+'\n')
    bounded = [r for r in records if r['geometry_all_caps']['outside'] == 0]
    print({'done':len(records), 'all_bank_bounded':len(bounded),
           'selected_bounded':sum(r['best_selection']['geometry']['outside']==0 for r in bounded)},flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--output',type=Path,default=OUT)
    p.add_argument('--max-systems',type=int,default=20000); p.add_argument('--seed',type=int,default=613103)
    p.add_argument('--labels',type=int,nargs='+',default=[2,3,4,5,6]); p.add_argument('--screen',action='store_true')
    args = p.parse_args()
    if args.screen: screen(args.output)
    else: run(args.output,args.max_systems,args.seed,args.labels)
