#!/usr/bin/env python3
"""Numerical discovery of seed metrics supporting three prescribed cap families."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import root

import fresh_h9_multiblock_seed_20261013 as base
import fresh_h9_three_parameter_seed_20261013 as helper

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'data/enumerations/fresh_h9_three_parameter_seed_20261013_scout.json'
OUT = ROOT/'data/enumerations/fresh_h9_three_norm_seed_20261013_scout.json'


def run(output, trials, rng_seed):
    start_time = time.monotonic(); rng = np.random.default_rng(rng_seed)
    base.VALS = np.array([0., 1., np.sqrt(2), base.PHI])
    seeds = json.loads(SOURCE.read_text())['records']
    tables = {}
    for ns in set(tuple(s['groups']) for s in seeds):
        P = base.profiles(ns); tv, dv, sq = P[2:5]
        U = -.5*np.array([helper.symprod(t,t) for t in tv])
        V = -np.array([helper.symprod(d,d) for d in dv])/32
        C = sq/8-2
        _, keep, inverse = np.unique(np.round(np.c_[U,V,C],10),axis=0,return_index=True,return_inverse=True)
        tables[ns] = (U[keep],V[keep],C[keep],keep,inverse)
    points = {}; seen = set(); attempts = 0; converged = 0; rank_deficient = 0
    for trial in range(trials):
        row = seeds[int(rng.integers(len(seeds)))]; ns = tuple(row['groups'])
        Uall, Vall, Call, original, inverse = tables[ns]
        available = np.unique(inverse[row['profile_indices']])
        if len(available) < 2:
            continue
        ij = rng.choice(available,2,replace=False); k = int(rng.integers(len(Call)))
        ids = tuple(sorted([int(ij[0]),int(ij[1]),k]))
        if len(set(ids)) < 3:
            continue
        key = (ns,ids)
        if key in seen:
            continue
        seen.add(key); U = Uall[list(ids)]; V = Vall[list(ids)]; C = Call[list(ids)]
        r0 = np.array(row['inverse_metric'])
        def fun(r):
            try:
                b, db = helper.matrices(r,ns)
                return U@r+V@b+C, U+V@db
            except ValueError:
                return np.full(3,1e12), np.eye(3)*1e12
        if np.linalg.matrix_rank(fun(r0)[1],tol=1e-9) < 3:
            rank_deficient += 1
            continue
        attempts += 1
        sol = root(fun,r0,jac=True,method='hybr',options={'xtol':1e-10,'maxfev':150})
        residual = max(abs(fun(sol.x)[0])); M = helper.metric(sol.x,ns)
        if residual > 1e-8 or M is None:
            continue
        converged += 1; mk = (ns,tuple(np.round(M.ravel(),7)))
        points.setdefault(mk,{'groups':list(ns),'metric':M.tolist(),'inverse_metric':sol.x.tolist(),
                               'anchor_profiles':original[list(ids)].tolist(),
                               'max_equation_residual':float(residual)})
        if attempts%5000 == 0:
            print({'trials':trial+1,'root_attempts':attempts,'converged':converged,
                   'distinct_metrics':len(points),'seconds':round(time.monotonic()-start_time,1)},flush=True)
    result = {'status':'NUMERICAL_THREE_NORM_METRIC_PROPOSALS','source':str(SOURCE.relative_to(ROOT)),
              'trials':trials,'rng_seed':rng_seed,'distinct_triples':len(seen),
              'rank_deficient_starts_skipped':rank_deficient,'root_attempts':attempts,
              'valid_solutions':converged,'records':list(points.values()),
              'scope':'Random triples of distinct prescribed cap-norm equations: two supported by one of the saved 4671 metrics and one from its complete four-label profile table. One starting point retained per triple. Singular initial Jacobians skipped. Numerical roots only, with metric entries in (2,60] and weighted negative eigenvalue in (-4,0). No exhaustive parameter or existence conclusion.'}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print({'done':len(points),'attempts':attempts,'seconds':round(time.monotonic()-start_time,1)},flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=OUT)
    p.add_argument('--trials',type=int,default=200000);p.add_argument('--seed',type=int,default=103613)
    p.add_argument('--screen',action='store_true');args=p.parse_args()
    if args.screen:helper.screen(args.output)
    else:run(args.output,args.trials,args.seed)
