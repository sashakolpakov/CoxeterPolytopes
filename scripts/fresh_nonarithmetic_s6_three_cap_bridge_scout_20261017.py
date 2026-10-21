#!/usr/bin/env python3
"""Numerical three-cap closure of the8/30/44/52 S6 patch.

A and B separately terminate30/44 and30/52; C terminates8/44 and8/52.
Finite contacts A-C and B-C determine the two endpoint parameters, then
A-B and all four new mutual corners are tested. All arithmetic is discovery.
"""
import argparse
import json
import math
from collections import Counter
from hashlib import sha256
from pathlib import Path

import numpy as np
import fresh_nonarithmetic_s6_split_bridge_scout_20261017 as F

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_s6_three_cap_bridge_scout_20261017'
CBANK = 'fresh_nonarithmetic_s6_component_b_bridge_scout_20261017'
S = F.S
M = F.METRIC
LABELS = list(range(2, 17))+[20,24,32,64]


def theta(m):
    return 0. if m == 2 else 2*math.cos(math.pi/m)


def pd(normals):
    return float(np.linalg.eigvalsh(normals @ M @ normals.T)[0]) > 2e-7


def run():
    d = F.read(F.TAG)
    oldbank = F.read(CBANK)
    normals = np.array(d['old_normals'])
    h = np.array([1.]+[0.]*9)
    groups = []
    for group in d['families']:
        groups.append([{**r, **{k:np.array(r[k]) for k in ['p','a','b']},
                        'intervals':[[lo,math.inf if hi is None else hi] for lo,hi in r['intervals']]}
                       for r in group])
    cframe = S+[8,44,52]
    nn = normals[cframe]
    inverse = np.linalg.inv(nn @ M @ nn.T)
    kernel = h-(h @ M @ nn.T) @ inverse @ nn
    assert kernel @ M @ kernel > 0
    kernel /= -math.sqrt(kernel @ M @ kernel)
    candidates, decisions = [], Counter()
    grid = np.array([theta(m) for m in range(2,257)])
    for labels in sorted({tuple(r['labels']) for r in oldbank['groups'][0]['survivors']}):
        contacts = np.array([0.]*4+[-theta(m) for m in labels])
        p = contacts @ inverse @ nn
        residual = 2-p @ M @ p
        if residual <= 0:
            decisions['nonpositive_residual'] += 1
            continue
        for sign in [-1,1]:
            c = p+sign*math.sqrt(residual)*kernel
            values = c @ M @ normals.T
            if max(values) > 3e-6:
                decisions['obtuse_C'] += 1
                continue
            finite = values[values > -2+1e-6]
            if max(np.min(abs(finite[:,None]+grid[None,:]),axis=1),default=0.) > 3e-6:
                decisions['C_non_Coxeter_old_contact'] += 1
                continue
            if not all(pd(np.vstack([normals[S+edge],c])) for edge in [[8,44],[8,52]]):
                decisions['C_endpoint_not_finite'] += 1
                continue
            candidates.append(dict(labels=list(labels),branch=sign,normal=c))
    print('THREE_CAP_C_DOMAIN',len(candidates),dict(decisions),flush=True)
    allhits, counts, stages = [], Counter(), []
    for ci, candidate in enumerate(candidates):
        c = candidate['normal']
        events = []
        for gi, side in enumerate([44,52]):
            rows = []
            for fi, family in enumerate(groups[gi]):
                a,b,p = [float(c @ M @ family[k]) for k in ['a','b','p']]
                for m in LABELS:
                    for t in F.roots(a,p+theta(m),b):
                        if not any(t>lo+1e-7*max(1.,lo) and
                                   (math.isinf(hi) or t<hi-1e-7*max(1.,hi))
                                   for lo,hi in family['intervals']):
                            continue
                        normal = F.evaluate(family,t)
                        if not pd(np.vstack([normals[S+[side]],normal,c])):
                            counts['nonfinite_outer_mutual_corner'] += 1
                            continue
                        other = [j for j in range(63) if j not in family['frame']]
                        if max(normal @ M @ normals[other].T)>-2+2e-6:
                            counts['numerical_old_gate_failed'] += 1
                            continue
                        rows.append(dict(family=fi,parameter=float(t),label_to_C=m,normal=normal))
            events.append(rows)
        stage = dict(C=ci,A_events=len(events[0]),B_events=len(events[1]))
        stages.append(stage)
        counts['A_C_events'] += len(events[0])
        counts['B_C_events'] += len(events[1])
        if not all(events):
            continue
        aa,bb = [np.array([r['normal'] for r in rows]) for rows in events]
        values = aa @ M @ bb.T
        valid = np.argwhere((values<3e-6)&(values>-2+1e-6))
        counts['finite_AB_range_pairs'] += len(valid)
        for ia,ib in valid:
            value = float(values[ia,ib])
            m = int(np.argmin(abs(grid+value)))+2
            error = abs(value+theta(m))
            if error>3e-6:
                counts['AB_Coxeter_gap'] += 1
                continue
            a,b = aa[ia],bb[ib]
            if not pd(np.vstack([normals[S+[30]],a,b])):
                counts['nonfinite_inner_old_corner'] += 1
                continue
            if not pd(np.vstack([normals[S],a,b,c])):
                counts['nonfinite_three_cap_corner'] += 1
                continue
            counts['three_cap_hits'] += 1
            if len(allhits)<100:
                allhits.append(dict(C=ci,C_labels=candidate['labels'],C_branch=candidate['branch'],
                                    A={k:v for k,v in events[0][ia].items() if k!='normal'},
                                    B={k:v for k,v in events[1][ib].items() if k!='normal'},
                                    AB_label=m,AB_product=value,AB_error=error,
                                    normals=[a.tolist(),b.tolist(),c.tolist()]))
        print('THREE_CAP_EVENTS',stage,dict(counts),flush=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{F.TAG}.py',f'results/{F.TAG}.json',f'results/{CBANK}.json']
    return dict(status='COMPLETED_NUMERICAL_THREE_CAP_PATCH_SCOUT',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                retained_walls=63,patch=[8,30,44,52],event_labels=LABELS,
                mutual_AB_recognition_through=256,C_domain_count=len(candidates),C_decisions=dict(decisions),
                counts=dict(counts),events=stages,hits=allhits,
                scope='Numerical discovery on bounded endpoint-label grids and a previously sampled C bank. A/B unprescribed old contacts are strictly dotted. Four new mutual corners are numerically positive definite. No exact chamber, full patch closure, all-label obstruction, or finite-volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],'saved hits',len(result['hits']),flush=True)
