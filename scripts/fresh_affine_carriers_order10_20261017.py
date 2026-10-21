#!/usr/bin/env python3
"""Check the genuinely new order10--17 affine-carrier cutoff ingredients.

The all-rank proof is in the accompanying report.  Frozen high-rank role
enumerations are bound as data, not rerun.  Literal small-rank probes are
corroboration of explicit leaf and classical-extension arguments.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import fresh_private_affine_tree_exclusion_20261017 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_affine_carriers_order10_20261017'
SOURCES=[
    'scripts/fresh_private_affine_tree_exclusion_20261017.py',
    'results/fresh_private_affine_tree_exclusion_20261017.md',
    'scripts/fresh_affine_tree_carrier_exclusion_20261017.py',
    'results/fresh_affine_tree_carrier_exclusion_20261017.json',
    'results/fresh_affine_tree_carrier_exclusion_20261017.md',
    'results/fresh_private_affine_cycle_exclusion_20261017.md',
    'results/fresh_affine_cycle_carrier_exclusion_20261017.md',
    'results/fresh_simple_cusp_vinberg_weight_frontier_20261017.md',
]

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def build():
    keys,catalogue=P.affine_catalogue();rows=[];counts=Counter()
    for n in range(10,18):
        for kind in 'BCD':
            es,nu=P.graph(kind,n)
            witnesses=[P.leaf_witness(n,es,i,label,keys)
                       for i in range(n)for label in (3,4)]
            for w in witnesses:counts[w['affine_type']]+=1
            rows.append(dict(kind=kind,order=n,witnesses=witnesses))
    # New all-rank dichotomy: available simple arms L,R have L+R>=7.
    # If both are >=3 use affine E7; otherwise min=2 and max>=5 gives E8.
    arm_checks=[]
    for L in range(2,8):
        for R in range(2,8):
            if L+R<7:continue
            witness='E7'if min(L,R)>=3 else'E8'
            assert (L>=3 and R>=3)or(min(L,R)==2 and max(L,R)>=5)
            arm_checks.append([L,R,witness])
    # Above7 an arm can only make these fixed witnesses easier to embed.
    # The old resistance role range is independent of ambient tree order.
    role=json.loads((ROOT/'results/fresh_affine_tree_carrier_exclusion_20261017.json').read_text())
    assert role['role_census']['trials']==908 and len(role['role_census']['rows'])==11
    assert F(9,4)>2
    # New boundary-order continuation checks; no old18/19/23 replay.
    continuations=[]
    for kind in 'BC':
        for n in range(10,18):
            es,nu=P.graph(kind,n);adj=P.adjacency(n,es)
            x=next(i for i,a in enumerate(nu)if a==F(1,2))
            p=adj[x][0][0];y=next(i for i,_ in adj[p]if i!=x)
            z,w=n,n+1;ee=es+[(x,z,3),(y,z,3),(z,w,3)]
            keep=[i for i in range(n+2)if i!=p]
            target='B'if kind=='C'else'D';te,_=P.graph(target,n+1)
            assert P.tree_key(n+1,P.induced(ee,keep))==P.tree_key(n+1,te)
            continuations.append(dict(kind=kind,order=n,deleted=p,result=target,result_order=n+1))
    cycles=[]
    for n in range(10,18):
        es,_=P.graph('A',n)
        for label in(3,4):cycles.append(dict(order=n,label=label,witness=P.leaf_witness(n,es,0,label,keys)))
    return dict(status='PASS_AFFINE_CARRIER_ORDER10_CUTOFF_INGREDIENTS',
        sources={p:sha(p)for p in SOURCES},source_sha256=sha('scripts/'+TAG+'.py'),
        minimum_affine_order=10,forbidden_affine_order_through=9,
        all_rank_new_bounds=dict(simple_arm_sum_lower_bound=7,opposite_end_resistance_lower_bound=9,
                                opposite_end_cost_lower_bound='9/4',finite_deletion_order_lower_bound=9),
        affine_witness_catalogue=catalogue,leaf_probes=rows,leaf_type_counts=dict(sorted(counts.items())),
        finite_arm_dichotomy_checks=arm_checks,continuations=continuations,cycle_leaf_probes=cycles,
        inherited_role_trials=908,inherited_surviving_roles=11,
        counts=dict(tree_leaf_witnesses=sum(len(x['witnesses'])for x in rows),cycle_leaf_witnesses=len(cycles),
                    boundary_continuations=len(continuations)),
        scope='New cutoff ingredients and literal10--17 corroboration; the report gives the all-rank private/carrier transfer. No ideal-quadrilateral payment or ambient dimension conclusion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=json.loads(json.dumps(build()));out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],d['counts'],flush=True)
