#!/usr/bin/env python3
"""SOURCE ONLY: current independent affineD13 centre-zero output theorem.

All17 finite/ideal/short4 alternatives are bound. The complete all-native
profile and energy-one classification remain explicit accepted premises.
The source never invokes S8 or calls the zero corner an intrinsic unit.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import types

ROOT=Path(__file__).resolve().parents[1]
LOADER='scripts/check_finite57_current_adaptive9_path_binding_20261025.py'
LOADER_SHA='76022d95fcc5ec6e61bf1ccaea65635adda1ab8593d0bf921d75754f14bd139c'
HALF='fresh_order13_d_center_pair_half_transfer_20261018'
ZERO='fresh_order13_d_center_pair_zero_transfer_20261018'
EXTRA='results/finite57_current_d13_center_zero_extra_inputs_20261025.json'
EXTRA_SHA='4bf29435f19fac4c352e3fd183cd2e484cf92947e45c74cab4001d71c8b7d1f3'
MAX_OUTPUT=500*1024


def need(ok,message):
    if not ok:raise ValueError(message)


def raw(name):
    p=(ROOT/name).resolve()
    need(p.is_relative_to(ROOT) and p.is_file(),'bounded repository source')
    with p.open('rb') as f:b=f.read(4*1024**2+1)
    need(0<len(b)<=4*1024**2,'bounded source bytes')
    return b


def prepare():
    encoded=raw(LOADER);need(sha256(encoded).hexdigest()==LOADER_SHA,'exact shared current loader')
    a=types.ModuleType('finite57_current_D13zero_loader');a.__file__=str(ROOT/LOADER)
    exec(compile(encoded,a.__file__,'exec'),a.__dict__)
    r,c,first_extra=a.prepare()
    extra=raw(EXTRA);need(sha256(extra).hexdigest()==EXTRA_SHA,'exact D13zero extra inputs')
    entries=json.loads(extra,object_pairs_hook=c.unique)['inputs']
    expected={f'{folder}/{stem}.{ext}' for stem in (HALF,ZERO)
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))}
    need(len(entries)==6 and {p['path'] for p in entries}==expected,'exact six D13zero proof/data inputs')
    for p in entries:
        if p['path'] in c.PINS:need(c.PINS[p['path']]==p['sha256'],'extra/base identity collision')
        c.PINS[p['path']]=p['sha256'];need(c.digest(p['path'])==p['sha256'],'D13zero source input identity')
        c.READ.add(p['path'])
    return a,r,c,first_extra,extra


def domain(rows,key,expected,message):
    need(type(rows) is list and len(rows)==len(expected),message)
    actual=[key(p) for p in rows]
    need(len(actual)==len(set(actual)) and set(actual)==set(expected),message)


def short4(c):
    counts={'finite':0,'affine':0};digest=sha256()
    for ideal,start,kinds in ((False,12,('A','B','D')),(True,10,('A','B','C','D'))):
        for kind in kinds:
            for order in range(start,59):
                g=c.W.classical_graph(kind,order,affine=ideal)
                ds=[c.W._distances(g,i) for i in range(order)]
                for p,q in combinations(range(order),2):
                    c.E.charge(c.FUEL)
                    if ds[p][q]>4:continue
                    value=c.W.affine_weight(kind,order,p,q) if ideal else c.W.finite_same_weight((kind,order),p,q)
                    if not ideal and (kind,order)==('D',13) and (p,q) in ((0,11),(0,12)):
                        value=min(value,F(1,2))
                    need(value==1,'current whole short4 pair lost its intrinsic unit')
                    counts['affine' if ideal else 'finite']+=1
                    digest.update(c.canonical([ideal,kind,order,p,q,str(value)]).encode()+b'\n')
    need(all(counts.values()),'complete nonempty short4 domains')
    # Every small finite type is unit through rank11. Exceptional affine
    # factors are unit through order9; no unsupported large exceptional
    # family is suppressed by the classical loops.
    for typ in (('A',1),('A',2),('I',2),('H',3),('H',4),('F',4),('E',6),('E',7),('E',8)):
        for p,q in combinations(range(typ[1]),2):
            need(c.W.finite_same_weight(typ,p,q)==1,'small/exceptional current unit interface')
    return dict(counts=counts,ordered_domain_sha256=digest.hexdigest(),ambient=58,
        finite_small_rule='All whole finite ranks<=11 are unit, including every I2(m) class.',
        affine_small_rule='All whole affine orders<=9 are unit, including exceptional types.',
        geometric_meaning='Any actual adjacent corner containing a white-interior path of length<=4 has whole marked-factor graph distance<=4. No S8 theorem is used.',
        unmarked_positive_witnesses_credited=False,negative_D13_guard_accounted=True)


def base_gram():
    n=13;g=[['2' if i==j else '0' for j in range(n)] for i in range(n)]
    for i,j in [(i,i+1) for i in range(8)]+[(0,9),(0,10),(8,11),(8,12)]:
        g[i][j]=g[j][i]='-1'
    return g


def physical_coupon(side,leaf,extra):
    n=14+extra;g=[['2' if i==j else '0' for j in range(n)] for i in range(n)]
    old=base_gram()
    for i in range(13):
        for j in range(13):g[i][j]=old[i][j]
    g[13][leaf]=g[leaf][13]='-1'
    if extra:g[13][14]=g[14][13]='-1'
    own,far=(0,8) if side==0 else (8,0)
    marks=[far,13];white=[i for i in range(n) if i not in (own,far,13)]
    return g,white,marks


def bind(c):
    formulas=c.field(HALF,('formulas',))
    need(formulas['resistance']==8 and formulas['extreme_cost']==6 and
         formulas['remaining_nonextreme_maximum']==3 and
         formulas['two_nonextreme_strict_gap']=='1/4' and
         formulas['one_extreme_weighted_energy_maximum']==1 and
         formulas['nonzero_extreme_support_energy_lower_bound']==2 and
         formulas['both_extreme_unscaled_energy']==1,
         'exact universal profile/energy-domain propositions used by zero composition')
    immediate=c.field(HALF,('immediate_affine_unit_corners',))
    finite=c.field(HALF,('finite_coupon_corners',))
    mixed=c.field(HALF,('mixed_energy_one_cases',))
    branches=c.field(ZERO,('inherited_branches_with_independent_units',))
    coupon_updates=c.field(ZERO,('upgraded_coupons',))
    need(c.field(HALF,('order',))==c.field(ZERO,('order',))==13 and
         c.field(HALF,('affine_type',))==c.field(ZERO,('affine_type',))=='D12' and
         c.field(ZERO,('changed_pairs',))==[[0,8]] and c.field(ZERO,('new_weight',))=='0' and
         c.W.affine_weight('D',13,0,8)==0,'exact CURRENT D13 center-zero role')
    domain(immediate,lambda p:(p['side'],p['profile']),
           set(product((0,1),('inner4','next3'))),'complete four immediate keys')
    domain(finite,lambda p:(p['side'],p['attached_old_leaf'],p['spectator_A1']),
           {(s,leaf,bool(e)) for s in (0,1) for leaf in ((9,10) if s==0 else (11,12)) for e in (0,1)},
           'complete eight finite coupon keys')
    domain(mixed,lambda p:p['type'],{'A3','D4','D5','D6','D7'},'complete five energy-one types')
    need(type(branches) is list and len(branches)==17,'complete17 zero-composition interfaces')
    immediate_zero=[p for p in branches if p['branch']=='immediate_affine']
    finite_zero=[p for p in branches if p['branch']=='whole_finite_coupon']
    mixed_zero=[p for p in branches if p['branch']=='mixed_energy_one']
    domain(immediate_zero,lambda p:(p['side'],p['profile']),
           {(p['side'],p['profile']) for p in immediate},'same four immediate zero-composition keys')
    domain(finite_zero,lambda p:(p['side'],p['old_leaf'],p['small']),
           {(p['side'],p['attached_old_leaf'],p['whole_small_type']) for p in finite},
           'same eight finite zero-composition keys')
    domain(mixed_zero,lambda p:p['type'],{p['type'] for p in mixed},'same five mixed zero-composition keys')
    out=[]
    for p in immediate:
        c.E.charge(c.FUEL)
        q=next(v for v in immediate_zero if (v['side'],v['profile'])==(p['side'],p['profile']))
        nodes=p['nodes'];marks=[nodes.index(m) for m in p['marks']]
        need(len(nodes)==len(set(nodes))==11 and len(marks)==2 and p['retained_payment']=='1',
             'complete physical immediate whole factor')
        white=[i for i in range(11) if i not in marks]
        value=c.E.corner(p['gram'],white,marks,c.W,c.FUEL)
        kind='B' if p['profile']=='inner4' else 'D'
        need(value['kind']=='ideal' and value['forest']==[[kind,11]] and F(value['weight'])==1,
             'CURRENT immediate B11/D11 whole-factor unit')
        need(p['canonical_pair'] in value['marked_images'] and
             q['type']==p['whole_factor'] and q['canonical_pair']==p['canonical_pair'] and
             q['distinct_payment']=='1','physical immediate/source/composition image association')
        out.append(dict(key=['immediate',p['side'],p['profile']],current=value,
                        original_zero_contribution_used=False))
    for p in finite:
        c.E.charge(c.FUEL)
        need(type(p['spectator_A1']) is bool,'finite coupon variant has exact boolean type')
        side=p['side'];leaf=p['attached_old_leaf'];extra=int(p['spectator_A1'])
        g,white,marks=physical_coupon(side,leaf,extra)
        small='A3' if extra else 'A2';role='middle' if extra else 'endpoint'
        need(p['marks']==marks and p['whole_small_type']==small and p['small_mark_role']==role and
             p['whole_large_type']=='D10' and p['large_mark_role']=='fork_center',
             'exact physical finite coupon/role association')
        q=next(v for v in finite_zero if (v['side'],v['old_leaf'],v['small'])==(side,leaf,small))
        need(q['small_role']==role and q['large']=='D10' and q['large_role']=='fork_center' and
             q['original_ideal_payment']=='0' and q['distinct_payment']==q['new_total']=='1',
             'same distinct coupon after zero composition')
        value=c.E.corner(g,white,marks,c.W,c.FUEL)
        physical=c.E.matrix(g);parts=c.E.components(physical,sorted(set(white)|set(marks)),c.FUEL)
        need(sorted(map(tuple,parts))==sorted(tuple(sorted(z)) for z in p['components']) and
             p['signature']==[len(g)-1,0,0],'same complete physical source corner partition')
        small_nodes=next(z for z in parts if 13 in z);large_nodes=next(z for z in parts if marks[0] in z)
        st,si,sa=c.E.classify(physical,small_nodes,[13],c.W,c.FUEL)
        lt,li,la=c.E.classify(physical,large_nodes,[marks[0]],c.W,c.FUEL)
        need(st==('A',2+extra) and not sa and
             set(si)==({(1,)} if extra else {(0,),(1,)}) and
             lt==('D',10) and not la and set(li)=={(7,)},
             'complete native A2/A3-middle and D10-center images')
        need(value['kind']=='finite' and sorted(map(tuple,value['forest']))==sorted([st,lt]) and
             F(value['weight'])==1,'CURRENT whole A2/A3-middle D10-center unit')
        out.append(dict(key=['finite',side,leaf,small],current=value,
                        source_whole_factor_premise='Accepted extreme-side vanishing/energy-one equality seals both factors; graph matching alone does not seal them.',
                        original_zero_contribution_used=False,legacy_half_payment_used=False))
    for p in mixed:
        rank=3 if p['type']=='A3' else int(p['type'][1:])
        q=next(v for v in mixed_zero if v['type']==p['type'])
        need(F(p['projection'])==F(rank,4)<2 and p['energy']=='1' and
             p['opposite_black_path']==4 and p['payment']=='1' and
             q['path_length']==4 and q['distinct_payment']=='1',
             'exact whole energy-one type/current short4 branch')
        out.append(dict(key=['mixed',p['type']],required_current_basis='finite57.current.short4_intrinsic_unit.v1',
            actual_path_length_bound=4,current_lower='1',
            geometry_premise='Accepted complete energy-one whole spectator classification and actual white-interior c-to-d path; no native/rank sampling.',
            original_zero_contribution_used=False))
    domain(coupon_updates,lambda p:(p['small'],p['small_role']),
           {('whole_A2','either_endpoint'),('whole_A3','unique_middle')},'exact two upgraded coupon roles')
    for p in coupon_updates:
        need(p['large']=='whole_D10' and p['large_role']=='fork_center' and p['new_weight']=='1',
             'same unconditional whole finite coupon upgrade')
    need(len(out)==17,'every distinct universal branch is current-bound')
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=('d13-center-zero',),required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args=parser.parse_args();own=str(Path(__file__).resolve().relative_to(ROOT));source=raw(own)
    a,r,c,first_extra,extra=prepare();basis=short4(c);branches=bind(c)
    for p in c.READ:need(c.digest(p)==c.PINS[p],'consumed universal source/input changed')
    need(sha256(raw(LOADER)).hexdigest()==LOADER_SHA and sha256(raw(a.ROUTE)).hexdigest()==a.ROUTE_SHA and
         sha256(raw(r.COLLECTOR)).hexdigest()==r.COLLECTOR_SHA,'shared universal loader/controller changed')
    need(raw(a.EXTRA)==first_extra and raw(EXTRA)==extra and raw(own)==source and
         sha256(c.raw(c.MANIFEST)).hexdigest()==c.MANIFEST_SHA,'universal source/manifest changed')
    payload=json.dumps(dict(status='PASS_CURRENT_D13_CENTER_ZERO_OUTPUT_BINDING_AT_EXPLICIT_PREMISES',
        batch=args.batch,source_sha256=sha256(source).hexdigest(),short4_current_interface=basis,
        universal_provider=dict(theorem='finite57.current.D13_center_zero.whole_face.v1',
            affine_kind='D',affine_rank=12,whole_order=13,canonical_marked_images=[[0,8]],intrinsic_current='0',
            distinct_current_alternatives=branches,requires_independent_S8=False,
            premises=['Complete actual all-native R8 profile/determinant theorem',
                      'Complete whole energy-one supporting-forest classification',
                      'Exact saved source matrix/kernel/role interfaces, not regenerated']),
        consumed_inputs={p:c.PINS[p] for p in sorted(c.READ)},fuel_used=c.FUEL[0],
        independent_S8_verified=False,current_terminal_bridge_complete=False,global_finite57_verified=False),sort_keys=True)
    need(len(payload.encode())+1<=MAX_OUTPUT,'complete universal output cap; never truncate branches')
    print(payload)


if __name__=='__main__':main()
