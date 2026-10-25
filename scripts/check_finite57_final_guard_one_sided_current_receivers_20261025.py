#!/usr/bin/env python3
"""SOURCE ONLY: the eight unstored rank58 one-sided guard receivers.

The complete native/rank15 energy-one classification and whole inactive
witnesses are explicit accepted proof inputs. This binds only the two
physical30-root models with four ordered inactive A15/D15 choices each.
It neither checks selected G5 nor closes the24 quiet branches or S8.
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
from itertools import product
import json
from pathlib import Path
import types

ROOT=Path(__file__).resolve().parents[1]
ROUTE='scripts/index_finite57_current_terminal_routes_syntax_safe_20261025.py'
ROUTE_SHA='43a0c4e986cbdb92d112771bf0edc1d128ce2d0feced30fffd9d9825dda7c29d'
FINAL='fresh_d13_three_rank15_context_transfer_20261018'
INPUTS={
 f'scripts/{FINAL}.py':'a95ea679854ba907ead8640964a258c3baeb3e998052a3e2e3ff20f230cce486',
 f'results/{FINAL}.md':'91e93c3f7ff5eb4e64b435eea2892bf20af08711aed131110f40e898a9322014',
 f'results/{FINAL}.json':'2c7d615cd38123bef20cffad68de8ab30baa91826c147a6bb084892590e8ca86',
}
PAIRS={(0,11),(0,12)}
MAX_OUTPUT=500*1024


def need(ok,message):
    if not ok:raise ValueError(message)


def raw(name):
    p=(ROOT/name).resolve()
    need(p.is_relative_to(ROOT) and p.is_file(),'bounded source path')
    with p.open('rb') as stream:data=stream.read(2*1024**2+1)
    need(0<len(data)<=2*1024**2,'bounded source read')
    return data


def prepare():
    encoded=raw(ROUTE);need(sha256(encoded).hexdigest()==ROUTE_SHA,'exact current route loader')
    r=types.ModuleType('finite57_final_guard_one_sided_route');r.__file__=str(ROOT/ROUTE)
    exec(compile(encoded,r.__file__,'exec'),r.__dict__)
    c=r.prepare()  # Enforces the already-held shared parent proof lock.
    need(c.FINAL==FINAL,'exact finalguard provider stem')
    for name,digest in INPUTS.items():
        if name in c.PINS:need(c.PINS[name]==digest,'finalguard/base pin collision')
        c.PINS[name]=digest
        need(c.digest(name)==digest,'exact finalguard proof/source/data identity')
        c.READ.add(name)
    return r,c


def physical(attachment):
    need(type(attachment) is int and attachment in (9,10),'exact active attachment')
    g=[['2' if i==j else '0' for j in range(30)] for i in range(30)]
    edges=[(i,i+1) for i in range(8)]+[(0,9),(0,10),(8,11),(8,12)]
    edges += [(15+i,16+i) for i in range(12)]+[(27,28),(27,29)]
    edges += [(13,attachment),(14,12),(13,15)]
    for i,j in edges:g[i][j]=g[j][i]='-1'
    return g


def source_model(c,p):
    attachment=p['active_attachment'];expected=physical(attachment)
    g=c.E.matrix(p['gram'])
    need(g==c.E.matrix(expected),'complete original physical30-root Gram identity')
    white=sorted(set(range(30))-{9,11,13,14})
    need(p['white']==white and len(p['corners'])==4 and
         [q['marks'] for q in p['corners']]==[[9,11],[11,13],[13,14],[14,9]],
         'complete physical original white and ordered actual-corner domain')
    need(p['quiet_corner']==[14,9] and p['new_paying_corner']==[13,14] and
         set(p['quiet_corner'])!=set(p['new_paying_corner']),
         'exact distinct quiet/opposite actual receiver association')
    need(p['quiet_half']==p['new_half']=='1/2' and p['total']=='1',
         'saved receiver claim is the exact half-plus-half statement')
    quiet=p['corners'][3];pay=p['corners'][2]
    parts={}
    for name,q in (('quiet',quiet),('pay',pay)):
        ids=sorted(set(white)|set(q['marks']))
        blocks=c.E.components(g,ids,c.FUEL)
        need(q['signature']==[28,0,0] and
             sorted(map(tuple,blocks))==sorted(tuple(b['nodes']) for b in q['components']),
             'same complete physical receiver partition and saved finite signature')
        parts[name]=blocks
    small=next(b for b in parts['quiet'] if set(quiet['marks'])<=set(b))
    rest=[b for b in parts['quiet'] if b!=small]
    need(len(parts['quiet'])==2 and len(small)==13 and len(rest)==1 and len(rest[0])==15,
         'quiet corner contains precisely marked D13 and active unmarked rank15')
    typ,images,affine=c.E.classify(g,small,quiet['marks'],c.W,c.FUEL)
    witness,_,wa=c.E.classify(g,rest[0],[],c.W,c.FUEL)
    need(typ==('D',13) and not affine and set(images)==PAIRS and
         witness==('D',15) and not wa,'exact quiet marked orbit and active D15 factor')
    roles=[]
    for b in parts['pay']:
        marks=[m for m in pay['marks'] if m in b]
        need(marks,'no unsealed unmarked component is credited at the paying corner')
        typ,images,affine=c.E.classify(g,b,marks,c.W,c.FUEL)
        need(not affine,'paying factors are finite')
        roles.append(dict(type=typ[0]+str(typ[1]),images=[list(v) for v in images]))
    canon=lambda values:sorted(values,key=c.canonical)
    required=[dict(type='A12',images=[[0],[11]]),dict(type='D16',images=[[0]])] if attachment==9 else [dict(type='D28',images=[[0,12]])]
    need(canon(roles)==canon(required)==canon(p['paying_images']),
         'complete actual native paying roles and saved provider association')
    return white,rest[0],roles


def bind(c):
    rows=c.field(FINAL,('one_sided_whole_models',))
    need(type(rows) is list and len(rows)==2 and
         [p['active_attachment'] for p in rows]==[9,10],'complete two physical one-sided model keys')
    need(c.W.finite_same_weight(('D',28),0,12)==0 and
         c.W.finite_cross_weight(('A',12),0,('D',16),0)==0 and
         c.W.finite_cross_weight(('A',12),11,('D',16),0)==0,
         'both new receiver roles have ordinary current value zero')
    out=[]
    for index,p in enumerate(rows):
        white,active,roles=source_model(c,p)
        for kinds in product(('A','D'),repeat=2):
            inactive=tuple((kind,15) for kind in kinds)
            quiet=c.E.corner(p['gram'],white,p['quiet_corner'],c.W,c.FUEL,
                extra=inactive,complete=True,sealed_witnesses=[active])
            pay=c.E.corner(p['gram'],white,p['new_paying_corner'],c.W,c.FUEL,
                extra=inactive,complete=True)
            need(quiet['kind']==pay['kind']=='finite' and quiet['complete'] and pay['complete'] and
                 not quiet['ignored_unsealed_components'] and not pay['ignored_unsealed_components'],
                 'all actual rank58 factors are supplied and no induced witness is credited')
            need(sum(rank for _,rank in quiet['forest'])==sum(rank for _,rank in pay['forest'])==58,
                 'both exact full actual finite ranks are58')
            need(F(quiet['weight'])==F(pay['weight'])==F(1,2),
                 'CURRENT quiet negative guard and opposite positive guard each pay exactly half')
            total=F(quiet['weight'])+F(pay['weight'])
            need(total==1,'CURRENT distinct actual receiver total is one')
            out.append(dict(key=[p['active_attachment'],*kinds],
                source_selector=['one_sided_whole_models',index],
                physical_model_sha256=sha256(c.canonical([p['gram'],white]).encode()).hexdigest(),
                inactive_whole_factors=[list(t) for t in inactive],
                active_whole_D15_nodes=active,paying_native_roles=roles,
                quiet_marks=p['quiet_corner'],opposite_marks=p['new_paying_corner'],
                quiet_current=quiet,opposite_current=pay,distinct_current_total=str(total),
                whole_witness_premise='Accepted rank15 energy-one theorem forces exactly one active whole D15; the other two whole A15/D15 factors are orthogonal to the entire star and remain whole at BOTH actual corners. No graph list or ambient leftover rank is substituted for this proof.'))
    need(len(out)==8 and len({tuple(p['key']) for p in out})==8,'exact eight labelled receiver keys')
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=('one-sided-eight',),required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args=parser.parse_args();own=str(Path(__file__).resolve().relative_to(ROOT));before=raw(own)
    r,c=prepare();receivers=bind(c)
    for name in c.READ:need(c.digest(name)==c.PINS[name],'consumed finalguard source/input changed')
    need(raw(own)==before and sha256(raw(ROUTE)).hexdigest()==ROUTE_SHA,'finalguard adapter/route changed')
    need(sha256(raw(r.COLLECTOR)).hexdigest()==r.COLLECTOR_SHA and
         sha256(c.raw(c.MANIFEST)).hexdigest()==c.MANIFEST_SHA,'collector/controller manifest changed')
    payload=json.dumps(dict(status='PASS_NAMED_FINAL_GUARD_EIGHT_ONE_SIDED_CURRENT_RECEIVERS',
        batch=args.batch,source_sha256=sha256(before).hexdigest(),receivers=receivers,
        consumed_inputs={p:c.PINS[p] for p in sorted(c.READ)},fuel_used=c.FUEL[0],
        scope='Only the eight unstored one-sided positive-context physical receivers, at accepted native/rank15 classification, saved inertia and whole-witness premises. Selected guards,24quiet branches,S8 and the global theorem remain separate.',
        rank15_classification_regenerated=False,whole_witness_geometry_proved_by_code=False,
        quiet24_complete=False,independent_S8_verified=False,
        current_terminal_bridge_complete=False,global_finite57_verified=False),sort_keys=True)
    need(len(payload.encode())+1<=MAX_OUTPUT,'complete output byte cap; no truncated receivers')
    print(payload)


if __name__=='__main__':main()
