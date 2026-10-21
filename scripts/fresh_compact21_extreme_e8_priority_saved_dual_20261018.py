#!/usr/bin/env python3
"""Stream the latest saved dual and exactly repair the proposed two-row change."""
import argparse
from collections import defaultdict
import fcntl
from fractions import Fraction as F
from hashlib import file_digest
from itertools import combinations
import json
from pathlib import Path
import resource
import signal
import sys
import fresh_compact21_sealed_e8_quad_triangle_activation_20261018 as E

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact21_extreme_e8_priority_saved_dual_20261018'
MODEL = 'results/fresh_sqrt5_dimension21_a4_source_orbits_20261016_cap2.json'


def objects(key, indent=2):
    active,index,item,size = False,-1,[],0
    start = ' '*indent+'"'+key+'": ['
    begin,end = ' '*(indent+2)+'{',' '*(indent+2)+'}'
    with (ROOT/MODEL).open() as stream:
        for line in stream:
            if not active:
                if line.startswith(start):
                    active=True
                continue
            if line.startswith(' '*indent+'],'):
                return
            if line.startswith(begin):
                index+=1
                item,size=[],0
            item.append(line)
            size+=len(line)
            assert size < 2*1024**2
            if line.startswith(end):
                yield index,json.loads(''.join(item).rstrip().rstrip(','))
                item,size=[],0


def run():
    extreme=[]
    for _,record in objects('connected_orbits'):
        if record['type']!='E8':
            continue
        ids=record['nodes']
        edges=dict(zip(combinations(ids,2),record['labels']))
        naming=E.canonical_e8_map(ids,edges)
        if sorted(naming[i] for i in record['marked'])==[0,6]:
            extreme.append(record['key'])
    assert len(extreme)==1
    columns=[]
    variables=0
    for index,record in objects('variable_order'):
        variables=index+1
        if record['orbit']==extreme[0] and all(record[k]==0 for k in
                ('A1_count','I25_context','H3_context','A4_context')):
            columns.append(index)
    assert variables==93045 and len(columns)==1
    beta=columns[0]
    loads=defaultdict(F)
    old_value=F()
    terms=0
    target=[]
    for index,term in objects('terms',4):
        y=F(term['multiplier'])
        assert y<0
        old_value+=y*F(term['rhs'])
        for j,a in term['coefficients']:
            assert 0<=j<=variables
            loads[j]+=y*F(a)
        if term['metadata']['kind']=='forest' and sorted(term['metadata']['types'])==['D5','E8','E8']:
            target.append((index,term))
        terms+=1
    assert terms==2804 and len(target)==1
    assert all(load<=F(j==variables) for j,load in loads.items())
    assert old_value==F(121786916150321,3043000000000)>40
    index,term=target[0]
    row={j:F(a)for j,a in term['coefficients']}
    mu=-F(term['multiplier'])
    assert row[beta]==2 and row[variables]==-1 and F(term['rhs'])==-36
    old_beta_load=loads[beta]
    active=max(F(),mu+old_beta_load/2)
    inactive=mu-active
    assert 0<=active<=mu and 0<=inactive<=mu
    loads[beta]+=2*inactive
    assert all(load<=F(j==variables) for j,load in loads.items())
    value=old_value-active
    # Active row has the old coefficients and rhs+1. Inactive row has
    # the old rhs and removes beta's multiplicity2. All other rows stay.
    direct=old_value-F(term['multiplier'])*F(term['rhs'])
    direct-=active*(F(term['rhs'])+1)+inactive*F(term['rhs'])
    assert direct==value
    sources=[f'scripts/{TAG}.py',MODEL,
             'scripts/fresh_compact21_sealed_e8_quad_triangle_activation_20261018.py',
             'results/fresh_compact21_extreme_e8_quad_consumer_completeness_20261018.json',
             'results/fresh_compact21_d5_mixed_fork_quad_exclusion_20261018.json']
    pins={}
    for path in sources:
        with (ROOT/path).open('rb') as stream:
            pins[path]=file_digest(stream,'sha256').hexdigest()
    return {'status':'PASS_EXACT_SAVED_DUAL_REPAIR_FOR_EXTREME_E8_PRIORITY_ROWS',
            'sources':pins,'variables':variables,'literal_dual_terms':terms,
            'touched_columns':len(loads),'extreme_E8_bonus_column':beta,
            'target_forest_term':index,'old_dual_value':str(old_value),
            'D5_2E8_mass':str(mu),'old_extreme_column_load':str(old_beta_load),
            'minimal_active_repair_mass':str(active),'inactive_repair_mass':str(inactive),
            'repaired_dual_value':str(value),'repaired_excess_above_40':str(value-40),
            'repaired_dual_remains_above_40':value>40,
            'every_literal_column_inequality_checked':True,
            'scope':'Exact repair of one saved dual for the stated two-row successor; '
                    'no optimization, no face-model rebuild, no new dimension theorem.',
            'geometry_review_required':True,'new_dimension_bound':None}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    signal.alarm(59)
    rss=lambda:resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
    assert rss()<1024**3
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        result=run()
        assert rss()<1024**3
        target=ROOT/f'results/{TAG}.json'
        if args.verify:
            assert json.loads(target.read_text())==result
        else:
            assert not target.exists()
            target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'repaired_dual':result['repaired_dual_value'],
                      'still_above_40':result['repaired_dual_remains_above_40'],
                      'peak_rss_bytes':rss(),'new_dimension_bound':None}))
