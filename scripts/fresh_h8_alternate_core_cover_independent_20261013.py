#!/usr/bin/env python3
"""Replay complete 493-node escape cutters and the finite-volume core proof."""
import argparse
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
from scripts.fresh_h8_alternate_core_angles_independent_20261013 import load,qsign
from scripts.fresh_h4_blocks_independent_verify import Z,parse,enc,sign,dot,pair,matrix

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results/fresh_h8_alternate_core_cover_independent_20261013.json'

def dpll(clauses,adjacency,include_parallel):
    n=len(adjacency);hit=[sum(1<<j for j,c in enumerate(clauses) if (c>>i)&1) for i in range(n)]
    states=[];ids={}
    @lru_cache(None)
    def solve(available,need):
        key=available,need;state_id=len(states);ids[key]=state_id
        row={'id':state_id,'available_mask':str(available),'remaining_clauses_mask':str(need),'branches':[]};states.append(row)
        if not need:
            row['satisfiable']=True;return True
        clause=min((j for j in range(len(clauses)) if (need>>j)&1),key=lambda j:(available&clauses[j]).bit_count())
        row['branch_clause']=clause;options=available&clauses[clause]
        while options:
            bit=options&-options;options-=bit;i=bit.bit_length()-1
            child=(available&adjacency[i],need&~hit[i]);answer=solve(*child)
            row['branches'].append({'selected_local_node':i,'child_state':ids[child]})
            if answer:row['satisfiable']=True;return True
        row['satisfiable']=False;return False
    answer=solve((1<<n)-1,(1<<len(clauses))-1)
    assert not answer
    # Replay every recorded branch without trusting the recursive solver.
    for row in states:
        available=int(row['available_mask']);need=int(row['remaining_clauses_mask']);j=row['branch_clause']
        assert need and (need>>j)&1 and not row['satisfiable']
        options=available&clauses[j]
        assert [b['selected_local_node'] for b in row['branches']]==[i for i in range(n) if (options>>i)&1]
        for branch in row['branches']:
            i=branch['selected_local_node'];child=states[branch['child_state']]
            assert int(child['available_mask'])==available&adjacency[i]
            assert int(child['remaining_clauses_mask'])==need&~hit[i]
            assert int(child['remaining_clauses_mask'])!=need and not child['satisfiable']
    return {'status':'UNSAT','parallel_pairs_allowed':include_parallel,'distinct_states':len(states),'proof_states':states}

def run(verify=False):
    graph,core,nodes=load()
    angle_path='results/fresh_h8_alternate_core_angles_independent_20261013.json'
    angle=json.loads((ROOT/angle_path).read_text())
    assert angle['status']=='COMPLETE_CORE_COMPACT_COMPATIBILITY_GRAPH' and not angle['missing_compatible_edges']
    assert angle['core_nodes']==core['union_nodes'] and angle['checked_pairs']==1431
    expected=set(combinations(core['union_nodes'],2))
    assert {(r['i'],r['j']) for r in angle['pair_checks']}==expected
    base=json.loads((ROOT/graph['base_source']).read_text());I=matrix(base['seed_gram_inverse']);center=tuple(map(parse,base['strict_center_y']))
    assert all(sign(x)>0 for x in center) and sign(pair(center,I,center))<0
    assert all(qsign((dot(n['a'],center),dot(n['b'],center)),n['d'])>0 for n in nodes)
    old=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_sat_20261013.json').read_text())
    new=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_graph_sat_20261013.json').read_text())
    rows=old['escape_clauses']+new['new_escape_clauses'];checks=[]
    for index in core['clause_indices']:
        row=rows[index];certificate=row['certificate'];y=tuple(map(parse,certificate['y']))
        q=pair(y,I,y);cp=pair(y,I,center)
        assert all(sign(x)>0 for x in y) and sign(q)>0 and sign(cp)<0
        assert q==parse(certificate['norm']) and cp==parse(certificate['center_pairing'])
        cutters=[]
        for j,node in enumerate(nodes):
            value=(dot(node['a'],y),dot(node['b'],y))
            if qsign(value,node['d'])<0:cutters.append(j)
        assert cutters==row['cutters']
        checks.append({'clause_index':index,'source':row['source'],'escape_y':certificate['y'],
            'escape_norm':enc(q),'escape_center_pairing':enc(cp),'complete_cutter_indices':cutters})
    union=sorted(set(j for row in checks for j in row['complete_cutter_indices']))
    assert union==core['union_nodes'] and [len(row['complete_cutter_indices']) for row in checks]==core['cutter_counts']
    local={j:i for i,j in enumerate(union)}
    clauses=[sum(1<<local[j] for j in row['complete_cutter_indices']) for row in checks]
    solutions=[]
    for parallel in (False,True):
        adjacency=[0]*len(union)
        for row in angle['pair_checks']:
            if row['kind'] not in (('finite','dotted','parallel') if parallel else ('finite','dotted')):continue
            i,j=local[row['i']],local[row['j']];adjacency[i]|=1<<j;adjacency[j]|=1<<i
        solutions.append(dpll(clauses,adjacency,parallel))
    result={'status':'EXACT_FINITE_BANK_FINITE_VOLUME_OBSTRUCTION',
        'angle_audit':angle_path,'dictionary_nodes':len(nodes),'core_nodes':union,
        'common_strict_interior_y':base['strict_center_y'],'common_interior_norm':enc(pair(center,I,center)),
        'five_complete_escape_clauses':checks,'parallel_pairs':[[r['i'],r['j']] for r in angle['pair_checks'] if r['kind']=='parallel'],
        'clause_masks_local_order':list(map(str,clauses)),'proofs':solutions,
        'scope':'No selection from the complete saved493-cap bank, retaining the ten fixed seed walls, can be a finite-volume Coxeter polyhedron. All possible finite Coxeter angles and parallel pairs on the54-node witness core are included. The finite numerical proposal bank is not an exhaustive search of unrestricted caps or H9 polytopes.'}
    if verify:assert json.loads(OUTPUT.read_text())==result
    else:OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k in ('status','dictionary_nodes','parallel_pairs','scope')},indent=2))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
