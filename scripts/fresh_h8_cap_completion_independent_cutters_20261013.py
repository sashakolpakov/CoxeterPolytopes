#!/usr/bin/env python3
"""Independent exact replay of the optional-polar core's full-bank cutters.

This checks witness geometry and every1035-node cutter decision. It does not
reclassify the core compatibility graph or establish the Boolean conclusion.
"""
import argparse
import json
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import Z,O,parse,enc,add,sub,scale,mul,sign,dot,mv,pair,matrix,inverse
from scripts.fresh_h8_alternate_core_angles_independent_20261013 import load,qsign
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results/fresh_h8_cap_completion_independent_cutters_20261013.json'

def run(verify=False):
    core_path='data/enumerations/fresh_h8_cap_completion_20261013_optional_core.json'
    core=json.loads((ROOT/core_path).read_text());source=json.loads((ROOT/core['source']).read_text())
    raw=json.loads((ROOT/source['new_cap_source']).read_text());base=json.loads((ROOT/raw['base_source']).read_text())
    G=matrix(base['seed_gram']);I=inverse(G);assert I==matrix(base['seed_gram_inverse'])
    center=tuple(map(parse,base['strict_center_y']));assert all(sign(x)>0 for x in center) and sign(pair(center,I,center))<0
    _,_,nodes=load()
    for row in raw['caps']:
        a=tuple(map(parse,row['coefficients_rational']));b=tuple(map(parse,row['coefficients_radical']));d=parse(row['discriminant'])
        nodes.append({'a':a,'b':b,'d':d,'q':scale(O,2),'pa':mv(G,a),'pb':mv(G,b)})
    assert len(nodes)==source['node_count']==1035
    for n in nodes:
        q=(add(dot(n['a'],n['pa']),mul(n['d'],dot(n['b'],n['pb']))),scale(dot(n['a'],n['pb']),2))
        assert qsign((sub(q[0],n['q']),q[1]),n['d'])==0 and sign(n['q'])>0
        assert qsign((dot(n['a'],center),dot(n['b'],center)),n['d'])>0
    rows=[]
    assert core['clauses']==[source['escape_clauses'][i]for i in core['clause_indices']]
    for index,row in zip(core['clause_indices'],core['clauses']):
        y=tuple(map(parse,row['y']));norm=pair(y,I,y);cp=pair(y,I,center)
        assert all(sign(t)>0 for t in y) and sign(norm)>0 and sign(cp)<0
        assert norm==parse(row['norm']) and cp==parse(row['center_pairing'])
        cutters=[]
        for j,n in enumerate(nodes):
            p=(dot(n['a'],y),dot(n['b'],y))
            if qsign(p,n['d'])<0:cutters.append(j)
        assert cutters==row['complete_dictionary_cutters']
        rows.append({'source_clause':index,'strict_spacelike_norm':enc(norm),'future_center_pairing':enc(cp),
            'complete_dictionary_cutters':cutters})
    union=sorted(set(i for row in rows for i in row['complete_dictionary_cutters']))
    assert union==core['core_nodes'] and len(union)==147 and len(rows)==13
    result={'status':'PASS_COMPLETE_CORE_WITNESS_CUTTERS','source':core_path,
        'dictionary_size':len(nodes),'core_witnesses':len(rows),'complete_cutter_sign_checks':len(nodes)*len(rows),
        'union_nodes':union,'witness_checks':rows,
        'scope':'Exact strict spacelike witness geometry, common future orientation, and all1035 cutter decisions for each of the13 core clauses. This independent replay does not certify graph completeness, compact UNSAT, or finite-volume exclusion.'}
    if verify:assert json.loads(OUTPUT.read_text())==result
    else:OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k]for k in ('status','dictionary_size','core_witnesses','complete_cutter_sign_checks')}));return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
