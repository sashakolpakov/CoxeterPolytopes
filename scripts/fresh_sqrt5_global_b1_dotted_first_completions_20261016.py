#!/usr/bin/env python3
"""Exact zero-further-contact dotted completions of unpaid first extensions."""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_dotted_completion_20261015 as D

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_global_b1_dotted_first_completions_20261016'
INPUT='results/fresh_sqrt5_global_b1_first_extensions_20261016.json'

def run():
    doc=json.loads((ROOT/INPUT).read_text())
    rows=[r for r in doc['failures']if r['parameter_type']=='dotted']
    assert len(rows)==122 and all(r['row']['labels'].count(0)==1 for r in rows)
    records=[];counts=Counter();outcomes=Counter();admissible=[]
    for i,r in enumerate(rows):
        z=D.evaluate(r['row']);counts[z['status']]+=1
        for j,c in enumerate(z.get('candidates',[])):
            if c['status']=='PASS_ARITHMETIC_NECESSARY_COMPLETION':
                G=[[D.K(*x)for x in rr]for rr in c['gram']]
                kernel=D.Q.kernel(G);assert len(kernel)==1
                signs=[x.sign()for x in kernel[0]]
                c['kernel']=[x.enc()for x in kernel[0]];c['kernel_signs']=signs
                if all(v>=0 for v in signs)or all(v<=0 for v in signs):
                    c['status']='REJECT_NONNEGATIVE_ROOT_DEPENDENCE'
                elif not(signs[0]==signs[2]and signs[1]==signs[3]and signs[0]*signs[1]==-1):
                    c['status']='UNRESOLVED_QUADRILATERAL_BLACK_KERNEL_PATTERN'
                else:
                    c['status']='PASS_NECESSARY_ARITHMETIC_QUADRILATERAL_COMPLETION'
                    admissible.append(dict(record=i,candidate=j,parent=r['index'],contacts=r['contacts'],
                                           nodes=r['nodes'],payment=r['minimum'],
                                           external_counts=r['worst_external_counts'],ratio=c['squared_dotted_ratio']))
            outcomes[c['status']]+=1
        records.append(dict(index=i,parent=r['index'],contacts=r['contacts'],row=r['row'],
                            worst_external_counts=r['worst_external_counts'],minimum_payment=r['minimum'],result=z))
        if i%10==0:print('DOTTED_B1',i,len(rows),dict(counts),dict(outcomes),flush=True)
    src=[f'scripts/{TAG}.py',INPUT,'scripts/fresh_sqrt5_dotted_completion_20261015.py',
         'scripts/fresh_sqrt5_singular_extension_frontier_20261015.py',
         'scripts/fresh_sqrt5_quad_contacts_20261015_review.py']
    return dict(status='PASS_EXACT_UNPAID_DOTTED_FIRST_EXTENSION_COMPLETION_SCREEN',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in src},
                input_count=len(rows),counts=dict(counts),candidate_counts=dict(outcomes),
                records=records,admissible_unpaid_completions=admissible,
                scope='Exact one-missing-entry corank1 completions with no further touching white roots, restricted to122 unpaid symbolic-dotted first extensions. Physical/conjugate inertia, cyclic integrality and black kernel alternation checked on surviving candidates. No arbitrary attachment closure or global compact realization; singular common bases remain explicitly unresolved.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['input_count'],len(d['admissible_unpaid_completions']),flush=True)
