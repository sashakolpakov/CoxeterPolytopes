#!/usr/bin/env python3
"""Global40 reoptimization with all admissible finite and unpaid-dotted first stops.
The full white-extension closure remains outstanding; this is a proposal model.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_global_discount_20261016 as G
M,X,ROOT=G.M,G.X,G.ROOT
TAG='fresh_sqrt5_global_b1_first_repair_20261016';N=G.N;TARGET=40
FINITE='results/fresh_sqrt5_global_b1_finite_first_20261016.json'
DOTTED='results/fresh_sqrt5_global_b1_dotted_first_completions_20261016.json'
OLD_RESULT='results/fresh_sqrt5_global_discount_20261016.json'
OLD_CHECKPOINT=G.CHECKPOINT
INITIAL=f'data/enumerations/{TAG}_initial.json'
CHECKPOINT=Path('/private/tmp/'+TAG+'_checkpoint.json')
OLD_FACES=G.faces

def install():
    G.install();G.TAG=TAG;G.INITIAL=INITIAL;G.CHECKPOINT=CHECKPOINT

def prepare_initial():
    path=ROOT/INITIAL
    if path.exists():return
    old=json.loads(OLD_CHECKPOINT.read_text());vertices=[]
    for rec in old['vertex_cuts']:
        front=rec['metadata'];p=G.literal_profile(front);assert p==G.assembled_profile(front)
        row={j:z for j,z in p.items()if j>=0};row[N]=F(-1)
        assert row=={int(j):F(z)for j,z in rec['row']}and F(rec['rhs'])==-p.get(-1,F(0))
        vertices.append(rec)
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_global_discount_20261016.py',OLD_RESULT]
    path.write_text(json.dumps(dict(status='LITERAL_RECOUNTED_GLOBAL_VERTEX_HANDOFF_NOT_A_CERTIFICATE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},vertex_cuts=vertices,
        steps=[],latest_numeric_weights=None),separators=(',',':'))+'\n')
    print('FIRST_REPAIR_INITIAL_LITERAL_FORESTS',len(vertices),flush=True)

def faces():
    prior,counts=OLD_FACES();obligations={(M.key(p),t):meta for p,t,meta in prior}
    finite=json.loads((ROOT/FINITE).read_text());dotted=json.loads((ROOT/DOTTED).read_text())
    frows=[(r['index'],r['row'])for r in finite['records']if r['classification']['status']=='ACCEPT_SINGULAR_STOP']
    assert len(frows)==196
    drows=[]
    for rec in dotted['admissible_unpaid_completions']:
        full=dotted['records'][rec['record']];assert full['index']==rec['record']
        candidate=full['result']['candidates'][rec['candidate']]
        assert candidate['status']=='PASS_NECESSARY_ARITHMETIC_QUADRILATERAL_COMPLETION'
        drows.append((rec['record'],full['row']))
    assert len(drows)==14
    counts.update(finite_first_singular_rows=len(frows),dotted_first_singular_rows=len(drows),
        finite_first_common_contexts=0,dotted_first_common_contexts=0)
    oldweights=list(map(F,json.loads((ROOT/OLD_RESULT).read_text())['certificate']['primal']))
    old_failures=[]
    for family,rows in [('finite_first',frows),('dotted_first',drows)]:
        for pos,(index,row)in enumerate(rows):
            es=[(i,j,m)for(i,j),m in zip(combinations(range(row['nodes']),2),row['labels'])]
            compiled=[X.compile_corner(row['white']+[a,b],es,[a,b],[None])for a,b in product((0,2),(1,3))]
            minimum=None;worst=None;failed=0
            for ext in X.external_contexts(23-row['nodes']):
                p,tables=X.add_bonus(compiled,ext)
                for table in tables:p=M.plus(p,table[None])
                obligations.setdefault((M.key(p),F(1)),dict(kind='quad',family=family,index=index,
                    producer_status='ADMISSIBLE_SINGULAR_FIRST_EXTENSION_STOP',external_counts=list(ext)))
                counts[family+'_common_contexts']+=1
                value=M.value(p,oldweights);failed+=value<1
                if minimum is None or value<minimum:minimum=value;worst=ext
            if failed:old_failures.append(dict(family=family,index=index,minimum=str(minimum),worst_external_counts=list(worst),failed_contexts=failed))
            if pos%25==0:print('FIRST_REPAIR_FACE',family,pos,len(rows),'profiles',len(obligations),flush=True)
    counts['old_rule_unpaid_stops']=len(old_failures)
    return [(dict(p),t,meta)for(p,t),meta in obligations.items()],counts,old_failures

def run(saved=None):
    install();obligations,counts,failures=faces()
    print('GLOBAL_FIRST_REPAIR_MODEL',N,'variables',len(obligations),'profiles',counts,flush=True)
    cert,steps=G.solve(obligations)if saved is None else(saved['certificate'],saved['steps'])
    proof=G.verify(obligations,cert)
    paths=[f'scripts/{TAG}.py',INITIAL,'scripts/fresh_sqrt5_global_discount_20261016.py',OLD_RESULT,FINITE,DOTTED,
           f'results/{X.A.TRI}_rank22.json',f'results/{X.A.QUAD}.json']
    old=json.loads((ROOT/OLD_RESULT).read_text())
    return dict(status='PASS_EXACT_GLOBAL_FIRST_STOP_REPAIR_MODEL',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},dimension=21,target='40',variables=N,
        variable_order=old['variable_order'],face_counts=counts,unique_face_profiles=len(obligations),
        original_rule_first_stop_failures=failures,certificate=cert,verification=proof,steps=steps,
        vertex_context_order=old['vertex_context_order'],triangle_external_context_order=old['triangle_external_context_order'],
        quadrilateral_external_context_order=old['quadrilateral_external_context_order'],
        scope='Same global vertex payment family, now includes all196 admissible finite first-extension singular stopping classes and all14 admissible singular completions among unpaid dotted first extensions. Every added row uses all common external spectator contexts, in addition to all996 triangles,37 original paid stops and1979 zero-payment quad rows. No first-extension stopping row is removed on account of payment or guessed geometry. Full deeper finite/dotted white-extension closure and all remaining B1 quadrilateral coverage are still outstanding. This is a restricted model result only, not a dimension bound.',
        optimizer_scope='Fixed budget40 and minimum sum of nonnegative bonus coefficients; if finite cuts are infeasible, an exact rational weak-duality lower bound is extracted. Every inherited vertex cut is independently reassembled and literally recounted before use; exact feasible results require full global finite-forest DP and rational face checks.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();install()
    if not a.verify:
        prepare_initial()
        if not CHECKPOINT.exists():CHECKPOINT.write_bytes((ROOT/INITIAL).read_bytes())
    p=ROOT/f'results/{TAG}.json';saved=json.loads(p.read_text())if a.verify else None
    result=json.loads(json.dumps(run(saved),default=str))
    if a.verify:assert result==saved
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    number=result['verification'].get('maximum',result['verification'].get('lower_bound'))
    print(result['status'],result['certificate']['kind'],number,flush=True)
