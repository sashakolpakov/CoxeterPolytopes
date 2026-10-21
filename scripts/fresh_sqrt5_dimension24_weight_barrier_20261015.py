#!/usr/bin/env python3
"""Exact dual obstruction to the frozen normalized-A1 weight model in rank24.

This is a restricted-model obstruction, not a dimension theorem. The existing
producer is imported read-only; its model is rebuilt with dimension24.
"""
import argparse,json,sys
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_complete_face_weights_20261015 as W

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_dimension24_weight_barrier_20261015'
OUT=ROOT/f'results/{TAG}.json'
W.BOUND=24;W.TARGET=F(54)

def sources():
    paths={f'scripts/{TAG}.py','results/fresh_sqrt5_complete_face_weights_20261015_n25.json'}
    for name,module in tuple(sys.modules.items()):
        if name.startswith('fresh_') and getattr(module,'__file__',None):
            p=Path(module.__file__).resolve()
            if ROOT in p.parents:paths.add(str(p.relative_to(ROOT)))
    queue=[p for p in paths if p.endswith('.json')]
    while queue:
        p=queue.pop();data=json.loads((ROOT/p).read_text())
        for q,h in data.get('sources',{}).items():
            assert sha256((ROOT/q).read_bytes()).hexdigest()==h,(p,q)
            if q not in paths:
                paths.add(q)
                if q.endswith('.json'):queue.append(q)
    return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}

def forest_meta(data,front):
    return dict(kind='forest',context=front['context'],forest=list(front['forest']),
                types=[data['components'][i]['type'] for i in front['forest']])

def dual_check(data,rows,rhs,multipliers):
    total=Counter();bound=F(0)
    assert len(rows)==len(rhs)==len(multipliers)
    for row,b,y in zip(rows,rhs,multipliers):
        assert y<=0
        bound+=b*y
        for j,v in row.items():total[j]+=v*y
    slack={j:F(j==data['zi'])-total[j] for j in range(data['zi']+1)}
    assert min(slack.values())>=0
    return bound,{j:v for j,v in slack.items() if v}

def solve(data):
    rows=list(data['rows']);rhs=list(data['rhs']);meta=list(data['meta']);cuts=set();steps=[]
    def add(front):
        key=front['context'],tuple(front['forest'])
        if key in cuts:return
        cuts.add(key);r,b=W.forestrow(data,*key)
        rows.append(r);rhs.append(b);meta.append(forest_meta(data,front))
    for front in W.maxforests(data,[F(0)]*data['zi']):add(front)
    for iteration in range(200):
        rr=[];cc=[];vv=[]
        for i,row in enumerate(rows):
            for j,v in row.items():rr.append(i);cc.append(j);vv.append(float(v))
        A=W.C.coo_matrix((vv,(rr,cc)),shape=(len(rows),data['zi']+1)).tocsr()
        objective=W.C.np.zeros(data['zi']+1);objective[-1]=1
        lp=W.C.linprog(objective,A_ub=A,b_ub=W.C.np.array(list(map(float,rhs))),bounds=(0,None),method='highs')
        assert lp.success,lp.message
        steps.append(dict(iteration=iteration,numerical_lower=float(lp.fun),forest_cuts=len(cuts)))
        print('dimension24 model LP',steps[-1],flush=True)
        if lp.fun>54+1e-7:
            for denominator in (1000,1000000,1000000000):
                y=[F(float(v)).limit_denominator(denominator) if abs(v)>1e-9 else F(0) for v in lp.ineqlin.marginals]
                try:bound,slack=dual_check(data,rows,rhs,y)
                except AssertionError:continue
                if bound>54:
                    selected=[]
                    for i,v in enumerate(y):
                        if not v:continue
                        selected.append(dict(multiplier=str(v),rhs=str(rhs[i]),
                            coefficients=[[j,str(a)]for j,a in sorted(rows[i].items())],metadata=meta[i]))
                    return dict(lower_bound=str(bound),terms=selected,
                        dual_slack=[[j,str(v)]for j,v in sorted(slack.items())],
                        numerical_steps=steps,optimality_claim=False)
            raise AssertionError('Numerical obstruction found, exact dual recovery remains open')
        x=[max(F(0),F(float(v)).limit_denominator(10**8))for v in lp.x]
        old=len(cuts)
        for front in W.maxforests(data,x):
            if float(front['maximum'])>lp.fun+1e-7:add(front)
        assert len(cuts)>old,'Model feasible at target; no dual obstruction'
    raise AssertionError('Bounded LP scout did not terminate')

def verify(data,certificate):
    rows=[];rhs=[];multipliers=[];face_lookup={
        (tuple(sorted(r.items())),b):i for i,(r,b)in enumerate(zip(data['rows'],data['rhs']))}
    for term in certificate['terms']:
        row={int(j):F(v)for j,v in term['coefficients']};b=F(term['rhs']);meta=term['metadata']
        if meta['kind']=='forest':
            expected,eb=W.forestrow(data,meta['context'],tuple(meta['forest']))
            assert row==expected and b==eb
            assert meta['types']==[data['components'][i]['type']for i in meta['forest']]
        else:
            i=face_lookup[tuple(sorted(row.items())),b];assert meta==data['meta'][i]
        rows.append(row);rhs.append(b);multipliers.append(F(term['multiplier']))
    lower,slack=dual_check(data,rows,rhs,multipliers)
    assert lower==F(certificate['lower_bound'])>54
    assert [[j,str(v)]for j,v in sorted(slack.items())]==certificate['dual_slack']
    return dict(lower_bound=str(lower),threshold='54',terms=len(rows),
        nonzero_slack_columns=len(slack),exact_weak_duality=True)

def run(saved=None):
    data=W.build();print('built',len(data['faces']),'faces',len(data['rows']),'constraints',data['zi'],'variables',flush=True)
    cert=solve(data)if saved is None else saved['certificate'];proof=verify(data,cert)
    used=sorted({j for t in cert['terms']for j,v in t['coefficients']if j!=data['zi']})
    faceids=sorted({w['face']for t in cert['terms']if t['metadata']['kind']=='face'for w in t['metadata']['witnesses']})
    return dict(status='PASS_EXACT_RESTRICTED_NORMALIZED_A1_MODEL_OBSTRUCTION_DIMENSION24',
        sources=sources(),dimension=24,faces=len(data['faces']),face_cases=len(data['faceconstraints']),
        unique_face_constraints=len(data['rows']),variables=data['zi'],
        used_variables=[dict(index=j,orbit=data['variables'][j][0],exact_A1_count=data['variables'][j][1])for j in used],
        referenced_faces=[dict(model_face_index=i,**data['faces'][i])for i in faceids],
        certificate=cert,verification=proof,
        scope='Necessary obstruction only to the frozen additive q2/H4 baseline model with nonnegative connected-pair bonuses indexed by the exact A1 count and normalized A1-to-component-root bonuses. All other disconnected pairs have weight zero. This does not exclude richer forest contexts, non-A1 cross weights, baseline reductions, or compact arithmetic polytopes in dimension24. Some face constraints are conservative necessary principal-star obligations, not independently asserted actual quadrilateral realizations.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    saved=json.loads(OUT.read_text())if args.verify else None
    result=json.loads(json.dumps(run(saved),default=str))
    if args.verify:assert result==saved
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['verification'])
