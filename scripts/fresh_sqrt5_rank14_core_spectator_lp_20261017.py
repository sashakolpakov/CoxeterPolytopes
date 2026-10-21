#!/usr/bin/env python3
"""One bounded direct rank14 simplex LP with complete finite vertex budgets."""
import argparse, hashlib, json
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
import fresh_sqrt5_rank14_core_spectator_domain_20261017 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank14_core_spectator_lp_20261017'


def domain():
    prior=B.read(P.TAG)
    models=[(name,g) for name,g in B.C.models() if len(g)<=14]
    model=dict(models)
    codes={name:json.loads(B.dump(B.C.code(g))) for name,g in models}
    @lru_cache(None)
    def forests(rank,first=0):
        if rank==0:return ((),)
        return tuple(sorted(tuple(sorted((name,)+tail))
            for j,(name,g) in enumerate(models[first:],first) if len(g)<=rank
            for tail in forests(rank-len(g),j)))
    all_forests=forests(14)
    assert len(all_forests)==prior['counts']['all_finite_rank14_forests']==11751
    variables=[];vertices=[];lookup={}
    for forest in all_forests:
        g={}
        for name in forest:
            h=model[name];offset=len(g)
            for i,row in h.items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==14
        pieces=B.C.pieces(g);orbits=defaultdict(list)
        for pair in combinations(range(14),2):
            orbits[B.dump(sorted(B.C.code(p,pair) for p in pieces))].append(list(pair))
        vid=len(vertices);vv=[]
        for code,pairs in sorted(orbits.items()):
            index=len(variables);lookup[forest,code]=index;vv.append([index,len(pairs)])
            variables.append(dict(index=index,vertex=vid,forest=list(forest),marked_forest_code=json.loads(code),
                                  multiplicity=len(pairs),literal_pairs=pairs))
        assert sum(m for i,m in vv)==91
        vertices.append(dict(index=vid,forest=list(forest),variables=vv))
        if len(vertices)%2000==0:print('RANK14_VERTEX_BUDGETS',len(vertices),len(variables),flush=True)
    faces=[];face_lookup={}
    for ctx in prior['contexts']:
        core=prior['cores'][ctx['core']];spectator=ctx['spectator_forest'];cc=Counter()
        for corner in core['corner_occurrences']:
            forest=tuple(sorted(corner['forest']+spectator))
            code=B.dump(sorted(corner['marked_forest_code']+[codes[name] for name in spectator]))
            cc[lookup[forest,code]]+=1
        coeff=tuple(sorted(cc.items()));key=(core['required'],coeff)
        if key not in face_lookup:
            face_lookup[key]=len(faces)
            faces.append(dict(index=len(faces),required=core['required'],coefficients=[list(x) for x in coeff],contexts=[]))
        faces[face_lookup[key]]['contexts'].append(ctx['index'])
    active=sorted({i for f in faces for i,m in f['coefficients']})
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',
           f'scripts/{P.TAG}.py',f'results/{P.TAG}.json']
    return dict(sources={p:P.sha(p) for p in paths},dimension=14,target=14,
        variables=variables,vertices=vertices,faces=faces,active_variables=active,
        counts=dict(vertices=len(vertices),variables=len(variables),active_variables=len(active),
            zero_only_variables=len(variables)-len(active),distinct_face_constraints=len(faces),
            literal_simplex_stars=len(prior['contexts'])),
        definition='All178723 inherited-core/finite-spectator simplex stars in rank14, with their literal aggregate demands and all91-pair budgets for every11751 finite rank14 forest. No determinant partition. Variables absent from all face rows are fixed to zero during optimization and retained explicitly in the complete budget/certificate domain.')


def primal(d,x):
    best=None
    for denominator in (10**6,10**8,10**10):
        w={j:F(max(0,round(float(a)*denominator)),denominator) for j,a in x.items()}
        w={j:a for j,a in w.items() if a};scale=F(1)
        for f in d['faces']:
            total=sum((m*w.get(j,F(0)) for j,m in f['coefficients']),F(0))
            if not total:break
            scale=max(scale,F(f['required'])/total)
        else:
            w={j:a*scale for j,a in w.items()}
            budgets=[sum((m*w.get(j,F(0)) for j,m in v['variables']),F(0)) for v in d['vertices']]
            c=dict(kind='EXACT_RANK14_PARTIAL_PRIMAL',maximum=str(max(budgets)),scale=str(scale),
                weights=[[j,str(a)] for j,a in sorted(w.items())],vertex_budgets=list(map(str,budgets)),
                reaches_target14=max(budgets)<=14)
            if best is None or F(c['maximum'])<F(best['maximum']):best=c
    return best


def dual(d,values):
    c=D.recover_unminimized(d,values)
    c['kind']='EXACT_RANK14_PARTIAL_DUAL';c['excludes_target14']=F(c['bound'])>14
    del c['excludes_target12']
    return c


def verify(d,c):
    if c['kind']=='EXACT_RANK14_PARTIAL_PRIMAL':
        w={i:F(a) for i,a in c['weights']}
        assert len(w)==len(c['weights']) and all(a>0 for a in w.values())
        assert set(w)<=set(d['active_variables'])
        for f in d['faces']:
            assert sum((a*w.get(i,F(0)) for i,a in f['coefficients']),F(0))>=f['required']
        budgets=[sum((m*w.get(j,F(0)) for j,m in v['variables']),F(0)) for v in d['vertices']]
        assert list(map(str,budgets))==c['vertex_budgets'] and str(max(budgets))==c['maximum']
        assert c['reaches_target14']==(max(budgets)<=14)
    elif c['kind']=='EXACT_RANK14_PARTIAL_DUAL':
        a={i:F(x) for i,x in c['face_multipliers']}
        bound,beta,tight,num,den=D.sparse_dual_bound(d,a)
        assert c['bound']==str(bound) and c['objective_numerator']==str(num) and c['normalizer']==str(den)
        assert c['vertex_multipliers']==[[i,str(x)] for i,x in sorted(beta.items())]
        assert c['tight_orbits']=={str(i):x for i,x in tight.items()}
        assert c['excludes_target14']==(bound>14)
    else:raise AssertionError(c['kind'])


def optimize(d):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    active=d['active_variables'];lookup={i:j for j,i in enumerate(active)};N=len(active)
    rows=[];cols=[];values=[];rhs=[]
    for f in d['faces']:
        r=len(rhs);rhs.append(-f['required'])
        for j,m in f['coefficients']:rows.append(r);cols.append(lookup[j]);values.append(-m)
    for v in d['vertices']:
        r=len(rhs);rhs.append(0)
        for j,m in v['variables']:
            if j in lookup:rows.append(r);cols.append(lookup[j]);values.append(m)
        rows.append(r);cols.append(N);values.append(-1)
    matrix=coo_matrix((values,(rows,cols)),shape=(len(rhs),N+1)).tocsr()
    objective=np.zeros(N+1);objective[-1]=1
    print('RANK14_SINGLE_LP_START',d['counts'],flush=True)
    result=linprog(objective,A_ub=matrix,b_ub=np.array(rhs,dtype=float),bounds=(0,None),
                   method='highs',options={'time_limit':300})
    solver=dict(status=int(result.status),message=str(result.message),success=bool(result.success),
                objective=None if result.fun is None else float(result.fun),time_limit_seconds=300)
    certificates=[]
    if result.success:
        p=primal(d,{i:result.x[j] for j,i in enumerate(active)})
        if p is not None:certificates.append(p)
        certificates.append(dual(d,result.ineqlin.marginals[:len(d['faces'])]))
    for c in certificates:verify(d,c)
    return solver,certificates


def main():
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--solve-prepared',action='store_true');p.add_argument('--verify',action='store_true');args=p.parse_args()
    prepared=ROOT/f'results/{TAG}_domain.json';out=ROOT/f'results/{TAG}.json'
    if args.solve_prepared:
        d=json.loads(prepared.read_text())
        for path,h in d['sources'].items():assert P.sha(path)==h,path
    else:d=domain()
    if args.prepare:
        prepared.write_text(json.dumps(d,indent=2)+'\n');print('PASS_RANK14_LP_DOMAIN_PREPARED',d['counts'],flush=True);return
    if args.verify:
        saved=json.loads(out.read_text())
        for k,v in d.items():assert saved[k]==v,k
        for c in saved['certificates']:verify(d,c)
        d=saved
    else:
        solver,certificates=optimize(d)
        feasible=any(c.get('reaches_target14',False) for c in certificates)
        excluded=any(c.get('excludes_target14',False) for c in certificates)
        assert not(feasible and excluded)
        d.update(status='PASS_EXACT_RANK14_PARTIAL_RULE_AT14' if feasible else 'PASS_EXACT_RANK14_PARTIAL_MODEL_OBSTRUCTION' if excluded else 'UNRESOLVED_RANK14_PARTIAL_MODEL',
            solver=solver,certificates=certificates,
            scope='One300-second partial simplex model. All finite rank14 vertex budgets are present, but connected order15 cores and other missing cores, quadrilaterals and general faces are absent. A feasible rule is not a dimension theorem; an exact dual obstructs only this nonnegative whole-forest corner-payment model.')
        out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound'))) for c in d['certificates']],flush=True)


if __name__=='__main__':main()
