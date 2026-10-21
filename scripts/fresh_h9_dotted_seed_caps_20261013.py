#!/usr/bin/env python3
"""Discovery of cap banks including a fixed dotted seed pairing magnitude.

Choose D^2=4+4cos(pi/m). Two orbit members differing by one transposition
then have pairing 2-D^2/2=-2cos(pi/m), removing a self-orbit obstruction
of the finite-seed-pairing alphabets. Numerical proposals only.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import root
import z3

import fresh_h9_multiblock_seed_20261013 as base
import fresh_h9_three_parameter_seed_20261013 as helper

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h9_three_parameter_seed_20261013_scout.json'


def cover_search(G,R,round_limit=120):
    """Clique constraints plus separating spacelike-point coverage clauses."""
    C=base.compat(R@G@R.T);n=len(R);variables=[z3.Bool(f'cap{i}')for i in range(n)]
    solver=z3.Solver();solver.set(timeout=2000)
    for i,j in zip(*np.where(np.triu(~C,1)),strict=True):
        solver.add(z3.Or(z3.Not(variables[int(i)]),z3.Not(variables[int(j)])))
    val,V=np.linalg.eigh(G);idx=list(np.flatnonzero(val < -1e-7))+list(np.flatnonzero(val > 1e-7))
    N=V[:,idx]*np.sqrt(abs(val[idx]));N[:,0]*=np.sign(sum(N[:,0]));caps=R@N
    clauses=[];best=None
    for attempt in range(round_limit):
        result=solver.check()
        if result==z3.unsat:
            return {'status':'NUMERICAL_COVER_CONSTRAINTS_UNSAT','rounds':attempt,
                    'clauses':clauses,'best_selection':best,'compatibility_pairs':int(C.sum()//2)}
        if result!=z3.sat:
            return {'status':'SOLVER_TIMEOUT','rounds':attempt,'clauses':clauses,'best_selection':best}
        model=solver.model();selected=[i for i,v in enumerate(variables)if z3.is_true(model.eval(v,model_completion=True))]
        available=np.ones(n,dtype=bool)
        if selected:available=np.all(C[selected],axis=0)
        while np.any(available):
            ids=np.flatnonzero(available);j=int(ids[np.argmax(C[np.ix_(ids,ids)].sum(axis=1))])
            selected.append(j);available&=C[j]
        geo=base.geometry(G,R,selected)
        if best is None or geo['max_klein_square']<best['geometry']['max_klein_square']:
            best={'cap_indices':selected,'geometry':geo}
        if geo['outside']==0:
            return {'status':'NUMERICAL_FINITE_VOLUME_CANDIDATE','rounds':attempt+1,
                    'clauses':clauses,'best_selection':best,'compatibility_pairs':int(C.sum()//2)}
        point=np.array(geo['worst_vertex']);cuts=np.flatnonzero(caps[:,1:]@point-caps[:,0]>1e-7).tolist()
        clauses.append({'spacelike_point':point.tolist(),'klein_square':geo['max_klein_square'],'cutting_caps':cuts})
        if not cuts:
            return {'status':'NUMERICAL_SHARED_ESCAPE','rounds':attempt+1,'clauses':clauses,'best_selection':best}
        solver.add(z3.Or([variables[i]for i in cuts]))
    return {'status':'ROUND_LIMIT','rounds':round_limit,'clauses':clauses,'best_selection':best}


def screen(output):
    data=json.loads(output.read_text());base.VALS=np.array(data['pairing_magnitudes']);start=time.monotonic()
    bounded=0
    for i,row in enumerate(data['records']):
        if 'geometry_all_caps'in row:continue
        ns=row['groups'];G,R,covs,ii=helper.complete_bank(ns,np.array(row['metric']),base.profiles(ns))
        geo=helper.reduced_geometry(ns,G,R);geo.pop('N')
        row.update({'cap_count':len(R),'profile_indices':ii.tolist(),'geometry_all_caps':geo})
        if geo['outside']==0:
            bounded+=1;result=cover_search(G,R);row['cover_search']=result
            print({'record':i,'caps':len(R),'cover_status':result['status'],'rounds':result['rounds']},flush=True)
            if result['status']=='NUMERICAL_FINITE_VOLUME_CANDIDATE':
                result['selected_covectors']=covs[result['best_selection']['cap_indices']].tolist()
                data['status']='NUMERICAL_POSITIVE_CANDIDATE_PENDING_EXACT_RECONSTRUCTION'
                output.write_text(json.dumps(data,indent=2)+'\n')
                return
        if (i+1)%500==0:
            output.write_text(json.dumps(data,indent=2)+'\n')
            print({'screened':i+1,'total':len(data['records']),'bounded':bounded,'seconds':round(time.monotonic()-start,1)},flush=True)
    data['status']='NUMERICAL_DOTTED_SEED_CAP_SCREEN_COMPLETE'
    output.write_text(json.dumps(data,indent=2)+'\n')
    print({'done':len(data['records']),'all_bank_bounded':sum(x['geometry_all_caps']['outside']==0 for x in data['records'])},flush=True)


def run(output,trials,m,rng_seed):
    start=time.monotonic();rng=np.random.default_rng(rng_seed)
    D=np.sqrt(4+4*np.cos(np.pi/m));base.VALS=np.array([0.,1.,base.PHI,D])
    seeds=json.loads(SOURCE.read_text())['records'];tables={};starts=[]
    for ns in set(tuple(s['groups'])for s in seeds):
        P=base.profiles(ns);types,inds,tv,dv,sq,mult=P
        U=-.5*np.array([helper.symprod(t,t)for t in tv])
        V=-np.array([helper.symprod(d,d)for d in dv])/32;C=sq/8-2
        _,keep=np.unique(np.round(np.c_[U,V,C],10),axis=0,return_index=True)
        uses=np.array([any(types[g][i][1]==3 for g,i in enumerate(inds[k]))for k in keep])
        tables[ns]=(U[keep],V[keep],C[keep],keep,uses)
    for s in seeds:
        ns=tuple(s['groups']);U,V,C,keep,uses=tables[ns];r0=np.array(s['inverse_metric']);b,_=helper.matrices(r0,ns)
        errors=abs(U@r0+V@b+C);close=np.argsort(errors)[:32]
        dotted=np.flatnonzero(uses);dotted=dotted[np.argsort(errors[dotted])[:24]]
        starts.append((ns,r0,close,dotted))
    seen=set();points={};attempts=0;converged=0;singular=0
    for trial in range(trials):
        ns,r0,close,dotted=starts[int(rng.integers(len(starts)))];U0,V0,C0,keep,uses=tables[ns]
        ij=rng.choice(close if trial%10 else len(C0),2,replace=False)
        k=int(rng.choice(dotted));ids=tuple(sorted([int(ij[0]),int(ij[1]),k]))
        if len(set(ids))<3 or(ns,ids)in seen:continue
        seen.add((ns,ids));U,V,C=U0[list(ids)],V0[list(ids)],C0[list(ids)]
        def fun(r):
            try:
                b,db=helper.matrices(r,ns);return U@r+V@b+C,U+V@db
            except ValueError:return np.full(3,1e12),np.eye(3)*1e12
        if np.linalg.matrix_rank(fun(r0)[1],tol=1e-9)<3:singular+=1;continue
        attempts+=1;sol=root(fun,r0,jac=True,method='hybr',options={'xtol':1e-10,'maxfev':150})
        error=max(abs(fun(sol.x)[0]));M=helper.metric(sol.x,ns)
        if error>1e-8 or M is None:continue
        converged+=1;key=(ns,tuple(np.round(M.ravel(),7)))
        points.setdefault(key,{'groups':list(ns),'metric':M.tolist(),'inverse_metric':sol.x.tolist(),
                               'anchor_profiles':keep[list(ids)].tolist(),'max_equation_residual':float(error)})
        if attempts%5000==0:print({'m':m,'attempts':attempts,'metrics':len(points),'seconds':round(time.monotonic()-start,1)},flush=True)
    data={'status':'NUMERICAL_DOTTED_SEED_CAP_METRICS','source':str(SOURCE.relative_to(ROOT)),
          'pairing_magnitudes':base.VALS.tolist(),'transposition_angle_label':m,'trials':trials,'rng_seed':rng_seed,
          'distinct_triples':len(seen),'singular_starts_skipped':singular,'root_attempts':attempts,
          'valid_solutions':converged,'records':list(points.values()),
          'scope':'Random three-profile norm systems with at least one dotted-seed profile. Ninety percent use profiles numerically near norm2 at a saved old start; ten percent choose two unrestricted profiles. One start per triple, singular initial Jacobians skipped. Metric entries in(2,60], negative weighted-core eigenvalue in(-4,0). The dotted magnitude is fixed by its displayed transposition angle. Numerical only; no exhaustion or finite-volume claim.'}
    output.write_text(json.dumps(data,indent=2)+'\n');print({'done':len(points),'m':m,'seconds':round(time.monotonic()-start,1)},flush=True)


def symmetric(output,m,select):
    import fresh_h9_asymmetric_intersections_20261013 as intersections
    base.VALS=np.array([0.,1.,base.PHI,np.sqrt(4+4*np.cos(np.pi/m))])
    if not select:
        intersections.run(output,300,[[4,4]])
        return
    data=json.loads(output.read_text());base.VALS=np.array(data['finite_pairing_magnitudes']);records=[]
    for i,row in enumerate(data['records']):
        if row['geometry_all_caps']['outside']:continue
        G,R,_=base.bank(row['groups'],row['a'],row['b'],base.profiles(row['groups']),row['profile_indices'])
        result=cover_search(G,R);result.update({'source_record':i,'cap_count':len(R)});records.append(result)
        print(i,len(R),result['status'],result['rounds'],flush=True)
    target=output.with_name(output.name.replace('_scout.json','_cover.json'))
    target.write_text(json.dumps({'status':'NUMERICAL_COMPATIBLE_COVER_SEARCH','source':str(output.relative_to(ROOT)),
      'records':records,'scope':'Numerical compatibility and separating-point clauses for every whole-bank-bounded entry. UNSAT only excludes the discretized graphs and coverage clauses; no exact parameter or volume claim.'},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--m',type=int,default=4);p.add_argument('--trials',type=int,default=100000)
    p.add_argument('--seed',type=int,default=610134);p.add_argument('--output',type=Path);p.add_argument('--screen',action='store_true')
    p.add_argument('--symmetric',action='store_true')
    a=p.parse_args();prefix='dotted_symmetric'if a.symmetric else'dotted_seed_caps'
    output=a.output or ROOT/f'data/enumerations/fresh_h9_{prefix}_m{a.m}_20261013_scout.json'
    if a.symmetric:symmetric(output,a.m,a.screen)
    elif a.screen:screen(output)
    else:run(output,a.trials,a.m,a.seed)
