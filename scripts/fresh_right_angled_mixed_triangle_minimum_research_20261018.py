#!/usr/bin/env python3
"""One capped LP, then exact primal/dual verification of its small support."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse
import fcntl
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
import resource
import signal
import sys
import time
from check_right_angled_mixed_weighted_triangles_research_20261018 import ds_rows
from check_right_angled_mixed_opposition_moments_research_20261018 import partitions,clique_vector
import run_lean_batches_20261018 as R

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_triangle_minimum_research_20261018.json'


class L:
    def __init__(self,values):self.v=list(map(Q,values))
    def __add__(self,other):
        if not isinstance(other,L):other=L([0]*(len(self.v)-1)+[other])
        return L(a+b for a,b in zip(self.v,other.v))
    __radd__=__add__
    def __neg__(self):return L(-a for a in self.v)
    def __sub__(self,other):return self+-other if isinstance(other,L) else self+(-Q(other))
    def __rsub__(self,other):return -self+other
    def __mul__(self,scalar):return L(Q(scalar)*a for a in self.v)
    __rmul__=__mul__
    def __truediv__(self,scalar):return self*(1/Q(scalar))


def linear_system(rows,n,free_values=None):
    rows=[list(map(Q,row)) for row in rows];where=0
    pivots=[]
    for col in range(n):
        pivot=next((i for i in range(where,len(rows)) if rows[i][col]),None)
        if pivot is None:continue
        rows[where],rows[pivot]=rows[pivot],rows[where]
        divisor=rows[where][col];rows[where]=[a/divisor for a in rows[where]]
        for i,row in enumerate(rows):
            if i!=where and row[col]:
                factor=row[col];rows[i]=[a-factor*b for a,b in zip(row,rows[where])]
        pivots.append(col);where+=1
    assert all(not any(row[:n]) and not row[-1] for row in rows[where:]),'Inconsistent exact support'
    if len(pivots)<n and free_values is None:raise ValueError('Exact support system is rank deficient')
    result=[Q(0)]*n
    for j in set(range(n))-set(pivots):result[j]=Q(str(free_values[j])).limit_denominator(10**6)
    for i,j in enumerate(pivots):result[j]=rows[i][-1]-sum(rows[i][k]*result[k] for k in range(n) if k!=j)
    return result


def model():
    parts=list(partitions(12));names=['c']+[f'a{k}' for k in range(6,12)]+['D']+['w:'+','.join(map(str,p)) for p in parts]
    names+=['polygon:0,3','polygon:2,2','polygon:3,1','polygon:5,0']
    count=len(names);var=[L([int(i==j) for i in range(count)]+[0]) for j in range(count)]
    c=var[0];free=var[1:7];D=var[7];weights=var[8:85];polygons=var[85:]
    zero=L([0]*(count+1));known=[zero]*6+free+[c,zero]
    a=[sum((x*y for x,y in zip(row,known)),zero) for row in ds_rows(12)]+free
    f=a[0];vectors=[clique_vector(p) for p in parts]
    T=sum((v[2]*w for v,w in zip(vectors,weights)),zero)
    twoB=sum((v[3]*w for v,w in zip(vectors,weights)),zero);B=twoB/2
    C=sum((v[4]*w for v,w in zip(vectors,weights)),zero)/5
    K=sum((v[5]*w for v,w in zip(vectors,weights)),zero)
    rows=[];eq=[('finite_normalization',(f-1).v),('corner_normalization',(sum(weights,zero)-1).v),
                ('noncounted_polygon_count',(sum(polygons,zero)-a[2]+T).v)]
    add=lambda name,value:rows.append((name,value.v))
    beta=lambda k:comb(11,k-1)*2**(12-k)
    for k in range(6):add(f'face_positive_{k}',a[k])
    # Each noncounted polygon dominates one of these four minimal pairs;
    # arbitrary extra finite/ideal corners are represented by the slacks.
    add('polygon_finite_residual',66*f-T-2*polygons[1]-3*polygons[2]-5*polygons[3])
    add('polygon_ideal_residual',11264*c-2*T-3*polygons[0]-2*polygons[1]-polygons[2])
    minima={3:6,4:10,5:16,6:27,7:46,8:82,9:218,10:1923,11:183968}
    nus={3:9,4:15,5:26,6:45,7:81,8:217,9:1922,10:183967,11:1674688396}
    for k,v in minima.items():add(f'facet_min_{k}',(13-k)*a[k-1]-v*a[k])
    for k,v in nus.items():add(f'nu_min_{k}',(13-k)*a[k-1]+beta(k)*c-v*a[k])
    for k in range(2,7):
        add(f'nk_{k}_0',Q(comb(12,k),comb(6,k))*a[k]-comb(12,k)*f-beta(k)*c)
        for l in range(1,k):add(f'nk_{k}_{l}',Q(comb(6,l),comb(6,k))*a[k]-a[l])
    for label,value in (
        ('weighted',Q(22,5)*a[2]-66*f-Q(2123,1280)*11264*c),
        ('ratio_lower',f-Q(1474,21)*c),('ratio_upper',Q(1298,9)*c-f),
        ('T_faces',T-Q(331,437)*a[2]),('T_finite',T-Q(1986,59)*f),
        ('B_density',B-Q(2515,59)*f),('C_density',C-Q(627,1475)*f),
        ('B_facet_gap',B-7*a[3]+10*a[2]),('B_nu_gap',B-10*a[3]+10*a[2]+beta(3)*c),
        ('C_facet_gap',C-11*a[4]+9*a[3]),('C_nu_gap',C-16*a[4]+9*a[3]+beta(4)*c),
        ('D_facet_gap',D-17*a[5]+8*a[4]),('K_canonical',K-16*D),
        ('incoming_K5',792*f+160*D-11*K),('moment_cap',Q(516992,3125)*f-twoB)):
        add(label,value)
    for k,value,v,w in ((3,B,2,3),(4,C,5,5),(5,D,16,10)):
        add(f'canonical_faces_{k}',a[k]-value)
        add(f'canonical_finite_{k}',comb(12,k)*f-v*value)
        add(f'canonical_ideal_{k}',beta(k)*c-w*value)
    for label,value in (('BT',10*T-6*B),('CB',9*B-10*C),('DC',8*C-16*D),
                        ('CT',45*T-30*C),('DT',120*T-160*D),('DB',36*B-80*D)):
        add(label,value)
    return names,rows,eq,twoB.v


def verify(record):
    names,rows,eq,obj=model();n=len(names)
    x=[Q(record['primal'].get(name,'0')) for name in names]
    multipliers={name:Q(value) for name,value in record['dual_nonnegative_rows'].items()}
    equality={name:Q(value) for name,value in record['dual_equalities'].items()}
    assert all(value>=0 for value in x) and all(value>=0 for value in multipliers.values())
    evaluate=lambda row:sum(a*b for a,b in zip(row[:-1],x))+row[-1]
    assert all(evaluate(row)>=0 for _,row in rows)
    assert all(evaluate(row)==0 for _,row in eq)
    combined=[Q(0)]*(n+1)
    for name,row in rows:
        for j,value in enumerate(row):combined[j]+=multipliers.get(name,0)*value
    for name,row in eq:
        for j,value in enumerate(row):combined[j]+=equality.get(name,0)*value
    reduced=[obj[j]-combined[j] for j in range(n)]
    assert all(value>=0 for value in reduced)
    primal=evaluate(obj);dual=obj[-1]-combined[-1]
    assert primal==dual==Q(record['minimum_two_B_over_f'])
    assert (Q(6620,59)-64)/50==Q(1422,1475)
    assert (Q(6620,59)-64)/(220-64)==Q(237,767)
    assert 36*Q(21230,11)-69120==360
    assert 2*(Q(75,2)+360*Q(3,236))==Q(4965,59)
    assert (Q(4965,59)-49)/45==Q(2074,2655)
    assert (Q(4965,59)-49)/(165-49)==Q(1037,3422)
    return primal,min(reduced),min(evaluate(row) for _,row in rows)


def optimize():
    from scipy.optimize import linprog
    names,rows,eq,obj=model();n=len(names)
    # Column scaling only conditions high-dimensional tiny face counts.
    scales=[Q(1)]*n;scales[4]=Q(1,100);scales[5]=Q(1,10**6);scales[6]=Q(1,10**12)
    normal=lambda row:max([abs(row[j]*scales[j]) for j in range(n)]+[abs(row[-1]),Q(1)])
    norms=[normal(row) for _,row in rows];enorms=[normal(row) for _,row in eq]
    result=linprog([float(obj[j]*scales[j]) for j in range(n)],
        A_ub=[[-float(row[j]*scales[j]/norm) for j in range(n)] for (_,row),norm in zip(rows,norms)],
        b_ub=[float(row[-1]/norm) for (_,row),norm in zip(rows,norms)],
        A_eq=[[float(row[j]*scales[j]/norm) for j in range(n)] for (_,row),norm in zip(eq,enorms)],
        b_eq=[-float(row[-1]/norm) for (_,row),norm in zip(eq,enorms)],
        bounds=[(0,None)]*n,method='highs-ds',options={'time_limit':35,'threads':1,'parallel':False})
    assert result.success,result.message
    snapshot=dict(objective=float(result.fun),primal_scaled=result.x.tolist(),
        inequality_residual=result.ineqlin.residual.tolist(),inequality_marginal=result.ineqlin.marginals.tolist(),
        equality_marginal=result.eqlin.marginals.tolist(),lower_marginal=result.lower.marginals.tolist(),
        status=result.message,one_optimizer_call=True)
    CERT.with_name(CERT.stem+'_numeric.json').write_text(json.dumps(snapshot,indent=2)+'\n')
    return recover(snapshot)


def recover(snapshot):
    names,rows,eq,obj=model();n=len(names)
    scales=[Q(1)]*n;scales[4]=Q(1,100);scales[5]=Q(1,10**6);scales[6]=Q(1,10**12)
    normal=lambda row:max([abs(row[j]*scales[j]) for j in range(n)]+[abs(row[-1]),Q(1)])
    norms=[normal(row) for _,row in rows];enorms=[normal(row) for _,row in eq]
    positive=[j for j,value in enumerate(snapshot['primal_scaled']) if value>1e-7]
    active=[i for i,value in enumerate(snapshot['inequality_residual']) if abs(value)<1e-8]
    chosen=[i for i,value in enumerate(snapshot['inequality_marginal']) if abs(value)>1e-8]
    primal_rows=[[row[j] for j in positive]+[-row[-1]] for _,row in eq]
    primal_rows += [[rows[i][1][j] for j in positive]+[-rows[i][1][-1]] for i in active]
    px=linear_system(primal_rows,len(positive));x=[Q(0)]*n
    for j,value in zip(positive,px):x[j]=value
    zero_columns=[j for j,value in enumerate(snapshot['lower_marginal']) if abs(value)<1e-8]
    dual_rows=[[rows[i][1][j] for i in chosen]+[row[j] for _,row in eq]+[obj[j]] for j in zero_columns]
    approximate=[-snapshot['inequality_marginal'][i]/float(norms[i]) for i in chosen]
    approximate += [value/float(norm) for value,norm in zip(snapshot['equality_marginal'],enorms)]
    dy=linear_system(dual_rows,len(chosen)+len(eq),approximate)
    record=dict(status='EXACT_PRIMAL_DUAL_FOR_LINEAR_RELAXATION',dimension=12,normalization='finite vertices f=1',
        minimum_two_B_over_f=str(sum(obj[j]*x[j] for j in range(n))+obj[-1]),
        primal={name:str(value) for name,value in zip(names,x) if value},
        dual_nonnegative_rows={rows[i][0]:str(value) for i,value in zip(chosen,dy)},
        dual_equalities={name:str(value) for (name,_),value in zip(eq,dy[len(chosen):])},
        variables=n,nonnegative_inequalities=len(rows),equalities=len(eq),
        optimizer_calls=1,optimizer_status=snapshot['status'],
        scope=dict(strict_rows_relaxed_to_nonstrict=True,nonlinear_macaulay_conditions_included=False,
            actual_corner_graph_realization_included=False,clique_vector_hull_77_partitions=True,
            closed_cluster_defect_not_quantified=True,dimension_exclusion=False,lean_compilation_performed=False))
    verify(record);CERT.write_text(json.dumps(record,indent=2)+'\n')
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--optimize',action='store_true')
    parser.add_argument('--recover',action='store_true');args=parser.parse_args()
    if args.optimize:
        lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        R.resources();started=time.monotonic()
        # macOS rejects finite setrlimit values on this host.  Poll the actual
        # peak RSS and wall time instead; this remains one foreground process.
        def guard(signum,frame):
            peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
            if peak>=1024**3:raise MemoryError('One-GiB RSS cap reached')
            if time.monotonic()-started>=59:raise TimeoutError('59-second cap reached')
        signal.signal(signal.SIGALRM,guard);signal.setitimer(signal.ITIMER_REAL,.1,.1)
        record=optimize();R.resources(start=False)
        signal.setitimer(signal.ITIMER_REAL,0)
        peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
        assert peak<1024**3;print('Peak RSS MiB:',round(peak/1024**2,2))
    elif args.recover:
        with CERT.with_name(CERT.stem+'_numeric.json').open() as stream:record=recover(json.load(stream))
    else:
        with CERT.open() as stream:record=json.load(stream)
    value,_,_=verify(record)
    print('PASS exact primal/dual minimum 2B/f =',value,'=',float(value),'; no dimension exclusion')


if __name__=='__main__':main()
