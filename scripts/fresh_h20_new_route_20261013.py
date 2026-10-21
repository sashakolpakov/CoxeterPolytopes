#!/usr/bin/env python3
"""Cusp-adapted H20 sections preserving a complete Euclidean Coxeter link.

Cut a bridge of an ambient affine cusp diagram, then translate the normal
by the cusp vector. Local projected angles stay exactly constant. Global
angle equations propose algebraic sections; actual ridge certificates decide
their status. This is a finite list of literal cusp pencils, not all cusps.
"""
from __future__ import annotations
import argparse
from collections import Counter
from itertools import combinations
import json
import numpy as np
import sympy as sp
import fresh_h20_global_interior_sections_20261013 as base
import fresh_h20_triple_angle_targets_20261013 as old
import fresh_h20_bivariate_nonelliptic_20261013 as ridge
import fresh_h20_fourwall_sections_20261013 as four
from fresh_h20_global_section_parameters_20261013 import Field,coeffs

ROOT=base.ROOT;P=old.P;t=old.t
OUTPUT=ROOT/'data/enumerations/fresh_h20_new_route_20261013.json'
CUSP_SOURCE=ROOT/'data/enumerations/h20_gluing_independent_cusps_20261013.json'


def components(M,ids=None):
    remaining=set(range(len(M)) if ids is None else ids);out=[]
    while remaining:
        first=min(remaining);remaining.remove(first);seen={first};todo=[first]
        while todo:
            a=todo.pop();neighbors={b for b in remaining if M[a,b]!=0}
            remaining-=neighbors;seen|=neighbors;todo.extend(neighbors)
        out.append(sorted(seen))
    return out


def positive_kernel(M):
    ns=sp.Matrix(M.tolist()).nullspace()
    if len(ns)!=1:return None
    v=ns[0]
    if all(x<0 for x in v):v=-v
    if not all(x>0 for x in v):return None
    denominator=sp.ilcm(*[x.q for x in v]);return list(map(int,base.primitive(np.asarray(v,dtype=object).reshape(-1)*int(denominator))))


def link(Q,W,G,h,u,e):
    ids=[i for i,p in enumerate(W@Q@u) if p==0];q=e@Q@e;v=W@Q@e
    H=np.array([[sp.Rational(G[i,j])-v[i]*v[j]/q for j in ids] for i in ids],dtype=object)
    live=[k for k in range(len(ids)) if H[k,k]>0]
    if any(H[k,k]<0 for k in range(len(ids))):return None
    H=H[np.ix_(live,live)];ids=[ids[k] for k in live]
    if sp.Matrix(H.tolist()).rank()!=19:return None
    for a,b in combinations(range(len(ids)),2):
        if H[a,b]>0:return None
        ratio=H[a,b]**2/(H[a,a]*H[b,b])
        if ratio<1 and ratio not in (0,sp.Rational(1,4),sp.Rational(1,2),sp.Rational(3,4)):return None
        if ratio>1:return None
    comps=[]
    for js in components(H):
        M=H[np.ix_(js,js)];marks=positive_kernel(M)
        if marks is None:return None
        projected=[W[ids[k]]-v[ids[k]]*e/q for k in js]
        null=sum((c*w for c,w in zip(marks,projected)),np.zeros(22,dtype=object))
        pos=next(i for i in range(22) if u[i]);multiple=null[pos]/u[pos]
        if multiple<=0 or any(null[d]!=multiple*u[d] for d in range(22)):return None
        comps.append({'walls_one_based':[ids[k]+1 for k in js],'gram':[[str(x) for x in row] for row in M],
                      'positive_kernel_marks':marks,'null_sum_cusp_multiple':str(multiple),'rank':len(js)-1})
    return {'active_incident_walls_one_based':[i+1 for i in ids],
            'inactive_incident_walls_one_based':[i+1 for i,p in enumerate(W@Q@u) if p==0 and i not in ids],
            'components':comps,'rank':19}


def pencils(Q,W,G,h):
    out=[];seen=set();source=json.loads(CUSP_SOURCE.read_text());counts=Counter()
    for cusp_index,raw in enumerate(source['cusps']):
        u=np.array(raw['vector'],dtype=object);up=W@Q@u
        if max(up)>0:continue
        assert u@Q@u==0 and h@Q@u<0;counts['ambient_cusps']+=1
        ids=[i for i,p in enumerate(up) if p==0]
        for comp in components(G,ids):
            M=G[np.ix_(comp,comp)];marks=positive_kernel(M)
            if marks is None:continue
            for a,b in combinations(range(len(comp)),2):
                if M[a,b]==0:continue
                cut=M.copy();cut[a,b]=cut[b,a]=0
                if len(components(cut))!=2:continue
                counts['ambient_affine_bridge_cuts']+=1
                values=sp.zeros(len(comp),1);values[a]=marks[b];values[b]=-marks[a]
                solution,parameters=sp.Matrix(M.tolist()).gauss_jordan_solve(values)
                solution=solution.subs({x:0 for x in parameters})
                e=sum((c*W[i] for c,i in zip(solution,comp)),np.zeros(22,dtype=object))
                e-=sp.Rational(h@Q@e,h@Q@u)*u
                assert e@Q@h==0 and e@Q@u==0 and e@Q@e>0
                key=tuple(sp.Matrix([list(e),list(u)]).rref()[0])
                if key in seen:continue
                proof=link(Q,W,G,h,u,e)
                if proof is None:continue
                seen.add(key)
                out.append({'cusp_source_index':cusp_index,'cusp_vector':list(map(int,u)),
                            'ambient_component_walls_one_based':[i+1 for i in comp],
                            'bridge_walls_one_based':[comp[a]+1,comp[b]+1],
                            'base_normal':[str(x) for x in e],'normal_square':str(e@Q@e),
                            'central_cusp_link':proof})
    return out,dict(counts)


def section(Q,W,G,h,row):
    e0=np.asarray(list(map(sp.Rational,row['base_normal'])),dtype=object)
    u=np.asarray(row['cusp_vector'],dtype=object)
    e=[P(e0[d]+t*u[d]) for d in range(22)];q=P(e0@Q@e0)
    vals=[P(w@Q@e0+t*(w@Q@u)) for w in W]
    he=P(t*(h@Q@u))
    def gr(a,b):return q*int(G[a,b])-vals[a]*vals[b]
    C=[P(q*int(h[d])-he*e[d]) for d in range(22)]
    return (e,q,vals,gr),C,he


def choose_target(Q,W,G,h,row):
    sd,C,he=section(Q,W,G,h,row);field=type('Sample',(),{'lo':sp.Rational(0),'hi':sp.Rational(0)})()
    e,q,v,gr=sd;targets=[];seen=set()
    for a,b in old.candidates(Q,W,h,sd,C,field):
        if gr(a,b).eval(0)>=0:continue
        N=gr(a,b)**2;den=gr(a,a)*gr(b,b);ratio=sp.cancel(N.as_expr()/den.as_expr())
        if t not in ratio.free_symbols or ratio in seen:continue
        seen.add(ratio);targets.append((a,b))
        if len(targets)==8:break
    best=None
    Qf=np.asarray(Q,float);Wf=np.asarray(W,float);hf=np.asarray(h,float)
    for a,b in targets:
        N=gr(a,b)**2;den=gr(a,a)*gr(b,b)
        for m,eq in ((3,4*N-den),(4,2*N-den),(5,16*N*N-12*N*den+den*den),(6,4*N-3*den)):
            if eq.is_zero:continue
            for root in np.roots([float(x) for x in eq.all_coeffs()]):
                if abs(root.imag)>1e-7 or abs(root.real)<1e-7:continue
                value=float(root.real);ef=np.array([float(p.eval(value)) for p in e]);qf=float(q.eval(value));vf=Wf@Qf@ef
                H=np.asarray(G,float)-np.outer(vf,vf)/qf
                if min(H[a,a],H[b,b])<1e-8 or H[a,b]>=0:continue
                if abs(H[a,b]**2/(H[a,a]*H[b,b])-np.cos(np.pi/m)**2)>1e-6:continue
                Cf=qf*hf-(hf@Qf@ef)*ef;cp=Wf@Qf@Cf
                if max(cp)>=-1e-8:continue
                if next(four.numerical_ridges(H,cp,[(a,b)]),None) is None:continue
                # Prefer a root near the certified central section.
                key=(abs(value),m,a,b)
                if best is None or key<best:best=key
    return None if best is None else (best[2],best[3])


def exact_candidates(Q,W,G,h,row,target):
    sd,rawC,he=section(Q,W,G,h,row);e,q,v,gr=sd;a,b=target
    N=gr(a,b)**2;den=gr(a,a)*gr(b,b)
    equations={2:N,3:4*N-den,4:2*N-den,5:16*N*N-12*N*den+den*den,6:4*N-3*den}
    factors={}
    for eq in equations.values():
        assert not eq.is_zero
        for f,power in sp.factor_list(eq)[1]:factors[tuple(coeffs(f.monic()))]=f.monic()
    rows=[]
    for key,poly in sorted(factors.items()):
        for (lo,hi),mult in poly.intervals(eps=sp.Rational(1,10**18)):
            if lo==hi:lo-=sp.Rational(1,10**20);hi+=sp.Rational(1,10**20)
            F=Field(poly,lo,hi);c={'minimal_polynomial':list(key),'isolating_interval':[str(lo),str(hi)],'approximate_parameter':float((lo+hi)/2)}
            if F.sign(P(t))==0:c['status']='CENTRAL_PARAMETER';rows.append(c);continue
            if min(F.sign(gr(a,a)),F.sign(gr(b,b)))<=0 or F.sign(gr(a,b))>0:
                c['status']='TARGET_NOT_ACUTE_ELLIPTIC_PAIR';rows.append(c);continue
            labels=[m for m,eq in equations.items() if F.sign(eq)==0 and (m!=5 or F.sign(2*N-den)>0)]
            if not labels:c['status']='CONJUGATE_ANGLE_BRANCH';rows.append(c);continue
            obstruction=None
            for component_index,component in enumerate(row['central_cusp_link']['components']):
                js=[i-1 for i in component['walls_one_based']];marks=component['positive_kernel_marks']
                null=[F.rem(sum((mark*(q*int(W[i,d])-v[i]*e[d]) for mark,i in zip(marks,js)),P(0))) for d in range(22)]
                index=next(d for d,x in enumerate(row['cusp_vector']) if x)
                multiple=F.rem(null[index]/row['cusp_vector'][index])
                assert all(F.rem(null[d]-multiple*row['cusp_vector'][d]).is_zero for d in range(22))
                if F.sign(multiple)<=0:
                    obstruction={'component_index':component_index,'cusp_multiple':coeffs(multiple)};break
            if obstruction is not None:
                c.update(status='NO_INTERIOR_BY_CUSP_CONE_IDENTITY',target_labels=labels,cusp_cone_obstruction=obstruction)
                rows.append(c);continue
            C=[F.rem(p) for p in rawC];cp=[old.pdot(Q,C,[P(int(x)) for x in w],F) for w in W];method='projected_original_center'
            if any(F.sign(p)>=0 for p in cp):
                C=four.alternate_center(Q,W,h,sd,F)
                if C is None:c.update(status='UNRESOLVED_INTERIOR_SCOUT',target_labels=labels);rows.append(c);continue
                cp=[old.pdot(Q,C,[P(int(x)) for x in w],F) for w in W];method='rational_LP_point_projected_exactly'
            target_proof=ridge.actual_ridge(Q,W,h,sd,C,F,a,b,cp);bad=None
            for aa,bb in old.candidates(Q,W,h,sd,C,F):
                bad=old.witness(Q,W,h,sd,C,F,aa,bb,cp)
                if bad is not None:break
            # Positive multiples of u certify all affine link components at this parameter.
            link_proofs=[]
            for component in row['central_cusp_link']['components']:
                js=[i-1 for i in component['walls_one_based']];marks=component['positive_kernel_marks']
                normals=[[F.rem(q*int(W[i,d])-v[i]*e[d]) for d in range(22)] for i in js]
                null=[F.rem(sum((c*n[d] for c,n in zip(marks,normals)),P(0))) for d in range(22)]
                index=next(d for d,x in enumerate(row['cusp_vector']) if x)
                multiple=F.rem(null[index]/row['cusp_vector'][index])
                assert F.sign(multiple)>0 and all(F.rem(null[d]-multiple*row['cusp_vector'][d]).is_zero for d in range(22))
                link_proofs.append(coeffs(multiple))
            c.update(isolating_interval=[str(F.lo),str(F.hi)],target_labels=labels,
                     strict_interior=[coeffs(p) for p in C],interior_method=method,
                     repaired_target_ridge=target_proof,actual_forbidden_ridge=bad,
                     cusp_component_null_multiples=link_proofs,
                     status='NON_COXETER_ACTUAL_RIDGE' if bad else 'OPEN_COXETER_CANDIDATE')
            rows.append(c);print('candidate',row['cusp_source_index'],row['bridge_walls_one_based'],labels,c['approximate_parameter'],c['status'],flush=True)
    return {**row,'target_pair_one_based':[a+1,b+1],'target_labels':[2,3,4,5,6],
            'target_equations':{str(m):coeffs(p) for m,p in equations.items()},'candidates':rows}


def build(limit,replay=False):
    Q,W,G,h=base.load()
    if replay:
        prior=json.loads(OUTPUT.read_text());assert prior['complete'];census=prior['cusp_cut_census'];selected=prior['pencils']
    else:
        options,census=pencils(Q,W,G,h);census['distinct_Coxeter_cusp_normal_planes']=len(options);selected=[]
        for row in options:
            target=choose_target(Q,W,G,h,row)
            print('pencil',row['cusp_source_index'],row['bridge_walls_one_based'],'target',target,flush=True)
            if target is not None:selected.append({**row,'target_pair_one_based':[i+1 for i in target]})
            if len(selected)==limit:break
    rows=[]
    for row in selected:
        target=tuple(i-1 for i in row['target_pair_one_based']);rows.append(exact_candidates(Q,W,G,h,row,target))
        OUTPUT.write_text(json.dumps({'complete':False,'cusp_cut_census':census,'pencils':rows},indent=2)+'\n')
    out={'schema':'cusp-adapted-H20-sections-v1','complete':True,'source':str(base.SOURCE.relative_to(ROOT)),
         'cusp_source':str(CUSP_SOURCE.relative_to(ROOT)),'cusp_cut_census':census,'pencil_count':len(rows),
         'pencils':rows,'statuses':dict(Counter(c['status'] for p in rows for c in p['candidates'])),
         'scope':'Bounded literal affine-cusp bridge cuts with a full rank19 Coxeter cusp link and an exact global angle equation. Complete real roots only for the archived pencils/labels2..6; no all-cusp or all-section exclusion.'}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print('DONE',json.dumps(out['statuses']),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--limit',type=int,default=8);ap.add_argument('--replay-pencils',action='store_true')
    args=ap.parse_args();build(args.limit,args.replay_pencils)
