#!/usr/bin/env python3
"""New ideal Z-corner recipes on the frozen 95 J prefixes; bounded discovery."""
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import numpy as np
import sympy as sp
import fresh_wall66_seam_subdivision_20261018 as C

T=C.T; ROOT=C.ROOT; TAG='fresh_wall66_ideal_junction_pilot_20261018'
LABELS=[2,3,4,5,6,0]  # 0 means an explicitly prescribed parallel product -2.


def theta(m):
    return sp.Integer(2) if m==0 else 2*sp.cos(sp.pi/m)


@lru_cache(None)
def corner(profiles,labels):
    """Exact PSD principal-minor check of A2 plus the three norm-two roots."""
    g=sp.eye(5)*2;g[0,1]=g[1,0]=-1
    for j,pi in enumerate(profiles,2):
        if pi:
            node=0 if pi<=3 else 1;m=pi+2 if pi<=3 else pi-1
            g[node,j]=g[j,node]=-theta(m)
    for (a,b),m in zip(combinations(range(2,5),2),labels):g[a,b]=g[b,a]=-theta(m)
    ev=np.linalg.eigvalsh(np.array(g,dtype=float))
    if ev[0]<-1e-8:return None
    if ev[0]>1e-8:
        assert T.P.corner3(*profiles,*labels)is not None
        return dict(kind='finite')
    # Numerical filtering is only acceleration; acceptance is exact.
    minors=[]
    for n in range(1,6):
        for ids in combinations(range(5),n):
            d=sp.simplify(g.extract(ids,ids).det())
            assert d.is_nonnegative is not None,(profiles,labels,ids,d)
            if not d.is_nonnegative:return None
            minors.append(str(d))
    rank=g.rank();assert rank<5
    return dict(kind='affine_necessary',rank=rank,principal_minors=minors,
                gram=[[str(x)for x in row]for row in g.tolist()])


def run():
    p=ROOT/f'results/{C.TAG}.json';old=json.loads(p.read_text())
    for name,h in old['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h
    inp=ROOT/'results/fresh_wall59_low_bond_companions_20261018.json';prior=json.loads(inp.read_text())
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[]
    base={i:v for i,v in s.old.items()if i not in[59,60,61]}
    ids=[i for i in s.cur['actual_sides']if i not in[59,60,61]]
    debit=T.H.deletion(s.primary,s.cur,[59,60,61])
    assert debit['lost_finite_count']==7 and debit['compact_faces_by_dimension']=={2:6,3:1}
    assert all(8 not in r for r in debit['lost_finite_supports'])
    tests=[];native_caps=[];counts=Counter();blockers=Counter();domain=[];unresolved=[]
    diagnostics=[];diagnostic_index={}
    # The literal finite/ideal triangle domain is fixed by each old X, not by J.
    for ri,rec in enumerate(old['records']):
        pair=rec['selection']['pair'];seed=prior['seeds'][pair['seed']]
        x=seed['banks'][1]['caps'][pair['x']];y=seed['banks'][3]['caps'][pair['y']]
        nd={**base,67:np.array(seed['normal']),68:np.array(x['normal']),69:np.array(y['normal'])}
        types={67:0,68:x['profile'],69:y['profile']};frame=[30,67,68]
        m01=T.G.relation(s.dot(nd[30],s.ps(30,types),nd[67],0))
        m12=T.G.relation(s.dot(nd[67],0,nd[68],x['profile']))
        assert isinstance(m01,int)and isinstance(m12,int)
        recipes=[]
        for ms in product(LABELS,repeat=3):
            a=corner((0,0,0),(m01,ms[0],ms[1]))
            b=corner((0,x['profile'],0),(m12,ms[1],ms[2]))
            if a and b and any(z['kind']=='affine_necessary'for z in(a,b)):
                recipes.append(dict(labels=list(ms),corner30_67_Z=a,corner67_X_Z=b))
        domain.append(dict(prefix=ri,X_profile=x['profile'],old_labels=[m01,m12],recipes=recipes))
        # Solve once per X recipe, then screen separately against every frozen J.
        for di,r in enumerate(recipes):
            ms=r['labels'];rhs=[-float(theta(m)) for m in ms]
            rr,info=T.N.norm_roots([nd[i]for i in frame],rhs,2.)
            row=dict(prefix=ri,recipe=di,norm=info,branches=[]);tests.append(row);counts['norm_equations']+=1
            if info['status'].startswith('unresolved'):
                unresolved.append(dict(prefix=ri,recipe=di,norm=info));continue
            for sign,v in rr:
                counts['norm_branches']+=1
                err=max(abs(float(v@T.M@v)-2),max(abs(float(v@T.M@nd[i])-q)for i,q in zip(frame,rhs)))
                vals,bad,cp=s.native(v,0,nd,types)
                assigned_parallel={i for i,m in zip(frame,ms)if m==0}
                bad=[i for i in bad if not(i in assigned_parallel and abs(vals[i]+2)<2e-7)]
                branch=dict(sign=sign,normal=v.tolist(),norm_error=err,center_product=cp,
                    blockers=[[i,vals[i]]for i in bad],J_tests=[]);row['branches'].append(branch)
                if err>2e-5:
                    unresolved.append(dict(prefix=ri,recipe=di,branch=sign,error=err));continue
                if bad or cp>=-1e-8:
                    blockers[tuple(bad)]+=1;continue
                counts['native_before_J']+=1
                for ji,j in enumerate(rec['J_bank']['caps']):
                    value=s.dot(v,0,np.array(j['normal']),j['profile']);rel=T.G.relation(value)
                    jr=dict(J=ji,native_JZ=value,relation=rel);branch['J_tests'].append(jr);counts['J_tests']+=1
                    if rel in(None,'parallel_or_unresolved'):continue
                    nn={**nd,70:np.array(j['normal']),71:v};tt={**types,70:j['profile'],71:0}
                    # Old separated-PD scorer deliberately does not accept near-affine supports.
                    # Save the complete projective cone and independent support diagnostic;
                    # any affine candidate stays unresolved pending exact geometry.
                    bd=T.S.C.detailed_cone(nn,ids+[67,68,69,70,71])
                    squares=[s.profiles[s.ps(i,tt)]['square']if i in tt or i in s.types else float(s.old[i]@T.M@s.old[i])for i in ids+[67,68,69,70,71]]
                    ic=T.H.incidence_score([nn[i]for i in ids+[67,68,69,70,71]],ids+[67,68,69,70,71],squares)
                    diagnostic=dict(boundary=bd,PD_incidence=ic)
                    key=json.dumps(diagnostic,sort_keys=True)
                    if key not in diagnostic_index:
                        diagnostic_index[key]=len(diagnostics);diagnostics.append(diagnostic)
                    missing=[name for name,g in [('Y/Z',vals[69]),('J/Z',value)]if g < -2-2e-7]
                    native_caps.append(dict(prefix=ri,recipe=di,branch=sign,J=ji,native_products=sorted(vals.items()),
                        native_JZ=value,native_YZ=vals[69],diagnostic=diagnostic_index[key],mandatory_dotted_seams=missing,
                        status='UNPROMOTED_IDEAL_CANDIDATE_REQUIRES_AFFINE_INCIDENCE_CERTIFICATE'))
        print('IDEAL_PREFIX',ri,'J',len(rec['J_bank']['caps']),'recipes',len(recipes),'native',len(native_caps),flush=True)
    paths=[Path(__file__),Path(C.__file__),Path(T.__file__),Path(T.N.__file__),Path(T.G.__file__),Path(T.H.__file__),p,inp]
    return dict(status='BOUNDED_IDEAL_JUNCTION_NECESSARY_DOMAIN_PILOT',sources={str(q.relative_to(ROOT)):sha256(q.read_bytes()).hexdigest()for q in paths},
        exact_old_deletion_debit=debit,old_control=dict(walls=63,finite=199,open_ends=10,finite_volume=False),
        profile=0,label_alphabet=LABELS,parallel_symbol=0,prefix_count=16,frozen_J_count=95,
        topology=C.topology(),domains=domain,tests=tests,counts=dict(counts),
        blockers=[dict(walls=list(k),count=v)for k,v in blockers.most_common()],native_candidates=native_caps,
        complete_cone_diagnostics=diagnostics,unresolved=unresolved,
        scope='Only frozen16 X/Y prefixes and95 J roots; new profile0 Z domain requires at least one exact PSD affine boundary-corner recipe with labels2..6 or explicit parallel. Norm roots and full native/common-center screens are numerical. Any surviving affine cone remains explicitly unpromoted pending exact affine incidence/ambient certification. No all-profile/all-label/all-parameter exclusion; unchanged7/6/1 debit; no new finite-volume example.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'native',len(out['native_candidates']),'unresolved',len(out['unresolved']),flush=True)
