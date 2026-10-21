#!/usr/bin/env python3
"""Remaining high-RM/high-RS duties over the frozen low-S physical domain.

No new normal bank: exact physical intervals on all continuous high-N
parameters, with all S candidates inherited from the complete low-S proof.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_oblique_low_s_continuous_n_interval_20261018 as A

ROOT=A.ROOT
TAG='fresh_oblique_low_s_high_r_polar_interval_20261018'
V,add,sub,mul,div,neg,sqrt,dot=A.V,A.add,A.sub,A.mul,A.div,A.neg,A.sqrt,A.dot


def restore(g,proof):
    wall=proof['wall'];out={int(j):V(z)for j,z in proof['new_products'].items()}
    assert set(out)==set(g)
    for j,z in out.items():g[j][wall]=z
    g[wall]={**out,wall:V(proof['projected_square'])}


def parents(prior):
    base={int(i):{int(j):V(z)for j,z in row.items()}for i,row in prior['refined_primary']['projected_Gram'].items()}
    result=[]
    for row in prior['records']:
        g={i:r.copy()for i,r in base.items()}
        for step in row['refined_parent_steps']:restore(g,step)
        assert 70 in g
        del g[70]
        for r in g.values():del r[70]
        selected=[s for s in row['S']if s['blocker']is None and s['labels'][2]==2]
        assert len(selected)==12 and sum(s['profile']==0 for s in selected)==7
        ss=[]
        for s in selected:
            gs={i:r.copy()for i,r in g.items()};restore(gs,s['extension'])
            # L/M and X/M are exact defining rights; M has profile0.
            assert gs[69][71]==gs[72][71]==V(0)
            assert gs[69][69]==gs[71][71]==V(2)
            ss.append((s,gs))
        result.append((row,g,ss))
    return result


def n_extension(g,n,a):
    frame=[10,67,69];inv,det=A.P.inverse([[g[i][j]for j in frame]for i in frame])
    assert det[1]<0
    cv=[g[i][64]for i in frame];rhs=[V(-a),V(0),neg(n)]
    sig=sub(g[64][64],dot(cv,inv,cv));rho=sub(V(2),dot(rhs,inv,rhs))
    assert sig[0]>0 and rho[0]>0
    factor=sqrt(div(rho,sig));out={};opposite8=None
    for j in g:
        cj=[g[i][j]for i in frame];linear=dot(rhs,inv,cj)
        residual=mul(factor,sub(g[64][j],dot(cv,inv,cj)))
        out[j]=rhs[frame.index(j)]if j in frame else add(linear,residual)
        if j==8:opposite8=sub(linear,residual)
    assert opposite8[0]>0
    gn={i:{**row,73:out[i]}for i,row in g.items()};gn[73]={**out,73:V(2)}
    return gn,dict(frame=frame,inverse=inv,determinant=det,domain=n,N10_magnitude=a,
        witness=64,witness_residual=sig,norm_remainder=rho,opposite_N_sheet_old8=opposite8,new_products=out)


def high_rm(g):
    frame=[69,73,72];gram=[[g[i][j]for j in frame]for i in frame]
    inv,det=A.P.inverse(gram)
    # A positive diagonal and negative determinant certify signature(2,1).
    assert gram[0][0][0]>0 and det[1]<0
    covector=[g[i][71]for i in frame]
    assert covector[0]==covector[2]==V(0) and g[71][71]==V(2)
    minor=sub(mul(gram[0][0],gram[2][2]),mul(gram[0][2],gram[2][0]))
    # H=(-det G)*(3-2*(M^2-v^T G^-1 v)). No division by4-n^2.
    h=sub(mul(sub(V(3),mul(V(2),g[71][71])),neg(det)),
          mul(V(2),mul(mul(covector[1],covector[1]),minor)))
    assert h[0]>0,('RM_margin',h)
    square=mul(V(2),sub(g[71][71],dot(covector,inv,covector)))
    return dict(frame=frame,Gram=gram,inverse=inv,determinant=det,M_covector=covector,
        middle_cofactor=minor,strict_adjugate_margin=h,RM_square=square)


def high_rs(g):
    assert g[72][72]==V(2) and g[72][71]==V(0)
    frame=[69,71,73];gram=[[g[i][j]for j in frame]for i in frame]
    inv,det=A.P.inverse(gram);assert gram[0][0][0]>0 and det[1]<0
    w=[g[i][8]for i in frame];s=[g[i][72]for i in frame]
    residual8=sub(g[8][8],dot(w,inv,w));assert residual8[0]>0
    old8=sqrt(mul(V(2),residual8));assert old8[0]>0
    product=neg(mul(sqrt(div(V(2),residual8)),sub(g[8][72],dot(w,inv,s))))
    if product[0]>0:kind='positive'
    else:assert product[1]<-2,('RS_unresolved',product);kind='strictly_dotted'
    return dict(frame=frame,Gram=gram,inverse=inv,determinant=det,old8_covector=w,
        S_covector=s,old8_residual_square=residual8,positive_R_sheet_old8=old8,
        negative_R_sheet_SR=product,rejection=kind)


def check(ss,n,a):
    rows=[]
    for s,g in ss:
        gn,np=n_extension(g,n,a)
        rows.append(dict(S_recipe=s['recipe'],S_branch=s['branch'],S_profile=s['profile'],S_labels=s['labels'],
            N=np,high_RM=high_rm(gn),high_RS=high_rs(gn)if s['profile']==0 else None))
    return rows


def run(saved=None):
    prior=json.loads((ROOT/f'results/{A.TAG}.json').read_text())
    for path,h in prior['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert prior['summary']['S_recipes']==84 and prior['summary']['S_branches']==168
    assert prior['summary']['S_strict_rejections']==92 and prior['summary']['S_overdomain_survivors']==76
    assert prior['summary']['unresolved']==0
    rows=[];attempts=0;domains=[]
    for parent,g,ss in parents(prior):
        domains.append(dict(parent=parent['index'],S=[[s['recipe'],s['branch'],s['profile'],s['labels']]for s,_ in ss]))
        for a in(0,1):
            todo=[(A.theta(6)[0],F(2),0)]
            while todo:
                lo,hi,depth=todo.pop();attempts+=1
                try:proofs=check(ss,(lo,hi),a)
                except(AssertionError,ZeroDivisionError):
                    assert depth<12,('unresolved',parent['index'],a,lo,hi)
                    mid=(lo+hi)/2;todo.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
                else:
                    rows.append(dict(index=len(rows),parent=parent['index'],N10_magnitude=a,domain=(lo,hi),depth=depth,S=proofs))
                    print('POLAR_BOX',parent['index'],a,depth,len(rows),flush=True)
    for p in range(2):
        for a in(0,1):
            parts=sorted(r['domain']for r in rows if r['parent']==p and r['N10_magnitude']==a)
            assert parts[0][0]==A.theta(6)[0] and parts[-1][1]==2
            assert all(x[1]==y[0]for x,y in zip(parts,parts[1:]))
    allproofs=[s for r in rows for s in r['S']];rm=[s['high_RM']for s in allproofs];rs=[s['high_RS']for s in allproofs if s['high_RS']]
    assert len(rm)==12*len(rows) and len(rs)==7*len(rows)
    paths={Path(__file__),Path(A.__file__),Path(A.P.__file__),Path(A.P.I.__file__)}
    paths.update(ROOT/path for path in prior['sources'])
    paths.update(ROOT/f'{folder}/{A.TAG}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    result=A.P.I.enc(dict(status='PASS_PHYSICAL_LOW_S_ALL_HIGH_RM_RS_CONTINUOUS_N_POLAR_EXCLUSION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},interval_scale=A.P.I.PREC,
        parameter_domain=[A.theta(6)[0],F(2)],N_families=[[2,2,'m>=6'],[3,2,'m>=6']],selected_S=domains,
        records=rows,exact_deletion_debit=prior['exact_deletion_debit'],deleted_walls=prior['deleted_walls'],topology=prior['topology'],
        summary=dict(parents=2,parent_origins=4,low_S_recipes=84,low_S_branches=168,inherited_strict_S_rejections=92,
            inherited_S_overdomain=76,SM_right_survivors=24,SM_right_profile0_survivors=14,
            high_RM_families=48,high_RS_families=28,attempted_boxes=attempts,covering_boxes=len(rows),
            high_RM_box_proofs=len(rm),high_RS_box_proofs=len(rs),
            minimum_high_RM_adjugate_margin=min(r['strict_adjugate_margin'][0]for r in rm),
            maximum_high_RM_determinant=max(r['determinant'][1]for r in rm),
            maximum_high_RS_determinant=max(r['determinant'][1]for r in rs),
            minimum_positive_R_old8=min(r['positive_R_sheet_old8'][0]for r in rs),
            minimum_positive_SR=min(r['negative_R_sheet_SR'][0]for r in rs if r['rejection']=='positive'),
            maximum_dotted_SR=max(r['negative_R_sheet_SR'][1]for r in rs if r['rejection']=='strictly_dotted'),
            unresolved=0,accepted=0,new_cones=0),
        scope='All high RM or RS finite labels in the complete frozen low-S cube, for all high-N labels on two fixed physical parent chains. Necessary finite duties force the stated right products and profile0 restrictions. No new labels sampled, normal banks or cones; historical scaffold/profile completeness and earlier low-S native exclusions are inherited. No new chamber or volume claim.'))
    if saved:assert result==saved
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    assert args.verify or not path.exists();out=run(json.loads(path.read_text())if args.verify else None)
    if not args.verify:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
