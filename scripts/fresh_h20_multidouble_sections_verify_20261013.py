#!/usr/bin/env python3
"""Standalone Fraction replay of successive ambient doubles and H20 ridges."""
from copy import deepcopy
from fractions import Fraction as F
import json
import fresh_h20_reflected_sections_double_verify_20261013 as D

A=D.A;B=D.B;K=D.K;V=D.V;U=D.U;prior=D.prior;ROOT=D.ROOT
SOURCE=ROOT/'data/enumerations/fresh_h20_multidouble_sections_20261013.json'
OUT=ROOT/'data/enumerations/fresh_h20_multidouble_sections_verify_20261013.json'


class Gram:
    def __init__(self,R,Q):
        self.R=R;self.dual=[[sum(Q[i][j]*r[j] for j in range(22))for i in range(22)]for r in R]
    def __getitem__(self,i):
        outer=self
        class Row:
            def __getitem__(self,j):return sum(x*y for x,y in zip(outer.R[i],outer.dual[j]))
        return Row()


def double(R,Q,seam):
    k=R.index(seam);q=A.dot(seam,seam,Q)
    p=[A.dot(r,seam,Q)for r in R]
    assert q>0 and all(x<=0 for j,x in enumerate(p)if j!=k)
    rows=[];seen=set()
    for reflect in (False,True):
        for j,r in enumerate(R):
            if j==k:continue
            if reflect:
                c=F(-2*p[j],q);assert c>=0
                root=[x+c*y for x,y in zip(r,seam)]
            else:root=r[:]
            assert all(F(x).denominator==1 for x in root)
            key=tuple(map(int,root))
            if key not in seen:seen.add(key);rows.append(list(key))
    return rows,{'seam_index_zero_based':k,'seam_square':q,'before':len(R),'after':len(rows)}


def bad_ridge(proof,Q,R,H,h,sd,field):
    aa,ab,bb=B.actual(proof,Q,R,H,h,sd,field)
    N=field.rem(A.mul(ab,ab));den=field.rem(A.mul(aa,bb))
    if proof['angle_rejection']=='obtuse':assert field.sign(ab)>0
    else:
        assert field.sign(ab)<0 and field.sign(N)>0
        assert field.sign(A.sub(A.scale(den,3),A.scale(N,4)))>0
        for p in (A.sub(A.scale(N,4),den),A.sub(A.scale(N,2),den),
                  A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den))):
            assert field.sign(p)!=0


def row_check(row,case,Q,W,G,R,H,h,field):
    u=list(map(F,case['u']));base=list(map(A.readp,case['base_normal']));tau=A.readp(row['translation'])
    e=[field.rem(A.add(x,A.scale(tau,y)))for x,y in zip(base,u)]
    assert e==list(map(A.readp,row['normal']))
    sd=prior.products(e,Q,R,H,field);assert sd[1]==A.readp(row['normal_square']) and field.sign(sd[1])>0
    C=list(map(A.readp,row['strict_center']));dual,cp,square=U.field_products(C,Q,R,field)
    assert field.sign(square)<0 and all(field.sign(x)<0 for x in cp)
    assert field.sign(V.sum_polys(A.scale(x,y)for x,y in zip(dual,h)))<0
    assert field.rem(V.sum_polys(A.mul(x,y)for x,y in zip(dual,e)))==A.Z
    D.cusp(row,case,Q,W,G,e,field)
    if row['target_ridge'] is not None:
        a=row['orthogonal_new_boundary_index']-1;assert sd[2][a]==A.Z
        aa,ab,bb=B.actual(row['target_ridge'],Q,R,H,h,sd,field);assert ab==A.Z
    assert row['status']=='NON_COXETER_ACTUAL_RIDGE'
    bad_ridge(row['actual_forbidden_ridge'],Q,R,H,h,sd,field)
    repair=row.get('local_doubling_repair')
    if repair is not None:
        seam=repair['seam_normal'];r=W[0];q=A.dot(seam,seam,Q);p=A.dot(r,seam,Q)
        reflected=[x-F(2*p,q)*y for x,y in zip(r,seam)]
        a,b=[i-1 for i in repair['new_pair']];assert R[a]==r and R[b]==reflected
        assert p==-1 and q==2 and A.dot(r,r,Q)==2
        assert field.rem(V.sum_polys(A.scale(e[d],sum(Q[d][j]*r[j]for j in range(22)))for d in range(22)))==A.Z
        assert field.sign(sd[3](a,b))>0
        if repair['actual_obtuse_ridge'] is not None:
            bad_ridge(repair['actual_obtuse_ridge'],Q,R,H,h,sd,field)
    return row['target_ridge'] is not None,repair is not None and repair['actual_obtuse_ridge'] is not None


def main():
    Q,W,G,h=A.source();data=json.loads(SOURCE.read_text());source=json.loads((ROOT/data['source']).read_text());case=source['cases'][0]
    assert data['source_case']==0 and data['seams_one_based']==[169,194,71]
    ids=[168,193,70];assert data['seam_gram']==[[G[i][j]for j in ids]for i in ids]==[[4,-4,-4],[-4,4,-4],[-4,-4,4]]
    # On span(r169,r194), s169*s194=I+N with nonzero N and N²=0.
    N=[[2,-2],[2,-2]];assert all(sum(N[i][k]*N[k][j]for k in range(2))==0 for i in range(2)for j in range(2))
    assert N!=[[0,0],[0,0]]
    u=case['u'];R=[r[:]for r in W];states={}
    for n,(label,record)in enumerate(zip(data['seams_one_based'],data['stages']),1):
        seam=W[label-1];assert A.dot(seam,u,Q)<0
        R,step=double(R,Q,seam);step.update(original_seam_one_based=label,ambient_chambers=2**n)
        assert record==step
        if n>=2:states[2**n]=[r[:]for r in R]
    assert len(data['runs'])==2
    records=data['runs']+data['continuation']
    for record in data['continuation']:
        size=record['ambient_chambers'];assert size in (16,32)
        seam=states[8][395] if size==16 else W[21]
        assert A.dot(seam,u,Q)<0
        R,step=double(R,Q,seam);step.update(seam_normal=seam,ambient_chambers=size)
        assert step==record['stage'];states[size]=[r[:]for r in R]
    assert A.dot(states[8][395],W[21],Q)==0
    field=K.Field(A.readp(data['field_polynomial']),*data['field_interval'])
    assert field.p==A.trim((-2,0,1)) and field.sign(A.T)>0
    targets=bad=repairs=0;counts=[]
    for record in records:
        size=record['ambient_chambers'];R=states[size]
        assert R==record['boundary_normals']
        assert {tuple(r)for r in R if A.dot(r,u,Q)==0}=={tuple(w)for w in W if A.dot(w,u,Q)==0}
        assert all(A.dot(r,u,Q)<=0 for r in R)
        H=Gram(R,Q)
        for row in record['candidates']:
            target,repair=row_check(row,case,Q,W,G,R,H,h,field)
            targets+=target;repairs+=repair;bad+=1
        counts.append({'ambient_chambers':size,'facets':len(R),'section_checks':len(record['candidates'])})
        print('verified',counts[-1],flush=True)
    forged=deepcopy(records[0]['candidates'][0]);p=forged['actual_forbidden_ridge']['point_coordinate_polynomials'][0];p[-1]=str(F(p[-1])+1)
    try:row_check(forged,case,Q,W,G,states[4],Gram(states[4],Q),h,field)
    except AssertionError:pass
    else:raise AssertionError('Forged ridge accepted')
    result={'status':'PASS','standard_library_only':True,'infinite_generating_subgroup':True,
            'ambient_unions':counts,'actual_new_target_right_ridges':targets,
            'actual_forbidden_ridges':bad,'actual_obtuse_local_repairs':repairs,
            'complete_cusp_rank':19,'forged_ridge_rejected':True,'scope':data['scope']}
    OUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
