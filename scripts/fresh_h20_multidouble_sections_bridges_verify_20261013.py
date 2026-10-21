#!/usr/bin/env python3
"""Independent exact audit of the changed cusp quotient directions."""
from copy import deepcopy
from fractions import Fraction as F
import json
import fresh_h20_multidouble_sections_verify_20261013 as M

A=M.A;B=M.B;K=M.K;V=M.V;U=M.U;ROOT=M.ROOT;prior=M.prior
SOURCE=ROOT/'data/enumerations/fresh_h20_multidouble_sections_bridges_20261013.json'
OUT=ROOT/'data/enumerations/fresh_h20_multidouble_sections_bridges_verify_20261013.json'


def original_cusp(pencil,Q,W,G,h):
    # The prior helper's section constructor also checks the exact bridge
    # construction. Supply an unused original pair solely for that interface.
    row=deepcopy(pencil);row['target_pair_one_based']=[1,2];row['target_labels']=[2,3,4,5,6]
    e0=list(map(F,row['base_normal']));u=list(map(F,row['cusp_vector']));q=(A.dot(e0,e0,Q),)
    vals=[A.trim((A.dot(w,e0,Q),A.dot(w,u,Q)))for w in W]
    gr=lambda a,b:A.sub(A.scale(q,G[a][b]),A.mul(vals[a],vals[b]))
    N=A.mul(gr(0,1),gr(0,1));den=A.mul(gr(0,0),gr(1,1))
    equations={2:N,3:A.sub(A.scale(N,4),den),4:A.sub(A.scale(N,2),den),
               5:A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den)),
               6:A.sub(A.scale(N,4),A.scale(den,3))}
    row['target_equations']={str(k):list(map(str,reversed(v)))for k,v in equations.items()}
    sd,he,u,e0,incident=K.section(row,Q,W,G,h)
    multipliers=K.cusp_link(row,Q,W,G,h,sd,u,e0,incident)
    assert A.dot(e0,W[0],Q)!=0
    return sd,u,e0,multipliers


def candidate(c,pencil,Q,W,G,R,H,h,original,multipliers):
    field=K.Field(A.readp(c['polynomial']),*c['interval'])
    assert len(field.p)>2
    e=[field.rem(x)for x in original[0]]
    assert e==list(map(A.readp,c['normal']))
    sd=prior.products(e,Q,R,H,field);q=sd[1]
    a,b=[i-1 for i in c['pair']];aa,ab,bb=sd[3](a,a),sd[3](a,b),sd[3](b,b)
    N=field.rem(A.mul(ab,ab));den=field.rem(A.mul(aa,bb));label=c['label']
    equations={2:ab,3:A.sub(A.scale(N,4),den),4:A.sub(A.scale(N,2),den),
               5:A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den)),
               6:A.sub(A.scale(N,4),A.scale(den,3))}
    assert field.sign(equations[label])==0
    assert field.sign(ab)<=0 and min(field.sign(aa),field.sign(bb))>0
    if label==5:assert field.sign(A.sub(A.scale(N,2),den))>0
    assert all(field.sign(mu)>0 for mu in multipliers)
    proof=c['cusp_center'];assert [field.rem(mu)for mu in multipliers]==list(map(A.readp,proof['cusp_component_null_multiples']))
    C0=list(map(A.readp,proof['unshifted_point']));C=list(map(A.readp,c['strict_center']));u=list(map(F,pencil['cusp_vector']));shift=proof['shift']
    assert isinstance(shift,int) and shift>0
    assert C==[field.rem(A.add(x,A.scale((y,),shift)))for x,y in zip(C0,u)]
    dual0,_,_=U.field_products(C0,Q,[],field)
    assert field.rem(V.sum_polys(A.scale(x,y)for x,y in zip(dual0,u)))==A.trim((-1,))
    assert field.rem(V.sum_polys(A.scale(x,y)for x,y in zip(dual0,h)))==A.Z
    dual,cp,square=U.field_products(C,Q,R,field)
    assert field.sign(square)<0 and all(field.sign(x)<0 for x in cp)
    assert field.sign(V.sum_polys(A.scale(x,y)for x,y in zip(dual,h)))<0
    assert field.rem(V.sum_polys(A.mul(x,y)for x,y in zip(dual,e)))==A.Z
    if c['target_ridge'] is not None:
        assert c['target_ridge']['ridge_walls_one_based']==c['pair']
        B.actual(c['target_ridge'],Q,R,H,h,sd,field)
    assert c['status']=='NON_COXETER_ACTUAL_RIDGE'
    M.bad_ridge(c['actual_forbidden_ridge'],Q,R,H,h,sd,field)
    return c['target_ridge'] is not None


def main():
    Q,W,G,h=A.source();data=json.loads(SOURCE.read_text());assert data['complete'] and data['ambient_chambers']==8
    ambient=json.loads((ROOT/data['ambient_source']).read_text());R=[r[:]for r in W]
    for label in(169,194,71):R,_=M.double(R,Q,W[label-1])
    assert R==ambient['runs'][1]['boundary_normals'] and len(R)==1456
    H=M.Gram(R,Q);cusps=json.loads((ROOT/data['cusp_source']).read_text())
    count=targets=0;cases=[]
    assert [r['option_index_in_prior_census']for r in data['groups']]==[5,6,17,18]
    for group in data['groups']:
        pencil=group['pencil'];assert pencil['cusp_vector']==cusps['cusps'][pencil['cusp_source_index']]['vector']
        sd,u,e0,multipliers=original_cusp(pencil,Q,W,G,h)
        assert all(A.dot(W[j-1],u,Q)<0 for j in(169,194,71))
        assert {tuple(r)for r in R if A.dot(r,u,Q)==0}=={tuple(w)for w in W if A.dot(w,u,Q)==0}
        assert all(A.dot(r,u,Q)<=0 for r in R)
        assert len(group['candidates'])==2
        for c in group['candidates']:
            targets+=candidate(c,pencil,Q,W,G,R,H,h,sd,multipliers);count+=1
            cases.append((c,pencil,sd,multipliers))
        print('verified',pencil['bridge_walls_one_based'],flush=True)
    c,pencil,sd,multipliers=cases[0];bad=deepcopy(c)
    bad['strict_center'][0][-1]=str(F(bad['strict_center'][0][-1])+1)
    try:candidate(bad,pencil,Q,W,G,R,H,h,sd,multipliers)
    except AssertionError:pass
    else:raise AssertionError('Forged center accepted')
    result={'status':'PASS','standard_library_only':True,'new_quotient_directions':4,
            'exact_finite_volume_sections':count,'complete_Coxeter_cusp_rank':19,
            'actual_target_ridges':targets,'actual_forbidden_ridges':count,
            'forged_center_rejected':True,'scope':data['scope']}
    OUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
