#!/usr/bin/env python3
"""Independent stdlib/Fraction replay of simultaneous adapted H20 sections.

No producer, numpy, sympy, numeric search, or floating-point geometry imports.
Rebuilds every ambient union and validates the complete parameter argument.
"""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from math import lcm
import fresh_h20_multidouble_sections_bridges_verify_20261013 as V
import fresh_h20_triple_verify_20261013 as Sturm

M=V.M;A=V.A;B=V.B;K=V.K;U=V.U;ROOT=V.ROOT
SOURCE=ROOT/'data/enumerations/fresh_h20_adapted_sections_20261013.json'
OUT=ROOT/'data/enumerations/fresh_h20_adapted_sections_verify_20261013.json'


def digest(rows):
    return hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()


def rational_pairings(Q,R,point):
    den=lcm(*(F(x).denominator for x in point))
    integer=[int(x*den)for x in point]
    dual=[sum(Q[i][j]*integer[j]for j in range(22))for i in range(22)]
    return [sum(x*y for x,y in zip(r,dual))for r in R]


def field_pairings(Q,R,point,field):
    dual,_,_=U.field_products(point,Q,[],field)
    degree=len(field.p)-1;columns=[]
    for j in range(degree):
        column=[p[j]if j<len(p)else F(0)for p in dual]
        den=lcm(*(x.denominator for x in column))
        columns.append(([int(x*den)for x in column],den))
    return [A.trim(F(sum(x*y for x,y in zip(r,col)),den)for col,den in columns)for r in R]


def endpoint(row,e0,u,Q,W,h):
    end=F(row['parameter']);point=list(map(F,row['point']))
    ids=[i-1 for i in row['zero_original_walls_one_based']]
    coefficients=list(map(F,row['projection_coefficients']))
    assert len(ids)==len(coefficients)
    assert point==[F(h[d])-sum(c*W[i][d]for c,i in zip(coefficients,ids))for d in range(22)]
    assert A.dot(point,point,Q)==F(row['square'])<0
    assert A.dot(point,h,Q)==F(row['h_pairing'])<0
    assert A.dot(point,u,Q)==F(row['u_pairing'])<0
    assert A.dot(point,[a+end*b for a,b in zip(e0,u)],Q)==0
    vals=rational_pairings(Q,W,point)
    assert max(vals)==0 and [i for i,v in enumerate(vals)if v==0]==ids
    return end,point


def actual_residual(proof,Q,W,G,h,sd,field):
    assert proof['angle_rejection']=='obtuse'
    M.bad_ridge(proof,Q,W,G,h,sd,field)
    assert proof['zero_walls_one_based']==proof['ridge_walls_one_based']
    return list(map(A.readp,proof['point_coordinate_polynomials']))


def main():
    data=json.loads(SOURCE.read_text());assert data['complete']
    assert data['source']=='data/enumerations/fresh_h20_multidouble_sections_bridges_20261013.json'
    prior_data=json.loads((ROOT/data['source']).read_text())
    prior_pencils={g['option_index_in_prior_census']:g['pencil']for g in prior_data['groups']}
    assert data['base_seams_one_based']==[169,194,71]
    Q,W,G,h=A.source();base=[r[:]for r in W]
    for wall in data['base_seams_one_based']:
        base,_=M.double(base,Q,W[wall-1])
    assert len(base)==1456
    states={};ambient_counts={}
    assert [r['available_seams_one_based']for r in data['families']]==[[23,49,50],[22,23,24]]
    for family in data['families']:
        seams=family['available_seams_one_based']
        assert all(G[i-1][j-1]==(2 if i==j else 0)for i in seams for j in seams)
        assert [r['mask']for r in family['assemblies']]==list(range(8))
        for raw in family['assemblies']:
            selected=[seams[k]for k in range(3)if raw['mask']>>k&1]
            assert selected==raw['seams_one_based']
            R=[r[:]for r in base]
            for wall,old in zip(selected,raw['stages']):
                R,record=M.double(R,Q,W[wall-1]);record['original_seam_one_based']=wall
                assert record==old
            assert len(raw['stages'])==len(selected)
            assert raw['ambient_chambers']==8*2**len(selected)
            assert raw['boundary_count']==len(R)and raw['boundary_digest']==digest(R)
            states[(tuple(seams),raw['mask'])]=R
            ambient_counts[(tuple(seams),raw['mask'])]=len(R)
        print('verified ambient family',seams,flush=True)
    assert [g['option_index_in_prior_census']for g in data['groups']]==[5,6,17,18]
    count=0;regression_cases=[]
    for group in data['groups']:
        pencil=group['pencil'];seams=tuple(group['available_seams_one_based'])
        assert pencil==prior_pencils[group['option_index_in_prior_census']]
        sd,u,e0,multipliers=V.original_cusp(pencil,Q,W,G,h)
        e,q,vals,gr=sd[:4]
        assert q==(F(4),)and vals[0]==(F(2),)
        assert multipliers==list(map(A.readp,group['component_multipliers']))
        lo,hi=map(F,group['admissible_parameter_interval'])
        assert hi-lo==1 and vals[22]==A.trim((2*lo,-2))
        assert all(A.dot(W[j-1],u,Q)<0 for j in (169,194,71)+seams)
        # These affine positive relations make the interval necessary, not
        # merely a sufficient range inherited from a particular center.
        left=[];right=[]
        for mu in multipliers:
            assert len(mu)<=2 and Sturm.evalp(mu,(lo+hi)/2)>0
            if len(mu)==1:continue
            root=-mu[0]/mu[1]
            if mu[1]>0:left.append(root)
            else:right.append(root)
        assert max(left)==lo and min(right)==hi
        (l,pl),(r,pr)=[endpoint(raw,e0,u,Q,W,h)for raw in group['endpoint_ridge_points']]
        assert (l,r)==(lo,hi)
        zl=set(group['endpoint_ridge_points'][0]['zero_original_walls_one_based'])
        zr=set(group['endpoint_ridge_points'][1]['zero_original_walls_one_based'])
        assert zl&zr=={1,23}
        ap=-A.dot(pl,u,Q);bp=-A.dot(pr,u,Q)
        point=[A.add(A.scale((hi,-1),bp*x),A.scale((-lo,1),ap*y))for x,y in zip(pl,pr)]
        assert point==list(map(A.readp,group['all_parameter_ridge_point']))
        # Both coefficients are positive for lo<t<hi; future timelike
        # convexity and the exact endpoint signs prove actual ridges there.
        pairing=A.Z
        for i in range(22):
            for j in range(22):
                pairing=A.add(pairing,A.scale(A.mul(point[i],e[j]),Q[i][j]))
        assert pairing==A.Z
        x=vals[22]
        assert gr(0,0)==(F(4),)
        assert gr(0,22)==A.sub((-F(4),),A.scale(x,2))
        assert gr(22,22)==A.sub((F(8),),A.mul(x,x))
        # x ranges exactly over (-2,0). For the undoubled pair:
        # cos²=(x+2)²/(8-x²), strictly between0 and1/2, so only m=3.
        # Its necessary equation is 5x²+16x+8=0.
        polynomial=A.add(A.add(A.scale(A.mul(x,x),5),A.scale(x,16)),(F(8),))
        polynomial=A.scale(polynomial,1/polynomial[-1])
        c=group['candidate']
        assert polynomial==A.readp(c['polynomial'])
        assert Sturm.rootcount(polynomial,lo,hi)==1
        field=K.Field(polynomial,*c['interval'])
        assert lo<field.lo<field.hi<hi
        assert field.sign(A.add(x,(F(2),)))>0 and field.sign(x)<0
        efield=[field.rem(z)for z in e]
        assert efield==list(map(A.readp,c['normal']))
        original=V.prior.products(efield,Q,W,G,field)
        assert all(field.sign(mu)>0 for mu in multipliers)
        # The exact center proves full-dimensional finite-volume sections.
        C=list(map(A.readp,c['strict_center']))
        dual,cp,square=U.field_products(C,Q,W,field)
        assert field.sign(square)<0 and all(field.sign(v)<0 for v in cp)
        assert field.sign(M.V.sum_polys(A.scale(v,z)for v,z in zip(dual,h)))<0
        assert field.rem(M.V.sum_polys(A.mul(v,z)for v,z in zip(dual,efield)))==A.Z
        centerproof=c['cusp_center'];C0=list(map(A.readp,centerproof['unshifted_point']));shift=centerproof['shift']
        assert isinstance(shift,int)and shift>0
        assert C==[field.rem(A.add(v,(shift*z,)))for v,z in zip(C0,u)]
        assert list(map(A.readp,centerproof['cusp_component_null_multiples']))==[field.rem(mu)for mu in multipliers]
        dual0,_,_=U.field_products(C0,Q,[],field)
        assert field.rem(M.V.sum_polys(A.scale(v,z)for v,z in zip(dual0,u)))==(-F(1),)
        assert field.rem(M.V.sum_polys(A.scale(v,z)for v,z in zip(dual0,h)))==A.Z
        aa,ab,bb=B.actual(c['actual_pi_over_3_ridge'],Q,W,G,h,original,field)
        assert c['actual_pi_over_3_ridge']['ridge_walls_one_based']==[1,23]
        assert field.sign(ab)<0 and field.rem(A.sub(A.scale(A.mul(ab,ab),4),A.mul(aa,bb)))==A.Z
        bad=c['actual_obtuse_ridge'];badpoint=actual_residual(bad,Q,W,G,h,original,field)
        pair=[i-1 for i in bad['ridge_walls_one_based']]
        assert pair[0]==0 and all(i+1 not in seams+(169,194,71)for i in pair)
        assert [s['mask']for s in group['assembly_checks']]==list(range(8))
        sign_cache={}
        for check in group['assembly_checks']:
            R=states[(seams,check['mask'])];lookup={tuple(v):i for i,v in enumerate(R)}
            selected=[seams[k]for k in range(3)if check['mask']>>k&1]
            doubled=23 in selected
            reflected=[a+b for a,b in zip(W[0],W[22])]
            primary=[lookup[tuple(W[0])],lookup[tuple(reflected if doubled else W[22])]]
            assert [i+1 for i in primary]==check['primary_boundary_pair_one_based']
            vl=rational_pairings(Q,R,pl);vr=rational_pairings(Q,R,pr)
            assert max(vl)<=0 and max(vr)<=0
            assert [i for i,(a,b)in enumerate(zip(vl,vr))if a==b==0]==sorted(primary)
            if doubled:
                assert check['primary_type']=='always_obtuse_after_doubling'
                assert A.dot(W[0],W[22],Q)==-1 and A.dot(W[22],W[22],Q)==2
                # Projected r1/reflected-r1 Gram (scaled by q=4):
                # [4,-2x; -2x,8-(2+x)^2]. Since -2<x<0,
                # offdiagonal>0 and determinant16-16x-8x²>0.
                assert A.dot(W[0],reflected,Q)==1 and A.dot(reflected,reflected,Q)==2
                vrpoly=A.add(vals[0],vals[22])
                assert A.sub((F(4),),A.mul(vals[0],vrpoly))==A.scale(x,-2)
            else:
                assert check['primary_type']=='only_pi_over_3_parameter'
            residual=[lookup[tuple(W[i])]for i in pair]
            assert [i+1 for i in residual]==check['candidate_obtuse_boundary_pair_one_based']
            signs=[]
            for v in field_pairings(Q,R,badpoint,field):
                if v not in sign_cache:sign_cache[v]=field.sign(v)
                signs.append(sign_cache[v])
            assert max(signs)==0 and [i for i,s in enumerate(signs)if s==0]==sorted(residual)
            inc={tuple(r)for r,p in zip(R,rational_pairings(Q,R,u))if p==0}
            assert inc=={tuple(w)for w in W if A.dot(w,u,Q)==0}
            count+=1
        regression_cases.append((group,e0,u,original,field))
        print('verified all parameters in eight assemblies, option',group['option_index_in_prior_census'],flush=True)
    group,e0,u,original,field=regression_cases[0]
    forged=deepcopy(group['endpoint_ridge_points'][0]);forged['point'][0]=str(F(forged['point'][0])+1)
    try:endpoint(forged,e0,u,Q,W,h)
    except AssertionError:pass
    else:raise AssertionError('Forged endpoint accepted')
    forged=deepcopy(group['candidate']['actual_obtuse_ridge']);forged['point_coordinate_polynomials']=[['0']]*22
    try:actual_residual(forged,Q,W,G,h,original,field)
    except AssertionError:pass
    else:raise AssertionError('Forged actual ridge accepted')
    assert count==32 and data['status']=='NO_COXETER_SECTION_IN_THE_32_DECLARED_CONTINUOUS_FAMILIES'
    result={'status':'PASS','standard_library_only':True,'continuous_cusp_pencils':4,
            'assembly_families_per_pencil':8,'all_parameter_families_excluded':count,
            'declared_ambient_states':16,
            'distinct_boundary_lists':len({digest(r)for r in states.values()}),
            'ambient_chamber_counts':[8,16,32,64],
            'boundary_counts_by_family':[list(ambient_counts[(tuple(f['available_seams_one_based']),i)]for i in range(8))for f in data['families']],
            'quadratic_candidates_with_actual_pi_over_3_repair':4,
            'candidates_with_actual_obtuse_residual':4,
            'forged_endpoint_rejected':True,'forged_actual_ridge_rejected':True,
            'sturm_regressions':Sturm.sturm_regressions(),'scope':data['scope']}
    OUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
