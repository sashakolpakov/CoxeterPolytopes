#!/usr/bin/env python3
"""Simultaneous neighbor-seam repairs for four changed Coxeter-cusp pencils.

Construct all subsets of three commuting noncusp doublings after the earlier
eight-chamber union. A persistent actual ridge reduces every real parameter
to a single quadratic candidate, unless the ridge's seam is doubled, when
its replacement is always obtuse. Exact actual ridges reject the candidates.
"""
import hashlib
import json
from itertools import combinations
import numpy as np
import sympy as sp
import fresh_h20_multidouble_sections_bridges_20261013 as B

M=B.M;O=B.O;D=B.D;P=B.P;t=B.t;old=B.old;ROOT=B.ROOT
OUT=ROOT/'data/enumerations/fresh_h20_adapted_sections_20261013.json'


def digest(rows):
    return hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()


def project_face(Q,W,G,h,ids):
    A=sp.Matrix(W[ids].tolist());GG=sp.Matrix(G[np.ix_(ids,ids)].tolist())
    coefficients=GG.inv()*A*sp.Matrix(Q.tolist())*sp.Matrix(list(map(int,h)))
    point=np.array(sp.Matrix(list(map(int,h)))-A.T*coefficients,dtype=object).reshape(-1)
    values=W@Q@point
    assert max(values)==0 and point@Q@point<0 and point@Q@h<0
    assert [i for i,v in enumerate(values)if v==0]==ids
    return point,{'zero_original_walls_one_based':[i+1 for i in ids],
                  'projection_coefficients':list(map(str,coefficients)),
                  'point':list(map(str,point)),
                  'square':str(point@Q@point),'h_pairing':str(point@Q@h)}


def field_pairings(Q,R,point,F):
    # Reduce once, then pair coefficient vectors with the integral normals.
    degree=F.poly.degree() if hasattr(F,'poly') else F.p.degree()
    coeff=np.array([[F.rem(p).nth(j)for j in range(degree)]for p in point],dtype=object)
    values=(R@Q)@coeff
    return [F.rem(P(sum((v*t**j for j,v in enumerate(row)),sp.S.Zero)))for row in values]


def run():
    Q,W,G,h=O.base.load();Q=np.asarray(Q,dtype=np.int64);W=np.asarray(W,dtype=np.int64);G=np.asarray(G,dtype=np.int64)
    source=json.loads(B.OUT.read_text());base=M.build(Q,W,G)[1][8]
    families=[];states={}
    for seams in ((23,49,50),(22,23,24)):
        assert all(G[i-1,j-1]==0 for i,j in combinations(seams,2))
        rows=[]
        for mask in range(8):
            R=base.copy();stages=[];chosen=[seams[k]for k in range(3)if mask>>k&1]
            for wall in chosen:
                R,record=M.double(R,Q,W[wall-1]);record['original_seam_one_based']=wall;stages.append(record)
            states[(seams,mask)]=R
            rows.append({'mask':mask,'seams_one_based':chosen,'stages':stages,
                         'ambient_chambers':8*2**len(chosen),'boundary_count':len(R),
                         'boundary_digest':digest(R.tolist())})
        families.append({'available_seams_one_based':list(seams),'assemblies':rows})
        print('built simultaneous seams',seams,[r['boundary_count']for r in rows],flush=True)
    groups=[]
    for oldgroup,secondary in zip(source['groups'],(45,11,3,5)):
        pencil=oldgroup['pencil'];option=oldgroup['option_index_in_prior_census']
        sd,e0,u=B.section(Q,W,G,h,pencil);e,q,v,gr=sd
        assert q==P(4) and v[0]==P(2)
        seams=(23,49,50)if option<10 else(22,23,24)
        assert all(W[i-1]@Q@u<0 for i in (169,194,71)+seams)
        roots=[];components=[]
        for comp in pencil['central_cusp_link']['components']:
            ids=[i-1 for i in comp['walls_one_based']];marks=comp['positive_kernel_marks']
            mu=P(4*sp.Rational(comp['null_sum_cusp_multiple'])-t*sum(m*(W[i]@Q@e0)for m,i in zip(marks,ids)))
            components.append(O.coeffs(mu))
            if mu.degree()==1:
                end=-mu.nth(0)/mu.nth(1)
                ids=sorted(set(ids+[0,22]));point,record=project_face(Q,W,G,h,ids)
                assert (e0+end*u)@Q@point==0 and point@Q@u<0
                record.update(parameter=str(end),u_pairing=str(point@Q@u))
                roots.append((end,point,record))
        assert len(roots)==2
        roots.sort(key=lambda r:r[0]);(lo,pl,lr),(hi,pr,rr)=roots
        assert hi-lo==1 and v[22]==P(-2*(t-lo))
        for mu in components:
            poly=sp.Poly.from_list(list(map(sp.Rational,mu)),t,domain=sp.QQ)
            assert poly.eval((lo+hi)/2)>0
        assert set(lr['zero_original_walls_one_based'])&set(rr['zero_original_walls_one_based'])=={1,23}
        # Positive weights give a future timelike section ridge for all lo<t<hi.
        ap=-pl@Q@u;bp=-pr@Q@u
        point=[P((hi-t)*bp*x+(t-lo)*ap*y)for x,y in zip(pl,pr)]
        assert old.pdot(Q,point,e,type('Identity',(),{'rem':staticmethod(lambda p:p)})()).is_zero
        x=v[22];equation=P(5*x*x+16*x+8).monic()
        intervals=[(a,b)for(a,b),m in equation.intervals(eps=sp.Rational(1,10**18))if lo<(a+b)/2<hi]
        assert len(intervals)==1
        F=O.Field(equation,*intervals[0]);C,centerproof=B.cusp_center(Q,W,W,h,pencil,sd,F)
        target=D.ridge(Q,W,h,sd,point,F,0,22);assert target is not None
        assert F.rem(4*gr(0,22)**2-gr(0,0)*gr(22,22)).is_zero and F.sign(gr(0,22))<0
        bad=D.ridge(Q,W,h,sd,C,F,0,secondary-1,True)
        assert bad is not None and bad['angle_rejection']=='obtuse'
        badpoint=[sp.Poly.from_list(list(map(sp.Rational,p)),t,domain=sp.QQ)for p in bad['point_coordinate_polynomials']]
        checks=[]
        for mask in range(8):
            R=states[(seams,mask)];lookup={tuple(r):i for i,r in enumerate(R)}
            doubled=23 in [seams[k]for k in range(3)if mask>>k&1]
            primary=[lookup[tuple(W[0])],lookup[tuple(W[0]+W[22]if doubled else W[22])]]
            vl=R@Q@pl;vr=R@Q@pr
            assert max(vl)<=0 and max(vr)<=0
            assert [i for i in range(len(R))if vl[i]==vr[i]==0]==sorted(primary)
            assert {tuple(r)for r in R if r@Q@u==0}=={tuple(w)for w in W if w@Q@u==0}
            residual=[lookup[tuple(W[0])],lookup[tuple(W[secondary-1])]]
            vals=field_pairings(Q,R,badpoint,F)
            assert all(F.sign(v)<=0 for v in vals)
            assert [i for i,v in enumerate(vals)if F.sign(v)==0]==sorted(residual)
            checks.append({'mask':mask,'primary_boundary_pair_one_based':[i+1 for i in primary],
                           'primary_type':'always_obtuse_after_doubling'if doubled else'only_pi_over_3_parameter',
                           'candidate_obtuse_boundary_pair_one_based':[i+1 for i in residual]})
        groups.append({'option_index_in_prior_census':option,'pencil':pencil,
                       'available_seams_one_based':list(seams),'component_multipliers':components,
                       'admissible_parameter_interval':[str(lo),str(hi)],
                       'endpoint_ridge_points':[lr,rr],
                       'all_parameter_ridge_point':[O.coeffs(p)for p in point],
                       'candidate':{'polynomial':O.coeffs(equation),
                           'interval':[str(F.lo),str(F.hi)],'approximate_parameter':float(sum(intervals[0])/2),
                           'normal':[O.coeffs(F.rem(p))for p in e],
                           'strict_center':[O.coeffs(p)for p in C],'cusp_center':centerproof,
                           'actual_pi_over_3_ridge':target,'actual_obtuse_ridge':bad},
                       'assembly_checks':checks})
        print('option',option,'all real parameters rejected in eight assemblies; candidate',float(sum(intervals[0])/2),flush=True)
    data={'schema':'simultaneous-neighbor-doubles-continuous-cusp-sections-v1','complete':True,
          'source':str(B.OUT.relative_to(ROOT)),'base_seams_one_based':[169,194,71],
          'families':families,'groups':groups,
          'status':'NO_COXETER_SECTION_IN_THE_32_DECLARED_CONTINUOUS_FAMILIES',
          'scope':'Four specified cusp pencils, each with every subset of its three orthogonal noncusp short-neighbor seams, after the specified eight-chamber union. All real parameters are covered; no exclusion of other cusp directions, other ambient unions, further doublings, or H20 itself.'}
    OUT.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],flush=True)


if __name__=='__main__':run()
