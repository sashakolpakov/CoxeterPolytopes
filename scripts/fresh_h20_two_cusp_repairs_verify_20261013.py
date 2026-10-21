#!/usr/bin/env python3
"""Independent stdlib/Fraction audit of the new two-parameter sections."""
from fractions import Fraction as F
from copy import deepcopy
from math import gcd,lcm
import json
import fresh_h20_multidouble_sections_verify_20261013 as M
import fresh_h20_new_route_20261013_verify as N

A=M.A;ROOT=M.ROOT
SOURCE=ROOT/'data/enumerations/fresh_h20_two_cusp_repairs_exact_20261013.json'
OUT=ROOT/'data/enumerations/fresh_h20_two_cusp_repairs_verify_20261013.json'


def sumps(xs):return M.V.sum_polys(xs)
def field_dot(x,y,Q,f):
    return f.rem(sumps(A.scale(A.mul(x[i],y[j]),Q[i][j])for i in range(22)for j in range(22)if Q[i][j]))


def target_check(record,Q,W,G,h,sd,field):
    a,b=[i-1 for i in record['pair']];aa,bb,ab=sd[3](a,a),sd[3](b,b),sd[3](a,b)
    if record['actual_ridge'] is not None:
        got=M.B.actual(record['actual_ridge'],Q,W,G,h,sd,field)
        assert got==tuple(map(field.rem,(aa,ab,bb)))
    numerator=A.mul(ab,ab);den=A.mul(aa,bb);label=record['label']
    eq={3:A.sub(A.scale(numerator,4),den),4:A.sub(A.scale(numerator,2),den),
        5:A.add(A.sub(A.scale(A.mul(numerator,numerator),16),A.scale(A.mul(numerator,den),12)),A.mul(den,den)),
        6:A.sub(A.scale(numerator,4),A.scale(den,3))}[label]
    assert field.rem(eq)==A.Z and field.sign(ab)<0
    if label==5:assert field.sign(A.sub(A.scale(numerator,2),den))>0


def no_rational_root_mod(p):
    den=lcm(*(x.denominator for x in p));v=[int(x*den)for x in p];d=gcd(*v);v=[x//d for x in v]
    for prime in(17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101):
        if v[-1]%prime and all(sum(c*pow(z,i,prime)for i,c in enumerate(v))%prime for z in range(prime)):
            return prime
    assert N.C.no_rational_root(p)
    return None


def star_check(record,Q,W,G,h):
    def reflect(v,r):
        d=A.dot(v,r,Q);return tuple(x-d*y for x,y in zip(v,r))
    def apply(columns,v):return tuple(sum(v[j]*columns[j][i]for j in range(22))for i in range(22))
    identity=tuple(tuple(int(i==j)for i in range(22))for j in range(22))
    group=[identity];seen={identity}
    for columns in group:
        for r in(W[0],W[21]):
            nxt=tuple(reflect(c,r)for c in columns)
            if nxt not in seen:group.append(nxt);seen.add(nxt)
        assert len(group)<=6
    normals=set()
    for columns in group:
        for i,w in enumerate(W):
            if i in(0,21):continue
            root=apply(columns,w);delta=[a-b for a,b in zip(root,w)]
            a=F(2*A.dot(delta,W[0],Q)+A.dot(delta,W[21],Q),3)
            b=F(A.dot(delta,W[0],Q)+2*A.dot(delta,W[21],Q),3)
            assert min(a,b)>=0 and delta==[a*x+b*y for x,y in zip(W[0],W[21])]
            normals.add(root)
    assert len(group)==record['group_order']==6 and len(normals)==record['boundary_count']==416
    u=record['cusp_vector'];assert A.dot(u,u,Q)==0
    assert A.dot(W[0],u,Q)==0 and A.dot(W[21],u,Q)==-1
    j=record['neighbor_original_wall']-1
    roots=[W[j],[x+y for x,y in zip(W[j],W[0])]]
    assert roots==record['obtuse_boundary_normals']
    assert [[A.dot(a,b,Q)for b in roots]for a in roots]==record['obtuse_gram']==[[2,1],[1,2]]
    assert all(G[k][i]==0 for k in(0,j)for i in(59,98,124))
    pp,pm=[list(map(F,record[key]))for key in('positive_face_point','negative_face_point')]
    for point,zeros in((pp,{0,j,59,98}),(pm,{0,j,124})):
        assert A.dot(point,point,Q)<0 and A.dot(point,h,Q)<0
        vals=[A.dot(point,w,Q)for w in W]
        assert max(vals)==0 and {i for i,x in enumerate(vals)if x==0}==zeros
    vp={r:A.dot(pp,r,Q)for r in normals};vm={r:A.dot(pm,r,Q)for r in normals}
    assert max(vp.values())<=0 and max(vm.values())<=0
    assert {r for r in normals if vp[r]==vm[r]==0}==set(map(tuple,roots))
    # e.pp=-t*(r125.pp)>0 and e.pm=r60.pm+s*r99.pm<0.
    # Positive interpolation is future timelike, e-orthogonal, and strict on
    # every star wall except these two. Their Gram is unchanged by projection.
    assert A.dot(W[124],pp,Q)<0
    assert A.dot(W[59],pm,Q)<0 and A.dot(W[98],pm,Q)<0


def run():
    data=json.loads(SOURCE.read_text());Q,W,G,h=A.source()
    old=json.loads((ROOT/data['cusp_source']).read_text())['cusps']
    U=[old[i]['vector']for i in data['cusp_indices']]
    assert U==data['cusp_vectors']
    assert [[A.dot(a,b,Q)for b in U]for a in U]==data['cusp_gram']==[[0,-24,-2],[-24,0,-14],[-2,-14,0]]
    for u in U:assert A.dot(u,h,Q)<0 and all(A.dot(w,u,Q)<=0 for w in W)
    common=[i+1 for i,w in enumerate(W)if all(A.dot(w,u,Q)==0 for u in U)]
    assert common==data['common_wall_labels'] and len(common)==17
    count=targets=both=irrational_five=0;degrees={};support_records=[];mutant=None
    for row in data['cases']:
        poly=A.readp(row['field_polynomial']);field=M.K.Field(poly,*row['field_interval'])
        x,y=map(A.readp,row['parameters']);assert field.sign(x)>0 and field.sign(y)>0
        e=[field.rem(A.sub(A.sub((F(7*u0),),A.scale(x,u1)),A.scale(y,12*u2)))for u0,u1,u2 in zip(*U)]
        assert e==list(map(A.readp,row['normal']))
        sd=M.prior.products(e,Q,W,G,field);q=sd[1]
        assert field.sign(q)>0 and q==A.readp(row['normal_square'])
        assert field.rem(A.sub(q,A.scale(A.sub(A.add(x,y),A.mul(x,y)),336)))==A.Z
        C=list(map(A.readp,row['strict_center']));dual,cp,square=M.U.field_products(C,Q,W,field)
        assert field.sign(square)<0 and all(field.sign(p)<0 for p in cp)
        assert field_dot(C,e,Q,field)==A.Z
        assert field.sign(field_dot(C,[(F(z),)for z in h],Q,field))<0
        for target in row['repaired_actual_ridges']:
            target_check(target,Q,W,G,h,sd,field)
            targets+=target['actual_ridge']is not None
        both+=all(t['actual_ridge']is not None for t in row['repaired_actual_ridges'])
        M.bad_ridge(row['actual_forbidden_ridge'],Q,W,G,h,sd,field)
        if mutant is None:mutant=(row,Q,W,G,h,sd,field)
        degree=len(poly)-1;degrees[str(degree)]=degrees.get(str(degree),0)+1
        if degree>1 and all(len(p)<=2 for p in e):
            E=[[p[0]if p else F(0)for p in e],[p[1]if len(p)>1 else F(0)for p in e]]
            assert len(N.C.rref(E)[1])==2
            prime=no_rational_root_mod(poly);test=N.plane_four_support(E,W)
            gram=[[A.dot(a,b,Q)for b in E]for a in E]
            record={'scout_proposal':row['scout_proposal'],'rational_root_exclusion_prime':prime,
                    'normal_plane_gram':[[str(x)for x in r]for r in gram],
                    'normal_plane_determinant':str(gram[0][0]*gram[1][1]-gram[0][1]**2),**test}
            support_records.append(record)
            irrational_five+=test['spanning_walls_one_based']is None
        count+=1;print('verified section',row['scout_proposal'],'degree',degree,flush=True)
    star_check(data['initial_star_diagnostic'],Q,W,G,h)
    row,Q,W,G,h,sd,field=mutant
    forged=deepcopy(row['actual_forbidden_ridge']);p=forged['point_coordinate_polynomials'][0];p[-1]=str(F(p[-1])+1)
    try:M.bad_ridge(forged,Q,W,G,h,sd,field)
    except AssertionError:pass
    else:raise AssertionError('forged ridge point accepted')
    forged=deepcopy(row['repaired_actual_ridges'][0]);forged['label']=3
    try:target_check(forged,Q,W,G,h,sd,field)
    except AssertionError:pass
    else:raise AssertionError('forged target angle accepted')
    out={'status':'PASS_EXACT_TWO_PARAMETER_SECTION_REPLAY','standard_library_only':True,
         'finite_volume_sections':count,'actual_repaired_target_ridges':targets,
         'sections_with_both_actual_targets':both,'actual_forbidden_ridges':count,
         'field_degrees':degrees,'irrational_normals_requiring_at_least_five_original_walls':irrational_five,
         'normal_plane_support_checks':support_records,'initial_six_star_all_positive_parameters_obtuse':True,
         'forged_ridge_point_rejected':True,'forged_target_angle_rejected':True,
         'scope':data['scope']}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':run()
