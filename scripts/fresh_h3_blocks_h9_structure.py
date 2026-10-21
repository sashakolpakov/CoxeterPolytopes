#!/usr/bin/env python3
"""Exact duplication checks and dyadic structure for the H3-cubed hosts."""
import json
from fractions import Fraction as F
from itertools import product,permutations
from scripts.fresh_h3_blocks_h9_search import ROOT,G,PHI,O,Z,H3,HI,host,norm_pool,zero_census,inverse,mv,pair,enc,determinant

def tr(a):return tuple(zip(*a))
def mm(a,b):return tuple(tuple(sum((x*y for x,y in zip(row,col)),Z) for col in tr(b)) for row in a)
def eye(n):return tuple(tuple(G(i==j) for j in range(n)) for i in range(n))
def matenc(a):return [enc(r) for r in a]
def refl(v,B):
    cov=mv(B,v);q=pair(v,B,v)
    return tuple(tuple(G(i==j)-2*v[i]*cov[j]/q for j in range(len(v))) for i in range(len(v)))
def mod2(v):return int(v.a)%2,int(v.b)%2
def rmul(x,y):return ((x[0]*y[0]+x[1]*y[1])%2,(x[0]*y[1]+x[1]*y[0]+x[1]*y[1])%2)
def radd(x,y):return ((x[0]+y[0])%2,(x[1]+y[1])%2)
RES=((0,0),(1,0),(0,1),(1,1))
def line(v):
    first=next(x for x in v if x!=(0,0));inverse=next(x for x in RES if rmul(x,first)==(1,0))
    return tuple(rmul(inverse,x) for x in v)
def run():
    A=((O,-PHI/2,Z),(Z,(PHI-1)/2,Z),(Z,O/2,-O))
    assert mm(tr(A),mm(tuple(tuple(G(2*(i==j)) for j in range(3)) for i in range(3)),A))==H3
    assert all(x.is_integral for row in inverse(A) for x in row)
    assert determinant(H3)==2*(2-PHI)
    dual=tuple((PHI-1)*HI[i][0] for i in range(3));assert pair(dual,H3,dual)==G(F(3,2))
    assert all((2*x).is_integral for x in dual) and any(not x.is_integral for x in dual)
    # All projective isotropic lines of the three-dimensional discriminant
    # pairing diag(3/2), over O/(2)=F4.
    lines={line(v) for v in product(RES,repeat=3) if any(x!=(0,0) for x in v)
        and tuple(sum(x[k] for x in v)%2 for k in range(2))==(0,0)}
    assert len(lines)==5
    orbit_list=[];unseen=set(lines)
    while unseen:
        v=next(iter(unseen));orbit={line(tuple(v[i] for i in p)) for p in permutations(range(3))}
        orbit_list.append(sorted(orbit));unseen-=orbit
    assert sorted(map(len,orbit_list))==[2,3]
    # The even time coefficient adds a fourth discriminant coordinate.
    # Normalize its generator to -phi/2 in the time line: its doubled
    # discriminant norm is -phi^3 == 1 mod2, matching the three H3 factors.
    timedual=(-PHI/2,);assert -2*PHI*timedual[0]**2==-(1+2*PHI)/2
    def residue_dot(v,w):
        value=(0,0)
        for a,b in zip(v,w):value=radd(value,rmul(a,b))
        return value
    assert all(residue_dot(v,w)!=(0,0) for v in lines for w in lines if v!=w)
    lines4={line(v) for v in product(RES,repeat=4) if any(x!=(0,0) for x in v)
        and tuple(sum(x[k] for x in v)%2 for k in range(2))==(0,0)}
    planes=set()
    for v,w in product(lines4,repeat=2):
        if v==w or residue_dot(v,w)!=(0,0):continue
        vectors={tuple(radd(rmul(a,x),rmul(b,y)) for x,y in zip(v,w)) for a,b in product(RES,repeat=2)}
        planes.add(tuple(sorted({line(z) for z in vectors if any(x!=(0,0) for x in z)})))
    assert len(lines4)==21 and len(planes)==5
    common=set.intersection(*(set(p) for p in planes));assert common=={((1,0),)*4}
    unseen=set(planes);plane_orbits=[]
    while unseen:
        plane=next(iter(unseen));orbit={tuple(sorted(line((v[0],*(v[1+i] for i in perm))) for v in plane)) for perm in permutations(range(3))}
        plane_orbits.append(len(orbit));unseen-=orbit
    assert sorted(plane_orbits)==[2,3]
    # Standard coordinate model of E8, with Bourbaki/Bugaenko simple roots.
    cols=[tuple(G(F(s,2)) for s in (1,-1,-1,-1,-1,-1,-1,1)),(O,O,Z,Z,Z,Z,Z,Z),(-O,O,Z,Z,Z,Z,Z,Z)]
    for i in range(1,6):cols.append(tuple(-O if k==i else O if k==i+1 else Z for k in range(8)))
    R=tr(tuple(cols));assert abs(determinant(R).a)==1 and determinant(R).b==0
    E=mm(tr(R),R);assert all(E[i][i]==2 for i in range(8))
    R9=tuple(tuple(R[i][j] if i<8 and j<8 else G(i==j) for j in range(9)) for i in range(9))
    A9=tuple(tuple(A[i%3][j%3] if i//3==j//3 else Z for j in range(9)) for i in range(9))
    spatial_e8=mm(inverse(R9),A9)
    rows=[]
    for alpha,time_scale,target_alpha in ((PHI,G(F(1,2)),2*PHI),(2*PHI,O,PHI)):
        B=host(alpha);norms=norm_pool(alpha)
        assert all((2*x).is_integral for row in inverse(B) for x in row)
        diagonal=tuple(tuple(-target_alpha if i==j==0 else O if i==j else Z for j in range(10)) for i in range(10))
        target=tuple(tuple(-target_alpha if i==j==0 else E[i-1][j-1] if 1<=i<=8 and 1<=j<=8 else O if i==j==9 else Z for j in range(10)) for i in range(10))
        Ts=[]
        for spatial in (A9,spatial_e8):
            T=tuple(tuple(time_scale if i==j==0 else spatial[i-1][j-1] if i*j else Z for j in range(10)) for i in range(10));Ts.append(T)
        assert tuple(tuple(2*x for x in row) for row in mm(tr(Ts[0]),mm(diagonal,Ts[0])) )==B
        assert tuple(tuple(2*x for x in row) for row in mm(tr(Ts[1]),mm(target,Ts[1])) )==B
        distinctions=[]
        # First H3 simple reflection fails the transported E8 lattice; the
        # middle H3 simple reflection fails the diagonal integral lattice.
        for name,T,Btarget,j in (('E8_plus_scalar',Ts[1],target,1),('diagonal',Ts[0],diagonal,2)):
            source_root=tuple(G(i==j) for i in range(10));r=mv(T,source_root);reflection=refl(r,Btarget)
            failures=[{'row':i,'column':j,'value':reflection[i][j].to_json()} for i in range(10) for j in range(10) if not reflection[i][j].is_integral]
            assert failures and all(x.is_integral for row in refl(source_root,B) for x in row)
            distinctions.append({'comparison':name,'source_simple_root_index':j,'transported_root':enc(r),'nonintegral_reflection_entries':failures})
        rows.append({'alpha':alpha.to_json(),'gram':matenc(B),'determinant':determinant(B).to_json(),'norm_pool':enc(norms),
            'zero_height_dominant_census':zero_census(alpha),'comparison_target_alpha':target_alpha.to_json(),
            'field_similarity_factor':2,'to_diagonal':matenc(Ts[0]),'to_E8_plus_scalar':matenc(Ts[1]),
            'target_E8_gram':matenc(target),'explicit_group_distinctions':distinctions})
    payload={'status':'EXACT_H3_STRUCTURAL_CHECKS_PASS','H3_gram':matenc(H3),'H3_inverse':matenc(HI),'H3_to_2I3':matenc(A),
        'H3_discriminant_generator':enc(dual),'H3_discriminant_generator_norm':G(F(3,2)).to_json(),
        'spatial_isotropic_discriminant_lines':sorted(lines),'block_permutation_orbit_sizes':sorted(map(len,orbit_list)),
        'even_host_isotropic_lines':sorted(lines4),'even_host_maximal_isotropic_planes':sorted(planes),'even_host_common_isotropic_line':list(common)[0],'even_host_maximal_plane_orbit_sizes':sorted(plane_orbits),
        'hosts':rows,'scope':'No new field quadratic spaces. Reflections differ from old E8/scalar and diagonal groups under the explicit displayed similarities. No abstract conjugacy exclusion or worldwide group novelty is claimed; the previous unique-overclosure duplicate argument does not apply.'}
    path=ROOT/'results/fresh_h3_blocks_h9_structure_20261013.json';path.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':payload['status'],'spatial_isotropic_lines':5,'block_orbits':[2,3],'zero_censuses':[len(r['zero_height_dominant_census']) for r in rows]}))
    return payload
if __name__=='__main__':run()
