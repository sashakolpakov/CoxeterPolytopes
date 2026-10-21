#!/usr/bin/env python3
"""All-radius endpoint-cap obstruction at the golden triangle scale.

Uses the frozen complete endpoint census, not a bounded angle/radius scout.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.fresh_h20_intrinsic_caps_20261014_endpoint import ZERO,c,add,scale,mul,sign,enc,approx
TAG='fresh_h20_intrinsic_caps_20261014_endpoint_bounds'
SOURCE='results/fresh_h20_intrinsic_caps_20261014_endpoint.json'
OUT=ROOT/f'results/{TAG}.json'
PHI=(F(1,2),F(0),F(0),F(0),F(1,2),F(0),F(0),F(0))
SQRT2=(F(0),F(1),F(0),F(0),F(0),F(0),F(0),F(0))
sub=lambda a,b:add(a,scale(b,-1))
square=lambda a:mul(a,a)
read=lambda a:tuple(map(F,a))


def sums(xs):
    out=ZERO
    for x in xs:out=add(out,x)
    return out
def less(r,s):
    return sign(sub(mul(r[0],s[1]),mul(s[0],r[1])))<0
def radius_profile(row,typ):
    ratios=row['omitted_contact_ratios']
    A=[scale(PHI if rr['numerator_has_triangle_inradius_s']else c(1),F(rr['numerator_over_radius']))for rr in ratios]
    B=[read(rr['subtract'])for rr in ratios]
    assert all(sign(a)>0 for a in A)and all(sign(b)>=0 for b in B)
    base=[(a,add(b,c(F(1,2))))for a,b in zip(A,B,strict=True)]
    bound=base[0]
    for r in base[1:]:
        if less(r,bound):bound=r
    accepted=[];rejected=[]
    for i in range(3):
        if B[i]==ZERO:continue
        # At R=A_i/B_i, each other ratio is zero or at least 1/2.
        nums=[sub(mul(A[j],B[i]),mul(B[j],A[i]))for j in range(3)]
        if all(n==ZERO or sign(sub(scale(n,2),A[i]))>=0 for n in nums):
            rr=(A[i],B[i]);accepted.append(i)
            if less(bound,rr):bound=rr
        else:rejected.append(i)
    tightened=False
    if typ==2 and B[1]==B[2]==ZERO:
        # E6 ratio4/R cannot lie in (0,1/2), so R<=8. Then the
        # triangle ratio is >=3phi/8>1/2, forcing >=cos(pi/4).
        assert not less((c(8),c(1)),bound)
        assert sign(sub(scale(PHI,3),c(4)))>0
        golden=(scale(mul(PHI,SQRT2),3),c(1))
        if less(golden,bound):bound=golden;tightened=True
    assert sign(bound[0])>0 and sign(bound[1])>0
    return dict(A=A,B=B,bound=bound,accepted_right_factors=accepted,
                rejected_right_factors=rejected,golden_gap_tightened=tightened)


def run():
    p=json.loads((ROOT/SOURCE).read_text())
    for path,h in p['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert p['status']=='PASS_COMPLETE_INTRINSIC_CUSP_EDGE_PSD_CONTACT_PROFILES'
    ids=p['seed_order'];index={w:i for i,w in enumerate(ids)}
    b=[c(1)]*12+[PHI]*3+[c(1)]*7
    reference_omitted=[4,23,22]
    reference_heights=[c(12),scale(PHI,3),c(4)]
    gaps=[scale(v,F(1,2))for v in b]
    for w,a in zip(reference_omitted,reference_heights,strict=True):
        gaps[index[w]]=add(gaps[index[w]],scale(a,F(1,2)))
    assert all(sign(v)>0 for v in gaps)
    cached=[];radius_rows=[]
    for typ in p['types']:
        I=[[F(x)for x in row]for row in typ['inverse_Gram']]
        rows=[]
        for row in typ['profiles']:
            rr=radius_profile(row,typ['id']);co=list(map(read,row['cosine_vector']))
            cost=read(row['normalized_spatial_cost']);assert sign(sub(c(1),cost))>=0
            rows.append((co,cost,rr['bound']))
            radius_rows.append(dict(type=typ['id'],profile=row['id'],
                bound_numerator=enc(rr['bound'][0]),bound_denominator=enc(rr['bound'][1]),
                accepted_right_factors=rr['accepted_right_factors'],
                rejected_right_factors=rr['rejected_right_factors'],
                golden_gap_tightened=rr['golden_gap_tightened']))
        cached.append((I,rows))
    counts=Counter();edge_rows=[];sample_margins=[];min_float=(float('inf'),None)
    def check_point(de,typ,margin_lower,edge_tag):
        nonlocal min_float
        I,rows=cached[typ]
        V=[sums(scale(v,2*x)for v,x in zip(de,row,strict=True))for row in I]
        norm=sums(mul(x,v)for x,v in zip(de,V,strict=True))
        assert sign(sub(norm,c(margin_lower)))>0
        for j,(co,cost,(A,C))in enumerate(rows):
            t=scale(sums(mul(x,v)for x,v in zip(co,V,strict=True)),-1)
            # C^2 times the power of the point at the largest relaxed ball.
            numerator=add(add(mul(norm,square(C)),scale(mul(mul(A,C),t),2)),
                          mul(sub(cost,c(1)),square(A)))
            assert sign(sub(numerator,scale(square(C),margin_lower)))>0,(edge_tag,j,approx(numerator)/approx(C)**2)
            counts['half_vertex_profile_bounds'if edge_tag is not None else'incenter_profile_bounds']+=1
            if edge_tag is not None:
                ap=approx(numerator)/(approx(C)**2)
                if ap<min_float[0]:min_float=(ap,(list(edge_tag),j,enc(numerator),enc(C)))
                if typ==2 and list(edge_tag)==reference_omitted and j in[1,2,3,5,39]:
                    sample_margins.append(dict(omitted=list(edge_tag),profile=j,
                        positive_power_numerator=enc(numerator),positive_denominator=enc(square(C))))
        return norm
    for tr in p['literal_edge_transports']:
        typ=tr['type'];edge=tr['edge_support_in_profile_order']
        de=[gaps[index[w]]for w in edge]
        norm=check_point(de,typ,F(1,4),tr['omitted'])
        if typ==1:
            assert sign(sub(norm,c(100)))>0
            isolated=p['types'][1]['dihedral_tails'][0]['wall']
            at=p['types'][1]['edge_support'].index(isolated)
            assert de[at]==c(F(1,2))
        edge_rows.append(dict(omitted=tr['omitted'],type=typ,
            point_to_vertex_square=enc(norm),finite_pool_profile_count=len(cached[typ][1])))
    # Incenter is fixed by the transport group, so three representatives suffice.
    for typ in p['types']:
        de=[b[index[w]]for w in typ['edge_support']]
        norm=check_point(de,typ['id'],F(100),None)
        if typ['id']==1:assert sign(sub(norm,c(300)))>0
    # Tail c=cos(pi/m), m>=7 (including parallel limit c=1). It is
    # safe to enlarge to c in[1/2,1]. E6-right R=12/c>=12 is impossible
    # because the untouched triangle gives R<=6phi<12. Thus R<=12/(1+c).
    assert sign(sub(c(12),scale(PHI,6)))>0
    # At the half-vertex, squared distance F>100 and isolated-wall gap1/2:
    # (1+c)^2*power > 232c²+188c−44.
    # Substitute c=(1+t)/2, t in[0,1], obtaining108+210t+58t².
    assert F(232,4)+F(188,2)-44==108
    assert F(232,2)+F(188,2)==210 and F(232,4)==58
    # Dividing by(1+c)^2<=4 gives power>27. The incenter has F>300,
    # gap1: subtract100(1+c)^2, leaving320c²+376c+56>0.
    assert 300-24+144-100==320 and 600-24-200==376 and 300-144-100==56
    assert counts['half_vertex_profile_bounds']==16272
    assert counts['incenter_profile_bounds']==260
    # A point halfway to the reference vertex has quarter of its squared
    # norm. In half-space projective coordinates a formal height square
    # −1/8 gives Lorentz norm+1/8 and still leaves every cap strictly satisfied.
    first=next(r for r in edge_rows if r['omitted']==reference_omitted)
    xnorm=read(first['point_to_vertex_square'])
    expected=(F(76),F(0),F(0),F(0),F(1,2),F(0),F(0),F(0))
    assert xnorm==expected
    paths=[SOURCE,'scripts/fresh_h20_intrinsic_caps_20261014_endpoint.py',f'scripts/{TAG}.py']
    return dict(status='PASS_ALL_INTRINSIC_ENDPOINT_CAPS_MISS_TWO_STRICT_ESCAPES_AT_GOLDEN_SCALE',
        sources={x:sha256((ROOT/x).read_bytes()).hexdigest()for x in paths},
        triangle_inradius='phi=(1+sqrt5)/2',other_inradii=[1,1],
        exact_number_basis=p['exact_number_basis'],radius_relaxation=radius_rows,
        counts=dict(counts),literal_edge_count=252,isolated_tail_edge_count=108,
        edge_distance_squares=edge_rows,critical_five_margins=sample_margins,
        smallest_approximate_power_for_diagnostics_only=min_float[0],
        smallest_row_for_diagnostics_only=list(min_float[1]),
        half_vertex_witness=dict(reference_omitted=reference_omitted,
            point='x=one half of the Euclidean product vertex with these omitted facets',
            all22_seed_gaps=[enc(v)for v in gaps],euclidean_norm_square=enc(xnorm),
            uniform_ball_power_strict_lower_bound='1/4',formal_height_square='-1/8',
            resulting_Lorentz_norm='1/8'),
        incenter_witness=dict(point='x=0',uniform_ball_power_strict_lower_bound=100),
        infinite_tail_proof=dict(cosine_superset=['1/2','1'],radius_bound='12/(1+c)',
            E6_right_radius_rejected='12/c>=12>6phi',half_vertex_distance_square_lower=100,
            half_vertex_isolated_wall_gap='1/2',substituted_power_numerator_lower=[108,210,58],
            resulting_strict_power_lower=27,incenter_distance_square_lower=300,
            incenter_excess_over_power100_polynomial=[56,376,320]),
        scope='For the intrinsic A11~ x A2~ x E6~ seed with inradii(1,phi,1), every acute Coxeter cap that supplies a finite or ideal opposite endpoint of any original cusp edge misses both saved points. All real radii and every integer endpoint label are included, with isolated I2 tails and parallel limits. Any finite completion using only such caps has a strict spacelike escape. Interior caps not incident to any original cusp edge, changed seed facets or other scale ratios are outside scope.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run()
    if args.verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print(out['counts']);print(out['smallest_approximate_power_for_diagnostics_only'])
