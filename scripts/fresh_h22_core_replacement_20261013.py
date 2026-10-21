#!/usr/bin/env python3
"""Single-wall replacement at a D20 cusp: exact 813-wall H22 chamber.

Standalone stdlib arithmetic. Reconstruct all 812 parent normals, check their
complete Gram gate, delete wall 1, and insert the two prescribed radical caps.
Every radical comparison is reduced to a rational sign/square identity.
No producer or project arithmetic is imported.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path

P = Path(__file__).resolve().parents[1]
PRE = 'fresh_h22_core_replacement_20261013'
OUT = 'results/' + PRE + '.json'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
CORE = 'results/fresh_h22_incidence_core_20261013.json'
BANK = 'data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json'
WALK = 'results/fresh_h21_source_corner_walk_20261013.json'

def read(path): return json.loads((P/path).read_text())
def sha(path): return hashlib.sha256((P/path).read_bytes()).hexdigest()
def dot(a,b): return sum(x*y for x,y in zip(a,b,strict=True))
def mv(M,v): return [dot(row,v) for row in M]
def pair(M,a,b): return dot(a,mv(M,b))
def add(a,b): return [x+y for x,y in zip(a,b,strict=True)]
def scale(a,k): return [k*x for x in a]
def enc(a): return [str(x) for x in a]

def solve(M,b):
    a=[list(map(F,row))+[F(x)] for row,x in zip(M,b,strict=True)]
    for j in range(len(a)):
        i=next(i for i in range(j,len(a)) if a[i][j])
        a[i],a[j]=a[j],a[i]
        q=a[j][j]; a[j]=[x/q for x in a[j]]
        for i in range(len(a)):
            if i!=j and a[i][j]:
                q=a[i][j]; a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return [row[-1] for row in a]

def positive_pivots(M):
    a=[list(map(F,row)) for row in M]; pivots=[]
    for j in range(len(a)):
        q=a[j][j]; assert q>0; pivots.append(q)
        for i in range(j+1,len(a)):
            for k in range(j+1,len(a)):
                a[i][k]-=a[i][j]*a[j][k]/q
    return pivots

def parent_dotted(A,U):
    # A/sqrt(U)+sqrt(2)<-2 iff -A>(2+sqrt2)sqrt(U).
    # Both sides are positive before each squaring.
    B=A*A-6*U
    return A<0 and B>0 and B*B>32*U*U

def setup():
    raw=read(RAW)['coxeter_input']
    Q=[[-int(x) for x in row] for row in raw['basis_gram']]
    W=[[int(x) for x in row] for row in raw['wall_coordinates_as_columns_transposed']]
    core=read(CORE); h=core['center']; points=core['incidence_parts'][0]
    vp=core['old_projections'][1]; C=[]
    for row in read(BANK)['candidates']:
        if row['twice_height_square']!=4: continue
        A=row['A']; mask=row['mask']; bits=[mask>>i&1 for i in range(21)]
        assert sum(bits)==row['weight'] and A*A-sum(bits)+4==4
        twice=[A*vp[j]-sum(bits[i]*W[points[i]][j] for i in range(21)) for j in range(22)]
        assert twice==row['twice_old_projection'] and all(x%2==0 for x in twice)
        C.append([x//2 for x in twice])
    assert len(C)==301 and len(set(map(tuple,C)))==301
    M=[row+[0] for row in Q]+[[0]*22+[2]]
    R=[w+[0] for w in W]+[c+[1] for c in C]+[c+[-1] for c in C]
    return Q,W,C,h,M,R

def run():
    Q,W,C,h,M,R=setup(); hh=h+[0]
    assert len(W)==210 and pair(Q,h,h)==-14
    cov=[mv(M,r) for r in R]; norms=[dot(r,c) for r,c in zip(R,cov)]
    assert Counter(norms)=={2:644,4:168}
    assert all(dot(hh,c)<0 for c in cov)
    base_counts=Counter()
    for i in range(812):
        for j in range(i):
            g=dot(R[i],cov[j]); assert g<=0
            q=F(g*g,norms[i]*norms[j])
            assert q in (0,F(1,4),F(1,2),F(3,4)) or q>=1
            base_counts['dotted' if q>1 else 'parallel' if q==1 else str(q)]+=1
    row=read(WALK)['corner_orbits'][14]
    assert row['id']==14 and row['type']=='ideal'
    delta=row['primitive_ray']; assert pair(Q,delta,delta)==0
    ds=[pair(Q,delta,w) for w in W]
    active=[i for i,x in enumerate(ds) if x==0]
    assert active==row['active_walls'] and len(active)==21 and max(ds)==0
    assert pair(Q,delta,h)==-46
    comps=row['affine_components']; assert len(comps)==1
    weights=list(map(F,comps[0]['positive_kernel']))
    assert comps[0]['walls']==active and all(x>0 for x in weights)
    relation=[sum(a*W[i][j] for a,i in zip(weights,active)) for j in range(22)]
    assert relation==scale(delta,F(comps[0]['ray_multiple']))
    assert all(pair(Q,delta,c)<0 for c in C)
    deleted=1; I=[i for i in active if i!=deleted]; assert len(I)==20
    GI=[[pair(Q,W[i],W[j]) for j in I] for i in I]
    pivots=positive_pivots(GI)
    coeff=solve(GI,[pair(Q,h,W[i]) for i in I])
    p=[F(h[j])-sum(a*W[i][j] for a,i in zip(coeff,I)) for j in range(22)]
    assert all(pair(Q,p,W[i])==0 for i in I)
    A=pair(Q,p,p); d=pair(Q,p,delta)
    assert A==-2484 and d==-46 and pair(Q,p,h)==A
    t=F(350); u=add(scale(p,-1),scale(delta,t)); U=pair(Q,u,u)
    assert U==29716 and U==A-2*t*d and pair(Q,u,h)<0
    v=add(p,scale(delta,t-A/d))
    assert pair(Q,v,v)==-U and pair(Q,v,u)==0 and pair(Q,v,h)==t*d<0
    assert all(pair(Q,u,W[i])==pair(Q,v,W[i])==0 for i in I)
    # PD20 plus its mutually orthogonal positive u and negative v proves
    # Q has signature(21,1), hence diag(Q,2) has signature(22,1).
    old_pair_rows=[]
    for i,w in enumerate(W):
        if i==deleted: continue
        a=pair(Q,u,w); n=pair(Q,w,w)
        if i in I: assert a==0; label='2'
        else: assert a<0 and a*a>2*n*U; label='dotted'
        old_pair_rows.append({'old_wall':i,'unnormalized_pairing':str(a),'label':label})
    parent_A=[pair(Q,u,c) for c in C]
    assert all(parent_dotted(a,U) for a in parent_A)
    assert pair(Q,u,delta)==46>0
    # Each new cap has square U/U + 2/2 =2; their pairing is1-1=0.
    retained=[i for i in range(812) if i!=deleted]
    vp=[dot(v+[0],cov[i]) for i in retained]
    assert {i for i,x in zip(retained,vp) if x==0}==set(I)
    assert all(x<=0 for x in vp)
    # Thus(v,0) is an actual finite vertex, with exactly the20 oldI walls
    # and the two caps incident. Its Gram is GI direct-sum 2I2, PD22.
    delta0=read(WALK)['corner_orbits'][0]['primitive_ray']
    escape=add(scale(delta0,100),h)+[40]
    eq=pair(M,escape,escape); eh=pair(M,escape,hh)
    assert eq==1586 and eh==-814
    old_escape=[dot(escape,cov[i]) for i in retained]; assert max(old_escape)<0
    ex=pair(Q,u,escape[:22]); z=escape[-1]
    assert ex<0 and ex*ex>2*z*z*U
    # Therefore ex/sqrt(U) +/- sqrt2*z is strictly negative for both caps.
    # A strict future spacelike point certifies an open ideal escape patch.
    shorts=[i for i,n in enumerate(norms[:210]) if n==2 and i!=deleted]
    shortsum=[sum(W[i][j] for i in shorts) for j in range(22)]
    negative=pair(Q,shortsum,shortsum); assert negative<0
    b=parent_A[0]; assert b!=0
    radical=2*U; assert radical.denominator==1
    assert isqrt(radical.numerator)**2!=radical.numerator
    # Norm-two pair cyclic product is(pair)^2 = b²/U+2+2b/U sqrt(2U).
    cyc=[b*b/U+2,2*b/U]
    return {
        'status':'PASS_EXACT_813_FACET_CHANGED_CORE_WITH_FINITE_VERTEX_AND_STRICT_ESCAPE',
        'sources':{x:sha(x) for x in (RAW,CORE,BANK,WALK)},
        'coordinate_convention':'H22 form diag(Q,2); last coordinate multiplies physical sqrt2',
        'base_pair_counts':dict(base_counts),'base_facets':812,
        'retained_base_indices':retained,'deleted_old_wall':deleted,'cusp_orbit':14,
        'cusp_ray':delta,'affine_kernel':enc(weights),'twenty_retained_incident':I,
        'positive_gram_pivots':enc(pivots),'projection_coefficients':enc(coeff),
        'projection_p':enc(p),'p_square':str(A),'p_delta':str(d),'parameter_t':str(t),
        'old_projection_numerator_u':enc(u),'u_square':str(U),
        'new_normal_formula':'r_plus=(u/sqrt(U),1/sqrt2), r_minus=(u/sqrt(U),-1/sqrt2)',
        'new_norms':['2','2'],'new_mutual_pairing':'0',
        'new_center_numerator':str(pair(Q,u,h)),'old_cusp_cut_numerator':'46',
        'old_wall_relations':old_pair_rows,'parent_cap_unnormalized_pairings':enc(parent_A),
        'new_pair_counts':{'right':41,'dotted':1582},'facets':813,
        'finite_vertex':{'old_coordinates':enc(v),'last_coordinate':'0','square':str(-U),
            'center_pairing':str(pair(Q,v,h)),'old_incident':I,'new_incident':[0,1],
            'max_nonincident_pairing':str(max(x for x in vp if x<0)),'signature':[22,0,0]},
        'strict_escape':{'coordinates':escape,'square':str(eq),'center_pairing':str(eh),
            'max_retained_pairing':str(max(old_escape)),'new_old_part_numerator':str(ex),
            'last_coordinate':z,'strict_squared_margin':str(ex*ex-2*z*z*U)},
        'conditional_nonquasi_arithmetic':{'cycle_parent_index':210,'radicand':str(radical),
            'cycle_coefficients_in_1_sqrt_radicand':enc(cyc),
            'rational_indefinite_short_block':shorts,'ones_square':str(negative)},
        'scope':'Replaces one of210 vertical old walls in the812 branch, not the860 branch. Full acute Coxeter-pair gate, strict common center and actual finite vertex. Spacelike old projections cut the specified inherited cusp. Exact strict escape proves infinite volume. A cofinite completion retaining this facet core would be non-quasi-arithmetic. No finite-volume example, exhaustive changed-core search or literature novelty claim.'
    }

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--verify',action='store_true'); args=ap.parse_args()
    output=run()
    if args.verify: assert output==read(OUT)
    else: (P/OUT).write_text(json.dumps(output,indent=2)+'\n')
    print(output['status'])
