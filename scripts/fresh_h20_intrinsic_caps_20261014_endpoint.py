#!/usr/bin/env python3
"""Complete PSD endpoint profiles for the normalized A11~ x A2~ x E6~ cusp.

The finite pool is derived from Schur/three-node PSD inequalities, not an
arithmetic angle restriction. An isolated I2(m) tail is retained explicitly.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import isqrt, sqrt
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_intrinsic_caps_20261014_endpoint'
SEED='results/fresh_h20_quartic_orbit_cap_gate_20261014.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUT=ROOT/f'results/{TAG}.json'
RAD=[1,2,3,6,5,10,15,30]
SCALE=1<<120
ZERO=(F(0),)*8


def c(x):return (F(x),)+(F(0),)*7
def add(a,b):return tuple(x+y for x,y in zip(a,b,strict=True))
def scale(a,x):return tuple(x*v for v in a)
def mul(a,b):
    out=[F(0)]*8
    for i,x in enumerate(a):
        if not x:continue
        for j,y in enumerate(b):
            if not y:continue
            mask=i&j
            out[i^j]+=x*y*((2 if mask&1 else 1)*(3 if mask&2 else 1)*(5 if mask&4 else 1))
    return tuple(out)
def sign(a):
    if a==ZERO:return 0
    lo=hi=F(0)
    for n,x in zip(RAD,a,strict=True):
        if not x:continue
        b=isqrt(n*SCALE*SCALE);l=F(b,SCALE)
        h=l if b*b==n*SCALE*SCALE else F(b+1,SCALE)
        lo+=x*(l if x>0 else h);hi+=x*(h if x>0 else l)
    assert lo>0 or hi<0,'unresolved radical sign'
    return 1 if lo>0 else -1
def approx(a):return sum(float(x)*sqrt(n)for x,n in zip(a,RAD,strict=True))
def enc(a):return list(map(str,a))
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def inverse(A):
    n=len(A);M=[list(map(F,row))+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
    for j in range(n):
        p=next(i for i in range(j,n)if M[i][j]);M[j],M[p]=M[p],M[j]
        v=M[j][j];M[j]=[x/v for x in M[j]]
        for i in range(n):
            if i==j:continue
            v=M[i][j]
            if v:M[i]=[x-v*y for x,y in zip(M[i],M[j],strict=True)]
    return [row[n:]for row in M]
def pd(A):
    M=[list(map(F,row))for row in A];p=[]
    for i in range(len(M)):
        v=M[i][i];assert v>0;p.append(str(v))
        for j in range(i+1,len(M)):
            for k in range(j,len(M)):
                M[j][k]=M[k][j]=M[j][k]-M[j][i]*M[k][i]/v
    return p
def components(A):
    left=set(range(len(A)));out=[]
    while left:
        seen={min(left)};todo=list(seen)
        for i in todo:
            for j in sorted(left-seen):
                if A[i][j]:seen.add(j);todo.append(j)
        out.append(sorted(seen));left-=seen
    return out


def run():
    packet=json.loads((ROOT/SEED).read_text());seed=packet['seed_cusp']
    for p,h in packet['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
    Q=[[-F(x)for x in row]for row in raw['basis_gram']]
    W=[list(map(F,row))for row in raw['wall_coordinates_as_columns_transposed']]
    C=[[dot(row,w)for row in Q]for w in W]
    G=[[dot(w,cw)for cw in C]for w in W]
    factors=seed['unchanged_affine_factors']
    A=next(x['walls']for x in factors if len(x['walls'])==12)
    E=next(x['walls']for x in factors if len(x['walls'])==7)
    T=[23,50,200];IDS=A+T+E
    assert len(IDS)==len(set(IDS))==22 and sorted(IDS)==seed['seed_walls']
    assert all(G[i][i]==2 for i in A+E)
    H=[[2 if i==j else (-1 if i in T and j in T else G[i][j])for j in IDS]for i in IDS]
    # The active projected triangle has normalized Gram A2~. This source
    # identity is independently checked in the bound seed certificate.
    assert seed['projected_Gram_rank']==19 and seed['projected_span_rank']==20
    assert all(x in[-1,0,2]for row in H for x in row)
    index={w:i for i,w in enumerate(IDS)}
    gram=lambda i,j:H[index[i]][index[j]]
    a_adj={i:[j for j in A if i!=j and gram(i,j)]for i in A}
    assert all(len(v)==2 for v in a_adj.values())
    cycle=[4,min(a_adj[4])]
    while len(cycle)<12:
        cycle.append(next(j for j in a_adj[cycle[-1]]if j!=cycle[-2]))
    assert set(cycle)==set(A) and cycle[0]in a_adj[cycle[-1]]
    branches=[(33,26),(40,34),(45,24)]
    assert {(i,j)for i,j in combinations(E,2)if gram(i,j)}=={
        tuple(sorted(p))for mid,leaf in branches for p in[(22,mid),(mid,leaf)]}
    marks={i:1 for i in A+T}
    marks.update({22:3,33:2,40:2,45:2,26:1,34:1,24:1})
    assert all(sum(marks[j]*gram(i,j)for j in A)==0 for i in A)
    assert all(sum(marks[j]*gram(i,j)for j in E)==0 for i in E)
    assert sum(marks[i]for i in A)==sum(marks[i]for i in E)==12
    assert [sum(marks[i]*W[i][d]for i in A)for d in range(22)]==list(map(F,seed['null_ray']))
    assert [sum(marks[i]*W[i][d]for i in E)for d in range(22)]==list(map(F,seed['null_ray']))
    mags={3:c(1),4:(F(0),F(1),F(0),F(0),F(0),F(0),F(0),F(0)),
          5:(F(1,2),F(0),F(0),F(0),F(1,2),F(0),F(0),F(0)),
          6:(F(0),F(0),F(1),F(0),F(0),F(0),F(0),F(0))}
    cases=[]
    for typ,omit,multiplicity in [('outer',24,108),('middle',33,108),('center',22,36)]:
        omitted=[4,23,omit];ids=[i for i in IDS if i not in omitted]
        B=[[gram(i,j)for j in ids]for i in ids]
        piv=pd(B);I=inverse(B)
        assert all(I[i][j]>=0 for i in range(19)for j in range(19))
        assert all(I[i][i]>=F(1,2)for i in range(19))
        comps=components(B);singleton={comp[0]for comp in comps if len(comp)==1}
        rows=[];counts=Counter();tested=Counter();min_support_counts={}
        def add_profile(inds,labels):
            cc=[ZERO]*19
            for i,m in zip(inds,labels,strict=True):cc[i]=mags[m]
            cost=ZERO
            for i in inds:
                for j in inds:cost=add(cost,scale(mul(cc[i],cc[j]),I[i][j]))
            residual=add(c(2),scale(cost,-1));sg=sign(residual)
            tested[str(len(inds))]+=1
            if sg<0:return
            unit_cos=[scale(x,F(1,2))for x in cc]
            omitted_ratios=[]
            for name,fac,o in [('A11',A,4),('A2',T,23),('E6',E,omit)]:
                b=ZERO
                for i,w in enumerate(ids):
                    if w in fac:b=add(b,scale(unit_cos[i],F(marks[w],marks[o])))
                aval=F(sum(marks[w]for w in fac),marks[o])
                omitted_ratios.append(dict(factor=name,omitted=o,
                    numerator_over_radius=str(aval),numerator_has_triangle_inradius_s=name=='A2',
                    subtract=enc(b),subtract_approx=approx(b)))
            kind='finite'if sg>0 else'ideal_PSD'
            counts[kind]+=1
            rows.append(dict(id=len(rows),contacts=[[ids[i],m]for i,m in zip(inds,labels,strict=True)],
                cosine_vector=[enc(x)for x in unit_cos],cosine_vector_approx=list(map(approx,unit_cos)),
                normalized_spatial_cost=enc(scale(cost,F(1,2))),
                normalized_spatial_cost_approx=approx(cost)/2,
                norm_two_projection_square=enc(cost),norm_two_schur_residual=enc(residual),
                endpoint=kind,omitted_contact_ratios=omitted_ratios))
        add_profile((),())
        for n in range(1,5):
            good=[]
            for inds in combinations(range(19),n):
                lower=sum(I[i][j]for i in inds for j in inds)
                if lower<=2:good.append(inds)
            min_support_counts[str(n)]=len(good)
            for inds in good:
                for labels in product([3,4,5,6],repeat=n):add_profile(inds,labels)
        tails=[]
        for i in sorted(singleton):
            assert I[i][i]==F(1,2)and all(I[i][j]==0 for j in range(19)if j!=i)
            tails.append(dict(wall=ids[i],finite_labels='every integer m >= 7',
                cosine_vector='single entry cos(pi/m), all others zero',
                normalized_spatial_cost='cos(pi/m)^2',norm_two_schur_residual='2 sin(pi/m)^2',
                parallel_limit=dict(label='infinity',contact_magnitude=2,endpoint='ideal_PSD',
                                    normalized_spatial_cost=1)))
        cases.append(dict(id=len(cases),E6_deletion_type=typ,omitted=omitted,edge_support=ids,
            literal_edge_multiplicity=multiplicity,PD19_Gram=[[int(x)for x in r]for r in B],
            inverse_Gram=[[str(x)for x in r]for r in I],LDL_pivots=piv,
            component_walls=[[ids[i]for i in comp]for comp in comps],
            minimum_cost_support_counts=min_support_counts,finite_pool_trials=dict(tested),
            counts=dict(counts),profiles=rows,dihedral_tails=tails))
        print(typ,dict(counts),'tails',len(tails),flush=True)
    # Explicit automorphisms transport these three complete labelled lists
    # to all 252 literal edge supports. No unverified orbit equivalence.
    transports=[]
    for da,dt,de in product(A,T,E):
        typ=next(t for t,o in [('outer',24),('middle',33),('center',22)]
                 if marks[o]==marks[de])
        base=next(x for x in cases if x['E6_deletion_type']==typ)
        shift=cycle.index(da)
        perm={x:cycle[(i+shift)%12]for i,x in enumerate(cycle)}
        shift=T.index(dt);perm.update({x:T[(i+shift)%3]for i,x in enumerate(T)})
        perm[22]=22
        if typ=='center':order=list(range(3))
        else:
            coord=1 if typ=='outer'else 0
            old=next(i for i,b in enumerate(branches)if b[coord]==base['omitted'][2])
            new=next(i for i,b in enumerate(branches)if b[coord]==de)
            rest=[i for i in range(3)if i!=new];order=[]
            for i in range(3):order.append(new if i==old else rest.pop(0))
        for i,j in enumerate(order):
            for x,y in zip(branches[i],branches[j],strict=True):perm[x]=y
        assert set(perm)==set(IDS)and set(perm.values())==set(IDS)
        assert all(gram(i,j)==gram(perm[i],perm[j])for i in IDS for j in IDS)
        assert all(marks[i]==marks[perm[i]]for i in IDS)
        assert [perm[i]for i in base['omitted']]==[da,dt,de]
        transports.append(dict(omitted=[da,dt,de],type=base['id'],
            edge_support_in_profile_order=[perm[i]for i in base['edge_support']],
            seed_permutation=[perm[i]for i in IDS]))
    assert len(transports)==252 and len({tuple(x['omitted'])for x in transports})==252
    assert Counter(x['type']for x in transports)=={0:108,1:108,2:36}
    sources=[SEED,RAW,f'scripts/{TAG}.py']
    return dict(status='PASS_COMPLETE_INTRINSIC_CUSP_EDGE_PSD_CONTACT_PROFILES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sources},
        exact_number_basis=['1','sqrt2','sqrt3','sqrt6','sqrt5','sqrt10','sqrt15','sqrt30'],
        seed_type='affine A11 x affine A2 x affine E6',seed_norm=2,
        source_factor_correction='The literal rank-six seven-wall factor is affine E6, not D6.',
        seed_order=IDS,normalized_seed_Gram=[[int(x)for x in row]for row in H],
        positive_affine_marks=[marks[i]for i in IDS],
        upper_half_space_inradii=['1','s','1'],
        omitted_ratio_formula='distance_to_omitted_wall / R = A/R - B; A for the triangle includes the free positive factor s.',
        seed_compatibility='Each omitted ratio must be cos(pi/m), m>=2, or at least1 (parallel/dotted). The retained ratios are the saved cosines.',
        completeness='Support<=4 follows inverse Gram entrywise nonnegative and diagonal>=1/2. A label>=7 violates a three-wall PSD minor unless its old wall is isolated and it is the sole new contact. Exhaust all labels3,4,5,6 with Schur>=0, retaining the isolated I2 tail and parallel limit separately.',
        types=cases,literal_edge_transports=transports,
        scope='All necessary finite or ideal Coxeter endpoint contact profiles against the nineteen retained walls of each of the252 cusp edges, with three omitted-seed scalar gates. No rational ambient restriction. Profiles are not mutually compatible cap selections or a finite-volume completion.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    z=run()
    if args.verify:assert z==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(z,indent=2)+'\n')
    print(z['status'])
