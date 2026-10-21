#!/usr/bin/env python3
"""Independent rank-three profile/graph proof and rational51-wall replay."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product,permutations
from pathlib import Path
import sympy as s
import fresh_h20_a2_h18_coupling_20261014_review as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_h17_rank_three_coupling_20261015_review'
PRIMARY='results/fresh_h20_h17_rank_three_coupling_20261015.json'
P19='data/known_examples/k3_high_rank_nef_chambers.json'
P48='data/enumerations/coble_type_mii_full_ns_h17.json'
K,Z,O=A.K,A.Z,A.O

def dot(v,w):return sum((a*b for a,b in zip(v[2:],w[2:])),F(0))-v[0]*w[1]-v[1]*w[0]
def gram(N):return [[dot(v,w)for w in N]for v in N]
def vec(v):return list(map(F,v))
def comp(G):
    remaining=set(range(len(G)));out=[]
    while remaining:
        current={min(remaining)}
        while True:
            nxt=current|{j for i in current for j in remaining if G[i][j]<0 and 4*G[i][j]**2==G[i][i]*G[j][j]}
            if nxt==current:break
            current=nxt
        remaining-=current;out.append(sorted(current))
    return out
def label(g,q):
    assert g<=0
    if not g:return'right'
    c=4*g*g/q
    if c in(F(1),F(2),F(3)):return{F(1):'pi3',F(2):'pi4',F(3):'pi6'}[c]
    assert c>=4
    return'parallel'if c==4 else'dotted'
def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,list):return[encode(v)for v in x]
    if isinstance(x,dict):return{k:encode(v)for k,v in x.items()}
    return x

def run():
    d=A.read(PRIMARY)
    raw19,raw48=[json.loads((ROOT/p).read_text())for p in(P19,P48)]
    source19=next(row for row in raw19['cases']if row['rho']==18)
    G19=[[2*int(i==j)for j in range(19)]for i in range(19)]
    for value,edges in source19['intersections'].items():
        for i,j in edges:G19[i-1][j-1]=G19[j-1][i-1]=-int(value)
    assert G19==d['H17_19']['gram']and comp(G19)==[list(range(19))]
    assert all(not G19[i][j]for i,j in combinations([0,2,5,7],2))
    N=list(map(vec,raw48['facets']['physical_NS_wall_vectors']));G=gram(N)
    assert len(N)==48 and all(len(v)==18 for v in N)
    assert G==[[-F(x)for x in row]for row in raw48['facets']['physical_wall_gram']]
    assert s.Matrix(N).rank()==18
    comps=comp(G);assert list(map(len,comps))==[24,12,12]
    assert comps==[row['vertices']for row in d['H17_48']['components']]
    assert all(not G[i][j]for i,j in combinations([0,1,2,3],2))
    groups=[]
    for component in comps[1:]:
        remaining=set(component);parts=[]
        while remaining:
            i=min(remaining);part=sorted(j for j in remaining if j==i or G[i][j]==-2)
            assert len(part)==4;parts.append(part);remaining-=set(part)
        assert len(parts)==3
        for i,j in combinations(component,2):
            assert G[i][j]==(-2 if any(i in p and j in p for p in parts)else-1)
        groups.append(parts)
        for i in component:
            anchors=[j for j in range(8,24)if G[i][j]]
            assert len(anchors)==4 and G[i][i]==2
            assert all(G[j][j]==1 and G[i][j]==-1 for j in anchors)
    assert groups==d['H17_48']['parallel_cosets']
    crossrights=[(i,j)for i in comps[1]for j in comps[2]if not G[i][j]]
    assert crossrights and all(G[i][j]in(0,-2)for i in comps[1]for j in comps[2])
    # Finite connected rank3 diagrams are paths. The determinant budget
    # excludes4-4 and3-6, leaving3-3,3-4,3-5; a triangle is not spherical.
    assert F(1,2)+F(1,2)==1 and F(1,4)+F(3,4)==1
    rt2,phi=A.rt2,A.phi
    support={
        'A3':[[2*O,-O,Z],[-O,2*O,-O],[Z,-O,2*O]],
        'B3':[[2*O,-O,Z],[-O,2*O,-rt2],[Z,-rt2,2*O]],
        'H3':[[2*O,-phi,Z],[-phi,2*O,-O],[Z,-O,2*O]],
        'A2+A1':[[2*O,-O,Z],[-O,2*O,Z],[Z,Z,2*O]],
        'I2(4)+A1':[[2*O,-rt2,Z],[-rt2,2*O,Z],[Z,Z,2*O]],
        'A1^3':[[2*O,Z,Z],[Z,2*O,Z],[Z,Z,2*O]]}
    tables={}
    for name,B in support.items():
        J=A.inverse(B);assert all(A.sign(x)>=0 for row in J for x in row)
        assert all(A.sign(J[i][i]-O/2)>=0 for i in range(3))
        # Labels>=5 cannot have cost1 because phi²/2>1.
        assert A.sign(phi**2/2-O)>0
        profiles=[]
        for c in product((Z,O,rt2),repeat=3):
            if A.dot(c,A.mv(J,c))==O:profiles.append(list(c))
        saved=d['rank_three_cost_one'][name]
        saved_profiles=[[K.from_sympy(s.sympify(x))for x in row]for row in saved['cost_one_profiles']]
        assert set(map(tuple,profiles))==set(map(tuple,saved_profiles))
        gram_profiles=[[A.dot(c,A.mv(J,e))for e in saved_profiles]for c in saved_profiles]
        assert gram_profiles==[[K.from_sympy(s.sympify(x))for x in row]for row in saved['pair_products']]
        half=[[i,j]for i,j in combinations(range(len(saved_profiles)),2)if gram_profiles[i][j]==O/2]
        assert half==saved['half_product_edges']
        assert half==([[3,4],[3,5],[4,5]]if name=='A1^3'else[])
        tables[name]=dict(profile_count=len(profiles),inverse_diagonal=[A.enc(J[i][i])for i in range(3)],half_product_edges=half)
    # The remaining disconnected I2(m)+A1 tail has I2 inverse diagonals
    # >=1+1/sqrt5>1 for m>=5, so only its isolated A1 pi4 survives.
    tail=O/(2*((5*O-A.rt5)/8))
    assert tail==O+O/A.rt5 and A.sign(tail-O)>0
    # On each complete tripartite12-component, only the three doublepi3
    # profiles have half-product edges. Parallel-within-part forbids
    # different double profiles. The two parity components cannot both
    # be special, because all double-profile mutual products are positive.
    DP=[[O,O,Z],[O,Z,O],[Z,O,O]];J3=support['A1^3'];J3=A.inverse(J3)
    dots=[[A.dot(v,A.mv(J3,w))for w in DP]for v in DP]
    assert dots==[[O if i==j else O/2 for j in range(3)]for i in range(3)]
    rep=d['nonordinary_representative'];NN=[v+[F(0)]*3 for v in N]
    for p,color in zip(groups[0],[(0,1),(0,2),(1,2)]):
        for i in p:
            for c in color:NN[i][18+c]=-1
    for i in range(3):v=[F(0)]*21;v[18+i]=1;NN.append(v)
    assert NN==list(map(vec,rep['normals']))and len(NN)==51
    GG=gram(NN);assert GG==[vec(row)for row in rep['gram']]
    hist=Counter(label(GG[i][j],GG[i][i]*GG[j][j])for i,j in combinations(range(51),2))
    assert dict(hist)==rep['pair_histogram']==dict(right=985,pi3=152,pi4=120,parallel=18)
    H=vec(rep['strict_center']);E=vec(rep['strict_escape']);X=vec(rep['boundary_spacelike'])
    assert H==vec([5,5]+[-1]*19)and E==[x+h/F(1000)for x,h in zip(X,H)]
    assert dot(H,H)==-31 and dot(X,X)==7 and dot(E,E)==F(6949969,1000000)
    assert dot(E,H)<0 and all(dot(H,v)<0 and dot(E,v)<0 and dot(X,v)<=0 for v in NN)
    assert [dot(E,v)for v in NN]==vec(rep['escape_pairings'])
    for i,n in enumerate(NN):
        P=[h-dot(H,n)/dot(n,n)*x for h,x in zip(H,n)]
        assert P==vec(rep['facet_points'][i])and dot(P,P)<0 and dot(P,H)<0
        assert all(dot(P,m)==0 if i==j else dot(P,m)<0 for j,m in enumerate(NN))
    scaled=[list(map(F,v))for v in NN]
    for i in list(range(8))+comps[1]:scaled[i]=[x/2 for x in scaled[i]]
    assert scaled==list(map(vec,rep['scaled_integral_generators']))and s.Matrix(scaled).rank()==21
    GS=gram(scaled);cartan=[[2*GS[i][j]/GS[j][j]for j in range(51)]for i in range(51)]
    assert cartan==[vec(row)for row in rep['cartan']]and all(x.denominator==1 for row in cartan for x in row)
    variant_keys=set()
    for row in d['exact_isometric_variants']:
        perm=row['coordinate_permutation'];assert sorted(perm)==list(range(21))
        def transform(v):
            out=[F(0)]*21
            for i,x in enumerate(v):out[perm[i]]=x
            return out
        T=list(map(transform,NN));TE=transform(E)
        assert transform(H)==H and gram(T)==GG and TE==vec(row['escape'])
        assert dot(TE,TE)==dot(E,E)and all(dot(TE,n)<0 for n in T)
        variant_keys.add(tuple(sorted(map(tuple,T))))
        changed={N.index(n[:18])for n in T[:48]if any(n[18:])}
        assert changed==set(comps[1+row['changed_component']])
    assert len(variant_keys)==len(d['exact_isometric_variants'])==12
    paths=[PRIMARY,P19,P48,f'scripts/{TAG}.py','scripts/fresh_h20_a2_h18_coupling_20261014_review.py']
    return encode(dict(status='PASS_INDEPENDENT_H17_RANK_THREE_DIRECT_COUPLINGS',sources={p:A.sha(p)for p in paths},
        arithmetic='Independent Fraction physical51-wall calculations and separate SymPy multiquadratic support inverses; no rankthree producer imports.',
        H17_19_connected=True,H17_19_orthogonal_four=[0,2,5,7],H17_48_component_sizes=list(map(len,comps)),
        H17_48_orthogonal_four=[0,1,2,3],parallel_parts=groups,cross_parity_right_pairs=len(crossrights),
        cost_one_tables=tables,all_I2_tail_checked=True,coherent_variants=12,
        physical_wall_count=51,pair_count=1275,pair_histogram=dict(hist),facet_sign_checks=2601,
        center_norm=dot(H,H),escape_norm=dot(E,E),escape_sign_checks=51,integral_cartan_entries=2601,
        scope='Independent literal source-graph/profile proof checks and complete51-wall representative/12isometries, all pairs/facets/escapes. The completeness/finitevolume of the two original H17 sources remains inherited. The direct inducing family is infinite-volume; additional noninducing walls and changed sources remain outside. No ambient vertex census or reflection-lattice reflectivity theorem claimed.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();d=run()
    out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('PASS saved independent rankthree coupling replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
