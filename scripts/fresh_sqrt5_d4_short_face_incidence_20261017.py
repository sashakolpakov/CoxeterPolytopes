#!/usr/bin/env python3
"""Exact local identities for D4-copy vertices, with a scoped incidence bound."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as A
K,Z,O=A.K,A.Z,A.O;ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_d4_short_face_incidence_20261017';PHI=K(F(1,2),F(1,2))

def d4(copies=1):
    n=4*copies;g=[[K(2*(i==j))for j in range(n)]for i in range(n)]
    for c in range(copies):
        for j in range(1,4):g[4*c][4*c+j]=g[4*c+j][4*c]=-O
    return g

def is_d4(g):
    return len(g)==4 and all(g[i][i]==K(2)for i in range(4)) and all(g[i][j]in(Z,-O)for i,j in combinations(range(4),2)) and sorted(sum(bool(x)for j,x in enumerate(r)if i!=j)for i,r in enumerate(g))==[1,1,1,3]

def extend(g,contacts):
    return [r+[contacts[i]]for i,r in enumerate(g)]+[contacts+[K(2)]]

def copy_endpoint(g,i,framing):
    return extend(g,[framing if j==i else g[i][j]for j in range(len(g))])

def reattachments():
    g=d4();out=[]
    for port in range(4):
        ids=[i for i in range(4)if i!=port];support=A.sub(g,ids)
        accepted=[]
        for bits in product((0,-1),repeat=3):
            if is_d4(extend(support,list(map(K,bits)))):accepted.append(bits)
        expected=tuple(g[port][j].a for j in ids)
        assert accepted==[expected]
        out.append(dict(port=port,retained=ids,unique_contacts=list(map(str,expected))))
    return out

# Two-variable rational polynomials; determinant verification is symbolic.
def pc(x):return {(0,0):F(x)}if x else{}
def plus(a,b):
    c=a.copy()
    for k,v in b.items():
        c[k]=c.get(k,F(0))+v
        if not c[k]:del c[k]
    return c
def neg(a):return{k:-v for k,v in a.items()}
def mul(a,b):
    c={}
    for k,v in a.items():
        for l,w in b.items():
            key=(k[0]+l[0],k[1]+l[1]);c[key]=c.get(key,F(0))+v*w
    return{k:v for k,v in c.items()if v}
def determinant(m):
    out={};n=len(m)
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j]for i in range(n)for j in range(i+1,n));term=pc(sign)
        for i in range(n):term=mul(term,m[i][p[i]])
        out=plus(out,term)
    return out

def schur_identities():
    out=[];li={(1,0):F(1)};lj={(0,1):F(1)}
    for copies in(1,2):
        g=d4(copies);gi=A.inverse(g)
        pairs=list(combinations(range(4),2))if copies==1 else list(product(range(4),range(4,8)))
        for i,j in pairs:
            h=A.inverse(A.sub(gi,[i,j]));assert all(not v.b for r in h for v in r)
            a,b,c=h[0][0].a,h[1][1].a,h[0][1].a
            aa,bb,cc=pc(a),pc(b),pc(c);al=plus(aa,neg(li));bl=plus(bb,neg(lj))
            block=[[aa,cc,al,cc],[cc,bb,cc,bl],[al,cc,aa,cc],[cc,bl,cc,bb]]
            wanted=mul(mul(li,lj),plus(mul(plus(pc(2*a),neg(li)),plus(pc(2*b),neg(lj))),pc(-4*c*c)))
            assert determinant(block)==wanted
            kind='cross_component'if copies==2 else('center_leaf'if 0 in(i,j)else'leaf_leaf')
            if kind=='leaf_leaf':assert(a,b,c)==(F(4,3),F(4,3),F(-2,3))
            if kind=='cross_component':assert c==0
            out.append(dict(ports=[i,j],copies=copies,kind=kind,schur=[[str(x.a)for x in r]for r in h],
                determinant=[[list(k),str(v)]for k,v in sorted(wanted.items())]))
    assert Counter(r['kind']for r in out)=={'center_leaf':3,'leaf_leaf':3,'cross_component':16}
    return out

def norm_obstruction():
    # O_K=Z[phi], N(x+y phi)=x^2+xy-y^2. Mod2 it is anisotropic.
    norm=lambda x,y:x*x+x*y-y*y
    zero_mod2=[(x,y)for x,y in product(range(2),repeat=2)if norm(x,y)%2==0]
    assert zero_mod2==[(0,0)]
    assert norm(1,0)%3==1
    possible_negative_norms=[-(4**r)for r in range(5)]
    assert all(n%3==2 for n in possible_negative_norms)
    return dict(ring='Z[phi], phi^2=phi+1',norm='x^2+x*y-y^2',zero_norm_classes_mod2=zero_mod2,
        leaf_pair_parameters='a=-3*g_i-2, b=-3*g_j-2; a*b=16',
        residue='a=b=1 modulo3, hence Norm(a)=1 modulo3',
        signs='a,b>0 physically and a_conjugate,b_conjugate<0, hence Norm(a)<0',
        possible_negative_norms=possible_negative_norms,
        proof='Norm(a) divides Norm(16)=256. An even norm forces both coordinates even, so its2-adic valuation is even. Thus the negative possibilities are exactly -4^r,0<=r<=4, all2 modulo3. This contradicts the required residue1.')

def triangles():
    g=d4();tests=[];allowed=[]
    for port in range(4):
        for label,value in[(2,Z),(3,-O),(5,-PHI)]:
            big=copy_endpoint(g,port,value);pi=A.inertia(big);co=A.inertia([[x.conj()for x in r]for r in big])
            for other in range(4):
                if other==port:continue
                third=A.sub(big,[i for i in range(5)if i!=other]);finite=A.inertia(third)==[4,0,0]
                ok=pi==[4,1,0]and co==[5,0,0]and finite
                record=dict(port=port,other=other,label=label,physical=pi,conjugate=co,third_corner_finite=finite,accepted=ok)
                tests.append(record)
                if ok:
                    assert port in(1,2,3)and other==0 and label==5
                    assert not is_d4(third);allowed.append([port,other,label])
    assert allowed==[[1,0,5],[2,0,5],[3,0,5]]
    # If the other port lies in a different copy, the old center survives.
    # Its three-root subdiagram with the old/new leaf is3,3,5 and indefinite.
    triangle=[[K(2),-O,-PHI],[-O,K(2),-O],[-PHI,-O,K(2)]]
    assert A.inertia(triangle)==[2,1,0]
    return dict(tests=tests,accepted=allowed,cross_component_obstruction=A.enc(triangle),
        finite_label_reduction='The shared old neighbor gives a triangular cyclic product equal to the framing g, so g lies integrally in K. If its pair occurs in a finite corner, the real cyclotomic degree restriction gives m=2,3,4,5,6; sqrt2 and sqrt3 are not in Q(sqrt5), leaving2,3,5.',
        consequence='No triangle has three D4-copy vertices. A triangle with exactly two can use only a center-leaf pair at each of these vertices; the third vertex replaces that D4 by2A1+I2(5).')

def geometry(quad):
    g=d4()
    if quad:
        g=copy_endpoint(g,0,-2*PHI)
        g=extend(g,[-O,-2*PHI,Z,Z,-O]);black=[0,1,4,5];expected=[[0,1],[0,5],[1,4],[4,5]]
        pi=[4,1,1];co=[5,0,1];ids=list(range(5))
    else:
        g=copy_endpoint(g,1,-PHI);black=[0,1,4];expected=[[0,1],[0,4],[1,4]]
        pi=[4,1,0];co=[5,0,0];ids=list(range(5))
    assert A.inertia(g)==pi and A.inertia([[x.conj()for x in r]for r in g])==co
    assert all(x.integral()for r in g for x in r)
    q=A.sub(g,ids);qi=A.inverse(q);normals=[A.mv(qi,[g[j][i]for j in ids])for i in range(len(g))]
    assert [[A.pair(q,x,y)for y in normals]for x in normals]==g
    h=A.mv(qi,[-O]*5);assert A.pair(q,h,h).sign()<0
    hp=[A.pair(q,h,n)for n in normals];assert max(x.sign()for x in hp)<0
    vertices=[];census=[]
    for pair in combinations(black,2):
        active=[2,3]+list(pair);ker=A.kernel([A.mv(q,normals[i])for i in active]);assert len(ker)==1
        signs=[]
        for sign in(-1,1):
            v=[sign*x for x in ker[0]];products=[A.pair(q,v,n)for n in normals]
            if max(x.sign()for x in products)<=0:
                assert A.pair(q,v,v).sign()<0 and A.pair(q,v,h).sign()<0
                assert all(products[j].sign()<0 for j in range(len(g))if j not in active)
                assert A.inertia(A.sub(g,active))==[4,0,0]
                if quad:assert is_d4(A.sub(g,active))
                vertices.append(dict(black=list(pair),coordinates=[x.enc()for x in v],norm=A.pair(q,v,v).enc(),products=[x.enc()for x in products]));signs.append(sign)
        census.append(dict(black=list(pair),feasible_signs=signs))
    assert sorted(v['black']for v in vertices)==expected
    # Every wall is actual: project the strict timelike center onto that wall.
    facets=[]
    for i,n in enumerate(normals):
        v=[h[j]-hp[i]*n[j]/2 for j in range(5)];pp=[A.pair(q,v,m)for m in normals]
        assert not pp[i]and all(z.sign()<0 for j,z in enumerate(pp)if j!=i)and A.pair(q,v,v).sign()<0
        facets.append(A.pair(q,v,v).enc())
    return dict(gram=A.enc(g),physical_inertia=pi,conjugate_inertia=co,white=[2,3],black=black,
        basis_indices=ids,root_coordinates=A.enc(normals),center=[x.enc()for x in h],
        center_norm=A.pair(q,h,h).enc(),vertices=vertices,complete_face_census=census,actual_facet_projection_norms=facets,
        scope='Exact compact local triangle/quadrilateral in a partial arithmetic chamber. No compact ambient polytope or finite-covolume group is asserted.')

def incidence():
    out=[]
    for k in range(2,8):
        d=4*k;c=d*(d-1)//2;t=F(5*((d-1)**2//4),2)-c;alpha=F(3*k,2)
        assert t==2*k*k-3*k and 2*(t-alpha)==k*(4*k-9)and c-t==k*(6*k+1)
        out.append(dict(dimension=d,copies=k,corners=c,threshold=str(t),special_vertex_charge=str(alpha),
            special_coefficient=str(t-alpha),other_coefficient=str(c-t),
            inequality=f'{t-alpha} V_D < {c-t} V_other'if t>alpha else'No useful positive frequency bound from this estimate'))
    return dict(corner_capacity='2*T_two_D+4*Q_all_D <=3*k*V_D',
        payment='2*T+Q <=(3*k/2)*V_D+binom(d,2)*V_other',
        charging='Each triangle with two D4-copy vertices receives1 from its other vertex and1/2 from each special corner. Every other triangle has at least two other vertices and is paid at most1 per such corner. Every quadrilateral with another vertex receives1 from one other corner; an all-special quadrilateral receives1/4 per special corner. All special charged corners are among the3k center-leaf pairs.',
        strict_h_vector='T(d)*(V_D+V_other) <2*T+Q, where T(d)=5*floor((d-1)^2/4)/2-binom(d,2)',
        dimensions=out,general_frequency_inequality='(4*k-9)*V_D < (12*k+2)*V_other for k>=3',
        consequence='For d=4k>=12, a compact arithmetic polytope over Q(sqrt5) cannot have every vertex forest equal to kD4. This is a type-frequency obstruction, not a general dimension bound.')

def run():
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_adjacent_discount_witness_20261016.py',
           'scripts/fresh_sqrt5_discount_structural_floor_20261016.py','results/fresh_sqrt5_discount_structural_floor_20261016.json']
    return dict(status='PASS_EXACT_D4_SHORT_FACE_IDENTITIES_AND_SCOPED_INCIDENCE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},reattachments=reattachments(),
        schur_cases=schur_identities(),leaf_pair_arithmetic_obstruction=norm_obstruction(),triangle_cases=triangles(),
        sharp_local_triangle=geometry(False),allowed_center_leaf_quadrilateral=geometry(True),incidence=incidence(),
        assumptions='Compact arithmetic Coxeter polytope of simplest type over K=Q(sqrt5), dimension d=4k; special vertices have exactly kD4 as their full finite diagram. No narrowness/minimality bound is used.',
        quadrilateral_argument='Same-type edge attachments are forced copies. For ports in different D4 copies, the two edge Schur inequalities give lambda_i>2a and lambda_j>2b, while c=0, so the four-root determinant cannot vanish. For leaf-leaf ports in one copy, the exact determinant and the two embedding signatures give a*b=16 with opposite embedding signs and a=b=1 modulo3, contradicted by the norm argument. Only center-leaf ports remain.',
        limitations='The D4-copy condition is essential. Mixed D5/D4 forests can move components during endpoint replacement and are not classified here. No universal lowering of a corner payment, no completed new weight rule and no unconditional dimension bound follows.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=json.loads(json.dumps(run()))
    p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
