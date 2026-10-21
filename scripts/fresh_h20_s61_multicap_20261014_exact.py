#!/usr/bin/env python3
"""Exact algebraic right-chain completion of the retained S61 plane.

The one new algebraic parameter is defined by a rational isolating interval,
opposite endpoint signs, and a strictly positive derivative.  No floating
point acceptance or expanded compositum is used.  Orthogonality identities
come from the displayed rational rotation, not from interval overlap.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from collections import Counter
from fresh_h20_coupled_roles_20261014_field import Field,parse,enc,pd,B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_s61_multicap_20261014_exact'
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json',
       'results/fresh_h20_e7_bridge_face_20261014_exact.json',
       'results/fresh_h20_bridge_edge_closure_20261014.json',
       'results/fresh_h20_polar_pair_20261014.json',
       'results/fresh_h20_coupled_roles_20261014_exact.json']
BITS=448;DEN=1<<BITS
LO=F('16.76673201844734545420499946573595570803689974995962125597990920435661799831290443011705')
HI=LO+F(1,10**80)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def down(x):return F((x*DEN).__floor__(),DEN)
def up(x):return F((x*DEN).__ceil__(),DEN)
class I:
    def __init__(self,x,y=None):
        if isinstance(x,I):self.l,self.h=x.l,x.h;return
        self.l=F(x);self.h=self.l if y is None else F(y);assert self.l<=self.h
    def __add__(self,o):
        o=I(o);return I(down(self.l+o.l),up(self.h+o.h))
    __radd__=__add__
    def __neg__(self):return I(-self.h,-self.l)
    def __sub__(self,o):return self+-I(o)
    def __rsub__(self,o):return I(o)+-self
    def __mul__(self,o):
        o=I(o);v=[a*b for a in[self.l,self.h]for b in[o.l,o.h]]
        return I(down(min(v)),up(max(v)))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=I(o);assert o.l>0 or o.h<0,(o.l,o.h)
        return self*I(down(1/o.h),up(1/o.l))
    def __rtruediv__(self,o):return I(o)/self
    def sqrt(self):
        assert self.l>0
        lo=isqrt(self.l.numerator*DEN*DEN//self.l.denominator)
        hi=isqrt(self.h.numerator*DEN*DEN//self.h.denominator)+1
        return I(F(lo,DEN),F(hi,DEN))
    def out(self):return[str(self.l),str(self.h)]
class D:
    """Value and derivative, both with outward rational enclosures."""
    def __init__(self,x,d=0):self.v=I(x);self.d=I(d)
    @staticmethod
    def get(x):return x if isinstance(x,D)else D(x)
    def __add__(self,o):
        o=D.get(o);return D(self.v+o.v,self.d+o.d)
    __radd__=__add__
    def __neg__(self):return D(-self.v,-self.d)
    def __sub__(self,o):return self+-D.get(o)
    def __rsub__(self,o):return D.get(o)+-self
    def __mul__(self,o):
        o=D.get(o);return D(self.v*o.v,self.d*o.v+self.v*o.d)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=D.get(o);return D(self.v/o.v,(self.d*o.v-self.v*o.d)/(o.v*o.v))
    def sqrt(self):
        r=self.v.sqrt();return D(r,self.d/(2*r))
def chain(a):
    """Algebraic formula, upper endpoint first; 26 norm-two normals."""
    z=(a*a+1).sqrt();cc=a*a/(a*a+2);ss=2*z/(a*a+2)
    b,c=-z,D(1);half=[(a,b,c)]
    for _ in range(12):b,c=cc*b+ss*c,cc*c-ss*b;half.append((a,b,c))
    ans=list(reversed(half))+[(x,y,-z)for x,y,z in half]
    assert len(ans)==26
    return ans,cc,ss
def fdot(v,w):return -v[0]*w[0]+v[1]*w[1]+v[2]*w[2]
def pairing(v,p):return sum((a*b for a,b in zip(v,p)),I(0))
def out(v):return[x.out()for x in v]
def run():
    old,bridge,edge,polar,cur=[json.loads((ROOT/p).read_text())for p in FILES]
    for data in[old,bridge,edge,polar,cur]:
        for p,h in data['sources'].items():assert sha(p)==h,p
    A=Field(old,bridge,edge);z=A.z
    M=[tuple(parse(x,1)for x in v)for v in old['normals']]
    N=[A.vlift(v,1)for v in M]+[A.vlift(tuple(parse(x,2)for x in bridge['new_normal']),2)]
    N +=[tuple(parse(x,3)for x in v)for v in edge['new_normals']]
    N +=[tuple(A.c(x)for x in v)for v in polar['new_normals']]
    norms=[A.c(q)for q in old['norms']]+[A.c(2)]+[parse(edge['new_normal_norm'],3)]*2+[A.c(q)for q in polar['new_norms']]
    N.append(A.vlift(tuple(parse(x,1)for x in cur['cap_R']['normal']),1))
    norms.append(A.lift(parse(cur['cap_R']['norm'],1),1))
    assert len(N)==101
    S=cur['S61_limitation']['support'];assert S==edge['recipe']['support']
    roots=[M[i]for i in S]
    G=[[A.rational(A.dot(v,w,1))for w in roots]for v in roots]
    assert pd(G);J=B.inv(G)
    rec=edge['recipe'];assert int(rec['boost_integer'])==1
    X=tuple(parse(x,2)for x in rec['X']);W=tuple(parse(x,2)for x in rec['W'])
    rr=tuple(parse(x,2)for x in bridge['new_normal'])
    AA=parse(rec['A'],2);ww=parse(rec['W_norm'],2)
    om=tuple(parse(x,2)for x in rec['dual_weight']);part=tuple(parse(x,2)for x in rec['particular'])
    assert part==A.lin([om,rr],[A.c(-1,2),A.c(F(-1,2),2)],2)
    assert A.dot(om,om,2)==A.c(F(2,3),2)
    frame=[X,rr,W];fsq=[AA,A.c(2,2),ww]
    assert [[A.dot(v,w,2)for w in frame]for v in frame]==[
        [A.scale(AA,-1,2),z[2],z[2]],[z[2],A.c(2,2),z[2]],[z[2],z[2],ww]]
    assert all(A.dot(v,A.vlift(M[i],1,2),2)==z[2]for v in frame for i in S)
    assert all(A.dot(om,v,2)==z[2]for v in frame)
    delta=parse(edge['field']['new_positive_sqrt_radical'],2)
    assert delta==A.mul(AA,A.add(ww,A.c(F(-5,6),2),2),2)
    assert A.sign(AA,2)>0
    for raw,sg in zip(edge['new_normals'],[1,-1]):
        assert tuple(parse(x,3)for x in raw)==tuple((A.mul(AA,A.add(part[j],A.scale(W[j],sg,2),2),2),X[j])for j in range(21))
    def bound(v,k=3):return I(*A.bounds(v,k,BITS))
    fden=[bound(q,2).sqrt()for q in fsq]
    frame3=[A.vlift(v,2)for v in frame]
    # Every row is the exact physical pairing with each normalized frame vector.
    # Old normals are positively normalized to square two for uniform pair gates.
    P=[]
    for i,(n,q)in enumerate(zip(N,norms)):
        den=(bound(q)/2).sqrt()
        row=[bound(A.dot(v,n))/(r*den)for v,r in zip(frame3,fden)]
        if i in S:
            assert all(A.dot(v,n)==z[3]for v in frame3);row=[I(0)]*3
        P.append(row)
    hb=[tuple(parse(x,1)for x in cur['hex_basis'][key])for key in['C','U','V']]
    hq=[parse(x,1)for x in cur['hex_basis']['squared_lengths']]
    for co in cur['hexagon']['frame_coefficients']:
        row=[]
        for v,den in zip(frame3,fden):
            row.append(sum((bound(A.lift(parse(c,0),0))*bound(A.dot(v,A.vlift(u,1)))/(den*bound(A.lift(q,1)).sqrt())for c,u,q in zip(co,hb,hq)),I(0)))
        P.append(row)
    assert len(P)==107
    wb=bound(ww,2);rootw=wb.sqrt();kk=(wb-F(5,6)).sqrt();root2=I(2).sqrt()
    # Literal predecessor identity supplies the exact side rows, replacing
    # dependency-inflated enclosures by this independently checked formula.
    P[95]=[I(0),root2,I(0)]
    P[96]=[-kk,-1/root2,rootw];P[97]=[-kk,-1/root2,-rootw]
    def endpoint(a):
        co,_,_=chain(a);a,b,c=co[0]
        return a*D(kk)+b/D(root2)-c*D(rootw)
    flo=endpoint(D(LO));fhi=endpoint(D(HI));whole=endpoint(D(I(LO,HI),1))
    assert flo.v.h<0<fhi.v.l and whole.d.l>0
    co,rotc,rots=chain(D(I(LO,HI),1));coeff=[[v.v for v in row]for row in co]
    assert rotc.v.l>0 and rots.v.l>0
    assert all(c[1].h<0 for c in coeff)
    assert all(c[2].l>0 for c in coeff[:13])and all(c[2].h<0 for c in coeff[13:])
    assert all(coeff[i][2].l>coeff[i+1][2].h for i in range(25))
    # Formal exact identities: C²+S²=1, z²=a²+1, every norm=2,
    # adjacent product0. Reflection in the third coordinate gives the second
    # endpoint equation from the first. The derivative above isolates its root.
    rows=[];counts=Counter();allnew=[]
    for j,c in enumerate(coeff):
        row=[]
        for i,p in enumerate(P):
            g=pairing(c,p)
            right=i in S or(j==0 and i==96)or(j==25 and i==97)
            if right:assert g.l<=0<=g.h or i in S;rel='right'
            else:assert g.h<-2,(j,i,g.out());rel='dotted'
            counts[rel]+=1;row.append(dict(wall=i,relation=rel,bounds=g.out()))
        rows.append(row)
        for k in range(j):
            g=fdot(c,coeff[k]);rel='right'if j-k==1 else'dotted'
            if rel=='right':assert g.l<=0<=g.h
            else:assert g.h<-2,(j,k,g.out())
            counts[rel]+=1;allnew.append(dict(walls=[107+k,107+j],relation=rel,bounds=g.out()))
    assert sum(counts.values())==26*107+26*25//2
    # Cyclic S61 halfplanes: old95, upper96, chain, lower97.
    facecoeff=[[I(0),root2,I(0)],[kk,-1/root2,rootw]]+coeff+[[kk,-1/root2,-rootw]]
    faceids=[95,96]+list(range(107,133))+[97]
    n=len(faceids);assert n==29
    vertices=[];excluded=[]
    for i in range(n):
        for j in range(i+1,n):
            u,v=facecoeff[i],facecoeff[j];det=u[1]*v[2]-u[2]*v[1]
            assert det.l>0 or det.h<0,(i,j,'parallel spatial rows')
            x=(u[0]*v[2]-u[2]*v[0])/det
            y=(u[1]*v[0]-u[0]*v[1])/det
            pt=[I(1),x,y];facepairs=[fdot(pt,c)for c in facecoeff]
            adjacent=j==i+1 or(i==0 and j==n-1)
            if not adjacent:
                blockers=[k for k,g in enumerate(facepairs)if g.l>0]
                assert blockers,(i,j,'nonconsecutive feasible point')
                k=max(blockers,key=lambda k:facepairs[k].l)
                excluded.append(dict(walls=[faceids[i],faceids[j]],positive_blocker=faceids[k],pair_bound=facepairs[k].out()))
                continue
            q=fdot(pt,pt);assert q.h<0
            # Positive time coordinate fixes the same future sheet as X.
            assert all(g.h<0 for k,g in enumerate(facepairs)if k not in[i,j])
            own={faceids[i],faceids[j]};oldpairs=[]
            for k,p in enumerate(P):
                g=pairing(pt,p)
                if k in S or k in own:assert g.l<=0<=g.h or k in S
                else:assert g.h<0,(i,j,k,g.out())
                oldpairs.append(g.out())
            vertices.append(dict(walls=sorted(own),frame_point_bounds=out(pt),norm_bound=q.out(),old107_pair_bounds=oldpairs,face_pair_bounds=out(facepairs),support=S+sorted(own),new=bool(own&set(range(107,133)))))
    assert len(vertices)==29 and sum(v['new']for v in vertices)==27 and len(excluded)==377
    # Boundedness is exact before imposing the hyperbolic unit disk:
    # n_middle+ + n_middle- + sqrt(2)*z*r95 = (2a,0,0).
    # These three spatial normals span R² and have a strictly positive relation.
    # Thus the spatial halfplane intersection is bounded. The exhaustive pair
    # census above proves it has exactly29 vertices, all strictly inside disk.
    # Ambient strict center: C=x-e/2 + eps*H_old/sqrt(-H_old²).
    Hold=tuple(parse(x,1)for x in cur['strict_center']['point'])
    qH=A.dot(Hold,Hold,1);assert A.sign(qH,1)<0
    hden=bound(A.scale(qH,-1,1),1).sqrt()
    hf=[bound(A.dot(v,A.vlift(Hold,1)))/(d*hden)for v,d in zip(frame3,fden)]
    Cface=[I(1),I(F(-1,2)),I(0)];eps=F(1,2**30)
    assert pairing(Cface,hf).h<0
    centerpairs=[]
    for c in coeff:
        g=fdot(c,Cface)+eps*pairing(c,hf);assert g.h<0;centerpairs.append(g.out())
    for i,p in enumerate(P):
        g=pairing(Cface,p)
        if i in S:assert g.l==g.h==0
        else:assert g.h<0
    cq=I(F(-3,4))+2*eps*pairing(Cface,hf)-eps*eps;assert cq.h<0
    # Exactly27 finite vertices on the new walls: their complete finite
    # neighbor lists contain S and only the cyclic adjacent boundary walls.
    # The two additional neighbors are mutually dotted, so only either
    # consecutive20-support is spherical. Their actual points were checked.
    neighbors=[S+([96]if j==0 else[106+j])+([97]if j==25 else[108+j])for j in range(26)]
    # Proper bounded old cusps remain strict, checked directly as well.
    cuspdata=[]
    for raw in[old['completed_cusp']['point']]+[r['point']for r in old['preserved_old_bounded_cusps']]:
        nu=tuple(A.lift(parse(x,0),0)for x in raw)
        pp=[bound(A.dot(v,nu))/d for v,d in zip(frame3,fden)]
        vals=[pairing(c,pp)for c in coeff];assert all(g.h<0 for g in vals)
        cuspdata.append(dict(point=raw,frame_pair_bounds=out(pp),new_pair_bounds=out(vals)))
    # Current107 rational escape survives. Its norm is rechecked from Q.
    eraw=cur['strict_escape']['rational_point'];E=tuple(A.c(x)for x in eraw)
    qe=A.dot(E,E);assert qe==A.c(F(38017800175279,500000000000))
    ep=[bound(A.dot(v,E))/d for v,d in zip(frame3,fden)]
    eg=[pairing(c,ep)for c in coeff];assert all(g.h<0 for g in eg)
    efuture=ep[0]-ep[1]/2+eps*bound(A.dot(E,A.vlift(Hold,1)))/hden
    assert efuture.h<0
    # New caps are redundant on the already compact S75 plane. Only their75
    # pairing can be nonzero on S75; the saved inverse diagonal there is5.
    S75=cur['hex_basis']['support'];assert set(S75)-set(S)=={75}
    assert norms[75]==A.c(2)
    J75=[[F(x)for x in row]for row in cur['hex_basis']['support_inverse']]
    assert J75[S75.index(75)][S75.index(75)]==5
    s75norms=[I(2)-5*pairing(c,P[75])*pairing(c,P[75])for c in coeff]
    assert all(g.h<0 for g in s75norms)
    paths=FILES+[f'scripts/{TAG}.py','scripts/fresh_h20_coupled_roles_20261014_field.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py','results/acute_extension_preservation_20261014.md']
    return enc(dict(status='PASS_133_WALL_CHAMBER_WITH_COMPACT_S61_29_GON',sources={p:sha(p)for p in paths},interval_bits=BITS,
        support=S,support_gram=G,support_inverse=J,frame=dict(X=X,A=AA,r95=rr,W=W,W_norm=ww,normalization='x=X/sqrt(A), e=r95/sqrt(2), f=W/sqrt(W_norm)',old107_unit_normal_pair_bounds=[out(row)for row in P]),
        parameter=dict(name='a',isolating_interval=[LO,HI],endpoint_lower_sign=flo.v.out(),endpoint_upper_sign=fhi.v.out(),derivative_bound=whole.d.out(),equation='a*sqrt(w-5/6)+b_endpoint/sqrt(2)-c_endpoint*sqrt(w)=0',definition='The unique real root in the displayed rational interval; existence by opposite signs and uniqueness by strictly positive derivative.'),
        chain=dict(walls=list(range(107,133)),norms=[2]*26,recurrence='z=sqrt(a²+1), C=a²/(a²+2), S=2z/(a²+2). Starting(b,c)=(-z,1), rotate(b,c)=(Cb+Sc,Cc-Sb) twelve times; list the13 upper rows in reverse order and then their reflected lower rows in forward order. Every time coefficient is a.',frame_coefficient_bounds=[out(c)for c in coeff],rotation_cosine_bound=rotc.v.out(),rotation_sine_bound=rots.v.out(),exact_Gram_identities='C²+S²=1; initial norm -a²+z²+1=2; adjacent inner product0, including the central reflected pair. The root equation and reflection give contacts96/107 and97/132 exactly0.',old107_pairs=rows,new_new_pairs=allnew,finite_neighbors=neighbors),
        new_pair_counts=dict(counts),new_pair_count=sum(counts.values()),face=dict(cyclic_walls=faceids,vertices=vertices,nonconsecutive_intersections_excluded=excluded,boundedness_identity='n119+n120+sqrt(2)*sqrt(a²+1)*r95=(2a,0,0) in the unit Lorentz frame; the three spatial normals span R² and their coefficients are positive.',compactness='The spatial halfplane polygon is bounded. All406 boundary-pair intersections are classified:29 consecutive intersections are feasible future timelike, and377 others violate a retained polygon halfspace. Hence the whole S61 face is a compact29-gon. The two old corners have their inherited projected angles; no claim that the entire intrinsic polygon is Coxeter.'),
        strict_center=dict(recipe='x-e/2+epsilon*H107/sqrt(-H107²)',epsilon=eps,norm_bound=cq.out(),new_pair_bounds=centerpairs,proof='The face point is strict on every old wall outsideS61; positive inward old-center perturbation makes allS61 inequalities strict. The displayed bounds check every new inequality.'),
        complete_new_vertex_count=27,known_total_vertices_at_least=168,local_census_proof='Each new wall has exactly S61 plus its two cyclic neighbors as finite-angle neighbors. Those two neighbors are mutually dotted. Thus every possible finite20-support uses all18 S61 roots and one adjacent boundary; precisely the27 displayed new points result. All141 known old vertices are retained by acute extension.',
        preserved_bounded_cusps=cuspdata,S75_projected_norm_bounds=out(s75norms),strict_escape=dict(rational_point=eraw,norm=F(38017800175279,500000000000),frame_pair_bounds=out(ep),new_pair_bounds=out(eg),new_strict_center_pair_bound=efuture.out(),old107_signs_dependency=FILES[-1],meaning='The same future-oriented strict spacelike escape survives all133 walls, proving infinite hyperbolic volume.'),
        scope='Exact26-wall coherent addition retaining107. Both named S61/S75 faces are now compact, with27 new finite vertices and at least168 overall, but the chamber still has infinite volume. No complete ambient vertex or edge census is asserted.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();data=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['new_pair_counts'],data['complete_new_vertex_count'],flush=True)
