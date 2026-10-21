#!/usr/bin/env python3
"""Exact global escape for radial and single-factor face-cap replacements.

The metric candidate has only a numerical all-pair audit.  This checker
instead certifies its norm/contact identities and an escape for the larger
continuous family, without assuming compatibility of that relaxation.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_endpoint_replacement_20261014'
LINKS=list(range(5))+list(range(6,11))+list(range(12,18))

def dot(a,b):
    return sum(x*y for x,y in zip(a,b))-sum((a[i]*b[i+1]+a[i+1]*b[i])/2 for i in LINKS)

def inverse_A(n,rhs):
    return [sum(F(2*min(i,j)*(n+1-max(i,j)),n+1)*rhs[j-1] for j in range(1,n+1)) for i in range(1,n+1)]

def add(a,b):return [x+y for x,y in zip(a,b)]
def scale(a,t):return [t*x for x in a]
def sub(a,b):return add(a,scale(b,-1))
def sign5(a,b):
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if a*b>0:return (a>0)-(a<0)
    d=a*a-5*b*b
    assert d
    return ((d>0)-(d<0))*((a>0)-(a<0))

def run():
    # Reconstruct all 392 unit-inradius product vertices without imports.
    factors=[];seed=[];start=0
    for n in (6,6,7):
        vv=[]
        for o in range(n+1):
            v=inverse_A(n,[F(1)-(n+1 if i==o else 0) for i in range(1,n+1)])
            vv.append([F(0)]*start+v+[F(0)]*(19-start-n))
            w=[F(0)]*19
            if o==0:w[start:start+n]=[-F(1)]*n
            else:w[start+o-1]=F(1)
            seed.append(w)
        factors.append(vv);start+=n
    local=list(product(range(7),range(7),range(8)))
    Z=[add(add(factors[0][a],factors[1][b]),factors[2][c]) for a,b,c in local]
    assert all(dot(z,z)==196 for z in Z)
    index=lambda a:a[0]*56+a[1]*8+a[2]
    lo,hi=F(16),F(49,3)
    # Adjacent A6 triangle copies have cosine c(r); its positive pi/5 root
    # lies strictly inside [16,49/3].  The polynomial has one positive root.
    # 4f(r)=(4-4sqrt5)r²+(49sqrt5-539)r+8232.
    for r,want in [(lo,1),(hi,-1)]:
        assert sign5(4*r*r-539*r+8232,-4*r*r+49*r)==want
    assert sign5(F(4),-F(4))<0 and sign5(-F(539),F(49))<0

    # Check norm, gap, and right-contact identities as polynomial identities.
    # Every expression has degree at most two; three rational values suffice.
    for f,n in enumerate((7,7,8)):
        ids=[]
        for j in range(3):
            a=[0,0,0];a[f]=j;ids.append(index(a))
        z0,z1,z2=[Z[i] for i in ids];m01=scale(add(z0,z1),F(1,2));m12=scale(add(z1,z2),F(1,2))
        h=F(n*(n-1),2)
        assert dot(sub(z0,z1),sub(z0,z1))==4*h
        for r in (F(15),F(16),F(17)):
            a=r/(n*n);w=add(scale(add(z0,z2),a),scale(z1,1-2*a))
            q=4*r*r/(n*n)-r
            assert dot(w,w)==196-4*F(n-1,n)*r+4*r*r/(n*n)
            assert dot(sub(w,z1),sub(w,z1))==q+r
            assert dot(sub(w,m01),sub(w,m01))==q+h-r
            assert dot(sub(w,m12),sub(w,m12))==q+h-r
            gaps=[1-dot(w,s) for s in seed]
            base=(0,7,14)[f]
            assert gaps[base:base+3]==[r/n,n-2*r/n,r/n]
            # Shifted neighboring triangle center and its mutual cosine.
            aa=[0,0,0];aa[f]=3;z3=Z[index(aa)]
            w2=add(scale(add(z1,z3),a),scale(z2,1-2*a))
            dd=dot(sub(w,w2),sub(w,w2))
            assert dd==12*r*r/(n*n)-8*r+2*n*(n-1)
            if n==7:
                assert (dd-2*q)*49/2==2*r*r-147*r+2058

    # A cross-factor square cannot be covered by the entire parameter family.
    support=[0,56,8,64]
    m=scale([sum(Z[i][j] for i in support) for j in range(19)],F(1,4))
    assert dot(m,m)==154
    t=F(99,100);x=scale(m,t);xq=dot(x,x)
    assert min(1-dot(x,s) for s in seed)==F(1,100)
    minima={};attainers={}
    def record(kind,p,recipe):
        if kind not in minima or p<minima[kind]:minima[kind]=p;attainers[kind]=recipe
    record('central',xq-4,None)
    for i,z in enumerate(Z):
        record('endpoint_r_in_interval',dot(sub(x,z),sub(x,z))-hi,dict(vertex=i,r=str(hi)))
        # This minimum is over ALL REAL gamma, a deliberate relaxation of
        # gamma<1 and R<=2(1-gamma) forced by the retained positive seed gaps.
        d=dot(x,z)-4
        record('all_real_vertex_radial',xq-4-d*d/192,dict(vertex=i))
    for f,n in enumerate((7,7,8)):
        h=F(n*(n-1),2)
        for abc in local:
            b=list(abc);b[f]=(b[f]+1)%n
            c=list(abc);c[f]=(c[f]+2)%n
            ids=[index(abc),index(b),index(c)]
            zz=[Z[i] for i in ids]
            center=scale(add(zz[0],zz[1]),F(1,2))
            # Edge power is affine increasing in r.
            record('all_nearest_edge_caps',dot(sub(x,center),sub(x,center))-(h-lo),dict(factor=f,vertices=ids[:2],r=str(lo)))
            # The triangle center is affine and its radius² quadratic in r.
            # Quadratic coefficients cancel, leaving an affine ball power.
            vals=[]
            for r in (lo,(lo+hi)/2,hi):
                a=r/(n*n);w=add(scale(add(zz[0],zz[2]),a),scale(zz[1],1-2*a));q=4*r*r/(n*n)-r
                vals.append(dot(sub(x,w),sub(x,w))-q)
            assert 2*vals[1]==vals[0]+vals[2]
            for r,p in ((lo,vals[0]),(hi,vals[2])):
                record('all_consecutive_triangle_caps',p,dict(factor=f,vertices=ids,r=str(r)))
    oldpath='results/fresh_h20_balanced_integer_exact_20261014.json'
    old=json.loads((ROOT/oldpath).read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    for i,row in enumerate(old['walls'][415:],415):
        w=list(map(F,row));q=dot(w[:19],w[:19])-2*w[19]
        record('old_three_motif_and_392_middle',dot(sub(x,w[:19]),sub(x,w[:19]))-q,dict(wall=i))
    assert len(old['walls'][415:])==395
    assert min(minima.values())>1
    D=(xq-1)/2
    assert xq-2*D==1 and -200-D<0

    # Separate simpler radial obstruction near one nearest-edge midpoint,
    # even allowing every replacement endpoint radius²<=21.
    em=scale(add(Z[0],Z[56]),F(1,2));ex=scale(em,t);eq=dot(ex,ex)
    assert dot(em,em)==175
    er=min(eq-4-(dot(ex,z)-4)**2/192 for z in Z)
    ep=min(dot(sub(ex,z),sub(ex,z))-21 for z in Z)
    assert er==F(1407119,76800) and ep==F(7,400)>F(1,100)
    paths=[f'scripts/{TAG}.py',oldpath]
    return dict(status='PASS_EXACT_CONTINUOUS_RADIAL_AND_FACE_HIERARCHY_ESCAPE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        arithmetic='stdlib Fraction; no project imports',
        parameter_interval=['16','49/3'],
        metric_root_polynomial='(4-4sqrt5)r²+(49sqrt5-539)r+8232=0, unique positive root',
        endpoint_radius_squared='r',edge_radius_squared='n(n-1)/2-r',
        triangle_center='(r*z0+(n²-2r)*z1+r*z2)/n²',triangle_radius_squared='4r²/n²-r',
        literal_counts=dict(vertices=392,nearest_edges=1176,consecutive_triangles=1176,old_middle_and_motif=395),
        exact_contact_identities='Triangle cap right to middle endpoint and both neighboring edge caps; A6 adjacent triangle cosine is the pi/5 polynomial ratio.',
        common_escape=dict(spatial= list(map(str,x)),spatial_square=str(xq),D=str(D),U=1,norm=1,
            source_face_vertices=support,inward_scale='99/100',minimum_seed_gap='1/100',
            minimum_ball_powers={k:str(v) for k,v in minima.items()},attainers=attainers),
        simpler_radial_escape=dict(midpoint_vertices=[0,56],inward_scale='99/100',
            endpoint_radius_squared_upper_bound=21,radial_power_lower_bound=str(er),
            endpoint_power_lower_bound=str(ep),usable_spacelike_norm='1/100'),
        scope='A common strict escape excludes every finite compatible selection from this continuous family, including optional deletions. It does not exclude arbitrary mixed-factor face caps or differently positioned single-factor caps. Numerical full-family Coxeter compatibility at the metric root is a separate scout, not certified here.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print(out['common_escape']['minimum_ball_powers'])
