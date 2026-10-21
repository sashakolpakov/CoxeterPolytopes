#!/usr/bin/env python3
"""Small independent Fraction checks supporting the source/formula review.

No primary imports. This is not a second all-family transport enumeration.
The full exact primary saved replay was run separately during review.
"""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_endpoint_replacement_20261014_review'
SOURCE='results/fresh_h20_endpoint_replacement_20261014.json'
def solve(A,b):
    B=[list(r)+[v] for r,v in zip(A,b)];n=len(B)
    for i in range(n):
        j=next(j for j in range(i,n) if B[j][i]);B[i],B[j]=B[j],B[i]
        p=B[i][i];B[i]=[x/p for x in B[i]]
        for j in range(n):
            if j!=i:
                p=B[j][i];B[j]=[x-p*y for x,y in zip(B[j],B[i])]
    return [r[-1] for r in B]
def run():
    data=json.loads((ROOT/SOURCE).read_text());paths={SOURCE}
    for p,h in data['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p;paths.add(p)
    B=[[F(i==j) for j in range(19)] for i in range(19)]
    offsets=[0,6,12];dims=[6,6,7]
    for a,n in zip(offsets,dims):
        for i in range(a,a+n-1): B[i][i+1]=B[i+1][i]=-F(1,2)
    def dot(x,y):return sum(x[i]*B[i][j]*y[j] for i in range(19) for j in range(19))
    def vertex(ids):
        rhs=[]
        for o,n in zip(ids,dims):rhs += [F(1-(n+1)*(i==o)) for i in range(1,n+1)]
        return solve(B,rhs)
    def sub(a,b):return [x-y for x,y in zip(a,b)]
    x=list(map(F,data['common_escape']['spatial']));q=dot(x,x)
    assert q==F(754677,5000)==F(data['common_escape']['spatial_square'])
    D=F(data['common_escape']['D']);assert q-2*D==1 and -200-D<0
    seed=[]
    for a,n in zip(offsets,dims):
        for j in range(n):
            w=[F(0)]*19;w[a+j]=1;seed.append(w)
        w=[F(0)]*19;w[a:a+n]=[-F(1)]*n;seed.append(w)
    assert min(1-dot(x,s) for s in seed)==F(1,100)
    vertices=[vertex(t) for t in [(0,0,0),(1,0,0),(0,1,0),(1,1,0)]]
    assert x==[F(99,400)*sum(v[i] for v in vertices) for i in range(19)]
    # Polynomial coefficients of |x-gamma*z|²-4(1-gamma)²:
    # (q-4), -2(xz-4), 192. Completing the square gives saved minimum.
    z=vertices[0];assert dot(z,z)==196
    radial=q-4-(dot(x,z)-4)**2/F(192)
    assert radial==F(data['common_escape']['minimum_ball_powers']['all_real_vertex_radial'])
    triangle=[]
    for f,n in enumerate([7,7,8]):
        vv=[]
        for j in range(3):
            a=[0,0,0];a[f]=j;vv.append(vertex(a))
        z0,z1,z2=vv
        direction=[(a-2*b+c)/F(n*n) for a,b,c in zip(z0,z1,z2)]
        assert dot(direction,direction)==F(4,n*n)
        assert dot(z1,direction)==-F(2*(n-1),n)
        # w=z1+r*direction and R²=4r²/n²-r. Therefore the r²
        # term cancels identically, not only at sampled parameter values.
        linear=1-2*dot(sub(x,z1),direction)
        constant=dot(sub(x,z1),sub(x,z1))
        values=[constant+linear*r for r in [F(16),F(49,3)]]
        assert min(values)>1
        triangle.append(dict(nodes=n,squared_direction=str(dot(direction,direction)),
            power_constant=str(constant),power_linear=str(linear),
            endpoint_powers=list(map(str,values))))
    powers={k:F(v) for k,v in data['common_escape']['minimum_ball_powers'].items()}
    assert len(powers)==6 and min(powers.values())>1
    # Only the two extremal signs are required to isolate the positive
    # metric root; leading and linear coefficients are both negative.
    def sign5(a,b):
        if a*b>=0:return 1 if a+b>0 else -1
        return (1 if a>0 else -1)*(1 if a*a>5*b*b else -1)
    assert sign5(F(4),-F(4))==-1 and sign5(-F(539),F(49))==-1
    assert sign5(F(632),-F(240))==1  # polynomial at16
    r=F(49,3)
    assert sign5(4*r*r-539*r+8232,-4*r*r+49*r)==-1
    paths.add(str(Path(__file__).resolve().relative_to(ROOT)))
    return dict(status='PASS_ENDPOINT_REPLACEMENT_SOURCE_AND_SMALL_FRACTION_REVIEW',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)},
        escape_norm=1,minimum_seed_gap='1/100',radial_attainer_power=str(radial),
        exact_triangle_quadratic_cancellation=triangle,
        scope='Source/proof review plus fresh primary saved replay. Additional separate Fraction '
              'checks reconstruct the literal escape, all22 seed gaps, radial minimum attainer '
              'and exact triangle quadratic cancellation for both factor sizes. No second '
              'all-transport census, numerical metric audit, or chamber compatibility claim.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
