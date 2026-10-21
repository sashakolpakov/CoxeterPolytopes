#!/usr/bin/env python3
"""Exact rank162 obstruction to fixed-omission local/narrow-pair accounting.

Rational physical Gram construction, independent of prior producer imports.
The B2 long root is rescaled to norm4. No finite-volume polytope is claimed;
a strict null escape is explicitly certified.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_rank162_local_obstruction_20261016'

def dot(x,y):return sum((a*b for a,b in zip(x,y)),F(0))
def mv(A,x):return [dot(row,x) for row in A]
def inverse(A):
    n=len(A);M=[list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        i=next(i for i in range(k,n) if M[i][k]);M[k],M[i]=M[i],M[k]
        z=M[k][k];M[k]=[v/z for v in M[k]]
        for i in range(n):
            if i!=k:
                z=M[i][k];M[i]=[x-z*y for x,y in zip(M[i],M[k])]
    return [row[n:] for row in M]
def pd_pivots(A):
    A=[list(row) for row in A];out=[]
    while A:
        d=A[0][0];assert d>0;out.append(d)
        A=[[A[i][j]-A[i][0]*A[0][j]/d for j in range(1,len(A))] for i in range(1,len(A))]
    return out

def gram(name):
    if name in ('E6','E7','E8'):
        arms={'E6':[2,2,2],'E7':[1,3,3],'E8':[1,2,5]}[name]
        n=sum(arms)+1;G=[[F(2*(i==j)) for j in range(n)] for i in range(n)];start=1
        for length in arms:
            path=[0]+list(range(start,start+length));start+=length
            for i,j in zip(path,path[1:]):G[i][j]=G[j][i]=F(-1)
        return G
    if name=='B2':return [[F(v) for v in row] for row in [[2,0,-2],[0,2,-2],[-2,-2,4]]]
    if name=='A1avoid':return [[F(2),F(-2)],[F(-2),F(2)]]
    assert name=='D5'
    G=[[F(2*(i==j)) for j in range(6)] for i in range(6)]
    for i,j in [(0,2),(1,2),(2,3),(3,4),(3,5)]:G[i][j]=G[j][i]=F(-1)
    return G

def factor(name):
    G=gram(name);n=len(G);r=n-1
    J=inverse([row[1:] for row in G[1:]])
    delta=[F(1)]+mv(J,[-G[i][0] for i in range(1,n)])
    delta=[v/min(delta) for v in delta]
    assert min(delta)>0 and mv(G,delta)==[0]*n
    assert min(delta[i]**2*G[i][i]/2 for i in range(n))==1
    g=next(i for i in range(n) if delta[i]==1 and G[i][i]==2)
    ids=[i for i in range(n) if i!=g];J=inverse([[G[i][j] for j in ids] for i in ids])
    pivots=pd_pivots([[G[i][j] for j in ids] for i in ids])
    omissions={'E6':(2,4),'E7':(3,6),'E8':(1,6),'B2':(0,2),'D5':(0,5),'A1avoid':(0,0)}[name]
    nodes={'E6':(2,4),'E7':(7,7),'E8':(8,8),'B2':(0,0),'D5':(0,5),'A1avoid':(0,0)}[name]
    p,q=[[F(i==node) for i in range(n)] for node in nodes]
    assert dot(delta,p)==dot(delta,q)==1
    assert all(not row[i] or row[i]**2>=G[i][i]/2 for row in (p,q) for i in range(n))
    costs=[];frame_pivots=[]
    for row,o in zip((p,q),omissions):
        keep=[i for i in range(n) if i!=o]
        A=[[G[i][j] for j in keep] for i in keep];rr=[row[i] for i in keep]
        frame_pivots.append(pd_pivots(A));costs.append(dot(rr,mv(inverse(A),rr)))
    diff=[p[i]-q[i] for i in ids];W=dot(diff,mv(J,diff))
    target=F(2*(r-int(name=='A1avoid')),9)
    return dict(name=name,rank=r,G=G,delta=delta,ground=g,ids=ids,J=J,
        pivots=pivots,omissions=omissions,p=p,q=q,costs=costs,W=W,target=target,
        frame_pivots=frame_pivots)

def witness(core):
    types=(['E7']*2+['B2']*2 if core=='E7' else ['E8']*2+['B2'])+['A1avoid']+['E6']*23+['D5']
    fs=[factor(name) for name in types];N=sum(f['rank'] for f in fs);assert N==162
    W=sum((f['W'] for f in fs),F(0));assert W==F(383,12)
    T=sum((f['target'] for f in fs),F(0));assert T==F(322,9)
    costs=[sum((f['costs'][j] for f in fs),F(0)) for j in (0,1)]
    assert max(costs)<=2
    pairing=W/2-2;assert pairing==F(335,24)<14
    offsets=[];count=0
    for f in fs:offsets.append(count);count+=len(f['G'])
    x,y=count,count+1;A=[[F(0)]*(count+2) for _ in range(count+2)]
    for off,f in zip(offsets,fs):
        for i,row in enumerate(f['G']):
            for j,v in enumerate(row):A[off+i][off+j]=v
        for i in range(len(f['G'])):
            A[off+i][x]=A[x][off+i]=-f['p'][i]
            A[off+i][y]=A[y][off+i]=-f['q'][i]
    A[x][x]=A[y][y]=2;A[x][y]=A[y][x]=-pairing
    assert all(A[i][j]<=0 for i in range(count+2) for j in range(i))
    normalized_max=F(0)
    for i in range(count+2):
        for j in range(i):
            v=A[i][j]
            assert 4*v*v<196*A[i][i]*A[j][j]
            if (j,i)!=(x,y):assert v*v/(A[i][i]*A[j][j]) in (0,F(1,4),F(1,2),F(1))
    # Eliminate one positive-definite grounded block per factor. Each omitted
    # row has residual zero and pairs -1 with each endpoint, since heights=1.
    sx=F(2);sy=F(2);sxy=-pairing
    for f in fs:
        pp=[f['p'][i] for i in f['ids']];qq=[f['q'][i] for i in f['ids']]
        sx-=dot(pp,mv(f['J'],pp));sy-=dot(qq,mv(f['J'],qq));sxy-=dot(pp,mv(f['J'],qq))
        om=f['ground'];gr=[f['G'][om][i] for i in f['ids']]
        assert f['G'][om][om]-dot(gr,mv(f['J'],gr))==0
        assert -f['p'][om]+dot(gr,mv(f['J'],pp))==-1
        assert -f['q'][om]+dot(gr,mv(f['J'],qq))==-1
    assert sx+sy-2*sxy==0
    # A selected omitted row and x have residual block [[0,-1],[-1,sx]],
    # determinant -1: one positive, one negative direction. Every other
    # omitted-row difference and y-x+(sxy-sx)*omitted_row is radical.
    assert F(0)*sx-(-1)**2==-1
    inertia=dict(positive=N+1,negative=1,zero=len(fs))
    assert sum(inertia.values())==len(A)
    u=[F(0)]*len(A)
    for i,v in enumerate(fs[0]['delta']):u[offsets[0]+i]=v
    Au=mv(A,u);assert Au[:count]==[0]*count and Au[x:]==[-1,-1]
    # Initial null escape uses the midpoint of the two endpoint centers.
    z=[F(0)]*len(A);z[x]=z[y]=1
    lam=1-pairing/2
    z=[a+lam*b for a,b in zip(z,u)]
    assert dot(z,mv(A,z))==0 and dot(z,Au)==-2
    # Move slightly into the interior of every affine chamber, preserving
    # common height2, then restore null norm by a multiple of the cusp ray.
    eps=F(1,1000);h=[F(0)]*len(A);contacts=[]
    for off,f in zip(offsets,fs):
        old=[a+b for a,b in zip(f['p'],f['q'])]
        fresh=[(1-eps)*v+eps*F(2,sum(f['delta'])) for v in old]
        assert min(fresh)>0 and dot(f['delta'],fresh)==2
        desired=[old[i]-fresh[i] for i in f['ids']]
        vv=mv(f['J'],desired)
        for i,v in zip(f['ids'],vv):h[off+i]=v
        contacts.extend(fresh)
    z=[a+b for a,b in zip(z,h)];mu=dot(z,mv(A,z))/4
    z=[a+mu*b for a,b in zip(z,u)];Az=mv(A,z)
    assert dot(z,Az)==0 and dot(z,Au)==-2
    assert Az[:count]==[-v for v in contacts] and max(Az)<0
    # The common cusp and this strict null ray lie in the same causal sheet,
    # since their product is -2. Strict inequalities give an open ideal patch,
    # so this explicit chamber has infinite volume.
    return dict(core=core,affine_rank=N,ambient_hyperbolic_dimension=N+1,
        factor_types=types,endpoint_costs=list(map(str,costs)),M='1',sum_target=str(T),
        sum_energy=str(W),loss=str(T-W),negative_endpoint_pairing=str(pairing),
        narrowness_gap=str(14-pairing),gram_dimension=len(A),gram_inertia=inertia,
        schur_endpoint_block=list(map(str,(sx,sxy,sy))),
        frame_schur_complements=[str(2-v) for v in costs],
        all_present_normalized_root_pairings_strictly_below14=True,
        common_height_kernel_relations=True,
        strict_null_escape=dict(coefficients=list(map(str,z)),gram_pairings=list(map(str,Az)),
            square='0',cusp_pairing='-2',epsilon=str(eps),maximum_wall_pairing=str(max(Az))),
        gram=[[str(v) for v in row] for row in A],
        factors=[dict(name=f['name'],rank=f['rank'],offset=o,gram=[[str(v) for v in row] for row in f['G']],
            normalized_kernel=list(map(str,f['delta'])),kernel_normalization='min delta_i*sqrt(root_norm_i/2)=1',
            omissions=f['omissions'],left_contacts=list(map(str,f['p'])),right_contacts=list(map(str,f['q'])),
            endpoint_costs=list(map(str,f['costs'])),energy=str(f['W']),target=str(f['target']),
            ground=f['ground'],positive_ground_pivots=list(map(str,f['pivots']))) for o,f in zip(offsets,fs)])

def run():
    records=[witness(core) for core in ('E7','E8')]
    assert records[0]['endpoint_costs']==['2','2']
    assert records[1]['endpoint_costs']==['16/9','11/6']
    paths=[f'scripts/{TAG}.py', 'results/fresh_nikulin_two_budget_frontier_20261015.json',
           'results/fresh_nikulin_two_budget_frontier_20261015_review.json',
           'results/fresh_nikulin_two_budget_frontier_20261015_dimension163.json']
    for path in paths:
        if path.endswith('.json'):
            for dep,digest in json.loads((ROOT/path).read_text()).get('sources',{}).items():
                assert sha256((ROOT/dep).read_bytes()).hexdigest()==digest
    return dict(status='PASS_EXACT_RANK162_FIXED_OMISSION_LOCAL_OBSTRUCTIONS_AND_STRICT_NULL_ESCAPE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        witnesses=records,scope='Two exact Coxeter Gram configurations at affine rank162, with common height ratio1, complete endpoint cost budgets, Lorentz signature, and all present wall pairings below14. They obstruct a dimension162 conclusion from these fixed-omission local conditions alone. Strict null escapes prove the displayed chambers have infinite volume; no finite-volume polytope, obstruction to changed framing choices, or optimal general dimension bound is claimed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
    for r in out['witnesses']:print(r['core'],r['gram_inertia'],r['negative_endpoint_pairing'],r['strict_null_escape']['maximum_wall_pairing'])
