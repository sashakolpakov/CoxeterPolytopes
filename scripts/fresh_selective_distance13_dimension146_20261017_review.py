#!/usr/bin/env python3
"""Independent Fraction-only review of selective-distance13 dimension146.

Imports no producer arithmetic, SymPy, NumPy, or optimization package.
"""
import argparse
from collections import deque, Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import comb, prod
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_selective_distance13_dimension146_20261017_review'
BASE='fresh_selective_distance13_dimension146_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'


def sha(p):
    h=sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def matrix(n):return [[F(2 if i==j else 0)for j in range(n)]for i in range(n)]
def extract(a,rows,cols=None):
    cols=rows if cols is None else cols
    return [[a[i][j]for j in cols]for i in rows]
def transpose(a):return list(map(list,zip(*a)))
def multiply(a,b):return [[sum((x*y for x,y in zip(r,c)),F())for c in zip(*b)]for r in a]
def square(a,x):return sum((x[i]*a[i][j]*x[j]for i in range(len(x))for j in range(len(x))),F())
def encode(a):return [[str(x)for x in r]for r in a]


def reduced(a):
    a=[list(map(F,r))for r in a];pivots=[];i=0
    for j in range(len(a[0])):
        k=next((k for k in range(i,len(a))if a[k][j]),None)
        if k is None:continue
        a[i],a[k]=a[k],a[i];x=a[i][j];a[i]=[v/x for v in a[i]]
        for k in range(len(a)):
            if k!=i:
                x=a[k][j];a[k]=[v-x*w for v,w in zip(a[k],a[i])]
        pivots.append(j);i+=1
        if i==len(a):break
    return a,pivots


def inverse(a):
    n=len(a)
    b,p=reduced([r+[F(i==j)for j in range(n)]for i,r in enumerate(a)])
    assert p[:n]==list(range(n))
    return [r[n:]for r in b]


def nullspace(a):
    b,p=reduced(a);n=len(a[0]);free=[i for i in range(n)if i not in p];out=[]
    for j in free:
        x=[F(0)]*n;x[j]=F(1)
        for k,i in enumerate(p):x[i]=-b[k][j]
        out.append(x)
    return out


def psd(a,strict=False):
    a=[r[:]for r in a];p=[]
    while a:
        q=a[0][0];assert q>=0
        if q==0:
            assert not strict and all(v==0 for v in a[0]);p.append(q)
            a=[r[1:]for r in a[1:]]
        else:
            p.append(q)
            a=[[a[i][j]-a[i][0]*a[0][j]/q for j in range(1,len(a))]for i in range(1,len(a))]
    return p


def graph_dist(a,marks,connected=True):
    out={i:0 for i in marks};todo=deque(marks)
    while todo:
        i=todo.popleft()
        for j,v in enumerate(a[i]):
            if j!=i and v and j not in out:out[j]=out[i]+1;todo.append(j)
    if connected:assert len(out)==len(a)
    return [out.get(i)for i in range(len(a))]


def rational_gram(d):
    n=d['order'];adj=[[]for _ in range(n)]
    for i,j,m in d['edges']:
        assert m in(3,4);adj[i].append((j,m==4));adj[j].append((i,m==4))
    bits={0:0};todo=deque([0])
    while todo:
        i=todo.popleft()
        for j,flip in adj[i]:
            wanted=bits[i]^flip
            if j in bits:assert bits[j]==wanted
            else:bits[j]=wanted;todo.append(j)
    assert len(bits)==n
    norms=[1<<bits[i]for i in range(n)]
    a=[[F(2*norms[i]if i==j else 0)for j in range(n)]for i in range(n)]
    for i,j,m in d['edges']:
        a[i][j]=a[j][i]=-F(2 if m==4 or bits[i]else 1)
    return a,norms


def port_checks(ds,producer):
    rows=[];ports=[]
    for d in ds:
        adjacency=matrix(d['order'])
        for i,j,_ in d['edges']:adjacency[i][j]=adjacency[j][i]=F(-1)
        for marks in combinations(d['finite_deletions'],2):
            distances=graph_dist(adjacency,marks)
            rows.append(dict(diagram=d['index'],black=list(marks),distances=distances))
            if max(distances)<6:continue
            a,norms=rational_gram(d)
            for j,dist in enumerate(distances):
                if dist<6:continue
                vals=[]
                for black in marks:
                    ids=[i for i in range(len(a))if i!=black];b=extract(a,ids)
                    psd(b,True);inv=inverse(b);vals.append(norms[j]*inv[ids.index(j)][ids.index(j)])
                ports.append(dict(diagram=d['index'],black=list(marks),node=j,distance=dist,
                    inverse_diagonals=list(map(str,vals)),may_contact_external_white=max(vals)<2))
    assert len(rows)==384 and len(ports)==21
    pflat=[dict(diagram=q['diagram'],black=q['black'],**v)for q in producer['radius_contexts']for v in q['ports']]
    assert ports==pflat and len(producer['radius_contexts'])==16
    assert sum(x['may_contact_external_white']for x in ports)==7
    return dict(marked_pairs=len(rows),high_distance_contexts=16,high_distance_ports=ports)


def core_checks(d,producer):
    a,_=rational_gram(d);assert all(x.denominator==1 for r in a for x in r)
    g=matrix(19)
    for mp in [list(range(10)),list(range(10,19))+[9]]:
        for i,j,_ in d['edges']:g[mp[i]][mp[j]]=g[mp[j]][mp[i]]=F(-1)
    z=[F(x)for x in[-6,-3,-4,-2,-5,-4,-3,-2,-1,0,6,3,4,2,5,4,3,2,1]]
    assert multiply(g,[[x]for x in z])==[[0]for _ in z]
    assert len(nullspace(g))==1
    ids=[i for i in range(19)if i not in(0,10)];b=extract(g,ids);psd(b,True)
    x=extract(g,[0,10],ids);corr=multiply(multiply(x,inverse(b)),transpose(x))
    schur=[[g[i][j]-corr[k][l]for l,j in enumerate([0,10])]for k,i in enumerate([0,10])]
    assert schur==[[F(-1,12)]*2 for _ in range(2)]
    # Independent disjoint far-far and white far-near calculations.
    h=matrix(20)
    for off in [0,10]:
        for i,j,_ in d['edges']:h[i+off][j+off]=h[j+off][i+off]=F(-1)
    h[9][19]=h[19][9]=F(-1)
    ids2=[i for i in range(20)if i not in(0,10)];b=extract(h,ids2);psd(b,True)
    x=extract(h,[0,10],ids2);corr=multiply(multiply(x,inverse(b)),transpose(x))
    ss=[[h[i][j]-corr[k][l]for l,j in enumerate([0,10])]for k,i in enumerate([0,10])]
    assert ss==[[F(-7,78),F(-1,13)],[F(-1,13),F(-7,78)]]
    assert ss[0][0]+ss[0][1]==-F(1,6)and ss[0][0]-ss[0][1]==-F(1,78)
    h[9][19]=h[19][9]=F(0);h[9][18]=h[18][9]=F(-1)
    white=[18,19,17,16,9,8,7,6,5];b=extract(h,white);psd(b)
    ker=nullspace(b);assert len(ker)==1
    k=ker[0];k=k if all(t>0 for t in k)else[-t for t in k]
    assert all(t>0 for t in k)
    a3=[[F(2),-F(1),F(0)],[-F(1),F(2),-F(1)],[F(0),-F(1),F(2)]]
    assert 1/inverse(a3)[2][2]==F(4,3)and F(16,9)<2
    # Reconstruct all9 preferred corners, including every allowed extra contact.
    out=[]
    for left,right0 in product([(1,2),(1,3),(2,3)],repeat=2):
        right=tuple(i+10 for i in right0)
        p=2 if left!=(1,3)else 1;q=(2 if right0!=(1,3)else 1)+10
        omitted=sorted((set(left)|set(right))-{p,q});ids=[i for i in range(19)if i not in omitted]
        b=extract(g,ids);kernel=nullspace(b);psd(b,not kernel)
        if kernel:
            assert len(kernel)==1;v=kernel[0];assert all(t>0 for t in v)or all(t<0 for t in v)
            name='affine_D16';allowed=[]
        else:
            inv=inverse(b);assert all(t>=0 for r in inv for t in r)
            allowed=[ids[i]for i in range(17)if inv[i][i]<=2]
            name='A17'if left==right0==(1,2)else'D17'
        allowed=sorted(allowed+omitted);du=graph_dist(g,left);dv=graph_dist(g,right)
        assert all(z[i]and(du[i]<=1 if z[i]<0 else dv[i]<=1)for i in allowed)
        assert graph_dist(b,[ids.index(p)])[ids.index(q)]==14
        out.append(dict(left_black=list(left),right_black=list(right),preferred_corner=[p,q],
            possible_core_contacts=allowed,type=name))
    assert Counter(r['type']for r in out)=={'A17':1,'D17':4,'affine_D16':4}
    for x,y in zip(out,producer['supporting_white_isolation']['cases']):
        assert all(x[k]==y[k]for k in x)
    return g,z,dict(shared_radical=list(map(str,z)),shared_schur=encode(schur),
        disjoint_schur=encode(ss),far_near_positive_kernel=list(map(str,k)),preferred_corners=out)


def pyramid_checks(g,z,producer):
    ids=list(range(1,19));j=inverse(extract(g,ids));rays=[]
    for i in range(19):
        if not z[i]:
            p=[F(0)]*19;p[i]=-1;rays.append(p)
    for a in range(19):
        for b in range(19):
            if z[a]<0<z[b]:
                p=[F(0)]*19;p[a]=-z[b];p[b]=z[a];rays.append(p)
    assert len(rays)==82
    c=[sum(p[i]for p in rays)for i in range(19)];cc=[c[i]for i in ids]
    assert all(t<0 for t in c)and square(j,cc)<0
    counts=Counter();cusps=[]
    for p in rays:
        q=[p[i]for i in ids];norm=square(j,q);assert norm<=0
        assert sum((q[i]*j[i][k]*cc[k]for i in range(18)for k in range(18)),F())<0
        active=[i for i,t in enumerate(p)if not t];a=extract(g,active)
        psd(a,norm<0)
        if norm<0:counts['finite']+=1;assert len(active)==17
        else:
            counts['ideal']+=1
            todo=set(active);comps=[]
            while todo:
                comp={min(todo)};queue=deque(comp);todo-=comp
                while queue:
                    i=queue.popleft()
                    for k in list(todo):
                        if g[i][k]:comp.add(k);queue.append(k);todo.remove(k)
                block=extract(g,sorted(comp));kern=nullspace(block)
                assert len(kern)==1 and (all(t>0 for t in kern[0])or all(t<0 for t in kern[0]))
                comps.append(sorted(comp))
            assert sum(len(x)-1 for x in comps)==16
            cusps.append(dict(active=active,affine_component_nodes=comps))
    assert counts=={'finite':80,'ideal':2}
    assert sorted(sorted(len(c)for c in q['affine_component_nodes'])for q in cusps)==[[9,9],[17]]
    for i in range(19):
        mat=[[p[k]for k in ids]for p in rays if p[i]==0]
        assert len(reduced(mat)[1])==17
    assert dict(counts)==producer['known_sharp_distance14_pyramid']['vertices']
    quadrilaterals=[]
    for left,right0 in product([(1,2),(1,3),(2,3)],repeat=2):
        right=tuple(k+10 for k in right0);black=set(left+right)
        white=[k for k in range(19)if k not in black];psd(extract(g,white),True)
        face=[p for p in rays if all(p[k]==0 for k in white)]
        assert len(face)==4 and len(reduced([[p[k]for k in ids]for p in face])[1])==3
        distances=[]
        for u,v in product(left,right):
            corner=white+[u,v];b=extract(g,corner);psd(b)
            distance=graph_dist(b,[corner.index(u)],False)[corner.index(v)]
            assert distance is None or distance>=14;distances.append(distance)
        quadrilaterals.append(dict(left_black=list(left),right_black=list(right),
            vertices=4,corner_distances=distances))
    return dict(extreme_rays=82,finite=80,ideal=2,proper_cusp_links=cusps,
        actual_distance13_counterexample_faces=quadrilaterals,
        description='Complete cone p<=0,z.p=0; its rays are the zero-coefficient axis and every opposite-sign pair. Each ideal link has affine rank16 and positive relations, giving bounded products of alcoves.')


def polynomial(ranks):
    p=[1]
    for r in ranks:
        q=[0]*(len(p)+r)
        for i,t in enumerate(p):
            for j in range(r+1):q[i+j]+=t
        p=q
    return p


def profile(ranks):
    d=sum(ranks);h=polynomial(ranks);v=prod(r+1 for r in ranks)
    assert h==h[::-1]and sum(h)==v
    G=sum(t*min(j,d-j)for j,t in enumerate(h));curv=F(5*G-d*v,2*v)
    moment=0;s=len(ranks)
    for k in range(s+1):
        for subset in combinations(ranks,k):
            shift=sum(x+1 for x in subset);m=(d-1)//2-shift
            if m>=0:moment+=(-1)**k*((d-2*shift)*comb(m+s,s)-2*s*comb(m+s,s+1))
    assert curv==F(3*d,4)-F(5*moment,2*v)
    cost=sum((min(F(13),F(r,2))for r in ranks),F())
    return curv,cost,curv-cost


def budget_checks(producer):
    output=[];lookup={tuple(sorted(r['ranks'])):r for r in producer['reduced_profiles']}
    for d in range(146,156):
        shapes=set()
        # Independent normal-form enumeration: low in0..25, high0or>=27.
        for twentysixes in range(d//26+1):
            for low in range(26):
                high=d-26*twentysixes-low
                if high==0 or high>=27:
                    shapes.add(tuple(sorted([26]*twentysixes+([low]if low else[])+([high]if high else[]))))
        for ranks in sorted(shapes):
            a,b,c=profile(ranks);assert c>0
            q=lookup[ranks]
            assert (a,b,c)==tuple(F(q[k])for k in['normalized_curvature','normalized_cost','margin'])
            output.append(dict(ranks=list(ranks),margin=str(c)))
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>13*n
    assert len(output)==len(lookup)==1255
    assert min(F(x['margin'])for x in output)==F(830407,9565938)
    assert profile((26,)*5+(15,))[2]==-F(515185,19131876)
    # Exact parity interpolation, then degree2 algebraic verification at four points.
    polys=[]
    for initial,expected in[(111,[12,224,4]),(112,[120,228,4])]:
        vals=[]
        for t in range(4):
            d=initial+2*t;beta=F(d*(d-6),8)if d%2==0 else F((d-1)*(d-5),8)
            vals.append(8*(beta-13*(d+1)))
        coeff=[vals[0],vals[1]-vals[0]-(vals[2]-2*vals[1]+vals[0])/2,(vals[2]-2*vals[1]+vals[0])/2]
        assert coeff==expected and all(vals[t]==sum(coeff[i]*t**i for i in range(3))for t in range(4))
        polys.append([str(x)for x in coeff])
    return dict(profiles=output,minimum_margin='830407/9565938',compact_facet_parity_slack_times8=polys)


def run():
    producer=json.loads((ROOT/f'results/{BASE}.json').read_text())
    assert producer['status']=='PASS_SELECTIVE_DISTANCE13_SMALL_COMPONENT_EXCEPTION_DIMENSION146'
    for path,h in producer['sources'].items():assert sha(ROOT/path)==h,path
    ds=json.loads((ROOT/f'results/{CENSUS}.json').read_text())['diagrams']
    ports=port_checks(ds,producer);d=next(x for x in ds if x['index']==25)
    g,z,core=core_checks(d,producer);pyramid=pyramid_checks(g,z,producer)
    budgets=budget_checks(producer)
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',f'results/{BASE}.md',
        f'scripts/{CENSUS}.py',f'results/{CENSUS}.json',
        'results/fresh_cusp_factor_payment_framework_20261017.md',
        'results/fresh_cusp_link_shape_curvature_20261017.md',
        'results/fresh_distance14_unrestricted_dimension156_20261017.md']
    return dict(status='PASS_INDEPENDENT_FRACTION_SELECTIVE_DISTANCE13_DIMENSION146_REVIEW',
        sources={p:sha(ROOT/p)for p in paths},producer_source_bindings=len(producer['sources']),
        ports=ports,core=core,known_pyramid=pyramid,budgets=budgets,
        scope='Independent exact Fraction reconstruction plus separate mathematical proof/source read; inherits the26-type complete minimal-core classification, all-rank graph budgets, convex-order proof and earlier156 upper tail. No SymPy, shared arithmetic, optimization, formal verification or external review.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'384/21/9/82/1255',flush=True)
