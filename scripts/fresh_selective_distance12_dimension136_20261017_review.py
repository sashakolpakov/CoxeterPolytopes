#!/usr/bin/env python3
"""Independent Fraction-only review of selective-distance12 dimension136.

Imports no producer arithmetic, SymPy, NumPy, or optimization package.
"""
import argparse
from collections import deque, Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product, permutations
from math import comb, prod
from functools import lru_cache
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_selective_distance12_dimension136_20261017_review'
BASE='fresh_selective_distance12_dimension136_20261017'
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


def bond_graph(d):
    a=matrix(d['order'])
    for i,j,m in d['edges']:a[i][j]=a[j][i]=F(-(m-2))
    return a


def all_contexts(diagrams):
    out=[]
    for d in diagrams:
        g=bond_graph(d)
        for black in combinations(d['finite_deletions'],2):
            out.append(dict(diagram=d['index'],black=tuple(black),distances=graph_dist(g,black)))
    assert len(out)==384
    return out


def overlap_key(row):
    return (row['left'][0],tuple(row['left'][1]),row['right'][0],tuple(row['right'][1]),
            tuple(tuple(x)for x in row['match']))


def unite(a,b,match):
    mapping={j:i for i,j in match};n=len(a)
    for j in range(len(b)):
        if j not in mapping:mapping[j]=n;n+=1
    out=matrix(n)
    for i in range(len(a)):
        for j in range(len(a)):out[i][j]=a[i][j]
    for i in range(len(b)):
        for j in range(len(b)):
            if mapping[i]<len(a)and mapping[j]<len(a):assert out[mapping[i]][mapping[j]]==b[i][j]
            out[mapping[i]][mapping[j]]=b[i][j]
    return out,[mapping[j]for j in range(len(b))]


def positive_kernel(a):
    psd(a);ks=nullspace(a);assert len(ks)==1
    k=ks[0]
    if all(x<0 for x in k):k=[-x for x in k]
    assert all(x>0 for x in k)
    return k


def classification(diagrams,primary):
    ds={d['index']:d for d in diagrams};graphs={i:bond_graph(d)for i,d in ds.items()}
    contexts=all_contexts(diagrams)
    @lru_cache(None)
    def diagonals(di,black,node):
        a,scale=rational_gram(ds[di]);out=[]
        for b in black:
            ids=[i for i in range(len(a))if i!=b];g=extract(a,ids)
            psd(g,True);j=inverse(g);out.append(scale[node]*j[ids.index(node)][ids.index(node)])
        return tuple(out)
    ports=[]
    for r in contexts:
        for node,distance in enumerate(r['distances']):
            if distance<5:continue
            vals=diagonals(r['diagram'],r['black'],node)
            ports.append(dict(diagram=r['diagram'],black=list(r['black']),node=node,distance=distance,
                              inverse_diagonals=list(map(str,vals)),allowed=max(vals)<2))
    assert ports==primary['ports']
    allowed=[p for p in ports if p['allowed']];assert len(allowed)==18
    joins=[]
    for a,b in product(allowed,repeat=2):
        if a['distance']+b['distance']>=12:
            joins.append(dict(left=[a['diagram'],a['black'],a['node']],right=[b['diagram'],b['black'],b['node']]))
    assert joins==primary['disjoint_joins']and len(joins)==115
    types=Counter(tuple(sorted(((r['left'][0],r['left'][2]),(r['right'][0],r['right'][2]))))for r in joins)
    assert sorted(types.values())==[1,6,6,18,24,60]
    high=[r for r in contexts if max(r['distances'])>=6]
    assert len(high)==16 and all(sum(x>=6 for x in r['distances'])<=2 for r in high)
    saved={overlap_key(r):r for r in primary['overlap_records']}
    records=[];candidate_count=0
    for left,right in product(high,repeat=2):
        g,h=graphs[left['diagram']],graphs[right['diagram']]
        candidates=[(i,j)for i,d in enumerate(left['distances'])for j,e in enumerate(right['distances'])if d+e>=13]
        candidate_count+=len(candidates)+comb(len(candidates),2)
        # Independent matching enumeration: subsets and bijections, instead
        # of filtering combinations of the producer's candidate pairs.
        left_nodes=sorted({i for i,j in candidates});right_nodes=sorted({j for i,j in candidates})
        for size in (1,2):
            for ls,rs in product(combinations(left_nodes,size),combinations(right_nodes,size)):
                for pp in permutations(rs):
                    match=tuple(zip(ls,pp))
                    if any(x not in candidates for x in match):continue
                    if any(g[i][k]!=h[j][l]for i,j in match for k,l in match):continue
                    key=(left['diagram'],left['black'],right['diagram'],right['black'],match)
                    assert key in saved
                    old=saved[key];rejects=[]
                    for i,j in match:
                        for side,x,y,rr,other,shared in [('left',i,j,left,h,set(pp)),('right',j,i,right,g,set(ls))]:
                            external=[k for k in range(len(other))if k!=y and other[y][k]and k not in shared]
                            if external and max(diagonals(rr['diagram'],rr['black'],x))>=2:
                                assert all(k not in (right['black']if side=='left'else left['black'])for k in external)
                                rejects.append([side,x,external])
                    if rejects:
                        kind='inverse_rejection';assert kind in old
                    else:
                        assert left['diagram']==right['diagram']==25
                        a,_=rational_gram(ds[25]);u,mp=unite(a,a,match)
                        if match in (((8,9),),((9,8),)):
                            kind='white_affine_rejection';assert kind in old
                            white=old[kind]['nodes']
                            assert not set(white)&set(left['black']+tuple(mp[i]for i in right['black']))
                            k=positive_kernel(extract(u,white))
                            assert list(map(str,k))==old[kind]['positive_kernel']
                        else:
                            kind='shared_vertex'if len(match)==1 else'shared_edge'
                            assert match==(((9,9),)if kind=='shared_vertex'else((8,9),(9,8)))
                            assert old['surviving_type']==kind
                    records.append(dict(left=[left['diagram'],list(left['black'])],
                                        right=[right['diagram'],list(right['black'])],match=[list(x)for x in match],kind=kind))
    assert candidate_count==208 and len(records)==len(saved)==162
    counts=Counter(r['kind']for r in records)
    assert counts==primary['overlap_counts']=={'inverse_rejection':120,'white_affine_rejection':6,'shared_vertex':27,'shared_edge':9}
    good={((18,7),(25,9)),((22,4),(25,9)),((22,8),(25,9)),((25,9),(25,9))}
    unique={}
    for r in joins:
        a,b=r['left'],r['right'];typ=tuple(sorted(((a[0],a[2]),(b[0],b[2]))))
        if typ not in good:
            assert a[0]==b[0]==25
            native,_=rational_gram(ds[25]);g=disconnected_union(native,native)
            g[a[2]][10+b[2]]=g[10+b[2]][a[2]]=-F(1)
            if a[2]==b[2]==8:whites=[8,18,7,9,17,19]
            else:
                assert {a[2],b[2]}=={8,9}
                whites=[18,19,17,16,9,8,7,6,5]
                if a[2]==8:whites=[i+10 if i<10 else i-10 for i in whites]
            assert set(whites).isdisjoint(a[1]+[10+i for i in b[1]])
            positive_kernel(extract(g,whites))
            continue
        key=(a[0],tuple(a[1]),b[0],tuple(b[1]));unique.setdefault(key,[]).append((a[2],b[2]))
        for di,black,node in (a,b):
            gg=graphs[di];whites=[i for i in range(len(gg))if i not in black]
            assert any(len({node,x,y})==3 and gg[node][x]==gg[x][y]==-1 and not gg[node][y]
                       for x,y in product(whites,repeat=2))
    assert all(len(v)==1 for v in unique.values())
    return dict(marked_contexts=384,all_distance_at_least5_ports=ports,allowed_ports=18,
                disjoint_joins=115,disjoint_underlying_types=[dict(ports=k,count=v)for k,v in sorted(types.items())],
                candidate_matchings=208,consistent_matchings=162,overlap_counts=dict(counts),overlap_records=records)


def schur(g,centers):
    ids=[i for i in range(len(g))if i not in centers];a=extract(g,ids);psd(a,True)
    x=extract(g,centers,ids);c=multiply(multiply(x,inverse(a)),transpose(x))
    return [[g[i][j]-c[k][l]for l,j in enumerate(centers)]for k,i in enumerate(centers)]


def disconnected_union(a,b):
    g=matrix(len(a)+len(b))
    for off,c in [(0,a),(len(a),b)]:
        for i in range(len(c)):
            for j in range(len(c)):g[i+off][j+off]=c[i][j]
    return g


def join_checks(diagrams,primary):
    ds={d['index']:d for d in diagrams};a,_=rational_gram(ds[25]);rows=[]
    for di,v in [(18,7),(22,8),(22,4),(25,9)]:
        b,scales=rational_gram(ds[di]);assert set(scales)=={1}
        g=disconnected_union(a,b);g[9][10+v]=g[10+v][9]=-F(1)
        s=schur(g,[0,10]);det=s[0][0]*s[1][1]-s[0][1]*s[1][0]
        assert s[0][0]<0 and det>0
        old=next(r for r in primary['simple_leaf_joins']if r['right_diagram']==di and r['right_port']==v)
        assert encode(s)==old['schur']and str(det)==old['determinant']
        rows.append(dict(right_diagram=di,right_port=v,schur=encode(s),determinant=str(det)))
    kernels=[]
    for name,pair in [('far_near_white',(9,18)),('near_near_white',(8,18))]:
        g=disconnected_union(a,a);i,j=pair;g[i][j]=g[j][i]=-F(1)
        ids=primary[name]['nodes'];k=positive_kernel(extract(g,ids))
        assert list(map(str,k))==primary[name]['positive_kernel']
        kernels.append(dict(name=name,nodes=ids,positive_kernel=list(map(str,k))))
    a3=[[F(2),-F(1),F(0)],[-F(1),F(2),-F(1)],[F(0),-F(1),F(2)]]
    assert 1/inverse(a3)[0][0]==F(4,3)and F(16,9)<2
    return dict(negative_definite_schurs=rows,white_affine_kernels=kernels,
                cross_bond_squared_bound='16/9<2')


def pieces(g):
    rest=set(range(len(g)));out=[]
    while rest:
        seen={min(rest)};todo=list(seen);rest-=seen
        for i in todo:
            for j in sorted(rest):
                if g[i][j]:seen.add(j);todo.append(j);rest.remove(j)
        out.append(sorted(seen))
    return out


def preferred(black):
    return 0 if 0 in black else 2 if tuple(black)in((1,2),(2,3))else 1


def isolation_checks(diagrams,classification_data,primary):
    d=next(d for d in diagrams if d['index']==25);a,_=rational_gram(d)
    saved={(r['type'],tuple(r['left_black']),tuple(r['right_black'])):r for r in primary['cases']}
    cores=[];rows=[]
    for kind,match in [('shared_vertex',((9,9),)),('shared_edge',((8,9),(9,8)))]:
        g,mp=unite(a,a,match);assert encode(g)==primary['cores'][kind]['gram']
        s=schur(g,[0,10]);assert encode(s)==primary['cores'][kind]['schur']
        if kind=='shared_vertex':
            z=nullspace(g);assert len(z)==1;z=z[0]
            assert all(z[i]<0 for i in range(9))and z[9]==0 and all(z[i]>0 for i in range(10,19))
            assert s==[[-F(1,12)]*2 for _ in range(2)]
            assert list(map(str,z))==primary['cores'][kind]['radical']
            cores.append(dict(kind=kind,schur=encode(s),radical=list(map(str,z))))
        else:
            assert s==[[-F(5,66),-F(1,11)],[-F(1,11),-F(5,66)]]
            assert s[0][0]*s[1][1]-s[0][1]**2<0
            inv=inverse(g);end=[1,2,3];other=[mp[i]for i in end]
            ids=[i for i in range(len(g))if i not in(0,10)]
            determinant=prod(psd(extract(g,ids),True))*(s[0][0]*s[1][1]-s[0][1]*s[1][0])
            assert determinant==-1 and str(determinant)==primary['cores'][kind]['determinant']
            block=extract(inv,end);assert block==[[8,10,5],[10,14,7],[5,7,4]]
            assert extract(inv,other)==block and min(block[i][i]for i in range(3))==4
            cores.append(dict(kind=kind,schur=encode(s),determinant=str(determinant),end_inverse_block=encode(block),
                              nonzero_one_end_projection_cost_lower_bound='4'))
        contexts=[r for r in classification_data['overlap_records']if r['kind']==kind]
        for r in contexts:
            u=r['left'][1];v0=r['right'][1];v=[mp[i]for i in v0];p,q=preferred(u),mp[preferred(v0)]
            old=saved[kind,tuple(u),tuple(v)]
            omitted=sorted(set(u+v)-{p,q});ids=[i for i in range(len(g))if i not in omitted]
            h=extract(g,ids);allowed=omitted.copy();blocks=[];marked=None
            for indices in pieces(h):
                actual=[ids[i]for i in indices];b=extract(h,indices);ks=nullspace(b)
                if ks:
                    k=positive_kernel(b)
                    blocks.append(dict(nodes=actual,kind='affine',positive_kernel=list(map(str,k))))
                else:
                    psd(b,True);iv=inverse(b);assert all(x>=0 for rr in iv for x in rr)
                    vals=[iv[i][i]for i in range(len(b))]
                    allowed.extend(actual[i]for i,x in enumerate(vals)if x<=2)
                    blocks.append(dict(nodes=actual,kind='finite',inverse_diagonals=list(map(str,vals))))
                if p in actual:assert q in actual;marked=actual
            assert marked is not None and len(marked)<=17
            distance=graph_dist(h,[ids.index(p)],False)[ids.index(q)]
            assert distance>=13
            allowed.sort();assert allowed==old['possible_contacts']
            assert blocks==old['blocks']and marked==old['marked_component']and distance==old['marked_distance']
            du,dv=graph_dist(g,u),graph_dist(g,v)
            if kind=='shared_vertex':
                assert all(z[i]and(du[i]<=2 if z[i]<0 else dv[i]<=2)for i in allowed)
                bound=6
            else:
                assert set(allowed)<=set(end+other)
                assert all((du[i]if i in end else dv[i])<=1 for i in allowed)
                bound=4
            left_only=set(range(10))-{i for i,j in match}
            right_only={mp[i]for i in range(10)}-{mp[j]for i,j in match}
            extra=[]
            for i,j in product(sorted(left_only),sorted(right_only)):
                if g[i][j]:continue
                changed=[rr.copy()for rr in g];changed[i][j]=changed[j][i]=-F(1)
                if min(graph_dist(changed,u)[v0]for v0 in v)>12:
                    assert kind=='shared_vertex'and{i,j}=={8,mp[8]}
                    assert {8,9,mp[8]}.isdisjoint(u+v)
                    positive_kernel(extract(changed,[8,9,mp[8]]));extra.append([i,j])
            assert len(extra)==len(old['extra_cross_bonds_rejected_by_white_triangle'])
            rows.append(dict(type=kind,left_black=u,right_black=v,preferred_corner=[p,q],
                             blocks=blocks,marked_component=marked,marked_distance=distance,
                             allowed_core_contacts=allowed,shortcut_bound=bound,
                             extra_cross_bonds_excluded_by_affine_triangle=extra))
    assert len(rows)==36 and Counter(r['type']for r in rows)=={'shared_vertex':27,'shared_edge':9}
    return dict(cores=cores,preferred_corners=rows)


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
    cost=sum((min(F(12),F(r,2))for r in ranks),F())
    return curv,cost,curv-cost


def budget_checks(primary):
    output=[];lookup={tuple(sorted(r['ranks'])):r for r in primary['reduced_profiles']}
    for d in range(136,146):
        shapes=set()
        for capped in range(d//24+1):
            for low in range(24):
                high=d-24*capped-low
                if high==0 or high>=25:
                    shapes.add(tuple(sorted([24]*capped+([low]if low else[])+([high]if high else[]))))
        for ranks in sorted(shapes):
            a,b,c=profile(ranks);assert c>0;q=lookup[ranks]
            assert (a,b,c)==tuple(F(q[k])for k in['normalized_curvature','normalized_cost','margin'])
            output.append(dict(ranks=list(ranks),margin=str(c)))
        n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>12*n
    assert len(output)==len(lookup)==1175
    assert min(F(x['margin'])for x in output)==F(172096,1953125)
    assert profile((24,)*5+(15,))[2]==-F(59157,3906250)
    polynomials=[]
    for initial,expected in[(103,[12,208,4]),(104,[112,212,4])]:
        vals=[]
        for t in range(4):
            d=initial+2*t;beta=F(d*(d-6),8)if d%2==0 else F((d-1)*(d-5),8)
            vals.append(8*(beta-12*(d+1)))
        coefficients=[vals[0],vals[1]-vals[0]-(vals[2]-2*vals[1]+vals[0])/2,(vals[2]-2*vals[1]+vals[0])/2]
        assert coefficients==expected and all(vals[t]==sum(coefficients[i]*t**i for i in range(3))for t in range(4))
        polynomials.append(list(map(str,coefficients)))
    return dict(profiles=output,minimum_margin='172096/1953125',boundary_margin='-59157/3906250',
                compact_facet_parity_slack_times8=polynomials)


def run():
    primary=json.loads((ROOT/f'results/{BASE}.json').read_text())
    assert primary['status']=='PASS_SELECTIVE_DISTANCE12_SMALL_COMPONENT_EXCEPTION_DIMENSION136'
    for path,h in primary['sources'].items():assert sha(ROOT/path)==h,path
    diagrams=json.loads((ROOT/f'results/{CENSUS}.json').read_text())['diagrams']
    split=classification(diagrams,primary['census'])
    joins=join_checks(diagrams,primary['disjoint_obstructions'])
    iso=isolation_checks(diagrams,split,primary['isolation'])
    budgets=budget_checks(primary)
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',f'results/{BASE}.md',
           f'scripts/{CENSUS}.py',f'results/{CENSUS}.json',
           'results/fresh_cusp_factor_payment_framework_20261017.md',
           'results/fresh_cusp_link_shape_curvature_20261017.md',
           'results/fresh_selective_distance13_dimension146_20261017.md',
           'results/fresh_selective_distance13_dimension146_20261017_review.md']
    return dict(status='PASS_INDEPENDENT_FRACTION_SELECTIVE_DISTANCE12_DIMENSION136_REVIEW',
                sources={p:sha(ROOT/p)for p in paths},producer_bindings=len(primary['sources']),
                classification=split,disjoint_joins=joins,isolation=iso,budgets=budgets,
                scope='Standalone Fraction computations and independent matching enumeration, plus written geometric/source review. No producer arithmetic imports, SymPy, optimizer or external referee/formal proof. Inherited complete minimal-core classification, incidence, convex-order and all-rank budget theorems remain explicit proof inputs.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    out=json.loads(json.dumps(run()));target=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(target.read_text())
    else:target.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'384contexts/18ports/162overlaps/36corners/4Schurs/1175profiles',flush=True)

