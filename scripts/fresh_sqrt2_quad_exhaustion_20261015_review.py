#!/usr/bin/env python3
"""Independent finite-classifier and recursive-radical quadrilateral audit.

No producer arithmetic imports: finite Dynkin recognition, independently
ordered union DFS, marked permutation orbits, direct exact radical inertia,
and literal finite checks supporting the separately audited all-size proof.
"""
import argparse
from collections import Counter,deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_sqrt2_quad_exhaustion_20261015'
TAG=BASE+'_review'


class K:
    def __init__(self,v=0):
        self.v=v.v if isinstance(v,K)else tuple(map(F,v))if isinstance(v,(list,tuple))else(F(v),)+(F(0),)*7
        assert len(self.v)==8
    def __add__(self,y):return K(tuple(a+b for a,b in zip(self.v,K(y).v)))
    __radd__=__add__
    def __neg__(self):return K(tuple(-a for a in self.v))
    def __sub__(self,y):return self+-K(y)
    def __rsub__(self,y):return K(y)+-self
    def __mul__(self,y):
        out=[F(0)]*8;y=K(y)
        for i,a in enumerate(self.v):
            if not a:continue
            for j,b in enumerate(y.v):
                if not b:continue
                shared=i&j;ix=i^j;c=a*b*(2 if shared&1 else 1)*(3 if shared&2 else 1)
                if shared&4:
                    out[ix]+=2*c;out[ix^1]+=c*(2 if ix&1 else 1)
                else:out[ix]+=c
        return K(out)
    __rmul__=__mul__
    def __truediv__(self,y):
        y=K(y);assert y
        cu=K(tuple((-a if i&4 else a)for i,a in enumerate(y.v)));nu=y*cu
        ct=K(tuple((-a if i&2 else a)for i,a in enumerate(nu.v)));nt=nu*ct
        cs=K(tuple((-a if i&1 else a)for i,a in enumerate(nt.v)));norm=nt*cs
        assert all(not a for a in norm.v[1:]) and norm.v[0]
        out=self*cu*ct*cs
        return K(tuple(a/norm.v[0]for a in out.v))
    def __bool__(self):return any(self.v)
    def __eq__(self,y):return self.v==K(y).v
    def sign(self):return sign_level(self,2)
    def enc(self):return [str(a)for a in self.v]


S=K((0,1,0,0,0,0,0,0));T=K((0,0,1,0,0,0,0,0));U=K((0,0,0,0,1,0,0,0))


def sign_level(x,level):
    if not x:return 0
    if level<0:
        assert not any(x.v[1:]);return (x.v[0]>0)-(x.v[0]<0)
    bit=1<<level
    a=K(tuple(x.v[i]if i<bit else 0 for i in range(8)))
    b=K(tuple(x.v[i+bit]if i<bit else 0 for i in range(8)))
    sa,sb=sign_level(a,level-1),sign_level(b,level-1)
    if not sa:return sb
    if not sb or sa==sb:return sa
    rad=[K(2),K(3),2+S][level]
    return sa*sign_level(a*a-rad*b*b,level-1)


def inertia(M):
    A=[[K(z)for z in row]for row in M];out=[0,0,0]
    while A:
        n=len(A);k=next((i for i in range(n)if A[i][i]),None)
        if k is None:
            edge=next(((i,j)for i,j in combinations(range(n),2)if A[i][j]),None)
            if edge is None:out[2]+=n;break
            i,j=edge
            for h in range(n):A[i][h]=A[i][h]+A[j][h]
            for h in range(n):A[h][i]=A[h][i]+A[h][j]
            k=i
        d=A[k][k];out[0 if d.sign()>0 else 1]+=1
        I=[i for i in range(n)if i!=k]
        A=[[A[i][j]-A[i][k]*A[k][j]/d for j in I]for i in I]
    return out


def pair(a,b):return tuple(sorted((a,b)))


def matrix(n,g,conjugate=False):
    values={2:K(),3:K(F(1,2)),4:(-S if conjugate else S)*F(1,2),6:T*F(1,2),8:((S-1)*U if conjugate else U)*F(1,2)}
    M=[[K(i==j)for j in range(n)]for i in range(n)]
    for i,j in combinations(range(n),2):M[i][j]=M[j][i]=-values[g.get((i,j),2)]
    return M


def finite(nodes,g):
    """Recognize the finite Dynkin forests, rather than eliminate the matrix."""
    adj={i:[]for i in nodes}
    for i,j in combinations(sorted(nodes),2):
        m=g.get((i,j),2)
        if m!=2:
            if m==0:return False
            adj[i].append(j);adj[j].append(i)
    unseen=set(nodes)
    while unseen:
        root=min(unseen);todo=[root];component={root}
        while todo:
            i=todo.pop()
            for j in adj[i]:
                if j not in component:component.add(j);todo.append(j)
        unseen-=component;n=len(component)
        edges=[(i,j,g[pair(i,j)])for i in component for j in adj[i]if i<j]
        if len(edges)!=n-1:return False
        if n<=2:continue
        if any(m>=6 for i,j,m in edges):return False
        four=[(i,j)for i,j,m in edges if m==4]
        if four:
            if len(four)!=1 or any(len(adj[i])>2 for i in component):return False
            i,j=four[0]
            if len(adj[i])==1 or len(adj[j])==1:continue
            if n==4:continue
            return False
        branches=[i for i in component if len(adj[i])>2]
        if not branches:continue
        if len(branches)!=1 or len(adj[branches[0]])!=3:return False
        center=branches[0];lengths=[]
        for first in adj[center]:
            prev,cur=center,first;length=1
            while len(adj[cur])==2:
                nxt=next(x for x in adj[cur]if x!=prev);prev,cur=cur,nxt;length+=1
            lengths.append(length)
        lengths.sort()
        if lengths[:2]==[1,1]:continue
        if lengths in[[1,2,2],[1,2,3],[1,2,4]]:continue
        return False
    return True


def corners(n,g):
    W=list(range(4,n))
    return all(finite(W+[i,j],g)for i,j in product((0,2),(1,3)))


def no_short(n,g):
    return all(g.get(pair(i,j),2)==2 and not any(g.get(pair(i,w),2)!=2 and g.get(pair(j,w),2)!=2 for w in range(4,n))for i,j in product((0,2),(1,3)))


def white_key(n,g):
    return min(tuple(g.get(pair(p[i],p[j]),2)for i,j in combinations(range(n),2))for tail in permutations(range(2,n))for p in[(0,1)+tail])


def independent_pieces(saved):
    cat=json.loads((ROOT/'results/fresh_sqrt2_lanner_catalogue_20261015.json').read_text())
    keys={(2,(0,))}
    for row in cat['order_three']['accepted']:
        keys.update((3,p)for p in permutations(row['labels']))
    for n in(4,5):
        template={pair(i,(i+1)%n):4 if i==0 else 3 for i in range(n)}
        for p in permutations(range(n)):
            g={pair(i,j):template.get(pair(p[i],p[j]),2)for i,j in combinations(range(n),2)}
            keys.add((n,white_key(n,g)))
    literal=[(r['order'],tuple(r['labels']))for r in saved]
    assert len(keys)==len(literal)==51 and set(literal)==keys
    # Input ordering is only a deterministic naming convention.
    return [(n,dict(zip(combinations(range(n),2),labels)))for n,labels in literal]


def identifications(left,right):
    for k in range(min(left,right)+1):
        for R in combinations(range(right),k):
            for L in permutations(range(4,4+left),k):
                d=dict(zip(R,L));nxt=4+left;mapping=[]
                for j in range(right):
                    if j in d:mapping.append(d[j])
                    else:mapping.append(nxt);nxt+=1
                yield tuple(mapping)


def union_enumeration(base):
    bank=independent_pieces(base['piece_catalogue']);result=set();count=Counter()
    for li,(ln,lg)in enumerate(bank):
        for rn,rg in bank[li:]:
            for rmw in identifications(ln-2,rn-2):
                count['identifications']+=1;lm=(0,2)+tuple(range(4,ln+2));rm=(1,3)+rmw;g={};consistent=True
                for old,mp in[(lg,lm),(rg,rm)]:
                    for (i,j),m in old.items():
                        e=pair(mp[i],mp[j])
                        if e in g and g[e]!=m:consistent=False;break
                        g[e]=m
                    if not consistent:break
                if not consistent:continue
                n=max(lm+rm)+1
                if not no_short(n,g)or not corners(n,g):continue
                for i,j in product((0,2),(1,3)):g[pair(i,j)]=2
                free=list(reversed([e for e in combinations(range(n),2)if e not in g]))
                def rec(k):
                    count['dfs_nodes']+=1
                    if k==len(free):
                        if not(set(lm)&set(rm))and not any(g.get(pair(i,j),2)!=2 for i,j in product(lm,rm)):
                            return
                        count['completions']+=1;result.add((n,tuple(g.get(e,2)for e in combinations(range(n),2))));return
                    e=free[k]
                    for m in(8,6,4,3,2):
                        g[e]=m
                        if no_short(n,g)and corners(n,g):rec(k+1)
                    del g[e]
                rec(0)
    expected={(r['nodes'],tuple(r['labels']))for r in base['labelled_union_supports']}
    assert count['identifications']==6010 and count['completions']==115
    assert result==expected and len(result)==115
    return dict(counts=dict(count),exact_literal_support_set_matches=True)


def cycle_field(n,g):
    values={3:K(1),4:S,6:T,8:U};adj={i:[]for i in range(n)}
    for i,j in combinations(range(n),2):
        if g.get((i,j),2)!=2:adj[i].append(j);adj[j].append(i)
    paths={};forest=set()
    for start in range(n):
        if start in paths:continue
        paths[start]=[start];todo=[start]
        while todo:
            i=todo.pop()
            for j in adj[i]:
                if j not in paths:paths[j]=paths[i]+[j];forest.add(pair(i,j));todo.append(j)
    for i,j in combinations(range(n),2):
        if g.get((i,j),2)==2 or(i,j)in forest:continue
        a,b=paths[i],paths[j];k=0
        while k<min(len(a),len(b))and a[k]==b[k]:k+=1
        route=list(reversed(a[k-1:]))+b[k:]
        value=values[g[(i,j)]]
        for u,v in zip(route,route[1:]):value=value*values[g[pair(u,v)]]
        if any(value.v[2:]):return False
    return True


def filters(base,cores,extensions):
    results=[];extra_results=[];counts=Counter()
    for record in cores['records']:
        row=base['labelled_union_supports'][record['union_index']];n=row['nodes'];g=dict(zip(combinations(range(n),2),row['labels']))
        if 0 in row['labels']:
            assert record['status']=='DOTTED_PARAMETER_UNFILTERED';continue
        if not cycle_field(n,g):
            assert record['status']=='REJECT_CYCLE_FIELD';counts['bad_cycle']+=1;continue
        phys=inertia(matrix(n,g));conj=inertia(matrix(n,g,True))
        status='REJECT_TWO_NEGATIVE_DIRECTIONS'if phys[1]>1 else'REJECT_CONJUGATE_NOT_PSD'if conj[1]else'PASS_PRINCIPAL_UNION_NECESSARY_TESTS'
        assert status==record['status'];counts[status]+=1
        results.append(dict(class_index=record['class_index'],physical_inertia=phys,conjugate_inertia=conj,status=status))
        if not status.startswith('PASS'):continue
        ex=next(r for r in extensions['records']if r['class_index']==record['class_index'])
        profiles=[];current=[2]*n
        def rec(k):
            if k==n:
                if any(m!=2 for m in current):profiles.append(tuple(current))
                return
            for m in(8,6,4,3,2):
                g[(k,n)]=m;current[k]=m
                if no_short(n+1,g)and corners(n+1,g):rec(k+1)
            del g[(k,n)]
        rec(0)
        assert set(profiles)=={tuple(r['contacts'])for r in ex['trials']}
        accepts=[];trials=[]
        for profile in profiles:
            for i,m in enumerate(profile):g[(i,n)]=m
            field=cycle_field(n+1,g)
            if field:
                ep=inertia(matrix(n+1,g));ec=inertia(matrix(n+1,g,True));accepted=ep[1]<=1 and ec[1]==0
            else:ep=ec=None;accepted=False
            saved=next(r for r in ex['trials']if tuple(r['contacts'])==profile)
            assert accepted==(saved['status']=='PASS_ONE_EXTRA_WHITE')
            if saved['status']=='REJECT_CONJUGATE_KERNEL':assert ec[1]>0
            if saved['status']=='REJECT_TWO_NEGATIVE_DIRECTIONS':assert ep[1]>1
            if saved['status']=='REJECT_CONJUGATE_NOT_PSD':assert ec[1]>0
            if saved['status']=='REJECT_CYCLE_FIELD':assert not field
            if accepted:accepts.append(profile)
            trials.append(dict(contacts=profile,cycle_field=field,physical_inertia=ep,conjugate_inertia=ec,accepted=accepted))
        assert set(accepts)==set(map(tuple,ex['accepted']))
        extra_results.append(dict(class_index=record['class_index'],structural_profiles=len(profiles),accepted=accepts,trials=trials))
    assert counts['REJECT_TWO_NEGATIVE_DIRECTIONS']==40 and counts['bad_cycle']==1
    assert counts['PASS_PRINCIPAL_UNION_NECESSARY_TESTS']==9
    assert sum(r['structural_profiles']for r in extra_results)==44
    return dict(core_counts=dict(counts),core_inertias=results,extra_white_profiles=extra_results)


BLACK_PERMS=[p for p in permutations(range(4))if {frozenset((p[0],p[2])),frozenset((p[1],p[3]))}=={frozenset((0,2)),frozenset((1,3))}]


def orbit_key(row):
    n=row['nodes'];g=dict(zip(combinations(range(n),2),row['labels']));edges=list(combinations(range(n),2))
    return n,min(tuple(g[pair(mp[i],mp[j])]for i,j in edges)for b in BLACK_PERMS for w in permutations(range(4,n))for mp in[b+w])


def check_map(left,right,mp):
    n=left['nodes'];assert right['nodes']==n and sorted(mp)==list(range(n))
    assert all((i<4)==(mp[i]<4)for i in range(n))
    assert {frozenset((mp[0],mp[2])),frozenset((mp[1],mp[3]))}=={frozenset((0,2)),frozenset((1,3))}
    a=dict(zip(combinations(range(n),2),left['labels']));b=dict(zip(combinations(range(n),2),right['labels']))
    assert all(m==b[pair(mp[i],mp[j])]for(i,j),m in a.items())


def marked_orbits(base):
    rows=base['labelled_union_supports'];groups={};keys=[]
    for i,row in enumerate(rows):
        key=orbit_key(row);keys.append(key);groups.setdefault(key,set()).add(i)
    saved=[];seen=set()
    for rec in base['marked_isomorphism_classes']:
        ids={m['index']for m in rec['members']};assert not(ids&seen);seen|=ids
        assert ids==groups[keys[rec['representative']]]
        for m in rec['members']:check_map(rows[m['index']],rows[rec['representative']],m['to_representative'])
        saved.append(sorted(ids))
    assert seen==set(range(115))and len(groups)==len(saved)==52
    return dict(independent_classes=52,all_labelled_rows_covered=115,full_pair_bijections_checked=115,class_members=saved)


def determinant(A):
    A=[[K(x)for x in r]for r in A];out=K(1);n=len(A)
    for k in range(n):
        pivot=next((i for i in range(k,n)if A[i][k]),None)
        if pivot is None:return K()
        if pivot!=k:A[k],A[pivot]=A[pivot],A[k];out=-out
        v=A[k][k];out=out*v
        for i in range(k+1,n):
            q=A[i][k]/v
            for j in range(k+1,n):A[i][j]=A[i][j]-q*A[k][j]
    return out


def dotted_profiles(base,ext):
    output=[];counts=Counter()
    for record in ext['dotted_records']:
        row=base['labelled_union_supports'][record['union_index']];n=row['nodes'];g=dict(zip(combinations(range(n),2),row['labels']))
        profiles=[];current=[2]*n
        def rec(k):
            if k==n:
                if any(v!=2 for v in current):profiles.append(tuple(current))
                return
            for m in(8,6,4,3,2):
                g[(k,n)]=m;current[k]=m
                if no_short(n+1,g)and corners(n+1,g):rec(k+1)
            del g[(k,n)]
        rec(0)
        assert set(profiles)=={tuple(r['contacts'])for r in record['trials']}
        checks=[]
        for trial in record['trials']:
            for i,m in enumerate(trial['contacts']):g[(i,n)]=m
            touched={i for i,m in enumerate(trial['contacts'])if m!=2};status=trial['status']
            if status=='ONLY_DOTTED_ENDPOINT_CONTACTS':assert touched<=set((0,2))
            elif status=='REJECT_AFFINE_433_STAR':
                support=trial['support'];gg={(i,j):g[pair(support[i],support[j])]for i,j in combinations(range(4),2)}
                M=matrix(4,gg);v=[K(1),S/2,K(F(1,2)),K(F(1,2))]
                assert all(x.sign()>0 for x in v)and all(not sum((M[i][j]*v[j]for j in range(4)),K())for i in range(4))
                assert inertia(M)==[3,0,1]
            elif status=='REJECT_ORTHOGONAL_INDEFINITE_PIECES':
                a=trial['dotted_pair'];b=trial['triangle']
                assert g[pair(*a)]==0 and all(g[pair(i,j)]==2 for i in a for j in b)
                gg={(i,j):g[pair(b[i],b[j])]for i,j in combinations(range(3),2)}
                assert determinant(matrix(3,gg))==-S/4
            else:raise AssertionError(status)
            checks.append(dict(contacts=trial['contacts'],status=status));counts[status]+=1
        output.append(dict(class_index=record['class_index'],checks=checks))
    assert sum(counts.values())==18
    return dict(counts=dict(counts),records=output)


def add_path(row,length,ends):
    n=row['nodes'];N=n+length;g=dict(zip(combinations(range(n),2),row['labels']))
    route=[ends[0]]+list(range(n,N))+list(ends[1:])
    for i,j in zip(route,route[1:]):g[pair(i,j)]=3
    return dict(nodes=N,white=list(range(4,N)),labels=[g.get(p,2)for p in combinations(range(N),2)])


def closure_checks(base,ext,closure):
    classes={i:base['labelled_union_supports'][r['representative']]for i,r in enumerate(base['marked_isomorphism_classes'])}
    byid={r['class_index']:r for r in ext['records']};sealed=[5,9,17,32,41,47,50]
    assert all(not byid[i]['accepted']for i in sealed)
    assert byid[14]['accepted']==[[2,2,3,2,2,2,2]]
    assert set(map(tuple,byid[37]['accepted']))=={(2,2,3,2,2,2,2),(3,2,2,2,2,2,2),(3,2,3,2,2,2,2)}
    checked=[]
    for rec in closure['exact_core_extension_maps']:
        left=add_path(classes[rec['source_core']],rec['path_length'],rec['ends']);right=classes[rec['target_core']]
        check_map(left,right,rec['bijection']);checked.append(rec)
    assert {(r['source_core'],r['path_length'],tuple(r['ends']),r['target_core'])for r in checked}=={(14,1,(2,),17),(37,1,(0,2),41),(37,2,(0,2),47),(37,3,(0,2),50)}
    row=classes[37];g=dict(zip(combinations(range(row['nodes']),2),row['labels']));evaluations=[]
    # Multilinearity in the two diagonal variables makes four evaluations
    # an exact independent polynomial identity check, not numerical sampling.
    for p,q in product((F(0),F(1)),repeat=2):
        A=matrix(row['nodes'],g);A[0][0]=K(p);A[2][2]=K(q)
        expected=(1+S)/64*(p+q-4*p*q+(1-p-q)*S)
        assert determinant(A)==expected
        evaluations.append(dict(p=str(p),q=str(q),determinant=expected.enc()))
    tails=(F(1),F(3,4),F(2,3),F(5,8));tailcases=[]
    for p,q in product(tails,repeat=2):
        value=p+q-4*p*q+(1-p-q)*S;assert value and p+q>1
        tailcases.append(dict(p=str(p),q=str(q),nonzero_rank_polynomial=value.enc()))
    assert (4-F(3,4)*(2+S)).sign()>0 and (2-(F(1,2)+S)).sign()>0
    # Reconstruct D5 and its short-arm continuations directly. The first
    # forbidden continuation is positive semidefinite affine E8.
    arm_checks=[]
    for ell in range(5):
        arms=[1,2,1+ell];gg={};n=1
        for length in arms:
            prev=0
            for _ in range(length):gg[(prev,n)]=3;prev=n;n+=1
        I=inertia(matrix(n,gg));assert I==([n,0,0]if ell<4 else[n-1,0,1])
        arm_checks.append(dict(extra_nodes=ell,arms=arms,inertia=I))
    return dict(sealed_cores=sealed,literal_maps=checked,core37_multilinear_identity=evaluations,all16_pendant_cases=tailcases,
                dotted_single_squared_bound=(F(3,4)*(2+S)).enc(),dotted_double_bound=(F(1,2)+S).enc(),D5_short_arm_extensions=arm_checks,
                all_size_proof='Source/proof reviewed: sealed first-contact obstruction, finite D/E one-branch classification, and Schur bounds for arbitrary connecting white components. No bound on total external whites assumed.')


def manifests(paths):
    seen=set()
    def visit(path):
        if path in seen:return
        seen.add(path)
        obj=json.loads((ROOT/path).read_text())
        for name,expected in obj.get('sources',{}).items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==expected,(path,name)
            if name.endswith('.json'):visit(name)
    for p in paths:
        if p.endswith('.json'):visit(p)
    return len(seen)


def run():
    base=json.loads((ROOT/'results'/f'{BASE}.json').read_text())
    cores=json.loads((ROOT/'results'/f'{BASE}_cores.json').read_text())
    extensions=json.loads((ROOT/'results'/f'{BASE}_extensions.json').read_text())
    closure=json.loads((ROOT/'results'/f'{BASE}_closure.json').read_text())
    assert S*S==2 and T*T==3 and U*U==2+S and (S-1)*U*((S-1)*U)==2-S
    unions=union_enumeration(base);print('independent union set PASS',flush=True)
    exact=filters(base,cores,extensions);print('independent direct inertia/filter PASS',flush=True)
    orbits=marked_orbits(base);print('independent marked orbits PASS',flush=True)
    dotted=dotted_profiles(base,extensions);closed=closure_checks(base,extensions,closure)
    paths=[f'scripts/{TAG}.py']+[p for suffix in('', '_cores','_extensions','_closure')for p in(f'scripts/{BASE}{suffix}.py',f'results/{BASE}{suffix}.json')]
    paths += [f'{folder}/fresh_sqrt2_dotted_tail_20261015{suffix}.{extension}'for suffix in('','_review')for folder,extension in(('scripts','py'),('results','json'))]
    paths += [f'{folder}/fresh_sqrt2_weight_bound_20261015{suffix}.{extension}'for suffix in('','_review')for folder,extension in(('scripts','py'),('results','json'))]
    count=manifests(paths)
    return dict(status='PASS_INDEPENDENT_UNION_FIELD_PROFILE_REPLAY_AND_ALL_SIZE_CLOSURE_PROOF_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                recursive_json_manifests_checked=count,independent_union_enumeration=unions,independent_exact_filters=exact,
                independent_marked_orbits=orbits,independent_dotted_profiles=dotted,closure_checks=closed,
                scope='Complete independent finite union and contact arithmetic; hostile all-size geometric closure proof review. Global weight arithmetic belongs to the separate weight-bound review; no claim for noncompact polytopes or other ground fields.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text());print('PASS saved replay')
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
