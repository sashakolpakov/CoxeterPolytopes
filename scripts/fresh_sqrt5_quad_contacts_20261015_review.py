#!/usr/bin/env python3
"""Independent Dynkin-profile and full-inertia review of Qsqrt5 contacts.

Standalone Fraction arithmetic in basis1,sqrt5. No producer arithmetic or
graph imports. Uses tree-product scaling, not the producer's minimal colors.
"""
import argparse,json
from collections import Counter,deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_sqrt5_quad_exhaustion_20261015'
PRIMARY='fresh_sqrt5_quad_contacts_20261015'
TAG=PRIMARY+'_review'


class K:
    def __init__(self,a=0,b=0):
        if isinstance(a,K):self.a,self.b=a.a,a.b
        else:self.a,self.b=F(a),F(b)
    def __add__(self,y):y=K(y);return K(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return K(-self.a,-self.b)
    def __sub__(self,y):return self+-K(y)
    def __rsub__(self,y):return K(y)+-self
    def __mul__(self,y):y=K(y);return K(self.a*y.a+5*self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def __truediv__(self,y):
        y=K(y);n=y.a*y.a-5*y.b*y.b;assert n
        return self*K(y.a/n,-y.b/n)
    def __pow__(self,n):
        assert n>=0;z=K(1)
        while n:
            if n&1:z=z*self
            self=self*self;n//=2
        return z
    def __eq__(self,y):y=K(y);return self.a==y.a and self.b==y.b
    def __bool__(self):return bool(self.a or self.b)
    def conj(self):return K(self.a,-self.b)
    def sign(self):
        a,b=self.a,self.b;sg=lambda x:(x>0)-(x<0)
        if not a:return sg(b)
        if not b or sg(a)==sg(b):return sg(a)
        return sg(a)*sg(a*a-5*b*b)
    def enc(self):return [str(self.a),str(self.b)]


PHI=K(F(1,2),F(1,2));RADS=(K(2),K(3),2+PHI)


def pair(i,j):return tuple(sorted((i,j)))
def edges(row):return dict(zip(combinations(range(row['nodes']),2),row['labels']))
def restrict(M,ids):return [[M[i][j]for j in ids]for i in ids]


def inertia(M):
    A=[[K(x)for x in row]for row in M];out=[0,0,0]
    while A:
        n=len(A);pivot=next((i for i in range(n)if A[i][i]),None)
        if pivot is None:
            e=next(((i,j)for i,j in combinations(range(n),2)if A[i][j]),None)
            if e is None:out[2]+=n;break
            i,j=e
            for h in range(n):A[i][h]=A[i][h]+A[j][h]
            for h in range(n):A[h][i]=A[h][i]+A[h][j]
            pivot=i
        d=A[pivot][pivot];out[0 if d.sign()>0 else 1]+=1
        ids=[i for i in range(n)if i!=pivot]
        A=[[A[i][j]-A[i][pivot]*A[pivot][j]/d for j in ids]for i in ids]
    return out


def kernel(M):
    A=[[K(x)for x in r]for r in M];n=len(A);piv=[];r=0
    for col in range(n):
        p=next((i for i in range(r,n)if A[i][col]),None)
        if p is None:continue
        A[r],A[p]=A[p],A[r];z=A[r][col];A[r]=[v/z for v in A[r]]
        for i in range(n):
            if i!=r:
                q=A[i][col];A[i]=[v-q*w for v,w in zip(A[i],A[r])]
        piv.append(col);r+=1
    out=[]
    for f in range(n):
        if f in piv:continue
        v=[K(i==f)for i in range(n)]
        for i,p in enumerate(piv):v[p]=-A[i][f]
        assert all(not sum((a*b for a,b in zip(row,v)),K())for row in M)
        out.append(v)
    return out


def tree_scaled(n,g):
    """Formal products of edge magnitudes along a separately chosen DFS tree."""
    token={3:(K(1),(0,0,0)),4:(K(1),(1,0,0)),5:(PHI,(0,0,0)),6:(K(1),(0,1,0)),10:(K(1),(0,0,1))}
    paths={};adj={i:[]for i in range(n)}
    for i,j in combinations(range(n),2):
        lab=g.get((i,j),2);assert lab!=0
        if lab!=2:adj[i].append(j);adj[j].append(i)
    for start in reversed(range(n)):
        if start in paths:continue
        paths[start]=(K(1),(0,0,0));stack=[start]
        while stack:
            i=stack.pop();factor,exponents=paths[i]
            for j in sorted(adj[i],reverse=True):
                if j not in paths:
                    f,e=token[g[pair(i,j)]];paths[j]=(factor*f,tuple(a+b for a,b in zip(exponents,e)));stack.append(j)
    M=[[K()for _ in range(n)]for _ in range(n)]
    for i,(factor,e)in paths.items():
        z=2*factor*factor
        for d,power in zip(RADS,e):z=z*d**power
        assert z.sign()>0 and z.conj().sign()>0;M[i][i]=z
    for i,j in combinations(range(n),2):
        lab=g.get((i,j),2)
        if lab==2:continue
        f,e=token[lab];fi,ei=paths[i];fj,ej=paths[j]
        totals=[a+b+c for a,b,c in zip(ei,ej,e)]
        if any(v%2 for v in totals):return None
        z=-fi*fj*f
        for d,power in zip(RADS,totals):z=z*d**(power//2)
        M[i][j]=M[j][i]=z
        square={3:K(1),4:K(2),5:PHI*PHI,6:K(3),10:2+PHI}[lab]
        assert 4*z*z/(M[i][i]*M[j][j])==square
    return M


def graph(nodes,g):
    adj={i:[]for i in nodes}
    for i,j in combinations(sorted(nodes),2):
        if g.get((i,j),2)!=2:adj[i].append(j);adj[j].append(i)
    return adj


def components(adj):
    unseen=set(adj);out=[]
    while unseen:
        todo=[min(unseen)];found=set()
        while todo:
            i=todo.pop()
            if i in found:continue
            found.add(i);todo+=adj[i]
        unseen-=found;out.append(sorted(found))
    return out


def distance(adj,a,b):
    todo=deque([(a,0)]);seen={a}
    while todo:
        i,d=todo.popleft()
        if i==b:return d
        for j in adj[i]:
            if j not in seen:seen.add(j);todo.append((j,d+1))
    return None


@lru_cache(None)
def finite_word(n,labels):
    g=dict(zip(combinations(range(n),2),labels));adj=graph(range(n),g)
    for comp in components(adj):
        es=[(i,j,g[(i,j)])for i,j in combinations(comp,2)if g[(i,j)]!=2];r=len(comp)
        if len(es)!=r-1 or any(m==0 for _,_,m in es):return False
        if r<=2:continue
        if any(m>=6 for _,_,m in es):return False
        special=[(i,j,m)for i,j,m in es if m!=3]
        if special:
            if len(special)!=1 or any(len(adj[i])>2 for i in comp):return False
            i,j,m=special[0]
            if m==5:
                if r not in(3,4)or min(len(adj[i]),len(adj[j]))!=1:return False
            elif m==4:
                if min(len(adj[i]),len(adj[j]))!=1 and r!=4:return False
            else:return False
            continue
        branches=[i for i in comp if len(adj[i])>2]
        if not branches:continue
        if len(branches)!=1 or len(adj[branches[0]])!=3:return False
        lengths=[];c=branches[0]
        for j in adj[c]:
            prev,cur=c,j;length=1
            while len(adj[cur])==2:
                nxt=next(v for v in adj[cur]if v!=prev);prev,cur=cur,nxt;length+=1
            lengths.append(length)
        lengths.sort()
        if lengths[:2]!=[1,1]and lengths not in([1,2,2],[1,2,3],[1,2,4]):return False
    return True


def finite(ids,g):
    ids=sorted(ids);return finite_word(len(ids),tuple(g.get(pair(i,j),2)for i,j in combinations(ids,2)))


def no_short(n,g):
    return all(g.get(pair(a,b),2)==2 and not any(g.get(pair(a,w),2)!=2 and g.get(pair(b,w),2)!=2 for w in range(4,n))for a,b in product((0,2),(1,3)))


def domain(row):
    n=row['nodes'];g=edges(row);current=[2]*n;out=[]
    def rec(k):
        if k==n:
            if any(v!=2 for v in current):out.append(tuple(current))
            return
        for label in(10,6,5,4,3,2):
            g[(k,n)]=label;current[k]=label
            if no_short(n+1,g)and all(finite(list(range(4,n+1))+[a,b],g)for a,b in product((0,2),(1,3))):rec(k+1)
        del g[(k,n)]
    rec(0);return out


def extension(row,contact):
    n=row['nodes'];g=edges(row)
    for i,m in enumerate(contact):g[(i,n)]=m
    return dict(nodes=n+1,white=list(range(4,n+1)),labels=[g[p]for p in combinations(range(n+1),2)])


def paid(row):
    n=row['nodes'];g=edges(row);out=[]
    for a,b in product((0,2),(1,3)):
        adj=graph(row['white']+[a,b],g)
        for comp in components(adj):
            if a not in comp or b not in comp or len(comp)!=4:continue
            es=[(i,j,g[(i,j)])for i,j in combinations(comp,2)if g[(i,j)]!=2]
            if len(es)!=3 or sorted(m for _,_,m in es)!=[3,3,5]:continue
            if sorted(len(adj[i])for i in comp)!=[1,1,2,2]or len(adj[a])!=1 or len(adj[b])!=1:continue
            edge5=next((i,j)for i,j,m in es if m==5)
            assert min(len(adj[i])for i in edge5)==1
            out.append((tuple(sorted((a,b))),tuple(comp)))
    return set(out)


@lru_cache(None)
def affine_kernel(labels):
    n=next(n for n in range(3,20)if n*(n-1)//2==len(labels));g=dict(zip(combinations(range(n),2),labels))
    M=tree_scaled(n,g);assert M is not None and inertia(M)==[n-1,0,1]
    adj=graph(range(n),g);assert len(components(adj))==1
    vs=kernel(M);assert len(vs)==1;v=vs[0]
    if v[0].sign()<0:v=[-a for a in v]
    assert all(a.sign()>0 for a in v)
    return tuple(tuple(x.enc())for x in v)


def check_affine(row,witness):
    ids=witness['nodes'];g=edges(row);labels=tuple(g[pair(i,j)]for i,j in combinations(ids,2))
    assert list(labels)==witness['labels'];return affine_kernel(labels)


def check_paid(row,saved):
    actual=paid(row);reported={(tuple(sorted(r['marked'])),tuple(r['entire_H4']))for r in saved}
    assert actual==reported and actual
    return [dict(marked=list(p),entire_H4=list(c))for p,c in sorted(actual)]


def manifest_check(paths):
    seen=set()
    def visit(p):
        if p in seen:return
        seen.add(p);data=json.loads((ROOT/p).read_text())
        for name,h in data.get('sources',{}).items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(p,name)
            if name.endswith('.json'):visit(name)
    for p in paths:
        if p.endswith('.json'):visit(p)
    return len(seen)


def run():
    base=json.loads((ROOT/f'results/{BASE}.json').read_text());primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    filters=json.loads((ROOT/'results/fresh_sqrt5_quad_filters_20261015.json').read_text())
    expected={r['class_index']:r for r in filters['records']if r['status']in('PASS_PRINCIPAL_UNION_NECESSARY_TESTS','DOTTED_PARAMETER_UNFILTERED')}
    assert set(expected)=={r['class_index']for r in primary['records']}
    records=[];counts=Counter();trials=Counter();sealed=Counter();domain_count=0
    for rec in primary['records']:
        row=base['labelled_union_supports'][rec['union_index']];n=row['nodes'];g=edges(row);status=rec['status'];counts[status]+=1
        out=dict(class_index=rec['class_index'],union_index=rec['union_index'],status=status)
        assert rec['union_index']==expected[rec['class_index']]['union_index']
        if status=='ALREADY_PAID_BY_ENTIRE_H4':out['payments']=check_paid(row,rec['H4_payments'])
        elif status=='REJECT_AFFINE_PRINCIPAL_SUBDIAGRAM':out['positive_kernel']=check_affine(row,rec['affine_witness'])
        else:
            assert not paid(row)
            allrows=domain(row);assert set(allrows)=={tuple(r['contacts'])for r in rec['trials']};domain_count+=len(allrows)
            dotted=rec['dotted'];accepted=[];checks=[]
            if not dotted:
                M=tree_scaled(n,g);assert M is not None
                I=inertia(M);J=inertia([[z.conj()for z in r]for r in M])
                assert I[1]==1 and J[1]==0 and (I[2]>0)==rec['singular']
                out['core_physical_inertia']=I;out['core_conjugate_inertia']=J
            for trial in rec['trials']:
                contact=trial['contacts'];ext=extension(row,contact);label=trial['status'];check=dict(contacts=contact,status=label);trials[label]+=1
                if label=='BECOMES_PAID_BY_ENTIRE_H4':check['payments']=check_paid(ext,trial['H4_payments'])
                elif label=='REJECT_AFFINE_PRINCIPAL_SUBDIAGRAM':check['positive_kernel']=check_affine(ext,trial['affine_witness'])
                elif dotted:
                    assert label=='DOTTED_PARAMETER_UNFILTERED'and not paid(ext);accepted.append(tuple(contact))
                else:
                    assert not paid(ext);N=tree_scaled(n+1,edges(ext))
                    if label=='REJECT_CYCLE_FIELD':assert N is None
                    else:
                        assert N is not None
                        I=inertia(N);J=inertia([[x.conj()for x in r]for r in N])
                        check['physical_inertia']=I;check['conjugate_inertia']=J
                        if label=='REJECT_TWO_NEGATIVE_DIRECTIONS':assert I[1]>1
                        elif label in('REJECT_CONJUGATE_NOT_PSD','REJECT_CONJUGATE_KERNEL'):assert J[1]>0
                        elif label=='PASS_EXTRA_WHITE_NECESSARY_TESTS':assert I[1]<=1 and J[1]==0;accepted.append(tuple(contact))
                        else:raise AssertionError(label)
                checks.append(check)
            assert set(accepted)==set(map(tuple,rec['accepted']))
            if not dotted and not accepted:sealed['singular'if rec['singular']else'nonsingular']+=1
            if not dotted and accepted:sealed['nonsingular_open']+=1;assert not rec['singular']
            out.update(structural_profiles=len(allrows),accepted=[list(c)for c in accepted],trials=checks)
        records.append(out)
        if len(records)%50==0:print('independent records',len(records),'profiles',domain_count,flush=True)
    assert dict(counts)==primary['counts']and dict(trials)==primary['profile_counts']
    assert domain_count==4811 and sealed==dict(singular=112,nonsingular=44,nonsingular_open=130)
    paths=[f'scripts/{TAG}.py']
    for stem in(PRIMARY,BASE,'fresh_sqrt5_quad_filters_20261015',BASE+'_incidence_review','fresh_sqrt5_quad_filters_20261015_review'):
        paths += [f'scripts/{stem}.py',f'results/{stem}.json']
    manifests=manifest_check(paths)
    return json.loads(json.dumps(dict(status='PASS_INDEPENDENT_COMPLETE_ONE_WHITE_DOMAINS_AND_EXACT_QSQRT5_CONTACT_FILTERS',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},recursive_json_manifests=manifests,
                core_counts=dict(counts),profile_counts=dict(trials),total_profiles=domain_count,unpaid_sealing_counts=dict(sealed),records=records,
                arithmetic='Independent Fraction basis1,sqrt5; DFS tree-product radical scaling; direct full congruence inertia. No producer imports.',
                scope='Complete one-extra-white domains conditional on the separately reviewed full union list. Stable H4 payments and literal positive affine kernels. No complete dotted/all-external-white closure or global dimension bound.')))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text());print('PASS saved review replay')
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['core_counts'],out['unpaid_sealing_counts'])
