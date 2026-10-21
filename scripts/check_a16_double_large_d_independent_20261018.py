#!/usr/bin/env python3
"""Small exact check of coupled A16 certificates; no mathematical producer imports.

The check reconstructs the contact domain, core linear algebra, native roles,
all payment sums, and complete small-DP domain. The universal use of these
identities remains the argument in the accompanying written proof.
"""
import fcntl
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import resource
import signal
import run_lean_batches_20261018 as R

ROOT=Path(__file__).resolve().parents[1]
UP='fresh_a16_double_large_d_upper_scout_20261018'
LOCAL='fresh_a16_double_large_d_physical_scout_20261018'
TARGET={(0,13),(2,15),(1,13),(2,14),(1,12),(3,14),(2,13)}


def load(stem):
    p=ROOT/f'results/{stem}.json';assert p.stat().st_size<20*1024**2
    with p.open()as stream:return json.load(stream)


def digest(p):
    h=sha256()
    with p.open('rb')as stream:
        for chunk in iter(lambda:stream.read(1048576),b''):h.update(chunk)
    return h.hexdigest()


def inverse(g):
    n=len(g);a=[list(map(F,row))+[F(i==j)for j in range(n)]for i,row in enumerate(g)]
    for i in range(n):
        k=next(k for k in range(i,n)if a[k][i]);a[i],a[k]=a[k],a[i]
        p=a[i][i];a[i]=[v/p for v in a[i]]
        for j in range(n):
            if j!=i:
                p=a[j][i];a[j]=[v-p*w for v,w in zip(a[j],a[i])]
    return[r[n:]for r in a]


def inertia(g):
    a=[list(map(F,row))for row in g];p=m=z=0
    while a:
        n=len(a);i=next((i for i in range(n)if a[i][i]),None)
        if i is not None:
            ids=[j for j in range(n)if j!=i];q=a[i][i];p+=q>0;m+=q<0
            a=[[a[j][k]-a[j][i]*a[i][k]/q for k in ids]for j in ids]
        else:
            pair=next(((i,j)for i in range(n)for j in range(i)if a[i][j]),None)
            if pair is None:z+=n;break
            i,j=pair;q=a[i][j];ids=[k for k in range(n)if k not in pair]
            a=[[a[k][l]-(a[k][i]*a[j][l]+a[k][j]*a[i][l])/q for l in ids]for k in ids];p+=1;m+=1
    return[p,m,z]


def sub(g,ids):return[[g[i][j]for j in ids]for i in ids]


def blocks(g,ids):
    rest=set(ids);out=[]
    while rest:
        todo=[min(rest)];rest.remove(todo[0])
        for i in todo:
            for j in sorted(rest):
                if g[i][j]:rest.remove(j);todo.append(j)
        ns=sorted(todo);adj={i:[j for j in ns if j!=i and g[i][j]]for i in ns}
        assert all(g[i][j]==-1 for i in ns for j in adj[i])
        if max(map(len,adj.values()),default=0)<=2:
            assert sum(map(len,adj.values()))==2*(len(ns)-1);kind='A'
            start=min(i for i in ns if len(adj[i])<=1);path=[start]
            for i in path:
                path.extend(j for j in adj[i]if j not in path)
            maps=[path,path[::-1]]
        else:
            assert len(ns)>=4 and sorted(map(len,adj.values()))==[1,1,1]+[2]*(len(ns)-4)+[3]
            fork=next(i for i in ns if len(adj[i])==3);leaves=[i for i in adj[fork]if len(adj[i])==1]
            assert len(leaves)==2;kind='D';path=[fork];prev=None
            while len(path)<len(ns)-2:
                nxt=next(j for j in adj[path[-1]]if j!=prev and j not in leaves)
                prev=path[-1];path.append(nxt)
            maps=[path[::-1]+leaves,path[::-1]+leaves[::-1]]
        out.append(dict(type=kind+str(len(ns)),nodes=ns,maps=maps))
    return out


def images(block,marks):return sorted({tuple(sorted(mp.index(i)for i in marks))for mp in block['maps']})


def change(g,w,marks,ext):
    ids=sorted(set(w)|set(marks))
    if inertia(sub(g,ids))!=[len(ids),0,0]:return F(),None
    bs=blocks(g,ids);a,b=[next(x for x in bs if m in x['nodes'])for m in marks]
    co={t:sum(x['type']==t for x in bs)for t in('A1','A2')};has1=ext or co['A1']
    if a is b:
        im=images(a,marks)
        if a['type']=='A16'and set(im)&TARGET:
            assert set(im)<=TARGET
            old=F(5,8)if set(im)&{(1,12),(3,14)}and has1 else F()if set(im)&{(1,12),(3,14)}else F(1)if im==[(2,13)]else F(1,2)
            return-old,old
        if a['type']=='A15'and set(im)&{(0,12),(2,14)}:return F(1,2),F()
        return F(),None
    out=F()
    for sm,la,m,n in((a,b,*marks),(b,a,marks[1],marks[0])):
        st,lt=sm['type'],la['type'];si=images(sm,[m]);li=images(la,[n])
        if st=='A1'and lt=='A13'and(li==[(0,),(12,)]or li==[(2,),(10,)]and co['A2']):out+=F(1,2)
        if st=='A4'and lt=='A12'and si==[(0,),(3,)]and li==[(1,),(10,)]and has1:out+=F(5,8)
        if st=='A2'and lt=='A14'and li==[(2,),(11,)]and has1:out+=F(1,2)
        if st=='A4'and lt=='A13'and si==[(1,),(2,)]and li==[(2,),(10,)]:out+=F(1,2)
        if st=='A1'and lt=='D13'and li==[(11,),(12,)]and co['A2']:out+=F(1,2)
        if st=='A4'and lt=='D13'and si==[(1,),(2,)]and li==[(11,),(12,)]:out+=F(1,2)
    return out,None


def distances(g,starts):
    ds={i:0 for i in starts};todo=list(starts)
    for i in todo:
        for j in range(len(g)):
            if i!=j and g[i][j]and j not in ds:ds[j]=ds[i]+1;todo.append(j)
    return ds


def geometry(local):
    checks=0
    for p in local['local']['selected_and_later_seals']:
        g=[list(map(F,row))for row in p['gram']];w=p['white'];n=len(g);s=p['seal'];cs=p['corners'];sel=p['selected_marks']
        assert set(x for row in g for x in row)<={-1,0,2}
        black=sorted(set(range(n))-set(w));assert len(black)==4
        left=[sel[0]]+[i for i in black if i!=sel[0]and all(set([sel[0],i])!=set(c['marks'])for c in cs)]
        right=sorted(set(black)-set(left));assert len(left)==len(right)==2
        du,dv=distances(g,left),distances(g,right)
        for ids in[w]+[w+[b]for b in black]+[sorted(set(w)|set(c['marks']))for c in cs]:
            assert inertia(sub(g,ids))==[len(ids),0,0];checks+=1
        ids=sorted(set(w)|set(sel));bs=blocks(g,ids);allowed=sorted(set(black)-set(sel))
        for b in bs:
            iv=inverse(sub(g,b['nodes']));assert all(x>=0 for row in iv for x in row)
            allowed.extend(i for k,i in enumerate(b['nodes'])if iv[k][k]<=2)
        allowed=sorted(allowed)
        if s['method']=='full_radical':
            assert inertia(g)==[17,1,1]and allowed==s['allowed_contacts']
            z=list(map(F,s['radical']));assert all(sum(g[i][j]*z[j]for j in range(n))==0 for i in range(n))
            assert all(z[i]for i in allowed)
            expected={(i,j):min(du[i]+dv[j],du[j]+dv[i])+2 for i,j in combinations(allowed,2)if z[i]*z[j]<0}
            assert expected and max(expected.values())<=8
            assert{tuple(r['nodes']):r['length']for r in s['mixed_sign_shortcuts']}==expected
        else:
            assert inertia(g)==[17,1,0];iv=inverse(g);assert iv==[list(map(F,row))for row in s['inverse']]
            rem=sorted(set(allowed)-set(p['selected_nodes']));assert rem==s['remaining_contacts']
            assert all(iv[i][i]>2 for i in rem)
            end=s['end_proof'];assert allowed==end['allowed_contacts'];a,b=end['end_sets']
            assert not set(a)&set(b)and sorted(a+b)==allowed
            assert all(iv[i][j]>0 for side in(a,b)for i in side for j in side)
            bound=max(du[i]for i in a)+max(dv[i]for i in b)+2
            assert bound==end['both_ends_shortcut_bound']<=8
            if p['index']==1010:assert rem==[1,11]and iv[1][1]==iv[11][11]==8
        for var in p['variants']:
            total=F()
            for c in var['corners']:
                d,before=change(g,w,c['marks'],var['external_A1']);assert d==F(c['delta'])
                if before is None:before=F(c['saved_previous'])
                assert before==F(c['current_before'])and before+d==F(c['payment'])>=0;total+=before+d
            assert total==F(var['total'])>=1
    assert len(local['local']['selected_and_later_seals'])==17
    return checks


def domain():
    out=set()
    def rec(r,left,cs):
        if r==18:
            if not cs[7]:return
            used=sum(r*c for r,c in zip(range(9,18),cs))
            for a12 in range(cs[3]+1):
                if not(a12 or cs[4]or cs[5]or cs[6]):continue
                for k5 in range((57-used-18)//5+1):out.add((tuple(cs),k5,a12))
            return
        for c in range(left//r+1):rec(r+1,left-r*c,cs+[c])
    rec(9,39,[]);return out


def upper(doc):
    assert{(tuple(c['counts']),c['k5'],c['a12'])for c in doc['classes']}==domain()
    assert len(doc['classes'])==16;maximum=F();checks=preds=rows=0
    for cl in doc['classes']:
        cs=cl['counts'];k5=cl['k5'];a12=cl['a12'];rem=57-sum(r*c for r,c in zip(range(9,18),cs))-5*k5
        assert rem==cl['remaining_rank'];tabs=cl['tables'];weights=dict(tabs[0]['weights'])
        assert dict(tabs[1]['weights'])=={r:w for r,w in weights.items()if r!=2}
        for t in tabs:
            ws=dict(t['weights']);v=t['values'];last=t['last'];assert v[0]==0 and last[0]is None
            for n in range(1,58):
                possible=[r for r in ws if r<=n and v[n-r]is not None]
                for r in possible:assert v[n]>=v[n-r]+ws[r];checks+=1
                if possible:r=last[n];assert r in possible and v[n]==v[n-r]+ws[r];preds+=1
                else:assert v[n]is None and last[n]is None
        constant=None;mx=F();assert[r['reserved_D_rank']for r in cl['reserved_rows']]==list(range(18,rem+1))
        for row in cl['reserved_rows']:
            s=row['reserved_D_rank'];n=rem-s;assert n==row['residual_rank'];normal,no2=[t['values']for t in tabs];opts=[normal[n]]
            for a1 in(1,2):
                old=a1==2 and bool(cs[1]);new=bool(cs[0]and k5)
                if n<a1 or not(old or new)or(a1==2 and not old):continue
                for a2 in range((n-a1)//2+1):
                    v=no2[n-a1-2*a2]
                    if v is not None:opts.append(a1*weights[1]+a2*weights[2]+v+(36*a2*(a2-1)if old else 0)+(64*a2*k5 if new else 0))
            best=max(v for v in opts if v is not None);delta=16*k5*a12+(8*k5*n if cs[2]else 0)
            assert delta==64*F(row['order17_upper']);fc=64*F(row['ordinary_upper'])-best-delta-weights[s]
            if constant is None:constant=fc
            assert fc==constant;q=None
            if cs[1]and n>=7:
                for a1 in range(1,n-5):
                    for a4 in range(1,(n-a1-2)//4+1):
                        for a2 in range(1,(n-a1-4*a4)//2+1):
                            v=no2[n-a1-4*a4-2*a2]
                            if v is None:continue
                            z=fc+weights[s]+a1*weights[1]+a2*weights[2]+a4*weights[4]+v+delta+(36*a2*(a2-1)if a1>=2 else 0)+(64*a2*k5 if cs[0]else 0)+64*a2*a4
                            q=z if q is None else max(q,z)
            assert(q is None)==(row['quarter_upper']is None)
            if q is not None:assert q==64*F(row['quarter_upper'])
            z=max(64*F(row['ordinary_upper']),q if q is not None else -1)
            assert z==64*F(row['upper']);mx=max(mx,F(z,64));rows+=1
        assert mx==F(cl['upper']);maximum=max(maximum,mx)
    assert maximum==371 and checks==doc['Bellman_inequalities']and preds==doc['attaining_predecessors']
    assert rows==doc['reserved_D_options'];return checks,preds,rows


def main():
    up,local=load(UP),load(LOCAL);hashes=0
    for doc in(up,local):
        for p,h in doc['sources'].items():assert digest(ROOT/p)==h;hashes+=1
    geom=geometry(local);dp=upper(up)
    lit=local['current_literal_interface'];assert lit['configurations']==216 and lit['guarded_variants']==432
    for p in lit['rows']:
        for mode in p['modes']:
            assert sum(F(c['payment'])for c in mode['corners'])==F(mode['total'])>=1
            for c in mode['corners']:assert F(c['current_before'])+F(c['delta'])==F(c['payment'])>=0
    affected={r['case']for r in lit['negative_A16_outputs']};assert affected=={6,7,12,13,18,19,24,25}
    old=load('fresh_affine_d11_leaf_pair_transfer_20261017')['two_A1_A3_spectators']['rows']
    for i in affected:
        p=old[i];g=[list(map(F,row))for row in p['gram']];found=[]
        for c in p['corners']:
            bs=blocks(g,sorted(set(p['white'])|set(c['marks'])))
            b=next((b for b in bs if b['type']=='D13'and set(c['marks'])<=set(b['nodes'])),None)
            if b:assert images(b,c['marks'])==[(1,11),(1,12)];found.append(c['marks'])
        assert len(found)==1
    later=local['current_later_interfaces'];assert later['later_corners']==564 and len(later['target_hits'])==20
    assert all(p['index']in(1167,2015)for p in later['replaced_own_newer_outputs'])
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<1024**3
    print('PASS_INDEPENDENT_A16_DOUBLE_LARGE_D',hashes,geom,*dp,432,len(affected),flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);R.resources();main();R.resources(start=False)
