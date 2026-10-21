#!/usr/bin/env python3
"""Small standard-library independent certificate check; no producer imports."""
import json
import signal
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
UP='fresh_d17_node3_large_d_guard_upper_20261018'
LOCAL='fresh_d17_node3_large_d_guard_transfer_20261018'


def digest(path):
    h=sha256()
    with path.open('rb')as stream:
        for chunk in iter(lambda:stream.read(1048576),b''):h.update(chunk)
    return h.hexdigest()


def inertia(matrix):
    a=[list(map(F,row))for row in matrix];pos=neg=zero=0
    while a:
        n=len(a);i=next((i for i in range(n)if a[i][i]),None)
        if i is not None:
            ids=[j for j in range(n)if j!=i];p=a[i][i]
            pos+=int(p>0);neg+=int(p<0)
            a=[[a[j][k]-a[j][i]*a[i][k]/p for k in ids]for j in ids]
        else:
            pair=next(((i,j)for i in range(n)for j in range(i)if a[i][j]),None)
            if pair is None:
                zero+=n;break
            i,j=pair;b=a[i][j];ids=[k for k in range(n)if k not in pair]
            a=[[a[k][l]-(a[k][i]*a[j][l]+a[k][j]*a[i][l])/b for l in ids]for k in ids]
            pos+=1;neg+=1
    return[pos,neg,zero]


def domain():
    expected=set()
    def rec(rank,left,counts):
        if rank==18:
            if not(counts[4]or counts[5]):return
            used=sum(r*c for r,c in zip(range(9,18),counts))
            for k5 in range((57-used-18)//5+1):
                for a12 in range(counts[3]+1):expected.add((tuple(counts),k5,a12))
            return
        for c in range(left//rank+1):rec(rank+1,left-rank*c,counts+[c])
    rec(9,39,[])
    return expected


def main():
    upper=json.loads((ROOT/f'results/{UP}.json').read_text())
    local=json.loads((ROOT/f'results/{LOCAL}.json').read_text())
    sourcechecks=0
    for doc in(upper,local):
        for p,h in doc['sources'].items():assert digest(ROOT/p)==h;sourcechecks+=1
    observed={(tuple(r['counts']),r['k5'],r['a12'])for r in upper['classes']}
    assert observed==domain()and len(observed)==135==len(upper['classes'])
    checks=preds=rows=0;maximum=F()
    for cl in upper['classes']:
        cs=cl['counts'];k5=cl['k5'];a12=cl['a12'];rem=57-sum(r*c for r,c in zip(range(9,18),cs))-5*k5
        assert rem==cl['remaining_rank']
        tabs=cl['tables'];weights=dict(tabs[0]['weights'])
        assert dict(tabs[1]['weights'])=={r:w for r,w in weights.items()if r!=2}
        for tab in tabs:
            ws=dict(tab['weights']);v=tab['values'];last=tab['last']
            assert v[0]==0 and last[0]is None and len(v)==len(last)==58
            for n in range(1,58):
                possible=[r for r in ws if r<=n and v[n-r]is not None]
                for r in possible:assert v[n]>=v[n-r]+ws[r];checks+=1
                if possible:
                    r=last[n];assert r in possible and v[n]==v[n-r]+ws[r];preds+=1
                else:assert v[n]is None and last[n]is None
        assert [r['reserved_D_rank']for r in cl['reserved_rows']]==list(range(18,rem+1))
        constant=None;clmax=F()
        for row in cl['reserved_rows']:
            s=row['reserved_D_rank'];n=rem-s;assert n==row['residual_rank'];normal,no2=[t['values']for t in tabs]
            opts=[normal[n]]
            for a1 in(1,2):
                old=a1==2 and bool(cs[1]);new=bool(cs[0]and k5)
                if n<a1 or not(old or new)or(a1==2 and not old):continue
                for a2 in range((n-a1)//2+1):
                    v=no2[n-a1-2*a2]
                    if v is not None:opts.append(a1*weights[1]+a2*weights[2]+v+(36*a2*(a2-1)if old else 0)+(64*a2*k5 if new else 0))
            best=max(v for v in opts if v is not None)
            delta=16*k5*a12+(8*k5*n if cs[2]else 0)
            assert delta==64*F(row['order17_delta_upper'])
            # The reserved large D item has no interaction and no new coupon.
            fc=64*F(row['ordinary_upper'])-best-delta-weights[s]
            if constant is None:constant=fc
            assert fc==constant
            q=None
            if cs[1]and n>=7:
                for a1 in range(1,n-5):
                    for a4 in range(1,(n-a1-2)//4+1):
                        for a2 in range(1,(n-a1-4*a4)//2+1):
                            rest=n-a1-4*a4-2*a2;v=no2[rest]
                            if v is None:continue
                            z=fc+weights[s]+a1*weights[1]+a2*weights[2]+a4*weights[4]+v+delta+(36*a2*(a2-1)if a1>=2 else 0)+(64*a2*k5 if cs[0]else 0)+64*a2*a4
                            q=z if q is None else max(q,z)
            assert(q is None)==(row['quarter']is None)
            if q is not None:assert q==row['quarter'][0]
            z=max(64*F(row['ordinary_upper']),q if q is not None else -1)
            assert z==64*F(row['upper']);clmax=max(clmax,F(z,64));rows+=1
        assert clmax==F(cl['upper']);maximum=max(maximum,clmax)
    assert maximum==365 and checks==upper['Bellman_inequalities']and preds==upper['attaining_predecessors']
    corechecks=0
    for proof in local['local_proof']['proofs']:
        g=[list(map(F,row))for row in proof['gram']];w=proof['white'];z=list(map(F,proof['radical']))
        assert inertia(g)==[17,1,1]
        assert all(sum(g[i][j]*z[j]for j in range(19))==0 for i in range(19))
        black=sorted(set(range(19))-set(w))
        for ids in[w]+[w+[b]for b in black]+[sorted(set(w)|set(c['marks']))for c in proof['corners']]:
            assert inertia([[g[i][j]for j in ids]for i in ids])==[len(ids),0,0];corechecks+=1
        assert sum(F(c['payment'])for c in proof['corners'])==1
    assert len(local['current_literal_interface']['negative_outputs'])==8
    assert local['finite']['new_exact_global_formal_maximum']=='372'
    print('PASS_INDEPENDENT_D17_GUARD_CERTIFICATES',sourcechecks,checks,preds,rows,corechecks,maximum)


if __name__=='__main__':
    signal.alarm(59);main()
