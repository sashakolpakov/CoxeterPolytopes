#!/usr/bin/env python3
"""Independent stdlib scalar reconstruction; no producer arithmetic imports.

Optional low/high factors are chosen first, followed by multiplicity vectors
for ranks10..17. Sliding windows compute H; CDF tails compute G. The primary
may still change its separate global-constant proof: only the immutable
scalar rows and stated extrema are compared here.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from math import prod
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_ambient65_ideal_margin_closure_20261017'
TAG=PRIMARY+'_review'


def middle_vectors(total,rank=10,prefix=()):
    if rank==18:
        if total==0:yield prefix
        return
    for count in range(total//rank+1):
        yield from middle_vectors(total-count*rank,rank+1,prefix+(count,))


def terminal_profiles():
    out=[]
    for low in range(10):
        for high in [0,*range(18,65)]:
            rem=64-low-high
            if rem<0:continue
            for counts in middle_vectors(rem):
                r=([low]if low else[])+[k for k,c in zip(range(10,18),counts)for _ in range(c)]+([high]if high else[])
                assert sum(r)==64 and r==sorted(r)
                out.append(tuple(r))
    assert len(out)==len(set(out))==1650
    return sorted(out)


def costs(r):
    if r<=10:return F(r,2)
    if r>=17:return F(163,23)
    return {11:F(16,3),12:F(74,13),13:F(6),14:F(94,15),15:F(107,16),16:F(7)}[r]


def row(rs):
    h=[1]
    for r in rs:
        new=[];window=0
        for k in range(len(h)+r):
            if k<len(h):window+=h[k]
            if k-r-1>=0:window-=h[k-r-1]
            new.append(window)
        h=new
    V=prod(r+1 for r in rs);assert len(h)==65 and sum(h)==V and h==h[::-1]
    prefix=0;G=0
    for k in range(32):
        prefix+=h[k];G+=V-2*prefix
    assert 2*G==64*V-sum(n*abs(64-2*k)for k,n in enumerate(h))
    charge=sum(map(costs,rs),F());curv=F(5*G-64*V,2);margin=curv/V-charge
    assert margin>0
    return dict(ranks=list(rs),h_polynomial=h,link_vertices=V,G=G,normalized_cost=str(charge),curvature=str(curv),normalized_margin=str(margin),absolute_margin=str(curv-V*charge))


def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()


def run():
    primary=json.loads((ROOT/'results'/f'{PRIMARY}.json').read_text())
    own=[row(rs)for rs in terminal_profiles()]
    saved={tuple(x['ranks']):x for x in primary['rows']}
    assert len(saved)==len(primary['rows'])==1650
    assert set(saved)=={tuple(r['ranks'])for r in own}
    for r in own:assert r==saved[tuple(r['ranks'])],r['ranks']
    nmin=min(F(r['normalized_margin'])for r in own);amin=min(F(r['absolute_margin'])for r in own)
    nm=[r for r in own if F(r['normalized_margin'])==nmin];am=[r for r in own if F(r['absolute_margin'])==amin]
    assert nmin==F(61,1936)and amin==F(445,23)
    assert len(nm)==len(am)==1 and nm[0]['ranks']==[9,10,10,10,10,15]and am[0]['ranks']==[64]
    assert primary['minimum_normalized_margin']==nm[0]and primary['minimum_absolute_margin']==am[0]
    assert primary['ambient_dimension']==65 and primary['complete_terminal_profiles']==1650
    return dict(status='PASS_INDEPENDENT_AMBIENT65_ALL1650_SCALAR_ROWS',sources={f'scripts/{TAG}.py':sha256(Path(__file__).read_bytes()).hexdigest()},
        producer_scalar_rows_sha256=sha256(canon(sorted(primary['rows'],key=lambda r:r['ranks']))).hexdigest(),
        independent_rows_sha256=sha256(canon(own)).hexdigest(),profiles=1650,rank_sum=64,
        enumeration='Choose optional low1..9 and optional high18..64 first, then all nonnegative multiplicities at10..17 summing to the remainder.',
        calculation='Integer sliding-window convolution for all65 H coefficients; symmetric CDF-tail sum for G, checked by absolute-deviation identity; exact Fraction costs and margins.',
        minimum_normalized=nm[0],minimum_absolute=am[0],
        all_rows=[{k:v for k,v in r.items()if k!='h_polynomial'}for r in own],
        scope='Every primary scalar row and both unique extrema independently reconstructed, importing only the standard library. No producer arithmetic, geometry census, all-rank concentration theorem, finite-coupon budget or additive global-constant proof is rerun. The primary file/source hash is deliberately not frozen by this scalar-only certificate while its separate global-constant proof is being finalized.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    assert out['producer_scalar_rows_sha256']==out['independent_rows_sha256']
    path=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['profiles'],out['minimum_normalized']['normalized_margin'],out['minimum_absolute']['absolute_margin'],flush=True)
