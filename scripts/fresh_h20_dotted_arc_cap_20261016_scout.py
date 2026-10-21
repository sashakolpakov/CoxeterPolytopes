#!/usr/bin/env python3
"""Continuous one-weight / adjacent-chain right-corner numerical screen.

The saved bounds are proposals, not exact exclusion certificates. No labels or
coefficient fields are bounded: the scanned family has one nonzero support
contact and is right to two adjacent chain walls.
"""
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import fresh_h20_nonright_s61_20261014_scout as C

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_dotted_arc_cap_20261016_scout'

def run():
    d=C.load();base=d['base'];J=d['J'];S=d['S'];P=d['P'];W=d['W']
    bank=json.loads((ROOT/'results/fresh_h20_nonright_s61_20261014_bank_gate.json').read_text())
    E=mp.matrix([d['f'](x) for x in bank['escape']['rational_point']])
    ep=mp.matrix([d['dot'](v,E) for v in d['frame']])
    # Reconstruct support dual pairing from the rational physical escape.
    ae=mp.matrix([-sum(J[j,k]*d['dot'](d['N'][S[k]],E) for k in range(18)) for j in range(18)])
    def mid(v):return (d['f'](v[0])+d['f'](v[1]))/2
    def out(x):return mp.nstr(x,45)
    rows=[];events=[]
    for vertex in base['face']['vertices']:
        pair=vertex['walls']
        if not (min(pair)>=107 and pair[1]-pair[0]==1):continue
        v=mp.matrix([mid(x)for x in vertex['frame_point_bounds']])
        q=-v[0]**2+v[1]**2+v[2]**2
        assert q<0
        v/=mp.sqrt(-q)
        vp=P*v;ve=(ep.T*v)[0]
        for j in range(18):
            low=1/mp.sqrt(J[j,j]);cut=-ve/ae[j]
            upper=[]
            for i in range(133):
                if i in S or i in pair:continue
                assert W[j,i]<mp.mpf('1e-60')
                if W[j,i]<-mp.mpf('1e-60'):
                    upper.append((vp[i]/W[j,i],i))
            hi,block=min(upper)
            lb=max(low,cut)
            row=dict(pair=pair,wall=S[j],dual_index=j,dual_norm=out(J[j,j]),
                     norm_lower=out(low),escape_lower=out(cut),acute_upper=out(hi),
                     blocker=block,gap=out(hi-lb),acute_possible=bool(hi>low),
                     escape_possible=bool(hi>lb))
            if hi>lb:
                lam=(lb+hi)/2;t=mp.sqrt(2/(J[j,j]*lam**2-1));c=lam*t
                g=[t*(vp[i]-lam*W[j,i])for i in range(133)]
                row.update(sample_lambda=out(lam),sample_t=out(t),sample_contact=out(c),
                           sample_escape=out(t*(ve+lam*ae[j])),sample_pairs=list(map(out,g)),
                           forbidden=[i for i,x in enumerate(g)if C.classify(float(x))in['obtuse','forbidden']])
            rows.append(row)
            outer=pair[0]-1 if pair[0]>107 else pair[1]+1
            ell=-vp[outer]
            for m in [3,4,5,6,7,8,9,10,12,16,32,100,0]:
                k=2 if m==0 else 2*mp.cos(mp.pi/m)
                t=k/ell;c=mp.sqrt((2+t*t)/J[j,j]);lam=c/t
                cutpair=t*ve+c*ae[j]
                if cutpair<=0 or lam>hi:continue
                gs=[t*vp[i]-c*W[j,i]for i in range(133)]
                bad=[dict(wall=i,pair=out(g),kind=C.classify(float(g)))for i,g in enumerate(gs)
                     if C.classify(float(g))in['obtuse','forbidden']]
                events.append(dict(pair=pair,wall=S[j],outer_label=m or 'parallel',
                                   t=out(t),contact=out(c),escape=out(cutpair),bad=bad,
                                   all_pairs=list(map(out,gs))))
    assert len(rows)==450
    paths=['scripts/'+TAG+'.py','scripts/fresh_h20_nonright_s61_20261014_scout.py',
           'results/fresh_h20_s61_multicap_20261014_exact.json',
           'results/fresh_h20_nonright_s61_20261014_bank_gate.json']
    return dict(status='NUMERICAL_CONTINUOUS_ADJACENT_RIGHT_CHAIN_SINGLE_SUPPORT_SCREEN',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                cases=len(rows),acute_possible=sum(r['acute_possible']for r in rows),
                escape_possible=sum(r['escape_possible']for r in rows),rows=rows,events=events,
                inputs=dict(P=[[out(x)for x in row]for row in P.tolist()],
                            W=[[out(x)for x in row]for row in W.tolist()],
                            J=[[out(x)for x in row]for row in J.tolist()],
                            ep=list(map(out,ep)),AE=list(map(out,ae)),support=S))

if __name__=='__main__':
    z=run();(ROOT/'results'/f'{TAG}.json').write_text(json.dumps(z,indent=2)+'\n')
    print(z['status'],z['cases'],z['acute_possible'],z['escape_possible'],flush=True)
