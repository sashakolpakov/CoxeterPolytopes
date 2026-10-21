#!/usr/bin/env python3
"""Certify a strict escape for the transported height141 chamber prefix."""
import argparse,json,sys
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import search_fresh_cubic7_h9_20261013 as V
import fresh_h9_cubic_corner_20261014 as P
TAG='fresh_h9_cubic_corner_20261014_escape'
SOURCE='results/fresh_h9_cubic_corner_20261014_h141_orbit.json'
def dec(x):return V.C(x['a'],x['b'],x['c'])
def run(saved=None):
    d=json.loads((ROOT/SOURCE).read_text());N=[tuple(map(dec,r))for r in d['normals']]
    a=dec(d['timelike_coefficient']);G=V.gram_for(a);H=tuple(map(dec,d['strict_center']))
    if saved is None:
        proposal=V.exact_escape(N,G)
        assert proposal['verdict']=='EXACT_SPACELIKE_ESCAPE',proposal
        X=tuple(map(dec,proposal['vector']));epsilon=Fraction(1,2)
        while True:
            Y=tuple(x+epsilon*h for x,h in zip(X,H))
            if V.pair(Y,G,Y).sign()>0:break
            epsilon/=2
    else:Y=tuple(map(dec,saved['strict_escape']))
    ips=[V.pair(Y,G,r)for r in N];q=V.pair(Y,G,Y)
    assert q.sign()>0 and Y[0].sign()>0 and all(x.sign()<0 for x in ips)
    assert V.pair(H,G,H).sign()<0 and all(V.pair(H,G,r).sign()<0 for r in N)
    # Repeat only the decisive geometry using the separate Fraction-polynomial engine.
    K=P.K;ak=K.decode(d['timelike_coefficient']);yk=[K.decode(x.to_json())for x in Y]
    nk=[[K.decode(x)for x in r]for r in d['normals']]
    assert K.sign(K.pairing(yk,yk,ak))>0
    assert all(K.sign(K.pairing(yk,r,ak))<0 for r in nk)
    paths=[SOURCE,f'scripts/{TAG}.py','scripts/search_fresh_cubic7_h9_20261013.py',
           'scripts/fresh_h9_cubic_corner_20261014.py','scripts/audit_fresh_cubic7_h9_no_import_20261013.py',
           'src/coxeter_polytopes/exact/cubic7.py']
    return dict(status='EXACT_STRICT_ESCAPE_OF_32_WALL_PREFIX',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        strict_escape=V.encode_vector(Y),norm=q.to_json(),wall_pairings=V.encode_vector(ips),
        scope='The 32-wall prefix has infinite volume. This does not establish nonreflectivity of the full arithmetic lattice; undiscovered mirrors may cut this escape. Decisive signs are repeated in a separate field implementation.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();p=ROOT/f'results/{TAG}.json'
    saved=json.loads(p.read_text())if args.verify else None;out=run(saved)
    if args.verify:assert out==saved
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
