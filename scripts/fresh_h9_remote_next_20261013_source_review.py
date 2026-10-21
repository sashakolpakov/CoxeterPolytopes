#!/usr/bin/env python3
"""Narrow formula review and higher-precision replay, not a new pair engine."""
import argparse,json,hashlib
from pathlib import Path
from scripts import fresh_h9_remote_next_20261013_verify as N
V=N.V;K=V.K
OUT=V.ROOT/'results/fresh_h9_remote_next_20261013_source_review.json'

def run(verify=False):
    # Twice the constructor verifier's 80-decimal scale, with cached radicals
    # discarded before any reconstruction. Arithmetic implementation is reused.
    V.SCALE=10**160
    V.PHI=V.div(V.add(V.rat(1),V.sqrt(V.rat(5))),V.rat(2));V.radical.cache_clear()
    source=N.PREFIX+'_caps.json';data=V.read(source);replay=N.check(data)
    # Independent symbolic frame calculation in variables a,b,c,u,L:
    # b=-(1+L)/2,c=(L-1)/2,2u²=2+a²-b²-c².
    import sympy as sp
    aa,L=sp.symbols('a L');b=-(1+L)/2;c=(L-1)/2;u2=(2+aa*aa-b*b-c*c)/2
    assert sp.expand(b+c)==-1 and sp.expand(b-c)==-L
    assert sp.expand(-aa*aa+b*b+c*c+2*u2)==2
    assert sp.expand(b*b+c*c)==(1+L*L)/2
    # H4 chain [5,3,3] has leading minors 2,3-phi,4-2phi,5-3phi.
    M=[[K.scale(K.O,2)if i==j else K.neg(K.P)if{i,j}=={0,1}else K.neg(K.O)if abs(i-j)==1 else K.Z for j in range(4)]for i in range(4)]
    pivots=[]
    while M:
        p=M[0][0];assert K.sign(p)>0;pivots.append(K.enc(p))
        M=[[K.sub(M[i][j],K.div(K.mul(M[i][0],M[0][j]),p))for j in range(1,len(M))]for i in range(1,len(M))]
    result={'status':'PASS_NARROW_SOURCE_REVIEW_AND_160_DIGIT_INTERVAL_REPLAY','source':source,'source_sha256':V.sha(source),
      'verifier_source_sha256':V.sha('scripts/fresh_h9_remote_next_20261013_verify.py'),
      'facets':replay['facets'],'new_pair_counts':replay['new_pair_counts'],
      'actual_vertex_type':replay['actual_vertex']['type'],'nonincident_vertex_checks':len(replay['actual_vertex']['nonincident_checks']),
      'norm_identity':True,'pi3_and_pi10_pairing_identities':True,'H4_positive_pivots':pivots,
      'remaining_escape_norm':replay['remaining_escape']['norm'],'strict_escape_inequalities':len(replay['remaining_escape']['strict_pairings']),
      'scope':'Narrow source-level proof review, separate symbolic frame identities, and full constructor-verifier replay at twice its decimal interval scale. This reuses the existing interval arithmetic and source reconstruction; it is not a separate implementation of all86 root or pair formulas. No finite-volume claim.'}
    if verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],result['facets'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
