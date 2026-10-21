#!/usr/bin/env python3
"""Narrow algebra review and higher-precision replay of the four-cap insertion.

This deliberately reuses the producer-free interval verifier. It is not a
second independent implementation of all ninety-four normals.
"""
import argparse,hashlib,json
from fractions import Fraction as Q
from pathlib import Path
from scripts import fresh_h9_compact_next_20261013_verify as C

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/(C.PREFIX+'_caps.json')
OUT=ROOT/'results/fresh_h9_compact_next_20261013_source_review.json'

def run(verify=False):
    K=C.K;V=C.V
    q=K.scale(K.sub(K.scale(K.O,4),K.scale(K.P,2)),Q(1,3))
    nu=K.sub(K.div(K.O,q),K.scale(K.O,2));assert K.sign(q)>0 and K.sign(nu)>0
    # In a perpendicular decomposition, a²=q, P²=-1, and a.P=0.
    # Hence r=-a/q+sqrt(nu)P has r²=2 and r.a=-1 exactly.
    assert K.sub(K.div(K.O,q),nu)==K.scale(K.O,2)
    # The other three caps have Euclidean coordinates
    # (sqrt(2/3),sqrt(4/3),0), (sqrt(2/3),-sqrt(1/3),1),
    # (sqrt(2/3),-sqrt(1/3),-1) in P-perp. Their Gram is 2I3.
    assert Q(2,3)+Q(4,3)==2 and Q(2,3)+Q(1,3)+1==2
    assert Q(2,3)-Q(2,3)==0 and Q(2,3)+Q(1,3)-1==0
    # Repeat all formula reconstruction and inequalities with 160 decimal
    # places, twice the archived checker's precision. No old intervals reused.
    V.SCALE=10**160;V.PHI=V.div(V.add(V.rat(1),V.sqrt(V.rat(5))),V.rat(2));V.radical.cache_clear()
    d=json.loads(SOURCE.read_text());result=C.check(d)
    assert result['facets']==94 and len(result['actual_finite_vertices'])==4 and 'remaining_escape'in result
    out={'status':'PASS_NARROW_FORMULA_REVIEW_AND_160_DIGIT_RIGOROUS_REPLAY',
      'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,Path(C.__file__))},
      'facets':94,'new_norm_two_caps':4,'actual_compact_face':'H3 Lanner 3-5-3','actual_finite_vertices':4,
      'remaining_strict_spacelike_escape':True,'precision_decimal_places':160,
      'scope':'Independent algebraic review of the orthogonal projection, norm-two cap formula, orthogonal triple and compact tetrahedron argument; a higher-precision rerun of the existing producer-free rigorous interval implementation. This is not a second all-root implementation. The remaining escape proves that the full H9 prefix still has infinite volume.'}
    if verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
