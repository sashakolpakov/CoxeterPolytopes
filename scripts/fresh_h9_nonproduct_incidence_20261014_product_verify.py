#!/usr/bin/env python3
"""Exact classification of all sixteen archived 216-vertex proposals.

The geometric exclusion invokes Alexandrov, arXiv2203.07248v2, Theorems
A/B. This script proves the literal abstract combinatorial identification.
"""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
import argparse, hashlib, json
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_h9_incidence_repair_20261014_released_sat.json'
OUT=ROOT/'results/fresh_h9_nonproduct_incidence_20261014_product_verify.json'

def verify():
    d=json.loads(SOURCE.read_text());out=[]
    for ci,case in enumerate(d['cases']):
        for mi,model in enumerate(case['models']):
            vertices=set(map(tuple,model['vertices']))
            assert len(vertices)==216 and all(len(V)==len(set(V))==9 for V in vertices)
            pairs={B for V in vertices for B in combinations(V,2)}
            triples={B for V in vertices for B in combinations(V,3)}
            M=[B for B in combinations(range(15),2)if B not in pairs]
            M +=[B for B in combinations(range(15),3)if B not in triples and all(C in pairs for C in combinations(B,2))]
            assert Counter(map(len,M))=={2:3,3:3}
            assert set().union(*map(set,M))==set(range(15))
            assert all(not set(A)&set(B)for A,B in combinations(M,2))
            expected={tuple(sorted(set(range(15))-set(omitted)))for omitted in product(*M)}
            assert expected==vertices
            # Equality of all maximal faces proves equality of the entire
            # downward-closed complexes, so there are no larger missing faces.
            out.append({'source_case':ci,'source_model':mi,'minimal_nonfaces':[list(B)for B in M],
                        'factor_dimensions':[len(B)-1 for B in M],
                        'complete_maximal_face_equality':True,'vertices':len(vertices)})
    assert len(out)==16
    return {'status':'EXACT_ALL_SIXTEEN_COMBINATORIAL_TYPES_IDENTIFIED',
            'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[SOURCE,Path(__file__).resolve()]},
            'cases':out,'type':'Delta2 x Delta2 x Delta2 x I x I x I',
            'external_theorem':{'author':'Stepan Alexandrov','arxiv':'2203.07248v2','url':'https://arxiv.org/html/2203.07248v2','statements':'Theorems A and B'},
            'conclusion':'No compact hyperbolic Coxeter realization of this combinatorial type, for any finite labels or dotted distances. Six simplex factors include noninterval factors; this is excluded by Theorem B, and the corresponding six disjoint minimal nonfaces also meet Theorem A.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=verify()
    if a.verify:assert d==json.loads(OUT.read_text());print('PASS exact product-type replay')
    else:OUT.write_text(json.dumps(d,indent=2)+'\n');print('WROTE',OUT)
if __name__=='__main__':main()
