#!/usr/bin/env python3
"""Independent rational elimination replay of the H18 seam-face witnesses.

Does not import the search or its field/matrix arithmetic. The mixed field
occurs only in the pairing between two projected walls, never in their
pairings with a seam face. Thus rational elimination suffices.
"""
from __future__ import annotations

import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


def pair(value):
    return (F(value['a']),F(value['b']))


def sign(value):
    a,b=value
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if (a>0)==(b>0):return (a>0)-(a<0)
    return ((a>0)-(a<0)) if a*a>2*b*b else ((b>0)-(b<0))


def inverse_positive(matrix):
    n=len(matrix)
    work=[row[:]+[F(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    for k in range(n):
        pivot=work[k][k]
        if pivot<=0:return None
        work[k]=[x/pivot for x in work[k]]
        for i in range(n):
            if i==k:continue
            factor=work[i][k]
            work[i]=[a-factor*b for a,b in zip(work[i],work[k])]
    return [row[n:] for row in work]


def main():
    path=Path('data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json')
    payload=json.loads(path.read_text())
    labels=payload['labels']; gram=[[pair(x) for x in row] for row in payload['gram']]
    seam=sorted(labels.index(f'L{int(x)+1}') for x in payload['construction']['matching'])
    allowed={(F(0),F(0)),(F(1,4),F(0)),(F(1,2),F(0)),(F(3,4),F(0)),(F(1,2),F(1,4))}
    seam_witnesses=[]
    for index in seam:
        facet=next(x for x in payload['facets'] if x['facet']==labels[index])
        assert not facet['coxeter_angle_gate']
        witness=facet['first_forbidden_pair']
        seam_witnesses.append({'face':[labels[index]],'forbidden_pair':witness['walls'], 'pairing':witness['pairing'], 'normalized_pairing_square':witness['normalized_pairing_square']})
    facet_row={'codimension':1,'dimension':17,'elliptic_subsets':len(seam),'exact_rejections':len(seam),'survivors':[],'witnesses':seam_witnesses}
    checked=[]
    for row in [facet_row,*payload['seam_faces']]:
        witnesses={tuple(labels.index(x) for x in item['face']):item for item in row['witnesses']}
        eligible=set()
        for selection in combinations(seam,row['codimension']):
            assert all(gram[i][j][1]==0 for i in selection for j in range(len(gram)))
            inverse=inverse_positive([[gram[i][j][0] for j in selection] for i in selection])
            if inverse is None:continue
            eligible.add(selection)
            item=witnesses[selection]
            i,j=[labels.index(x) for x in item['forbidden_pair']]
            assert i not in selection and j not in selection and i!=j
            def schur(k,l):
                subtract=sum(gram[k][s][0]*inverse[a][b]*gram[t][l][0] for a,s in enumerate(selection) for b,t in enumerate(selection))
                return gram[k][l][0]-subtract,gram[k][l][1]
            ni=schur(i,i);nj=schur(j,j);p=schur(i,j)
            assert ni[1]==nj[1]==0 and ni[0]>0 and nj[0]>0
            product=ni[0]*nj[0]
            q=((p[0]*p[0]+2*p[1]*p[1])/product,2*p[0]*p[1]/product)
            assert p==pair(item['pairing']) and q==pair(item['normalized_pairing_square'])
            assert sign(p)>0 or (sign((q[0]-1,q[1]))<0 and q not in allowed)
        assert eligible==set(witnesses)
        assert len(eligible)==row['elliptic_subsets']==row['exact_rejections']
        assert not row['survivors']
        checked.append({'dimension':row['dimension'],'exact_witnesses_verified':len(eligible)})
    print(json.dumps({'status':'PASS_INDEPENDENT_RATIONAL_REPLAY','checks':checked},indent=2))


if __name__=='__main__':main()
