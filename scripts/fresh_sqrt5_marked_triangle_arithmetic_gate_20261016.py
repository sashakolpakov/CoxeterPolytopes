#!/usr/bin/env python3
"""Exact arithmetic gate for the fixed196 selected triangle families only."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_marked_a1_facet_mixed_cross_20261016 as E
import fresh_sqrt5_quad_filters_20261015 as G
ROOT=E.ROOT;TAG='fresh_sqrt5_marked_triangle_arithmetic_gate_20261016'

def literal_triangles():
    tri=E.H.read('fresh_sqrt5_triangle_weights_20261015');literal=[]
    for i,row in enumerate(tri['markings']):
        literal.append((f'core{i}',row['order'],row['white'],row['black'],[(a,b,m)for(a,b),m in zip(combinations(range(row['order']),2),row['labels'])]))
    for i in(46,47,51,55):
        row=tri['markings'][i]
        for tail in range(1,18):
            es=[(a,b,m)for(a,b),m in zip(combinations(range(5),2),row['labels'])]+[(0,5,3)]+[(u,u+1,3)for u in range(5,4+tail)]
            literal.append((f'core{i}_tail{tail}',5+tail,row['white']+list(range(5,5+tail)),row['black'],es))
    for tail in range(1,20):
        es=[(0,2,3),(1,2,10),(0,3,3)]+[(u,u+1,3)for u in range(3,2+tail)]
        literal.append((f'I210_tail{tail}',3+tail,list(range(3,3+tail)),[0,1,2],es))
    for row in E.pending_triangle_attachment_scout()['records']:
        name=row['attached_component'];root=row['attachment_root']
        if name.startswith('A')and len(E.MODELS[name][root])<=1:continue
        literal.append((f'core51_attach_{name}_root{root}',row['nodes'],row['white'],row['black'],row['edges']))
    assert len(literal)==196
    return literal

def analyze(label,n,white,black,edges):
    contacts={tuple(sorted((u,v))):m for u,v,m in edges}
    row=dict(nodes=n,white=white,black=black,labels=[contacts.get(pair,2)for pair in combinations(range(n),2)])
    coloring=G.coloring(row)
    out=dict(family=label,literal_row=row,colors=coloring)
    if not coloring['passes']:
        out['status']='EXCLUDE_CYCLIC_FIELD';return out
    matrix=G.scaled_matrix(row,coloring['colors']);conjugate=[[x.conj()for x in r]for r in matrix]
    # Check all three actual finite corners at both embeddings independently
    # by Schur pivots, then the one-dimensional full-star Schur remainder.
    corners=[]
    for omitted in black:
        ids=white+[b for b in black if b!=omitted];record=dict(omitted_black=omitted,corner_ids=ids)
        for name,mat in [('physical',matrix),('conjugate',conjugate)]:
            residual,pivots=G.schur(mat,ids);assert len(residual)==1
            record[name]=dict(corner_pivots=[x.enc()for x in pivots],schur=residual[0][0].enc(),sign=residual[0][0].sign())
        corners.append(record)
    assert len({r['physical']['sign']for r in corners})==1
    assert len({r['conjugate']['sign']for r in corners})==1
    physical=corners[0]['physical']['sign'];other=corners[0]['conjugate']['sign']
    status='PASS_LOCAL_V_ARITHMETIC_GRAM_NECESSARY_TESTS'if physical==-1 and other>=0 else'EXCLUDE_PHYSICAL_NON_LORENTZ'if physical!=-1 else'EXCLUDE_CONJUGATE_NEGATIVE_DIRECTION'
    out.update(status=status,corners=corners,physical_inertia=[n-1,int(physical<0),int(physical==0)]if physical<=0 else[n,0,0],
        conjugate_inertia=[n-1,int(other<0),int(other==0)]if other<=0 else[n,0,0])
    if status.startswith('EXCLUDE'):
        out['scaled_gram']=G.encode(matrix)
        witnesses=[]
        for removed in range(n):
            ids=[u for u in range(n)if u!=removed];sub=G.sub(conjugate,ids);inertia=G.inertia(sub)
            if inertia==[n-2,0,1]:
                vectors=G.nullspace(sub);assert len(vectors)==1
                null=vectors[0];coupling=G.dot([conjugate[removed][u]for u in ids],null)
                if coupling:
                    witnesses.append(dict(removed_node=removed,principal_ids=ids,inertia=inertia,
                        null_vector=[v.enc()for v in null],nonzero_coupling=coupling.enc()))
        out['affine_principal_coupling_witnesses']=witnesses
    return out

def run():
    rows=[]
    for i,row in enumerate(literal_triangles()):
        rows.append(analyze(*row))
        if(i+1)%20==0:print('TRIANGLES',i+1,dict(Counter(r['status']for r in rows)),flush=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{E.TAG}.py',f'results/{E.TAG}.json',f'scripts/{G.TAG}.py',
        'results/fresh_sqrt5_quad_filters_20261015_review.json','results/fresh_sqrt5_triangle_weights_20261015.json']
    return dict(status='PASS_EXACT_SELECTED_TRIANGLE_ARITHMETIC_GATE',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        triangle_families=len(rows),counts=dict(Counter(r['status']for r in rows)),records=rows,
        retained_families=[r['family']for r in rows if r['status'].startswith('PASS_')],
        excluded_families=[r['family']for r in rows if r['status'].startswith('EXCLUDE_')],
        arithmetic='Exact inherited Q(phi) squareclass-rescaled Gram engine; no numerical eigenvalue tests. Every spherical corner is pivoted at both embeddings; the full-star residual is one-dimensional.',
        scope='An exact gate for precisely the196 selected literal triangle families. Retained rows satisfy local necessary Gram tests, not a complete chamber/incidence theorem. Excluded rows cannot be arithmetic with these literal contacts, including any orthogonal spectator extension. This does not enumerate all white extensions or complete the discounted short-face coverage.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(p.read_text())
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
