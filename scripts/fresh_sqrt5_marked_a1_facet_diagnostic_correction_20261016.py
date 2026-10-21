#!/usr/bin/env python3
"""Correct two prose identifications without altering the frozen certificate."""
import argparse,json
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_marked_a1_facet_diagnostic_correction_20261016'
PRIOR='fresh_sqrt5_marked_a1_facet_mixed_cross_20261016'

def graph(vertices,edges):
    g={u:{}for u in vertices}
    for u,v,m in edges:
        if u in g and v in g and m!=2:g[u][v]=m;g[v][u]=m
    return g

def components(g):
    unseen=set(g);out=[]
    while unseen:
        todo=[min(unseen)];seen=set()
        while todo:
            u=todo.pop()
            if u in seen:continue
            seen.add(u);todo.extend(v for v in g[u]if v not in seen)
        unseen-=seen;out.append({u:{v:m for v,m in g[u].items()if v in seen}for u in seen})
    return out

def arms(g):
    centers=[u for u in g if len(g[u])==3];assert len(centers)==1
    center=centers[0];lengths=[]
    for first in g[center]:
        previous=center;u=first;n=1
        while len(g[u])!=1:
            nxt=[v for v in g[u]if v!=previous];assert len(nxt)==1
            previous,u=u,nxt[0];n+=1
        lengths.append(n)
    return center,sorted(lengths)

def classify(g):
    n=len(g);labels=[m for u in g for v,m in g[u].items()if u<v]
    assert len(labels)==n-1
    if n==1:return'A1'
    if max(map(len,g.values()))==2 and 5 in labels:
        assert n in(3,4)and sorted(labels)==[3]*(n-2)+[5]
        return'H'+str(n)
    assert set(labels)=={3}
    _,a=arms(g)
    if a[:2]==[1,1]:return'D'+str(n)
    assert a in([1,2,2],[1,2,4]);return'E'+str(n)

def run():
    p=ROOT/f'results/{PRIOR}.json';doc=json.loads(p.read_text());out=[]
    for row in doc['single_attachment_triangle_domain']['records']:
        if row['attached_component']!='D5':continue
        root=row['attachment_root'];attached=graph(list(range(5,row['nodes'])),row['edges'])
        center,_=arms(attached);front={root+5};distance=0
        while center not in front:
            front={v for u in front for v in attached[u]};distance+=1;assert distance<=2
        corner_types=[]
        for omitted in row['black']:
            pair=[u for u in row['black']if u!=omitted]
            names=sorted(classify(g)for g in components(graph(row['white']+pair,row['edges'])))
            wanted=next(c['components']for c in row['corners']if c['pair']==pair)
            assert names==wanted;corner_types.append(dict(pair=pair,components=names))
        expected={2:([['D5','H4'],['A1','D8'],['D6','H3']],2),3:([['D5','H4'],['A1','E8'],['E6','H3']],1)}
        assert ([r['components']for r in corner_types],distance)==expected[root]
        out.append(dict(attachment_root=root,D5_endpoint_kind='long'if distance==2 else'short',distance_to_branch=distance,corners=corner_types))
    terms=doc['certificate']['terms']
    dominant_triangle=max((t for t in terms if t['metadata']['kind']=='triangle'),key=lambda t:F(t['multiplier']))
    assert dominant_triangle['metadata']['family']=='core51_attach_D5_root2'
    dominant_vertex=max((t for t in terms if t['metadata']['kind']=='vertex'),key=lambda t:F(t['multiplier']))
    assert dominant_vertex['metadata']['marked_component']=='D5'and dominant_vertex['metadata']['facet_root']==2
    paths=[f'scripts/{TAG}.py',f'results/{PRIOR}.json',f'results/{PRIOR}.md']
    return dict(status='PASS_LITERAL_D5_LONG_SHORT_PROSE_CORRECTION',sources={x:sha256((ROOT/x).read_bytes()).hexdigest()for x in paths},
        records=out,dominant_triangle=dominant_triangle['metadata'],dominant_vertex=dominant_vertex['metadata'],
        correction='The dominant attachment root2 is the D5 long endpoint and has an H3+D6 corner. Root3 is a short leaf and has the H3+E6 corner. The dominant D5 marked vertex likewise uses the long endpoint. The frozen source and exact certificate already contain the correct literal data; only the report prose had these identifications reversed.',
        scope='A graph-theoretic diagnostic correction only; no change to the exact lower bound, model, or face-coverage scope.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();r=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert r==json.loads(p.read_text())
    else:p.write_text(json.dumps(r,indent=2)+'\n')
    print(r['status'])
