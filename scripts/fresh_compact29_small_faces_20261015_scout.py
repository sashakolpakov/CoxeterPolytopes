#!/usr/bin/env python3
"""Bounded constructive scouting of triangle stars; no exclusion claim."""
import json
from itertools import combinations, product
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact29_small_faces_20261015_scout'


def diagram(kind, n):
    e = {(i, i+1): 3 for i in range(n-1)}
    if kind == 'B': e[0, 1] = 4
    if kind == 'F': e[1, 2] = 4
    if kind == 'H': e[0, 1] = 5
    if kind == 'D':
        e = {(i, i+1): 3 for i in range(n-3)}
        e[n-3,n-2] = e[n-3,n-1] = 3
    if kind == 'E':
        e = {(i, i+1): 3 for i in range(n-2)}
        e[2,n-1] = 3
    return e


def components(nodes, e):
    left = set(nodes); out = []
    adj = {i: set() for i in nodes}
    for (i,j),m in e.items():
        if i in adj and j in adj and m > 2:
            adj[i].add(j); adj[j].add(i)
    while left:
        todo = [min(left)]; left.remove(todo[0])
        for i in todo:
            for j in sorted(adj[i] & left):
                todo.append(j); left.remove(j)
        out.append(todo)
    return out, adj


def finite(nodes, e):
    cc, adj = components(nodes, e)
    for c in cc:
        n = len(c)
        if n <= 2: continue
        ee = {p:m for p,m in e.items() if set(p) <= set(c) and m > 2}
        if len(ee) != n-1 or max(len(adj[i]) for i in c) > 3: return False
        forks = [i for i in c if len(adj[i]) == 3]
        if not forks:
            i = next(i for i in c if len(adj[i]) == 1); prev = None; seq=[]
            while True:
                nxt = adj[i] - ({prev} if prev is not None else set())
                if not nxt: break
                j = next(iter(nxt)); seq.append(ee[tuple(sorted((i,j)))]); prev,i=i,j
            heavy = [(i,m) for i,m in enumerate(seq) if m != 3]
            if not heavy: continue
            if len(heavy) != 1: return False
            i,m = heavy[0]
            if m == 4 and (i in [0,n-2] or (n == 4 and i == 1)): continue
            if m == 5 and n in [3,4] and i in [0,n-2]: continue
            return False
        if len(forks) != 1 or any(m != 3 for m in ee.values()): return False
        arms=[]; center=forks[0]
        for i in adj[center]:
            prev=center; length=1
            while adj[i]-{prev}:
                j=next(iter(adj[i]-{prev})); prev,i=i,j; length+=1
            arms.append(length)
        arms.sort()
        if arms[:2] == [1,1]: continue
        if tuple(arms) in [(1,2,2),(1,2,3),(1,2,4)]: continue
        return False
    return True


def extensions(S, e, new):
    cc,_ = components(S,e); where={i:k for k,c in enumerate(cc) for i in c}
    out=[]
    for k in range(4):
        for support in combinations(S,k):
            if len({where[i] for i in support}) != k: continue
            choices = [(3,3,3)] if k == 3 else product([3,4,5],repeat=k)
            for ms in choices:
                ee=e | {(min(i,new),max(i,new)):m for i,m in zip(support,ms)}
                if finite(S+[new],ee): out.append({i:m for i,m in zip(support,ms)})
    for c in cc:
        if len(c) == 1:
            for m in [6,7,8,10,12]: out.append({c[0]:m})
    return out


def gram(n,e):
    g=2*np.eye(n)
    for (i,j),m in e.items(): g[i,j]=g[j,i]=-2*np.cos(np.pi/m)
    return g


def run():
    out={}
    for kind,n in [('A',4),('B',4),('D',4),('F',4),('H',4),('E',7),('E',8)]:
        e=diagram(kind,n); assert finite(list(range(n)),e)
        rows=[]
        for a,b in combinations(range(n),2):
            W=[i for i in range(n) if i not in [a,b]]
            options=extensions(W+[a],e,n); candidates=[]
            for opt in options:
                for mb in [2,3,4,5,6,7,8,10,12]:
                    ee=e | {(i,n):m for i,m in opt.items()}
                    if mb>2: ee[b,n]=mb
                    if not finite(W+[b,n],ee): continue
                    g=gram(n+1,ee)
                    eig=np.linalg.eigvalsh(g)
                    if eig[0] >= -1e-9 or eig[1] <= 1e-9: continue
                    parabolic=False
                    for size in range(3,n+1):
                        for inds in combinations(range(n+1),size):
                            ev=np.linalg.eigvalsh(g[np.ix_(inds,inds)])
                            if abs(ev[0])<1e-9 and ev[1]>1e-9:
                                parabolic=True; break
                        if parabolic: break
                    if parabolic: continue
                    # Negative unit-sum dual provides a chamber center if all
                    # inverse row sums are negative. This is sufficient only.
                    I=np.linalg.inv(g); center=-I@np.ones(n+1)
                    centerq=float(center@g@center)
                    # Three vertex rays: perpendicular to W and adjacent pair.
                    face=[]
                    for missing in [a,b,n]:
                        x=-I[:,missing]; face.append(x)
                    dots=[[float(x@g@y) for y in face] for x in face]
                    if not all(dots[i][i]<-1e-9 for i in range(3)): continue
                    if not all(dots[i][j]<-1e-9 for i in range(3) for j in range(3)): continue
                    candidates.append(dict(new_contacts=sorted([[i,m] for i,m in opt.items()]+([[b,mb]] if mb>2 else [])),
                                           eigenvalues=eig.tolist(),center_square=centerq,vertex_pairings=dots))
            rows.append(dict(marked=[a,b],triangle_count=len(candidates),triangle=candidates[0] if candidates else None))
        out[kind+str(n)]=rows
        print(kind+str(n),sum(bool(r['triangle']) for r in rows),'/',len(rows),flush=True)
    (ROOT/f'results/{TAG}.json').write_text(json.dumps(dict(scope='Bounded numeric existence scout only; no negative classification.',diagrams=out),indent=2)+'\n')


if __name__ == '__main__': run()
