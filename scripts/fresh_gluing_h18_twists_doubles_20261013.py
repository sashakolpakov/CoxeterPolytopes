#!/usr/bin/env python3
"""Exact relative-twist and parent-double experiments for Vinberg's H18 mix."""
from __future__ import annotations
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import networkx as nx
from sympy import Matrix
from coxeter_polytopes.exact.silver import Silver as Q
from scripts.explore_vinberg_h18_mixture_faces import parents, graph

ALLOWED={Q(0),Q(Fraction(1,4)),Q(Fraction(1,2)),Q(Fraction(3,4)),Q(Fraction(1,2),Fraction(1,4))}

def auts(G,indices):
    g=graph(G,indices)
    return list(nx.algorithms.isomorphism.GraphMatcher(g,g,node_match=lambda a,b:a['norm']==b['norm'],edge_match=lambda a,b:a['pairing']==b['pairing']).isomorphisms_iter())

def even_walls(G):
    # These two parents have only finite labels 2,3,4.
    return [i for i in range(len(G)) if all(not G[i][j] or 2*G[i][j]**2>=G[i][i]*G[j][j] for j in range(len(G)) if i!=j)]

def restricted(G,S):
    M=Matrix([[G[i][j].a for j in S] for i in S])
    assert all(not G[i][j].b for i in S for j in S)
    if any(M[:k,:k].det()<=0 for k in range(1,len(S)+1)):return None,None
    R=M.inv();inv=[[Q(Fraction(x)) for x in row] for row in R.tolist()]
    cache={}
    def p(i,j):
        if (i,j) not in cache:cache[i,j]=G[i][j]-sum((G[i][s]*inv[a][b]*G[t][j] for a,s in enumerate(S) for b,t in enumerate(S)),Q(0))
        return cache[i,j]
    ids=[i for i in range(len(G)) if i not in S and p(i,i)>0]
    return ids,p

def face_gate(G,S,labels):
    ids,p=restricted(G,S)
    if ids is None:return {'elliptic':False}
    for i,j in combinations(ids,2):
        value=p(i,j);sq=value**2/(p(i,i)*p(j,j))
        if value>0 or (sq<1 and sq not in ALLOWED):
            return {'elliptic':True,'coxeter':False,'facet_count':len(ids),'witness':[labels[i],labels[j]],'pairing':value.to_json(),'square':sq.to_json()}
    nonlinear=[]
    for i,j in combinations(ids,2):
        value=4*p(i,j)**2/(p(i,i)*p(j,j))
        if value.b:nonlinear.append({'walls':[labels[i],labels[j]],'two_cycle':value.to_json()})
    return {'elliptic':True,'coxeter':True,'facet_count':len(ids),'nonrational_two_cycle':nonlinear[:1]}

def main():
    original=json.loads(Path('data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json').read_text())
    G=[[Q.from_json(x) for x in row] for row in original['gram']]
    labels=original['labels'];L,R=parents();a,b=original['construction']['gluing_facets_zero_based']
    mapping={int(i):j for i,j in original['construction']['matching'].items()}
    C=sorted(mapping)
    AG=auts(L,C);AL=auts(L,list(range(len(L))));AR=auts(R,list(range(len(R))))
    SL=[x for x in AL if x[a]==a];SR=[x for x in AR if x[b]==b]
    left_images={tuple(x[i] for i in C) for x in SL}
    right_inverse={j:i for i,j in mapping.items()}
    right_images={tuple(right_inverse[x[mapping[i]]] for i in C) for x in SR}
    full={tuple(x[i] for i in C) for x in AG}
    assert left_images==right_images==full
    twist={'common_facet_automorphism_order':len(AG),'left_automorphism_order':len(AL),'right_automorphism_order':len(AR),'left_seam_stabilizer_order':len(SL),'right_seam_stabilizer_order':len(SR),'both_restriction_images_full':True,'relative_twist_double_cosets':1}
    print(json.dumps(twist),flush=True)
    # Add the removed seam normal g, of norm one. The right seam normal is -sqrt(2)g.
    reservoir=[row[:] for row in G]
    seam=[]
    for label in labels:
        if label.startswith('L'):seam.append(Q(L[int(label[1:])-1][a]))
        else:seam.append(Q(0,-Fraction(R[int(label[1:])-1][b],2)))
    for i,row in enumerate(reservoir):row.append(seam[i])
    reservoir.append([*seam,Q(1)])
    n=len(reservoir)
    def unit(i,c=Q(1)):return {i:c}
    def bil(v,w):return sum((x*y*reservoir[i][j] for i,x in v.items() for j,y in w.items()),Q(0))
    base={label:unit(i) for i,label in enumerate(labels)}
    pl={i:(unit(40) if i==a else base[f'L{i+1}']) for i in range(len(L))}
    pr={j:(unit(40,Q(0,-1)) if j==b else pl[right_inverse[j]] if j in right_inverse else base[f'R{j+1}']) for j in range(len(R))}
    target_rows=[]
    for i in C:
        f=next(x for x in original['facets'] if x['facet']==f'L{i+1}')
        w=f['first_forbidden_pair']
        target_rows.append({'face':[f'L{i+1}'],'forbidden_pair':w['walls'],'dimension':17})
    for row in original['seam_faces']:
        target_rows.extend({**w,'dimension':row['dimension']} for w in row['witnesses'])
    output=[]
    for name,P,k,stabilizer,pvectors in [('left',L,a,SL,pl),('right',R,b,SR,pr)]:
        eligible=[i for i in even_walls(P) if i!=k]
        assert all(P[k][i]**2>=P[k][k]*P[i][i] for i in eligible)
        todo=set(eligible);orbits=[]
        while todo:
            chosen=min(todo);orbit={g[chosen] for g in stabilizer};todo-=orbit;orbits.append(sorted(orbit))
        for orbit in orbits:
            d=orbit[0];dlabel=f'{"L" if name=="left" else "R"}{d+1}';dvec=pvectors[d];dnorm=bil(dvec,dvec)
            wall_vectors={lab:v for lab,v in base.items() if lab!=dlabel}
            for i,v in pvectors.items():
                if i==d:continue
                factor=2*bil(v,dvec)/dnorm
                reflected=v.copy()
                for j,c in dvec.items():reflected[j]=reflected.get(j,Q(0))-factor*c
                reflected={j:c for j,c in reflected.items() if c}
                if not factor:continue
                wall_vectors[f'{name}_s{d+1}_{i+1}']=reflected
            # Reflection in a disjoint wall leaves the original common facet intact.
            names=list(wall_vectors);vectors=list(wall_vectors.values())
            D=[[bil(v,w) for w in vectors] for v in vectors]
            assert all(D[i][j]<=0 for i in range(len(D)) for j in range(i))
            assert all(D[i][j]**2>=D[i][i]*D[j][j] or D[i][j]**2/(D[i][i]*D[j][j]) in ALLOWED for i in range(len(D)) for j in range(i))
            retained=Counter();fresh=[]
            for row in target_rows:
                if dlabel not in row['forbidden_pair']:
                    retained[row['dimension']]+=1
                    continue
                selection=[names.index(x) for x in row['face']]
                fresh.append({'face':row['face'],'dimension':row['dimension'],'gate':face_gate(D,selection,names)})
            record={'parent':name,'doubling_facet':dlabel,'seam_stabilizer_orbit_zero_based':orbit,'chamber_facets':len(D),'old_witnesses_retained':dict(retained),'fresh_angle_tests':fresh,'coxeter_survivors':[x for x in fresh if x['gate'].get('coxeter')]}
            output.append(record)
            print(json.dumps(record),flush=True)
    result={'status':'EXACT_TWIST_AND_SINGLE_PARENT_DOUBLE_EXPERIMENT','twists':twist,'doubles':output,'source_parent_facets':[37,24],'scope':'All relative direct twists and all even single-parent doubles preserving the original seam, modulo the original seam stabilizers. Only original seam-crossing faces in dimensions17,16,15,13 are audited.'}
    Path('data/enumerations/fresh_gluing_h18_twists_doubles_20261013.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':main()
