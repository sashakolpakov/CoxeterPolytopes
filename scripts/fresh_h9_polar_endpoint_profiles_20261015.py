#!/usr/bin/env python3
"""Complete finite endpoint profiles for supports with isolated I2(10)."""
from itertools import combinations,product
from collections import Counter
from fractions import Fraction as F
import fresh_h9_mixed_d4_endpoint_profiles_20261014 as P
Q,Z,O=P.Q,P.ZERO,P.ONE

def profiles(S,labels):
    edges={tuple(sorted((i,j))):labels[tuple(sorted((i,j)))]for i,j in combinations(S,2)}
    assert len(S)==8 and all(m in(2,3,5,10)for m in edges.values())
    comps,_=P.components(S,edges);types=[P.component_type(c,edges)for c in comps]
    inert=[c for c,t in zip(comps,types)if t=='I2(10)'];blocked={i for c in inert for i in c}
    assert all(len(c)==2 and edges[tuple(c)]==10 for c in inert)
    # Normalized norm-two inverse diagonal of I2(10) is 3+sqrt5>2.
    heavy_diagonal=Q(3,1);assert(heavy_diagonal-2).sign()>0
    J=[[Z for _ in S]for _ in S]
    for c,t in zip(comps,types):
        if t=='I2(10)':
            for i in c:J[S.index(i)][S.index(i)]=heavy_diagonal
            continue
        M=[[Q(2)if i==j else{2:Z,3:-O,5:-P.PHI}[edges[tuple(sorted((i,j)))]]for j in c]for i in c]
        P.pd(M);B=P.inverse(M)
        for a,i in enumerate(c):
            for b,j in enumerate(c):J[S.index(i)][S.index(j)]=B[a][b]
    assert all(v.sign()>=0 for r in J for v in r)
    where={i:k for k,c in enumerate(comps)for i in c};rows=[];rejected=Counter();tested=0
    allowed=[i for i in S if i not in blocked]
    for n in range(4):
        for touched in combinations(allowed,n):
            if len({where[i]for i in touched})!=n:continue
            for ms in(product((3,4,5),repeat=n)if n<3 else[(3,3,3)]):
                tested+=1;cost=sum((P.SQUARES[m]*J[S.index(i)][S.index(i)]for i,m in zip(touched,ms)),Z);rho=2-cost
                if rho.sign()<=0:rejected['affine'if not rho else'indefinite']+=1;continue
                es=dict(edges);es.update({(i,12):m for i,m in zip(touched,ms)})
                cs,_=P.components(S+[12],es);ts=[P.component_type(c,es)for c in cs]
                rows.append(dict(contacts=[[i,m]for i,m in zip(touched,ms)],projection_square=cost,residual_norm=rho,
                                 endpoint_components=ts,new_node_component=ts[next(k for k,c in enumerate(cs)if 12 in c)]))
    tails=[]
    for c in comps:
        if len(c)!=1:continue
        i=c[0];assert J[S.index(i)][S.index(i)]==F(1,2)
        tails.append(dict(wall=i,integer_labels='every integer m>=6',other_contacts='all right',
                          projection_square='2cos(pi/m)^2',residual_norm='2sin(pi/m)^2',new_node_component='I2(m)'))
    assert tested==len(rows)+sum(rejected.values())
    return dict(support=S,support_components=comps,support_component_types=types,
                inert_I2_10_components=inert,forbidden_attachment_nodes=sorted(blocked),
                effective_projection_inverse=J,
                inverse_scope='The I2(10) off-diagonal entries are omitted, not zero in the true inverse. Every admissible contact vector vanishes on those components, so this matrix gives the exact projection for every recorded profile.',
                finite_profiles=rows,dihedral_tail_profiles=tails,finite_trial_count=tested,rejected_trial_counts=dict(rejected),
                completeness='I2(10) inverse diagonal3+sqrt5>2 forbids every nonright attachment. On all other components, cycles and four-leaf stars are nonpositive; a new label>=6 requires an isolated A1 and no other contact. All remaining3/4/5 profiles are checked by exact Schur cost.')
