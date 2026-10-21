#!/usr/bin/env python3
"""Exact finite spherical-extension profiles, with isolated I2 tails explicit.

Uses Fraction Gram inversion and separate squarefree-radical signs; imports
the frozen exact860 reconstruction, never the numerical search producer.
"""
import argparse,json
from itertools import combinations
from collections import Counter
from fractions import Fraction as F
import fresh_h22_contact_bridge_20261013_verify as E
R=E.R;TAG='fresh_h22_nonright_bridge_20261013';OUT='results/'+TAG+'_profiles.json'

def cosine(m):
    return {3:R.C(F(1,2)),4:R.sqrtq(F(1,2)),
            5:R.scale(R.add(R.C(1),R.sqrtq(5)),F(1,4)),
            6:R.scale(R.sqrtq(3),F(1,2))}[m]

def run(verify=False):
    source,M,W,cov,n,h,points=E.setup();rows=[];counts=Counter()
    for vi,src in enumerate(source['compact_H2_vertices']):
        inc=src['incident'];A=[]
        for j in inc:
            row=[]
            for k in inc:
                a=R.pair(W[j],cov[k]);assert set(a)<={1};row.append(a.get(1,F(0)))
            A.append(row)
        R.pd(A);AI=E.inverse(A)
        for k,released in enumerate(inc):
            ix=[i for i in range(22)if i!=k];edge=[inc[i]for i in ix]
            GI=[[AI[i][j]-AI[i][k]*AI[k][j]/AI[k][k]for j in ix]for i in ix]
            assert all(x>=0 for row in GI for x in row)
            choices=[];profile=[]
            for j in range(21):
                for m in (3,4,5,6):
                    p=R.scale(R.mul(R.sqrtq(2*n[edge[j]]),cosine(m)),-1)
                    q=R.scale(R.mul(p,p),GI[j][j])
                    if R.sign(E.sub(R.C(2),q))>0:choices.append((j,m,p))
            for links in (1,2,3):
                selected=choices if links<3 else [x for x in choices if x[1]==3]
                for pr in combinations(selected,links):
                    if len({x[0]for x in pr})!=links:continue
                    q=R.psum(R.scale(R.mul(p,s),GI[j][i])for j,m,p in pr for i,l,s in pr)
                    residual=E.sub(R.C(2),q)
                    if R.sign(residual)<=0:continue
                    counts[str(links)+'link']+=1
                    profile.append({'contacts':[[edge[j],m]for j,m,p in pr],
                                    'schur_residual':R.enc(residual)})
            isolated=[edge[j]for j,i in enumerate(ix)if all(A[i][l]==0 for l in ix if l!=i)]
            # Every isolated existing A1 admits a one-link I2(m) extension for
            # every m>=7: residual=2(1-cos(pi/m)^2)>0. No bound is imposed here.
            counts['unbounded_dihedral_families']+=len(isolated)
            rows.append({'vertex':vi,'released':released,'edge':edge,
                         'finite_profiles':profile,'I2_labels_at_least7':isolated})
    assert counts=={'1link':1951,'2link':3706,'3link':823,'unbounded_dihedral_families':217},counts
    out={'status':'PASS_EXACT_FINITE_SPHERICAL_PROFILE_CENSUS_WITH_I2_TAILS_EXPLICIT',
         'sources':{p:R.sha(p)for p in(E.SOURCE,'scripts/'+TAG+'_profiles.py','scripts/'+E.TAG+'_verify.py','scripts/fresh_h22_multiquadratic_20261013_review.py')},
         'counts':dict(counts),'rows':rows,
         'scope':'Positive-definite rank22 Coxeter extensions of the88 specified PD21 supports, with a new norm2 root. Finite irreducible Coxeter classification bounds the new-node valence by3; valence3 has labels3, valence2 labels at most5, and valence1 label>=7 is possible only on an isolated oldA1. Thus the6480 explicit profiles plus217 parametric I2(m>=7) families exhaust abstract nonright spherical extension profiles. This is not existence of compatible geometric caps, nor exclusion of these parametric tails.'}
    if verify:assert R.read(OUT)==out
    else:(R.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],dict(counts),flush=True);return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');run(ap.parse_args().verify)
