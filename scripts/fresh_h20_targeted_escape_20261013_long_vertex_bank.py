#!/usr/bin/env python3
"""The complementary6300 original finite vertices with one norm-six wall.

Exact support reduction: a long wall cannot have two pi/6 neighbors in an
elliptic diagram; its one neighbor must be isolated in the short subdiagram.
The excluded3x3 principal matrices have nonpositive determinant, so this
reduction needs no floating inertia or external finite-diagram classifier.
"""
import argparse,hashlib,itertools,json
from fractions import Fraction as F
from pathlib import Path
import fresh_h20_targeted_escape_20261013_vertex_bank as V
ROOT=Path(__file__).resolve().parents[1]
VERTICES='data/enumerations/fresh_h20_targeted_escape_20261013_long_vertices.json'
OUT='data/enumerations/fresh_h20_targeted_escape_20261013_long_vertex_bank.json'
BASE=V.BASE

def vertices():
 b=V.read(BASE);Q=[[F(x[0])for x in row[:20]]for row in b['form'][:20]]
 R=[[F(x[0])for x in r[:20]]for r in b['old_normals']];h=[F(x[0])for x in b['strict_center'][:20]]
 dot=lambda a,b:sum((x*y for x,y in zip(a,b,strict=True)),F(0))
 mv=lambda M,x:[dot(row,x)for row in M]
 cov=[mv(Q,r)for r in R];G=[[dot(a,c)for c in cov]for a in R];hc=mv(Q,h);rhs=[dot(r,hc)for r in R]
 assert all(G[i][i]==2 for i in range(24))and all(G[i][i]==6 for i in range(24,36))
 assert all(G[i][j]in(0,-1)for i in range(24)for j in range(i))
 assert all(G[i][j]in(0,-3)for i in range(24,36)for j in range(24))
 # Any two long walls already give a non-PD2x2principal matrix.
 assert all(G[i][j]<=-6 for i in range(24,36)for j in range(24,i))
 # Two neighbors: det[[6,-3,-3],[-3,2,p],[-3,p,2]]<0 for p=0,-1.
 # One neighbor connected to another short node: principaldet0 if
 # thatnode has zero longpairing, andnegative if it too is a neighbor.
 assert [24-36-6*p*p+18*p for p in(0,-1)]==[-12,-36]
 assert 6*(4-1)-9*2==0
 count=0;out=[];bylong={}
 for L in range(24,36):
  N=[j for j in range(24)if G[L][j]<0];assert len(N)==3
  orth=[j for j in range(24)if G[L][j]==0]
  supports=list(itertools.combinations(orth,18))
  for a in N:
   permitted=[j for j in orth if G[a][j]==0]
   supports.extend(tuple(sorted((a,)+S))for S in itertools.combinations(permitted,17))
  assert len(set(supports))==len(supports);bylong[str(L)]={'candidate_supports':len(supports),'accepted':0}
  for S in supports:
   count+=1;ids=tuple(S)+(L,);sol=V.positive_solve([[G[i][j]for j in ids]for i in ids],[rhs[i]for i in ids])
   if sol is None:continue
   x=[h[k]-sum((s*R[i][k]for s,i in zip(sol,ids,strict=True)),F(0))for k in range(20)]
   xx=dot(x,mv(Q,x));assert xx<0;ps=[dot(x,c)for c in cov]
   assert all(p==0 if i in ids else p<0 for i,p in enumerate(ps))
   bylong[str(L)]['accepted']+=1
   out.append({'incident_walls':list(ids),'long_wall':L,'projection':list(map(str,x)),
    'norm':str(xx),'gram_solution':list(map(str,sol))})
 assert len(out)==6300
 d={'status':'EXACT_COMPLETE_ONE_LONG_FINITE_VERTEX_DIRECTION_CENSUS','source':BASE,'source_sha256':V.sha(BASE),
  'reduced_candidate_supports':count,'positive_definite_subsets':len(out),'by_long_wall':bylong,'vertices':out,
  'support_reduction':'Two pi6neighbors give negative3x3determinant. Onepi6neighbor joinedtoanothershortwall gives zero/negative3x3determinant. Thus thelongwall has zero neighbors, oroneisolatedshortneighbor. Every remaining18short+longsupport is enumerated exactly.',
  'scope':'Complete finitevertexdirection class withoneoriginalnorm6wall; togetherwiththe2046short-onlydirections this matchesall8346archivedfinitevertices. No cap or volumeconclusion follows fromthecensus.'}
 (ROOT/VERTICES).write_text(json.dumps(d,indent=2)+'\n');print({'vertices':len(out),'supports':count},flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--vertices',action='store_true');p.add_argument('--bank',action='store_true');a=p.parse_args()
 if a.vertices:vertices()
 if a.bank:
  V.VERTICES=VERTICES;V.OUT=OUT;V.bank()
