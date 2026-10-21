#!/usr/bin/env python3
"""Rigorous nine-plane obstruction to another mutually acute cap pair."""
import argparse
import json
from itertools import combinations
from scripts import fresh_h9_remote_next_20261013_bridge_verify as B
V=B.V;K=V.K
OUT=B.PREFIX+'_projection_obstruction.json'
INCIDENT=['seed1','seed2','seed4','seed6','seed8','wall16','wall17','wall35','nextH4']
WITNESSES={'nextH4':['wall46','nextH4'],'wall35':['wall35','wall37'],'wall17':['seed3','wall17'],
           'wall16':['wall16','seed5'],'seed8':['seed8','seed3'],'seed6':['seed5','seed6'],
           'seed4':['seed4','wall46'],'seed2':['wall60','seed2'],'seed1':['wall60','seed1']}

def inverse(M):
    n=len(M);A=[list(row)+[V.rat(int(i==j))for j in range(n)]for i,row in enumerate(M)];pivots=[]
    for j in range(n):
        candidates=[i for i in range(j,n)if A[i][j][0]>0 or A[i][j][1]<0]
        assert candidates,'No certified nonzero Gaussian pivot'
        i=max(candidates,key=lambda i:min(abs(x)for x in A[i][j]));A[i],A[j]=A[j],A[i]
        p=A[j][j];pivots.append(V.enc(p));A[j]=[V.div(x,p)for x in A[j]];A[j][j]=V.rat(1)
        for i in range(n):
            if i==j:continue
            t=A[i][j];A[i]=[V.sub(x,V.mul(t,y))for x,y in zip(A[i],A[j],strict=True)];A[i][j]=V.rat(0)
    # Gaussian identities justify exact zeros/ones in the eliminated block.
    return[row[n:]for row in A],pivots

def run():
    source=B.PREFIX+'_caps.json';data=V.read(source);B.check(data)
    G,I,h,g,roots,names,caps,center,_=B.reconstruct(data);roots+=caps;names+=['bridge0','bridge1']
    target=tuple(map(K.parse,data['escape']));yi=list(map(V.golden,target));rows=[]
    # These named sets came from a numerical vertex. The obstruction below
    # applies to their exact spans, without requiring exact incidence of that
    # numerical proposal or of the deliberately strict saved escape.
    for omitted in INCIDENT:
        block=[n for n in INCIDENT if n!=omitted];A=[roots[names.index(n)]for n in block]
        M=[[V.pair(a,g,b)for b in A]for a in A];Mi,pivots=inverse(M)
        neg=next((i,j)for i,j in combinations(range(8),2)if M[i][j][1]<-2*V.SCALE)
        nq=V.add(V.rat(4),V.mul(V.rat(2),M[neg[0]][neg[1]]));assert nq[1]<0
        # A nonsingular 8D span with a negative direction inside signature(9,1)
        # has signature(7,1). Its orthogonal complement is positive definite2D.
        def project(x):
            rhs=[V.pair(a,g,x)for a in A];coeff=V.mv(Mi,rhs)
            return[V.sub(x[j],V.total(V.mul(c,a[j])for c,a in zip(coeff,A,strict=True)))for j in range(10)]
        witness=WITNESSES[omitted];u,v=[project(roots[names.index(n)])for n in witness]
        d=V.pair(u,g,u);e=V.pair(v,g,v);p=V.pair(u,g,v)
        assert d[0]>0 and e[0]>0 and p[1]<0
        det=V.sub(V.mul(d,e),V.mul(p,p));assert det[0]>0
        rows.append({'omitted':omitted,'span_walls':block,'span_signature':[7,1],
                     'inverse_nonzero_pivots':pivots,'negative_direction_walls':[block[i]for i in neg],
                     'negative_direction_norm':V.enc(nq),'projected_witness_walls':witness,
                     'projected_gram':[[V.enc(d),V.enc(p)],[V.enc(p),V.enc(e)]],
                     'projected_determinant':V.enc(det)})
    return {'status':'PASS_NINE_POSITIVE_PLANES_HAVE_NO_MUTUALLY_ACUTE_PAIR',
            'source':source,'source_sha256':V.sha(source),'ambient_signature':[9,1],
            'selected_nine_walls':INCIDENT,'strict_escape_pairings':[V.enc(V.dot(roots[names.index(n)],yi))for n in INCIDENT],
            'planes':rows,'proof':'Each exact eight-wall span is nonsingular and contains a negative direction, so its orthogonal complement is a positive2-plane. In that plane two projected retained walls u,v are independent and have u.v<0. Their positive-definite2x2 Gram matrix has an entrywise strictly positive inverse. For any two nonzero vectors r,s satisfying (r,u),(r,v),(s,u),(s,v)<=0, expressing their pairing via that inverse makes (r,s)>0. Equivalently, the two nonpositive halfplanes meet in a sector of opening strictly smaller than pi/2. Thus no mutually acute-compatible pair of nonzero new normals exists in each of these nine exact planes, with no label bound.',
            'scope':'Only these nine specified orthogonal complements. Does not exclude one cap in a plane, caps in different planes, lower-rank orthogonality, replacement of a retained wall, or different incident spans. The saved strict escape is not asserted to be incident to the nine named walls.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert V.read(OUT)==out
    else:(V.ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['planes']),'exact planes')
