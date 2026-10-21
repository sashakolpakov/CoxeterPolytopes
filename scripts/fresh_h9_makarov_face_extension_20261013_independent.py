#!/usr/bin/env python3
"""Independent squarefree-radical replay of the 12-wall Makarov-face seed.

No arithmetic, matrix, or geometry routine from the producer is imported.
"""
import argparse
import json
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h22_multiquadratic_20261013_review as R

TAG='fresh_h9_makarov_face_extension_20261013'
SOURCE='results/'+TAG+'.json'
SCRIPT='scripts/'+TAG+'.py'
OUT='results/'+TAG+'_independent.json'
RADROOT={11:F(16,11),55:F(4,11)}
def k(x):return R.add(R.C(F(x[0])),{5:F(x[1])} if F(x[1]) else {})
def l(x):return R.add(k(x[0]),R.mul(k(x[1]),RADROOT))
def sub(x,y):return R.add(x,R.scale(y,-1))
def invk(x):
    assert set(x)<=set((1,5))
    a,b=x.get(1,F(0)),x.get(5,F(0));q=a*a-5*b*b;assert q
    return R.add(R.C(a/q),{5:-b/q} if b else {})
def divk(x,y):return R.mul(x,invk(y))
def mv(A,v):return [R.psum(R.mul(a,b) for a,b in zip(row,v,strict=True)) for row in A]
def pair(a,Q,b):return R.pair(a,mv(Q,b))
def inertia(A):
    A=[[dict(x) for x in row]for row in A];p=n=0;piv=[]
    for j in range(len(A)):
        d=A[j][j];s=R.sign(d);assert s;p+=s>0;n+=s<0;piv.append(R.enc(d))
        for u in range(j+1,len(A)):
            for v in range(j+1,len(A)):
                A[u][v]=sub(A[u][v],divk(R.mul(A[u][j],A[j][v]),d))
    return [p,n],piv

def run(verify=False):
    d=R.read(SOURCE);cap=d['changed_center_cap']
    assert d['status']=='EXACT_COMPACT_H5_FACE_H9_SEED_AND_CHANGED_CENTER_CUTTER'
    assert R.sign(RADROOT)>0 and R.mul(RADROOT,RADROOT)==k(cap['extension_radicand'])
    Q=[[k(z) for z in row]for row in d['form']];assert len(Q)==10 and all(len(row)==10 for row in Q)
    assert Q==list(map(list,zip(*Q)));signature,piv=inertia(Q);assert signature==[9,1]
    roots=[[k(z) for z in row]for row in d['seed_normals']]
    assert len(roots)==11;roots.append([l(z) for z in cap['normal']]);assert all(len(r)==10 for r in roots)
    norms=[pair(r,Q,r)for r in roots]
    assert all(n==R.C(2)for i,n in enumerate(norms)if i!=6)
    assert norms[6]=={1:F(7,11),5:F(-1,11)} and R.sign(norms[6])>0
    assert norms[-1]==l(cap['norm'])
    cosines={'m2':{},'m3':R.C(F(1,4)),'m5':{1:F(3,8),5:F(1,8)}}
    types=[];counts=Counter()
    for i,r in enumerate(roots):
        for j in range(i):
            p=pair(r,Q,roots[j]);assert R.sign(p)<=0
            sq=R.mul(p,p);den=R.mul(norms[i],norms[j])
            if not p:kind='m2'
            elif R.sign(sub(sq,den))>0:kind='dotted'
            else:
                kind=next((name for name,c in cosines.items()if sq==R.mul(den,c)),None)
                assert kind,(i,j,R.enc(p))
            counts[kind]+=1
            if i==11:types.append(kind)
    assert types==cap['old_pair_types'] and sum(counts.values())==66
    center=list(map(k,d['strict_center']));cn=pair(center,Q,center)
    assert cn==k(d['strict_center_norm']) and R.sign(cn)<0
    centerpairs=[pair(r,Q,center)for r in roots];assert all(R.sign(z)<0 for z in centerpairs)
    assert centerpairs[-1]==l(cap['center_pair'])
    vertices=[]
    assert len(d['compact_H5_vertices'])==len(d['vertex_active_sets'])==10
    for row,oldactive in zip(d['compact_H5_vertices'],d['vertex_active_sets'],strict=True):
        v=list(map(k,row))+[{} for _ in range(4)];qn=pair(v,Q,v)
        ps=[pair(v,Q,r)for r in roots];active=[j for j,p in enumerate(ps)if not p]
        assert active==oldactive+list(range(7,11)) and len(active)==9
        assert R.sign(qn)<0 and R.sign(pair(v,Q,center))<0 and all(R.sign(p)<=0 for p in ps)
        assert R.sign(ps[-1])<0
        gram=[[pair(roots[i],Q,roots[j])for j in active]for i in active]
        assert inertia(gram)[0]==[9,0]
        vertices.append({'active':active,'norm':R.enc(qn),'new_cap_pairing':R.enc(ps[-1])})
    # The strict center projection exhibits a relative-interior point of the
    # added facet. Acute pairings similarly preserve every old true facet.
    f=list(map(l,cap['actual_facet_point']));fn=pair(f,Q,f)
    assert R.sign(fn)<0 and not pair(f,Q,roots[-1]) and R.sign(pair(f,Q,center))<0
    assert all(R.sign(pair(f,Q,w))<0 for w in roots[:-1])
    expected=[sub(a,R.scale(R.mul(centerpairs[-1],b),F(1,2)))for a,b in zip(center,roots[-1],strict=True)]
    assert f==expected

    old=d['all_label_vertex_center_family_obstruction'];E=list(map(R.C,old['strict_integer_escape']))
    eq=pair(E,Q,E);assert eq==k(old['norm']) and R.sign(eq)>0
    assert pair(E,Q,center)==k(old['center_pair']) and R.sign(pair(E,Q,center))<0
    assert all(R.sign(pair(E,Q,w))<0 for w in roots[:-1])
    cut=pair(E,Q,roots[-1]);assert cut==l(cap['old_escape_cut']) and R.sign(cut)>0
    EE=list(map(R.C,cap['new_strict_integer_escape']));newq=pair(EE,Q,EE)
    assert newq==k(cap['new_escape_norm']) and R.sign(newq)>0 and R.sign(pair(EE,Q,center))<0
    eps=[pair(EE,Q,r)for r in roots];assert all(R.sign(p)<0 for p in eps)
    assert eps[-1]==l(cap['new_escape_cap_pair'])

    H=[row[6:]for row in Q[6:]];J=[[k(z)for z in row]for row in d['H4_inverse']]
    assert len(J)==4 and all(len(row)==4 for row in J)
    for i in range(4):
        for j in range(4):assert R.psum(R.mul(H[i][t],J[t][j])for t in range(4))==R.C(int(i==j))
    assert inertia(H)[0]==[4,0]
    G=[row[:6]for row in Q[:6]];proofs={(p['vertex'],p['h4_index']):p for p in old['forty_strict_m3_checks']}
    assert len(proofs)==40;gaps=[]
    for vi,row in enumerate(d['compact_H5_vertices']):
        v=list(map(k,row));qn=pair(v,G,v);a=pair(v,G,E[:6]);assert R.sign(a)<0
        for j in range(4):
            kap=divk(sub(J[j][j],R.C(2)),R.scale(qn,-1));assert R.sign(kap)>0 and R.sign(E[6+j])<0
            gap=sub(R.mul(kap,R.mul(a,a)),R.mul(E[6+j],E[6+j]));p=proofs[(vi,j)]
            assert a==k(p['old_pair']) and gap==k(p['squared_strict_gap']) and R.sign(gap)>0
            gaps.append({'vertex':vi,'h4_index':j,'gap':R.enc(gap)})
    out={'status':'PASS_INDEPENDENT_12_WALL_MAKAROV_FACE_GEOMETRY',
         'sources':{p:R.sha(p)for p in [SOURCE,SCRIPT,'scripts/fresh_h22_multiquadratic_20261013_review.py','scripts/'+Path(__file__).name]},
         'radical_square_root':R.enc(RADROOT),'signature':signature,'LDL_pivots':piv,
         'walls':12,'norms':list(map(R.enc,norms)),'pair_counts':dict(counts),
         'center_norm':R.enc(cn),'center_pairings':list(map(R.enc,centerpairs)),
         'actual_finite_vertices':vertices,'new_facet_point_norm':R.enc(fn),
         'old_escape_norm':R.enc(eq),'old_escape_new_cap_pairing':R.enc(cut),
         'new_escape_norm':R.enc(newq),'new_escape_center_pairing':R.enc(pair(EE,Q,center)),
         'all_12_new_escape_pairings':list(map(R.enc,eps)),
         'all_40_vertex_weight_family_strict_gaps':gaps,
         'scope':'Separate sparse squarefree-radical arithmetic verifies all 12 norms, all 66 pairs, strict center, ten actual rank-nine finite vertices, a relative-interior new-facet point, the old escape cut and new strict escape. Also independently verifies the forty squared inequalities used in the all-label vertex/weight-family obstruction. It does not repeat the compact H5 extreme-ray census or the 200-profile classification and makes no finite-volume or novelty claim.'}
    if verify:assert out==R.read(OUT)
    else:(R.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
