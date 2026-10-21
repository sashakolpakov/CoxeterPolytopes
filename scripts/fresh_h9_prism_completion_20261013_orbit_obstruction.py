#!/usr/bin/env python3
"""Exact common escape for all2304 seed images of the six source caps."""
import argparse,json,sys,itertools
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_prism_completion_20261013_verify as T
V=T.V;K=T.K
PREFIX=T.PREFIX+'_orbit_obstruction'

def check(d):
    assert d['source']==T.PREFIX+'_caps.json'and d['source_sha256']==V.sha(d['source'])
    original=V.read(d['source']);certificate=T.check(original)
    G,I,h,g,old,names,new,center,_,_=T.reconstruct(original)
    seed=[]
    for i in range(16):
        s=[0]*10
        if i<8:s[i]=1
        else:
            j=i-8;block=j//4;s[8+block]=1
            if j%4:s[4*block]=1;s[j]=-1
        seed.append(s)
        assert old[i]==list(map(V.rat,s))
    sg=[[K.pair(tuple(K.scale(K.O,x)for x in a),G,tuple(K.scale(K.O,x)for x in b))for b in seed]for a in seed]
    A,B,C=sg[0][8],sg[1][8],sg[4][8]
    assert len({A,B,C})==3 and all(K.sign(K.add(x,K.scale(K.O,2)))<0 for x in(A,B,C))
    for i in range(16):
        for j in range(16):
            if i//8==j//8:expected=K.scale(K.O,2)if i==j else K.Z
            else:
                u,v=i%8,j%8;expected=A if u==v else B if u//4==v//4 else C
            assert sg[i][j]==expected
    # The two maximal orthogonal8-cliques determine the layers. The unique
    # A-colored matching fixes a common permutation of their indices. The
    # B-colored pairs determine the4+4 partition. Hence every weighted-Gram
    # automorphism is precisely one of the2304 enumerated permutations.
    y=tuple(map(K.parse,d['point']));assert all(x[1]==0 for x in y)
    norm=K.pair(y,I,y);future=K.pair(y,I,h);assert K.sign(norm)>0 and K.sign(future)<0
    yi=list(map(V.golden,y));oldpairs=[V.dot(r,yi)for r in old+new];assert min(x[0]for x in oldpairs)>0
    basis=list(range(8))+[8,12];sources=[old[names.index(n)]for n in T.SOURCE_NAMES]
    seedcenter=[sum(row[j]for row in seed)for j in range(10)];cq=tuple(K.scale(K.O,x)for x in seedcenter)
    center_norm=K.pair(cq,G,cq);assert K.sign(center_norm)<0 and K.sign(K.dot(cq,h))>0
    center_future=K.neg(K.dot(cq,y));assert K.sign(center_future)<0
    ci=list(map(V.rat,seedcenter));center_old=[V.pair(r,g,ci)for r in old+new];assert max(x[1]for x in center_old)<0
    center_source=[V.pair(r,g,ci)for r in sources];assert max(x[1]for x in center_source)<0
    seen=set();minima=[None]*6;worst=[None]*6;signcount=0
    for p,q,bs,uv in itertools.product(itertools.permutations(range(4)),itertools.permutations(range(4)),range(2),range(2)):
        perm=tuple(((i//8)^uv)*8+(((i%8)//4)^bs)*4+(p if (i%8)//4==0 else q)[i%4]for i in range(16))
        assert perm not in seen;seen.add(perm)
        assert all(sg[perm[i]][perm[j]]==sg[i][j]for i in range(16)for j in range(16))
        M=[seed[perm[i]]for i in basis]
        assert all([sum(a*M[k][j]for k,a in enumerate(seed[i]))for j in range(10)]==seed[perm[i]]for i in range(16))
        assert [sum(a*M[k][j]for k,a in enumerate(seedcenter))for j in range(10)]==seedcenter
        # r transforms as the row rM; (rM).y = r.(My), with My rational.
        transformed=[V.rat(sum(Q(a)*y[j][0]for j,a in enumerate(row)))for row in M]
        for i,r in enumerate(sources):
            pp=V.dot(r,transformed);assert pp[0]>0,(perm,T.SOURCE_NAMES[i],pp)
            signcount+=1
            if minima[i]is None or pp[0]<minima[i][0]:minima[i]=pp;worst[i]=list(perm)
    assert len(seen)==2304 and signcount==13824
    return {'status':'PASS_EXACT_FULL_SEED_ORBIT_COMMON_ESCAPE','source_sha256':V.sha(T.PREFIX+'_caps.json'),
            'seed_weighted_gram_group_order':2304,'literal_orbit_entries':13824,'source_normals':T.SOURCE_NAMES,
            'point':d['point'],'norm':K.enc(norm),'future_center_pairing':K.enc(future),
            'fixed_strict_timelike_center':{'integer_seed_sum':seedcenter,'norm':K.enc(center_norm),'escape_pairing':K.enc(center_future),
                'old102_pairings':list(map(V.enc,center_old)),'six_source_pairings':list(map(V.enc,center_source)),
                'proof':'Every group matrix fixes the seed sum. Its strict six-source inequalities transfer to all13824 images; all old102 inequalities are also strict.'},
            'old102_strict_pairings':list(map(V.enc,oldpairs)),
            'orbit_minima':[{'source':n,'smallest_lower_bound_pairing':V.enc(pp),'attaining_seed_permutation':perm}for n,pp,perm in zip(T.SOURCE_NAMES,minima,worst,strict=True)],
            'seed_integer_normals':seed,'seed_cross_gram_values':list(map(K.enc,(A,B,C))),
            'scope':'The same strict future spacelike point satisfies every old102 inequality and every one of the13824 literal images of compact0..3,facePrism,faceGap under the complete weighted16-seed Gram automorphism group. Any selection from that union, including arbitrary old-wall deletions, has infinite volume if it defines a full-dimensional hyperbolic polyhedron. No claim about caps outside this finite image family.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=V.read(PREFIX+'_point.json');out=check(d)
    if a.verify:assert out==V.read(PREFIX+'.json')
    else:(V.ROOT/(PREFIX+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['literal_orbit_entries'])
