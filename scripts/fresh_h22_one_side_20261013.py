#!/usr/bin/env python3
"""Exact symmetrization of the 25 cap types of the 839-facet branch.

Norms and three right-contact identities are inherited from frozen recipes.
Every newly compared strict relation is checked with outward intervals.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
import fresh_h22_incidence_prism_20261013_verify as V
from fresh_h22_pi7_coupled_caps_20261013 import ipair,ivec

TAG='fresh_h22_one_side_20261013'
OUT='results/'+TAG+'.json'
COUPLED='results/fresh_h22_pi7_coupled_caps_20261013.json'
CUSPS='results/fresh_h22_pi7_cusp_pair_20261013.json'

def enc(x):return list(map(str,x))
def flip(r):return r[:-1]+[V.neg(r[-1])]

def run(verify=False):
    cert=V.read(V.OUT);coupled=V.read(COUPLED);cusps=V.read(CUSPS)
    for obj in (cert,coupled,cusps):
        for p,s in obj['sources'].items():assert V.sha(p)==s
    M,R,h=V.setup();recipes=V.read(V.REC)['recipes'];raw=[];norms=[]
    for rec in recipes:
        u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm'])
        qu=V.pair(M,u,u);qv=V.pair(M,v,v)
        assert qv>0 and V.pair(M,u,v)==0
        d=(q/qv,F(0),-q*q*qu/(4*qv))
        root=V.scale(V.sqrt(V.ki(d)),rec['radical_branch'])
        raw.append([V.add(V.scale(V.T,-q*a/2),V.scale(root,b)) for a,b in zip(u,v,strict=True)])
        norms.append(q)
    for r in coupled['new_caps']:
        raw.append([tuple(map(int,x))for x in r['normal_intervals']]);norms.append(F(r['norm']))
    assert len(raw)==25
    upper=[];signs=[];projection_norms=[]
    for i,(r,q)in enumerate(zip(raw,norms,strict=True)):
        s=1 if r[-1][0]>0 else -1
        assert r[-1][1]<0 if s<0 else r[-1][0]>0
        rr=r if s>0 else flip(r);upper.append(rr);signs.append(s)
        assert ipair(M,rr,ivec(h))[1]<0
        oldq=V.add(V.rat(q),V.neg(V.scale(V.mul(rr[-1],rr[-1]),2)))
        if i in (4,5):
            assert q==2 and rr[-1]==V.rat(1)
            oldq=(0,0) # exact inherited norm and exact last coordinate
        else:assert oldq[1]<0
        projection_norms.append(enc(oldq))

    # Reflection of the added coordinate is an exact form isometry and
    # permutes every rational parent wall. Thus all signed parent relations
    # inherit their exact angles from the 835/839 certificates.
    lookup={tuple(r):i for i,r in enumerate(R)}
    perm=[lookup[tuple(r[:-1]+[-r[-1]])]for r in R]
    assert all(perm[perm[i]]==i for i in range(812)) and h[-1]==0
    assert M[-1][-1]==2 and all(M[-1][i]==M[i][-1]==0 for i in range(22))
    assert signs[17]==signs[23]==signs[24]==1
    assert 829 in coupled['new_caps'][0]['right_neighbors']
    assert 829 in coupled['new_caps'][1]['right_neighbors']
    assert F(coupled['new_mutual_pairing'])==0
    right={(17,23),(17,24),(23,24)}
    signed=[];types=[];heights=[]
    for i,r in enumerate(upper):
        signed.extend([r,flip(r)]);types.extend([i,i]);heights.extend([1,-1])
    counts=Counter();pairs=[]
    for i,r in enumerate(signed):
        for j in range(i):
            ti,tj=types[i],types[j];p=ipair(M,r,signed[j])
            if heights[i]==heights[j] and tuple(sorted((ti,tj)))in right:
                assert p[0]<=0<=p[1]
                # Both originals are upper; the two right contacts with
                # cap17 and the new mutual right contact were proved by
                # exact projection/frame identities in the predecessor.
                kind='right_by_inherited_identity'
            elif ti==tj and ti in (4,5):
                assert heights[i]!=heights[j]
                # A null old projection and heights ±1 give pairing -2.
                assert p[0]<=-2*V.S<=p[1];kind='parallel_by_exact_null_projection'
            else:
                V.dotted(p,norms[ti],norms[tj]);kind='dotted'
            counts[kind]+=1
            pairs.append({'i':i,'j':j,'kind':kind,'pairing':enc(p)})
    assert counts=={'right_by_inherited_identity':6,'parallel_by_exact_null_projection':2,'dotted':1217}
    # All old 839 facets are retained in this symmetric set.
    oldmap=list(range(812))+[812+2*i+(signs[i]<0) for i in range(23)]
    oldmap += [812+2*4,812+2*5,812+2*23,812+2*24]
    assert len(oldmap)==839 and len(set(oldmap))==839
    assert signs[23:]==[1,1]
    # The independent integer escape is checked against every signed cap.
    x=list(map(F,coupled['small_integer_strict_escape']['vector']))
    assert V.pair(M,x,x)==1586 and V.pair(M,x,h)==-814
    oldpairs=[V.pair(M,x,r)for r in R];assert max(oldpairs)<0
    newpairs=[ipair(M,ivec(x),r)for r in signed];assert max(p[1]for p in newpairs)<0
    delta=list(map(F,coupled['untouched_old_cusp']['ray']))
    assert len(delta)==23 and delta[-1]==0 and V.pair(M,delta,delta)==0
    dp=[ipair(M,ivec(delta),r)for r in signed];assert max(p[1]for p in dp)<0
    links=[]
    for k,link in zip((4,5),cusps['bounded_cusp_links'],strict=True):
        ray=list(map(F,link['primitive_integer_null_ray']));assert ray[-1]==0
        pp=[ipair(M,ivec(ray),r)for r in signed]
        for j,p in enumerate(pp):
            if types[j]==k:assert p[0]<=0<=p[1] # inherited exact old-cusp incidence
            else:assert p[1]<0
        links.append({'old_cusp_type':k,'ray':list(map(str,ray)),
                      'new_active_caps':[812+2*k,813+2*k],
                      'bounded_link_inherited':True,'all_signed_cap_pairings':[enc(p)for p in pp]})
    sources=(V.OUT,V.REC,COUPLED,CUSPS,'scripts/'+TAG+'.py',
             'scripts/'+V.PRE+'_verify.py','scripts/fresh_h22_pi7_coupled_caps_20261013.py')
    out={'status':'PASS_EXACT_862_FACET_SYMMETRIC_ONE_SIDE_BRANCH',
         'sources':{p:V.sha(p)for p in sources},'base_facets':812,'cap_types':25,
         'signed_caps':50,'facets':862,'interval_denominator':str(V.S),
         'original_cap_height_signs':signs,'old_projection_norm_intervals':projection_norms,
         'parent_reflection_permutation':perm,'old839_inclusion':oldmap,
         'cap_pair_counts':dict(counts),'all_signed_cap_pairs':pairs,
         'norms_by_type':list(map(str,norms)),
         'strict_escape':{'vector':list(map(str,x)),'norm':1586,'center_pair':-814,
                          'maximum_parent_pair':str(max(oldpairs)),
                          'all_signed_cap_pairings':[enc(p)for p in newpairs]},
         'untouched_cusp':{'ray':list(map(str,delta)),
                           'all_signed_cap_pairings':[enc(p)for p in dp]},
         'preserved_bounded_links':links,
         'scope':'All839 old facets are retained. Exact reflection inherits every parent relation; complete interval checks and inherited identities prove all1225 signed-cap relations. A strict center gives862 actual facets. The two old bounded links and compact quadrilateral survive, and reflection gives its second copy. The exact norm1586 escape and untouched old cusp prove infinite volume. This is not a finite-volume or novelty claim.'}
    if verify:assert out==V.read(OUT)
    else:(V.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],dict(counts))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true')
    run(p.parse_args().verify)
