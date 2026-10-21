#!/usr/bin/env python3
"""Exact 860-facet replacement using pi/6 and rational square roots only."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
import sympy as sp
import fresh_h22_incidence_prism_20261013_verify as V
from fresh_h22_pi7_coupled_caps_20261013 import ipair,ivec,iadd,iscale,imul,div,pd

TAG='fresh_h22_multiquadratic_20261013';OUT='results/'+TAG+'.json'
COUPLED='results/fresh_h22_pi7_coupled_caps_20261013.json'
SELECT=[i for i in range(23)if i!=12]
def enc(x):return list(map(str,x))
def flip(r):return r[:-1]+[V.neg(r[-1])]

def run(verify=False):
    old=V.read(V.OUT);c=V.read(COUPLED)
    for obj in(old,c):
        for p,s in obj['sources'].items():assert V.sha(p)==s
    M,R,h=V.setup();recipes=V.read(V.REC)['recipes'];rn=[V.pair(M,r,r)for r in R]
    T=V.sqrt(V.rat(3));norms={};rad={};original={};upper={};signs={};records=[]
    lookup={tuple(r):i for i,r in enumerate(R)}
    perm=[lookup[tuple(r[:-1]+[-r[-1]])]for r in R]
    assert all(perm[perm[j]]==j for j in range(812)) and h[-1]==0
    assert M[-1][-1]==2 and all(M[i][-1]==M[-1][i]==0 for i in range(22))
    for k,rec in enumerate(recipes):
        u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm'])
        qu,qv=V.pair(M,u,u),V.pair(M,v,v)
        assert qv>0 and qu<=0 and V.pair(M,u,v)==0
        A=q/qv;C=-q*q*qu/(4*qv);d=A+3*C;assert d>0
        assert 3*q*q*qu/4+d*qv==q
        rad[k]=d;norms[k]=q
        normal=[V.add(V.scale(T,-q*a/2),V.scale(V.sqrt(V.rat(d)),rec['radical_branch']*b))for a,b in zip(u,v,strict=True)]
        original[k]=normal;s=1 if normal[-1][0]>0 else -1
        assert normal[-1][1]<0 if s<0 else normal[-1][0]>0
        signs[k]=s;upper[k]=normal if s>0 else flip(normal)
        # Exact identities for a free positive T, not reduction modulo its
        # former cubic polynomial. Thus labels7 become6 when T=sqrt3.
        for fc in rec['finite_contacts']:
            j,m=fc['wall'],fc['label'];L=-q*V.pair(M,u,R[j])/2
            B=rec['radical_branch']*V.pair(M,v,R[j])
            if m==2:assert L==B==0
            elif m==7:assert B==0 and L<0 and L*L==q*rn[j]/4
            else:assert m==4 and L==0 and C==0 and B<0 and B*B*A==q*rn[j]/2
        if k not in SELECT:continue
        assert ipair(M,upper[k],ivec(h))[1]<0
        ff={x['wall']:x['label']for x in rec['finite_contacts']};ps=[];labels={}
        for j,w in enumerate(R):
            original_j=j if s>0 else perm[j];p=ipair(M,upper[k],ivec(w))
            if original_j in ff:
                m=6 if ff[original_j]==7 else ff[original_j];labels[str(j)]=m
                if m==2:assert p[0]<=0<=p[1]
                else:assert p[1]<0
            else:V.dotted(p,q,rn[j])
            ps.append(enc(p))
        records.append({'type':k,'norm':str(q),'radicand':str(d),'original_height_sign':s,
                        'finite_parent_labels':labels,'all_parent_pairings':ps})

    # The omitted type12 has cos²=48/49 with parent40. If it equalled
    # cos²(pi/m), the rational number 2cos(2pi/m)=94/49 would be an
    # algebraic integer, impossible. This excludes every integer label.
    rr=recipes[12];u12=list(map(F,rr['u']));v12=list(map(F,rr['v']))
    assert norms[12]==rn[40]==2 and V.pair(M,u12,R[40])==F(8,7)
    assert V.pair(M,v12,R[40])==0
    assert F(3*64,49)/4==F(48,49) and F(4*48,49)-2==F(94,49)

    # A rational Lorentz-plane involution supplies the normalized future
    # vector without introducing a nested square-root field.
    r17=recipes[17];a=list(map(F,r17['u']));b=list(map(F,r17['v']))
    assert(V.pair(M,a,a),V.pair(M,b,b),V.pair(M,a,b))==(F(-9,2),2,0)
    assert r17['radical_branch']==1 and rad[17]==F(31,4) and norms[17]==2
    assert signs[17]==1
    S=c['old_incident'][:-1];assert c['old_pi7_cap']==829 and len(S)==21
    source_pivots=pd([[V.pair(M,R[i],R[j])for j in S]for i in S])
    assert all(V.pair(M,a,R[j])==V.pair(M,b,R[j])==0 for j in S)
    Z=iadd(imul(ivec(a),V.scale(V.sqrt(V.rat(31)),F(1,3))),
           imul(ivec(b),V.scale(T,F(-3,2))))
    orientation=1
    if ipair(M,Z,ivec(h))[0]>0:Z=iscale(Z,-1);orientation=-1
    assert ipair(M,Z,ivec(h))[1]<0
    assert F(31,9)*F(-9,2)+F(27,4)*2==-2
    # Coefficient of sqrt(93) in Z.cap17, up to its orientation sign.
    assert F(-1,3)*F(-9,2)+F(-3,4)*2==0
    du,dv=[list(map(F,r))for r in c['rational_duals']]
    assert(V.pair(M,du,du),V.pair(M,dv,dv),V.pair(M,du,dv))==(2,1,1)
    assert all(V.pair(M,du,R[j])==int(j==27)and V.pair(M,dv,R[j])==int(j==9)for j in S)
    # du,dv belong to span(S); hence they are orthogonal to a,b,Z.
    assert all(V.pair(M,z,t)==0 for z in(du,dv)for t in(a,b))
    rA=iscale(iadd(iscale(Z,209),ivec([-241*x for x in du])),F(1,120))
    rB=imul(iadd(iscale(Z,241),ivec([-418*x for x in dv])),div(V.rat(1),V.sqrt(V.rat(29281))))
    assert F(-2*209**2+2*241**2,120**2)==2
    assert F(-2*241**2+418**2,29281)==2 and -2*209*241+241*418==0
    for k,r,release in((23,rA,27),(24,rB,9)):
        assert r[-1][0]>0 and ipair(M,r,ivec(h))[1]<0
        upper[k]=r;norms[k]=F(2);signs[k]=1;ps=[];ff={str(j):2 for j in S if j!=release}
        for j,w in enumerate(R):
            p=ipair(M,r,ivec(w))
            if str(j)in ff:assert p[0]<=0<=p[1]
            else:V.dotted(p,2,rn[j])
            ps.append(enc(p))
        records.append({'type':k,'norm':'2','original_height_sign':1,'finite_parent_labels':ff,'all_parent_pairings':ps})

    ids=SELECT+[23,24];right={(17,23),(17,24),(23,24)}
    signed=[];types=[];heights=[]
    for k in ids:
        signed.extend((upper[k],flip(upper[k])));types.extend((k,k));heights.extend((1,-1))
    counts=Counter();pairdata=[]
    for i,r in enumerate(signed):
        for j in range(i):
            ti,tj=types[i],types[j];p=ipair(M,r,signed[j])
            if heights[i]==heights[j]and tuple(sorted((ti,tj)))in right:
                assert p[0]<=0<=p[1];kind='right_by_exact_frame_identity'
            elif ti==tj and ti in(4,5):
                assert upper[ti][-1]==V.rat(1) and norms[ti]==2
                assert p[0]<=-2*V.S<=p[1];kind='parallel_by_exact_null_projection'
            else:V.dotted(p,norms[ti],norms[tj]);kind='dotted'
            counts[kind]+=1;pairdata.append({'i':i,'j':j,'kind':kind,'pairing':enc(p)})
    assert counts=={'dotted':1120,'right_by_exact_frame_identity':6,'parallel_by_exact_null_projection':2}

    # The same exact Lorentz-plane construction still bounds an actual
    # compact quadrilateral; all its vertices are tested against every wall.
    cap17=812+2*ids.index(17);capA=812+2*ids.index(23);capB=812+2*ids.index(24)
    I=S+[cap17];edges=[[i for i in I if i!=j]for j in(27,9)]
    common=[i for i in I if i not in(27,9)];assert len(common)==20
    triples=[(1,0,0),(241,-209,0),(50369,-29281,-28800),(209,0,-241)]
    active=[I,edges[0]+[capA],common+[capA,capB],edges[1]+[capB]]
    walls=[ivec(r)for r in R]+signed;vertices=[]
    for co,inc in zip(triples,active,strict=True):
        aa,bb,cc=co
        point=iadd(iscale(Z,aa),ivec([bb*x+cc*y for x,y in zip(du,dv,strict=True)]))
        q=-2*aa*aa+2*bb*bb+cc*cc+2*bb*cc;assert q<0
        future=ipair(M,point,ivec(h));assert future[1]<0
        assert len(inc)==22
        if capA in inc:assert -418*aa-241*(2*bb+cc)==0
        if capB in inc:assert -482*aa-418*(bb+cc)==0
        ps=[]
        for j,w in enumerate(walls):
            val=ipair(M,point,w)
            if j in inc:assert val[0]<=0<=val[1]
            else:assert val[1]<0,(co,j,val)
            ps.append(enc(val))
        vertices.append({'frame_coefficients':list(co),'norm':q,'incident':inc,
                         'future_center_pairing':enc(future),'all_wall_pairings':ps})
    x=list(map(F,c['small_integer_strict_escape']['vector']));assert V.pair(M,x,x)==1586 and V.pair(M,x,h)==-814
    assert max(V.pair(M,x,r)for r in R)<0
    escape=[ipair(M,ivec(x),r)for r in signed];assert max(p[1]for p in escape)<0
    delta=list(map(F,c['untouched_old_cusp']['ray']));assert V.pair(M,delta,delta)==0 and delta[-1]==0
    cusp=[ipair(M,ivec(delta),r)for r in signed];assert max(p[1]for p in cusp)<0

    # A concrete irrational length-two cyclic product generates Q(sqrt3).
    rec=recipes[1];uu=list(map(F,rec['u']));vv=list(map(F,rec['v']));assert rad[1]==4
    for j,w in enumerate(R):
        L=-norms[1]*V.pair(M,uu,w)/2;B=rec['radical_branch']*V.pair(M,vv,w)
        if L*B:
            ca=4*(3*L*L+4*B*B)/(norms[1]*rn[j]);cb=16*L*B/(norms[1]*rn[j]);break
    assert cb!=0
    oldshort=[i for i in range(210)if rn[i]==2]
    assert len(oldshort)==42 and sum(V.pair(M,R[i],R[j])for i in oldshort for j in oldshort)==-126
    generators=sorted(set([F(3),F(31),F(29281)]+[rad[k]for k in SELECT]))
    sqclasses=[]
    for d in generators:
        f=sp.factorint(d.numerator*d.denominator);sf=1
        for prime,e in f.items():
            if e%2:sf*=int(prime)
        if sf!=1 and sf not in sqclasses:sqclasses.append(sf)
    sources=(V.OUT,V.REC,COUPLED,'scripts/'+TAG+'.py','scripts/'+V.PRE+'_verify.py')
    out={'status':'PASS_EXACT_860_FACET_MULTIQUADRATIC_PI6_REPLACEMENT',
         'sources':{p:V.sha(p)for p in sources},'base_facets':812,'selected_original_types':SELECT,
         'omitted_type':{'index':12,'parent':40,'forbidden_cosine_square':'48/49','rational_noninteger_2cos2angle':'94/49'},
         'signed_caps':48,'facets':860,'interval_denominator':str(V.S),'cap_records':records,
         'cap_pair_counts':dict(counts),'all_signed_cap_pairs':pairdata,
         'parent_reflection_permutation':perm,'normalized_frame_orientation':orientation,
         'old21_positive_pivots':enc(source_pivots),'common20_walls':common,
         'compact_H2_vertices':vertices,'intrinsic_quadrilateral_angles':['pi/4','pi/2','pi/2','pi/2'],
         'field':{'description':'All normal coordinates belong to a compositum of square roots of positive rational numbers; the displayed rational form is preserved.',
                  'rational_radicands':list(map(str,generators)),'squarefree_generators':sorted(sqclasses),
                  'quadratic_irrational_cycle':{'original_type':1,'parent_wall':j,'rational_part':str(ca),'sqrt3_coefficient':str(cb)},
                  'old_rational_indefinite_ones_value':-126},
         'strict_escape':{'vector':list(map(str,x)),'norm':1586,'center_pair':-814,'all_cap_pairings':[enc(p)for p in escape]},
         'untouched_cusp':{'ray':list(map(str,delta)),'all_cap_pairings':[enc(p)for p in cusp]},
         'scope':'A new construction replacing pi7 by pi6 and omitting one incompatible type. All860 facets pass exact norm, angle, center and compatibility checks, with four verified finite vertices bounding an actual compact H2 quadrilateral and its reflected copy. Their coordinates lie in an explicitly multiquadratic extension ofQ. An irrational quadratic cycle and retained rational indefinite block preclude quasi-arithmetic finite-volume completion; the cubic/non-pseudo obstruction from pi7 is absent. Exact strict escape proves this chamber still has infinite volume. No complete vertex/cusp census or finite-volume existence is claimed.'}
    if verify:assert out==V.read(OUT)
    else:(V.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],dict(counts),'squareclasses',sqclasses)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
