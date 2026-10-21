#!/usr/bin/env python3
"""Exact135-wall extension: old bank recipe44 cuts the new weak escape.

All contact identities are defining algebraic equations. Strict signs use
outward rational intervals. No complete ambient cone census is claimed.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from collections import Counter
from fresh_h20_coupled_roles_20261014_field import Field,parse,B
from fresh_h20_s61_multicap_20261014_exact import I,FILES,BITS,fdot,pairing

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_complementary_cap_20261016'
BASE='results/fresh_h20_s61_multicap_20261014_exact.json'
CAP='results/fresh_h20_dotted_arc_cap_20261016.json'
CAPREVIEW='results/fresh_h20_dotted_arc_cap_20261016_review.json'
WEAK='results/fresh_h20_asymmetric_arc_escape_gate_20261016.json'
BANK='results/fresh_h20_nonright_s61_20261014_scout.json'
BANKGATE='results/fresh_h20_nonright_s61_20261014_bank_gate.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def saved(x):
    den=1<<192
    return [str(F((x.l*den).__floor__(),den)),str(F((x.h*den).__ceil__(),den))]
def outs(v):return [saved(x)for x in v]
def contains(x,t=0):return x.l<=t<=x.h

def run():
    old,bridge,edge,polar,cur=[json.loads((ROOT/p).read_text())for p in FILES]
    base,cap,capreview,weak,bank,bankgate=[json.loads((ROOT/p).read_text())for p in [BASE,CAP,CAPREVIEW,WEAK,BANK,BANKGATE]]
    for doc in [old,bridge,edge,polar,cur,base,cap,capreview,weak,bank,bankgate]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert capreview['sources'][CAP]==sha(CAP)
    prior=bank['proposals'][44]
    assert prior['profile']==16 and prior['anchors']==[95,0] and prior['anchor_labels']==[2,6]
    assert prior['contacts']==[float(i==4)for i in range(18)]
    A=Field(old,bridge,edge);z=A.z
    def bound(x,k=3):return I(*A.bounds(x,k,BITS))
    M=[tuple(parse(x,1)for x in v)for v in old['normals']]
    N=[A.vlift(v,1)for v in M]+[A.vlift(tuple(parse(x,2)for x in bridge['new_normal']),2)]
    N +=[tuple(parse(x,3)for x in v)for v in edge['new_normals']]
    N +=[tuple(A.c(x)for x in v)for v in polar['new_normals']]
    norms=[A.c(q)for q in old['norms']]+[A.c(2)]+[parse(edge['new_normal_norm'],3)]*2+[A.c(q)for q in polar['new_norms']]
    N.append(A.vlift(tuple(parse(x,1)for x in cur['cap_R']['normal']),1))
    norms.append(A.lift(parse(cur['cap_R']['norm'],1),1));assert len(N)==101
    nd=[(bound(q)/2).sqrt()for q in norms]
    S=base['support'];J=[[F(q)for q in row]for row in base['support_inverse']]
    G=[[F(q)for q in row]for row in base['support_gram']]
    assert B.inv(G)==J and len(S)==18
    j=S.index(5);assert J[j][j]==F(30,11)
    omega=A.lin([M[i]for i in S],[A.c(q,1)for q in J[j]],1)
    assert A.dot(omega,omega,1)==A.c(F(30,11),1)
    assert A.dot(omega,M[0],1)==A.c(F(-6,11),1)
    assert A.dot(A.vlift(omega,1),N[95])==z[3]
    frame=[tuple(parse(x,2)for x in base['frame'][key])for key in ['X','r95','W']]
    fsq=[parse(base['frame']['A'],2),A.c(2,2),parse(base['frame']['W_norm'],2)]
    den=[bound(q,2).sqrt()for q in fsq];frame3=[A.vlift(v,2)for v in frame]
    assert all(A.dot(v,A.vlift(omega,1,2),2)==z[2]for v in frame)
    assert [[A.dot(v,w,2)for w in frame]for v in frame]==[
        [A.scale(fsq[0],-1,2),z[2],z[2]],[z[2],fsq[1],z[2]],[z[2],z[2],fsq[2]]]
    hb=[tuple(parse(x,1)for x in cur['hex_basis'][key])for key in ['C','U','V']]
    hq=[parse(x,1)for x in cur['hex_basis']['squared_lengths']]
    def old_pairs(v):
        v3=A.vlift(v,1)
        ans=[bound(A.dot(v3,n))/dd for n,dd in zip(N,nd)]
        for co in cur['hexagon']['frame_coefficients']:
            ans.append(sum((bound(A.lift(parse(c,0),0))*bound(A.dot(v3,A.vlift(w,1)))/bound(A.lift(q,1)).sqrt()
                            for c,w,q in zip(co,hb,hq)),I(0)))
        assert len(ans)==107
        return ans
    wp=old_pairs(omega)
    P=[[I(*x)for x in row]for row in base['frame']['old107_unit_normal_pair_bounds']]
    coeff=[[I(*x)for x in row]for row in base['chain']['frame_coefficient_bounds']]
    P +=[[-v[0],v[1],v[2]]for v in coeff];wp +=[I(0)]*26
    for k,i in enumerate(S):
        assert contains(wp[i],int(i==5));wp[i]=I(int(i==5))
        assert all(x.l==x.h==0 for x in P[i])
    assert contains(wp[0],F(-6,11));wp[0]=I(F(-6,11))
    P[95]=[I(0),I(2).sqrt(),I(0)];wp[95]=I(0)
    # Let p=P0_t, q=P0_f. y_e=0 fixes contact95. Solve p*t+q*z=h,
    # then t^2-z^2=8/11. Positive z branch is the literal old recipe44.
    p,q=P[0][0],P[0][2];h=I(F(-6,11))-I(3).sqrt()
    assert p.h<0 and q.h<0
    alpha=h/p;beta=-q/p
    qa=1-beta*beta;qb=-2*alpha*beta;qc=I(F(8,11))-alpha*alpha
    disc=qb*qb-4*qa*qc;assert qa.l>0 and disc.l>0
    zz=(-qb+disc.sqrt())/(2*qa);tt=alpha+beta*zz
    yy=[tt,I(0),zz];assert tt.l>0 and zz.l>0
    assert contains(fdot(yy,yy),F(-8,11))
    assert contains(pairing(yy,P[0])-wp[0]+I(3).sqrt())
    box=[(F('0.9491'),F('0.9493')),(F('-0.0001'),F('0.0001')),(F('0.4167'),F('0.4169'))]
    assert all(lo<x.l<=x.h<hi for x,(lo,hi)in zip(yy,box))
    alt=(-qb-disc.sqrt())/(2*qa);assert alt.h<0
    # Independent source algebra for finite-corner exclusions below.
    assert A.dot(M[0],M[1],1)==A.c(-1,1)
    assert A.sign(A.dot(M[0],M[87],1)+A.c(2,1),1)<0
    vals=[pairing(yy,v)-w for v,w in zip(P,wp)]
    y133=[I(*q)for q in cap['definition']['frame_bounds']]
    c133=I(*cap['definition']['dual_coefficient_bounds'])
    assert J[j][S.index(31)]==0
    mutual=fdot(yy,y133);assert mutual.h<-2
    vals.append(mutual)
    finite={i:2 for i in S};finite[5]=3;finite.update({0:6,95:2})
    counts=Counter();pairs=[]
    for i,g in enumerate(vals):
        if i in finite:
            m=finite[i];expected=I(0)if m==2 else -I(1)if m==3 else -I(3).sqrt()
            assert contains(g-expected);relation='right'if m==2 else f'pi/{m}'
        else:assert g.h<-2,(i,saved(g));relation='dotted'
        counts[relation]+=1;pairs.append(dict(wall=i,relation=relation,pair_bounds=saved(g)))
    assert counts=={'right':18,'pi/3':1,'pi/6':1,'dotted':114}
    print('PASS all134 new contact gates',dict(counts),flush=True)
    # Inherited strict center, actual-facet projections, oldfinite preservation.
    H=tuple(parse(x,1)for x in cur['strict_center']['point']);qH=A.dot(H,H,1);assert A.sign(qH,1)<0
    hd=bound(A.scale(qH,-1,1),1).sqrt()
    hp=[bound(A.dot(v,A.vlift(H,1)))/(d*hd)for v,d in zip(frame3,den)]
    wh=bound(A.dot(omega,H,1),1)/hd
    eps=F(base['strict_center']['epsilon']);center=[I(1),I(F(-1,2)),I(0)]
    newcenter=fdot(yy,center)+eps*(pairing(yy,hp)-wh);assert newcenter.h<0
    centerpairs=[I(*v)for v in cap['common_center']['all134_pair_bounds']]+[newcenter]
    cq=I(*cap['common_center']['norm_bounds']);assert cq.h<0 and all(g.h<0 for g in centerpairs)
    facetnorms=[cq-g*g/2 for g in centerpairs];assert all(v.h<0 for v in facetnorms)
    # Cut the E0 weakescape and its previously strict epsilon version.
    e0=fdot(yy,[I(1),I(F(-1,100)),I(0)])+J[j][j];assert e0.l>0
    oldeb=list(map(F,cap['strict_escape']['dual_coefficients']))
    oldwt=sum(J[j][k]*oldeb[k]for k in range(18))
    eeps=fdot(yy,[I(1),I(F(-1,100)),I(0)])-oldwt;assert eeps.l>0
    # A simple remaining strict escape discovered in the complete135 halfspaces.
    ez=[F(1),F(-1,100),F(-49,50)];eb=[F(-1,100000)]*18;eb[S.index(17)]=F(-1)
    wt=[sum(J[i][k]*eb[k]for k in range(18))for i in range(18)]
    V=A.lin([M[i]for i in S],[A.c(x,1)for x in wt],1)
    assert all(A.dot(A.vlift(V,1,2),v,2)==z[2]for v in frame)
    ve=old_pairs(V)+[I(0)]*26
    for k,i in enumerate(S):assert contains(ve[i],eb[k]);ve[i]=I(eb[k])
    eg=[pairing(ez,p0)+v for p0,v in zip(P,ve)]
    eg +=[fdot(y133,[I(x)for x in ez])-c133*wt[S.index(31)],fdot(yy,[I(x)for x in ez])-wt[j]]
    assert len(eg)==135 and all(g.h<0 for g in eg)
    qe=-ez[0]*ez[0]+ez[1]*ez[1]+ez[2]*ez[2]+sum(eb[i]*wt[i]for i in range(18));assert qe>0
    eh=pairing(ez,hp)+bound(A.dot(V,H,1),1)/hd
    ec=fdot([I(x)for x in ez],center)+eps*eh;assert ec.h<0
    cusps=[]
    for raw in [old['completed_cusp']['point']]+[r['point']for r in old['preserved_old_bounded_cusps']]:
        nu=tuple(A.lift(parse(x,0),0)for x in raw)
        fp=[bound(A.dot(v,nu))/d for v,d in zip(frame3,den)]
        g=pairing(yy,fp)-bound(A.dot(A.vlift(omega,1),nu))
        assert g.h<0;cusps.append(saved(g))
    print('PASS coherence,135 facets, escapes and3 preserved cusps',flush=True)
    paths=[f'scripts/{TAG}.py',BASE,CAP,CAPREVIEW,WEAK,BANK,BANKGATE]+FILES+[
        'scripts/fresh_h20_coupled_roles_20261014_field.py','scripts/fresh_h20_s61_multicap_20261014_exact.py',
        'results/acute_extension_preservation_20261014.md']
    return dict(status='PASS_EXACT_135_WALL_COMPLEMENTARY_E0_CUTTER',sources={p:sha(p)for p in paths},
        interval_bits=BITS,retained_walls=134,new_wall=134,actual_facets=135,
        provenance=dict(old_numerical_proposal=44,old_profile=16,prior_exact_recipe=BANKGATE,
            genuinely_new_verification='All134 Coxeter gates, compatibility with newarc cap133, actual135facets, E0cut and remainingescape; oldrecipe was not previously promoted as a full Coxeter wall.'),
        definition=dict(support=S,contact_wall=5,contact_magnitude=1,dual_norm='30/11',
            anchors=[95,0],anchor_labels=[2,6],frame='(x,e,f) frozenS61 unitLorentzframe',
            formula='p=P[0,0];q=P[0,2];h=-6/11-sqrt3;alpha=h/p;beta=-q/p; qa=1-beta^2;qb=-2alpha beta;qc=8/11-alpha^2;z=(-qb+sqrt(qb^2-4qa qc))/(2qa);y=(alpha+beta*z,0,z);r=y-omega5',
            exact_residual_norm='-8/11',exact_cap_norm=2,quadratic_coefficients=outs([qa,qb,qc]),
            discriminant_bounds=saved(disc),frame_bounds=outs(yy),isolating_box=[[str(a),str(b)]for a,b in box]),
        pair_counts=dict(counts),pairs=pairs,finite_neighbors=sorted(finite),mutual_arc_cap_pair_bounds=saved(mutual),
        common_center=dict(inherited_recipe=cap['common_center']['inherited_recipe'],norm_bounds=saved(cq),
            all135_pair_bounds=outs(centerpairs),all135_facet_projection_norm_bounds=outs(facetnorms)),
        complete_new_finite_vertices=0,known_total_finite_vertices_at_least=170,
        finite_census_proof='Exactly20 oldfinite neighbors:S plus0,95. A20-root finite support uses19 ofthese. Dotted0/87 forces omission0or87. Omit0 retains S+r, whose Schurvalue is-8/11. Omit87 retains the affine triple(r,0,1) with bonds(pi6,pi3,right). Hence nonePD. Cap133 is dotted to r and cannotoccur.',
        finite_obstructions=dict(dotted_pair=[0,87],support_S_plus_cap_Schur='-8/11',
            affine_triple=[134,0,1],affine_Gram=[['2','-sqrt3','0'],['-sqrt3','2','-1'],['0','-1','2']],
            positive_kernel=['sqrt3','2','1']),
        weak_escape_cut=dict(source=WEAK,pair_bounds=saved(e0)),previous_strict_escape_cut=saved(eeps),
        preserved_three_cusp_pair_bounds=cusps,
        strict_escape=dict(frame_coefficients=list(map(str,ez)),support=S,dual_coefficients=list(map(str,eb)),
            norm=str(qe),all135_pair_bounds=outs(eg),strict_center_pair_bounds=saved(ec)),
        scope='Exact coherent135-wall chamber and allactualfacets; retains alloldfinitevertices andthreepropercusps. Complete zero newfinitevertex census only. Atleast170finitevertices, not an ambientcompletecensus. Strict spacelike escape proves infinitevolume. No finite-volume construction or allcap-family exhaustion claimed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['actual_facets'],flush=True)
