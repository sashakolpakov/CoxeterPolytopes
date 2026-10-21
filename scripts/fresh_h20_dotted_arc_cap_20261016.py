#!/usr/bin/env python3
"""One exact escape-cutting cap with a three-wall finite rotational arc.

The defining rotation parameter is inherited by its isolating interval.
Norm and prescribed contacts are algebraic identities; strict signs use
outward rational intervals. The former full ambient vertex census is partial.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from collections import Counter
from fresh_h20_coupled_roles_20261014_field import Field,parse,B
from fresh_h20_s61_multicap_20261014_exact import I,FILES,BITS,fdot,pairing

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_dotted_arc_cap_20261016'
BASE='results/fresh_h20_s61_multicap_20261014_exact.json'
BASE_REVIEW='results/fresh_h20_s61_multicap_20261014_exact_review.json'
BANK='results/fresh_h20_nonright_s61_20261014_bank_gate.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def saved(x):
    den=1<<192
    return [str(F((x.l*den).__floor__(),den)),str(F((x.h*den).__ceil__(),den))]
def outs(v):return [saved(x)for x in v]
def contains(x,t=0):return x.l<=t<=x.h

def run():
    old,bridge,edge,polar,cur=[json.loads((ROOT/p).read_text())for p in FILES]
    base=json.loads((ROOT/BASE).read_text());review=json.loads((ROOT/BASE_REVIEW).read_text())
    bank=json.loads((ROOT/BANK).read_text())
    for doc in [old,bridge,edge,polar,cur,base,review,bank]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
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
    j=S.index(31);assert all(J[j][k]==(F(1,2)if j==k else 0)for k in range(18))
    assert A.dot(M[31],M[31],1)==A.c(2,1)
    omega=tuple(A.scale(x,F(1,2),1)for x in M[31])
    assert A.dot(omega,omega,1)==A.c(F(1,2),1)
    frame=[tuple(parse(x,2)for x in base['frame'][key])for key in ['X','r95','W']]
    fsq=[parse(base['frame']['A'],2),A.c(2,2),parse(base['frame']['W_norm'],2)]
    den=[bound(q,2).sqrt()for q in fsq];frame3=[A.vlift(v,2)for v in frame]
    assert all(A.dot(v,A.vlift(omega,1,2),2)==z[2]for v in frame)
    assert [[A.dot(v,w,2)for w in frame]for v in frame]==[
        [A.scale(fsq[0],-1,2),z[2],z[2]], [z[2],fsq[1],z[2]], [z[2],z[2],fsq[2]]]
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
    assert base['chain']['walls']==list(range(107,133))
    P +=[[-v[0],v[1],v[2]]for v in coeff];wp +=[I(0)]*26
    for k,i in enumerate(S):
        assert contains(wp[i],int(i==31))
        wp[i]=I(int(i==31));assert all(x.l==x.h==0 for x in P[i])
    aa=I(*base['parameter']['isolating_interval']);C=aa*aa/(aa*aa+2)
    root2=I(2).sqrt();h=root2-1
    # n129 has spatial square a²+2. These formulas impose the central label3
    # and the two neighboring label4 contacts identically.
    n=coeff[129-107]
    yy=[(h*aa*aa+2*root2)/(2*aa),h*n[1]/2,h*n[2]/2]
    residual=-I(F(1,2))-2/(aa*aa)
    cc=(5+4/(aa*aa)).sqrt()
    assert contains(fdot(yy,yy)-residual) and contains(residual+cc*cc/2,2)
    assert yy[0].l>0 and residual.h<0 and cc.l>2
    vals=[pairing(yy,p)-cc*w for p,w in zip(P,wp)]
    finite={i:2 for i in S if i!=31};finite.update({128:4,129:3,130:4})
    counts=Counter();pairs=[]
    for i,g in enumerate(vals):
        if i in finite:
            m=finite[i];expected=I(0)if m==2 else -I(1)if m==3 else -root2
            assert contains(g-expected)
            relation='right'if m==2 else f'pi/{m}'
        else:assert g.h<-2,(i,saved(g));relation='dotted'
        counts[relation]+=1;pairs.append(dict(wall=i,relation=relation,pair_bounds=saved(g)))
    assert counts=={'right':17,'pi/4':2,'pi/3':1,'dotted':113}
    print('PASS all133 new pair gates',dict(counts),flush=True)
    # The same strict center makes every one of the134 acute walls actual.
    Hold=tuple(parse(x,1)for x in cur['strict_center']['point'])
    qH=A.dot(Hold,Hold,1);assert A.sign(qH,1)<0
    hd=bound(A.scale(qH,-1,1),1).sqrt()
    hf=[bound(A.dot(v,A.vlift(Hold,1)))/(d*hd)for v,d in zip(frame3,den)]
    hp=[x/hd for x in old_pairs(Hold)]+[pairing(v,hf)for v in coeff]
    wh=bound(A.dot(omega,Hold,1),1)/hd
    eps=F(base['strict_center']['epsilon']);center=[I(1),I(F(-1,2)),I(0)]
    centerpairs=[pairing(center,p)+eps*g for p,g in zip(P,hp)]
    newcenter=fdot(yy,center)+eps*(pairing(yy,hf)-cc*wh)
    centerpairs.append(newcenter);assert all(g.h<0 for g in centerpairs)
    cq=I(F(-3,4))+2*eps*pairing(center,hf)-eps*eps;assert cq.h<0
    facet_norms=[cq-g*g/2 for g in centerpairs];assert all(q.h<0 for q in facet_norms)
    # Test both previously saved strict escapes. The847-bank witness is cut;
    # the older witness is recorded separately rather than assumed cut.
    cuts=[]
    for name,raw in [('old133',base['strict_escape']['rational_point']),('bank847',bank['escape']['rational_point'])]:
        E=tuple(A.c(x)for x in raw)
        ep=[bound(A.dot(v,E))/d for v,d in zip(frame3,den)]
        ew=bound(A.dot(A.vlift(omega,1),E))
        g=pairing(yy,ep)-cc*ew;assert g.l>0 or g.h<0
        if name=='bank847':assert g.l>0
        cuts.append(dict(escape=name,pair_bounds=saved(g),cut=bool(g.l>0)))
    # Complete new finite-corner census: exactly two of the20 old finite
    # neighbors can supply19 roots without containing dotted128/130.
    vertices=[]
    for ab in [[128,129],[129,130]]:
        ov=next(v for v in base['face']['vertices']if v['walls']==ab)
        point=[I(*x)for x in ov['frame_point_bounds']]
        zy=fdot(point,yy);dualcoef=2*zy/cc
        q=fdot(point,point)+dualcoef*dualcoef/2;assert q.h<0 and dualcoef.h<0
        active=sorted(set(S)-{31}|set(ab))
        vp=[pairing(point,p)+dualcoef*w for p,w in zip(P,wp)]
        assert all(contains(vp[i])if i in active else vp[i].h<0 for i in range(133))
        assert contains(zy-cc*dualcoef/2)
        bonds=['-sqrt(2)'if i in [128,130]else'-1'for i in ab]
        vertices.append(dict(old_active=active,active=active+[133],frame_point_bounds=outs(point),
                             omega31_coefficient_bounds=saved(dualcoef),norm_bounds=saved(q),
                             old133_pair_bounds=outs(vp),future_time_coordinate=1,
                             spherical_last3_Gram=[['2','0',bonds[0]],['0','2',bonds[1]],[bonds[0],bonds[1],'2']],
                             cap_Schur_remainder='1/2'))
    # Exact remaining escape with rational coordinates in the orthogonal
    # Lorentz3 plus support-dual decomposition.
    ez=[F(1),F(-1,100),F(0)];eb=[F(-1,100000)]*18;eb[S.index(5)]=F(-1)
    wt=[sum(J[i][k]*eb[k]for k in range(18))for i in range(18)]
    V=A.lin([M[i]for i in S],[A.c(x,1)for x in wt],1)
    assert all(A.dot(A.vlift(V,1,2),v,2)==z[2]for v in frame)
    ve=old_pairs(V)+[I(0)]*26
    for k,i in enumerate(S):assert contains(ve[i],eb[k]);ve[i]=I(eb[k])
    eg=[pairing(ez,p)+v for p,v in zip(P,ve)]
    nr=fdot(yy,[I(x)for x in ez])-cc*wt[j];eg.append(nr)
    assert all(g.h<0 for g in eg)
    qe=-ez[0]*ez[0]+ez[1]*ez[1]+ez[2]*ez[2]+sum(eb[i]*wt[i]for i in range(18))
    assert qe>0
    eh=pairing(ez,hf)+bound(A.dot(V,Hold,1),1)/hd
    ec=fdot([I(x)for x in ez],center)+eps*eh;assert ec.h<0
    # Three inherited proper cusps: give their literal new signs; no assertion
    # of complete ideal-vertex census is made.
    cusps=[]
    for raw in [old['completed_cusp']['point']]+[r['point']for r in old['preserved_old_bounded_cusps']]:
        nu=tuple(A.lift(parse(x,0),0)for x in raw)
        fp=[bound(A.dot(v,nu))/d for v,d in zip(frame3,den)]
        g=pairing(yy,fp)-cc*bound(A.dot(A.vlift(omega,1),nu))
        assert g.h<0;cusps.append(saved(g))
    paths=[f'scripts/{TAG}.py',BASE,BASE_REVIEW,BANK]+FILES+[
        'scripts/fresh_h20_coupled_roles_20261014_field.py',
        'scripts/fresh_h20_s61_multicap_20261014_exact.py',
        'results/acute_extension_preservation_20261014.md']
    return dict(status='PASS_EXACT_134_WALL_DOTTED_ARC_ESCAPE_CUTTER',sources={p:sha(p)for p in paths},
        interval_bits=BITS,retained_walls=133,new_wall=133,actual_facets=134,
        definition=dict(parameter_dependency=BASE,parameter_interval=base['parameter']['isolating_interval'],
            frame='unit Lorentz frame (x,e,f) from frozen S61 package',support=S,isolated_support_wall=31,
            dual_weight='omega31=n31/2',central_chain_wall=129,
            formula='h=sqrt(2)-1; y=((h*a^2+2*sqrt(2))/(2*a), h*n129_e/2, h*n129_f/2); c=sqrt(5+4/a^2); r=y-c*omega31',
            exact_residual_norm='-1/2-2/a^2',exact_cap_norm=2,frame_bounds=outs(yy),
            dual_coefficient_bounds=saved(cc),rotation_cosine_bounds=saved(C)),
        pair_counts=dict(counts),pairs=pairs,finite_neighbors=sorted(finite),
        common_center=dict(inherited_recipe=base['strict_center']['recipe'],norm_bounds=saved(cq),
            all134_pair_bounds=outs(centerpairs),all134_facet_projection_norm_bounds=outs(facet_norms),
            actual_facet_proof='For unit-norm2 wall n_i, H_i=H-(H.n_i)n_i/2 is future timelike, orthogonal to n_i and strictly negative on every other wall because all cross-pairs are nonpositive.'),
        former_escapes_tested=cuts,complete_new_finite_vertices=vertices,complete_new_finite_vertex_count=2,
        finite_census_proof='Finite neighbors are17 retained support roots plus128/129/130. Every19-subset except omission128 or130 contains the dotted pair128/130. Each remaining support has cap Schur remainder1/2 and the displayed feasible future timelike point.',
        known_total_finite_vertices_at_least=170,old_known_finite_vertices_preserved=168,
        three_new_finite_edges=[dict(old_endpoint=S+[128,129],new_endpoint=vertices[0]['active']),
                               dict(new_endpoints=[v['active']for v in vertices]),
                               dict(old_endpoint=S+[129,130],new_endpoint=vertices[1]['active'])],
        preserved_three_cusp_pair_bounds=cusps,
        strict_escape=dict(frame_coefficients=list(map(str,ez)),dual_coefficients=list(map(str,eb)),support=S,
            norm=str(qe),all134_pair_bounds=outs(eg),strict_center_pair_bounds=saved(ec),
            conclusion='A strict future-oriented spacelike escape remains, so the134-wall chamber has infinite volume.'),
        scope='Exact one-wall coherent extension retaining all133 old walls. Cuts the stronger847-bank escape, adds exactly two finite vertices on the new facet and at least three finite edges. The older saved escape is separately tested. The inherited ambient finite-vertex count was only a lower bound; no complete ambient vertex/ray census or finite-volume/new-dimension claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    r=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert r==json.loads(p.read_text())
    else:p.write_text(json.dumps(r,indent=2)+'\n')
    print(r['status'],r['pair_counts'],'finite+',r['complete_new_finite_vertex_count'],flush=True)
