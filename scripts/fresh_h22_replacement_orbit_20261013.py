#!/usr/bin/env python3
"""Deleted-wall stabilizer orbit and its full-plane common-escape obstruction.

The complete old group is a hash-bound independently certified dependency.
All new geometry is checked with integers; native array arithmetic has explicit
pre-evaluation size bounds. No floating-point predicate is used.
"""
import argparse,json,math
from collections import Counter
import numpy as np
from fresh_h22_core_replacement_20261013 import P,RAW,CORE,BANK,WALK,OUT as SINGLE,read,sha,setup,pair,parent_dotted,positive_pivots
from fractions import Fraction as F

PRE='fresh_h22_replacement_orbit_20261013'
OUT='results/'+PRE+'.json'
WALK_AUDIT='results/fresh_h21_source_corner_walk_20261013_independent.json'

def maxabs(a):return int(np.max(np.abs(a)))
def pair_bound(Q,A,B):
    assert 22**2*maxabs(Q)*maxabs(A)*maxabs(B)<2**63
def seq(v):return list(map(int,v))

def run():
    single=read(SINGLE)
    assert single['status']=='PASS_EXACT_813_FACET_CHANGED_CORE_WITH_FINITE_VERTEX_AND_STRICT_ESCAPE'
    for path,digest in single['sources'].items():assert sha(path)==digest
    walk=read(WALK);audit=read(WALK_AUDIT)
    assert audit['producer_sha256']==sha(WALK) and audit['queue_closed']
    groupfile=walk['group_array'];assert sha(groupfile)==walk['group_sha256']
    full=np.load(P/groupfile)['permutations'];assert full.shape==(241920,210)
    Q,W,C,h,M,R=setup()
    Q=np.array(Q,dtype=np.int64);W=np.array(W,dtype=np.int64);C=np.array(C,dtype=np.int64)
    h=np.array(h,dtype=np.int64);M=np.array(M,dtype=np.int64);R=np.array(R,dtype=np.int64)
    assert 23**2*maxabs(M)*maxabs(R)**2<2**63
    G=R@M@R.T;norms=np.diag(G)
    assert Counter(seq(norms))=={2:644,4:168}
    ii,jj=np.tril_indices(812,-1);ps=G[ii,jj];den=norms[ii]*norms[jj]
    assert np.all(ps<=0)
    assert 4*maxabs(ps)**2<2**63
    sq=ps*ps
    assert np.all((sq==0)|(4*sq==den)|(2*sq==den)|(4*sq==3*den)|(sq>=den))
    h23=np.append(h,0);assert int(h23@M@h23)==-14 and np.all(R@(M@h23)<0)
    deleted=single['deleted_old_wall'];assert deleted==1
    ids=np.flatnonzero(full[:,deleted]==deleted);g=full[ids]
    assert len(ids)==5760 and int(ids[0])==0
    assert len({row.tobytes() for row in g})==5760
    points,lines=read(CORE)['incidence_parts']
    assert np.all(np.sort(g[:,points],axis=1)==np.array(sorted(points)))
    assert np.all(np.sort(g[:,lines],axis=1)==np.array(sorted(lines)))
    vp=np.array(read(CORE)['old_projections'][1],dtype=np.int64)
    assert np.array_equal(21*vp,5*W[points].sum(axis=0)+2*W[lines].sum(axis=0))
    assert np.array_equal(3*h,W[points].sum(axis=0)+W[lines].sum(axis=0))
    # Explicit parent-bank preservation: vprime is determined by its old
    # short-wall pairings, and the stabilizer preserves both incidence parts.
    # Transporting a binary recipe therefore just permutes its point mask.
    pp=np.full(210,-1,dtype=np.int64);pp[points]=np.arange(21)
    pg=pp[g[:,points]];assert np.all(pg>=0)
    recipes=[r for r in read(BANK)['candidates'] if r['twice_height_square']==4]
    masks=np.array([r['mask'] for r in recipes],dtype=np.int64)
    As=np.array([r['A'] for r in recipes],dtype=np.int64)
    images=np.zeros((len(g),301),dtype=np.int64)
    for i in range(21):images|=((masks>>i)&1)[None,:]<<pg[:,i,None]
    for A in (3,4):
        choose=np.flatnonzero(As==A);allowed=masks[choose]
        assert np.all(np.isin(images[:,choose],allowed))
    assert Counter(seq(As))=={3:280,4:21}
    I=single['twenty_retained_incident'];active=walk['corner_orbits'][14]['active_walls']
    weights=np.array(list(map(int,single['affine_kernel'])),dtype=np.int64)
    coeff=np.array(list(map(int,single['projection_coefficients'])),dtype=np.int64)
    assert sum(abs(int(x)) for x in weights)*maxabs(W)<2**63
    assert sum(abs(int(x)) for x in coeff)*maxabs(W)+maxabs(h)<2**63
    D=np.einsum('i,gij->gj',weights,W[g[:,active]])
    PP=h-np.einsum('i,gij->gj',coeff,W[g[:,I]])
    p=PP[0];d=D[0]
    assert seq(p)==list(map(int,single['projection_p'])) and seq(d)==single['cusp_ray']
    assert len({row.tobytes() for row in D})==5760
    pair_bound(Q,PP,PP);pair_bound(Q,D,D);pair_bound(Q,PP,D)
    assert np.all(np.einsum('ij,ij->i',PP@Q,PP)==-2484)
    assert np.all(np.einsum('ij,ij->i',D@Q,D)==0)
    assert np.all(np.einsum('ij,ij->i',PP@Q,D)==-46)
    assert np.all(PP@(Q@h)==-2484) and np.all(D@(Q@h)==-46)
    A=PP@(Q@p);X=D@(Q@p);Y=PP@(Q@d);E=D@(Q@d)
    assert int(A[0])==-2484 and int(X[0])==int(Y[0])==-46 and int(E[0])==0
    assert np.all(E[1:]<0)
    # f_g(t)=(u(t),g*u(t))+3U(t). Negative means the larger,
    # same-height pairing is strictly dotted; opposite heights are smaller.
    t=350;U=92*t-2484;assert U==29716
    assert maxabs(A)+t*(maxabs(X)+maxabs(Y))+t*t*maxabs(E)+3*U<2**63
    f=A-t*(X+Y)+t*t*E+3*U
    df=-(X+Y)+2*t*E+276
    assert np.all(f[1:]<0) and np.all(df[1:]<0)
    # E<0 makes every derivative decrease for t>=350, proving the full tail.
    assert maxabs(PP)+350*maxabs(D)<2**63
    u=-p+t*d;V=p+296*d;Ug=-PP+t*D;Vg=PP+296*D
    assert len({row.tobytes() for row in Ug})==5760
    assert len({row.tobytes() for row in Vg})==5760
    pair_bound(Q,Ug,V);vp=Ug@(Q@V)
    assert int(vp[0])==0 and np.all(vp[1:]<0)
    # These strict vertex pairings also persist for every common t>=350.
    # v(t) dot g*u(t)=E*t²+(X-Y-54E)t+54Y-A.
    vd=2*t*E+X-Y-54*E
    assert np.array_equal(vp,t*t*E+t*(X-Y-54*E)+54*Y-A)
    assert np.all(vd[1:]<0)
    # Replay the representative old-wall and parent gates, and monotonicity
    # of -N(t)/sqrt(U(t)) on the whole t>=350 tail.
    WQ=W@Q;CQ=C@Q
    wa=-(WQ@p);wb=WQ@d;ca=-(CQ@p);cb=CQ@d
    assert np.all(wb[I]==0)
    for i in range(210):
        if i==deleted:continue
        n=int(wa[i]+t*wb[i])
        if i in I:assert n==0;continue
        assert n<0 and n*n>2*int(norms[i])*U
        assert int(wb[i])<0 and -int(wb[i])*U+46*n>0
    for a,b in zip(ca,cb):
        n=int(a)+t*int(b)
        assert parent_dotted(n,U) and int(b)<0 and -int(b)*U+46*n>0
    assert np.all((WQ@V)[I]==0)
    assert all(int(WQ[i]@V)<0 for i in range(210) if i not in I and i!=deleted)
    assert np.all(CQ@V<0)
    # All-real universal escape for arbitrary, independently chosen heights.
    # Norm positivity gives t>27. Forward old188 acuteness forces t>=96;
    # reverse orientation's old5 acuteness would force t<=4, impossible.
    assert (int(wa[188]),int(wb[188]))==(192,-2)
    assert (int(wa[5]),int(wb[5]))==(4,-1)
    assert 92*27-2484==0 and 4<27<96
    esc=np.array(single['strict_escape']['coordinates'],dtype=np.int64)
    assert int(esc@M@esc)==1586 and int(esc@M@h23)==-814
    assert np.all(np.delete(R@(M@esc),deleted)<0)
    pair_bound(Q,PP,esc[:22]);pair_bound(Q,D,esc[:22])
    ea=-(PP@(Q@esc[:22]));eb=D@(Q@esc[:22])
    assert np.all(eb<0)
    # For z=40: ea+t*eb +/-40sqrt(2U) <0 iff N<0 and N²>3200U.
    t0=96;U0=92*t0-2484
    assert maxabs(ea)+t0*maxabs(eb)<2**63
    n0=ea+t0*eb;assert np.all(n0<0)
    assert maxabs(n0)**2+3200*U0<2**63
    margin=n0*n0-3200*U0
    assert 2*maxabs(n0)*maxabs(eb)+3200*92<2**63
    slope=2*n0*eb-3200*92
    assert np.all(margin>0) and np.all(slope>0)
    # The squared margin has positive quadratic coefficient eb², so grows
    # on[96,infinity). Covers independently varying real t, both heights,
    # and arbitrary finite selections, with no compatibility assumption.
    ex350=ea+350*eb;assert np.all(ex350<0)
    # Stronger: the escape functional is strictly inside the cone of three
    # retained-wall functionals on every entire3D residual plane, with no
    # norm, center, parameter or Coxeter-label hypothesis.
    lookup={(r['A'],r['mask']):i for i,r in enumerate(recipes)}
    parent19=np.array([lookup[(int(As[19]),int(mask))] for mask in images[:,19]],dtype=np.int64)
    witness5=W[g[:,5]];witness188=W[g[:,188]];witnessC=C[parent19]
    for old,expected in ((witness5,(-4,-1)),(witness188,(-192,-2)),(witnessC,(-423,-8))):
        pair_bound(Q,old,PP);pair_bound(Q,old,D)
        assert np.all(np.einsum('ij,ij->i',old@Q,PP)==expected[0])
        assert np.all(np.einsum('ij,ij->i',old@Q,D)==expected[1])
    eP=-ea;eD=eb
    # beta=(-eP+4eD-15640)/184, alpha=-eD-320-2beta.
    assert maxabs(eP)+4*maxabs(eD)+15640<2**63
    beta184=-eP+4*eD-15640
    assert 184*(maxabs(eD)+320)+2*maxabs(beta184)<2**63
    alpha184=184*(-eD-320)-2*beta184
    assert np.all(alpha184>0) and np.all(beta184>0)
    assert 4*maxabs(alpha184)+192*maxabs(beta184)+423*40*184<2**63
    assert 184*max(maxabs(eP),maxabs(eD))<2**63
    assert np.all(-4*alpha184-192*beta184-423*40*184==184*eP)
    assert np.all(-alpha184-2*beta184-8*40*184==184*eD)
    assert 2*40==80
    assert (-4)*(-2)*2-(-1)*(-192)*2==-368
    # Thus E(r)=alpha W_g5(r)+beta W_g188(r)+40 R_g229(r).
    # If r is acute to those retained walls this is nonpositive. Equality
    # forces all three functionals to vanish, and their nonzero determinant
    # forces r=0 within the3D plane. Every nonzero candidate is strict.
    cusp=walk['corner_orbits'][0];delta0=np.array(cusp['primitive_ray'],dtype=np.int64)
    assert np.array_equal(esc[:22],100*delta0+h)
    assert int(delta0@Q@delta0)==0 and int(delta0@Q@h)==-8
    incident=np.flatnonzero(WQ@delta0==0)
    assert seq(incident)==cusp['active_walls'] and deleted not in incident and len(incident)==24
    assert np.all(WQ@delta0<=0) and np.all(CQ@delta0<0)
    pair_bound(Q,Ug,delta0);cusp_cap=Ug@(Q@delta0)
    assert np.all(cusp_cap<0)
    ranks=[]
    for comp in cusp['affine_components']:
        ids0=comp['walls'];weights0=list(map(F,comp['positive_kernel']))
        assert len(ids0)==comp['rank']+1 and min(weights0)>0
        summed=[sum(a*int(W[i,j]) for a,i in zip(weights0,ids0)) for j in range(22)]
        assert summed==[F(comp['ray_multiple'])*int(x) for x in delta0]
        positive_pivots([[int(G[i,j]) for j in ids0[:-1]] for i in ids0[:-1]])
        ranks.append(comp['rank'])
    assert sorted(i for comp in cusp['affine_components'] for i in comp['walls'])==seq(incident)
    assert sum(ranks)==20
    assert all(int(G[i,j])==0 for a in cusp['affine_components'] for b in cusp['affine_components']
        if a is not b for i in a['walls'] for j in b['walls'])
    vertical=np.array([0]*22+[1],dtype=np.int64)
    assert int(vertical@M@vertical)==2 and np.all(R[incident]@(M@vertical)==0)
    n=len(ids);new=2*n;right=20*new+n
    total=new*811+new*(new-1)//2
    profiles=Counter((int(a),int(x+y),int(e)) for a,x,y,e in zip(A[1:],X[1:],Y[1:],E[1:]))
    return {'status':'PASS_COMPLETE_5760_REPLACEMENT_ORBIT_AND_FULL_PLANE_ESCAPE',
        'sources':{s:sha(s) for s in (RAW,CORE,BANK,WALK,WALK_AUDIT,groupfile,SINGLE,
            'scripts/fresh_h22_core_replacement_20261013.py')},
        'deleted_old_wall':1,'subgroup_order':n,'subgroup_full_group_indices':seq(ids),
        'parent_bank_invariance':{'A3_masks':280,'A4_masks':21,'verified_images':n*301},
        'distinct_cusp_rays':n,'distinct_old_cap_projections':n,'signed_caps':new,
        'retained_base_facets':811,'facets':811+new,'parameter_t':350,
        'uniform_coherent_parameter_tail':'every common real t>=350',
        'new_pair_counts':{'right':right,'dotted':total-right},
        'representative_orbit_pair_polynomials':[{'a':a,'b':b,'c':e,'multiplicity':m}
            for (a,b,e),m in sorted(profiles.items())],
        'pair_formula':'(u,g*u)=a-b*t+c*t²; same-height dotted iff (u,g*u)+3U<0',
        'largest_pair_margin_at350':int(f[1:].max()),
        'largest_pair_margin_derivative_at350':int(df[1:].max()),
        'actual_finite_vertices':n,'largest_other_orbit_pairing_numerator_at_vertex':int(vp[1:].max()),
        'largest_vertex_pair_derivative_at350':int(vd[1:].max()),
        'specified_old_cusps_cut':n,'old_cusp_cut_numerator':46,
        'all_real_escape':{'coordinates':seq(esc),'square':1586,'center_pairing':-814,
            'forward_acute_necessary_t':'t>=96 from old188:192-2t<=0',
            'reverse_orientation_impossible':'old5 forces t<=4, but norm positivity forces t>27',
            'largest_old_numerator_at96':int(n0.max()),
            'smallest_squared_margin_at96':int(margin.min()),
            'smallest_margin_derivative_at96':int(slope.min()),
            'largest_old_numerator_at350':int(ex350.max()),
            'scope':'All finite selections of orbit normals, independently varying real t, either vertical sign, whenever each normal is acute to the811 retained base walls. This fixed-right family cannot give finite volume.'},
        'full_plane_obstruction':{'planes':5760,'plane_formula':'span{g*p,g*delta,vertical}',
            'three_retained_walls':'g*old5, g*old188, g*parent229 (positive height)',
            'functional_matrix':[[-4,-1,0],[-192,-2,0],[-423,-8,2]],
            'determinant':-368,'fixed_third_weight':40,'weight_denominator':184,
            'minimum_alpha_numerator':int(alpha184.min()),'minimum_beta_numerator':int(beta184.min()),
            'weights':[{'group_index':int(gid),'first_wall':int(perm[5]),'second_wall':int(perm[188]),
                'third_wall':210+int(k),'alpha_numerator':int(al),'beta_numerator':int(be)}
                for gid,perm,k,al,be in zip(ids,g,parent19,alpha184,beta184)],
            'scope':'Every nonzero real normal in any one ofthese5760planes, acute to the811 retained base normals, is strictly negative at the same escape. Arbitrary finite additions selected from their union cannot give finite volume, irrespective ofnorm, mutual cap angle, heights or finite labels. Linear combinations mixing different planes are not covered.'},
        'surviving_deficient_cusp':{'old_cusp_orbit':0,'ray':seq(delta0)+[0],
            'old_incident_walls':seq(incident),'affine_component_ranks':ranks,'quotient_rank':20,
            'horospherical_dimension':21,'lineality_direction':seq(vertical),'lineality_square':2,
            'maximum_parent_pairing':int((CQ@delta0).max()),
            'maximum_orbit_cap_pairing_numerator':int(cusp_cap.max()),
            'result':'The unchanged old cusp link has one free vertical Euclidean direction; all parent and replacement caps are nonincident. The link remains unbounded.'},
        'scope':'Complete stabilizer of deleted oldwall1 within the independently certified241920-element source group; explicitly preserves all602 signed parent caps. Exact12331 actual facets and5760 actual finite vertices at t350 (and the common parameter tail); no complete changed-chamber vertex/cusp census. The full-plane common escape excludes arbitrary acute additions from the union of5760residualplanes, not general old projections, additional old-wall deletions or replacement cores.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run()
    if args.verify:assert out==read(OUT)
    else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
