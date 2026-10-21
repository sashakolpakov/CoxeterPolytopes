#!/usr/bin/env python3
"""Exact all-parameter dotted-anchor continuation of the cycle cap.

For every lambda>2 the archived radical formula gives a compatible all-right
endpoint cap. A single strict spacelike witness survives the whole family;
positive-cone identities also show domination by the lambda=1 cycle cap.
"""
import argparse
import json
from scripts.fresh_h9_split_cycle_constraints_20261013 import (
    ROOT,K,Quadratic,O,P,Q,source,read,enc,parse,scale,dot,matvec,pairing,
    lift_K,key,relation,angle_pool,ANGLE_SOURCE,strict_escape,sha,pair_with_old,PREFIX)
from scripts.fresh_h9_split_endpoint_walk_exact_20261013_base import pair_profiles

OUT=ROOT/(PREFIX+'_continuous.json')
REMOTE='data/enumerations/fresh_h9_split_multiple_edge_caps_20261013_escape_scout.json'
Y=[17553546,197944749,203785283,203785283,47178720,47178720,47178720,207849122,186231837,160670502]


def build():
    cycle=read(PREFIX+'.json');G,I,h,old,_,_=source()
    E0=Quadratic(K,parse(cycle['base_radicands'][0]))
    E=Quadratic(E0,E0.lift(parse(cycle['base_radicands'][1])))
    u=tuple(map(E.parse,cycle['line_constant']));v=tuple(map(E.parse,cycle['line_direction']))
    A,B,C=map(E.parse,cycle['norm_polynomial']);assert cycle['branch']==-1
    a=tuple(E.sub(x,E.mul(E.div(B,E.mul(E.rational(2),A)),z))for x,z in zip(u,v,strict=True))
    b=tuple(E.mul(E.div(E.rational(-1),E.mul(E.rational(2),A)),z)for z in v)
    kappa=E.sub(E.mul(B,B),E.mul(E.rational(4),E.mul(A,E.add(C,E.rational(2)))))
    delta=E.mul(E.rational(8),A)
    assert E.sign(A)>0 and E.sign(kappa)>0 and E.sign(delta)>0
    IE=tuple(tuple(lift_K(E,x)for x in rr)for rr in I)
    ac,bc=matvec(E,IE,a),matvec(E,IE,b)
    aa=pairing(E,a,IE,a);bb=pairing(E,b,IE,b);ab=pairing(E,a,IE,b)
    assert E.sign(ab)==0
    assert E.sign(E.add(aa,E.div(kappa,E.mul(E.rational(4),A))))==0
    assert E.sign(E.sub(bb,E.div(E.one,E.mul(E.rational(4),A))))==0
    assert E.sign(E.sub(E.mul(delta,bb),E.rational(2)))==0
    # Hence (lambda a+sqrt(kappa lambda²+delta)b)^2=2 identically.
    coeffs={};zero=[];decreasing=[]
    for name,row in old.items():
        f,x=pair_with_old(E,ac,row);ff,z=pair_with_old(E,bc,row);assert key(f)==key(ff)
        sx,sz=f.sign(x),f.sign(z)
        assert sx<=0 and sz<=0
        if sx==sz==0:zero.append(name)
        else:decreasing.append(name)
        if name=='wall22':assert f.sign(f.sub(x,f.rational(-1)))==0 and sz==0
        coeffs[name]={'linear_sign':sx,'radical_sign':sz}
    assert zero==[f'seed{i}'for i in(0,1,2,4,5,6,7)]+['wall17']
    assert len(decreasing)==64
    edge_incidence_signs={};old_center_signs={}
    for name,row in old.items():
        old_center_signs[name]=row['field'].sign(dot(row['field'],row['c'],tuple(lift_K(row['field'],x)for x in h)))
        assert old_center_signs[name]>0
        signs=[]
        for n in zero:
            f,z=pair_profiles(row,old[n]);signs.append(f.sign(z))
            if name in zero:
                assert f.sign(f.sub(z,f.rational(2 if name==n else 0)))==0
            else:assert signs[-1]<=0
        edge_incidence_signs[name]=signs
    ah=dot(E,ac,tuple(lift_K(E,x)for x in h));bh=dot(E,bc,tuple(lift_K(E,x)for x in h))
    assert E.sign(ah)>0 and E.sign(bh)>0
    D2=E.add(E.mul(E.rational(4),kappa),delta);F2=Quadratic(E,D2)
    p2=tuple((E.mul(E.rational(2),x),z)for x,z in zip(a,b,strict=True))
    c2=matvec(F2,tuple(tuple(lift_K(F2,x)for x in rr)for rr in I),p2)
    cap2={'name':'continuous_boundary_lambda2','field':F2,'p':p2,'c':c2,'q':scale(O,2)}
    pool=angle_pool(ANGLE_SOURCE);relations={n:relation(cap2,r,pool)for n,r in old.items()}
    assert all(r==({'kind':'finite','label':2}if n in zero else {'kind':'parallel'}if n=='wall22'else {'kind':'dotted'})for n,r in relations.items())
    # Strict monotonicity makes every off-block pairing strictly <−2 for lambda>2.
    y=tuple(scale(O,Q(v))for v in Y);escape=strict_escape(y,I,h,list(old.values()))
    ay=dot(E,ac,tuple(lift_K(E,x)for x in y));by=dot(E,bc,tuple(lift_K(E,x)for x in y))
    assert E.sign(ay)>0 and E.sign(by)>0
    # Exact domination: z=seed3+(phi/2)wall17 has nonnegative old-wall coefficients.
    c17=tuple(E.lift(x)for x in old['wall17']['p'])
    z=tuple(E.add(lift_K(E,x),E.mul(lift_K(E,scale(P,Q(1,2))),v))for x,v in zip(old['seed3']['p'],c17,strict=True))
    cc=E.div(pairing(E,z,IE,a),aa);dd=E.div(pairing(E,z,IE,b),bb)
    assert all(E.sign(E.sub(zz,E.add(E.mul(cc,x),E.mul(dd,w))))==0 for zz,x,w in zip(z,a,b,strict=True))
    assert E.sign(cc)>0 and E.sign(dd)<0
    D1=E.add(kappa,delta);F1=Quadratic(E,D1)
    assert E.sign(E.sub(D1,E.parse(cycle['norm_discriminant'])))==0
    assert all(F1.sign(F1.sub((x,w),F1.parse(v)))==0 for x,w,v in zip(a,b,cycle['normal_covector'],strict=True))
    denominator=F1.sub(F1.lift(dd),(E.zero,cc));assert F1.sign(denominator)<0
    # For lambda>1, sqrt(k lambda²+d)<lambda sqrt(k+d). Together with
    # c>0,d<0 this proves both displayed cone coefficients strictly positive.
    return {'status':'EXACT_CONTINUOUS_DOTTED_CAP_FAMILY_AND_UNIVERSAL_ESCAPE',
        'source_sha256':{p:sha(p)for p in(PREFIX+'.json',REMOTE)},
        'parameter_domain':'lambda>2',
        'formula':'normal_covector(lambda)=lambda*a+sqrt(kappa*lambda^2+delta)*b',
        'a_covector':list(map(E.encode,a)),'b_covector':list(map(E.encode,b)),
        'a_coefficients':list(map(E.encode,ac)),'b_coefficients':list(map(E.encode,bc)),
        'kappa':E.encode(kappa),'delta':E.encode(delta),'a_square':E.encode(aa),
        'b_square':E.encode(bb),'a_b_pairing':E.encode(ab),
        'orthogonal_edge_walls':zero,'old_wall_coefficient_signs':coeffs,
        'old_wall_to_edge_pairing_signs':edge_incidence_signs,
        'old_strict_center_pairing_signs':old_center_signs,
        'center_coefficient_pairings':[E.encode(ah),E.encode(bh)],
        'lambda2_relations':relations,'all_parameter_norm':2,
        'all_parameter_finite_vertex_gram':'2*I9: eight orthogonal old edge normals plus the new normal',
        'remote_integer_covector':Y,'remote_strict_escape':escape,
        'remote_coefficient_pairings':[E.encode(ay),E.encode(by)],
        'domination':{'old_positive_combination':'z=seed3+(phi/2)*wall17',
            'z_covector':list(map(E.encode,z)),'c':E.encode(cc),'d':E.encode(dd),
            'c_sign':1,'d_sign':-1,'denominator_sign':-1,
            'formula':'r(lambda)=alpha*r(1)+beta*z',
            'alpha':'(lambda*d-c*sqrt(kappa*lambda^2+delta))/(d-c*sqrt(kappa+delta))',
            'beta':'(sqrt(kappa*lambda^2+delta)-lambda*sqrt(kappa+delta))/(d-c*sqrt(kappa+delta))',
            'coefficient_domain':'alpha>0 and beta>0 for every lambda>1'},
        'scope':'A genuine continuous family with only right angles at its eight-wall endpoint block, all64 other old relations strictly dotted, and a full-dimensional actual all-right finite vertex. The fixed remote strict spacelike escape survives every member. Moreover each lambda>2 inequality is implied by the lambda=1 cycle cap and old seed3/wall17 inequalities, so this family cannot improve the cycle chamber. No statement about other edge profiles or all split-seed completions.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();r=build()
    if args.verify:assert r==read(PREFIX+'_continuous.json')
    else:OUT.write_text(json.dumps(r,indent=2)+'\n')
    print({'status':r['status'],'orthogonal':len(r['orthogonal_edge_walls']),'strictly_dotted':64,
           'universal_escape':True,'dominated_by_cycle_cap':True},flush=True)
