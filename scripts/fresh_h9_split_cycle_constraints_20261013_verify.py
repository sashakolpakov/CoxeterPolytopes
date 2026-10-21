#!/usr/bin/env python3
"""Producer-free exact replay of the cycle cap and its continuous family."""
import argparse
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from scripts.fresh_h9_split_endpoint_walk_exact_20261013_base import (
    ROOT,K,Quadratic,Z,O,P,source,read,parse,enc,scale,add,dot,matvec,pairing,
    lift_K,embed,lift_existing,key,relation,pair_profiles,angle_pool,ANGLE_SOURCE,strict_escape)
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K

PREFIX='data/enumerations/fresh_h9_split_cycle_constraints_20261013'
OUT=ROOT/'results/fresh_h9_split_cycle_constraints_20261013_verify.json'


def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def old_pair(F,c,row):
    E,er=embed(F,row['field'])
    return E,dot(E,tuple(lift_existing(E,F,x)for x in c),tuple(er(x)for x in row['p']))


def verify_cycle(data,ctx):
    G,I,h,old,E,F=ctx;pool=angle_pool(ANGLE_SOURCE)
    assert data['source_old_walls']==list(old)
    assert data['seed_gram']==[[enc(x)for x in r]for r in G]
    assert data['seed_gram_inverse']==[[enc(x)for x in r]for r in I]
    assert data['center_covector']==list(map(enc,h))
    c=tuple(map(F.parse,data['normal_coefficients']));p=tuple(map(F.parse,data['normal_covector']))
    assert len(c)==len(p)==10
    assert all(F.sign(F.sub(x,y))==0 for x,y in zip(matvec(F,[[lift_K(F,z)for z in r]for r in G],c),p,strict=True))
    assert F.sign(F.sub(dot(F,c,p),F.rational(2)))==0
    cp=dot(F,c,tuple(lift_K(F,x)for x in h));assert F.sign(cp)>0 and F.encode(cp)==data['strict_center_pairing']
    row={'name':'cycle2009','field':F,'c':c,'p':p,'q':scale(O,2)}
    checks={n:relation(row,r,pool)for n,r in old.items()}
    assert checks==data['pair_relations']and all(v['kind']in('finite','dotted')for v in checks.values())
    for n,want in zip(data['local_equation_walls'],[0]*8+[-1],strict=True):
        ff,z=old_pair(F,c,old[n]);assert ff.sign(ff.sub(z,ff.rational(want)))==0
    u=tuple(map(E.parse,data['line_constant']));v=tuple(map(E.parse,data['line_direction']))
    IE=[[lift_K(E,x)for x in rr]for rr in I]
    A=pairing(E,v,IE,v);B=E.mul(E.rational(2),pairing(E,u,IE,v));C=E.sub(pairing(E,u,IE,u),E.rational(2))
    assert list(map(E.encode,(A,B,C)))==data['norm_polynomial']
    assert E.sign(E.sub(E.sub(E.mul(B,B),E.mul(E.rational(4),E.mul(A,C))),F.d))==0
    t=F.parse(data['parameter']);assert data['branch']==-1
    expected=(E.div(E.neg(B),E.mul(E.rational(2),A)),E.div(E.rational(-1),E.mul(E.rational(2),A)))
    assert F.sign(F.sub(t,expected))==0
    assert all(F.sign(F.sub(z,F.add(F.lift(x),F.mul(t,F.lift(y)))))==0 for z,x,y in zip(p,u,v,strict=True))
    selected=old|{row['name']:row};IF=[[lift_K(F,x)for x in rr]for rr in I]
    for proof in data['finite_vertices']:
        y=tuple(map(F.parse,proof['y']));assert len(y)==10
        square=pairing(F,y,IF,y);cp=pairing(F,y,IF,tuple(lift_K(F,x)for x in h))
        assert F.sign(square)<0 and F.sign(cp)<0
        assert F.encode(square)==proof['square']and F.encode(cp)==proof['center_pairing']
        signs={}
        for n,r in selected.items():
            ff,er=embed(F,r['field'])
            z=dot(ff,tuple(er(x)for x in r['c']),tuple(lift_existing(ff,F,x)for x in y))
            signs[n]=ff.sign(z)
        incident=proof['incident'];assert len(incident)==len(set(incident))==9
        assert signs==proof['all_wall_signs']and all(s==(0 if n in incident else 1)for n,s in signs.items())
        gram=[[parse(x)for x in r]for r in proof['gram']];assert inertia_K(gram)==(9,0,0)
        for i,n in enumerate(incident):
            for j,m in enumerate(incident):
                ff,z=pair_profiles(selected[n],selected[m]);assert ff.sign(ff.sub(z,lift_K(ff,gram[i][j])))==0
    vertices=data['finite_vertices'];assert len(vertices)==4
    for proof,(i,j)in zip(data['four_cycle_edges'],((0,1),(1,2),(2,3),(3,0)),strict=True):
        assert proof['vertices']==[i,j]
        common=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']));assert len(common)==8 and common==proof['common_walls']
        mid=tuple(F.add(F.parse(x),F.parse(y))for x,y in zip(vertices[i]['y'],vertices[j]['y'],strict=True))
        q=pairing(F,mid,IF,mid);assert F.sign(q)<0 and F.encode(q)==proof['midpoint_square']
    escape=data['old_escape_survives'];assert strict_escape(tuple(map(parse,escape['y'])),I,h,list(selected.values()))==escape
    return row


def verify_family(data,cycle,ctx):
    G,I,h,old,E,F=ctx;IE=[[lift_K(E,x)for x in rr]for rr in I]
    a,b=[tuple(map(E.parse,data[n]))for n in('a_covector','b_covector')]
    ac,bc=matvec(E,IE,a),matvec(E,IE,b)
    assert list(map(E.encode,ac))==data['a_coefficients']and list(map(E.encode,bc))==data['b_coefficients']
    kap,delta=E.parse(data['kappa']),E.parse(data['delta']);assert E.sign(kap)>0 and E.sign(delta)>0
    aa,ab,bb=pairing(E,a,IE,a),pairing(E,a,IE,b),pairing(E,b,IE,b)
    assert E.sign(ab)==0 and E.sign(E.add(aa,E.mul(kap,bb)))==0
    assert E.sign(E.sub(E.mul(delta,bb),E.rational(2)))==0
    assert E.encode(aa)==data['a_square']and E.encode(bb)==data['b_square']and E.encode(ab)==data['a_b_pairing']
    zeros=data['orthogonal_edge_walls'];assert len(zeros)==8
    for n,r in old.items():
        ff,x=old_pair(E,ac,r);gg,y=old_pair(E,bc,r);assert key(ff)==key(gg)
        signs={'linear_sign':ff.sign(x),'radical_sign':ff.sign(y)}
        assert signs==data['old_wall_coefficient_signs'][n]
        assert signs==({'linear_sign':0,'radical_sign':0}if n in zeros else {'linear_sign':-1,'radical_sign':0}if n=='wall22'else {'linear_sign':-1,'radical_sign':-1})
        assert r['field'].sign(dot(r['field'],r['c'],tuple(lift_K(r['field'],z)for z in h)))==data['old_strict_center_pairing_signs'][n]==1
        for j,m in enumerate(zeros):
            fld,z=pair_profiles(r,old[m]);sgn=fld.sign(z);assert sgn==data['old_wall_to_edge_pairing_signs'][n][j]
            if n in zeros:assert fld.sign(fld.sub(z,fld.rational(2 if n==m else 0)))==0
            else:assert sgn<=0
    for c,saved in zip((ac,bc),data['center_coefficient_pairings'],strict=True):
        z=dot(E,c,tuple(lift_K(E,x)for x in h));assert E.sign(z)>0 and E.encode(z)==saved
    F2=Quadratic(E,E.add(E.mul(E.rational(4),kap),delta));p2=tuple((E.mul(E.rational(2),x),z)for x,z in zip(a,b,strict=True))
    c2=matvec(F2,[[lift_K(F2,x)for x in rr]for rr in I],p2)
    row2={'field':F2,'p':p2,'c':c2,'q':scale(O,2)};pool=angle_pool(ANGLE_SOURCE)
    relations={n:relation(row2,r,pool)for n,r in old.items()};assert relations==data['lambda2_relations']
    assert all(v==({'kind':'finite','label':2}if n in zeros else {'kind':'parallel'}if n=='wall22'else {'kind':'dotted'})for n,v in relations.items())
    y=tuple(scale(O,Q(v))for v in data['remote_integer_covector']);assert len(y)==10
    assert strict_escape(y,I,h,list(old.values()))==data['remote_strict_escape']
    for c,saved in zip((ac,bc),data['remote_coefficient_pairings'],strict=True):
        z=dot(E,c,tuple(lift_K(E,x)for x in y));assert E.sign(z)>0 and E.encode(z)==saved
    dom=data['domination'];c,d=E.parse(dom['c']),E.parse(dom['d']);assert E.sign(c)>0 and E.sign(d)<0
    z=tuple(map(E.parse,dom['z_covector']))
    assert all(E.sign(E.sub(zz,E.add(E.mul(c,x),E.mul(d,y))))==0 for zz,x,y in zip(z,a,b,strict=True))
    c17=tuple(E.lift(x)for x in old['wall17']['p'])
    expected=tuple(E.add(lift_K(E,x),E.mul(lift_K(E,scale(P,Q(1,2))),t))for x,t in zip(old['seed3']['p'],c17,strict=True))
    assert all(E.sign(E.sub(x,y))==0 for x,y in zip(z,expected,strict=True))
    assert E.sign(E.sub(E.add(kap,delta),F.d))==0
    assert all(F.sign(F.sub((x,y),F.parse(p)))==0 for x,y,p in zip(a,b,cycle['normal_covector'],strict=True))
    assert F.sign(F.sub(F.lift(d),(E.zero,c)))<0


def run():
    cycle=read(PREFIX+'.json');family=read(PREFIX+'_continuous.json')
    for data in(cycle,family):assert all(sha(p)==v for p,v in data['source_sha256'].items())
    G,I,h,old,_,_=source();assert inertia_K(G)==(9,1,0)
    E0=Quadratic(K,parse(cycle['base_radicands'][0]));E=Quadratic(E0,E0.lift(parse(cycle['base_radicands'][1])));F=Quadratic(E,E.parse(cycle['norm_discriminant']))
    assert E0.d==old['wall17']['field'].d and E.d==E0.lift(old['wall22']['field'].d)
    ctx=(G,I,h,old,E,F);verify_cycle(cycle,ctx);verify_family(family,cycle,ctx)
    bad=deepcopy(cycle);bad['normal_coefficients'][0]=F.encode(F.add(F.parse(bad['normal_coefficients'][0]),F.one))
    try:verify_cycle(bad,ctx)
    except AssertionError:pass
    else:raise AssertionError('Forged cap accepted')
    bad=deepcopy(family);bad['domination']['c']=E.encode(E.neg(E.parse(bad['domination']['c'])))
    try:verify_family(bad,cycle,ctx)
    except AssertionError:pass
    else:raise AssertionError('Forged domination coefficient accepted')
    return {'status':'PASS_PRODUCER_FREE_EXACT_CYCLE_AND_CONTINUOUS_FAMILY',
        'source_sha256':{p:sha(p)for p in(PREFIX+'.json',PREFIX+'_continuous.json')},
        'new_cap_norm_two':True,'retained_wall_relations':72,'actual_finite_vertices':4,
        'actual_cycle_edges':4,'old_open_edges_closed_simultaneously':2,
        'continuous_family_domain':'lambda>2','continuous_family_actual_all_right_vertex':True,
        'whole_family_strict_spacelike_escape':True,'whole_family_dominated_by_lambda1':True,
        'forged_cap_rejected':True,'forged_domination_coefficient_rejected':True,
        'scope':'Literal cycle and universal one-parameter continuation only. No complete chamber or unrestricted nonexistence claim.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();r=run()
    if args.verify:assert json.loads(OUT.read_text())==r
    else:OUT.write_text(json.dumps(r,indent=2)+'\n')
    print(r,flush=True)
