#!/usr/bin/env python3
"""Exact two coupled polar repairs retaining the certified96-wall H20 cone.

The numerical archives supply active-set indices only. All normals, target
points, angles, and the final strict escape are reconstructed exactly.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from fresh_h20_edge_family_completion_20261013_four_volume import setup as setup94,ROOT,read,sha,check
from fresh_dimension_bridge_20261013_verify import E,Z,O,parse,vec,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv,same,label,rank
from fresh_h20_edge_family_completion_20261013_cusp_repairs import base_root

PREFIX='fresh_h20_gap_next_20261013'
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
SOURCE='data/enumerations/fresh_h20_edge_family_completion_20261013_cusp_repairs.json'
FIRST_SCOUT='data/enumerations/'+PREFIX+'_coupled.json'
SECOND_SCOUT='data/enumerations/'+PREFIX+'_second_coupled.json'
ACTIVE_SOURCE='data/enumerations/fresh_h20_edge_family_completion_20261013_four_volume.json'
OUT='data/enumerations/'+PREFIX+'_exact.json'

def setup():
    _,_,_,Q,h,hc,roots=setup94()
    roots.extend((B,vec(c['normal']))for c in read(SOURCE)['new_caps'])
    assert len(roots)==96
    return Q,h,hc,roots

def kernel(rows):
    a=[list(row)for row in rows];piv=[];n=len(a[0]);k=0
    for j in range(n):
        v=next((i for i in range(k,len(a))if a[i][j].sign()),None)
        if v is None:continue
        a[k],a[v]=a[v],a[k];z=a[k][j];a[k]=[x/z for x in a[k]]
        for i in range(len(a)):
            if i==k:continue
            z=a[i][j]
            if z.sign():a[i]=[x-z*y for x,y in zip(a[i],a[k],strict=True)]
        piv.append(j);k+=1
        if k==len(a):break
    assert len(piv)==n-1
    free=next(j for j in range(n)if j not in piv);v=[Z]*n;v[free]=O
    for i,j in enumerate(piv):v[j]=-a[i][free]
    assert all(dot(row,v)==Z for row in rows)
    return tuple(v),piv

def base_point(Q,hc,roots,ids,anchor=None,pairing=None):
    assert len(ids)==len(set(ids))and all(0<=i<len(roots)for i in ids)
    rs=[base_root(*roots[i])for i in ids];p,piv=kernel([mv(Q,r)for r in rs])
    if anchor is not None:
        z=dot(p,mv(Q,base_root(*roots[anchor])));assert z.sign()
        p=tuple(x*pairing/z for x in p)
    elif dot(p,hc).sign()>0:p=tuple(-x for x in p)
    pc=mv(Q,p);q=dot(p,pc);assert q.sign()>0 and dot(p,hc).sign()<0
    pairs=[];zeros=[]
    for i,(f,r)in enumerate(roots):
        z=fdot(f,r,pc if f is B else tuple(map(f.lift,pc)));assert f.sign(z)<=0
        if not f.sign(z):zeros.append(i)
        pairs.append(f.encode(z))
    assert all(i in zeros for i in ids)
    return p,{'vector':[x.enc()for x in p],'norm':q.enc(),'center_pairing':dot(p,hc).enc(),'selected_active_indices':ids,'kernel_rank':len(piv),'all_wall_pairings':pairs,'zero_indices':zeros}

def root_checks(f,r,Q,hc,roots):
    Qf=Q if f is B else tuple(tuple(map(f.lift,row))for row in Q)
    rc=fmv(f,Qf,r);q=fdot(f,r,rc);assert f.sign(q)>0
    ch=fdot(f,r,hc if f is B else tuple(map(f.lift,hc)));assert f.sign(ch)<0
    rows=[]
    for i,(g,s)in enumerate(roots):
        Mg=Q if g is B else tuple(tuple(map(g.lift,row))for row in Q)
        qs=fdot(g,s,fmv(g,Mg,s))
        if g is B:
            p=fdot(f,rc,s if f is B else tuple(map(f.lift,s)));H=f
            qnew=q;qold=qs if f is B else f.lift(qs)
        elif f is B:
            p=fdot(g,tuple(map(g.lift,rc)),s);H=g;qnew=g.lift(q);qold=qs
        else:
            H=Extension(f,f.lift(g.d));ss=tuple((f.lift(a),f.lift(b))for a,b in s)
            p=fdot(H,tuple(map(H.lift,rc)),ss);qnew=H.lift(q);qold=(f.lift(qs[0]),f.lift(qs[1]))
        rows.append({'old':i,'pairing':H.encode(p),'kind':label(H,p,qnew,qold)})
    return q,ch,rows

def strictify(v,Q,h,hc,roots,new=None):
    for power in range(12,81):
        x=tuple(a+E(F(1,2**power))*b for a,b in zip(v,h,strict=True))
        try:c=check(x,Q,hc,roots)
        except AssertionError:continue
        if new is not None:
            f,r=new;xc=mv(Q,x);cut=fdot(f,r,xc if f is B else tuple(map(f.lift,xc)))
            if f.sign(cut)<=0:continue
            c['new_pairing']=f.encode(cut)
        c['perturbation_power']=power;return c
    raise AssertionError('no exact strict perturbation')

def cusp_signs(f,r,Q):
    out=[]
    for c in read(BASE)['cusp_rows']:
        cv=mv(Q,vec(c['ray']));p=fdot(f,r,cv if f is B else tuple(map(f.lift,cv)))
        assert f.sign(p)<0;out.append(f.encode(p))
    assert len(out)==107;return out

def run():
    Q,h,hc,roots=setup();a=base_root(*roots[3]);b=tuple(x/E(F(3))for x in base_root(*roots[42]))
    active1=read(ACTIVE_SOURCE)['certificate']['numeric_active_indices']
    p1,pcert1=base_point(Q,hc,roots,active1,3,E(F(-2)))
    assert dot(p1,mv(Q,p1))==E(F(3))and dot(p1,mv(Q,b))==E(F(-2))
    assert dot(a,mv(Q,a))==dot(b,mv(Q,b))==E(F(2))and dot(a,mv(Q,b))==Z
    f=Extension(B,E(F(10)));t=(E(F(4,3)),E(F(1,3)))
    r1=tuple(f.add(f.lift(x+y),f.mul(t,f.lift(p)))for x,y,p in zip(a,b,p1,strict=True))
    q,ch,pairs=root_checks(f,r1,Q,hc,roots);assert same(f,q,f.rational(2))
    target1=read(SOURCE)['strict_escape'];v1=vec(target1['vector']);actual=check(v1,Q,hc,roots)
    assert all(actual[k]==target1[k]for k in actual)
    cut1=fdot(f,r1,tuple(map(f.lift,mv(Q,v1))));assert f.sign(cut1)>0
    c1={'wall':96,'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':['10','0']},'normal':list(map(f.encode,r1)),'norm':f.encode(q),'strict_center_pairing':f.encode(ch),'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(x['kind']for x in pairs).items())),'strict_target':dict(actual,new_pairing=f.encode(cut1)),'inherited_cusp_strict_pairings':cusp_signs(f,r1,Q),'formula':'r3+r42/3+((4+sqrt10)/3)*p1','polar_certificate':pcert1}
    roots.append((f,r1))
    active2=read(FIRST_SCOUT)['escape']['active_indices']
    p2,pcert2=base_point(Q,hc,roots,active2,5,E(F(-1)))
    leaves=[base_root(*roots[i])for i in(5,14)]+[tuple(x/E(F(3))for x in base_root(*roots[36]))]
    assert dot(p2,mv(Q,p2))==O
    assert [[dot(x,mv(Q,y))for y in leaves]for x in leaves]==[[E(F(2))if i==j else Z for j in range(3)]for i in range(3)]
    assert all(dot(p2,mv(Q,x))==-O for x in leaves)
    r2=tuple(E(F(5))*p+sum((x[k]for x in leaves),Z)for k,p in enumerate(p2))
    q,ch,pairs=root_checks(B,r2,Q,hc,roots);assert q==O
    target2=strictify(p2,Q,h,hc,roots,(B,r2))
    c2={'wall':97,'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':['1','0']},'normal':[[x.enc(),Z.enc()]for x in r2],'norm':[q.enc(),Z.enc()],'strict_center_pairing':[ch.enc(),Z.enc()],'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(x['kind']for x in pairs).items())),'strict_target':target2,'inherited_cusp_strict_pairings':cusp_signs(B,r2,Q),'formula':'5*p2+r5+r14+r36/3','polar_certificate':pcert2}
    roots.append((B,r2))
    active3=read(SECOND_SCOUT)['escape']['active_indices']
    p3,pcert3=base_point(Q,hc,roots,active3)
    final=strictify(p3,Q,h,hc,roots)
    return {'status':'EXACT_COHERENT_98_WALL_H20_COUPLED_POLAR_REPAIRS_WITH_STRICT_ESCAPE','source_sha256':{p:sha(p)for p in(BASE,SOURCE,ACTIVE_SOURCE,FIRST_SCOUT,SECOND_SCOUT)},'retained_walls':96,'wall_count':98,'new_caps':[c1,c2],'final_hyperideal_point':pcert3,'strict_escape':final,'scope':'Two coupled local polar repairs retain all96 certified walls. Exact normals,193 new pairs,214 inherited cusp signs,two strict target cuts,and final98 strict escape. Numerical sources select active-set indices only; exactrank20 reconstruction supplies every polar and final escape. No finite-volume or complete vertex/cusp claim.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert out==read(OUT)
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print([r['pair_counts']for r in out['new_caps']]);print('escape',out['strict_escape']['norm'])
