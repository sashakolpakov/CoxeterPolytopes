#!/usr/bin/env python3
"""Separate Fraction replay of the new A2+B6 full first-edge census."""
import argparse,hashlib,json
from functools import cmp_to_key
import fresh_h9_changed_spatial_20261013_independent as old
a=old.a
ROOT=old.ROOT;PREFIX='fresh_h9_changed_spatial_next_20261013'
SOURCE=ROOT/'results'/f'{PREFIX}.json'
OUTPUT=ROOT/'results'/f'{PREFIX}_independent.json'
D=a.scale(a.ALPHA,3)
P=old.P


def E(u,v):
    return(a.add(a.mul(u[0],v[0]),a.mul(D,a.mul(u[1],v[1]))),
           a.add(a.mul(u[0],v[1]),a.mul(u[1],v[0])))


def power(x,n):
    if n<0:return power(a.inv(x),-n)
    y=a.O
    for _ in range(n):y=a.mul(y,x)
    return y


def check():
    d=json.loads(SOURCE.read_text())
    for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    pp=ROOT/'results'/f'{PREFIX}_profiles.json';pd=json.loads(pp.read_text())
    edge=old.matrix(pd['edge']);C=tuple(tuple(P(r,s)for s in edge)for r in edge);ci=a.minverse(C)
    assert old.matrix(pd['gram'])==C
    prior=json.loads((ROOT/'results/fresh_h9_changed_spatial_20261013.json').read_text())
    initial=old.matrix(prior['initial_walls']);assert edge==initial[:2]+initial[3:]
    found={}
    # Reversed profile recursion independently reconstitutes all finite
    # pairings; only the six fixed tail coordinates must be integral.
    for q in old.POOL:
        alph=[]
        for r in edge:
            choices={a.Z}
            for lam in a.ALLOWED:
                for p in old.quartic_roots(a.scale(a.mul(a.mul(q,P(r,r)),lam),a.Q(1,4))):
                    if a.sign(p)>=0 and a.integral(p)and a.integral(a.divide(a.scale(p,2),q)):choices.add(p)
            alph.append(choices)
        def visit(pre):
            n=len(pre);ids=tuple(range(7,7-n,-1))
            if n:
                sub=tuple(tuple(ci[i][j]for j in ids)for i in ids)
                if a.sign(a.sub(q,a.dot(pre,a.mv(sub,pre))))<=0:return
            if n==8:
                mags=tuple(reversed(pre));cs=a.mv(ci,tuple(map(a.neg,mags)))
                r0=tuple(old.summ(a.mul(c,r[j])for c,r in zip(cs,edge))for j in range(10))
                assert r0[0]==a.Z and old.summ(r0[1:4])==a.Z
                if not all(map(a.integral,r0[4:])):return
                N=a.scale(a.sub(q,P(r0,r0)),3)
                if any(a.sign(N,e)<=0 for e in range(4)):return
                assert a.integral(N)
                found[q,mags]=(r0,N);return
            for p in alph[7-n]:visit((*pre,p))
        visit(())
    saved={(a.dec(r['q']),a.vec(r['magnitudes'])):(a.vec(r['fixed']),a.dec(r['target']))for r in pd['profiles']}
    assert found==saved and len(found)==d['profiles']==28
    raw=(ROOT/'results'/f'{PREFIX}_raw.txt').read_text()
    for s in('ABSOLUTE x^8 + 12*x^6 + 18*x^4 - 108*x^2 - 81\n',
             'CERTIFY 1\n','BASE_CERTIFY 1\n','CLASS 1\n','SIGNATURE [2, 3]\n'):assert s in raw
    base=[];units=[];normrows=[];absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('BASE_UNIT '):base.append(a.dec(json.loads(line.split(' ',2)[2])))
        elif line.startswith('UNIT '):units.append(a.vec(json.loads(line.split(' ',2)[2])))
        elif line.startswith('UNIT_NORM '):
            n,e=line.split(' ',2)[2].split('] ',1);normrows.append((a.dec(json.loads(n+']')),json.loads(e)))
        elif line.startswith('NORM '):
            _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
        elif line.startswith('SOLUTION '):
            _,n,r=line.split(' ',2);absolute[int(n)].append(a.vec(json.loads(r)))
    assert base==list(a.vec(d['base_units']))and units==[a.vec(u)for u in d['fundamental_units']]
    assert [e for n,e in normrows]==d['relative_norm_map']==[[0,2,0,0],[-2,0,0,0],[0,0,2,0],[0,0,0,0]]
    for (z,x),(n,e)in zip(units,normrows):
        assert a.sub(a.sq(z),a.mul(D,a.sq(x)))==n
        v=a.scale(a.O,(-1)**e[-1])
        for b,k in zip(base,e):v=a.mul(v,power(b,k))
        assert v==n
    eps=units[3]
    if a.sign(eps[0])<0:eps=tuple(map(a.neg,eps))
    if a.sign(eps[1])<0:eps=eps[0],a.neg(eps[1])
    assert eps==a.vec(d['relative_unit'])and a.sign(eps[1])>0 and a.sign(a.sub(eps[0],a.O))>0
    assert a.sub(a.sq(eps[0]),a.mul(D,a.sq(eps[1])))==a.O
    # The first three norm exponent columns have determinant 8. Thus their
    # integer kernel is zero, and the full relative kernel is +/-unit4^Z.
    assert set(absolute)=={int(old.norm(N))for r0,N in saved.values()}
    assert all(len(v)==counts[n]for n,v in absolute.items())
    upper=a.dec(d['upper_parameter']);back=eps[0],a.neg(eps[1]);visited=[];valid=[];intervals=[]
    for i,row in enumerate(pd['profiles']):
        q,N=a.dec(row['q']),a.dec(row['target']);r0=a.vec(row['fixed']);reps=[]
        for z,x in absolute[int(old.norm(N))]:
            nu=a.sub(a.sq(z),a.mul(D,a.sq(x)));assert old.norm(nu)==old.norm(N)
            delta=a.divide(nu,N)
            if not a.integral(delta)or not a.integral(a.inv(delta)):continue
            assert all(a.sign(delta,e)>0 for e in range(4))
            roots=[s for s in old.quartic_roots(delta)if a.sign(s)>0];assert len(roots)==1;s=roots[0]
            assert a.integral(s)and a.integral(a.inv(s));z=a.divide(z,s);x=a.divide(x,s)
            if a.sign(z)<0:z=a.neg(z);x=a.neg(x)
            reps.append((z,x))
        for j,u in enumerate(reps):
            k=0
            while a.sign(u[1])<=0:u=E(u,eps);k+=1;assert abs(k)<50
            while a.sign(E(u,back)[1])>0:u=E(u,back);k-=1;assert abs(k)<50
            previous=E(u,back);assert a.sign(previous[1])<=0
            first=k
            while a.sign(a.sub(a.mul(upper,u[0]),a.mul(a.ALPHA,u[1])))>=0:
                z,x=u;r=(x,*(a.add(r0[t],a.scale(z,a.Q(1,3)))for t in range(1,4)),*r0[4:])
                why=old.reason(r,q);parameter=a.divide(a.mul(a.ALPHA,x),z)
                visited.append((i,j,k,r,why))
                if why is None:valid.append((parameter,r))
                u=E(u,eps);k+=1;assert k-first<=200
            intervals.append({'profile':i,'orbit':j,'first':first,'last':k-1,
                              'previous':list(map(a.enc,previous)),'next':list(map(a.enc,u))})
    assert visited==[(r['profile'],r['orbit'],r['power'],a.vec(r['root']),r['reason'])for r in d['visited_points']]
    assert valid and all(t==upper for t,r in valid)
    chosen=a.vec(d['chosen_root']);assert chosen in [r for t,r in valid]
    frame=edge+(chosen,);M=tuple(tuple(P(r,s)for s in frame)for r in frame);mi=a.minverse(M)
    for j,(_,r)in enumerate(valid):
        cs=a.mv(mi,tuple(P(r,s)for s in frame));assert all(a.sign(c)>=0 for c in cs)
        assert tuple(old.summ(a.mul(c,s[t])for c,s in zip(cs,frame))for t in range(10))==r
        assert cs==a.vec(d['simultaneous_root_cone_coefficients'][j])
    walls=old.matrix(d['all_walls']);inside=a.vec(d['strict_interior']);escape=a.vec(d['strict_escape']);point=a.vec(d['endpoint'])
    assert len(walls)==20 and walls[-1]==chosen
    for i,r in enumerate(walls):
        assert old.reason(r,P(r,r))is None
        for s in walls[:i]:old.gate(r,s)
    assert point==(a.O,upper,upper,upper,*(a.Z for _ in range(6)))
    assert a.sign(P(point,point))<0 and all(a.sign(P(point,r))<=0 for r in walls)
    assert [i for i,r in enumerate(walls)if P(point,r)==a.Z]==[0,1,3,4,5,6,7,8,19]
    assert P(chosen,chosen)==a.scale(a.O,4)
    assert tuple(P(chosen,r)for r in edge)==(a.Z,a.scale(a.O,-2),*(a.Z for _ in range(6)))
    assert d['endpoint_type']=='B3+B6'and d['endpoint_mirrors']==45
    z,x=eps
    assert chosen==(a.scale(x,2),a.scale(a.add(a.scale(z,2),a.scale(a.O,2)),a.Q(1,3)),
        a.scale(a.add(a.scale(z,2),a.scale(a.O,2)),a.Q(1,3)),
        a.scale(a.sub(a.scale(z,2),a.scale(a.O,4)),a.Q(1,3)),*(a.Z for _ in range(6)))
    assert a.sign(P(inside,inside))<0 and all(a.sign(P(inside,r))<0 for r in walls)
    assert P(escape,escape)==a.sub(a.scale(a.O,a.Q(2537,1250)),a.T)==a.dec(d['strict_escape_norm'])
    assert a.sign(P(escape,escape))>0 and a.sign(P(escape,inside))<0
    assert all(a.sign(P(escape,r))<0 for r in walls)
    cp=json.loads((ROOT/'results/fresh_h9_changed_spatial_20261013_checkpoint.json').read_text())
    assert a.sign(P(a.vec(cp['strict_escape']),chosen))>0
    return {'status':'PASS_INDEPENDENT_COMPLETE_D9_A2_B6_FIRST_EDGE',
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'profiles':len(found),
        'interval_points':len(visited),'first_mirrors':len(valid),'walls':len(walls),'pairs':190,
        'endpoint_type':'B3+B6','endpoint_mirrors':45,'strict_escape_norm':a.enc(P(escape,escape)),
        'intervals':intervals,
        'scope':'Separate Fraction arithmetic, reversed complete profile recursion, full relative-unit normalization and interval replay, all20 root/190 angle gates and exact finite corner/escape. Certified PARI absolute-norm lists and full fundamental units are the explicit shared completeness dependency.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=check()
    if args.verify:assert out==json.loads(OUTPUT.read_text())
    else:OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['profiles'],out['interval_points'],flush=True)
