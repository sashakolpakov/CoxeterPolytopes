#!/usr/bin/env python3
"""Separate Fraction profiles, unit intervals and D9 root/geometry replay."""
import argparse
import hashlib
import json
from functools import cmp_to_key
import fresh_h9_glued_host_20261013_independent as a
from fresh_h9_quartic_first_endpoint_20261013_independent import quartic_roots
from fresh_h9_glued_host_20261013_structure import norm

ROOT=a.ROOT;PREFIX='fresh_h9_changed_spatial_20261013'
G=tuple(tuple(a.neg(a.ALPHA)if i==j==0 else a.O if i==j else a.Z for j in range(10))for i in range(10))
POOL=(*a.POOL,a.scale(a.PI,2),a.scale(a.sq(a.T),2),a.divide(a.scale(a.O,4),a.PI),a.scale(a.O,4),
      a.scale(a.PI,4),a.scale(a.sq(a.T),4),a.divide(a.scale(a.O,8),a.PI),a.scale(a.O,8))


def P(x,y):return a.dot(x,a.mv(G,y))
def matrix(x):return tuple(a.vec(r)for r in x)
def emul(u,v):return a.add(a.mul(u[0],v[0]),a.mul(a.ALPHA,a.mul(u[1],v[1]))),a.add(a.mul(u[0],v[1]),a.mul(u[1],v[0]))
def summ(xs):
    out=a.Z
    for x in xs:out=a.add(out,x)
    return out
def lattice(r):return all(map(a.integral,r))and a.integral(a.scale(summ(r[1:]),a.Q(1,2)))
def basis():
    out=[(a.O,*(a.Z for _ in range(9)))]
    for j in range(1,9):
        r=[a.Z]*10;r[j]=a.O;r[j+1]=a.neg(a.O);out.append(tuple(r))
    r=[a.Z]*10;r[8]=r[9]=a.O;out.append(tuple(r));return out
BASIS=basis()
def reason(r,q):
    assert P(r,r)==q
    if not lattice(r):return 'NONINTEGRAL'
    if lattice(tuple(a.divide(x,a.T)for x in r)):return 'IMPRIMITIVE'
    if not all(a.integral(a.divide(a.scale(P(r,b),2),q))for b in BASIS):return 'NONREFLECTIVE'
    return None
def gate(r,s):
    p=P(r,s);assert a.sign(p)<=0
    ratio=a.divide(a.scale(a.sq(p),4),a.mul(P(r,r),P(s,s)))
    assert a.sign(a.sub(ratio,a.scale(a.O,4)))>=0 or ratio in a.ALLOWED


def independent_profiles(pd,integral_fixed=True):
    edge=matrix(pd['edge']);C=tuple(tuple(P(r,s)for s in edge)for r in edge);ci=a.minverse(C)
    assert matrix(pd['gram'])==C
    saved={(a.dec(r['q']),a.vec(r['magnitudes'])):(a.vec(r['fixed']),a.dec(r['target']))for r in pd['profiles']}
    found={}
    for q in POOL:
        alph=[]
        for r in edge:
            candidates={a.Z}
            for lam in a.ALLOWED:
                for p in quartic_roots(a.scale(a.mul(a.mul(q,P(r,r)),lam),a.Q(1,4))):
                    if a.sign(p)>=0 and a.integral(p)and a.integral(a.divide(a.scale(p,2),q)):candidates.add(p)
            alph.append(candidates)
        def walk(pre):
            n=len(pre);indices=list(range(7,7-n,-1))
            if n:
                Q=tuple(tuple(ci[i][j]for j in indices)for i in indices)
                if a.sign(a.sub(q,a.dot(pre,a.mv(Q,pre))))<=0:return
            if n==8:
                mags=tuple(reversed(pre));coef=a.mv(ci,tuple(map(a.neg,mags)))
                r0=tuple(summ(a.mul(c,r[j])for c,r in zip(coef,edge))for j in range(10))
                if integral_fixed and not all(map(a.integral,r0)):return
                N=a.sub(q,P(r0,r0))
                if not integral_fixed:N=a.scale(N,9)
                if any(a.sign(N,e)<=0 for e in range(4)):return
                found[q,mags]=(r0,N);return
            for x in alph[7-n]:walk((*pre,x))
        walk(())
    assert found==saved
    return len(found)


def run(stages,verify=False,diagonal=False):
    first=ROOT/'results'/f'{PREFIX}.json';fd=json.loads(first.read_text())
    assert a.vec(fd['norm_pool'])==POOL and matrix(fd['lattice_basis'])==tuple(BASIS)
    epsilon=a.vec(fd['relative_unit']);back=(epsilon[0],a.neg(epsilon[1]))
    assert a.sub(a.sq(epsilon[0]),a.mul(a.ALPHA,a.sq(epsilon[1])))==a.O
    assert a.sign(epsilon[1])>0 and a.sign(a.sub(epsilon[0],a.O))>0
    assert fd['relative_norm_map']==[[0,-1,0,0],[2,-1,0,0],[0,0,2,0],[0,0,0,0]]
    units=[a.vec(u)for u in fd['fundamental_units']];bases=a.vec(fd['base_units'])
    def powr(x,n):
        if n<0:return powr(a.inv(x),-n)
        y=a.O
        for _ in range(n):y=a.mul(y,x)
        return y
    for (z,x),es in zip(units,fd['relative_norm_map']):
        n=a.scale(a.O,(-1)**es[-1])
        for b,e in zip(bases,es):n=a.mul(n,powr(b,e))
        assert a.sub(a.sq(z),a.mul(a.ALPHA,a.sq(x)))==n
    assert units[3]==epsilon
    data=[];sources={};total_points=0;total_profiles=0
    for stage in [*range(1,stages+1),*([10]if diagonal else [])]:
        diag=stage==10
        path=first if stage==1 else ROOT/'results'/(f'{PREFIX}_diagonal.json'if diag else f'{PREFIX}_step{stage}.json')
        d=json.loads(path.read_text());sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
        for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
        if stage==1:
            ps=[]
            for i,q in enumerate(POOL):ps.append((q,q,(a.Z,)*10))
            r0=[a.Z]*10;r0[2]=a.O;ps.append((a.scale(a.O,2),a.O,tuple(r0)))
            assert len(d['profiles'])==14
            for row,(q,N,r0)in zip(d['profiles'],ps):assert a.dec(row['norm'])==q and a.dec(row['relative_target'])==N
            AI=tuple(tuple(a.O if i==j else a.Z for j in range(10))for i in range(10));free=1
        else:
            pd=json.loads((ROOT/'results'/(f'{PREFIX}_diagonal_profiles.json'if diag else f'{PREFIX}_step{stage}_profiles.json')).read_text())
            assert independent_profiles(pd,not diag)==d['profiles']
            ps=[(a.dec(r['q']),a.dec(r['target']),a.vec(r['fixed']))for r in pd['profiles']]
            if not diag:
                AI=matrix(d['frame_from_axis']);A=matrix(d['frame_to_axis']);free=pd['free_coordinate']
                assert a.minverse(A)==AI and all(a.integral(x)for row in A+AI for x in row)
                assert tuple(tuple(P(a.mv(A,e),a.mv(A,f))for f in BASIS)for e in BASIS)==tuple(tuple(P(e,f)for f in BASIS)for e in BASIS)
        rawpath=ROOT/'results'/(f'{PREFIX}_d9_raw.txt'if stage==1 else f'{PREFIX}_diagonal_raw.txt'if diag else f'{PREFIX}_step{stage}_raw.txt')
        raw=rawpath.read_text();assert 'CERTIFY 1\n'in raw and'BASE_CERTIFY 1\n'in raw
        assert 'ABSOLUTE x^8 + 4*x^6 + 2*x^4 - 4*x^2 - 1\n'in raw
        fundamental_lines=lambda text:[s for s in text.splitlines()if s.startswith(('UNIT ','UNIT_NORM ','BASE_UNIT '))]
        assert fundamental_lines(raw)==fundamental_lines((ROOT/'results'/f'{PREFIX}_d9_raw.txt').read_text())
        absolute={};counts={}
        for line in raw.splitlines():
            if line.startswith('NORM '):
                _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
            elif line.startswith('SOLUTION '):
                _,n,r=line.split(' ',2);absolute[int(n)].append(tuple(a.dec(x)for x in json.loads(r)))
        assert set(absolute)=={int(norm(N))for q,N,r0 in ps}
        assert all(len(v)==counts[n]for n,v in absolute.items())
        upper=a.dec(d['upper_parameter']);visited=[];valid=[]
        for i,(q,N,r0)in enumerate(ps):
            reps=[]
            for z,x in absolute[int(norm(N))]:
                nu=a.sub(a.sq(z),a.mul(a.ALPHA,a.sq(x)));assert norm(nu)==norm(N)
                delta=a.divide(nu,N)
                if not a.integral(delta)or not a.integral(a.inv(delta)):continue
                assert all(a.sign(delta,e)>0 for e in range(4))
                roots=[s for s in quartic_roots(delta)if a.sign(s)>0];assert len(roots)==1;s=roots[0]
                assert a.integral(s)and a.integral(a.inv(s));z=a.divide(z,s);x=a.divide(x,s)
                if a.sign(z)<0:z=a.neg(z);x=a.neg(x)
                reps.append((z,x))
            for j,u in enumerate(reps):
                power=0
                while a.sign(u[1])<=0:u=emul(u,epsilon);power+=1;assert abs(power)<30
                while a.sign(emul(u,back)[1])>0:u=emul(u,back);power-=1;assert abs(power)<30
                assert a.sign(emul(u,back)[1])<=0
                while a.sign(a.sub(a.mul(upper,u[0]),a.mul(a.ALPHA,u[1])))>=0:
                    z,x=u
                    if diag:r=(a.scale(x,a.Q(1,3)),*(a.add(y,a.scale(z,a.Q(1,9)))for y in r0[1:]))
                    else:
                        rr=list(r0);rr[0]=x;rr[free]=z;r=a.mv(AI,rr)
                    why=reason(r,q)
                    visited.append((i,j,power,r,why))
                    if why is None:valid.append((a.divide(a.mul(a.ALPHA,x),z),r))
                    u=emul(u,epsilon);power+=1;assert abs(power)<200
        expected=[(r['profile'],r['orbit'],r['power'],a.vec(r['root']),r['reason'])for r in d['visited_points']]
        assert visited==expected
        valid.sort(key=cmp_to_key(lambda x,y:a.sign(a.sub(x[0],y[0]))))
        assert valid
        nearest=[r for p,r in valid if p==valid[0][0]]
        chosen=a.vec(d['first_mirror']if stage==1 else d['chosen_root']);assert chosen in nearest
        walls=matrix(d['all_walls'])
        for i,r in enumerate(walls):
            q=P(r,r);assert reason(r,q)is None
            for s in walls[:i]:gate(r,s)
        point=a.vec(d['first_endpoint']if stage==1 else d['endpoint'])
        assert a.sign(P(point,point))<0 and all(a.sign(P(point,w))<=0 for w in walls)
        if stage==1:
            incident=walls[1:9]+(chosen,)
            assert len(nearest)==1
        else:
            incident=matrix(d['endpoint_incident']);ci=a.minverse(tuple(tuple(P(r,s)for s in incident)for r in incident))
            for r,encoded in zip(nearest,d['simultaneous_root_cone_coefficients']):
                coeff=a.mv(ci,tuple(P(r,e)for e in incident));assert coeff==a.vec(encoded)
                assert all(a.sign(c)>=0 for c in coeff)
                assert tuple(summ(a.mul(c,e[j])for c,e in zip(coeff,incident))for j in range(10))==r
        assert len(incident)==9 and all(P(point,r)==a.Z for r in incident)
        # Timelike orthogonality plus an invertible Gram proves spherical rank9.
        a.minverse(tuple(tuple(P(r,s)for s in incident)for r in incident))
        total_points+=len(visited);total_profiles+=len(ps)
        data.append({'stage':stage,'profiles':len(ps),'interval_points':len(visited),
                     'first_contact_mirrors':len(nearest),'wall_count':len(walls)})
    result={'status':'PASS_INDEPENDENT_D9_COMPLETE_EDGE_CHAIN','sources':sources,'stages':data,
            'total_profiles':total_profiles,'total_interval_points':total_points,
            'arithmetic':'Separate Fraction field operations, independent quadratic-tower square roots, reversed complete profile recursion, and independently normalized full unit intervals. Certified PARI absolute norm representatives/full units are an explicit shared CAS dependency.',
            'scope':'All saved first contacts, primitive roots, inherited finite-angle gates, simultaneous-mirror cone choices and actual finite endpoints. No finite-volume assertion.'}
    if diagonal:
        checkpoint=ROOT/'results'/f'{PREFIX}_checkpoint.json';ck=json.loads(checkpoint.read_text())
        for path,h in ck['sources'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
        inside=a.vec(ck['strict_interior']);escape=a.vec(ck['strict_escape'])
        assert a.sign(P(inside,inside))<0 and a.sign(P(escape,escape))>0 and a.sign(P(inside,escape))<0
        assert P(escape,escape)==a.dec(ck['strict_escape_norm'])
        assert all(a.sign(P(inside,r))<0 and a.sign(P(escape,r))<0 for r in walls)
        ev=a.vec(ck['next_initial_edge']['spacelike_continuation'])
        assert a.sign(P(ev,ev))>0 and all(a.sign(P(ev,r))<=0 for r in walls)
        assert [i for i in range(9)if P(ev,walls[i])==a.Z]==ck['next_initial_edge']['support']
        for rec in [ck['opposite_B9_transport'],*ck['B2_diagram_swap_tests']]:
            H=matrix(rec['matrix'])
            assert a.minverse(H)==H
            assert tuple(tuple(P(a.mv(H,e),a.mv(H,f))for f in BASIS)for e in BASIS)==tuple(tuple(P(e,f)for f in BASIS)for e in BASIS)
            assert any(not lattice(a.mv(H,b))for b in BASIS)
        result.update(checkpoint_source=str(checkpoint.relative_to(ROOT)),checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                      final_facets=len(walls),strict_escape_norm=ck['strict_escape_norm'],
                      final_geometry='Independent strict timelike center, strict future spacelike escape, open A2+B6 edge and three nonintegral involution checks PASS.')
    op=ROOT/'results'/f'{PREFIX}_independent.json'
    if verify:assert result==json.loads(op.read_text())
    else:op.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'stages',len(data),'profiles',total_profiles,'points',total_points,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stages',type=int,default=1);p.add_argument('--verify',action='store_true');p.add_argument('--diagonal',action='store_true');args=p.parse_args();run(args.stages,args.verify,args.diagonal)
