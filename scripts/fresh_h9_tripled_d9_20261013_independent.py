#!/usr/bin/env python3
"""Separate Fraction arithmetic for the tripled-D9 complete initial edges."""
import argparse,hashlib,json
import fresh_h9_changed_spatial_20261013_independent as old
a=old.a;ROOT=old.ROOT
PREFIX='fresh_h9_tripled_d9_20261013'
BETA=a.scale(a.ALPHA,3)
G=tuple(tuple(a.neg(BETA)if i==j==0 else a.O if i==j else a.Z for j in range(10))for i in range(10))
POOL=(*old.POOL,*(a.scale(q,3)for q in old.POOL))
OUT=ROOT/'results'/f'{PREFIX}_independent.json'


def P(r,s):return a.dot(r,a.mv(G,s))
def reason(r,q):
    assert P(r,r)==q
    if not old.lattice(r):return'NONINTEGRAL'
    if any(old.lattice(tuple(a.divide(x,p)for x in r))for p in(a.T,a.scale(a.O,3))):return'IMPRIMITIVE'
    if not all(a.integral(a.divide(a.scale(P(r,b),2),q))for b in old.BASIS):return'NONREFLECTIVE'
    return None
def gate(r,s):
    p=P(r,s);assert a.sign(p)<=0
    lam=a.divide(a.scale(a.sq(p),4),a.mul(P(r,r),P(s,s)))
    assert lam in a.ALLOWED or a.sign(a.sub(lam,a.scale(a.O,4)))>=0
def power(x,n):
    if n<0:return power(a.inv(x),-n)
    ans=a.O
    for _ in range(n):ans=a.mul(ans,x)
    return ans


def profiles(pd,k,axis_free=None):
    edge=old.matrix(pd['edge']);C=tuple(tuple(P(r,s)for s in edge)for r in edge);ci=a.minverse(C)
    assert C==old.matrix(pd['gram']);found={}
    for q in POOL:
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
                assert r0[0]==a.Z
                if axis_free is None:
                    assert old.summ(r0[1:k+1])==a.Z
                    if not all(map(a.integral,r0[k+1:])):return
                else:
                    assert r0[axis_free]==a.Z
                    if not all(map(a.integral,r0)):return
                N=a.scale(a.sub(q,P(r0,r0)),k)
                if any(a.sign(N,e)<=0 for e in range(4)):return
                assert a.integral(N);found[q,mags]=(r0,N);return
            for p in alph[7-n]:visit((*pre,p))
        visit(())
    saved={(a.dec(r['q']),a.vec(r['magnitudes'])):(a.vec(r['fixed']),a.dec(r['target']))for r in pd['profiles']}
    assert found==saved
    return len(found)


def run(verify=False):
    sp=ROOT/'results'/f'{PREFIX}_structure.json';st=json.loads(sp.read_text())
    for p,h in st['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    assert a.vec(st['norm_pool'])==POOL and old.matrix(st['ambient_gram'])==G
    assert old.matrix(st['lattice_basis'])==tuple(old.BASIS)
    gram=tuple(tuple(P(r,s)for s in old.BASIS)for r in old.BASIS);gi=a.minverse(gram)
    assert gram==old.matrix(st['lattice_gram'])
    assert all(a.integral(a.scale(x,12))for row in gi for x in row)
    assert any(not a.integral(a.scale(x,6))for row in gi for x in row)
    assert any(not a.integral(a.scale(x,4))for row in gi for x in row)
    assert all(a.integral(a.divide(a.scale(a.O,24),q))and all(a.sign(q,e)>0 for e in range(4))for q in POOL)
    assert [old.norm(q)for q in POOL]==[b**4*2**j for b in(1,3)for j in range(13)]
    sources={str(sp.relative_to(ROOT)):hashlib.sha256(sp.read_bytes()).hexdigest()};results=[];added=[];initial=None
    for k in(1,3):
        tag=f'{PREFIX}_edge{k}';path=ROOT/'results'/f'{tag}.json';d=json.loads(path.read_text())
        sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
        for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
        pd=json.loads((ROOT/'results'/f'{tag}_profiles.json').read_text());nprofiles=profiles(pd,k)
        assert nprofiles==d['profiles'];field=BETA if k==1 else a.ALPHA;timefactor=1 if k==1 else 3
        assert a.dec(d['field_parameter'])==field and d['time_divisor']==timefactor
        edge=old.matrix(pd['edge']);walls=old.matrix(d['all_walls'])
        if initial is None:initial=walls[:9]
        assert initial==walls[:9]and edge==initial[:k-1]+initial[k:]
        raw=(ROOT/'results'/f'{tag}_raw.txt').read_text()
        for s in('CERTIFY 1\n','BASE_CERTIFY 1\n','CLASS 1\n','SIGNATURE [2, 3]\n'):assert s in raw
        polynomial='x^8 + 12*x^6 + 18*x^4 - 108*x^2 - 81'if k==1 else'x^8 + 4*x^6 + 2*x^4 - 4*x^2 - 1'
        assert 'ABSOLUTE '+polynomial+'\n'in raw
        bases=[];units=[];norms=[];absolute={};counts={}
        for line in raw.splitlines():
            if line.startswith('BASE_UNIT '):bases.append(a.dec(json.loads(line.split(' ',2)[2])))
            elif line.startswith('UNIT '):units.append(a.vec(json.loads(line.split(' ',2)[2])))
            elif line.startswith('UNIT_NORM '):
                n,e=line.split(' ',2)[2].split('] ',1);norms.append((a.dec(json.loads(n+']')),json.loads(e)))
            elif line.startswith('NORM '):
                _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
            elif line.startswith('SOLUTION '):
                _,n,s=line.split(' ',2);absolute[int(n)].append(a.vec(json.loads(s)))
        assert bases==list(a.vec(d['base_units']))and units==[a.vec(u)for u in d['fundamental_units']]
        expected=[[0,2,0,0],[-2,0,0,0],[0,0,2,0],[0,0,0,0]]if k==1 else[[0,-1,0,0],[2,-1,0,0],[0,0,2,0],[0,0,0,0]]
        assert[e for n,e in norms]==expected==d['relative_norm_map']
        for(z,x),(n,e)in zip(units,norms):
            assert a.sub(a.sq(z),a.mul(field,a.sq(x)))==n
            ans=a.scale(a.O,(-1)**e[-1])
            for b,j in zip(bases,e):ans=a.mul(ans,power(b,j))
            assert n==ans
        eps=units[3]
        if a.sign(eps[0])<0:eps=tuple(map(a.neg,eps))
        if a.sign(eps[1])<0:eps=eps[0],a.neg(eps[1])
        assert eps==a.vec(d['relative_unit'])and a.sign(a.sub(eps[0],a.O))>0 and a.sign(eps[1])>0
        assert a.sub(a.sq(eps[0]),a.mul(field,a.sq(eps[1])))==a.O
        def emul(u,v):return a.add(a.mul(u[0],v[0]),a.mul(field,a.mul(u[1],v[1]))),a.add(a.mul(u[0],v[1]),a.mul(u[1],v[0]))
        def parameter(u):return a.divide(a.mul(BETA,u[1]),a.scale(u[0],timefactor))
        assert set(absolute)==set(pd['absolute_norms'])and all(len(v)==counts[n]for n,v in absolute.items())
        upper=a.dec(d['upper_parameter']);back=eps[0],a.neg(eps[1]);visited=[];valid=[]
        for i,row in enumerate(pd['profiles']):
            q,N=a.dec(row['q']),a.dec(row['target']);r0=a.vec(row['fixed']);orbits=[]
            for z,x in absolute[int(old.norm(N))]:
                nu=a.sub(a.sq(z),a.mul(field,a.sq(x)));assert old.norm(nu)==old.norm(N)
                ratio=a.divide(nu,N)
                if not a.integral(ratio)or not a.integral(a.inv(ratio)):continue
                assert all(a.sign(ratio,e)>0 for e in range(4))
                roots=[u for u in old.quartic_roots(ratio)if a.sign(u)>0];assert len(roots)==1;s=roots[0]
                assert a.integral(s)and a.integral(a.inv(s));z=a.divide(z,s);x=a.divide(x,s)
                if a.sign(z)<0:z=a.neg(z);x=a.neg(x)
                orbits.append((z,x))
            for j,u in enumerate(orbits):
                powern=0
                while a.sign(u[1])<=0:u=emul(u,eps);powern+=1;assert abs(powern)<50
                while a.sign(emul(u,back)[1])>0:u=emul(u,back);powern-=1;assert abs(powern)<50
                assert a.sign(emul(u,back)[1])<=0
                first=powern
                while a.sign(a.sub(parameter(u),upper))<=0:
                    z,x=u;r=(a.scale(x,a.Q(1,timefactor)),*(a.add(r0[t],a.scale(z,a.Q(1,k)))for t in range(1,k+1)),*r0[k+1:])
                    why=reason(r,q);visited.append((i,j,powern,r,why))
                    if why is None:valid.append((parameter(u),r))
                    u=emul(u,eps);powern+=1;assert powern-first<=200
        assert visited==[(r['profile'],r['orbit'],r['power'],a.vec(r['root']),r['reason'])for r in d['visited_points']]
        assert valid and all(t==upper for t,r in valid)
        chosen=a.vec(d['chosen_root']);assert chosen in[r for t,r in valid]and chosen==walls[-1]
        frame=edge+(chosen,);M=tuple(tuple(P(r,s)for s in frame)for r in frame);mi=a.minverse(M)
        for j,(_,r)in enumerate(valid):
            cs=a.mv(mi,tuple(P(r,s)for s in frame));assert all(a.sign(c)>=0 for c in cs)
            assert tuple(old.summ(a.mul(c,s[t])for c,s in zip(cs,frame))for t in range(10))==r
            assert cs==a.vec(d['simultaneous_root_cone_coefficients'][j])
        point=a.vec(d['endpoint']);inside=a.vec(d['strict_interior'])
        assert point==(a.O,*(upper for _ in range(k)),*(a.Z for _ in range(9-k)))and a.sign(P(point,point))<0
        assert all(P(point,r)==a.Z for r in frame)and a.sign(P(point,initial[k-1]))<0
        assert a.sign(P((a.O,*(a.Z for _ in range(9))),chosen))<0
        assert a.sign(P(inside,inside))<0 and all(a.sign(P(inside,r))<0 for r in walls)
        for i,r in enumerate(walls):
            assert reason(r,P(r,r))is None
            for s in walls[:i]:gate(r,s)
        added.append(chosen)
        results.append({'edge':k,'profiles':nprofiles,'interval_points':len(visited),'first_mirrors':len(valid),
                        'norm':a.enc(P(chosen,chosen))})
    walls=initial+tuple(added)
    for i,r in enumerate(walls):
        for s in walls[:i]:gate(r,s)
    inside=(a.O,*(a.scale(a.O,a.Q(9-i,32))for i in range(9)))
    escape=(a.O,*(a.scale(a.O,a.Q(x,100))for x in(145,75,50,40,30,20,15,10,5)))
    assert a.sign(P(inside,inside))<0 and all(a.sign(P(inside,r))<0 for r in walls)
    assert P(escape,escape)==a.sub(a.scale(a.O,a.Q(156,25)),a.scale(a.T,3))
    assert a.sign(P(escape,escape))>0 and a.sign(P(escape,inside))<0 and all(a.sign(P(escape,r))<0 for r in walls)
    out={'status':'PASS_INDEPENDENT_TRIPLED_D9_COMPLETE_INITIAL_EDGES','sources':sources,
        'norm_classes':26,'edges':results,'total_profiles':sum(r['profiles']for r in results),
        'total_interval_points':sum(r['interval_points']for r in results),'combined_walls':11,'combined_pairs':55,
        'strict_interior':list(map(a.enc,inside)),'strict_escape':list(map(a.enc,escape)),
        'strict_escape_norm':a.enc(P(escape,escape)),
        'scope':'Independent Fraction field arithmetic, reversed complete profile enumeration, full unit-orbit normalization and interval census, D9 primitivity atT and3, all root/angle gates and combined geometry. PARI full absolute-norm representatives and fundamental-unit completeness remain a shared certified dependency. Full-host reflectivity is not decided.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['total_profiles'],out['total_interval_points'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();run(args.verify)
