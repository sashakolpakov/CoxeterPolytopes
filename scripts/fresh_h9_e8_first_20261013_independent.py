#!/usr/bin/env python3
"""Separate Fraction norm-orbit and geometry replay for the E8-fixed edge."""
import argparse,hashlib,json
from pathlib import Path
import fresh_h9_changed_spatial_20261013_independent as old
a=old.a;ROOT=old.ROOT;PREFIX='fresh_h9_e8_first_20261013'
BETA=a.scale(a.ALPHA,3)
G=tuple(tuple(a.neg(BETA)if i==j==0 else a.O if i==j else a.Z for j in range(10))for i in range(10))
POOL=(*old.POOL[:5],*(a.scale(q,3)for q in old.POOL[:5]))
OUT=ROOT/'results'/f'{PREFIX}_independent.json'
def P(r,s):return a.dot(r,a.mv(G,s))
def power(x,n):
    if n<0:return power(a.inv(x),-n)
    ans=a.O
    for _ in range(n):ans=a.mul(ans,x)
    return ans
def run(verify=False):
    source=ROOT/'results'/f'{PREFIX}.json';d=json.loads(source.read_text())
    assert d['status']=='PASS_COMPLETE_E8_FIXED_FIRST_EDGE'
    for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    basis=old.matrix(d['lattice_basis']);B=tuple(zip(*basis));BI=a.minverse(B)
    gram=tuple(tuple(P(r,s)for s in basis)for r in basis);gi=a.minverse(gram)
    assert all(a.integral(a.scale(x,3))for row in gi for x in row)
    assert any(not a.integral(x)for row in gi for x in row)
    assert a.vec(d['safe_norm_pool'])==POOL and len(POOL)==10
    assert all(a.integral(a.divide(a.scale(a.O,6),q))and all(a.sign(q,e)>0 for e in range(4))for q in POOL)
    assert[old.norm(q)for q in POOL]==[b**4*2**j for b in(1,3)for j in range(5)]
    assert[(a.dec(r['q']),a.dec(r['target']),a.vec(r['magnitudes']))for r in d['profiles']]==[(q,q,(a.Z,)*8)for q in POOL]
    def integral(v):return all(map(a.integral,a.mv(BI,v)))
    def reason(r,q):
        assert P(r,r)==q
        if not integral(r):return'NONINTEGRAL'
        if any(integral(tuple(a.divide(x,p)for x in r))for p in(a.T,a.scale(a.O,3))):return'IMPRIMITIVE'
        if not all(a.integral(a.divide(a.scale(P(r,b),2),q))for b in basis):return'NONREFLECTIVE'
        return None
    raw=(ROOT/'results'/f'{PREFIX}_raw.txt').read_text()
    for s in('ABSOLUTE x^8 + 12*x^6 + 18*x^4 - 108*x^2 - 81\n','CERTIFY 1\n','BASE_CERTIFY 1\n','CLASS 1\n','SIGNATURE [2, 3]\n'):assert s in raw
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
    expected=[[0,2,0,0],[-2,0,0,0],[0,0,2,0],[0,0,0,0]]
    assert[e for n,e in norms]==expected==d['relative_norm_map']
    for(z,x),(n,e)in zip(units,norms):
        assert a.sub(a.sq(z),a.mul(BETA,a.sq(x)))==n
        ans=a.scale(a.O,(-1)**e[-1])
        for b,j in zip(bases,e):ans=a.mul(ans,power(b,j))
        assert n==ans
    eps=units[3]
    if a.sign(eps[0])<0:eps=tuple(map(a.neg,eps))
    if a.sign(eps[1])<0:eps=eps[0],a.neg(eps[1])
    assert eps==a.vec(d['relative_unit'])and a.sign(a.sub(eps[0],a.O))>0 and a.sign(eps[1])>0
    assert a.sub(a.sq(eps[0]),a.mul(BETA,a.sq(eps[1])))==a.O
    def emul(u,v):return a.add(a.mul(u[0],v[0]),a.mul(BETA,a.mul(u[1],v[1]))),a.add(a.mul(u[0],v[1]),a.mul(u[1],v[0]))
    def param(u):return a.divide(a.mul(BETA,u[1]),u[0])
    upper=a.dec(d['upper_parameter']);back=eps[0],a.neg(eps[1]);visited=[];valid=[]
    assert set(absolute)=={int(old.norm(q))for q in POOL}
    assert all(len(v)==counts[n]for n,v in absolute.items())
    assert{str(n):c for n,c in counts.items()}==d['absolute_norm_counts']
    # Reverse the norm classes; reconstruct every corrected absolute orbit.
    for i in range(9,-1,-1):
        q=POOL[i];orbits=[]
        for z,x in absolute[int(old.norm(q))]:
            nu=a.sub(a.sq(z),a.mul(BETA,a.sq(x)));assert old.norm(nu)==old.norm(q)
            ratio=a.divide(nu,q)
            if not a.integral(ratio)or not a.integral(a.inv(ratio)):continue
            assert all(a.sign(ratio,e)>0 for e in range(4))
            roots=[u for u in old.quartic_roots(ratio)if a.sign(u)>0];assert len(roots)==1;s=roots[0]
            assert a.integral(s)and a.integral(a.inv(s));z=a.divide(z,s);x=a.divide(x,s)
            if a.sign(z)<0:z=a.neg(z);x=a.neg(x)
            orbits.append((z,x))
        saved=d['targets'][i]
        assert a.dec(saved['q'])==q and orbits==[a.vec(v['representative'])for v in saved['orbits']]
        for j,u in enumerate(orbits):
            n=0
            while a.sign(u[1])<=0:u=emul(u,eps);n+=1;assert abs(n)<50
            while a.sign(emul(u,back)[1])>0:u=emul(u,back);n-=1;assert abs(n)<50
            previous=emul(u,back);assert a.sign(previous[1])<=0
            first=n;powers=[]
            while a.sign(a.sub(param(u),upper))<=0:
                r=(u[1],*(a.Z for _ in range(8)),u[0]);why=reason(r,q)
                visited.append((i,j,n,r,why));powers.append(n)
                if why is None:valid.append((i,j,n,r))
                u=emul(u,eps);n+=1;assert n-first<=200
            interval=saved['orbits'][j]['interval']
            assert(interval['first_power'],interval['last_power'],interval['visited_powers'])==(first,n-1,powers)
            assert a.vec(interval['previous'])==previous and a.vec(interval['next'])==u
    actual=sorted((v['profile'],v['orbit'],v['power'],a.vec(v['root']),v['reason'])for v in d['visited_points'])
    assert sorted(visited)==actual and len(visited)==2 and len(valid)==1
    root=valid[0][-1];assert root==a.vec(d['chosen_root'])and P(root,root)==a.dec(d['chosen_norm'])
    walls=old.matrix(d['all_walls']);initial=list(basis[1:9])+[tuple(map(a.neg,basis[9]))]
    assert walls==tuple(initial+[root]);point=a.vec(d['endpoint']);center=a.vec(d['strict_interior']);escape=a.vec(d['strict_escape'])
    assert point==(a.O,*(a.Z for _ in range(8)),upper)
    assert P(point,point)==a.dec(d['endpoint_norm'])and a.sign(P(point,point))<0
    assert all(P(point,r)==a.Z for r in walls[:8]+(root,))and a.sign(P(point,walls[8]))<0
    assert a.sign(P((a.O,*(a.Z for _ in range(9))),root))<0
    for i,r in enumerate(walls):
        assert reason(r,P(r,r))is None
        for s in walls[:i]:
            p=P(r,s);assert a.sign(p)<=0
            lam=a.divide(a.scale(a.sq(p),4),a.mul(P(r,r),P(s,s)))
            assert lam in a.ALLOWED or a.sign(a.sub(lam,a.scale(a.O,4)))>=0
    assert a.sign(P(center,center))<0 and a.sign(P(escape,escape))>0
    assert P(escape,escape)==a.dec(d['escape_norm'])and a.sign(P(center,escape))<0
    assert all(a.sign(P(v,r))<0 for v in(center,escape)for r in walls)
    # Independent integer closure in the eight-root E8 basis.
    eg=[[int(P(r,s)[0])for s in walls[:8]]for r in walls[:8]]
    assert all(P(r,s)[1:]==(0,0,0)for r in walls[:8]for s in walls[:8])
    seen={tuple(int(i==j)for j in range(8))for i in range(8)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(7,-1,-1):
            w=list(v);w[j]-=sum(v[k]*eg[k][j]for k in range(8));w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=240
    assert len(seen)==240 and d['finite_mirror_counts']==[121,121]
    out={'status':'PASS_INDEPENDENT_E8_FIXED_COMPLETE_FIRST_EDGE',
        'sources':{**{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(source,Path(__file__).resolve(),Path(old.__file__).resolve(),Path(a.__file__).resolve())},**d['sources']},
        'safe_profiles':10,'complete_interval_points':2,'primitive_first_mirrors':1,
        'chosen_norm':a.enc(P(root,root)),'full_corner_mirror_counts':[121,121],
        'facets':10,'pair_gates':45,'strict_escape_norm':a.enc(P(escape,escape)),
        'scope':'Separate Fraction arithmetic replays every relative-norm orbit and exact first interval, primitive reflection integrality, full E8 root closure, finite endpoint and ten-wall geometry. PARI certified absolute norm/unit completeness and the finite reflection classification are shared explicit mathematical dependencies. Indefinite host conjugacy remains open.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'profiles',10,'points',2,'mirrors',1,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a0=p.parse_args();run(a0.verify)
