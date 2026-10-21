#!/usr/bin/env python3
"""Independent complete profile/orbit replay of the initial E7+A1 edge."""
from fresh_h9_e8_first_20261013_independent import *
TAG=PREFIX+'_e7';OUT=ROOT/'results'/f'{TAG}_independent.json'

def run(verify=False):
    source=ROOT/'results'/f'{TAG}.json';d=json.loads(source.read_text())
    axispath=ROOT/'results'/f'{PREFIX}.json';axis=json.loads(axispath.read_text())
    auditpath=ROOT/'results'/f'{PREFIX}_independent.json';audit=json.loads(auditpath.read_text())
    assert audit['status']=='PASS_INDEPENDENT_E8_FIXED_COMPLETE_FIRST_EDGE'
    for obj in(d,audit):
        for p,h in obj['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    basis=old.matrix(axis['lattice_basis']);B=tuple(zip(*basis));BI=a.minverse(B)
    initial=list(basis[1:9])+[tuple(map(a.neg,basis[9]))];edge=initial[:1]+initial[2:]
    assert old.matrix(d['edge'])==tuple(edge)
    C=tuple(tuple(P(r,s)for s in edge)for r in edge);ci=a.minverse(C);assert C==old.matrix(d['edge_gram'])
    omega=a.vec(d['fundamental_weight']);assert P(omega,omega)==a.scale(a.O,2)
    assert all(P(omega,r)==a.scale(a.O,int(i==1))for i,r in enumerate(initial[:8]))
    def integral(r):return all(map(a.integral,a.mv(BI,r)))
    def reason(r,q):
        assert P(r,r)==q
        if not integral(r):return'NONINTEGRAL'
        if any(integral(tuple(a.divide(x,p)for x in r))for p in(a.T,a.scale(a.O,3))):return'IMPRIMITIVE'
        if not all(a.integral(a.divide(a.scale(P(r,b),2),q))for b in basis):return'NONREFLECTIVE'
        return None
    found={}
    for q in reversed(POOL):
        alph=[]
        for r in edge:
            aa={a.Z}
            for lam in a.ALLOWED:
                for p in old.quartic_roots(a.scale(a.mul(a.mul(q,P(r,r)),lam),a.Q(1,4))):
                    if a.sign(p)>=0 and a.integral(p)and a.integral(a.divide(a.scale(p,2),q)):aa.add(p)
            alph.append(aa)
        def walk(pre):
            n=len(pre);ids=tuple(range(7,7-n,-1))
            if n:
                cc=tuple(tuple(ci[i][j]for j in ids)for i in ids)
                if a.sign(a.sub(q,a.dot(pre,a.mv(cc,pre))))<=0:return
            if n==8:
                mags=tuple(reversed(pre));co=a.mv(ci,tuple(map(a.neg,mags)))
                r0=tuple(old.summ(a.mul(c,r[j])for c,r in zip(co,edge))for j in range(10))
                N=a.scale(a.sub(q,P(r0,r0)),2)
                if any(a.sign(N,e)<=0 for e in range(4)):return
                assert a.integral(N)and P(r0,omega)==a.Z and r0[0]==a.Z
                found[q,mags]=(r0,N);return
            for p in alph[7-n]:walk((*pre,p))
        walk(())
    saved={(a.dec(v['q']),a.vec(v['magnitudes'])):(a.vec(v['fixed']),a.dec(v['target']))for v in d['profiles']}
    assert found==saved and len(found)==12
    raw=(ROOT/'results'/f'{TAG}_raw.txt').read_text();prior=(ROOT/'results'/f'{PREFIX}_raw.txt').read_text()
    prefixes=('ABSOLUTE ','CERTIFY ','BASE_CERTIFY ','CLASS ','SIGNATURE ','BASE_UNIT ','UNIT ','UNIT_NORM ')
    assert[line for line in raw.splitlines()if line.startswith(prefixes)]==[line for line in prior.splitlines()if line.startswith(prefixes)]
    assert d['base_units']==axis['base_units']and d['fundamental_units']==axis['fundamental_units']and d['relative_norm_map']==axis['relative_norm_map']
    absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('NORM '):
            _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
        elif line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(a.vec(json.loads(s)))
    assert set(absolute)=={int(old.norm(N))for r,N in found.values()}
    assert all(len(v)==counts[n]for n,v in absolute.items())and{str(n):c for n,c in counts.items()}==d['absolute_norm_counts']
    eps=a.vec(axis['relative_unit']);assert eps==a.vec(d['relative_unit']);back=eps[0],a.neg(eps[1])
    sq2=a.sub(a.sq(a.T),a.scale(a.O,2));assert a.sq(sq2)==a.scale(a.O,2)
    def emul(u,v):return a.add(a.mul(u[0],v[0]),a.mul(BETA,a.mul(u[1],v[1]))),a.add(a.mul(u[0],v[1]),a.mul(u[1],v[0]))
    def param(u):return a.divide(a.mul(BETA,u[1]),a.mul(sq2,u[0]))
    upper=a.dec(d['upper_parameter']);visited=[];valid=[]
    targets={a.dec(v['target']):v for v in d['targets']}
    for i,row in enumerate(d['profiles']):
        q,N=a.dec(row['q']),a.dec(row['target']);r0=a.vec(row['fixed']);orbits=[]
        for z,x in absolute[int(old.norm(N))]:
            nu=a.sub(a.sq(z),a.mul(BETA,a.sq(x)));assert old.norm(nu)==old.norm(N)
            ratio=a.divide(nu,N)
            if not a.integral(ratio)or not a.integral(a.inv(ratio)):continue
            assert all(a.sign(ratio,e)>0 for e in range(4))
            rr=[s for s in old.quartic_roots(ratio)if a.sign(s)>0];assert len(rr)==1;s=rr[0]
            assert a.integral(s)and a.integral(a.inv(s));z=a.divide(z,s);x=a.divide(x,s)
            if a.sign(z)<0:z=a.neg(z);x=a.neg(x)
            orbits.append((z,x))
        assert orbits==[a.vec(v['representative'])for v in targets[N]['orbits']]
        for j,u in enumerate(orbits):
            n=0
            while a.sign(u[1])<=0:u=emul(u,eps);n+=1;assert abs(n)<50
            while a.sign(emul(u,back)[1])>0:u=emul(u,back);n-=1;assert abs(n)<50
            previous=emul(u,back);assert a.sign(previous[1])<=0;first=n;powers=[]
            while a.sign(a.sub(param(u),upper))<=0:
                r=(a.divide(u[1],sq2),*(a.sub(r0[k],a.scale(a.mul(u[0],omega[k]),a.Q(1,2)))for k in range(1,10)))
                why=reason(r,q);visited.append((i,j,n,r,why));powers.append(n)
                if why is None:valid.append(r)
                u=emul(u,eps);n+=1;assert n-first<=200
            inter=next(v for v in targets[N]['orbits'][j]['intervals']if v['profile']==i)
            assert(inter['first_power'],inter['last_power'],inter['visited_powers'])==(first,n-1,powers)
            assert a.vec(inter['previous'])==previous and a.vec(inter['next'])==u
    assert sorted(visited)==sorted((v['profile'],v['orbit'],v['power'],a.vec(v['root']),v['reason'])for v in d['visited_points'])
    assert len(visited)==7 and len(valid)==1;root=valid[0];assert root==a.vec(d['chosen_root'])
    point=a.vec(d['endpoint']);assert point==(a.O,*(a.neg(a.mul(upper,x))for x in omega[1:]))
    assert P(point,point)==a.dec(d['endpoint_norm'])and a.sign(P(point,point))<0
    frame=edge+[root];assert all(P(point,s)==a.Z for s in frame)and a.sign(P(point,initial[1]))<0
    assert a.sign(P((a.O,*(a.Z for _ in range(9))),root))<0
    walls=old.matrix(d['combined_walls']);assert walls==tuple(initial+[a.vec(axis['chosen_root']),root])
    for i,r in enumerate(walls):
        assert reason(r,P(r,r))is None
        for s in walls[:i]:
            p=P(r,s);assert a.sign(p)<=0;lam=a.divide(a.scale(a.sq(p),4),a.mul(P(r,r),P(s,s)))
            assert lam in a.ALLOWED or a.sign(a.sub(lam,a.scale(a.O,4)))>=0
    center,escape=map(a.vec,(d['strict_interior'],d['strict_escape']))
    assert a.sign(P(center,center))<0 and a.sign(P(escape,escape))>0 and a.sign(P(center,escape))<0
    assert all(a.sign(P(v,r))<0 for v in(center,escape)for r in walls)
    fg=[[P(r,s)for s in frame]for r in frame]
    seen={tuple(int(i==j)for j in range(9))for i in range(9)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(8,-1,-1):
            c=a.divide(a.scale(old.summ(a.scale(fg[k][j],v[k])for k in range(9)),2),fg[j][j]);assert a.integral(c)and c[1:]==(0,0,0)
            w=list(v);w[j]-=int(c[0]);w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=130
    assert len(seen)==130 and d['finite_mirrors']==65
    out={'status':'PASS_INDEPENDENT_COMPLETE_E7_A1_INITIAL_EDGE','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(source,auditpath,Path(__file__).resolve())},
         'profiles':12,'complete_interval_points':7,'primitive_first_mirrors':1,'new_corner_mirrors':65,
         'combined_facets':11,'pair_gates':55,'strict_escape_norm':a.enc(P(escape,escape)),
         'scope':'Independent Fraction enumeration regenerates all spherical profiles and norm-orbit intervals, then verifies integrality, endpoint, full finite closure, all eleven-wall angle gates and strict escape. Certified PARI completeness is shared explicitly; no full-volume or host-conjugacy conclusion.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'profiles',12,'points',7,'mirrors',1,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
