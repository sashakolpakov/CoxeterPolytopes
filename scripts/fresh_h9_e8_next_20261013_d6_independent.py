#!/usr/bin/env python3
"""Independent Fraction complete replay of the E7-corner D6+A1+A1 step."""
from fresh_h9_e8_first_20261013_independent import *
TAG='fresh_h9_e8_next_20261013_d6';OUT=ROOT/'results'/f'{TAG}_independent.json'

def run(verify=False):
    source=ROOT/'results'/f'{TAG}.json';d=json.loads(source.read_text())
    ap=ROOT/'results'/f'{PREFIX}.json';ep=ROOT/'results'/f'{PREFIX}_e7.json';axis=json.loads(ap.read_text());e7=json.loads(ep.read_text())
    prevpath=ROOT/'results/fresh_h9_e8_next_20261013.json';prev=json.loads(prevpath.read_text())
    prevcheck=ROOT/'results/fresh_h9_e8_next_20261013_independent.json';pc=json.loads(prevcheck.read_text());assert pc['status']=='PASS_INDEPENDENT_COMPLETE_E7_CORNER_SCALAR_EDGE'
    prior=ROOT/'results'/f'{PREFIX}_e7_independent.json';audit=json.loads(prior.read_text())
    assert audit['status']=='PASS_INDEPENDENT_COMPLETE_E7_A1_INITIAL_EDGE'
    for obj in(d,audit,prev,pc):
        for p,h in obj['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    basis=old.matrix(axis['lattice_basis']);BI=a.minverse(tuple(zip(*basis)))
    omega=a.vec(e7['fundamental_weight']);zz,xx=a.vec(e7['relative_unit']);sq2=a.sub(a.sq(a.T),a.scale(a.O,2))
    def boost(v,sgn=-1):
        p=P(v,omega);c=a.add(a.scale(a.mul(a.sub(zz,a.O),p),a.Q(1,2)),a.scale(a.divide(a.mul(BETA,a.mul(xx,v[0])),sq2),sgn))
        return(a.add(a.mul(zz,v[0]),a.scale(a.divide(a.mul(xx,p),sq2),sgn)),*(a.add(v[j],a.mul(c,omega[j]))for j in range(1,10)))
    standard=[tuple(a.O if i==j else a.Z for i in range(10))for j in range(10)]
    assert all(boost(boost(v),1)==v for v in standard)
    assert all(P(boost(r),boost(s))==P(r,s)for r in basis for s in basis)
    A=tuple(zip(*(boost(v)for v in standard)));AI=tuple(zip(*(boost(v,1)for v in standard)))
    assert A==old.matrix(d['coordinate_change'])and AI==old.matrix(d['coordinate_inverse'])
    pull=tuple(boost(b,1)[i]for b in basis for i in(0,9));assert all(map(a.integral,pull))and pull==a.vec(d['integral_time_scalar_pullbacks'])
    module=tuple(zip(*(a.mv(BI,boost(b))for b in basis)));bad=[[i,j]for i,row in enumerate(module)for j,x in enumerate(row)if not a.integral(x)]
    assert bad==d['coordinate_module_nonintegral_entries']and len(bad)==3
    eta=(a.Z,a.neg(a.O),a.neg(a.O),*(a.Z for _ in range(7)))
    assert eta==a.vec(d['fundamental_weight'])and P(eta,eta)==a.scale(a.O,2)and boost(eta)==eta
    assert all(a.integral(P(boost(b,1),eta))for b in basis)
    te=[basis[i]for i in[3,4,5,6,7,8]]+[tuple(map(a.neg,basis[9])),tuple(map(a.neg,omega))];edge=[boost(r)for r in te]
    oldwalls=old.matrix(prev['combined_walls']);assert tuple(edge)==tuple(oldwalls[i]for i in d['edge_wall_indices'])==old.matrix(d['edge'])
    assert tuple(te)==old.matrix(d['transformed_edge'])
    start=boost(standard[0]);assert tuple(a.divide(x,start[0])for x in start)==a.vec(d['start'])==a.vec(e7['endpoint'])
    C=tuple(tuple(P(r,s)for s in te)for r in te);ci=a.minverse(C);assert C==old.matrix(d['edge_gram'])
    def integral(r):return all(map(a.integral,a.mv(BI,r)))
    def reason(r,q):
        assert P(r,r)==q
        if not integral(r):return'NONINTEGRAL'
        if any(integral(tuple(a.divide(x,p)for x in r))for p in(a.T,a.scale(a.O,3))):return'IMPRIMITIVE'
        if not all(a.integral(a.divide(a.scale(P(r,b),2),q))for b in basis):return'NONREFLECTIVE'
        return None
    found={};nonintegral={}
    for q in reversed(POOL):
        alph=[]
        for r in te:
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
                r0=tuple(old.summ(a.mul(c,r[j])for c,r in zip(co,te))for j in range(10));N=a.scale(a.sub(q,P(r0,r0)),2)
                if any(a.sign(N,e)<=0 for e in range(4)):return
                assert r0[0]==a.Z and P(r0,eta)==a.Z
                (found if a.integral(N)else nonintegral)[q,mags]=(r0,N);return
            for p in alph[7-n]:walk((*pre,p))
        walk(())
    saved=lambda rows:{(a.dec(v['q']),a.vec(v['magnitudes'])):(a.vec(v['fixed']),a.dec(v['target']))for v in rows}
    assert found==saved(d['profiles'])and nonintegral==saved(d['nonintegral_target_profiles'])and len(found)==24
    raw=(ROOT/'results'/f'{PREFIX}_e7_raw.txt').read_text();raw2=(ROOT/'results'/f'{PREFIX}_raw.txt').read_text();absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('NORM '):
            _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
        elif line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(a.vec(json.loads(s)))
    extra={}
    for line in raw2.splitlines():
        if line.startswith('NORM '):
            _,n,c=line.split();extra[int(n)]=[]
            if int(n)not in counts:counts[int(n)]=int(c)
        elif line.startswith('SOLUTION '):
            _,n,ss=line.split(' ',2);extra[int(n)].append(a.vec(json.loads(ss)))
    for n,rows in extra.items():
        if n in absolute:assert absolute[n]==rows
        else:absolute[n]=rows
    assert{int(old.norm(N))for r,N in found.values()}<=set(absolute)
    assert all(len(v)==counts[n]for n,v in absolute.items())
    eps=a.vec(axis['relative_unit']);assert eps==a.vec(d['relative_unit']);back=eps[0],a.neg(eps[1])
    def emul(u,v):return a.add(a.mul(u[0],v[0]),a.mul(BETA,a.mul(u[1],v[1]))),a.add(a.mul(u[0],v[1]),a.mul(u[1],v[0]))
    def param(u):return a.divide(a.mul(BETA,u[1]),a.mul(sq2,u[0]))
    upper=a.dec(d['upper_parameter']);visited=[];valid=[];targets={a.dec(v['target']):v for v in d['targets']}
    for i,row in enumerate(d['profiles']):
        q,N=a.dec(row['q']),a.dec(row['target']);r0=a.vec(row['fixed']);orbits=[]
        for z,x in absolute[int(old.norm(N))]:
            nu=a.sub(a.sq(z),a.mul(BETA,a.sq(x)));assert old.norm(nu)==old.norm(N);ratio=a.divide(nu,N)
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
                r=boost((a.divide(u[1],sq2),*(a.sub(r0[k],a.scale(a.mul(u[0],eta[k]),a.Q(1,2)))for k in range(1,10))));why=reason(r,q)
                visited.append((i,j,n,r,why));powers.append(n)
                if why is None:valid.append(r)
                u=emul(u,eps);n+=1;assert n-first<=200
            inter=next(v for v in targets[N]['orbits'][j]['intervals']if v['profile']==i)
            assert(inter['first_power'],inter['last_power'],inter['visited_powers'])==(first,n-1,powers)
            assert a.vec(inter['previous'])==previous and a.vec(inter['next'])==u
    assert sorted(visited)==sorted((v['profile'],v['orbit'],v['power'],a.vec(v['root']),v['reason'])for v in d['visited_points'])
    assert len(visited)==17 and len(valid)==1;root=valid[0];assert root==a.vec(d['chosen_root'])
    point=boost((a.O,*(a.neg(a.mul(upper,eta[k]))for k in range(1,10))));point=tuple(a.divide(x,point[0])for x in point)
    assert point==a.vec(d['endpoint'])and P(point,point)==a.dec(d['endpoint_norm'])and a.sign(P(point,point))<0
    frame=edge+[root];assert all(P(point,s)==a.Z for s in frame)and a.sign(P(point,oldwalls[0]))<0 and a.sign(P(start,root))<0
    walls=old.matrix(d['combined_walls']);assert walls==(*oldwalls,root)
    for i,r in enumerate(walls):
        assert reason(r,P(r,r))is None
        for s in walls[:i]:
            p=P(r,s);assert a.sign(p)<=0;lam=a.divide(a.scale(a.sq(p),4),a.mul(P(r,r),P(s,s)))
            assert lam in a.ALLOWED or a.sign(a.sub(lam,a.scale(a.O,4)))>=0
    center,escape=map(a.vec,(d['strict_interior'],d['strict_escape']))
    assert a.sign(P(center,center))<0 and a.sign(P(escape,escape))>0 and a.sign(P(center,escape))<0
    assert all(a.sign(P(v,r))<0 for v in(center,escape)for r in walls)
    fg=[[P(r,s)for s in frame]for r in frame];seen={tuple(a.O if i==j else a.Z for j in range(9))for i in range(9)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(8,-1,-1):
            c=a.divide(a.scale(old.summ(a.mul(fg[k][j],v[k])for k in range(9)),2),fg[j][j]);w=list(v);w[j]=a.sub(w[j],c);w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=500
    assert len(seen)==d['finite_signed_root_closure']==66
    out={'status':'PASS_INDEPENDENT_COMPLETE_E7_CORNER_D6_EDGE','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(source,prior,prevcheck,Path(__file__).resolve(),Path(a.__file__).resolve(),Path(old.__file__).resolve())},
         'profiles':len(found),'nonintegral_target_profiles':len(nonintegral),'complete_interval_points':len(visited),'primitive_first_mirrors':len(valid),'new_corner_mirrors':len(seen)//2,
         'combined_facets':len(walls),'pair_gates':len(walls)*(len(walls)-1)//2,'strict_escape_norm':a.enc(P(escape,escape)),
         'scope':'Independent Fraction implementation reconstructs the rational boost and its integral time/scalar pullbacks, all profiles and all interval representatives, actual E8 module/reflection/primitivity gates, full finite closure and geometry. The original certified PARI absolute-norm and unit-completeness interface is a shared explicit dependency.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'profiles',len(found),'points',len(visited),'mirrors',len(valid),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
