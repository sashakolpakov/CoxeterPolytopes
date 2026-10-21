#!/usr/bin/env python3
"""Separate Fraction replay of all norm/profile first-edge data in the A6+D3 host."""
import argparse,hashlib,json,math
from pathlib import Path
import fresh_h9_changed_spatial_20261013_independent as old
a=old.a;ROOT=old.ROOT;TAG='fresh_h9_a6_d3_20261013';OUT=ROOT/'results'/f'{TAG}_independent.json'
BETA=a.mul(a.ALPHA,a.add(a.scale(a.O,2),a.T));PI=a.add(a.scale(a.O,2),a.T)
G=tuple(tuple(a.neg(BETA)if i==j==0 else a.scale(a.O,2)if i==j and 1<=i<=6 else a.neg(a.O)if 1<=i<=6 and 1<=j<=6 and abs(i-j)==1 else a.O if i==j else a.Z for j in range(10))for i in range(10))
STD=[tuple(a.O if i==j else a.Z for i in range(10))for j in range(10)]
def plus(v,w):return tuple(map(a.add,v,w))
def neg(v):return tuple(map(a.neg,v))
BASIS=STD[:7]+[plus(STD[7],neg(STD[8])),plus(STD[8],neg(STD[9])),plus(STD[8],STD[9])]
BI=a.minverse(tuple(zip(*BASIS)));INITIAL=[neg(v)for v in STD[1:7]]+[neg(BASIS[7]),neg(BASIS[8]),tuple(a.scale(x,-2)for x in STD[9])]
DYADIC=old.POOL;P7=(a.add(a.sq(a.T),a.O),a.sub(a.scale(a.O,5),a.sq(a.T)));POOL=tuple(a.mul(q,a.mul(x,y))for x in(a.O,P7[0])for y in(a.O,P7[1])for q in DYADIC)
def P(r,s):return a.dot(r,a.mv(G,s))
def integral(r):return all(map(a.integral,a.mv(BI,r)))
def reason(r,q):
    assert P(r,r)==q
    if not integral(r):return'NONINTEGRAL'
    if any(integral(tuple(a.divide(x,p)for x in r))for p in(a.T,*P7)):return'IMPRIMITIVE'
    if not all(a.integral(a.divide(a.scale(P(r,b),2),q))for b in BASIS):return'NONREFLECTIVE'
    return None

def qdet(M):
    aa=[list(r)for r in M];ans=a.Q(1)
    for i in range(len(aa)):
        ks=[j for j in range(i,len(aa))if aa[j][i]]
        if not ks:return a.Q(0)
        j=ks[0];aa[i],aa[j]=aa[j],aa[i]
        if i!=j:ans=-ans
        z=aa[i][i];ans*=z
        for k in range(i+1,len(aa)):
            c=aa[k][i]/z
            for t in range(i+1,len(aa)):aa[k][t]-=c*aa[i][t]
    return ans

def replay(suffix):
    source=ROOT/'results'/f'{TAG}{suffix}.json';d=json.loads(source.read_text())
    for path,h in d['sources'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
    assert old.matrix(d['ambient_gram'])==G and old.matrix(d['lattice_basis'])==tuple(BASIS)
    assert a.vec(d['safe_norm_pool'])==POOL and len(set(POOL))==52
    assert a.mul(*P7)==a.scale(a.O,7)and all(old.norm(p)==49 for p in P7)
    assert [old.norm(q)for q in POOL]==[49**j*2**k for j in(0,1,1,2)for k in range(13)]
    gi=a.minverse(tuple(tuple(P(r,s)for s in BASIS)for r in BASIS))
    assert all(a.integral(a.scale(x,28))for row in gi for x in row)
    assert all(a.integral(a.divide(a.scale(a.O,56),q))and all(a.sign(q,e)>0 for e in range(4))for q in POOL)
    for rr in d['nonsquare_comparison_ratios'].values():assert not old.quartic_roots(a.dec(rr))
    if suffix=='_a5':
        prior=json.loads((ROOT/'results'/f'{TAG}_third.json').read_text());oldwalls=old.matrix(prior['combined_walls']);D=a.scale(BETA,42);pp=STD[0];ai=a.minverse(tuple(row[1:7]for row in G[1:7]));qq=(a.Z,*(ai[i][0]for i in range(6)),a.Z,a.Z,a.Z);ids=list(range(1,9));released=0
        assert P(qq,qq)==a.scale(a.O,a.Q(6,7))
    elif suffix=='_third':
        prior=json.loads((ROOT/'results'/f'{TAG}_second.json').read_text());oldwalls=old.matrix(prior['combined_walls']);D=BETA;point0=a.vec(prior['endpoint']);ss=[v for v in old.quartic_roots(a.divide(a.neg(P(point0,point0)),D))if a.sign(v)>0];assert len(ss)==1;pp=tuple(a.divide(x,ss[0])for x in point0);qq=STD[9]
        assert pp==a.vec(d['complement_time_vector'])and qq==a.vec(d['complement_spatial_vector']);ids=[0,1,2,3,4,5,9,10];released=8
    elif suffix:
        first=json.loads((ROOT/'results'/f'{TAG}.json').read_text());oldwalls=old.matrix(first['combined_walls']);D=a.ALPHA
        cap=oldwalls[-1];pp=(a.scale(cap[7],a.Q(1,2)),*(a.Z for _ in range(6)),a.scale(a.mul(BETA,cap[0]),a.Q(1,2)),a.Z,a.Z);qq=STD[8]
        assert pp==a.vec(d['complement_time_vector'])and qq==a.vec(d['complement_spatial_vector'])
        ids=[0,1,2,3,4,5,8,9];released=7
    else:oldwalls=tuple(INITIAL);D=BETA;pp=STD[0];qq=STD[7];ids=[0,1,2,3,4,5,7,8];released=6
    assert P(pp,qq)==a.Z
    if suffix!='_a5':assert P(pp,pp)==a.neg(D)and P(qq,qq)==a.O
    scale=2 if suffix=='_third'else 1
    if suffix!='_a5':assert all(a.integral(a.divide(a.scale(P(b,pp),-scale),D))and a.integral(a.scale(P(b,qq),scale))for b in BASIS)
    else:assert all(a.integral(a.scale(P(b,qq),7))and a.integral(b[0])for b in BASIS)
    edge=tuple(oldwalls[i]for i in ids);assert edge==old.matrix(d['edge']);ci=a.minverse(tuple(tuple(P(r,s)for s in edge)for r in edge))
    found={};rejected={}
    for q in reversed(POOL):
        alph=[]
        for r in edge:
            aa={a.Z}
            for lam in a.ALLOWED:
                for p in old.quartic_roots(a.scale(a.mul(a.mul(q,P(r,r)),lam),a.Q(1,4))):
                    if a.sign(p)>=0 and a.integral(p)and a.integral(a.divide(a.scale(p,2),q)):aa.add(p)
            alph.append(aa)
        def walk(pre):
            n=len(pre);ii=tuple(range(7,7-n,-1))
            if n and a.sign(a.sub(q,a.dot(pre,a.mv(tuple(tuple(ci[i][j]for j in ii)for i in ii),pre))))<=0:return
            if n==8:
                mags=tuple(reversed(pre));co=a.mv(ci,tuple(map(a.neg,mags)));r0=tuple(old.summ(a.mul(c,r[j])for c,r in zip(co,edge))for j in range(10));N=a.scale(a.sub(q,P(r0,r0)),42 if suffix=='_a5'else scale**2)
                if any(a.sign(N,e)<=0 for e in range(4)):return
                assert P(r0,pp)==P(r0,qq)==a.Z
                fixedids=(7,8,9)if suffix=='_a5'else tuple(range(1,7))if suffix=='_third'else(*range(1,7),9)if suffix else(*range(1,7),8,9)
                why=None
                if any(not a.integral(r0[i])for i in fixedids):why='FIXED_COORDINATE_NONINTEGRAL'
                elif q not in DYADIC and suffix!='_a5':
                    assert all(x==a.Z for x in mags)and all(x==a.Z for x in r0)
                    # All odd-prime profiles have zero pairings; time and Euclidean
                    # spatial reflection coefficients force divisibility by p7.
                    assert all(a.integral(a.divide(BETA,p))is False for p in P7)
                    why='PRIME_SEVEN_FORCES_ALL_COORDINATES_DIVISIBLE'
                elif not a.integral(N):why='NONINTEGRAL_RESIDUAL_NORM'
                if why:rejected[q,mags]=(r0,N,why)
                else:found[q,mags]=(r0,N)
                return
            for p in alph[7-n]:walk((*pre,p))
        walk(())
    assert found=={(a.dec(r['q']),a.vec(r['magnitudes'])):(a.vec(r['fixed']),a.dec(r['target']))for r in d['profiles']}
    assert rejected=={(a.dec(r['q']),a.vec(r['magnitudes'])):(a.vec(r['fixed']),a.dec(r['target']),r['reason'])for r in d['rejected_profiles']}
    raw=(ROOT/'results'/f'{TAG}{suffix}_raw.txt').read_text()
    for st in('CERTIFY 1\n','BASE_CERTIFY 1\n','SIGNATURE [2, 3]\n'):assert st in raw
    assert ('CLASS 8\n'if suffix=='_a5'else'CLASS 1\n')in raw
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
    assert[e for n,e in norms]==d['relative_norm_map']
    def epower(u,n):
        if n<0:
            nn=a.sub(a.sq(u[0]),a.mul(D,a.sq(u[1])));u=a.divide(u[0],nn),a.divide(a.neg(u[1]),nn);n=-n
        out=(a.O,a.Z)
        for _ in range(n):out=emul(out,u)
        return out
    def emul(u,v):return a.add(a.mul(u[0],v[0]),a.mul(D,a.mul(u[1],v[1]))),a.add(a.mul(u[0],v[1]),a.mul(u[1],v[0]))
    for(z,x),(nu,e)in zip(units,norms):
        assert a.sub(a.sq(z),a.mul(D,a.sq(x)))==nu
        ans=a.scale(a.O,(-1)**e[-1])
        for b,n in zip(bases,e):
            if n<0:b=a.inv(b);n=-n
            for _ in range(n):ans=a.mul(ans,b)
        assert ans==nu
    mat=[[row[1][i]for row in norms]for i in range(3)];ker=[int((-1)**j*qdet([[a.Q(mat[i][k])for k in range(4)if k!=j]for i in range(3)]))for j in range(4)]
    gg=math.gcd(*ker);assert gg;ker=[n//gg for n in ker]
    if sum(n*row[1][3]for n,row in zip(ker,norms))%2:ker=[2*n for n in ker]
    assert ker==d['relative_unit_kernel'];eps=(a.O,a.Z)
    for u,n in zip(units,ker):eps=emul(eps,epower(u,n))
    if a.sign(eps[0])<0:eps=tuple(map(a.neg,eps))
    if a.sign(eps[1])<0:eps=eps[0],a.neg(eps[1])
    assert eps==a.vec(d['relative_unit'])and a.sub(a.sq(eps[0]),a.mul(D,a.sq(eps[1])))==a.O
    assert set(absolute)=={int(old.norm(N))for r,N in found.values()}
    assert all(len(v)==counts[n]for n,v in absolute.items())and{str(n):c for n,c in counts.items()}==d['absolute_norm_counts']
    def param(u):return a.divide(a.mul(D,u[1]),a.scale(u[0],6 if suffix=='_a5'else 1))
    upper=a.dec(d['upper_parameter']);back=eps[0],a.neg(eps[1]);targets={a.dec(t['target']):t for t in d['targets']};visited=[];valid=[]
    for i,row in enumerate(d['profiles']):
        q,N=a.dec(row['q']),a.dec(row['target']);r0=a.vec(row['fixed']);orbits=[]
        for z,x in absolute[int(old.norm(N))]:
            nu=a.sub(a.sq(z),a.mul(D,a.sq(x)));assert old.norm(nu)==old.norm(N);ratio=a.divide(nu,N)
            if not a.integral(ratio)or not a.integral(a.inv(ratio)):continue
            assert all(a.sign(ratio,e)>0 for e in range(4));rr=[s for s in old.quartic_roots(ratio)if a.sign(s)>0];assert len(rr)==1;s=rr[0]
            assert a.integral(s)and a.integral(a.inv(s));z=a.divide(z,s);x=a.divide(x,s)
            if a.sign(z)<0:z=a.neg(z);x=a.neg(x)
            orbits.append((z,x))
        assert orbits==[a.vec(v['representative'])for v in targets[N]['orbits']]
        for j,u in enumerate(orbits):
            n=0
            while a.sign(u[1])<=0:u=emul(u,eps);n+=1;assert abs(n)<50
            while a.sign(emul(u,back)[1])>0:u=emul(u,back);n-=1;assert abs(n)<50
            prev=emul(u,back);assert a.sign(prev[1])<=0;first=n;powers=[]
            while a.sign(a.sub(param(u),upper))<=0:
                r=tuple(a.add(r0[k],a.add(a.mul(u[1],pp[k]),a.scale(a.mul(u[0],qq[k]),a.Q(1,6))))for k in range(10))if suffix=='_a5'else tuple(a.add(r0[k],a.scale(a.add(a.mul(u[1],pp[k]),a.mul(u[0],qq[k])),a.Q(1,scale)))for k in range(10));why=reason(r,q);visited.append((i,j,n,r,why));powers.append(n)
                if why is None:valid.append(r)
                u=emul(u,eps);n+=1;assert n-first<100
            inter=next(v for v in targets[N]['orbits'][j]['intervals']if v['profile']==i)
            assert(inter['first_power'],inter['last_power'],inter['visited_powers'])==(first,n-1,powers)
            assert a.vec(inter['previous'])==prev and a.vec(inter['next'])==u
    assert sorted(visited)==sorted((v['profile'],v['orbit'],v['power'],a.vec(v['root']),v['reason'])for v in d['visited_points'])
    assert valid;root=a.vec(d['chosen_root'])
    if len(valid)>1:
        opts=[]
        for cand in valid:
            fr=(*edge,cand);fi=a.minverse(tuple(tuple(P(r,s)for s in fr)for r in fr));coefs=[]
            for rr in valid:
                cc=a.mv(fi,tuple(P(rr,s)for s in fr))
                if any(a.sign(x)<0 for x in cc):break
                assert tuple(old.summ(a.mul(c,s[i])for c,s in zip(cc,fr))for i in range(10))==rr;coefs.append(cc)
            else:opts.append((cand,coefs))
        assert len(opts)==1 and opts[0][0]==root and tuple(opts[0][1])==old.matrix(d['simultaneous_root_cone_coefficients'])
    else:assert valid[0]==root
    walls=(*oldwalls,root);assert walls==old.matrix(d['combined_walls'])
    point=tuple(a.add(pp[i],a.mul(upper,qq[i]))for i in range(10));point=tuple(a.divide(x,point[0])for x in point);assert point==a.vec(d['endpoint'])and a.sign(P(point,point))<0
    frame=(*edge,root);assert all(P(point,r)==a.Z for r in frame)and a.sign(P(point,oldwalls[released]))<0
    assert all(a.sign(P(point,r))<=0 for r in oldwalls)and a.sign(P(pp,root))<0
    for i,r in enumerate(walls):
        assert reason(r,P(r,r))is None
        for s in walls[:i]:
            p=P(r,s);assert a.sign(p)<=0;lam=a.divide(a.scale(a.sq(p),4),a.mul(P(r,r),P(s,s)));assert lam in a.ALLOWED or a.sign(a.sub(lam,a.scale(a.O,4)))>=0
    center,escape=map(a.vec,(d['strict_interior'],d['strict_escape']));assert a.sign(P(center,center))<0 and a.sign(P(escape,escape))>0 and a.sign(P(center,escape))<0
    assert P(escape,escape)==a.dec(d['escape_norm'])and all(a.sign(P(v,r))<0 for v in(center,escape)for r in walls)
    fg=tuple(tuple(P(r,s)for s in frame)for r in frame);seen={tuple(a.O if i==j else a.Z for i in range(9))for j in range(9)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(8,-1,-1):
            c=a.divide(a.scale(old.summ(a.mul(v[k],fg[k][j])for k in range(9)),2),fg[j][j]);w=list(v);w[j]=a.sub(w[j],c);w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=500
    assert len(seen)//2==d['endpoint_mirrors']and all(tuple(map(a.neg,v))in seen for v in seen)
    return {'edge':suffix or 'first','source':str(source.relative_to(ROOT)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'safe_norm_classes':52,'profiles':len(found),'rejected_profiles':len(rejected),'complete_interval_points':len(visited),'primitive_first_mirrors':len(valid),'endpoint_mirrors':len(seen)//2,'facets':len(walls),'pair_gates':len(walls)*(len(walls)-1)//2,'strict_escape_norm':a.enc(P(escape,escape))}

def run(verify=False):
    rows=[replay(s)for s in('','_second','_third','_a5')]
    out={'status':'PASS_INDEPENDENT_A6_D3_FOUR_COMPLETE_FIRST_EDGES','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(Path(__file__).resolve(),Path(old.__file__).resolve(),Path(a.__file__).resolve())},'edges':rows,'scope':'Separate Fraction arithmetic regenerates all 52-norm spherical profiles in reverse order, local rejections, full relative-unit kernels, corrected norm orbits and every first-segment representative. It checks the actual A6+D3 module, all wall angles, finite endpoint root closures and exact strict escapes. PARI certified absolute norm/unit completeness and finite reflection classification remain shared dependencies. Full-host reflectivity remains open.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],rows,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
