#!/usr/bin/env python3
"""Complete scalar edge out of the certified E7+A1+A1 corner."""
from fresh_h9_e8_first_20261013 import *
TAG='fresh_h9_e8_next_20261013';OUT=ROOT/'results'/f'{TAG}.json'

def run(verify=False):
    axispath=ROOT/'results'/f'{PREFIX}.json';e7path=ROOT/'results'/f'{PREFIX}_e7.json'
    auditpath=ROOT/'results'/f'{PREFIX}_e7_independent.json'
    axis=json.loads(axispath.read_text());e7=json.loads(e7path.read_text());audit=json.loads(auditpath.read_text())
    assert audit['status']=='PASS_INDEPENDENT_COMPLETE_E7_A1_INITIAL_EDGE'
    for obj in(axis,e7,audit):
        for p,h in obj['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    basis=[tuple(map(K.from_json,r))for r in axis['lattice_basis']];BI=inverse(tuple(zip(*basis)))
    omega=tuple(map(K.from_json,e7['fundamental_weight']));zz,xx=map(K.from_json,e7['relative_unit'])
    def boost(v,sgn=-1):
        p=pair(v,G,omega);c=(zz-1)*p/2+sgn*BETA*xx*v[0]/SQRT2
        return(zz*v[0]+sgn*xx*p/SQRT2,*(v[j]+c*omega[j]for j in range(1,10)))
    standard=[tuple(O if i==j else Z for i in range(10))for j in range(10)]
    A=tuple(zip(*(boost(v)for v in standard)));AI=tuple(zip(*(boost(v,1)for v in standard)))
    assert all(boost(boost(v),1)==v for v in standard)
    assert all(pair(boost(r),G,boost(s))==pair(r,G,s)for r in basis for s in basis)
    assert all(boost(b,1)[0].is_integral and boost(b,1)[9].is_integral for b in basis)
    module=tuple(zip(*(mv(BI,boost(b))for b in basis)))
    failures=[[i,j]for i,row in enumerate(module)for j,x in enumerate(row)if not x.is_integral]
    assert len(failures)==3
    transformed_edge=[basis[i]for i in[1,3,4,5,6,7,8]]+[tuple(-x for x in omega)]
    edge=[boost(r)for r in transformed_edge];oldwalls=[tuple(map(K.from_json,r))for r in e7['combined_walls']]
    assert edge==[oldwalls[i]for i in[0,2,3,4,5,6,7,10]]
    start=boost(standard[0]);savedstart=tuple(map(K.from_json,e7['endpoint']))
    assert tuple(x/start[0]for x in start)==savedstart
    C=tuple(tuple(pair(r,G,s)for s in transformed_edge)for r in transformed_edge);ci=inverse(C)
    def integral(r):return all(x.is_integral for x in mv(BI,r))
    def reason(r,q):
        assert pair(r,G,r)==q
        if not integral(r):return'NONINTEGRAL'
        if any(integral(tuple(x/p for x in r))for p in(T,K(3))):return'IMPRIMITIVE'
        if not all((2*pair(r,G,b)/q).is_integral for b in basis):return'NONREFLECTIVE'
        return None
    profiles=[];nonintegral=[]
    for q in POOL:
        alph=[]
        for r in transformed_edge:
            aa={Z}
            for lam in ALLOWED:
                p=square_root(q*pair(r,G,r)*lam/4)
                if p is not None and p.is_integral and(2*p/q).is_integral:aa.add(p)
            alph.append(sorted(aa))
        def walk(pre):
            n=len(pre)
            if n and pair(pre,tuple(row[:n]for row in ci[:n]),pre)>=q:return
            if n==8:
                co=mv(ci,tuple(-v for v in pre));r0=tuple(sum((c*r[j]for c,r in zip(co,transformed_edge)),Z)for j in range(10))
                N=q-pair(r0,G,r0)
                if any(N.sign(e)<=0 for e in range(4)):return
                assert r0[0]==r0[9]==0
                row={'q':q.to_json(),'magnitudes':enc(pre),'fixed':enc(r0),'target':N.to_json()}
                if not N.is_integral:nonintegral.append(row);return
                row['absolute_norm']=int(N.norm());profiles.append(row);return
            for p in alph[n]:walk((*pre,p))
        walk(())
    rawp=ROOT/'results'/f'{PREFIX}_raw.txt';raw=rawp.read_text();absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('NORM '):
            _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
        elif line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(tuple(map(K.from_json,json.loads(s))))
    assert{r['absolute_norm']for r in profiles}<=set(absolute)
    assert all(len(v)==counts[n]for n,v in absolute.items())
    eps=tuple(map(K.from_json,axis['relative_unit']));back=eps[0],-eps[1]
    assert eps[0]>1 and eps[1]>0 and eps[0]**2-BETA*eps[1]**2==1
    def mul(u,v):return u[0]*v[0]+BETA*u[1]*v[1],u[0]*v[1]+u[1]*v[0]
    def root(r0,u):return boost((u[1],*r0[1:9],u[0]))
    def param(u):return BETA*u[1]/u[0]
    targets={};starts=[];proposals=[]
    for row in profiles:
        N=K.from_json(row['target'])
        if N in targets:continue
        orbits=[];bad=[]
        for z,x in absolute[int(N.norm())]:
            ratio=(z*z-BETA*x*x)/N
            if not ratio.is_integral or not ratio.inverse().is_integral:
                bad.append({'representative':enc((z,x)),'ratio':ratio.to_json()});continue
            assert all(ratio.sign(e)>0 for e in range(4))
            t=square_root(ratio);assert t is not None and t.is_integral and t.inverse().is_integral
            z/=t;x/=t
            if z<0:z=-z;x=-x
            assert z*z-BETA*x*x==N
            orbits.append({'representative':enc((z,x)),'unit_correction':t.to_json()})
        targets[N]={'target':N.to_json(),'orbits':orbits,'nonunit_rejections':bad}
    for i,row in enumerate(profiles):
        q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        for j,orb in enumerate(targets[N]['orbits']):
            u=tuple(map(K.from_json,orb['representative']));power=0
            while u[1]<=0:u=mul(u,eps);power+=1;assert abs(power)<50
            while mul(u,back)[1]>0:u=mul(u,back);power-=1;assert abs(power)<50
            starts.append((i,j,power,u))
            for step in range(100):
                r=root(r0,u)
                if reason(r,q)is None:proposals.append((param(u),i,j,power,r));break
                u=mul(u,eps);power+=1
    assert proposals;proposals.sort(key=lambda x:x[0]);upper=proposals[0][0];visited=[];valid=[]
    for i,j,power,u in starts:
        row=profiles[i];q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        previous=mul(u,back);first=power;pts=[];assert previous[1]<=0
        while param(u)<=upper:
            r=root(r0,u);why=reason(r,q);rec={'profile':i,'orbit':j,'power':power,'root':enc(r),'q':q.to_json(),'parameter':param(u).to_json(),'reason':why}
            visited.append(rec);pts.append(power)
            if why is None:valid.append(rec)
            u=mul(u,eps);power+=1;assert power-first<=200
        targets[N]['orbits'][j].setdefault('intervals',[]).append({'profile':i,'first_power':first,'last_power':power-1,'previous':enc(previous),'next':enc(u),'visited_powers':pts})
    assert valid and all(K.from_json(v['parameter'])==upper for v in valid)
    choices=[]
    for v in valid:
        r=tuple(map(K.from_json,v['root']));frame=edge+[r];fi=inverse(tuple(tuple(pair(s,G,t)for t in frame)for s in frame));cone=[]
        for w in valid:
            s=tuple(map(K.from_json,w['root']));co=mv(fi,tuple(pair(s,G,t)for t in frame))
            if any(c<0 for c in co):break
            assert tuple(sum((c*t[i]for c,t in zip(co,frame)),Z)for i in range(10))==s;cone.append(enc(co))
        else:choices.append((v,cone))
    assert len(choices)==1;chosen,cone=choices[0];r=tuple(map(K.from_json,chosen['root']))
    point=boost((O,*(Z for _ in range(8)),upper));point=tuple(x/point[0]for x in point)
    assert pair(point,G,point)<0 and all(pair(point,G,s)==0 for s in edge+[r])
    assert pair(point,G,oldwalls[8])<0 and pair(start,G,r)<0
    walls=oldwalls+[r]
    for i,s in enumerate(walls):
        assert reason(s,pair(s,G,s))is None
        for t in walls[:i]:
            p=pair(s,G,t);assert p<=0,(i,walls.index(t),p)
            lam=4*p*p/(pair(s,G,s)*pair(t,G,t));assert lam in ALLOWED or lam>=4
    # Complete spherical mirror closure in exact K coordinates; no integral Cartan assumption.
    frame=edge+[r];fg=tuple(tuple(pair(s,G,t)for t in frame)for s in frame)
    seen={tuple(O if i==j else Z for j in range(9))for i in range(9)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(9):
            c=2*sum((v[k]*fg[k][j]for k in range(9)),Z)/fg[j][j]
            w=list(v);w[j]-=c;w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=500
    assert all(tuple(-x for x in v)in seen for v in seen)
    center=tuple(map(K.from_json,axis['strict_interior']));escape=tuple(map(K.from_json,axis['strict_escape']))
    assert pair(center,G,center)<0 and pair(escape,G,escape)>0 and pair(center,G,escape)<0
    assert all(pair(v,G,s)<0 for v in(center,escape)for s in walls)
    out={'status':'PASS_COMPLETE_SCALAR_EDGE_FROM_E7_A1_A1','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(axispath,e7path,auditpath,rawp,Path(__file__).resolve(),ROOT/'scripts/fresh_h9_glued_host_20261013_arithmetic.py',Path(square_root.__code__.co_filename).resolve())},
         'coefficient':BETA.to_json(),'coordinate_change':menc(A),'coordinate_inverse':menc(AI),'coordinate_module_nonintegral_entries':failures,
         'integral_time_scalar_pullbacks':enc(tuple(boost(b,1)[i]for b in basis for i in(0,9))),
         'start':enc(savedstart),'released_wall_index':8,'edge_wall_indices':[0,2,3,4,5,6,7,10],'edge':menc(edge),'transformed_edge':menc(transformed_edge),'edge_gram':menc(C),
         'profiles':profiles,'nonintegral_target_profiles':nonintegral,'relative_unit':enc(eps),'targets':list(targets.values()),'proposal_power_limit':100,
         'upper_parameter':upper.to_json(),'visited_points':visited,'primitive_reflective_survivors':valid,'chosen_root':enc(r),'chosen_norm':pair(r,G,r).to_json(),
         'simultaneous_root_cone_coefficients':cone,'endpoint':enc(point),'endpoint_norm':pair(point,G,point).to_json(),'endpoint_incident':menc(frame),'endpoint_gram':menc(fg),
         'finite_signed_root_closure':len(seen),'combined_walls':menc(walls),'strict_interior':enc(center),'strict_escape':enc(escape),'escape_norm':pair(escape,G,escape).to_json(),
         'scope':'Complete finite first contact on the scalar edge from the certified E7+A1+A1 corner. All ten safe primitive norms, every integral finite-angle profile and every relative norm-one orbit interval are included. The rational boost is not an E8 lattice automorphism; every candidate is pulled back and tested in the actual lattice. Certified PARI norm/unit completeness is inherited explicitly. No full-volume or conjugacy claim.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'profiles',len(profiles),'points',len(visited),'mirrors',len(valid),'closure',len(seen),'q',pair(r,G,r).to_json(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
