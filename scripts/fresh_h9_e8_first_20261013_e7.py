#!/usr/bin/env python3
"""Complete second initial edge: the E7+A1 face in the E8+1 host."""
from fresh_h9_e8_first_20261013 import *
TAG=PREFIX+'_e7';OUT=ROOT/'results'/f'{TAG}.json'

def run(verify=False):
    axis=json.loads((ROOT/'results'/f'{PREFIX}.json').read_text())
    basis=[tuple(map(K.from_json,r))for r in axis['lattice_basis']];e8=basis[1:9]
    B=tuple(zip(*basis));BI=inverse(B);initial=e8+[tuple(-x for x in basis[9])]
    edge=initial[:1]+initial[2:];C=tuple(tuple(pair(r,G,s)for s in edge)for r in edge);ci=inverse(C)
    ec=tuple(tuple(pair(r,G,s)for s in e8)for r in e8);ei=inverse(ec)
    omega=tuple(sum((ei[1][i]*r[j]for i,r in enumerate(e8)),Z)for j in range(10))
    assert pair(omega,G,omega)==2 and all(pair(omega,G,r)==int(i==1)for i,r in enumerate(e8))
    def integral(r):return all(x.is_integral for x in mv(BI,r))
    def reason(r,q):
        assert pair(r,G,r)==q
        if not integral(r):return'NONINTEGRAL'
        if any(integral(tuple(x/p for x in r))for p in(T,K(3))):return'IMPRIMITIVE'
        if not all((2*pair(r,G,b)/q).is_integral for b in basis):return'NONREFLECTIVE'
        return None
    profiles=[]
    for q in POOL:
        alph=[]
        for r in edge:
            aa={Z}
            for lam in ALLOWED:
                p=square_root(q*pair(r,G,r)*lam/4)
                if p is not None and p.is_integral and(2*p/q).is_integral:aa.add(p)
            alph.append(sorted(aa))
        def walk(pre):
            n=len(pre)
            if n and pair(pre,tuple(row[:n]for row in ci[:n]),pre)>=q:return
            if n==8:
                co=mv(ci,tuple(-v for v in pre));r0=tuple(sum((c*r[j]for c,r in zip(co,edge)),Z)for j in range(10))
                N=2*(q-pair(r0,G,r0))
                if any(N.sign(e)<=0 for e in range(4)):return
                assert N.is_integral and pair(r0,G,omega)==0 and r0[0]==0
                profiles.append({'q':q.to_json(),'magnitudes':enc(pre),'fixed':enc(r0),'target':N.to_json(),'absolute_norm':int(N.norm())});return
            for p in alph[n]:walk((*pre,p))
        walk(())
    gp=ROOT/'scripts'/f'{TAG}.gp';rawp=ROOT/'results'/f'{TAG}_raw.txt';ns=sorted({r['absolute_norm']for r in profiles})
    if not verify:
        text=(ROOT/'scripts'/f'{PREFIX}.gp').read_text();line=next(s for s in text.splitlines()if s.startswith('ns='))
        gp.write_text(text.replace(line,'ns='+str(ns)+';'))
        with rawp.open('w')as f:
            p=subprocess.run(['gp','-q',str(gp)],stdout=f,stderr=subprocess.PIPE,text=True,timeout=90);assert p.returncode==0,p.stderr
    raw=rawp.read_text()
    for s in('ABSOLUTE x^8 + 12*x^6 + 18*x^4 - 108*x^2 - 81\n','CERTIFY 1\n','BASE_CERTIFY 1\n','CLASS 1\n','SIGNATURE [2, 3]\n'):assert s in raw
    bases=[];units=[];norms=[];absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('BASE_UNIT '):bases.append(K.from_json(json.loads(line.split(' ',2)[2])))
        elif line.startswith('UNIT '):units.append(tuple(map(K.from_json,json.loads(line.split(' ',2)[2]))))
        elif line.startswith('UNIT_NORM '):
            n,e=line.split(' ',2)[2].split('] ',1);norms.append((K.from_json(json.loads(n+']')),json.loads(e)))
        elif line.startswith('NORM '):
            _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
        elif line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(tuple(map(K.from_json,json.loads(s))))
    assert set(absolute)==set(ns)and all(len(v)==counts[n]for n,v in absolute.items())
    expected=[[0,2,0,0],[-2,0,0,0],[0,0,2,0],[0,0,0,0]]
    assert[e for n,e in norms]==expected and len(units)==4 and len(bases)==3
    for(z,x),(n,e)in zip(units,norms):assert z*z-BETA*x*x==n==K((-1)**e[3])*bases[0]**e[0]*bases[1]**e[1]*bases[2]**e[2]
    eps=units[3]
    if eps[0]<0:eps=tuple(-x for x in eps)
    if eps[1]<0:eps=eps[0],-eps[1]
    assert eps[0]>1 and eps[1]>0 and eps[0]**2-BETA*eps[1]**2==1
    def mul(u,v):return u[0]*v[0]+BETA*u[1]*v[1],u[0]*v[1]+u[1]*v[0]
    def root(r0,u):return(u[1]/SQRT2,*(r0[j]-u[0]*omega[j]/2 for j in range(1,10)))
    def param(u):return BETA*u[1]/(SQRT2*u[0])
    back=eps[0],-eps[1];targets={};starts=[];proposals=[]
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
    assert proposals
    proposals.sort(key=lambda x:x[0]);upper=proposals[0][0];visited=[];valid=[]
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
    assert len(choices)==1;chosen,cone=choices[0];r=tuple(map(K.from_json,chosen['root']));point=(O,*(-upper*omega[j]for j in range(1,10)))
    assert pair(point,G,point)<0 and all(pair(point,G,s)==0 for s in edge+[r])
    assert pair(point,G,initial[1])<0 and pair((O,*(Z for _ in range(9))),G,r)<0
    axisroot=tuple(map(K.from_json,axis['chosen_root']));walls=initial+[axisroot,r]
    for i,s in enumerate(walls):
        assert reason(s,pair(s,G,s))is None
        for t in walls[:i]:
            p=pair(s,G,t);assert p<=0
            lam=4*p*p/(pair(s,G,s)*pair(t,G,t));assert lam in ALLOWED or lam>=4
    center=tuple(map(K.from_json,axis['strict_interior']));escape=tuple(map(K.from_json,axis['strict_escape']))
    assert pair(center,G,center)<0 and pair(escape,G,escape)>0 and pair(center,G,escape)<0
    assert all(pair(v,G,s)<0 for v in(center,escape)for s in walls)
    frame=edge+[r];fg=tuple(tuple(pair(s,G,t)for t in frame)for s in frame)
    seen={tuple(int(i==j)for j in range(9))for i in range(9)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(9):
            c=2*sum((v[k]*fg[k][j]for k in range(9)),Z)/fg[j][j]
            assert c.is_integral and c.coeffs[1:]==(0,0,0)
            w=list(v);w[j]-=int(c.a);w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=130
    assert len(seen)==130 and all(tuple(-x for x in v)in seen for v in seen)
    out={'status':'PASS_COMPLETE_E7_A1_INITIAL_EDGE','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(TRIAGE,ROOT/'results'/f'{PREFIX}.json',gp,rawp,Path(__file__).resolve())},
         'coefficient':BETA.to_json(),'edge':menc(edge),'edge_gram':menc(C),'released_root':enc(initial[1]),'fundamental_weight':enc(omega),
         'profiles':profiles,'base_units':enc(bases),'fundamental_units':[enc(u)for u in units],'relative_norm_map':expected,'relative_unit':enc(eps),
         'absolute_norm_counts':{str(n):c for n,c in counts.items()},'targets':list(targets.values()),'proposal_power_limit':100,
         'upper_parameter':upper.to_json(),'visited_points':visited,'primitive_reflective_survivors':valid,'chosen_root':enc(r),'chosen_norm':pair(r,G,r).to_json(),
         'simultaneous_root_cone_coefficients':cone,'endpoint':enc(point),'endpoint_norm':pair(point,G,point).to_json(),'endpoint_incident':menc(edge+[r]),'combined_walls':menc(walls),
         'endpoint_type':'E7+A1+A1','finite_mirrors':65,'strict_interior':enc(center),'strict_escape':enc(escape),'escape_norm':pair(escape,G,escape).to_json(),
         'scope':'Complete first mirror on the initial E7+A1 edge. All ten safe primitive norms, every integral finite-angle spherical profile and every relative norm-one orbit interval are included. The search uses no old-wall compatibility pruning before firstness. All eleven combined walls pass exact Coxeter pair gates. Finite volume and global host conjugacy are not established.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'profiles',len(profiles),'points',len(visited),'mirrors',len(valid),'q',pair(r,G,r).to_json(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
