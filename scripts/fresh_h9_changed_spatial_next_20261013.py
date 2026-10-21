#!/usr/bin/env python3
"""Complete initial A2+B6 edge in the quartic D9 host; new norm field."""
from fresh_h9_changed_spatial_20261013_steps import *

OLD_PREFIX=PREFIX
PREFIX='fresh_h9_changed_spatial_next_20261013'
PP=ROOT/'results'/f'{PREFIX}_profiles.json'
GP2=ROOT/'scripts'/f'{PREFIX}.gp'
RP=ROOT/'results'/f'{PREFIX}_raw.txt'
OP=ROOT/'results'/f'{PREFIX}.json'
PRIOR=ROOT/'results'/f'{OLD_PREFIX}_diagonal.json'
CHECKPOINT=ROOT/'results'/f'{OLD_PREFIX}_checkpoint.json'
D=3*ALPHA
EDGE=INITIAL[:2]+INITIAL[3:]


def mul(u,v):return u[0]*v[0]+D*u[1]*v[1],u[0]*v[1]+u[1]*v[0]


def profiles():
    C=tuple(tuple(pair(r,G,s)for s in EDGE)for r in EDGE);ci=inverse(C)
    rows=[];alphabets=[]
    for q in DYADIC:
        alph=[]
        for r in EDGE:
            aa={Z}
            for lam in ALLOWED:
                p=square_root(q*pair(r,G,r)*lam/4)
                if p is not None and p.is_integral and(2*p/q).is_integral:aa.add(p)
            alph.append(sorted(aa))
        alphabets.append({'q':q.to_json(),'magnitudes':[enc(a)for a in alph]})
        def visit(pre):
            n=len(pre)
            if n and pair(pre,tuple(row[:n]for row in ci[:n]),pre)>=q:return
            if n==8:
                coef=mv(ci,tuple(-v for v in pre))
                r0=tuple(sum((c*r[j]for c,r in zip(coef,EDGE)),Z)for j in range(10))
                assert r0[0]==0 and sum(r0[1:4],Z)==0
                if not all(v.is_integral for v in r0[4:]):return
                N=3*(q-pair(r0,G,r0))
                if any(N.sign(e)<=0 for e in range(4)):return
                assert N.is_integral and N.norm().denominator==1
                rows.append({'q':q.to_json(),'magnitudes':enc(pre),'fixed':enc(r0),
                             'target':N.to_json(),'absolute_norm':int(N.norm())});return
            for s in alph[n]:visit((*pre,s))
        visit(())
    return {'status':'COMPLETE_INITIAL_A2_B6_FINITE_PROFILES','edge':menc(EDGE),
            'gram':menc(C),'alphabets':alphabets,'profiles':rows,
            'absolute_norms':sorted({r['absolute_norm']for r in rows})}


def read_new_interface(pd):
    raw=RP.read_text()
    for s in('ABSOLUTE x^8 + 12*x^6 + 18*x^4 - 108*x^2 - 81\n',
             'CERTIFY 1\n','BASE_CERTIFY 1\n','CLASS 1\n','SIGNATURE [2, 3]\n'):assert s in raw
    bases=[];units=[];norms=[];absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('BASE_UNIT '):bases.append(K.from_json(json.loads(line.split(' ',2)[2])))
        elif line.startswith('UNIT '):units.append(tuple(map(K.from_json,json.loads(line.split(' ',2)[2]))))
        elif line.startswith('UNIT_NORM '):
            _,i,s=line.split(' ',2);n,e=s.split('] ',1)
            norms.append((K.from_json(json.loads(n+']')),json.loads(e)))
        elif line.startswith('NORM '):
            _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
        elif line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(tuple(map(K.from_json,json.loads(s))))
    assert len(bases)==3 and len(units)==4
    assert [e for n,e in norms]==[[0,2,0,0],[-2,0,0,0],[0,0,2,0],[0,0,0,0]]
    for (z,x),(n,e)in zip(units,norms):
        assert z*z-D*x*x==n
        assert n==K((-1)**e[3])*bases[0]**e[0]*bases[1]**e[1]*bases[2]**e[2]
    eps=units[3]
    if eps[0]<0:eps=(-eps[0],-eps[1])
    if eps[1]<0:eps=(eps[0],-eps[1])
    assert eps[0]>1 and eps[1]>0 and eps[0]**2-D*eps[1]**2==1
    assert set(absolute)==set(pd['absolute_norms'])and all(len(v)==counts[n]for n,v in absolute.items())
    return bases,units,norms,absolute,counts,eps


def run(verify=False,prepare=False):
    pd=profiles()
    if verify:assert pd==json.loads(PP.read_text())
    else:
        PP.write_text(json.dumps(pd,indent=2)+'\n')
        text=GP.read_text().replace('d=t-1;','d=3*(t-1);').replace('ns=vector(13,j,2^(j-1));',f'ns={pd["absolute_norms"]};')
        GP2.write_text(text)
        with RP.open('w')as f:
            job=subprocess.run(['gp','-q',str(GP2)],stdout=f,stderr=subprocess.PIPE,text=True,timeout=90)
            assert job.returncode==0,job.stderr
    base,units,norms,absolute,counts,eps=read_new_interface(pd)
    print('profiles',len(pd['profiles']),'absolute targets',counts,'unit',enc(eps),flush=True)
    if prepare:return
    back=(eps[0],-eps[1]);targets={};starts=[];proposals=[]
    for row in pd['profiles']:
        N=K.from_json(row['target'])
        if N in targets:continue
        orbits=[];bad=[]
        for z,x in absolute[int(N.norm())]:
            assert(z*z-D*x*x).norm()==N.norm()
            ratio=(z*z-D*x*x)/N
            if not ratio.is_integral or not ratio.inverse().is_integral:
                bad.append({'representative':enc((z,x)),'ratio':ratio.to_json()});continue
            assert all(ratio.sign(e)>0 for e in range(4))
            unit=square_root(ratio);assert unit is not None and unit.is_integral and unit.inverse().is_integral
            z/=unit;x/=unit
            if z<0:z=-z;x=-x
            assert z*z-D*x*x==N
            orbits.append({'representative':enc((z,x)),'unit_correction':unit.to_json()})
        targets[N]={'target':N.to_json(),'orbits':orbits,'nonunit_rejections':bad}
    def root(r0,u):
        z,x=u
        return(x,*(r0[i]+z/3 for i in range(1,4)),*r0[4:])
    for i,row in enumerate(pd['profiles']):
        q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        for j,orb in enumerate(targets[N]['orbits']):
            u=tuple(map(K.from_json,orb['representative']));power=0
            while u[1]<=0:u=mul(u,eps);power+=1;assert abs(power)<50
            while mul(u,back)[1]>0:u=mul(u,back);power-=1;assert abs(power)<50
            starts.append((i,j,power,u))
            for step in range(100):
                r=root(r0,u)
                if root_reason(r,q)is None:
                    proposals.append((ALPHA*u[1]/u[0],i,j,power,r));break
                u=mul(u,eps);power+=1
    assert proposals,'No finite upper endpoint in the stated 100-unit proposal scout.'
    proposals.sort(key=lambda r:r[0]);upper=proposals[0][0]
    print('upper proposal profile/orbit/power',proposals[0][1:4],flush=True)
    visited=[];valid=[]
    for i,j,power,u in starts:
        row=pd['profiles'][i];q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        previous=mul(u,back);first=power;pts=[]
        while ALPHA*u[1]<=upper*u[0]:
            r=root(r0,u);why=root_reason(r,q)
            rec={'profile':i,'orbit':j,'power':power,'root':enc(r),'q':q.to_json(),
                 'parameter':(ALPHA*u[1]/u[0]).to_json(),'reason':why}
            visited.append(rec);pts.append(power)
            if why is None:valid.append(rec)
            u=mul(u,eps);power+=1;assert power-first<=200
        targets[N]['orbits'][j].setdefault('intervals',[]).append({'profile':i,'first_power':first,
            'last_power':power-1,'previous':enc(previous),'next':enc(u),'visited_powers':pts})
    assert valid and all(K.from_json(r['parameter'])==upper for r in valid)
    simple=[]
    for row in valid:
        r=tuple(map(K.from_json,row['root']));frame=EDGE+[r]
        gi=inverse(tuple(tuple(pair(a,G,b)for b in frame)for a in frame));csall=[]
        for other in valid:
            s=tuple(map(K.from_json,other['root']));cs=mv(gi,tuple(pair(s,G,w)for w in frame))
            if any(c<0 for c in cs):break
            assert tuple(sum((c*w[t]for c,w in zip(cs,frame)),Z)for t in range(10))==s
            csall.append(enc(cs))
        else:simple.append((row,csall))
    assert len(simple)==1
    chosen,cone=simple[0];r=tuple(map(K.from_json,chosen['root']))
    old=json.loads(PRIOR.read_text());walls=[tuple(map(K.from_json,w))for w in old['all_walls']]
    assert all(compatible(r,w)for w in walls);walls.append(r)
    point=(O,upper,upper,upper,*(Z for _ in range(6)))
    assert pair(point,G,point)<0 and all(pair(point,G,w)<=0 for w in walls)
    assert all(pair(point,G,w)==0 for w in EDGE+[r])
    ck=json.loads(CHECKPOINT.read_text());inside=tuple(map(K.from_json,ck['strict_interior']))
    assert pair(inside,G,inside)<0 and all(pair(inside,G,w)<0 for w in walls)
    oldescape=tuple(map(K.from_json,ck['strict_escape']))
    assert pair(oldescape,G,r)>0
    escape=(O,*(K(F(v,100))for v in(72,52,33,28,20,10,5,3,1)))
    assert pair(escape,G,escape)==K(F(2537,1250))-T>0
    assert pair(escape,G,inside)<0 and all(pair(escape,G,w)<0 for w in walls)
    for i,w in enumerate(walls):
        assert root_reason(w,pair(w,G,w))is None
        assert all(compatible(w,s)for s in walls[:i])
    assert pair(r,G,r)==4 and [pair(r,G,w)for w in EDGE]==[Z,-2*O,Z,Z,Z,Z,Z,Z]
    z,x=eps
    assert r==(2*x,(2*z+2)/3,(2*z+2)/3,(2*z-4)/3,*(Z for _ in range(6)))
    result={'status':'PASS_COMPLETE_D9_INITIAL_A2_B6_EDGE_WITH_CERTIFIED_PARI_INTERFACE',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(PP,RP,GP2,PRIOR,CHECKPOINT,FIELD)},
        'new_field':'K(sqrt(3*(T-1)))','base_units':enc(base),'fundamental_units':[enc(u)for u in units],
        'relative_norm_map':[e for n,e in norms],'relative_unit':enc(eps),
        'profiles':len(pd['profiles']),'absolute_norm_counts':{str(n):c for n,c in counts.items()},
        'targets':list(targets.values()),'proposal_scout_powers':100,'upper_parameter':upper.to_json(),
        'visited_points':visited,'primitive_reflective_survivors':valid,'chosen_root':enc(r),'chosen_norm':pair(r,G,r).to_json(),
        'simultaneous_root_cone_coefficients':cone,'endpoint':enc(point),'endpoint_norm':pair(point,G,point).to_json(),
        'endpoint_incident':menc(EDGE+[r]),'all_walls':menc(walls),'wall_count':len(walls),
        'strict_interior':enc(inside),'previous_escape_pairing':pair(oldescape,G,r).to_json(),
        'endpoint_type':'B3+B6','endpoint_mirrors':45,'all_pairs_checked':len(walls)*(len(walls)-1)//2,
        'strict_escape':enc(escape),'strict_escape_norm':pair(escape,G,escape).to_json(),
        'strict_escape_pairings':enc(tuple(pair(escape,G,w)for w in walls)),
        'scope':'Complete first contact on the original A2+B6 edge. Every safe primitive norm and integral finite spherical profile is included; certified full absolute norm representatives and the full relative norm-one unit kernel supply every interval orbit. The proposal scout only supplies a finite upper endpoint. No prefix compatibility pruning occurs before first contact. Full-host reflectivity is not decided.'}
    if verify:assert result==json.loads(OP.read_text())
    else:OP.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'profiles',len(pd['profiles']),'points',len(visited),'first mirrors',len(valid),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--prepare',action='store_true')
    args=p.parse_args();run(args.verify,args.prepare)
