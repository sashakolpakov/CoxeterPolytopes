#!/usr/bin/env python3
"""Complete first contact on the initial A8 (all coordinates equal) edge."""
from fresh_h9_changed_spatial_20261013_steps import *

TAG=f'{PREFIX}_diagonal'
PP=ROOT/'results'/f'{TAG}_profiles.json';RP=ROOT/'results'/f'{TAG}_raw.txt'
GPP=ROOT/'scripts'/f'{TAG}.gp';OP=ROOT/'results'/f'{TAG}.json'


def profiles_diagonal():
    edge=INITIAL[:8];C=tuple(tuple(pair(r,G,s)for s in edge)for r in edge);ci=inverse(C)
    rows=[];alphabets=[]
    for q in DYADIC:
        alphabet={Z}
        for lam in ALLOWED:
            s=square_root(q*lam/2)
            if s is not None and s.is_integral and(2*s/q).is_integral:alphabet.add(s)
        alphabet=sorted(alphabet);alphabets.append({'q':q.to_json(),'magnitudes':enc(alphabet)})
        def visit(pre):
            n=len(pre)
            if n and pair(pre,tuple(row[:n]for row in ci[:n]),pre)>=q:return
            if n==8:
                coef=mv(ci,tuple(-v for v in pre))
                r0=tuple(sum((c*r[j]for c,r in zip(coef,edge)),Z)for j in range(10))
                assert r0[0]==0 and sum(r0[1:],Z)==0
                N=9*(q-pair(r0,G,r0))
                if any(N.sign(e)<=0 for e in range(4)):return
                assert N.is_integral
                rows.append({'q':q.to_json(),'magnitudes':enc(pre),'fixed':enc(r0),'target':N.to_json(),'absolute_norm':int(N.norm())});return
            for s in alphabet:visit((*pre,s))
        visit(())
    return {'status':'COMPLETE_A8_DIAGONAL_FINITE_PROFILES','edge':menc(edge),'gram':menc(C),
            'alphabets':alphabets,'profiles':rows,'absolute_norms':sorted({r['absolute_norm']for r in rows})}


def run(verify=False):
    pd=profiles_diagonal()
    if verify:assert pd==json.loads(PP.read_text())
    else:
        PP.write_text(json.dumps(pd,indent=2)+'\n')
        GPP.write_text(GP.read_text().replace('ns=vector(13,j,2^(j-1));',f'ns={pd["absolute_norms"]};'))
        with RP.open('w')as f:
            p=subprocess.run(['gp','-q',str(GPP)],stdout=f,stderr=subprocess.PIPE,text=True,timeout=60)
            assert p.returncode==0,p.stderr
    raw=RP.read_text();assert 'CERTIFY 1\n'in raw and'BASE_CERTIFY 1\n'in raw
    absr={};counts={}
    for line in raw.splitlines():
        if line.startswith('NORM '):
            _,n,c=line.split();absr[int(n)]=[];counts[int(n)]=int(c)
        if line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absr[int(n)].append(tuple(map(K.from_json,json.loads(s))))
    assert set(absr)==set(pd['absolute_norms'])and all(len(v)==counts[n]for n,v in absr.items())
    _,_,_,_,_,eps=read_interface();back=(eps[0],-eps[1]);targets={}
    for row in pd['profiles']:
        N=K.from_json(row['target'])
        if N in targets:continue
        orbits=[];bad=[]
        for z,x in absr[int(N.norm())]:
            assert (z*z-ALPHA*x*x).norm()==N.norm()
            ratio=(z*z-ALPHA*x*x)/N
            if not ratio.is_integral or not ratio.inverse().is_integral:
                bad.append({'representative':enc((z,x)),'ratio':ratio.to_json()});continue
            assert all(ratio.sign(e)>0 for e in range(4))
            unit=square_root(ratio);assert unit is not None and unit.is_integral and unit.inverse().is_integral
            z/=unit;x/=unit
            if z<0:z=-z;x=-x
            orbits.append({'representative':enc((z,x)),'unit_correction':unit.to_json()})
        targets[N]={'target':N.to_json(),'orbits':orbits,'nonunit_rejections':bad}
    starts=[];proposals=[]
    for i,row in enumerate(pd['profiles']):
        q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        for j,orb in enumerate(targets[N]['orbits']):
            u=tuple(map(K.from_json,orb['representative']));power=0
            while u[1]<=0:u=emul(u,eps);power+=1;assert abs(power)<30
            while emul(u,back)[1]>0:u=emul(u,back);power-=1;assert abs(power)<30
            starts.append((i,j,power,u))
            # This is only a proposal search for a finite upper endpoint.
            # Completeness is the subsequent exact all-orbit interval census.
            for step in range(100):
                z,x=u;r=(x/3,*(v+z/9 for v in r0[1:]))
                if root_reason(r,q)is None:
                    proposals.append((ALPHA*x/z,i,j,power,r));break
                u=emul(u,eps);power+=1
    assert proposals,'No proposed cap in the explicitly bounded unit scout.'
    proposals.sort(key=lambda x:x[0]);upper=proposals[0][0]
    visited=[];valid=[]
    for i,j,power,u in starts:
        row=pd['profiles'][i];q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        previous=emul(u,back);first=power;pts=[]
        while ALPHA*u[1]<=upper*u[0]:
            z,x=u;r=(x/3,*(v+z/9 for v in r0[1:]));why=root_reason(r,q)
            rec={'profile':i,'orbit':j,'power':power,'root':enc(r),'q':q.to_json(),'parameter':(ALPHA*x/z).to_json(),'reason':why}
            visited.append(rec);pts.append(power)
            if why is None:valid.append(rec)
            u=emul(u,eps);power+=1;assert power-first<=100
        targets[N]['orbits'][j].setdefault('intervals',[]).append({'profile':i,'first_power':first,'last_power':power-1,
              'previous':enc(previous),'next':enc(u),'visited_powers':pts})
    assert valid and all(K.from_json(r['parameter'])==upper for r in valid)
    edge=INITIAL[:8];simple=[]
    for row in valid:
        root=tuple(map(K.from_json,row['root']));frame=edge+[root];ci=inverse(tuple(tuple(pair(r,G,s)for s in frame)for r in frame))
        coeffs=[]
        for other in valid:
            s=tuple(map(K.from_json,other['root']));cs=mv(ci,tuple(pair(s,G,r)for r in frame))
            if any(c<0 for c in cs):break
            assert tuple(sum((c*r[t]for c,r in zip(cs,frame)),Z)for t in range(10))==s
            coeffs.append(enc(cs))
        else:simple.append((row,coeffs))
    assert len(simple)==1
    chosen,coeffs=simple[0];root=tuple(map(K.from_json,chosen['root']));point=(O,*(upper/3 for _ in range(9)))
    assert pair(point,G,point)<0 and all(pair(point,G,r)==0 for r in edge+[root])
    source=ROOT/'results'/f'{PREFIX}_step9.json';prior=json.loads(source.read_text());walls=[tuple(map(K.from_json,r))for r in prior['all_walls']]
    assert all(compatible(root,w)for w in walls);walls.append(root)
    assert all(pair(point,G,w)<=0 for w in walls)
    result={'status':'PASS_COMPLETE_D9_DIAGONAL_FIRST_ENDPOINT_WITH_CERTIFIED_PARI_INTERFACE',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(PP,RP,GPP,OUT,source)},
        'profiles':len(pd['profiles']),'absolute_norm_counts':{str(n):c for n,c in counts.items()},
        'relative_unit':enc(eps),'targets':list(targets.values()),'upper_parameter':upper.to_json(),
        'visited_points':visited,'primitive_reflective_survivors':valid,'chosen_root':enc(root),
        'chosen_norm':pair(root,G,root).to_json(),'simultaneous_root_cone_coefficients':coeffs,
        'endpoint':enc(point),'endpoint_norm':pair(point,G,point).to_json(),'endpoint_incident':menc(edge+[root]),
        'all_walls':menc(walls),'wall_count':len(walls),
        'scope':'Full first contact on the initial A8 edge, with all safe norms and all integral spherical profiles. The 100-power scout only supplies an upper endpoint; every norm orbit and every point before that endpoint is subsequently exhausted exactly. Completeness depends explicitly on the certified PARI norm/unit interface.'}
    if verify:assert result==json.loads(OP.read_text())
    else:OP.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'profiles',len(pd['profiles']),'points',len(visited),'first mirrors',len(valid),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();run(args.verify)
