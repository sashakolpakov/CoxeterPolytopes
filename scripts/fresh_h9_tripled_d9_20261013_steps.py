#!/usr/bin/env python3
"""Two full axis-chain steps in the tripled D9 host, with all q|24."""
from fresh_h9_tripled_d9_20261013 import *


def B(j,eps):
    z,x=eps;M=[list(r)for r in eye()]
    M[0][0]=M[j][j]=z;M[0][j]=-x;M[j][0]=-BETA*x
    return tuple(map(tuple,M))


def finite_profiles(k):
    roots=corner(k);edge=roots[:k]+roots[k+1:];free=k+1
    assert all(r[free]==0 for r in edge)
    C=tuple(tuple(pair(r,G,s)for s in edge)for r in edge);ci=inverse(C)
    rows=[];alphabets=[]
    for q in POOL:
        alph=[]
        for r in edge:
            aa={Z}
            for lam in ALLOWED:
                p=square_root(q*pair(r,G,r)*lam/4)
                if p is not None and p.is_integral and(2*p/q).is_integral:aa.add(p)
            alph.append(sorted(aa))
        alphabets.append({'q':q.to_json(),'magnitudes':[enc(a)for a in alph]})
        def walk(pre):
            n=len(pre)
            if n and pair(pre,tuple(row[:n]for row in ci[:n]),pre)>=q:return
            if n==8:
                cs=mv(ci,tuple(-v for v in pre))
                r0=tuple(sum((c*r[j]for c,r in zip(cs,edge)),Z)for j in range(10))
                assert r0[0]==r0[free]==0
                if not all(x.is_integral for x in r0):return
                N=q-pair(r0,G,r0)
                if any(N.sign(e)<=0 for e in range(4)):return
                assert N.is_integral
                rows.append({'q':q.to_json(),'magnitudes':enc(pre),'fixed':enc(r0),
                             'target':N.to_json(),'absolute_norm':int(N.norm())});return
            for p in alph[n]:walk((*pre,p))
        walk(())
    return {'status':'COMPLETE_TRIPLED_D9_CHAIN_PROFILES','stage':k+1,'free_coordinate':free,
        'edge':menc(edge),'gram':menc(C),'profiles':rows,'alphabets':alphabets,
        'absolute_norms':sorted({r['absolute_norm']for r in rows})}


def run(k,verify=False):
    assert k in(1,2)
    tag=f'{PREFIX}_step{k+1}';pp=ROOT/'results'/f'{tag}_profiles.json';gp=ROOT/'scripts'/f'{tag}.gp'
    rp=ROOT/'results'/f'{tag}_raw.txt';op=ROOT/'results'/f'{tag}.json'
    first=ROOT/'results'/f'{PREFIX}_edge1.json';fd=json.loads(first.read_text());eps=tuple(map(K.from_json,fd['relative_unit']))
    pd=finite_profiles(k)
    if verify:assert pd==json.loads(pp.read_text())
    else:
        pp.write_text(json.dumps(pd,indent=2)+'\n')
        gp.write_text(GP.read_text().replace('d=t-1;','d=3*(t-1);').replace('ns=vector(13,j,2^(j-1));',f'ns={pd["absolute_norms"]};'))
        with rp.open('w')as f:
            job=subprocess.run(['gp','-q',str(gp)],stdout=f,stderr=subprocess.PIPE,text=True,timeout=90)
            assert job.returncode==0,job.stderr
    raw=rp.read_text();base_raw=(ROOT/'results'/f'{PREFIX}_edge1_raw.txt').read_text()
    unitlines=lambda text:[s for s in text.splitlines()if s.startswith(('UNIT ','UNIT_NORM ','BASE_UNIT '))]
    assert unitlines(raw)==unitlines(base_raw)
    for s in('CERTIFY 1\n','BASE_CERTIFY 1\n','CLASS 1\n','SIGNATURE [2, 3]\n'):assert s in raw
    absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('NORM '):
            _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
        elif line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(tuple(map(K.from_json,json.loads(s))))
    assert set(absolute)==set(pd['absolute_norms'])and all(len(v)==counts[n]for n,v in absolute.items())
    A=eye()
    for j in range(1,k+1):A=mm(B(j,eps),A)
    AI=inverse(A);assert mm(mm(tuple(zip(*A)),G),A)==G
    assert all(z.is_integral for row in A+AI for z in row)
    free=k+1;upper=BETA*eps[1]/eps[0];back=eps[0],-eps[1]
    def mul(u,v):return u[0]*v[0]+BETA*u[1]*v[1],u[0]*v[1]+u[1]*v[0]
    targets={};visited=[];valid=[]
    for i,row in enumerate(pd['profiles']):
        q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        if N not in targets:
            orbits=[];bad=[]
            for z,x in absolute[int(N.norm())]:
                assert(z*z-BETA*x*x).norm()==N.norm()
                ratio=(z*z-BETA*x*x)/N
                if not ratio.is_integral or not ratio.inverse().is_integral:
                    bad.append({'representative':enc((z,x)),'ratio':ratio.to_json()});continue
                assert all(ratio.sign(e)>0 for e in range(4))
                u=square_root(ratio);assert u is not None and u.is_integral and u.inverse().is_integral
                z/=u;x/=u
                if z<0:z=-z;x=-x
                orbits.append({'representative':enc((z,x)),'unit_correction':u.to_json()})
            targets[N]={'target':N.to_json(),'orbits':orbits,'nonunit_rejections':bad}
        for j,orb in enumerate(targets[N]['orbits']):
            u=tuple(map(K.from_json,orb['representative']));power=0
            while u[1]<=0:u=mul(u,eps);power+=1;assert abs(power)<50
            while mul(u,back)[1]>0:u=mul(u,back);power-=1;assert abs(power)<50
            prev=mul(u,back);start=power;pts=[]
            while BETA*u[1]<=upper*u[0]:
                z,x=u;rr=list(r0);rr[0]=x;rr[free]=z;r=mv(AI,rr);why=reason(r,q)
                row={'profile':i,'orbit':j,'power':power,'root':enc(r),'transformed_root':enc(rr),
                     'q':q.to_json(),'parameter':(BETA*x/z).to_json(),'reason':why}
                visited.append(row);pts.append(power)
                if why is None:valid.append(row)
                u=mul(u,eps);power+=1;assert power-start<=100
            orb.setdefault('intervals',[]).append({'profile':i,'first_power':start,'last_power':power-1,
                'previous':enc(prev),'next':enc(u),'visited_powers':pts})
    valid.sort(key=lambda r:K.from_json(r['parameter']));assert valid
    nearest=[r for r in valid if r['parameter']==valid[0]['parameter']]
    edge=[mv(AI,tuple(map(K.from_json,r)))for r in pd['edge']];simple=[]
    for row in nearest:
        r=tuple(map(K.from_json,row['root']));frame=edge+[r]
        ci=inverse(tuple(tuple(pair(s,G,t)for t in frame)for s in frame));allcs=[]
        for other in nearest:
            s=tuple(map(K.from_json,other['root']));cs=mv(ci,tuple(pair(s,G,t)for t in frame))
            if any(c<0 for c in cs):break
            assert tuple(sum((c*t[i]for c,t in zip(cs,frame)),Z)for i in range(10))==s
            allcs.append(enc(cs))
        else:simple.append((row,allcs))
    assert len(simple)==1
    chosen,cone=simple[0];r=tuple(map(K.from_json,chosen['root']))
    prior=first if k==1 else ROOT/'results'/f'{PREFIX}_step{k}.json';old=json.loads(prior.read_text())
    if k==2:assert old['expected_staircase_continues']
    walls=[tuple(map(K.from_json,w))for w in old['all_walls']]
    assert all(gate(r,w)for w in walls);walls.append(r)
    p=[Z]*10;p[0]=O;p[free]=K.from_json(chosen['parameter']);point=mv(AI,p);point=tuple(x/point[0]for x in point)
    assert pair(point,G,point)<0 and all(pair(point,G,w)<=0 for w in walls)
    assert all(pair(point,G,w)==0 for w in edge+[r])
    expected=[Z]*10;expected[0]=eps[1];expected[k]=-O;expected[free]=eps[0];expected=mv(AI,expected)
    continues=r==expected and chosen['parameter']==upper.to_json()
    out={'status':'PASS_COMPLETE_TRIPLED_D9_CHAIN_EDGE','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(pp,gp,rp,prior,first,STRUCTURE)},
        'stage':k+1,'profiles':len(pd['profiles']),'frame_to_axis':menc(A),'frame_from_axis':menc(AI),
        'relative_unit':enc(eps),'absolute_norm_counts':{str(n):c for n,c in counts.items()},
        'upper_parameter':upper.to_json(),'targets':list(targets.values()),'visited_points':visited,
        'primitive_reflective_survivors':valid,'nearest_candidates':nearest,'chosen_root':enc(r),'chosen_norm':pair(r,G,r).to_json(),
        'simultaneous_root_cone_coefficients':cone,'endpoint':enc(point),'endpoint_norm':pair(point,G,point).to_json(),
        'endpoint_incident':menc(edge+[r]),'all_walls':menc(walls),'expected_staircase_continues':continues,
        'scope':'Complete q|24 spherical profiles and all relative-unit orbits on this compact edge segment. Boosts preserve the surrounding integral coordinate module only; all pulled-back roots are tested for actual D9 primitivity/reflection integrality at2 and3. No prefix compatibility pruning precedes first contact. Certified PARI interface remains an explicit dependency.'}
    if verify:assert out==json.loads(op.read_text())
    else:op.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'stage',k+1,'profiles',len(pd['profiles']),'points',len(visited),'nearest',len(nearest),'continues',continues,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',type=int,choices=(2,3),required=True);p.add_argument('--verify',action='store_true')
    args=p.parse_args();run(args.stage-1,args.verify)
