#!/usr/bin/env python3
"""Independent Fraction census of tripled-D9 axis-chain stages2 and3."""
from fresh_h9_tripled_d9_20261013_independent import *

OUT=ROOT/'results'/f'{PREFIX}_chain_independent.json'


def run(verify=False):
    firstp=ROOT/'results'/f'{PREFIX}_edge1.json';first=json.loads(firstp.read_text())
    initial_audit=ROOT/'results'/f'{PREFIX}_independent.json';ia=json.loads(initial_audit.read_text())
    assert ia['status']=='PASS_INDEPENDENT_TRIPLED_D9_COMPLETE_INITIAL_EDGES'
    for p,h in ia['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    sources={str(initial_audit.relative_to(ROOT)):hashlib.sha256(initial_audit.read_bytes()).hexdigest()}
    eps=a.vec(first['relative_unit']);back=eps[0],a.neg(eps[1]);rows=[]
    def emul(u,v):return a.add(a.mul(u[0],v[0]),a.mul(BETA,a.mul(u[1],v[1]))),a.add(a.mul(u[0],v[1]),a.mul(u[1],v[0]))
    base_raw=(ROOT/'results'/f'{PREFIX}_edge1_raw.txt').read_text()
    unitlines=lambda text:[s for s in text.splitlines()if s.startswith(('UNIT ','UNIT_NORM ','BASE_UNIT '))]
    for stage in(2,3):
        tag=f'{PREFIX}_step{stage}';path=ROOT/'results'/f'{tag}.json';d=json.loads(path.read_text())
        sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
        for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
        pd=json.loads((ROOT/'results'/f'{tag}_profiles.json').read_text());free=pd['free_coordinate']
        nprofiles=profiles(pd,1,free);assert nprofiles==d['profiles']
        A=old.matrix(d['frame_to_axis']);AI=old.matrix(d['frame_from_axis'])
        assert a.minverse(A)==AI and all(a.integral(x)for M in(A,AI)for row in M for x in row)
        assert all(P(a.mv(A,r),a.mv(A,s))==P(r,s)for r in old.BASIS for s in old.BASIS)
        raw=(ROOT/'results'/f'{tag}_raw.txt').read_text();assert unitlines(raw)==unitlines(base_raw)
        for s in('CERTIFY 1\n','BASE_CERTIFY 1\n','CLASS 1\n','SIGNATURE [2, 3]\n'):assert s in raw
        absolute={};counts={}
        for line in raw.splitlines():
            if line.startswith('NORM '):
                _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
            elif line.startswith('SOLUTION '):
                _,n,s=line.split(' ',2);absolute[int(n)].append(a.vec(json.loads(s)))
        assert set(absolute)==set(pd['absolute_norms'])and all(len(v)==counts[n]for n,v in absolute.items())
        upper=a.dec(d['upper_parameter']);visited=[];valid=[]
        for i,row in enumerate(pd['profiles']):
            q,N=a.dec(row['q']),a.dec(row['target']);r0=a.vec(row['fixed']);orbits=[]
            for z,x in absolute[int(old.norm(N))]:
                nu=a.sub(a.sq(z),a.mul(BETA,a.sq(x)));assert old.norm(nu)==old.norm(N)
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
                firstn=powern
                while a.sign(a.sub(a.mul(upper,u[0]),a.mul(BETA,u[1])))>=0:
                    z,x=u;rr=list(r0);rr[0]=x;rr[free]=z;r=a.mv(AI,rr);why=reason(r,q)
                    visited.append((i,j,powern,r,why))
                    if why is None:valid.append((a.divide(a.mul(BETA,x),z),r))
                    u=emul(u,eps);powern+=1;assert powern-firstn<=100
        assert visited==[(r['profile'],r['orbit'],r['power'],a.vec(r['root']),r['reason'])for r in d['visited_points']]
        assert valid
        closest=min((p for p,r in valid),key=__import__('functools').cmp_to_key(lambda x,y:a.sign(a.sub(x,y))))
        nearest=[r for p,r in valid if p==closest];chosen=a.vec(d['chosen_root']);assert chosen in nearest
        edge=tuple(a.mv(AI,r)for r in old.matrix(pd['edge']));frame=edge+(chosen,)
        ci=a.minverse(tuple(tuple(P(r,s)for s in frame)for r in frame))
        for i,r in enumerate(nearest):
            cs=a.mv(ci,tuple(P(r,s)for s in frame));assert all(a.sign(c)>=0 for c in cs)
            assert tuple(old.summ(a.mul(c,s[t])for c,s in zip(cs,frame))for t in range(10))==r
            assert cs==a.vec(d['simultaneous_root_cone_coefficients'][i])
        expected=[a.Z]*10;expected[0]=eps[1];expected[stage-1]=a.neg(a.O);expected[stage]=eps[0]
        assert chosen==a.mv(AI,expected)and closest==upper and d['expected_staircase_continues']
        walls=old.matrix(d['all_walls']);point=a.vec(d['endpoint'])
        pp=[a.Z]*10;pp[0]=a.O;pp[free]=closest;pp=a.mv(AI,pp);pp=tuple(a.divide(x,pp[0])for x in pp)
        assert point==pp and a.sign(P(point,point))<0
        assert all(P(point,r)==a.Z for r in frame)and all(a.sign(P(point,r))<=0 for r in walls)
        for i,r in enumerate(walls):
            assert reason(r,P(r,r))is None
            for s in walls[:i]:gate(r,s)
        rows.append({'stage':stage,'profiles':nprofiles,'interval_points':len(visited),'nearest_mirrors':len(nearest)})
    otherp=ROOT/'results'/f'{PREFIX}_edge3.json';other=json.loads(otherp.read_text())
    sources[str(otherp.relative_to(ROOT))]=hashlib.sha256(otherp.read_bytes()).hexdigest()
    joint=walls+(a.vec(other['chosen_root']),)
    for i,r in enumerate(joint):
        for s in joint[:i]:gate(r,s)
    inside=a.vec(ia['strict_interior']);escape=a.vec(ia['strict_escape'])
    assert a.sign(P(inside,inside))<0 and all(a.sign(P(inside,r))<0 for r in joint)
    escape_survives=all(a.sign(P(escape,r))<0 for r in joint)
    out={'status':'PASS_INDEPENDENT_TRIPLED_D9_CHAIN_EDGES','sources':sources,'stages':rows,
        'profiles':sum(x['profiles']for x in rows),'interval_points':sum(x['interval_points']for x in rows),
        'combined_walls':len(joint),'combined_pairs':len(joint)*(len(joint)-1)//2,
        'strict_interior':list(map(a.enc,inside)),'initial_escape_survives':escape_survives,
        'scope':'Separate Fraction arithmetic, complete reversed q|24 profiles and all relative-unit intervals for two new chain edges, plus actual D9 root tests and all combined pairings. Initial edge unit completeness is inherited by source hashes; PARI remains explicit. No rational-space nonreflectivity transfer is assumed.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['profiles'],out['interval_points'],'escape',escape_survives,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();run(args.verify)
