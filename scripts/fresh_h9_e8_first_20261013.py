#!/usr/bin/env python3
"""Complete E8-fixed first edge for <-3(T-1)>+E8+<1>."""
import argparse,hashlib,json,subprocess
from pathlib import Path
from fresh_h9_glued_host_20261013_arithmetic import *
from fresh_h9_changed_spatial_20261013 import square_root,ALLOWED

ROOT=Path(__file__).resolve().parents[1];PREFIX='fresh_h9_e8_first_20261013'
TRIAGE=ROOT/'results/fresh_h9_tripled_d9_20261013_e8_triage.json'
BETA=3*ALPHA
G=tuple(tuple(-BETA if i==j==0 else O if i==j else Z for j in range(10))for i in range(10))
POOL=tuple(b*q for b in(1,3)for q in NORMS)
OUT=ROOT/'results'/f'{PREFIX}.json'

def run(verify=False):
    st=json.loads(TRIAGE.read_text());basis=tuple(tuple(map(K.from_json,r))for r in st['full_lattice_basis'])
    e8=list(basis[1:9]);B=tuple(zip(*basis));BI=inverse(B)
    initial=e8+[tuple(-x for x in basis[9])]
    gram=tuple(tuple(pair(r,G,s)for s in basis)for r in basis);gi=inverse(gram)
    assert all((3*x).is_integral for row in gi for x in row)
    assert list(map(K.from_json,st['safe_norm_pool']))==list(POOL)
    def integral(r):return all(x.is_integral for x in mv(BI,r))
    def reason(r,q):
        assert pair(r,G,r)==q
        if not integral(r):return'NONINTEGRAL'
        if any(integral(tuple(x/p for x in r))for p in(T,K(3))):return'IMPRIMITIVE'
        if not all((2*pair(r,G,b)/q).is_integral for b in basis):return'NONREFLECTIVE'
        return None
    # Classification forces E8+A1 at a first finite contact. Thus the
    # complete finite-profile list has zero pairings on all eight E8 walls.
    profiles=[{'q':q.to_json(),'target':q.to_json(),'magnitudes':enc((Z,)*8)}for q in POOL]
    gp=ROOT/'scripts'/f'{PREFIX}.gp';rawp=ROOT/'results'/f'{PREFIX}_raw.txt'
    ns=sorted({int(q.norm())for q in POOL})
    if not verify:
        text=(ROOT/'scripts/fresh_h9_tripled_d9_20261013_edge1.gp').read_text()
        line=next(s for s in text.splitlines()if s.startswith('ns='))
        gp.write_text(text.replace(line,'ns='+str(ns)+';'))
        with rawp.open('w')as f:
            p=subprocess.run(['gp','-q',str(gp)],stdout=f,stderr=subprocess.PIPE,text=True,timeout=90)
            assert p.returncode==0,p.stderr
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
    expect=[[0,2,0,0],[-2,0,0,0],[0,0,2,0],[0,0,0,0]]
    assert[e for n,e in norms]==expect and len(bases)==3 and len(units)==4
    for(z,x),(n,e)in zip(units,norms):
        assert z*z-BETA*x*x==n==K((-1)**e[3])*bases[0]**e[0]*bases[1]**e[1]*bases[2]**e[2]
    eps=units[3]
    if eps[0]<0:eps=tuple(-x for x in eps)
    if eps[1]<0:eps=eps[0],-eps[1]
    assert eps[0]>1 and eps[1]>0 and eps[0]**2-BETA*eps[1]**2==1
    def mul(u,v):return u[0]*v[0]+BETA*u[1]*v[1],u[0]*v[1]+u[1]*v[0]
    def root(u):return(u[1],*(Z for _ in range(8)),u[0])
    def parameter(u):return BETA*u[1]/u[0]
    back=eps[0],-eps[1];targets=[];starts=[];proposals=[]
    for i,q in enumerate(POOL):
        orbits=[];bad=[]
        for z,x in absolute[int(q.norm())]:
            nu=z*z-BETA*x*x;assert nu.norm()==q.norm();ratio=nu/q
            if not ratio.is_integral or not ratio.inverse().is_integral:
                bad.append({'representative':enc((z,x)),'ratio':ratio.to_json()});continue
            assert all(ratio.sign(e)>0 for e in range(4))
            u=square_root(ratio);assert u is not None and u.is_integral and u.inverse().is_integral
            z/=u;x/=u
            if z<0:z=-z;x=-x
            assert z*z-BETA*x*x==q
            orbits.append({'representative':enc((z,x)),'unit_correction':u.to_json()})
        targets.append({'q':q.to_json(),'orbits':orbits,'nonunit_rejections':bad})
        for j,orb in enumerate(orbits):
            u=tuple(map(K.from_json,orb['representative']));power=0
            while u[1]<=0:u=mul(u,eps);power+=1;assert abs(power)<50
            while mul(u,back)[1]>0:u=mul(u,back);power-=1;assert abs(power)<50
            starts.append((i,j,power,u))
            for step in range(100):
                r=root(u)
                if reason(r,q)is None:proposals.append((parameter(u),i,j,power,r));break
                u=mul(u,eps);power+=1
    assert proposals
    proposals.sort(key=lambda x:x[0]);upper=proposals[0][0]
    visited=[];valid=[]
    for i,j,power,u in starts:
        q=POOL[i];previous=mul(u,back);first=power;pts=[]
        assert previous[1]<=0
        while parameter(u)<=upper:
            r=root(u);why=reason(r,q)
            rec={'profile':i,'orbit':j,'power':power,'root':enc(r),'q':q.to_json(),'parameter':parameter(u).to_json(),'reason':why}
            visited.append(rec);pts.append(power)
            if why is None:valid.append(rec)
            u=mul(u,eps);power+=1;assert power-first<=200
        targets[i]['orbits'][j]['interval']={'first_power':first,'last_power':power-1,'previous':enc(previous),'next':enc(u),'visited_powers':pts}
    assert len(valid)==1 and K.from_json(valid[0]['parameter'])==upper
    r=tuple(map(K.from_json,valid[0]['root']));q=pair(r,G,r)
    point=(O,*(Z for _ in range(8)),upper);origin=(O,*(Z for _ in range(9)))
    walls=initial+[r];incident=e8+[r]
    assert pair(point,G,point)<0 and all(pair(point,G,s)==0 for s in incident)
    assert pair(point,G,initial[-1])<0 and pair(origin,G,r)<0
    for s in initial:
        p=pair(s,G,r);assert p<=0 and(p==0 or 4*p*p/(pair(s,G,s)*q)>=4)
    # Reconstruct the E8 mirror system in its simple-root coordinates.
    eg=tuple(tuple(pair(s,G,t)for t in e8)for s in e8)
    seeds=[tuple(int(i==j)for j in range(8))for i in range(8)];seen=set(seeds);todo=list(seeds)
    while todo:
        v=todo.pop()
        for j in range(8):
            p=sum(v[i]*int(eg[i][j].a)for i in range(8));w=list(v);w[j]-=p;w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=240
    assert len(seen)==240 and all(tuple(-x for x in v)in seen for v in seen)
    # Full finite stabilizers are E8+A1 at both ends by classification;
    # the root closure is 120+1 unoriented mirrors at each.
    spatial=mv(inverse(eg),tuple(-O for _ in range(8)))
    dvec=tuple(sum((a*b[j]for a,b in zip(spatial,e8)),Z)for j in range(10))
    delta=O
    while True:
        center=(O,*(delta*dvec[i]for i in range(1,9)),delta)
        if pair(center,G,center)<0 and all(pair(center,G,s)<0 for s in walls):break
        delta/=2
    # Extend only the E8 spatial ray. Choosing f=(1+C)/2 gives
    # f^2>C but f<C for C>1, hence a spacelike escape future of center.
    S=delta*delta*pair(dvec,G,dvec);C=(BETA-delta*delta)/S;assert C>1
    f=(1+C)/2
    escape=(O,*(f*delta*dvec[i]for i in range(1,9)),delta)
    assert pair(escape,G,escape)>0 and pair(center,G,escape)<0
    assert all(pair(escape,G,s)<0 for s in walls)
    result={'status':'PASS_COMPLETE_E8_FIXED_FIRST_EDGE',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(TRIAGE,gp,rawp,Path(__file__).resolve(),ROOT/'scripts/fresh_h9_glued_host_20261013_arithmetic.py',Path(square_root.__code__.co_filename).resolve())},
        'coefficient':BETA.to_json(),'lattice_basis':menc(basis),'safe_norm_pool':enc(POOL),
        'base_units':enc(bases),'fundamental_units':[enc(u)for u in units],'relative_norm_map':expect,'relative_unit':enc(eps),
        'profiles':profiles,'absolute_norm_counts':{str(n):c for n,c in counts.items()},'targets':targets,
        'proposal_power_limit':100,'upper_parameter':upper.to_json(),'visited_points':visited,
        'primitive_reflective_survivors':valid,'chosen_root':enc(r),'chosen_norm':q.to_json(),
        'endpoint':enc(point),'endpoint_norm':pair(point,G,point).to_json(),'endpoint_incident':menc(incident),
        'all_walls':menc(walls),'strict_interior':enc(center),'strict_escape':enc(escape),'escape_norm':pair(escape,G,escape).to_json(),
        'full_corner_types':['E8+A1','E8+A1'],'finite_mirror_counts':[121,121],
        'scope':'Complete first contact on the E8-fixed initial edge, including all ten safe primitive norm classes and all relative norm-one unit orbits. First finite contact must have type E8+A1, hence zero E8 pairings. A 100-power scout only supplies an upper bound, followed by a complete interval census. Certified PARI norm/unit completeness is an explicit dependency. Full indefinite conjugacy with prior hosts remains unresolved; the ten-wall prefix has infinite volume.'}
    if verify:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'profiles',len(profiles),'points',len(visited),'q',q.to_json(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
