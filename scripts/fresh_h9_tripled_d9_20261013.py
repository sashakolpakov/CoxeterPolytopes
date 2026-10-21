#!/usr/bin/env python3
"""Complete initial edges for <−3(T−1)>+D9; safe q|24 and new targets."""
from fresh_h9_changed_spatial_20261013_steps import *

PREFIX='fresh_h9_tripled_d9_20261013'
BETA=3*ALPHA
G=tuple(tuple(-BETA if i==j==0 else O if i==j else Z for j in range(10))for i in range(10))
POOL=tuple(v*q for v in(1,3)for q in DYADIC)
STRUCTURE=ROOT/'results'/f'{PREFIX}_structure.json'


def reason(r,q):
    assert pair(r,G,r)==q
    if not lattice(r):return 'NONINTEGRAL'
    if any(lattice(tuple(x/p for x in r))for p in(T,K(3))):return 'IMPRIMITIVE'
    if not all((2*pair(r,G,b)/q).is_integral for b in DBASIS):return 'NONREFLECTIVE'
    return None


def gate(r,s):
    p=pair(r,G,s)
    return p<=0 and(4*p*p>=pair(r,G,r)*pair(s,G,s)or 4*p*p/(pair(r,G,r)*pair(s,G,s))in ALLOWED)


def structure(verify=False):
    base=json.loads(FIELD.read_text())
    assert base['signature_group_size']==16 and base['prime3_irreducibility_table']
    assert BETA>0 and all(BETA.sign(e)<0 for e in range(1,4))
    gram=tuple(tuple(pair(r,G,s)for s in DBASIS)for r in DBASIS);gi=inverse(gram)
    assert all((12*x).is_integral for row in gi for x in row)
    assert any(not(6*x).is_integral for row in gi for x in row)
    assert any(not(4*x).is_integral for row in gi for x in row)
    for b in(1,3):
        for j,q in enumerate(DYADIC):
            v=b*q
            assert all(v.sign(e)>0 for e in range(4))and(24/v).is_integral
            assert v.norm()==b**4*2**j
            unit=v/(b*PI**j);assert unit.is_integral and unit.inverse().is_integral
    out={'status':'PASS_TRIPLED_D9_STRUCTURE_AND_SAFE_NORMS',
        'sources':{str(FIELD.relative_to(ROOT)):hashlib.sha256(FIELD.read_bytes()).hexdigest()},
        'coefficient':BETA.to_json(),'ambient_gram':menc(G),'lattice_basis':menc(DBASIS),'lattice_gram':menc(gram),
        'determinant':'-12*(T-1)','dual_exponent':12,'safe_norm_divisor':24,'norm_pool':enc(POOL),
        'prime_classes':'Unique dyadic prime (T) with v_T(24)=12, and inert (3) with v_3(24)=1.',
        'primitivity':'Checked in Oe0+D9 at both T and3. Any prime dividing every coordinate of a vector divides its norm squared; the safe pool has no other prime support.',
        'norm_pool_proof':'Reflection integrality places 2r/q in the dual. Its exponent12 gives 24r/q in L; PID primitivity forces24/q integral. Totally positive units are squares, giving26 canonical classes.',
        'zero_chamber':'B9','zero_mirrors':81,
        'zero_completeness':'The totally definite spatial lattice is unchanged. Its signed permutations generate B9 and the finite rank-nine reflection classification excludes a larger irreducible finite reflection group containing B9.',
        'comparison':'The determinant ratio to <-alpha>+D9 is3, nonsquare by odd valuation at the inert prime3. In even rank10 this distinguishes rational quadratic spaces even up to scalar similarity. Repository search found no earlier first-edge search of this tripled D9 host. No nonreflectivity transfer or worldwide novelty claim is made.'}
    if verify:assert out==json.loads(STRUCTURE.read_text())
    else:STRUCTURE.write_text(json.dumps(out,indent=2)+'\n')
    return out


def profiles(k):
    edge=INITIAL[:k-1]+INITIAL[k:]
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
                assert r0[0]==0 and sum(r0[1:k+1],Z)==0
                if not all(x.is_integral for x in r0[k+1:]):return
                N=k*(q-pair(r0,G,r0))
                if any(N.sign(e)<=0 for e in range(4)):return
                assert N.is_integral and N.norm().denominator==1
                rows.append({'q':q.to_json(),'magnitudes':enc(pre),'fixed':enc(r0),
                             'target':N.to_json(),'absolute_norm':int(N.norm())});return
            for p in alph[n]:walk((*pre,p))
        walk(())
    return {'status':'COMPLETE_TRIPLED_D9_INITIAL_FINITE_PROFILES','active_coordinates':k,
        'edge':menc(edge),'gram':menc(C),'alphabets':alphabets,'profiles':rows,
        'absolute_norms':sorted({r['absolute_norm']for r in rows})}


def run(k,verify=False):
    assert k in(1,3);structure(verify)
    tag=f'{PREFIX}_edge{k}';pp=ROOT/'results'/f'{tag}_profiles.json';gp=ROOT/'scripts'/f'{tag}.gp'
    rp=ROOT/'results'/f'{tag}_raw.txt';op=ROOT/'results'/f'{tag}.json'
    field=3*ALPHA if k==1 else ALPHA;timefactor=1 if k==1 else 3
    pd=profiles(k)
    if verify:assert pd==json.loads(pp.read_text())
    else:
        pp.write_text(json.dumps(pd,indent=2)+'\n')
        text=GP.read_text().replace('d=t-1;',f'd={3 if k==1 else 1}*(t-1);').replace('ns=vector(13,j,2^(j-1));',f'ns={pd["absolute_norms"]};')
        gp.write_text(text)
        with rp.open('w')as f:
            job=subprocess.run(['gp','-q',str(gp)],stdout=f,stderr=subprocess.PIPE,text=True,timeout=90)
            assert job.returncode==0,job.stderr
    raw=rp.read_text();absolute={};counts={};bases=[];units=[];norms=[]
    for s in('CERTIFY 1\n','BASE_CERTIFY 1\n','CLASS 1\n','SIGNATURE [2, 3]\n'):assert s in raw
    polynomial='x^8 + 12*x^6 + 18*x^4 - 108*x^2 - 81'if k==1 else'x^8 + 4*x^6 + 2*x^4 - 4*x^2 - 1'
    assert 'ABSOLUTE '+polynomial+'\n'in raw
    for line in raw.splitlines():
        if line.startswith('BASE_UNIT '):bases.append(K.from_json(json.loads(line.split(' ',2)[2])))
        elif line.startswith('UNIT '):units.append(tuple(map(K.from_json,json.loads(line.split(' ',2)[2]))))
        elif line.startswith('UNIT_NORM '):
            n,e=line.split(' ',2)[2].split('] ',1);norms.append((K.from_json(json.loads(n+']')),json.loads(e)))
        elif line.startswith('NORM '):
            _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
        elif line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(tuple(map(K.from_json,json.loads(s))))
    expect=[[0,2,0,0],[-2,0,0,0],[0,0,2,0],[0,0,0,0]]if k==1 else[[0,-1,0,0],[2,-1,0,0],[0,0,2,0],[0,0,0,0]]
    assert [e for n,e in norms]==expect and len(bases)==3 and len(units)==4
    for(z,x),(n,e)in zip(units,norms):
        assert z*z-field*x*x==n==K((-1)**e[3])*bases[0]**e[0]*bases[1]**e[1]*bases[2]**e[2]
    eps=units[3]
    if eps[0]<0:eps=tuple(-x for x in eps)
    if eps[1]<0:eps=eps[0],-eps[1]
    assert eps[0]>1 and eps[1]>0 and eps[0]**2-field*eps[1]**2==1
    assert set(absolute)==set(pd['absolute_norms'])and all(len(v)==counts[n]for n,v in absolute.items())
    print('edge',k,'profiles',len(pd['profiles']),'absolute targets',len(counts),'representatives',sum(counts.values()),flush=True)
    def mul(u,v):return u[0]*v[0]+field*u[1]*v[1],u[0]*v[1]+u[1]*v[0]
    def root(r0,u):return(u[1]/timefactor,*(r0[i]+u[0]/k for i in range(1,k+1)),*r0[k+1:])
    def parameter(u):return BETA*u[1]/(timefactor*u[0])
    back=eps[0],-eps[1];targets={};starts=[];proposals=[]
    for row in pd['profiles']:
        N=K.from_json(row['target'])
        if N in targets:continue
        orbits=[];bad=[]
        for z,x in absolute[int(N.norm())]:
            assert(z*z-field*x*x).norm()==N.norm()
            ratio=(z*z-field*x*x)/N
            if not ratio.is_integral or not ratio.inverse().is_integral:
                bad.append({'representative':enc((z,x)),'ratio':ratio.to_json()});continue
            assert all(ratio.sign(e)>0 for e in range(4))
            u=square_root(ratio);assert u is not None and u.is_integral and u.inverse().is_integral
            z/=u;x/=u
            if z<0:z=-z;x=-x
            assert z*z-field*x*x==N
            orbits.append({'representative':enc((z,x)),'unit_correction':u.to_json()})
        targets[N]={'target':N.to_json(),'orbits':orbits,'nonunit_rejections':bad}
    for i,row in enumerate(pd['profiles']):
        q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        for j,orb in enumerate(targets[N]['orbits']):
            u=tuple(map(K.from_json,orb['representative']));power=0
            while u[1]<=0:u=mul(u,eps);power+=1;assert abs(power)<50
            while mul(u,back)[1]>0:u=mul(u,back);power-=1;assert abs(power)<50
            starts.append((i,j,power,u))
            for step in range(100):
                r=root(r0,u)
                if reason(r,q)is None:proposals.append((parameter(u),i,j,power,r));break
                u=mul(u,eps);power+=1
    assert proposals,'No endpoint proposed in100 unit powers; no completeness claim.'
    proposals.sort(key=lambda x:x[0]);upper=proposals[0][0]
    print('edge',k,'upper proposal',proposals[0][1:4],flush=True)
    visited=[];valid=[]
    for i,j,power,u in starts:
        row=pd['profiles'][i];q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        prev=mul(u,back);first=power;pts=[]
        while parameter(u)<=upper:
            r=root(r0,u);why=reason(r,q)
            rec={'profile':i,'orbit':j,'power':power,'root':enc(r),'q':q.to_json(),'parameter':parameter(u).to_json(),'reason':why}
            visited.append(rec);pts.append(power)
            if why is None:valid.append(rec)
            u=mul(u,eps);power+=1;assert power-first<=200
        targets[N]['orbits'][j].setdefault('intervals',[]).append({'profile':i,'first_power':first,'last_power':power-1,
            'previous':enc(prev),'next':enc(u),'visited_powers':pts})
    assert valid and all(K.from_json(row['parameter'])==upper for row in valid)
    edge=[tuple(map(K.from_json,r))for r in pd['edge']];simple=[]
    for row in valid:
        r=tuple(map(K.from_json,row['root']));frame=edge+[r];ci=inverse(tuple(tuple(pair(s,G,t)for t in frame)for s in frame));allcs=[]
        for other in valid:
            s=tuple(map(K.from_json,other['root']));cs=mv(ci,tuple(pair(s,G,t)for t in frame))
            if any(c<0 for c in cs):break
            assert tuple(sum((c*t[i]for c,t in zip(cs,frame)),Z)for i in range(10))==s
            allcs.append(enc(cs))
        else:simple.append((row,allcs))
    assert len(simple)==1
    chosen,cone=simple[0];r=tuple(map(K.from_json,chosen['root']));walls=INITIAL+[r]
    point=(O,*(upper for _ in range(k)),*(Z for _ in range(9-k)))
    assert pair(point,G,point)<0 and all(pair(point,G,s)==0 for s in edge+[r])
    assert pair(point,G,INITIAL[k-1])<0 and pair((O,*(Z for _ in range(9))),G,r)<0
    assert all(gate(r,s)for s in INITIAL)
    delta=O
    while True:
        inside=(O,*(delta*K(9-i)for i in range(9)))
        if pair(inside,G,inside)<0 and all(pair(inside,G,s)<0 for s in walls):break
        delta/=2
    result={'status':'PASS_COMPLETE_TRIPLED_D9_FIRST_EDGE','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(STRUCTURE,pp,gp,rp)},
        'active_coordinates':k,'coefficient':BETA.to_json(),'field_parameter':field.to_json(),'time_divisor':timefactor,
        'base_units':enc(bases),'fundamental_units':[enc(u)for u in units],'relative_norm_map':expect,'relative_unit':enc(eps),
        'profiles':len(pd['profiles']),'absolute_norm_counts':{str(n):c for n,c in counts.items()},
        'targets':list(targets.values()),'proposal_power_limit':100,'upper_parameter':upper.to_json(),
        'visited_points':visited,'primitive_reflective_survivors':valid,'chosen_root':enc(r),'chosen_norm':pair(r,G,r).to_json(),
        'simultaneous_root_cone_coefficients':cone,'endpoint':enc(point),'endpoint_norm':pair(point,G,point).to_json(),
        'endpoint_incident':menc(edge+[r]),'all_walls':menc(walls),'strict_interior':enc(inside),
        'scope':'Full first contact for this initial edge of the tripled D9 lattice. All26 safe primitive norm classes, all integral spherical profiles and all relative norm-one unit orbits are included. The100-power scout supplies an upper endpoint only; the finite interval is exhausted afterward. No old-wall compatibility pruning precedes firstness. Certified PARI norm/unit completeness is explicit.'}
    if verify:assert result==json.loads(op.read_text())
    else:op.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'edge',k,'profiles',len(pd['profiles']),'points',len(visited),'first mirrors',len(valid),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--edge',type=int,choices=(1,3),required=True);p.add_argument('--verify',action='store_true')
    args=p.parse_args();run(args.edge,args.verify)
