#!/usr/bin/env python3
"""Complete initial A6+B2 edge for a new prime-seven quartic H9 host."""
import argparse,hashlib,json,subprocess,math
from pathlib import Path
from fractions import Fraction
from fresh_h9_glued_host_20261013_arithmetic import *
from fresh_h9_changed_spatial_20261013 import square_root,ALLOWED,DYADIC
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_a6_d3_20261013_second';OUT=ROOT/'results'/f'{TAG}.json'
BETA=ALPHA*PI;D=ALPHA;P7=(T*T+1,5-T*T);POOL=tuple(q*a*b for a in(O,P7[0])for b in(O,P7[1])for q in DYADIC)
G=tuple(tuple(-BETA if i==j==0 else K(2)if i==j and 1<=i<=6 else -O if 1<=i<=6 and 1<=j<=6 and abs(i-j)==1 else O if i==j else Z for j in range(10))for i in range(10))
STD=[tuple(O if i==j else Z for i in range(10))for j in range(10)]
def vadd(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(v):return tuple(-x for x in v)
BASIS=STD[:7]+[vadd(STD[7],neg(STD[8])),vadd(STD[8],neg(STD[9])),vadd(STD[8],STD[9])]
INITIAL=[neg(v)for v in STD[1:7]]+[neg(BASIS[7]),neg(BASIS[8]),tuple(-2*x for x in STD[9])]
BI=inverse(tuple(zip(*BASIS)))
def integral(r):return all(x.is_integral for x in mv(BI,r))
def reason(r,q):
    assert pair(r,G,r)==q
    if not integral(r):return'NONINTEGRAL'
    if any(integral(tuple(x/p for x in r))for p in(T,*P7)):return'IMPRIMITIVE'
    if not all((2*pair(r,G,b)/q).is_integral for b in BASIS):return'NONREFLECTIVE'
    return None

def det(M):
    a=[list(r)for r in M];d=O
    for i in range(len(a)):
        ks=[j for j in range(i,len(a))if a[j][i]]
        if not ks:return Z
        k=ks[0];a[i],a[k]=a[k],a[i]
        if k!=i:d=-d
        p=a[i][i];d*=p
        for j in range(i+1,len(a)):
            c=a[j][i]/p
            for k in range(i+1,len(a)):a[j][k]-=c*a[i][k]
    return d

def run(verify=False):
    fieldpath=ROOT/'results/fresh_h9_glued_host_20261013_structure.json';field=json.loads(fieldpath.read_text())
    assert field['signature_group_size']==16 and field['class_number_argument']
    assert BETA>0 and all(BETA.sign(e)<0 for e in(1,2,3))
    assert P7[0]*P7[1]==7 and all(p.norm()==49 and all(p.sign(e)>0 for e in range(4))for p in P7)
    gram=tuple(tuple(pair(r,G,s)for s in BASIS)for r in BASIS);gi=inverse(gram);assert det(gram)==-28*BETA
    assert all((28*x).is_integral for row in gi for x in row)
    assert any(not(14*x).is_integral for row in gi for x in row)and any(not(4*x).is_integral for row in gi for x in row)
    assert len(set(POOL))==52 and all(all(q.sign(e)>0 for e in range(4))for q in POOL)
    # Exact trace Gram of the power basis, and positive-square decomposition.
    tracegram=((4,0,8,0),(0,8,0,24),(8,0,24,0),(0,24,0,80))
    assert all(sum((T**(i+j)).coeffs[k]*[4,0,8,0][k]for k in range(4))==tracegram[i][j]for i in range(4)for j in range(4))
    # Tr(v^2)=4|a+2c|²+8|c|²+8|b+3d|²+8|d|².
    assert tracegram==tuple(tuple(sum(w*row[i]*row[j]for w,row in[(4,(1,0,2,0)),(8,(0,0,1,0)),(8,(0,1,0,3)),(8,(0,0,0,1))])for j in range(4))for i in range(4))
    comparisons={}
    for name,d in [('unit_unimodular',-ALPHA),('ramified_unimodular',-BETA),('A9',-10*ALPHA),('D9',-4*ALPHA),('tripled_D9',-12*ALPHA),('tripled_E8_A1',-3*ALPHA)]:
        ratio=(-28*BETA)/d;assert square_root(ratio)is None;comparisons[name]=ratio.to_json()
    firstpath=ROOT/'results/fresh_h9_a6_d3_20261013.json';first=json.loads(firstpath.read_text())
    for path,h in first['sources'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
    oldwalls=[tuple(map(K.from_json,r))for r in first['combined_walls']];cap=oldwalls[-1]
    P=(cap[7]/2,*(Z for _ in range(6)),BETA*cap[0]/2,Z,Z);Q=STD[8]
    assert pair(P,G,P)==-D and pair(P,G,Q)==0 and pair(Q,G,Q)==1
    assert all((-pair(b,G,P)/D).is_integral and pair(b,G,Q).is_integral for b in BASIS)
    assert tuple(x/P[0]for x in P)==tuple(map(K.from_json,first['endpoint']))
    edge=[oldwalls[i]for i in[0,1,2,3,4,5,8,9]];C=tuple(tuple(pair(r,G,s)for s in edge)for r in edge);ci=inverse(C)
    profiles=[];rejected=[]
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
                co=mv(ci,tuple(-v for v in pre));r0=tuple(sum((c*r[j]for c,r in zip(co,edge)),Z)for j in range(10));N=q-pair(r0,G,r0)
                if any(N.sign(e)<=0 for e in range(4)):return
                assert pair(r0,G,P)==0 and pair(r0,G,Q)==0
                row={'q':q.to_json(),'magnitudes':enc(pre),'fixed':enc(r0),'target':N.to_json()}
                if any(not r0[i].is_integral for i in(*range(1,7),9)):
                    row['reason']='FIXED_COORDINATE_NONINTEGRAL';rejected.append(row);return
                if q not in DYADIC:
                    assert all(x==0 for x in pre)and all(x==0 for x in r0)
                    row['reason']='PRIME_SEVEN_FORCES_ALL_COORDINATES_DIVISIBLE';rejected.append(row);return
                if not N.is_integral:
                    row['reason']='NONINTEGRAL_RESIDUAL_NORM';rejected.append(row);return
                row['absolute_norm']=int(N.norm());profiles.append(row);return
            for p in alph[n]:walk((*pre,p))
        walk(())
    print('PROFILE_CENSUS',len(profiles),'rejected',len(rejected),flush=True)
    ns=sorted({r['absolute_norm']for r in profiles});gp=ROOT/'scripts'/f'{TAG}.gp';rawp=ROOT/'results'/f'{TAG}_raw.txt'
    if not verify:
        src=(ROOT/'scripts/fresh_h9_changed_spatial_20261013_d9.gp').read_text().replace('ns=vector(13,j,2^(j-1));','ns='+str(ns)+';');gp.write_text(src)
        with rawp.open('w')as f:
            proc=subprocess.run(['gp','-q',str(gp)],stdout=f,stderr=subprocess.PIPE,text=True,timeout=120);assert proc.returncode==0,proc.stderr
    raw=rawp.read_text()
    for s in('CERTIFY 1\n','BASE_CERTIFY 1\n','CLASS 1\n','SIGNATURE [2, 3]\n'):assert s in raw
    base=[];units=[];norms=[];absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('BASE_UNIT '):base.append(K.from_json(json.loads(line.split(' ',2)[2])))
        elif line.startswith('UNIT '):units.append(tuple(map(K.from_json,json.loads(line.split(' ',2)[2]))))
        elif line.startswith('UNIT_NORM '):
            n,e=line.split(' ',2)[2].split('] ',1);norms.append((K.from_json(json.loads(n+']')),json.loads(e)))
        elif line.startswith('NORM '):
            _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
        elif line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(tuple(map(K.from_json,json.loads(s))))
    assert len(base)==3 and len(units)==4 and set(absolute)==set(ns)and all(len(v)==counts[n]for n,v in absolute.items())
    for(z,x),(n,e)in zip(units,norms):assert z*z-D*x*x==n==K((-1)**e[3])*base[0]**e[0]*base[1]**e[1]*base[2]**e[2]
    m=[[row[1][i]for row in norms]for i in range(3)];kernel=[]
    for j in range(4):kernel.append((-1)**j*int(det([[K(m[i][k])for k in range(4)if k!=j]for i in range(3)]).a))
    gg=math.gcd(*kernel);assert gg;kernel=[n//gg for n in kernel]
    if sum(n*row[1][3]for n,row in zip(kernel,norms))%2:kernel=[2*n for n in kernel]
    def mul(u,v):return u[0]*v[0]+D*u[1]*v[1],u[0]*v[1]+u[1]*v[0]
    def power(u,n):
        if n<0:
            nu=u[0]*u[0]-D*u[1]*u[1];u=u[0]/nu,-u[1]/nu;n=-n
        out=(O,Z)
        while n:
            if n%2:out=mul(out,u)
            u=mul(u,u);n//=2
        return out
    eps=(O,Z)
    for u,n in zip(units,kernel):eps=mul(eps,power(u,n))
    if eps[0]<0:eps=tuple(-x for x in eps)
    if eps[1]<0:eps=eps[0],-eps[1]
    assert eps[0]>1 and eps[1]>0 and eps[0]**2-D*eps[1]**2==1
    back=eps[0],-eps[1]
    def root(r0,u):return tuple(r0[i]+u[1]*P[i]+u[0]*Q[i]for i in range(10))
    def param(u):return D*u[1]/u[0]
    targets={};starts=[];proposals=[]
    for row in profiles:
        N=K.from_json(row['target'])
        if N in targets:continue
        orbits=[];bad=[]
        for z,x in absolute[int(N.norm())]:
            ratio=(z*z-D*x*x)/N
            if not ratio.is_integral or not ratio.inverse().is_integral:bad.append(enc((z,x)));continue
            assert all(ratio.sign(e)>0 for e in range(4));t=square_root(ratio);assert t is not None and t.is_integral and t.inverse().is_integral
            z/=t;x/=t
            if z<0:z=-z;x=-x
            assert z*z-D*x*x==N;orbits.append({'representative':enc((z,x)),'unit_correction':t.to_json()})
        targets[N]={'target':N.to_json(),'orbits':orbits,'nonunit_rejections':bad}
    for i,row in enumerate(profiles):
        q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        for j,orb in enumerate(targets[N]['orbits']):
            u=tuple(map(K.from_json,orb['representative']));n=0
            while u[1]<=0:u=mul(u,eps);n+=1;assert abs(n)<50
            while mul(u,back)[1]>0:u=mul(u,back);n-=1;assert abs(n)<50
            starts.append((i,j,n,u))
            for step in range(40):
                r=root(r0,u)
                if reason(r,q)is None:proposals.append((param(u),r));break
                u=mul(u,eps)
    assert proposals;upper=min(p for p,r in proposals);visited=[];valid=[]
    for i,j,n,u in starts:
        row=profiles[i];q,N=K.from_json(row['q']),K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']));prev=mul(u,back);first=n;powers=[]
        while param(u)<=upper:
            r=root(r0,u);why=reason(r,q);rec={'profile':i,'orbit':j,'power':n,'root':enc(r),'q':q.to_json(),'parameter':param(u).to_json(),'reason':why};visited.append(rec);powers.append(n)
            if why is None:valid.append(rec)
            u=mul(u,eps);n+=1;assert n-first<100
        targets[N]['orbits'][j].setdefault('intervals',[]).append({'profile':i,'first_power':first,'last_power':n-1,'previous':enc(prev),'next':enc(u),'visited_powers':powers})
    assert len(valid)==1 and K.from_json(valid[0]['parameter'])==upper;r=tuple(map(K.from_json,valid[0]['root']));point=tuple(P[i]+upper*Q[i]for i in range(10));point=tuple(x/point[0]for x in point);walls=oldwalls+[r];frame=edge+[r]
    assert pair(point,G,point)<0 and all(pair(point,G,s)==0 for s in frame)and pair(point,G,INITIAL[7])<0
    for i,s in enumerate(walls):
        assert reason(s,pair(s,G,s))is None
        for t in walls[:i]:
            p=pair(s,G,t);assert p<=0;lam=4*p*p/(pair(s,G,s)*pair(t,G,t));assert lam in ALLOWED or lam>=4
    # Explicit strict center and an A6-dominated spacelike escape.
    a6=[row[1:7]for row in G[1:7]];rho=mv(inverse(a6),(O,)*6);delta=O
    while True:
        center=(O,*(delta*x for x in rho),3*delta,2*delta,delta)
        if pair(center,G,center)<0 and all(pair(center,G,s)<0 for s in walls):break
        delta/=2
    cbase=(O,*(Z for _ in range(6)),3*delta,2*delta,delta);S=delta*delta*pair((Z,*rho,Z,Z,Z),G,(Z,*rho,Z,Z,Z));C=-pair(cbase,G,cbase)/S;assert C>1;f=(1+C)/2
    escape=(O,*(f*delta*x for x in rho),3*delta,2*delta,delta)
    assert pair(escape,G,escape)>0 and pair(center,G,escape)<0 and all(pair(escape,G,s)<0 for s in walls)
    fg=tuple(tuple(pair(a,G,b)for b in frame)for a in frame);seen={tuple(O if i==j else Z for i in range(9))for j in range(9)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(9):
            c=2*sum((v[k]*fg[k][j]for k in range(9)),Z)/fg[j][j];w=list(v);w[j]-=c;w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=500
    assert all(tuple(-x for x in v)in seen for v in seen)
    out={'status':'PASS_COMPLETE_PRIME_SEVEN_A6_D3_SECOND_EDGE','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(fieldpath,firstpath,gp,rawp,Path(__file__).resolve(),ROOT/'scripts/fresh_h9_glued_host_20261013_arithmetic.py',Path(square_root.__code__.co_filename).resolve())},
         'host':'<-(T-1)(2+T)> + A6 + D3 over Z[T]','coefficient':BETA.to_json(),'ambient_gram':menc(G),'lattice_basis':menc(BASIS),'lattice_gram':menc(gram),'determinant':(-28*BETA).to_json(),'dual_exponent':28,'safe_norm_divisor':56,'safe_norm_pool':enc(POOL),'prime_seven_generators':enc(P7),'nonsquare_comparison_ratios':comparisons,
         'zero_chamber':'A6+B3','zero_mirrors':30,'trace_gram':tracegram,'zero_completeness':'Trace(v²)=4|a+2c|²+8|c|²+8|b+3d|²+8|d|². The rational lattice A6+D3 has minimum2, so every norm2 vector is rational. These roots span, hence every O-linear spatial isometry has rational integral entries and lies in Aut(A6) x Aut(D3). Its reflections are exactly W(A6) x W(B3).',
         'complement_time_vector':enc(P),'complement_spatial_vector':enc(Q),'residual_coefficient':D.to_json(),'profiles':profiles,'rejected_profiles':rejected,'edge':menc(edge),'base_units':enc(base),'fundamental_units':[enc(u)for u in units],'relative_norm_map':[e for n,e in norms],'relative_unit_kernel':kernel,'relative_unit':enc(eps),'absolute_norm_counts':counts,'targets':list(targets.values()),'proposal_power_limit':40,
         'upper_parameter':upper.to_json(),'visited_points':visited,'primitive_reflective_survivors':valid,'initial_walls':menc(INITIAL),'chosen_root':enc(r),'chosen_norm':pair(r,G,r).to_json(),'endpoint':enc(point),'endpoint_norm':pair(point,G,point).to_json(),'endpoint_mirrors':len(seen)//2,'endpoint_incident':menc(frame),'combined_walls':menc(walls),'strict_interior':enc(center),'strict_escape':enc(escape),'escape_norm':pair(escape,G,escape).to_json(),
         'scope':'Complete second finite edge in a determinant-squareclass-distinct quartic host; all 52 safe norms and finite profiles are screened. Nonintegral fixed coordinates and prime7 imprimitive profiles are rejected before the full certified relative-norm orbit interval. Full chamber reflectivity remains open; the prefix has an exact strict escape.'}
    if verify:assert json.loads(json.dumps(out))==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'profiles',len(profiles),'points',len(visited),'q',out['chosen_norm'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
