#!/usr/bin/env python3
"""Separate Fraction replay of all four next A9 first-endpoint interfaces."""
import argparse,hashlib,json
from fresh_h9_new_host_next_20261013_a9_independent import a,ALPHA,G,POOL,P,compatible
from fresh_h9_quartic_first_endpoint_20261013_independent import quartic_roots,closure
from fresh_h9_glued_host_20261013_structure import norm
ROOT=a.ROOT
PROFILE=ROOT/'results/fresh_h9_a9_next_20261013_profiles.json'
ORBIT=ROOT/'results/fresh_h9_a9_next_20261013_orbits.json'
OUT=ROOT/'results/fresh_h9_a9_next_20261013_independent.json'

def power(x,n):
    if n<0:return power(a.inv(x),-n)
    z=a.O
    while n:
        if n%2:z=a.mul(z,x)
        x=a.sq(x);n//=2
    return z

def emul(u,v,D):return(a.add(a.mul(u[0],v[0]),a.mul(D,a.mul(u[1],v[1]))),a.add(a.mul(u[0],v[1]),a.mul(u[1],v[0])))
def positive(rep,eps,D):
    u=rep;n=0;back=(eps[0],a.neg(eps[1]))
    while a.sign(u[1])<=0:u=emul(u,eps,D);n+=1;assert abs(n)<100
    while a.sign(emul(u,back,D)[1])>0:u=emul(u,back,D);n-=1;assert abs(n)<100
    return u,n,emul(u,back,D)

def load_interface(kind,pd,record):
    for p,h in record['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    raw=(ROOT/f'results/fresh_h9_a9_next_20261013_{kind}_raw.txt').read_text()
    assert 'BASE_CERTIFY 1\n'in raw and'CERTIFY 1\n'in raw and'CLASS 1\n'in raw and'SIGNATURE [2, 3]\n'in raw
    expected='x^8 + 4*x^6 + 2*x^4 - 4*x^2 - 1'if kind=='end'else'x^8 + 24*x^6 + 72*x^4 - 864*x^2 - 1296'
    assert 'ABSOLUTE '+expected+'\n'in raw
    D=ALPHA if kind=='end'else a.scale(ALPHA,6);base=[];units=[];maps=[];absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('BASE_UNIT '):base.append(a.dec(json.loads(line.split(' ',2)[2])))
        elif line.startswith('UNIT '):
            left,right=line.split(' MAP ');units.append(tuple(a.dec(v)for v in json.loads(left.split(' ',2)[2])));maps.append(json.loads(right))
        elif line.startswith('NORM '):
            _,n,c=line.split();counts[int(n)]=int(c);absolute[int(n)]=[]
        elif line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(tuple(a.dec(v)for v in json.loads(s)))
    assert len(base)==3 and len(units)==4 and all(abs(norm(u))==1 for u in base)
    for (u,v),exps in zip(units,maps,strict=True):
        expected=a.scale(a.O,(-1)**exps[3])
        for x,n in zip(base,exps[:3],strict=True):expected=a.mul(expected,power(x,n))
        assert a.sub(a.sq(u),a.mul(D,a.sq(v)))==expected
    m=[[maps[j][i]for j in range(3)]for i in range(3)];det=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
    assert det==record['map_minor_determinant']and det!=0 and maps[3]==[0,0,0,0]
    eps=units[3]
    if a.sign(eps[0])<0:eps=tuple(a.neg(v)for v in eps)
    if a.sign(eps[1])<0:eps=(eps[0],a.neg(eps[1]))
    assert eps==a.vec(record['relative_unit'])and a.sub(a.sq(eps[0]),a.mul(D,a.sq(eps[1])))==a.O
    assert a.sign(a.sub(eps[0],a.O))>0 and a.sign(eps[1])>0
    targets={a.dec(row['relative_target'])for case in pd['cases']if case['kind']==kind for row in case['profiles']}
    assert set(absolute)=={int(norm(N))for N in targets}and all(len(absolute[n])==counts[n]for n in absolute)
    reps={}
    for N in targets:
        reps[N]=[]
        for z,x in absolute[int(norm(N))]:
            n=a.sub(a.sq(z),a.mul(D,a.sq(x)));assert norm(n)==norm(N);delta=a.divide(n,N)
            if not a.integral(delta)or not a.integral(a.inv(delta)):continue
            assert all(a.sign(delta,e)>0 for e in range(4));ss=[s for s in quartic_roots(delta)if a.sign(s)>0];assert len(ss)==1;s=ss[0]
            assert a.integral(s)and a.integral(a.inv(s));z=a.divide(z,s);x=a.divide(x,s)
            if a.sign(z)<0:z=a.neg(z);x=a.neg(x)
            assert a.sub(a.sq(z),a.mul(D,a.sq(x)))==N;reps[N].append((z,x))
    saved={a.dec(row['target']):row for row in record['targets'].values()};assert targets==set(saved)
    for N,rs in reps.items():assert rs==[a.vec(row['relative_representative'])for row in saved[N]['relative_orbits']]
    return D,eps,reps

def run(verify=False):
    pd,od=[json.loads(p.read_text())for p in(PROFILE,ORBIT)]
    assert od['status']=='COMPLETE_FOUR_EDGE_NORM_INTERFACE_ROUND'
    for d in(pd,od):
        for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    source=ROOT/'results/fresh_h9_new_host_next_20261013_a9_endpoint.json';sd=json.loads(source.read_text())
    firstpath=ROOT/'results/fresh_h9_new_host_next_20261013_a9_first_independent.json';prior=json.loads(firstpath.read_text())
    assert prior['status']=='PASS_COMPLETE_A9_FIRST_ENDPOINT_WITH_CERTIFIED_PARI_INTERFACE'and prior['sources'][str(source.relative_to(ROOT))]==hashlib.sha256(source.read_bytes()).hexdigest()
    walls=[a.vec(r)for r in sd['all_walls']];inside=a.vec(sd['strict_interior']);corner=sd['actual_finite_vertex']['indices']
    start=a.vec(sd['endpoint']);z=a.dec(sd['z']);V=tuple(a.mul(z,v)for v in start);assert P(V,V)==a.scale(ALPHA,-20)
    assert a.vec(pd['norm_pool'])==POOL
    alphabets=[]
    for q,saved in zip(POOL,pd['alphabets'],strict=True):
        values={a.Z}
        for angle in a.ALLOWED:
            for v in quartic_roots(a.scale(a.mul(q,angle),a.Q(1,2))):
                if a.sign(v)>=0 and a.integral(v)and a.integral(a.divide(a.scale(v,2),q)):values.add(v)
        assert values==set(a.vec(saved));alphabets.append(values)
    tables={kind:load_interface(kind,pd,od['interfaces'][kind])for kind in('end','inner')};results=[]
    class Host:pair=staticmethod(P)
    host=Host();assert len(closure([walls[i]for i in corner],[walls[i]for i in corner],host))==21
    for case,saved in zip(pd['cases'],od['cases'],strict=True):
        ci=saved['case'];released=case['released_wall'];ids=case['edge_indices'];edge=[walls[i]for i in ids]
        assert ids==[i for i in corner if i!=released]and edge==[a.vec(r)for r in case['edge_roots']]
        assert start==a.vec(case['start'])and V==a.vec(case['timelike_vector'])
        gram=tuple(tuple(P(r,s)for s in edge)for r in edge);gi=a.minverse(gram);a.positive_definite(gram)
        assert gram==tuple(a.vec(row)for row in case['edge_gram'])
        co=a.mv(gi,tuple(P(walls[released],r)for r in edge));b=tuple(a.sub(walls[released][j],a.dot(co,tuple(r[j]for r in edge)))for j in range(10));B=P(b,b)
        assert b==a.vec(case['released_residual'])and B==a.dec(case['residual_square'])and P(V,b)==a.Z
        den=case['timelike_denominator'];mult=case['target_multiplier'];D,eps,reps=tables[case['kind']]
        assert D==a.dec(case['relative_discriminant'])and D==a.divide(a.scale(ALPHA,500),a.scale(B,den*den))
        assert a.divide(a.scale(a.O,25),B)==a.scale(a.O,mult)
        zcov=tuple(a.divide(a.scale(v,-5),B)for v in a.mv(G,b));xcov=tuple(a.divide(a.scale(v,-den),a.scale(ALPHA,20))for v in a.mv(G,V))
        assert zcov==a.vec(case['integral_Z_covector'])and xcov==a.vec(case['integral_X_covector'])and all(map(a.integral,zcov+xcov))
        # Independent fixed component: solve all ten covector equations.
        basis=edge+[b,V];eq_inverse=a.minverse(tuple(a.mv(G,r)for r in basis))
        profiles={(a.dec(row['norm']),a.vec(row['magnitudes'])):row for row in case['profiles']};found=set();order=list(reversed(range(8)));Q=tuple(tuple(gi[i][j]for j in order)for i in order)
        for q,alphabet in zip(POOL,alphabets,strict=True):
            def walk(pre):
                k=len(pre)
                if k and a.sign(a.sub(q,a.dot(pre,a.mv(tuple(row[:k]for row in Q[:k]),pre))))<=0:return
                if k==8:
                    mag=tuple(reversed(pre));r0=a.mv(eq_inverse,(*map(a.neg,mag),a.Z,a.Z));R=a.sub(q,P(r0,r0))
                    if any(a.sign(R,e)<=0 for e in range(4)):return
                    key=q,mag;found.add(key);row=profiles[key]
                    assert r0==a.vec(row['fixed_component'])and R==a.dec(row['residual'])and a.scale(R,mult)==a.dec(row['relative_target'])
                    assert norm(a.scale(R,mult))==row['absolute_norm'];return
                for v in alphabet:walk((*pre,v))
            walk(())
        assert found==set(profiles)
        assert saved['status']=='COMPLETE_FIRST_ENDPOINT_WITH_CERTIFIED_NORM_INTERFACE'
        upper=a.dec(saved['upper_parameter']);visited=[];intervals=[];accepted=[];primitives=[]
        for pi,row in enumerate(case['profiles']):
            q=a.dec(row['norm']);r0=a.vec(row['fixed_component']);N=a.dec(row['relative_target'])
            for oi,rep in enumerate(reps[N]):
                u,power,previous=positive(rep,eps,D);first=power;powers=[]
                while True:
                    Zc,Xc=u;s=a.divide(a.scale(a.mul(ALPHA,Xc),100),a.scale(a.mul(B,Zc),den))
                    if a.sign(a.sub(s,upper))>0:break
                    r=tuple(a.add(a.sub(v,a.scale(a.mul(Zc,bv),a.Q(1,5))),a.scale(a.mul(Xc,w),a.Q(1,den)))for v,bv,w in zip(r0,b,V))
                    assert P(r,r)==q and a.dot(zcov,r)==Zc and a.dot(xcov,r)==Xc
                    reason=None
                    if not all(map(a.integral,r)):reason='NONINTEGRAL'
                    elif all(a.integral(a.divide(v,a.T))for v in r)or all(a.integral(a.scale(v,a.Q(1,5)))for v in r):reason='IMPRIMITIVE'
                    elif not all(a.integral(a.divide(a.scale(v,2),q))for v in a.mv(G,r)):reason='NONREFLECTIVE'
                    elif not compatible(r,walls):reason='NOT_PREFIX_COMPATIBLE'
                    if reason is None:accepted.append((s,r,q))
                    if reason in(None,'NOT_PREFIX_COMPATIBLE'):primitives.append((s,r,q))
                    visited.append((pi,oi,power,q,r,u,s,reason));powers.append(power);u=emul(u,eps,D);power+=1;assert power-first<1000
                intervals.append({'profile':pi,'orbit':oi,'first_power':first,'last_power':power-1,'previous':[a.enc(v)for v in previous],'next':[a.enc(v)for v in u],'visited_powers':powers})
        expected=[(row['profile'],row['orbit'],row['power'],a.dec(row['norm']),a.vec(row['vector']),a.vec(row['Z_X']),a.dec(row['parameter']),row.get('rejection'))for row in saved['visited_points']]
        assert visited==expected and intervals==saved['intervals']
        minimum=a.dec(saved['nearest_parameter']);assert all(a.sign(a.sub(s,minimum))>=0 for s,r,q in accepted+primitives)
        earliest=[r for s,r,q in accepted if s==minimum];assert set(earliest)==set(a.vec(row['vector'])for row in saved['nearest_candidates'])
        allgroup=closure(edge+earliest,edge+earliest,host);chosen=a.vec(saved['selected']['vector']);assert chosen in earliest
        assert set(allgroup)==set(closure(edge+[chosen],edge+[chosen],host))and len(allgroup)==saved['full_endpoint_mirrors']
        point=tuple(a.sub(v,a.mul(minimum,bv))for v,bv in zip(V,b));point=tuple(a.divide(v,point[0])for v in point)
        assert point==a.vec(saved['selected']['point'])and P(point,point)==a.dec(saved['selected']['point_norm'])and a.sign(P(point,point))<0
        assert a.sign(P(start,chosen))<0 and a.sign(P(inside,chosen))<0 and all(a.sign(P(point,r))<=0 for r in walls)
        assert all(P(point,r)==a.Z for r in edge+[chosen])and compatible(chosen,walls)
        a.positive_definite(tuple(tuple(P(r,s)for s in edge+[chosen])for r in edge+[chosen]))
        results.append({'case':ci,'released_wall':released,'finite_profiles':len(profiles),'relative_targets':len({a.dec(row['relative_target'])for row in case['profiles']}),
          'interval_points':len(visited),'primitive_mirrors':len(primitives),'nearest_mirrors':len(earliest),'full_endpoint_mirrors':len(allgroup),
          'source_corner_mirrors':21,'exact_orientation_and_sector_join':True,'root':list(map(a.enc,chosen))})
        print('independent case',ci,'PASS',len(profiles),len(visited),len(primitives),flush=True)
    result={'status':'PASS_COMPLETE_FOUR_A9_FIRST_ENDPOINTS','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(PROFILE,ORBIT,source,firstpath)},
      'cases':results,'arithmetic':'Separate Fraction field arithmetic; reversed complete finite-profile enumeration; ten-equation fixed-component solution; independent unit-norm map and every relative orbit interval; complete finite reflection closures and oriented sector joins.',
      'completeness_dependency':'The archived certified PARI absolute norm representatives and full unit bases are a stated common CAS interface. The independent replay does not claim a second ideal/unit algorithm.'}
    if verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
