#!/usr/bin/env python3
"""Separate arithmetic replay of the complete A9 finite-endpoint interface.

The certified PARI absolute norm representatives and fundamental unit basis
are the sole CAS completeness dependency. No claim of a CAS-free exhaustive
coefficient search is made.
"""
import argparse,hashlib,json
from fresh_h9_new_host_next_20261013_a9_independent import a,ALPHA,A,G,POOL,P,compatible
from fresh_h9_quartic_first_endpoint_20261013_independent import quartic_roots,closure
from fresh_h9_glued_host_20261013_structure import norm

ROOT=a.ROOT
PREFIX='fresh_h9_new_host_next_20261013_a9_'
PROFILE=ROOT/'results'/f'{PREFIX}profiles.json'
ORBIT=ROOT/'results'/f'{PREFIX}orbits.json'
CAP=ROOT/'results'/f'{PREFIX}endpoint.json'
RAW=ROOT/'results'/f'{PREFIX}orbits_raw.txt'
OUT=ROOT/'results'/f'{PREFIX}first_independent.json'
D=a.scale(ALPHA,10)
emul=lambda u,v:(a.add(a.mul(u[0],v[0]),a.mul(D,a.mul(u[1],v[1]))),a.add(a.mul(u[0],v[1]),a.mul(u[1],v[0])))

def run(verify=False):
    pd,od,cd=[json.loads(p.read_text())for p in(PROFILE,ORBIT,CAP)];raw=RAW.read_text()
    for p,h in od['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    ids=[i for i in range(9)if i!=4];edge=[tuple(a.neg(a.O)if j==i+1 else a.Z for j in range(10))for i in ids]
    gram=tuple(tuple(P(r,s)for s in edge)for r in edge);gi=a.minverse(gram)
    assert [a.vec(r)for r in pd['edge_roots']]==edge and a.vec(pd['norm_pool'])==POOL
    saved={(a.dec(row['norm']),a.vec(row['magnitudes'])):i for i,row in enumerate(pd['profiles'])};found=set()
    for q,ap in zip(POOL,pd['alphabets'],strict=True):
        alphabet={a.Z}
        for value in a.ALLOWED:
            for v in quartic_roots(a.scale(a.mul(q,value),a.Q(1,2))):
                if a.sign(v)>=0 and a.integral(v)and a.integral(a.divide(a.scale(v,2),q)):alphabet.add(v)
        assert alphabet==set(a.vec(ap['magnitudes']))and a.dec(ap['norm'])==q
        # Reverse order, and invert the two Cartan blocks to reconstruct the
        # unique fixed spatial component whose middle coordinate is zero.
        rev=list(reversed(range(8)));Q=tuple(tuple(gi[i][j]for j in rev)for i in rev)
        def walk(pre):
            k=len(pre)
            if k and a.sign(a.sub(q,a.dot(pre,a.mv(tuple(row[:k]for row in Q[:k]),pre))))<=0:return
            if k==8:
                mags=tuple(reversed(pre));coeff=a.mv(gi,mags);r0=[a.Z]*10
                for j,v in zip(ids,coeff):r0[j+1]=v
                r0=tuple(r0);R=a.sub(q,P(r0,r0))
                if any(a.sign(R,e)<=0 for e in range(4)):return
                key=q,mags;found.add(key);row=pd['profiles'][saved[key]]
                assert r0==a.vec(row['fixed_component'])and R==a.dec(row['residual'])
                assert a.scale(R,10)==a.dec(row['relative_target'])and norm(a.scale(R,10))==row['absolute_norm'];return
            for v in alphabet:walk((*pre,v))
        walk(())
    assert found==set(saved)and len(found)==50
    assert 'CERTIFY 1\n'in raw and'CLASS 2\n'in raw and'SIGNATURE [2, 3]\n'in raw
    assert 'ABSOLUTE x^8 + 40*x^6 + 200*x^4 - 4000*x^2 - 10000\n'in raw
    units=[];absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('NORM '):
            _,n,c=line.split();counts[int(n)]=int(c);absolute[int(n)]=[]
        elif line.startswith('UNIT '):units.append(tuple(a.dec(v)for v in json.loads(line.split(' ',2)[2])))
        elif line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(tuple(a.dec(v)for v in json.loads(s)))
    assert set(absolute)==set(pd['absolute_norms'])and all(len(absolute[n])==counts[n]for n in absolute)
    for n,reps in absolute.items():
        for z,x in reps:assert norm(a.sub(a.sq(z),a.mul(D,a.sq(x))))==n
    assert len(units)==4 and units[0]==(ALPHA,a.Z)and all(v==a.Z for u,v in units[:3])
    assert all(abs(norm(u))==1 for u,v in units[:3])
    u,v=units[3];assert a.sub(a.sq(u),a.mul(D,a.sq(v)))==a.neg(ALPHA)
    epsilon=tuple(a.divide(v,ALPHA)for v in emul(units[3],units[3]))
    if a.sign(epsilon[0])<0:epsilon=tuple(a.neg(v)for v in epsilon)
    if a.sign(epsilon[1])<0:epsilon=(epsilon[0],a.neg(epsilon[1]))
    assert a.sub(a.sq(epsilon[0]),a.mul(D,a.sq(epsilon[1])))==a.O and a.sign(a.sub(epsilon[0],a.O))>0 and a.sign(epsilon[1])>0
    assert epsilon==a.vec(od['relative_unit']);back=(epsilon[0],a.neg(epsilon[1]))
    upper=a.dec(od['upper_parameter']);assert upper==a.divide(a.mul(ALPHA,a.dec(cd['time'])),a.dec(cd['z']))
    direction=(a.Z,*(a.scale(a.O,i)for i in(1,2,3,4,5,4,3,2,1)))
    visited=[];reflective=[];targets={}
    for i,row in enumerate(pd['profiles']):
        q=a.dec(row['norm']);r0=a.vec(row['fixed_component']);N=a.dec(row['relative_target']);reps=[]
        for z,x in absolute[int(norm(N))]:
            delta=a.divide(a.sub(a.sq(z),a.mul(D,a.sq(x))),N)
            if not a.integral(delta)or not a.integral(a.inv(delta)):continue
            assert all(a.sign(delta,e)>0 for e in range(4))
            roots=[s for s in quartic_roots(delta)if a.sign(s)>0];assert len(roots)==1;s=roots[0]
            assert a.integral(s)and a.integral(a.inv(s));z=a.divide(z,s);x=a.divide(x,s)
            if a.sign(z)<0:z=a.neg(z);x=a.neg(x)
            reps.append((z,x))
        targets[N]=reps
        for j,representative in enumerate(reps):
            u=representative;power=0
            while a.sign(u[1])<=0:u=emul(u,epsilon);power+=1;assert abs(power)<20
            while a.sign(emul(u,back)[1])>0:u=emul(u,back);power-=1;assert abs(power)<20
            assert a.sign(emul(u,back)[1])<=0 and a.sign(u[1])>0
            while a.sign(a.sub(a.mul(upper,u[0]),a.mul(ALPHA,u[1])))>=0:
                z,x=u;r=tuple(a.add(v,a.scale(a.mul(z,d),a.Q(1,10)))for v,d in zip(r0,direction));r=(x,*r[1:])
                assert P(r,r)==q;reason=None
                if not all(map(a.integral,r)):reason='NONINTEGRAL'
                elif all(a.integral(a.divide(v,a.T))for v in r)or all(a.integral(a.scale(v,a.Q(1,5)))for v in r):reason='IMPRIMITIVE'
                elif not all(a.integral(a.divide(a.scale(v,2),q))for v in a.mv(G,r)):reason='NONREFLECTIVE'
                else:reflective.append(r)
                visited.append((i,j,power,r,reason));u=emul(u,epsilon);power+=1;assert abs(power)<20
            assert a.sign(a.sub(a.mul(ALPHA,u[1]),a.mul(upper,u[0])))>0
    expected=[(row['profile'],row['orbit'],row['power'],a.vec(row['vector']),row.get('rejection'))for row in od['visited_points']]
    assert visited==expected and len(visited)==5 and reflective==[a.vec(cd['root'])]
    # No prefix-compatibility rejection is used to reach this uniqueness.
    initial=[tuple(a.neg(a.O)if i==j else a.Z for i in range(10))for j in range(1,10)]
    class Host:
        pair=staticmethod(P)
    host=Host();end=edge+reflective
    assert [len(closure(g,g,host))for g in(initial,edge,end)]==[45,20,21]
    assert all(compatible(r,end[:i])for i,r in enumerate(end))
    a.positive_definite(tuple(tuple(P(r,s)for s in end)for r in end))
    point=a.vec(cd['endpoint']);assert all(P(point,r)==a.Z for r in end)and a.sign(P(point,point))<0
    result={'status':'PASS_COMPLETE_A9_FIRST_ENDPOINT_WITH_CERTIFIED_PARI_INTERFACE',
      'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(PROFILE,ORBIT,CAP,RAW)},
      'complete_finite_profiles':50,'relative_targets':len(targets),'absolute_norm_targets':len(absolute),
      'orbit_points_in_interval':len(visited),'primitive_reflective_mirrors_in_interval':1,
      'mirror_counts_initial_edge_endpoint':[45,20,21],'endpoint_type':'A4+A4+A1',
      'nearest_parameter':a.enc(upper),'nearest_root':list(map(a.enc,reflective[0])),
      'arithmetic':'Separate Fraction operations, reversed finite-profile census and direct Cartan-block reconstruction; independent normalization and complete unit-interval replay. Absolute norm representatives and full units remain an explicitly certified PARI dependency.',
      'firstness_argument':'The initial full A9 chamber and its middle A4+A4 edge are complete. The explicit cap bounds a finite segment. Any earlier first mirror of the full root arrangement gives a finite vertex; local finiteness makes the spherical vertex chamber simplicial, so its ninth facet has Coxeter angles to all eight existing edge facets and belongs to the enumerated profiles. Its integral z+x*sqrt(10alpha) lies in one of the complete absolute norm orbits, hence in the replayed interval. The only primitive reflective mirror there is the displayed cap. Thus the cap is the first facet, and its nine simple incident walls generate the full A4+A4+A1 stabilizer. This relies on the certified CAS orbit interface, not a claimed independent coordinate census.'}
    if verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],len(visited),len(reflective),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
