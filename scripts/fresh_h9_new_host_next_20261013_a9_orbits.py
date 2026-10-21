#!/usr/bin/env python3
"""All relative norm orbits of the finite middle-edge coupling profiles.

Completeness is conditional on the explicitly archived certified PARI
bnfisintnorm/fundamental-unit interface. Root and geometric identities use
exact field arithmetic. This is not a bounded unit-grid search.
"""
import hashlib
from fresh_h9_new_host_next_20261013_a9 import *
from fresh_h9_quartic_mixed_20261013_endpoint import square_root
PROFILE=ROOT/'results/fresh_h9_new_host_next_20261013_a9_profiles.json'
RAW=ROOT/'results/fresh_h9_new_host_next_20261013_a9_orbits_raw.txt'
CAP=ROOT/'results/fresh_h9_new_host_next_20261013_a9_endpoint.json'
SOURCE=ROOT/'results/fresh_h9_new_host_next_20261013_a9_h650.json'
OUT=ROOT/'results/fresh_h9_new_host_next_20261013_a9_orbits.json'
GP=ROOT/'scripts/fresh_h9_new_host_next_20261013_a9_orbits.gp'
D=10*ALPHA
emul=lambda u,v:(u[0]*v[0]+D*u[1]*v[1],u[0]*v[1]+u[1]*v[0])

def run():
    pd=json.loads(PROFILE.read_text());cd=json.loads(CAP.read_text());sd=json.loads(SOURCE.read_text());raw=RAW.read_text()
    assert 'CERTIFY 1\n'in raw and'CLASS 2\n'in raw and'SIGNATURE [2, 3]\n'in raw
    assert 'ABSOLUTE x^8 + 40*x^6 + 200*x^4 - 4000*x^2 - 10000\n'in raw
    units=[];absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('NORM '):
            _,n,c=line.split();absolute[int(n)]=[];counts[int(n)]=int(c)
        if line.startswith('UNIT '):units.append(tuple(map(K.from_json,json.loads(line.split(' ',2)[2]))))
        if line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absolute[int(n)].append(tuple(map(K.from_json,json.loads(s))))
    assert set(absolute)==set(pd['absolute_norms'])and all(len(absolute[n])==counts[n]for n in absolute)
    assert len(units)==4 and all(v==0 for u,v in units[:3])and units[0]==(ALPHA,Z)
    assert all(abs(u.norm())==1 for u,v in units[:3])
    u,v=units[3];assert u*u-D*v*v==-ALPHA
    epsilon=tuple(z/ALPHA for z in emul(units[3],units[3]));assert epsilon[0]**2-D*epsilon[1]**2==1
    if epsilon[0]<0:epsilon=tuple(-z for z in epsilon)
    if epsilon[1]<0:epsilon=(epsilon[0],-epsilon[1])
    assert epsilon[0]>1 and epsilon[1]>0;epsinv=(epsilon[0],-epsilon[1])
    direction=(Z,*(K(i)for i in(1,2,3,4,5,4,3,2,1)));upper=ALPHA*K.from_json(cd['time'])/K.from_json(cd['z'])
    walls=[tuple(map(K.from_json,r))for r in sd['all_walls']];targets={};visited=[];accepted=[]
    for i,profile in enumerate(pd['profiles']):
        q=K.from_json(profile['norm']);r0=tuple(map(K.from_json,profile['fixed_component']));N=K.from_json(profile['relative_target'])
        key=str(N)
        if key not in targets:
            reps=[];bad=[]
            for z,x in absolute[int(N.norm())]:
                delta=(z*z-D*x*x)/N
                if not delta.is_integral or not delta.inverse().is_integral:
                    bad.append({'representative':enc((z,x)),'ratio':delta.to_json()});continue
                assert all(delta.sign(e)>0 for e in range(4));s=square_root(delta)
                assert s is not None and s.is_integral and s.inverse().is_integral
                z/=s;x/=s
                if z<0:z=-z;x=-x
                assert z*z-D*x*x==N and z>0
                reps.append({'representative':enc((z,x)),'correction_unit':s.to_json(),'ratio':delta.to_json()})
            targets[key]={'target':N.to_json(),'absolute_norm':int(N.norm()),'orbits':reps,'excluded_nonunit_ratios':bad}
        for j,rep in enumerate(targets[key]['orbits']):
            u=tuple(map(K.from_json,rep['representative']));power=0
            while u[1]<=0:u=emul(u,epsilon);power+=1;assert abs(power)<20
            while emul(u,epsinv)[1]>0:u=emul(u,epsinv);power-=1;assert abs(power)<20
            previous=emul(u,epsinv);assert previous[1]<=0 and u[1]>0 and u[0]>0
            first=power;points=[]
            while ALPHA*u[1]<=upper*u[0]:
                z,x=u;r=tuple(v+z*d/10 for v,d in zip(r0,direction));r=(x,*r[1:])
                assert pair(r,G,r)==q and z*z-D*x*x==N
                row={'profile':i,'orbit':j,'power':power,'vector':enc(r),'norm':q.to_json(),'z':z.to_json(),'time':x.to_json(),'parameter':(ALPHA*x/z).to_json()}
                if not all(v.is_integral for v in r):row['rejection']='NONINTEGRAL'
                elif all((v/T).is_integral for v in r)or all((v/5).is_integral for v in r):row['rejection']='IMPRIMITIVE'
                elif not all((2*v/q).is_integral for v in mv(G,r)):row['rejection']='NONREFLECTIVE'
                elif not compatible(r,q,walls,G):row['rejection']='NOT_PREFIX_COMPATIBLE'
                else:
                    point=(O,*(ALPHA*x/z*d for d in direction[1:]));assert pair(point,G,r)==0 and pair(point,G,point)<0
                    assert all(pair(point,G,w)<=0 for w in walls)
                    row.update(endpoint=enc(point),endpoint_norm=pair(point,G,point).to_json(),height=(ALPHA*x*x/q).to_json());accepted.append(row)
                visited.append(row);points.append(power);u=emul(u,epsilon);power+=1;assert power-first<20
            rep.setdefault('intervals',[]).append({'profile':i,'first_power':first,'last_power':power-1,'previous':enc(previous),'next':enc(u),'visited':points})
    accepted.sort(key=lambda row:K.from_json(row['parameter']));assert accepted
    earliest=[r for r in accepted if r['parameter']==accepted[0]['parameter']]
    out={'status':'COMPLETE_RELATIVE_NORM_FINITE_ENDPOINT_CENSUS','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(PROFILE,RAW,CAP,SOURCE,GP)},
      'relative_discriminant':D.to_json(),'relative_unit':enc(epsilon),'fundamental_units':[enc(u)for u in units],
      'unit_kernel_proof':'The first three fundamental units lie inK and u1=alpha. Norm(u4)=−alpha. A relative norm-one unit must have even u4 exponent by its primary sign. Writing that exponent2k leaves the base-field factor exactly±alpha^(−k), so the kernel is generated by±epsilon with epsilon=u4²/alpha.',
      'PARI_dependency':'bnfcertify=1; bnfisintnorm supplies all integral absolute norm representatives modulo norm-one units. Class number2 is handled by that solver; principality is not assumed. Completeness of this finite interface is an explicit CAS dependency.',
      'absolute_norm_counts':counts,'targets':list(targets.values()),'upper_parameter':upper.to_json(),
      'visited_points':visited,'accepted':accepted,'nearest_candidates':earliest,'nearest_parameter':earliest[0]['parameter'],
      'scope':'Complete relative norm orbit census for all saved finite Coxeter coupling profiles up to the explicit cap. Exact normalization, full unit kernel, interval end exclusions, integrality and all prefix angle gates are checked. This conditional first-facet interface is stronger than a bounded exponent scout; independent replay and complete incident-stabilizer identification are separate.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],'profiles',len(pd['profiles']),'points',len(visited),'accepted',len(accepted),'nearest',len(earliest),'nearest height',K.from_json(earliest[0]['height']).real(),flush=True)

if __name__=='__main__':run()
