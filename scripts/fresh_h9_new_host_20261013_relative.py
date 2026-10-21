#!/usr/bin/env python3
"""Certified-CAS relative norm orbits on the open initial edge.

PARI supplies a certified complete absolute integral norm orbit system and a
fundamental unit basis. Everything after that finite interface is exact K
arithmetic, including normalization, exponent bounds, roots and angle gates.
"""
import hashlib,json,time
from collections import Counter
from fresh_h9_new_host_20261013_edge import *

EDGE=OUT
RAW=ROOT/'results/fresh_h9_new_host_20261013_relative_raw.txt'
PELL=ROOT/'results/fresh_h9_new_host_20261013_pell.json'
OUT=ROOT/'results/fresh_h9_new_host_20261013_relative.json'


def emul(u,v,D):return(u[0]*v[0]+D*u[1]*v[1],u[0]*v[1]+u[1]*v[0])

def run():
    started=time.monotonic();edge=json.loads(EDGE.read_text());data=json.loads(SOURCE.read_text());face=json.loads(FACE.read_text());pd=json.loads(PELL.read_text())
    raw=RAW.read_text();assert 'CERTIFY 1\n'in raw and 'CLASS 1\n'in raw and 'SIGNATURE [2, 3]\n'in raw
    units=[];absolute={};counts={}
    for line in raw.splitlines():
        if line.startswith('NORM '):
            _,n,c=line.split();counts[int(n)]=int(c);absolute[int(n)]=[]
        if not line.startswith(('UNIT ','SOLUTION ')):continue
        tag,n,encoded=line.split(' ',2);u=tuple(map(K.from_json,json.loads(encoded)))
        if tag=='UNIT':units.append(u)
        else:absolute[int(n)].append(u)
    assert len(units)==4 and all(len(absolute[n])==c for n,c in counts.items())
    D=6*ALPHA
    for a,b in units:assert abs((a*a-D*b*b).norm())==1
    assert all(b==0 for a,b in units[:3])
    assert units[3][0]**2-D*units[3][1]**2==1 and units[3][1]!=0
    A,B=units[3]
    if A<0:A=-A;B=-B
    if B<0:B=-B
    assert A>1 and B>0 and A*A-D*B*B==1
    epsilon=(A,B);epsilon_inverse=(A,-B)
    # Every fundamental unit except epsilon is in K. Hence conjugation fixes
    # their subgroup, while it sends epsilon to epsilon^-1. In a relative-norm
    # one unit their K-part has square1, leaving exactly +/- epsilon^n.
    g,fg=Host(2).g,Host(1).g;walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    fw=[tuple(map(K.from_json,r))for r in face['all_walls']];direction=tuple(map(K.from_json,edge['direction']))
    upper=K.from_json(pd['compatible_candidates'][0]['endpoint_parameter'])
    target_records={};rows=[];compatible_rows=[];stats=Counter()
    for i,profile in enumerate(edge['profiles']):
        q=K.from_json(profile['norm']);r0=tuple(map(K.from_json,profile['fixed_component']));R=K.from_json(profile['residual_norm']);N=6*R
        assert N.is_integral and all(N.sign(e)>0 for e in range(4));nn=N.norm();assert nn.denominator==1 and int(nn)in absolute
        key=str(N)
        if key not in target_records:
            reps=[];excluded=[]
            for a,b in absolute[int(nn)]:
                delta=(a*a-D*b*b)/N
                if not delta.is_integral or not delta.inverse().is_integral:
                    excluded.append({'absolute_representative':enc((a,b)),'norm_ratio':delta.to_json(),'reason':'NOT_A_UNIT'});continue
                assert all(delta.sign(e)>0 for e in range(4))
                root=square_root(delta);assert root is not None and root.is_integral and root.inverse().is_integral
                aa,bb=a/root,b/root
                if aa<0:aa=-aa;bb=-bb
                assert aa*aa-D*bb*bb==N and aa>0
                reps.append({'absolute_representative':enc((a,b)),'norm_ratio':delta.to_json(),'correction_unit':root.to_json(),'relative_representative':enc((aa,bb))})
            target_records[key]={'relative_target':N.to_json(),'absolute_norm':int(nn),'absolute_representative_count':len(absolute[int(nn)]),'relative_orbits':reps,'excluded':excluded}
        target=target_records[key]
        for orbit_index,rep in enumerate(target['relative_orbits']):
            u=tuple(map(K.from_json,rep['relative_representative']));exponent=0
            while u[1]<=0:
                u=emul(u,epsilon,D);exponent+=1;assert abs(exponent)<20
            while emul(u,epsilon_inverse,D)[1]>0:
                u=emul(u,epsilon_inverse,D);exponent-=1;assert abs(exponent)<20
            before=emul(u,epsilon_inverse,D);assert before[0]>0 and before[1]<=0 and u[0]>0 and u[1]>0
            first_exponent=exponent;visited=[]
            while ALPHA*u[1]<=upper*u[0]:
                z,x=u;assert z*z-D*x*x==N and z>0 and x>0
                r=tuple(a+z*b/6 for a,b in zip(r0,direction));r=(x,*r[1:]);assert pair(r,fg,r)==q
                stats['orbit_points']+=1
                row={'profile':i,'target':key,'orbit':orbit_index,'exponent':exponent,'z':z.to_json(),'x':x.to_json(),'vector_H5':enc(r),
                     'norm':q.to_json(),'endpoint_parameter':(ALPHA*x/z).to_json()}
                visited.append(exponent)
                if not all(a.is_integral for a in r):row['rejection']='NONINTEGRAL'
                elif all((a/T).is_integral for a in r):row['rejection']='IMPRIMITIVE'
                elif not all((2*a/q).is_integral for a in mv(fg,r)):row['rejection']='NONREFLECTIVE'
                else:
                    stats['primitive_reflective']+=1;r9=(r[0],*r[1:5],*(Z,)*4,r[5]);row['vector']=enc(r9)
                    if not compatible(r9,q,walls,g):row['rejection']='NOT_COMPATIBLE_WITH_PREFIX'
                    else:
                        point=tuple(a+(ALPHA*x/z)*b for a,b in zip((O,*(Z,)*5),direction))
                        assert pair(point,fg,point)<0 and all(pair(point,fg,w)<=0 for w in fw)
                        row.update(endpoint=enc(point),endpoint_norm=pair(point,fg,point).to_json(),height=(ALPHA*x*x/q).to_json())
                        compatible_rows.append(row)
                rows.append(row);u=emul(u,epsilon,D);exponent+=1;assert exponent-first_exponent<20
            assert u[1]>0 and ALPHA*u[1]>upper*u[0]
            rep.setdefault('profile_exponent_intervals',[]).append({'profile':i,'first':first_exponent,'last':exponent-1,
                'previous_pair':enc(before),'next_pair':enc(u),'visited':visited})
    compatible_rows.sort(key=lambda r:K.from_json(r['endpoint_parameter']))
    out={'status':'COMPLETE_RELATIVE_NORM_ORBIT_CENSUS_TO_PROPOSED_ENDPOINT','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(EDGE,SOURCE,FACE,PELL,RAW)},
         'PARI':{'absolute_polynomial':[ -2592,0,-1728,0,-288,0,0,0,1],'class_number':1,'certified':True,'signature':[2,3],
             'fundamental_units':[enc(u)for u in units],'absolute_norm_counts':counts,
             'dependency':'bnfcertify=1 and bnfisintnorm complete integral solutions modulo positive absolute-norm units; this CAS completeness is a stated dependency, not replaced by a numerical norm solver.'},
         'discriminant':D.to_json(),'relative_unit':enc(epsilon),'upper_endpoint_parameter':upper.to_json(),
         'targets':list(target_records.values()),'orbit_points':rows,'compatible_candidates':compatible_rows,
         'nearest_parameter':compatible_rows[0]['endpoint_parameter']if compatible_rows else None,
         'nearest_candidates':[r for r in compatible_rows if r['endpoint_parameter']==compatible_rows[0]['endpoint_parameter']]if compatible_rows else[],
         'counts':dict(stats),'seconds':round(time.monotonic()-started,3),
         'scope':'Complete relative-norm solutions for all30 finite-coupling profiles on the chosen initial H5 edge, within the proposed endpoint interval, conditional only on the archived certified PARI ideal/unit interface. All residual square classes, primitive norms and unit exponents are included. Exact arithmetic checks integrality and every H9 prefix angle. A global first-mirror interpretation additionally uses the complete initial chamber and finite-angle necessity; a separate local-stabilizer check is not asserted here.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print({'status':out['status'],'targets':len(target_records),'orbit_points':len(rows),'primitive':stats['primitive_reflective'],'compatible':len(compatible_rows),'nearest_count':len(out['nearest_candidates']),'seconds':out['seconds']},flush=True)
    if compatible_rows:print({'nearest_height':K.from_json(compatible_rows[0]['height']).real(),'nearest':compatible_rows[0]},flush=True)

if __name__=='__main__':run()
