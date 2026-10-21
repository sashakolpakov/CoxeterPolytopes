#!/usr/bin/env python3
"""Finite edge-coupling profiles with explicit bounded relative-unit orbits."""
from fresh_h9_quartic_edge_walk_20261013 import *

OUT=ROOT/'results/fresh_h9_quartic_edge_walk_20261013_profiles.json'


def run_profiles():
    started=time.monotonic();data=json.loads(SOURCE.read_text());g=Host(2).g;fg=Host(1).g
    walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    face=[tuple(map(K.from_json,r))for r in data['forced_H5_face']['all_walls']]
    edge_ids=[1,2,4,12];edge=[face[i]for i in edge_ids]
    eg=tuple(tuple(pair(r,fg,s)for s in edge)for r in edge);ei=inverse(eg)
    h,hi,j,p,_,_=edge_translation(face[12],fg)
    base=tuple(x+SQRT2*y/2 for x,y in zip(face[6],face[12]))
    assert pair(base,fg,base)==1 and all(pair(base,fg,e)==0 for e in edge)
    y0=T/(2-T);s0=2/(2-T);y,s=O,Z;translations={0:tuple(tuple(K(i==k)for k in range(6))for i in range(6))}
    for power in range(1,5):
        y,s=y*y0+ALPHA*s*s0,y*s0+s*y0
        for sign in(-1,1):translations[power*sign]=tuple(tuple(K(i==k)+(y-1)*p[i][k]+sign*s*j[i][k]for k in range(6))for i in range(6))
    oldkeys={ray_key(r)for r in walls};cov=[mv(g,r)for r in walls];norms=[dot(r,c)for r,c in zip(walls,cov)]
    target=tuple(map(K.from_json,data['forced_H5_face']['geometry']['vector']))
    profiles=[];dictionary={};integral_candidates=[]
    for q in NORMS:
        alphabet=[sorted({r for ratio in ALLOWED if(r:=square_root(q*pair(e,fg,e)*ratio/4))is not None and r.is_integral})for e in edge]
        for fixed in product(*alphabet):
            pp=tuple(-z for z in fixed);coeff=mv(ei,pp)
            r0=tuple(sum((c*e[k]for c,e in zip(coeff,edge)),Z)for k in range(6))
            qp=q-pair(r0,fg,r0)
            if any(qp.sign(e)<=0 for e in range(4)):continue
            factor=square_root(qp);profile={'norm':q.to_json(),'fixed_magnitudes':enc(fixed),
                'fixed_component':enc(r0),'orthogonal_norm':qp.to_json(),'square_factor':None if factor is None else factor.to_json()}
            index=len(profiles);profiles.append(profile)
            if factor is None:continue
            for power,H in translations.items():
                translated=mv(H,base)
                for sign in(-1,1):
                    r=tuple(a+sign*factor*b for a,b in zip(r0,translated))
                    assert pair(r,fg,r)==q
                    if r[0]<=0 or not all(z.is_integral for z in r):continue
                    if all((z/T).is_integral for z in r):continue
                    if not all((2*z/q).is_integral for z in mv(fg,r)):continue
                    r9=(r[0],*r[1:5],*(Z,)*4,r[5]);key=ray_key(r9)
                    record={'profile_index':index,'unit_power':power,'orthogonal_sign':sign,'vector':enc(r9),
                            'norm':q.to_json(),'height':(ALPHA*r[0]*r[0]/q).to_json(),
                            'old_H5_escape_pairing':pair(r,fg,target).to_json()}
                    integral_candidates.append(record)
                    if key in oldkeys or key in dictionary or not coxeter(r9,q,cov,norms):continue
                    dictionary[key]=record
    assert len(profiles)==18 and sum(r['square_factor']is not None for r in profiles)==15
    rows=sorted(dictionary.values(),key=lambda row:K.from_json(row['height']))
    inside=tuple(map(K.from_json,data['strict_interior']));fi=(*inside[:5],inside[9]);cases=[]
    for row in rows:
        if row['unit_power']!=1:continue
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm']);new=sorted(orbit(r));wc=list(cov);wq=list(norms)
        for image in new:
            assert coxeter(image,q,wc,wq);wc.append(mv(g,image));wq.append(q)
        combined=walls+new;fc=[(*w[:5],w[9])for w in combined if not any(w[5:9])]
        fr=(*r[:5],r[9]);new_index=fc.index(fr)
        ids9=[walls.index((e[0],*e[1:5],*(Z,)*4,e[5]))for e in edge]
        ids9 += [i for i,w in enumerate(walls)if not w[0]and any(w[5:9])]
        oldv=tuple(map(K.from_json,data['geometry']['vector']))
        geo=certificate(combined,g,inside,oldv)if all(pair(w,g,oldv)<=0 for w in combined)else scaled_geometry(combined,g,inside,9)
        case={'dictionary_index':rows.index(row),'norm':q.to_json(),'height':row['height'],'new_roots':menc(new),
              'all_walls':menc(combined),'wall_count':len(combined),'strict_interior':enc(inside),
              'finite_H9_vertex':finite_vertex(combined,g,inside,ids9+[combined.index(r)]),'geometry':geo,
              'forced_H5_face':{'all_walls':menc(fc),'wall_count':len(fc),'strict_interior':enc(fi),
                                'finite_vertex':finite_vertex(fc,fg,fi,edge_ids+[new_index]),
                                'geometry':scaled_geometry(fc,fg,fi,5)}}
        cases.append(case)
        print({'case':len(cases)-1,'height':q.real()and K.from_json(row['height']).real(),'face_geometry':case['forced_H5_face']['geometry']['status']},flush=True)
    out={'status':'EXACT_FINITE_COUPLING_PROFILE_ORBIT_SEARCH','source':str(SOURCE.relative_to(ROOT)),
         'edge_indices':edge_ids,'edge_gram':menc(eg),'unit_orbit_exponents':list(range(-4,5)),
         'orthogonal_unit_norm_seed':enc(base),'profiles':profiles,'integral_primitive_candidates':integral_candidates,
         'compatible_dictionary':rows,'cases':cases,'seconds':round(time.monotonic()-started,3),
         'scope':'Complete finite-angle coupling profiles for the fixed B2+A1+A1 edge; only the 15 profiles with square orthogonal norm in K are explored through the explicitly stated nine unit powers and both orthogonal signs. This is not a complete Pell solution classification. Each unit-power-one orbit is an alternative cone from the original174 walls; conflicting endpoint caps are not silently merged.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print({'profiles':len(profiles),'dictionary':len(rows),'cases':len(cases),'seconds':out['seconds']},flush=True)


if __name__=='__main__':run_profiles()
