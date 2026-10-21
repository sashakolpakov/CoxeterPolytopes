#!/usr/bin/env python3
"""One new B3+A1 edge, using a generic exact binary complement."""
from fresh_h9_quartic_edge_walk_20261013_profiles import *

SOURCE=ROOT/'results/fresh_h9_quartic_edge_walk_20261013_profiles.json'
OUT=ROOT/'results/fresh_h9_quartic_edge_walk_20261013_next.json'


def nullspace(rows):
    a=[list(row)for row in rows];n=len(a[0]);pivots=[];level=0
    for j in range(n):
        k=next((i for i in range(level,len(a))if a[i][j]),None)
        if k is None:continue
        a[level],a[k]=a[k],a[level];p=a[level][j];a[level]=[x/p for x in a[level]]
        for i in range(len(a)):
            if i!=level:
                p=a[i][j];a[i]=[x-p*y for x,y in zip(a[i],a[level])]
        pivots.append(j);level+=1
        if level==len(a):break
    out=[]
    for j in range(n):
        if j in pivots:continue
        v=[Z]*n;v[j]=O
        for i,p in enumerate(pivots):v[p]=-a[i][j]
        out.append(tuple(v))
    assert all(dot(row,v)==0 for row in rows for v in out)
    return tuple(out)


def complement(edge,g):
    basis=nullspace(tuple(mv(g,e)for e in edge));assert len(basis)==2
    A=tuple(tuple(pair(r,g,s)for s in basis)for r in basis);aa,bb,cc=A[0][0],A[0][1],A[1][1]
    delta=bb*bb-aa*cc;k=square_root(delta/ALPHA)
    assert k is not None and k>0
    j2=((-bb/k,-cc/k),(aa/k,bb/k));B=tuple(zip(*basis))
    projection=matmul(matmul(matmul(B,inverse(A)),basis),g)
    generator=matmul(matmul(matmul(matmul(B,j2),inverse(A)),basis),g)
    assert matmul(projection,projection)==projection
    assert matmul(generator,generator)==tuple(tuple(ALPHA*x for x in row)for row in projection)
    return basis,A,k,projection,generator


def run_next(final_edge=False):
    source_path=OUT if final_edge else SOURCE
    output_path=OUT.with_name(OUT.stem.replace('_next','_f4')+'.json')if final_edge else OUT
    started=time.monotonic();source=json.loads(source_path.read_text());data=source['cases'][0]
    g=Host(2).g;fg=Host(1).g;walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    face=[tuple(map(K.from_json,r))for r in data['forced_H5_face']['all_walls']]
    edge_ids=[1,2,70,72]if final_edge else[1,2,4,70];seed_index=4 if final_edge else 12
    edge=[face[i]for i in edge_ids];eg=tuple(tuple(pair(r,fg,s)for s in edge)for r in edge);ei=inverse(eg)
    coeff=mv(ei,tuple(pair(face[seed_index],fg,e)for e in edge))
    base=tuple(face[seed_index][j]-sum((c*e[j]for c,e in zip(coeff,edge)),Z)for j in range(6))
    assert pair(base,fg,base)==1 and all(pair(base,fg,e)==0 for e in edge)
    basis,A,k,P,J=complement(edge,fg);y0=T/(2-T);s0=2/(2-T);y,s=O,Z
    translations={0:tuple(tuple(K(i==j)for j in range(6))for i in range(6))}
    for power in range(1,3):
        y,s=y*y0+ALPHA*s*s0,y*s0+s*y0
        for direction in(-1,1):
            H=tuple(tuple(K(i==j)+(y-1)*P[i][j]+direction*s*J[i][j]for j in range(6))for i in range(6))
            assert matmul(matmul(tuple(zip(*H)),fg),H)==fg
            translations[power*direction]=H
    cv=[mv(g,w)for w in walls];qn=[dot(w,c)for w,c in zip(walls,cv)];keys={ray_key(w)for w in walls}
    target=tuple(map(K.from_json,data['forced_H5_face']['geometry']['vector']))
    profiles=[];candidates=[];dictionary={}
    for q in NORMS:
        alph=[sorted({r for ratio in ALLOWED if(r:=square_root(q*pair(e,fg,e)*ratio/4))is not None and r.is_integral})for e in edge]
        for fixed in product(*alph):
            cf=mv(ei,tuple(-z for z in fixed));r0=tuple(sum((c*e[j]for c,e in zip(cf,edge)),Z)for j in range(6));qp=q-pair(r0,fg,r0)
            if any(qp.sign(e)<=0 for e in range(4)):continue
            factor=square_root(qp);idx=len(profiles)
            profiles.append({'norm':q.to_json(),'fixed_magnitudes':enc(fixed),'fixed_component':enc(r0),
                             'orthogonal_norm':qp.to_json(),'square_factor':None if factor is None else factor.to_json()})
            if factor is None:continue
            for power,H in translations.items():
                v=mv(H,base)
                for sign in(-1,1):
                    r=tuple(a+sign*factor*b for a,b in zip(r0,v));assert pair(r,fg,r)==q
                    if r[0]<=0 or not all(z.is_integral for z in r):continue
                    if all((z/T).is_integral for z in r)or not all((2*z/q).is_integral for z in mv(fg,r)):continue
                    r9=(r[0],*r[1:5],*(Z,)*4,r[5]);key=ray_key(r9)
                    record={'profile_index':idx,'unit_power':power,'orthogonal_sign':sign,'vector':enc(r9),
                            'norm':q.to_json(),'height':(ALPHA*r[0]*r[0]/q).to_json(),'old_H5_escape_pairing':pair(r,fg,target).to_json()}
                    candidates.append(record)
                    if key not in keys and key not in dictionary and coxeter(r9,q,cv,qn):dictionary[key]=record
    rows=sorted(dictionary.values(),key=lambda r:K.from_json(r['height']));inside=tuple(map(K.from_json,data['strict_interior']));fi=(*inside[:5],inside[9]);cases=[]
    for row in rows:
        if K.from_json(row['old_H5_escape_pairing'])<=0:continue
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm']);new=sorted(orbit(r));cc=list(cv);qq=list(qn)
        if not all(coxeter(a,q,[mv(g,b)],[q])for i,a in enumerate(new)for b in new[i+1:]):continue
        for rnew in new:assert coxeter(rnew,q,cc,qq);cc.append(mv(g,rnew));qq.append(q)
        combined=walls+new;fc=[(*w[:5],w[9])for w in combined if not any(w[5:9])];fr=(*r[:5],r[9])
        oldv=tuple(map(K.from_json,data['geometry']['vector']))
        geo=certificate(combined,g,inside,oldv)if all(pair(w,g,oldv)<=0 for w in combined)else scaled_geometry(combined,g,inside,9)
        ids9=[walls.index((e[0],*e[1:5],*(Z,)*4,e[5]))for e in edge]
        ids9 += [i for i,w in enumerate(walls)if not w[0]and any(w[5:9])]
        case={'dictionary_index':rows.index(row),'new_roots':menc(new),'all_walls':menc(combined),'wall_count':len(combined),
              'strict_interior':enc(inside),'finite_H9_vertex':finite_vertex(combined,g,inside,ids9+[combined.index(r)]),'geometry':geo,
              'forced_H5_face':{'all_walls':menc(fc),'wall_count':len(fc),'strict_interior':enc(fi),
                                'finite_vertex':finite_vertex(fc,fg,fi,edge_ids+[fc.index(fr)]),'geometry':scaled_geometry(fc,fg,fi,5)}}
        cases.append(case)
        print({'height':K.from_json(row['height']).real(),'walls':len(combined),'face_geometry':case['forced_H5_face']['geometry']['status']},flush=True)
        if len(cases)>=3:break
    out={'status':'EXACT_F4_EDGE_PROFILE_ORBITS'if final_edge else'EXACT_B3_A1_EDGE_PROFILE_ORBITS','source':str(source_path.relative_to(ROOT)),'source_case':0,
         'seed_wall_index':seed_index,
         'edge_indices':edge_ids,'edge_gram':menc(eg),'complement_basis':menc(basis),'complement_gram':menc(A),
         'discriminant_square_factor':k.to_json(),'plane_projection':menc(P),'generator':menc(J),'unit_norm_seed':enc(base),
         'unit_exponents':[-2,-1,0,1,2],'profiles':profiles,'integral_candidates':candidates,'compatible_dictionary':rows,
         'cases':cases,'seconds':round(time.monotonic()-started,3),
         'scope':'All finite coupling profiles of the stated new edge, but only square residual-norm classes and five powers of one explicit relative unit. This is a bounded orbit construction, not a complete binary norm classification.'}
    output_path.write_text(json.dumps(out,indent=2)+'\n');print({'profiles':len(profiles),'dictionary':len(rows),'cases':len(cases),'seconds':out['seconds']},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--f4',action='store_true');args=p.parse_args();run_next(args.f4)
