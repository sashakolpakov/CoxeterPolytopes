#!/usr/bin/env python3
"""An integral edge translation constructs a second finite H5 endpoint."""
from fresh_h9_quartic_mixed_20261013_endpoint import *

SOURCE=ROOT/'results/fresh_h9_quartic_mixed_20261013_endpoint.json'
OUT=ROOT/'results/fresh_h9_quartic_edge_walk_20261013.json'


def matmul(a,b):
    return tuple(tuple(sum((a[i][k]*b[k][j]for k in range(len(b))),Z)
                       for j in range(len(b[0])))for i in range(len(a)))


def edge_translation(root,g):
    cov=mv(g,root);xx,aa,bb=root[0],cov[1],cov[4]
    assert cov[2]==cov[3]==root[5]==0 and pair(root,g,root)==2
    skew=((Z,-bb,aa),(bb,Z,-xx),(-aa,xx,Z))
    small_gram=((-ALPHA,Z,Z),(Z,K(2),SQRT2),(Z,SQRT2,K(2)))
    j3=tuple(tuple(ALPHA*z for z in row)for row in matmul(inverse(small_gram),skew))
    def embed(v):return(v[0],*(row[0]*v[1]+row[3]*v[2]for row in FI),Z)
    def project(v):return(v[0],v[1]-SQRT2*v[4]/2,v[4]-SQRT2*v[1]/2)
    cols=[embed(mv(j3,project(tuple(K(i==j)for i in range(6)))))for j in range(6)]
    j=tuple(tuple(cols[k][i]for k in range(6))for i in range(6))
    p=tuple(tuple(z/ALPHA for z in row)for row in matmul(j,j))
    assert matmul(p,p)==p and matmul(j,p)==j
    y0=T/(2-T);s0=2/(2-T)
    assert y0*y0-ALPHA*s0*s0==1
    y,s=y0*y0+ALPHA*s0*s0,2*y0*s0
    assert y*y-ALPHA*s*s==1
    h=tuple(tuple(K(i==k)+(y-1)*p[i][k]+s*j[i][k]for k in range(6))for i in range(6))
    hi=tuple(tuple(K(i==k)+(y-1)*p[i][k]-s*j[i][k]for k in range(6))for i in range(6))
    assert all(z.is_integral for row in h+hi for z in row)
    assert matmul(h,hi)==tuple(tuple(K(i==k)for k in range(6))for i in range(6))
    ht=tuple(zip(*h));assert matmul(matmul(ht,g),h)==g
    return h,hi,j,p,y,s


def finite_vertex(walls,g,inside,indices):
    roots=[walls[i]for i in indices]
    gram=tuple(tuple(pair(r,g,s)for s in roots)for r in roots)
    current=gram
    while current:
        pivot=current[0][0];assert pivot>0
        current=tuple(tuple(current[i][j]-current[i][0]*current[0][j]/pivot
                    for j in range(1,len(current)))for i in range(1,len(current)))
    coefficients=mv(inverse(gram),tuple(pair(inside,g,r)for r in roots))
    point=tuple(inside[j]-sum((t*r[j]for t,r in zip(coefficients,roots)),Z)for j in range(len(inside)))
    values=tuple(pair(point,g,r)for r in walls)
    assert pair(point,g,point)<0 and pair(point,g,inside)<0 and all(z<=0 for z in values)
    assert all(values[i]==0 for i in indices)
    return {'indices':indices,'gram':menc(gram),'vector':enc(point),'norm':pair(point,g,point).to_json(),
            'center_pairing':pair(point,g,inside).to_json(),'all_wall_pairings':enc(values),
            'all_active_indices':[i for i,z in enumerate(values)if z==0]}


def run(halfstep=False):
    started=time.monotonic();data=json.loads(SOURCE.read_text());g=Host(2).g;fg=Host(1).g
    old=[tuple(map(K.from_json,r))for r in data['all_walls']]
    fw=[tuple(map(K.from_json,r))for r in data['forced_H5_face']['all_walls']]
    edge=[1,2,4,12];h,hi,j,p,y,s=edge_translation(fw[12],fg)
    integral_square=h
    if halfstep:
        y,s=T/(2-T),2/(2-T)
        h=tuple(tuple(K(i==k)+(y-1)*p[i][k]+s*j[i][k]for k in range(6))for i in range(6))
        hi=tuple(tuple(K(i==k)+(y-1)*p[i][k]-s*j[i][k]for k in range(6))for i in range(6))
        assert matmul(h,h)==integral_square
    assert all(mv(h,fw[i])==fw[i]for i in edge)
    transported=mv(h,fw[6]);assert transported[0]<0
    oriented=tuple(-z for z in transported)
    reflection_multiplier=2*pair(oriented,fg,fw[12])/pair(fw[12],fg,fw[12])
    r=tuple(x-reflection_multiplier*z for x,z in zip(oriented,fw[12]))
    assert r[0]>0 and pair(r,fg,r)==2
    root=(r[0],*r[1:5],*(Z,)*4,r[5]);new=sorted(orbit(root))
    cov=[mv(g,w)for w in old];norms=[dot(w,c)for w,c in zip(old,cov)]
    for w in new:
        assert all(z.is_integral for z in w)and not all((z/T).is_integral for z in w)
        assert pair(w,g,w)==2 and coxeter(w,K(2),cov,norms)
        cov.append(mv(g,w));norms.append(K(2))
    walls=old+new;inside=tuple(map(K.from_json,data['strict_interior']))
    assert all(pair(inside,g,w)<0 for w in walls)
    face=[(*w[:5],w[9])for w in walls if not any(w[5:9])];fi=(*inside[:5],inside[9])
    face_new=face.index(r);old_face_target=tuple(map(K.from_json,data['forced_H5_face']['geometry']['vector']))
    assert pair(r,fg,old_face_target)>0
    h5vertex=finite_vertex(face,fg,fi,edge+[face_new])
    h9edge=[old.index((w[0],*w[1:5],*(Z,)*4,w[5]))for w in(fw[i]for i in edge)]
    other_zero=[i for i,w in enumerate(old)if not w[0]and any(w[5:9])]
    h9vertex=finite_vertex(walls,g,inside,h9edge+other_zero+[walls.index(root)])
    oldv=tuple(map(K.from_json,data['geometry']['vector']))
    geom=certificate(walls,g,inside,oldv)if all(pair(w,g,oldv)<=0 for w in walls)else scaled_geometry(walls,g,inside,9)
    out={'status':'EXACT_EDGE_TRANSLATION_AND_FINITE_ENDPOINT','source':str(SOURCE.relative_to(ROOT)),
         'fixed_H5_edge_indices':edge,'transported_H5_wall':6,'reflection_H5_wall':12,
         'translation':menc(h),'translation_inverse':menc(hi),'generator':menc(j),'plane_projection':menc(p),
         'translation_integral':all(z.is_integral for row in h for z in row),'integral_translation_square':menc(integral_square),
         'relative_unit':{'y':y.to_json(),'s':s.to_json(),'identity':'y^2-(T-1)*s^2=1','power_of_small_unit':1 if halfstep else 2},
         'raw_transported_root':enc(transported),'reflection_multiplier':reflection_multiplier.to_json(),
         'constructed_H5_root':enc(r),'height':(ALPHA*r[0]*r[0]/2).to_json(),
         'old_H5_escape_pairing':pair(r,fg,old_face_target).to_json(),
         'new_roots':menc(new),'all_walls':menc(walls),'wall_count':len(walls),
         'strict_interior':enc(inside),'finite_H9_vertex':h9vertex,'geometry':geom,
         'forced_H5_face':{'all_walls':menc(face),'wall_count':len(face),'strict_interior':enc(fi),
                           'finite_vertex':h5vertex,'geometry':scaled_geometry(face,fg,fi,5)},
         'seconds':round(time.monotonic()-started,3),
         'scope':'One explicit relative-unit isometry, reversal of one transported root, and a reflection in a fixed edge wall. For the half step, the isometry is rational and its square is integral; the constructed root itself is separately verified integral and reflective. The four symmetry images retain all original174 walls and close a different edge at an actual finite vertex. No uniform height census or full-host classification is claimed.'}
    output=OUT.with_name(OUT.stem+'_halfstep.json')if halfstep else OUT
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'walls':len(walls),'face_walls':len(face),'new_roots':len(new),'height':K.from_json(out['height']).real(),
                      'geometry':geom['status'],'face_geometry':out['forced_H5_face']['geometry']['status'],'seconds':out['seconds']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--halfstep',action='store_true');args=p.parse_args();run(args.halfstep)
