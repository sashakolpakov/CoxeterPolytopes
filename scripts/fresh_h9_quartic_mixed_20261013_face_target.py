#!/usr/bin/env python3
"""Neighboring-weight completion of the forced H5 face, with H9 pair gates."""
from fresh_h9_quartic_mixed_20261013_target import *

BASE=ROOT/'results/fresh_h9_quartic_mixed_20261013_h1000.json'
PRIOR=ROOT/'results/fresh_h9_quartic_mixed_20261013_target.json'
FACE=ROOT/'results/fresh_h9_quartic_mixed_20261013_h5.json'


def run_face(height,output):
    started=time.monotonic();data=json.loads(BASE.read_text());prior=json.loads(PRIOR.read_text())
    face=json.loads(FACE.read_text());g=Host(2).g;fg=Host(1).g
    old=[tuple(map(K.from_json,r))for r in data['all_walls']]
    cov=[mv(g,r)for r in old];norms=[dot(r,c)for r,c in zip(old,cov)]
    keys={ray_key(r)for r in old};fv=tuple(map(K.from_json,face['geometry']['vector']))
    controls=sorted((ALPHA*x*x/q,q,x)for q in NORMS
                   for x in balanced_alphabet((height*q/ALPHA,-q/ALPHA,-q/ALPHA,-q/ALPHA))if x>0)
    dictionary={};levels=[]
    print(json.dumps({'scope':'three neighboring-weight H5 slices plus symmetries',
                      'height':height,'controls':len(controls)}),flush=True)
    for index,(h,q,x)in enumerate(controls):
        energy=q+ALPHA*x*x
        for j,k in((0,2),(1,2),(2,3)):
            aa,bb,cc=FI[j][j],FI[k][k],FI[j][k];det=aa*bb-cc*cc
            us=balanced_alphabet((energy/aa,energy*bb/det,energy*bb/det,energy*bb/det))
            square_roots={z*z:z for z in balanced_alphabet((bb*energy,)*4)}
            roots=[];accepted=0;cuts=0
            for u in us:
                square=bb*energy-det*u*u
                if square not in square_roots:continue
                v=(square_roots[square]-cc*u)/bb
                if u<=0 or v<=0 or not v.is_integral:continue
                r=(x,*(row[j]*u+row[k]*v for row in FI),Z)
                assert pair(r,fg,r)==q and all(z.is_integral for z in r)
                if all((z/T).is_integral for z in r):continue
                assert all((2*z/q).is_integral for z in mv(fg,r))
                roots.append(r)
                rr=(r[0],*r[1:5],*(Z,)*5)
                for image in orbit(rr):
                    key=ray_key(image)
                    if key in keys or key in dictionary:continue
                    if not coxeter(image,q,cov,norms):continue
                    accepted+=1
                    face_pair=pair((*image[:5],image[9]),fg,fv)
                    if not any(image[5:9]):cuts+=face_pair>0
                    dictionary[key]={'height':h.to_json(),'norm':q.to_json(),'vector':enc(image),
                                     'origin_slice':[j,k],'first_H5_escape_pairing':face_pair.to_json()}
            levels.append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
                           'covector_indices':[j,k],'u_alphabet':len(us),'square_alphabet':len(square_roots),
                           'roots':menc(roots),'accepted_images':accepted,'cutting_face_images':cuts})
            if accepted:
                print(json.dumps({'index':index,'height':h.real(),'slice':[j,k],
                                  'roots':len(roots),'accepted':accepted,'cuts_face':cuts}),flush=True)
        if index%50==0:print(json.dumps({'index':index,'elapsed':round(time.monotonic()-started,2)}),flush=True)
    new=sorted(dictionary.values(),key=lambda row:(K.from_json(row['height']),tuple(map(K.from_json,row['vector']))))
    merged={ray_key(tuple(map(K.from_json,row['vector']))):row for row in prior['compatible_dictionary']}
    merged.update(dictionary)
    rows=sorted(merged.values(),key=lambda row:(K.from_json(row['height']),tuple(map(K.from_json,row['vector']))))
    walls=list(old);wc=list(cov);wq=list(norms);chosen=[];rejected=[]
    for i,row in enumerate(rows):
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm'])
        if coxeter(r,q,wc,wq):chosen.append(i);walls.append(r);wc.append(mv(g,r));wq.append(q)
        else:rejected.append(i)
    inside=tuple(map(K.from_json,data['strict_interior']))
    while any(dot(inside,c)>=0 for c in wc):inside=(inside[0],*(z/2 for z in inside[1:]))
    fw=[(*r[:5],r[9])for r in walls if not any(r[5:9])];fi=(*inside[:5],inside[9])
    out={'status':'COMPLETE_BOUNDED_NEIGHBORING_WEIGHT_FACE_COMPLETION','source':str(BASE.relative_to(ROOT)),
         'prior_source':str(PRIOR.relative_to(ROOT)),'height_bound':height,'control_count':len(controls),
         'norm_pool':enc(NORMS),'levels':levels,'new_compatible_dictionary':new,
         'combined_dictionary':rows,'selected_indices':chosen,'rejected_indices':rejected,
         'all_walls':menc(walls),'wall_count':len(walls),'strict_interior':enc(inside),
         'geometry':geometry(walls,g,inside,9),
         'forced_H5_face':{'all_walls':menc(fw),'wall_count':len(fw),'strict_interior':enc(fi),
                           'geometry':geometry(fw,fg,fi,5)},
         'seconds':round(time.monotonic()-started,3),
         'scope':'Bounded complete slices with two distinct, strictly positive F4 fundamental covectors and scalar zero, all three pairs containing the escaping inner weight and their diagram/block symmetries. Boundary one-weight planes are not included. All original62-compatible roots retained before height-ordering their union with the previous mixed-block dictionary. This does not exclude other roots or decide full-host reflectivity.'}
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'new_dictionary':len(new),'selected':len(chosen),'rejected':len(rejected),
                      'walls':len(walls),'geometry':out['geometry']['status'],'H5_walls':len(fw),
                      'H5_geometry':out['forced_H5_face']['geometry']['status'],'seconds':out['seconds']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--height',type=int,default=100000)
    p.add_argument('--output',type=Path,default=ROOT/'results/fresh_h9_quartic_mixed_20261013_face_target.json')
    a=p.parse_args();run_face(a.height,a.output)
