#!/usr/bin/env python3
"""Complete one-F4-plus-scalar slice, targeted at the new H9 escape."""
import json
from itertools import groupby
import time
from fresh_h9_glued_host_20261013 import *

SOURCE=ROOT/'results/fresh_h9_glued_host_20261013_h9_h100.json'
OUT=ROOT/'results/fresh_h9_glued_host_20261013_target.json'

def run(height=300):
    started=time.monotonic();source=json.loads(SOURCE.read_text());host=Host(1);full=Host(2);g=full.g
    old=tuple(tuple(map(K.from_json,r))for r in source['all_walls'])
    target=tuple(map(K.from_json,source['geometry']['vector']));inside=tuple(map(K.from_json,source['strict_interior']))
    controls=sorted((ALPHA*x*x/q,q,x)for q in NORMS for x in alphabet((height*q/ALPHA,-q/ALPHA,-q/ALPHA,-q/ALPHA))if x>0)
    levels=[];dictionary={}
    for h,q,x in controls:
        roots,stats=host.candidates(x,q);cutters=[];accepted=[];reject=0
        for r in roots:
            for block in(0,1):
                rr=(r[0],*(r[1:5]if block==0 else(Z,)*4),*(r[1:5]if block==1 else(Z,)*4),r[5])
                if pair(rr,g,target)<=0:continue
                cutters.append((block,rr))
                if not compatible(rr,q,old,g):reject+=1;continue
                accepted.append((block,rr));dictionary[ray_key(rr)]={'height':h.to_json(),'norm':q.to_json(),'vector':enc(rr),
                    'target_pairing':pair(rr,g,target).to_json(),'block':block}
        levels.append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),'slice_roots':len(roots),
            'cutting_embeddings':len(cutters),'old_incompatible':reject,'accepted_embeddings':len(accepted),**stats})
        if cutters:print(json.dumps({'height':h.real(),'cuts':len(cutters),'accepted':len(accepted)}),flush=True)
    rows=sorted(dictionary.values(),key=lambda row:(K.from_json(row['height']),tuple(map(K.from_json,row['vector']))))
    walls=list(old);chosen=[];rejected=[]
    for i,row in enumerate(rows):
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm'])
        if compatible(r,q,walls,g):chosen.append(i);walls.append(r)
        else:rejected.append(i)
    while any(pair(inside,g,r)>=0 for r in walls):inside=(inside[0],*(x/2 for x in inside[1:]))
    out={'status':'COMPLETE_BOUNDED_FULL_F4_SCALAR_TARGET','source':str(SOURCE.relative_to(ROOT)),
        'height_bound':height,'control_count':len(controls),'norm_pool':enc(NORMS),'levels':levels,
        'compatible_dictionary':rows,'selected_indices':chosen,'rejected_indices':rejected,
        'all_walls':menc(walls),'wall_count':len(walls),'strict_interior':enc(inside),
        'geometry':geometry(walls,g,inside,9),'seconds':round(time.monotonic()-started,3),
        'scope':'Complete bounded search with one full F4 block and scalar free, in both block embeddings, retaining roots which cut the saved H9 escape and pass all old Coxeter gates. The other F4 block is zero. No full H9 height-300 census or unrestricted completion claim.'}
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':out['status'],'controls':len(controls),'dictionary':len(rows),'selected':len(chosen),
        'walls':len(walls),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)

def finish_symmetry():
    data=json.loads(OUT.read_text());g=Host(2).g
    walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    inside=tuple(map(K.from_json,data['strict_interior']));added=[]
    for row in data['compatible_dictionary']:
        r=tuple(map(K.from_json,row['vector']));s=(r[0],*r[5:9],*r[1:5],r[9]);q=pair(s,g,s)
        assert compatible(s,q,walls,g);walls.append(s);added.append(s)
    assert all(pair(inside,g,r)<0 for r in walls)
    symmetry={'status':'EXACT_BLOCK_SWAP_COMPLETION','source':str(OUT.relative_to(ROOT)),
        'added_walls':menc(added),'all_walls':menc(walls),'wall_count':len(walls),'strict_interior':enc(inside),
        'geometry':geometry(walls,g,inside,9),
        'scope':'Exact block swaps of all four accepted target cutters; every old and new angle checked. This is a symmetry completion, not an additional root search.'}
    path=ROOT/'results/fresh_h9_glued_host_20261013_symmetry.json'
    path.write_text(json.dumps(symmetry,indent=2,sort_keys=True)+'\n')
    face_walls=[(*r[:5],r[9])for r in walls if all(not x for x in r[5:9])]
    face_inside=(*inside[:5],inside[9]);face_g=Host(1).g
    face={'status':'EXACT_FORCED_F4_COMPLEMENT','source':str(path.relative_to(ROOT)),'blocks':1,
        'all_walls':menc(face_walls),'wall_count':len(face_walls),'strict_interior':enc(face_inside),
        'geometry':geometry(face_walls,face_g,face_inside,5)}
    (ROOT/'results/fresh_h9_glued_host_20261013_h5_completed.json').write_text(json.dumps(face,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'symmetry_walls':len(walls),'symmetry_geometry':symmetry['geometry']['status'],
                     'forced_face_walls':len(face_walls),'forced_face_geometry':face['geometry']['status']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--symmetry-only',action='store_true');a=p.parse_args()
    if not a.symmetry_only:run()
    finish_symmetry()
