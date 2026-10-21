#!/usr/bin/env python3
"""Independent full slice census, symmetry completion and escape replay."""
import json
import time
from fresh_h9_glued_host_20261013_independent import *

PATH=ROOT/'results/fresh_h9_glued_host_20261013_target.json'
OUT=ROOT/'results/fresh_h9_glued_host_20261013_target_independent.json'

def main():
    started=time.monotonic();data=json.loads(PATH.read_text());base=json.loads((ROOT/data['source']).read_text())
    full=Host(2);slice_host=Host(1);old=tuple(vec(r)for r in base['all_walls']);target=vec(base['geometry']['vector'])
    controls=slice_host.controls(data['height_bound']);saved={(dec(r['height']),dec(r['norm']),dec(r['time'])):r for r in data['levels']}
    assert set(controls)==set(saved)and len(controls)==data['control_count']
    all_compatible=set();target_cuts={};checked=[]
    for h,q,x in controls:
        roots,imprimitive,bs,ss=slice_host.candidates(x,q);row=saved[h,q,x]
        assert len(roots)==row['slice_roots']and imprimitive==row['imprimitive']and bs==row['block_vectors']and ss==row['scalar_vectors']
        cut=accepted=rejected=0
        for r in roots:
            for block in(0,1):
                rr=(r[0],*(r[1:5]if block==0 else(Z,)*4),*(r[1:5]if block==1 else(Z,)*4),r[5])
                compatible=full.compatible(rr,old)
                if compatible:all_compatible.add(rr)
                pairing=full.pair(rr,target)
                if sign(pairing)<=0:continue
                cut+=1
                if not compatible:rejected+=1;continue
                accepted+=1;target_cuts[rr]=(h,q,pairing,block)
        assert(cut,accepted,rejected)==(row['cutting_embeddings'],row['accepted_embeddings'],row['old_incompatible'])
        checked.append({'height':enc(h),'norm':enc(q),'slice_roots':len(roots),'target_cuts':cut,'accepted':accepted})
    recorded={vec(r['vector']):r for r in data['compatible_dictionary']};assert set(recorded)==set(target_cuts)
    for r,(h,q,pairing,block)in target_cuts.items():
        row=recorded[r];assert dec(row['height'])==h and dec(row['norm'])==q and dec(row['target_pairing'])==pairing and row['block']==block
    walls=list(old);chosen=[];rejected=[]
    for i,row in enumerate(data['compatible_dictionary']):
        r=vec(row['vector'])
        if full.compatible(r,walls):chosen.append(i);walls.append(r)
        else:rejected.append(i)
    assert chosen==data['selected_indices']and rejected==data['rejected_indices']
    assert tuple(walls)==tuple(vec(r)for r in data['all_walls'])
    target_geometry=full.check_geometry(data)
    symmetry=json.loads((ROOT/'results/fresh_h9_glued_host_20261013_symmetry.json').read_text())
    swaps={(r[0],*r[5:9],*r[1:5],r[9])for r in target_cuts}
    assert swaps=={vec(r)for r in symmetry['added_walls']}
    all_walls={vec(r)for r in symmetry['all_walls']}
    assert all_walls==set(walls)|swaps
    # Stronger than target-only coverage: the completed bank includes EVERY
    # old-compatible mirror in both full F4+scalar slices through this bound.
    assert all_compatible-set(old)==all_walls-set(old)
    symmetry_geometry=full.check_geometry(symmetry)
    face=json.loads((ROOT/'results/fresh_h9_glued_host_20261013_h5_completed.json').read_text())
    expected_face={(*r[:5],r[9])for r in all_walls if all(v==Z for v in r[5:9])}
    assert expected_face=={vec(r)for r in face['all_walls']}
    face_geometry=slice_host.check_geometry(face)
    out={'status':'PASS_COMPLETE_TWO_SLICE_CENSUS_AND_COMPLETIONS','controls':checked,
        'target_dictionary_size':len(target_cuts),'new_compatible_mirrors_in_full_two_slice_census':len(all_compatible-set(old)),
        'complete_slice_bank_in_final_cone':True,'target_geometry':target_geometry,
        'symmetry_geometry':symmetry_geometry,'forced_face_geometry':face_geometry,
        'seconds':round(time.monotonic()-started,3),
        'scope':'Every primitive candidate in either one-full-F4-plus-scalar slice through height300 was independently enumerated. All original14-compatible mirrors belong to the final22-wall cone, which still has a strict exact escape. This is a finite-height slice result, not a full H9 or global nonreflectivity theorem.'}
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in('controls','scope')}))

if __name__=='__main__':main()
