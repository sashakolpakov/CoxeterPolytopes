#!/usr/bin/env python3
"""Independent complete targeted-slice census and selected-cone certificates."""
import argparse
from fractions import Fraction as F
from itertools import combinations,product
import json
from pathlib import Path
import time
from fresh_cubic9_positive_20261013_independent import f,ROOT,Census,alphabet,Z,O,T,AI,A4,sq,SORT

POOL=tuple(f.scale(O,k) for k in (1,2,5,10))

def verify(source):
    data=json.loads(source.read_text());base=json.loads((ROOT/data['source']).read_text());c=Census(base)
    old=tuple(f.vecdecode(r) for r in base['all_walls']);with_scalar='rank4' in source.name
    controls=c.controls(data['height_bound'],POOL);roots={};counts={}
    for h,q,x in controls:
        energy=f.add(q,f.mul(T,sq(x)))
        weights=tuple(f.mul(f.scale(q,F(1,2)),w) for w in alphabet((f.divide(f.scale(energy,8),sq(q)),)*3)
            if f.sign(w)>=0 and f.integral(f.mul(f.scale(q,F(1,2)),w)))
        scalars=tuple(z for z in alphabet((f.scale(energy,F(1,2)),)*3)
            if f.sign(z)>=0 and f.integral(f.divide(f.scale(z,4),q))) if with_scalar else (Z,)
        ss={f.scale(sq(z),2):z for z in scalars};before=len(roots)
        residue=lambda w:tuple(int(t)%5 for t in w)
        for i,j in combinations(range(4),2):
            grouped={}
            for v in weights:grouped.setdefault(residue(f.scale(v,j+1)),[]).append(v)
            for u in weights:
                for v in grouped.get(residue(f.scale(u,-i-1)),()):
                    ww=tuple(u if k==i else v if k==j else Z for k in range(4));y=f.mv(AI,ww)
                    assert all(map(f.integral,y));e=f.dot(y,ww);z=ss.get(f.sub(energy,e))
                    if z is None:continue
                    r=(x,Z,Z,Z,Z,*y,z);assert c.norm(r)==q
                    assert all(f.integral(f.divide(f.scale(t,2),q)) for t in c.cov(r))
                    roots[f.ray_key(r)]=(q,r,c.compatible(r,old))
        counts[(h,q,x)]=(len(weights),len(roots)-before)
    recorded={(f.decode(row['height']),f.decode(row['norm']),f.decode(row['time'])):row for row in data['controls']}
    assert set(recorded)==set(counts)
    assert all(recorded[k]['alphabet_count']==v[0] and recorded[k]['new_candidates']==v[1] for k,v in counts.items())
    saved={f.ray_key(f.vecdecode(row['vector'])):row for row in data['all_slice_candidates']}
    assert set(saved)==set(roots)
    for key,(q,r,ok) in roots.items():
        row=saved[key];assert f.decode(row['norm'])==q and f.vecdecode(row['vector'])==r
        assert row['compatible_with_original_prefix']==ok
    # Selection is a recorded ordering choice; replay its exact height-ordered
    # policy rather than treating the dictionary as mutually compatible.
    ordered=sorted((row for row in data['all_slice_candidates'] if row['compatible_with_original_prefix']),key=lambda row:SORT(f.decode(row['height'])))
    # Independently verify the exact height order is monotone.
    assert all(f.sign(f.sub(f.decode(b['height']),f.decode(a['height'])))>=0 for a,b in zip(ordered,ordered[1:]))
    selected=list(old);kept=[];rejected=[]
    for row in ordered:
        r=f.vecdecode(row['vector'])
        if f.ray_key(r) in {f.ray_key(w)for w in selected}:continue
        if c.compatible(r,selected):selected.append(r);kept.append(row)
        else:rejected.append(row)
    assert kept==data['selected_new_roots'] and rejected==data['cross_incompatible_roots']
    assert tuple(selected)==tuple(f.vecdecode(r)for r in data['all_walls'])
    assert len(ordered)==data['compatible_candidates'];geometry=c.check_geometry(data)
    # A saved escape may bound the whole compatible dictionary, stronger than
    # a statement about the selected cone. Check this separately and exactly.
    escape=f.vecdecode(data['geometry']['hyperideal_witness']['vector'])
    cutters=[i for i,row in enumerate(data['all_slice_candidates']) if row['compatible_with_original_prefix']
        and f.sign(f.dot(escape,c.cov(f.vecdecode(row['vector']))))>0]
    return {'source':str(source.relative_to(ROOT)),'status':'PASS','complete_bounded_slice_census':True,
        'time_norm_controls':len(controls),'candidate_roots':len(roots),'compatible_roots':len(ordered),
        'selected_new_roots':len(kept),'geometry':geometry,'dictionary_roots_cutting_selected_escape':cutters,
        'escape_satisfies_entire_compatible_dictionary':not cutters}

def run(sources,output):
    start=time.monotonic();cases=[verify(p.resolve())for p in sources]
    structural=[]
    for raw in product(range(2),repeat=3):
        x=tuple(map(F,raw));value=f.neg(f.mul(T,sq(x)));residue=[int(t)%4 for t in value]
        assert residue not in ([1,0,0],[2,0,0])
        structural.append({'x_mod2':list(raw),'minus_T_x_squared_mod4':residue})
    ratio=f.scale(O,F(1,2));assert f.field_norm(ratio)==F(1,8)
    rank4=next((p for p in sources if 'rank4' in p.name),None)
    if rank4:
        final=json.loads(rank4.read_text());target=json.loads((ROOT/'results/fresh_cubic9_positive_20261013_target.json').read_text())
        c=Census(json.loads((ROOT/final['source']).read_text()));expected=set()
        for row in final['all_slice_candidates']:
            r=f.vecdecode(row['vector']);w=f.mv(A4,r[5:9]);active=[i for i,t in enumerate(w) if t!=Z]
            ray=(len(active)<=1 and all(f.integral(f.scale(t,F(1,5)))for t in w)) or (active==[1,2] and w[1]==w[2])
            if ray:expected.add(f.ray_key(r))
        saved={f.ray_key(f.vecdecode(row['vector']))for row in target['all_slice_candidates']}
        assert expected==saved
        cases.append({'source':'results/fresh_cubic9_positive_20261013_target.json','status':'PASS',
            'complete_bounded_slice_census':True,'method':'Exact subset of independently reenumerated saturated rank-four dictionary.',
            'candidate_roots':len(saved),'geometry':c.check_geometry(target)})
    result={'status':'PASS','arithmetic':'Independent Fraction field, unpruned coefficient alphabet, residue-indexed Cartesian two-covector enumeration.',
        'cases':cases,'symmetric_middle_binary_plane_mod4_certificate':structural,
        'new_to_previous_prime5_determinant_ratio':f.encoded(ratio),'determinant_ratio_field_norm':'1/8',
        'seconds':round(time.monotonic()-start,3)}
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','cases':len(cases),'seconds':result['seconds']}));return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sources',type=Path,nargs='+');p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.sources,a.output)
