#!/usr/bin/env python3
"""Promote the completed H650 portion of a stopped H1000 attempt."""
import argparse,hashlib
from fresh_h9_new_host_next_20261013_a9 import *
SOURCE=ROOT/'results/fresh_h9_new_host_next_20261013_a9_h1000.json'
OUT=ROOT/'results/fresh_h9_new_host_next_20261013_a9_h650.json'

def run(verify=False):
    d=json.loads(SOURCE.read_text());assert not d['completed']and d['height_bound']==1000
    if d['status']=='RUNNING':
        assert not verify
        d['status']='STOPPED_AT_COMPLETED_HEIGHT_GROUP';d['attempt_scope']='The requested H1000 run was stopped after completing the saved height groups. It is not a complete H1000 census. A separately checked H650 prefix is authoritative.'
        SOURCE.write_text(json.dumps(d,indent=2)+'\n')
    assert d['status']=='STOPPED_AT_COMPLETED_HEIGHT_GROUP'
    controls=[];reject=0
    for q in POOL:
        for x in balanced_alphabet((650*q/ALPHA,*(-q/ALPHA for _ in range(3)))):
            if x<=0:continue
            if not(2*ALPHA*x/q).is_integral:reject+=1;continue
            controls.append((ALPHA*x*x/q,q,x))
    rows=[r for r in d['levels']if K.from_json(r['height'])<=650]
    assert set(controls)=={tuple(K.from_json(r[k])for k in('height','norm','time'))for r in rows}
    walls=[tuple(map(K.from_json,r))for r in d['initial_walls']]+[tuple(map(K.from_json,r))for row in rows for r in row['accepted']]
    direction=mv(AI,(O,)*9);eps=O
    while True:
        inside=(O,*(eps*y for y in direction))
        if pair(inside,G,inside)<0 and all(pair(inside,G,r)<0 for r in walls):break
        eps/=2
    assert all(compatible(r,pair(r,G,r),walls[:i],G)for i,r in enumerate(walls))
    out={**d,'status':'COMPLETE_BOUNDED_A9_HOST','completed':True,'height_bound':650,'control_count':len(controls),
      'time_integrality_rejections':reject,'levels':rows,'all_walls':menc(walls),'wall_count':len(walls),
      'strict_interior':enc(inside),'geometry':geometry(walls,G,inside,9),
      'checkpoint_source':str(SOURCE.relative_to(ROOT)),'checkpoint_source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
      'scope':'Complete dominant primitive-root census through H650, extracted from completed height groups of the stopped H1000 attempt. All time controls through650 are regenerated exactly and matched. No H1000 completeness claim.'}
    del out['attempt_scope']
    if verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['control_count'],out['wall_count'],out['geometry']['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
