#!/usr/bin/env python3
"""Four inherited vertices and a common escape for the small exact cutter family."""
import argparse
import json
from scripts.fresh_h9_split_endpoint_walk_exact_20261013_base import *
from scripts.fresh_h9_vertex_walk_20261013_caps import finite_vertex
from scripts.fresh_h9_split_endpoint_walk_exact_20261013 import k_gram
from scripts.fresh_h9_vertex_walk_20261013_verify import point_check

def run():
    G,I,h,old,scout,_=source();bank=read(PREFIX+'_family.json');geometry=read(PREFIX+'_family_geometry.json')
    new=[load_cap(r)for r in bank['caps']];selected=old|{r['name']:r for r in new}
    records=[]
    for vertex in scout['source_snapshot']['finite_vertices']:
        ids=vertex['active'];extra=[i for i in ids if i>=16];assert len(extra)==1
        name=wall_name(extra[0]);block=[wall_name(i)for i in ids if i<16];assert len(block)==8
        gram=k_gram([selected[n]for n in block+[name]])
        record={'block_names':block,'spherical_vertex_gram':[[enc(x)for x in r]for r in gram]}
        proof=finite_vertex(old[name],record,I,h,selected)
        point_check(proof,record,selected,I,h)
        records.append({'record':record,'certificate':proof})
    common=None;origin=None
    for ci,case in enumerate(geometry['cases']):
        for ei,entry in enumerate(case['geometry']['escapes']):
            y=tuple(map(parse,entry['certificate']['y']))
            try:common=strict_escape(y,I,h,list(selected.values()))
            except AssertionError:continue
            origin={'case':ci,'escape':ei};break
        if common:break
    assert common is not None
    return {'status':'EXACT_FOUR_INHERITED_VERTICES_AND_COMMON_ESCAPE',
            'source_sha256':{p:sha(p)for p in(PREFIX+'_family.json',PREFIX+'_family_geometry.json')},
            'all_walls':list(selected),'old_finite_vertices':records,
            'common_strict_escape':common,'escape_origin':origin,
            'scope':'The four original finite vertices survive all ten proposed caps. One strict spacelike point satisfies all82 inequalities simultaneously, so every subset of this finite ten-cap completion still has infinite volume. The full82 set is not pairwise Coxeter-compatible.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();r=run()
    path=ROOT/(PREFIX+'_checkpoint.json')
    if a.verify:assert json.loads(path.read_text())==r
    else:path.write_text(json.dumps(r,indent=2)+'\n')
    print({'status':r['status'],'inherited_vertices':len(r['old_finite_vertices']),
           'common_escape_norm':r['common_strict_escape']['norm'],'origin':r['escape_origin']})
