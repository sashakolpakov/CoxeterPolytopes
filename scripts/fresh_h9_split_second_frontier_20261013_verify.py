#!/usr/bin/env python3
"""Producer-free exact audit of five coherent second-frontier additions."""
import argparse
import copy
import json
from itertools import combinations
from scripts.fresh_h9_split_second_frontier_20261013_base import *
from scripts.fresh_h9_vertex_walk_20261013_verify import point_check
from scripts.fresh_h9_replaced_core_20261013_verify import census_check
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K
from scripts.fresh_h9_split_second_frontier_20261013_frontier import build as frontier_check,run as frontier_run

def cap_check(record,G,I,h,old,pool):
    row=profile(record);F=row['field'];E=F.base
    assert key(E)==key(old['walk2720']['field'])==key(field_parse(record['base_field']))
    assert E.sign(F.d)>0
    u,v=[tuple(map(E.parse,record[k]))for k in('line_constant','line_direction')]
    assert len(u)==len(v)==len(row['p'])==len(row['c'])==10
    IE=tuple(tuple(lift_K(E,x)for x in r)for r in I)
    A=pairing(E,v,IE,v);B=E.mul(E.rational(2),pairing(E,u,IE,v));C=E.sub(pairing(E,u,IE,u),E.rational(2))
    assert E.sign(A)and list(map(E.encode,(A,B,C)))==record['norm_polynomial']
    D=E.sub(E.mul(B,B),E.mul(E.rational(4),E.mul(A,C)));assert E.sign(E.sub(D,F.d))==0
    branch=record['branch'];assert branch in(-1,1)
    t=(E.div(E.neg(B),E.mul(E.rational(2),A)),E.div(E.rational(branch),E.mul(E.rational(2),A)))
    assert F.encode(t)==record['parameter']
    p=tuple(F.add(F.lift(x),F.mul(t,F.lift(y)))for x,y in zip(u,v,strict=True))
    assert all(F.sign(F.sub(x,y))==0 for x,y in zip(p,row['p'],strict=True))
    assert all(F.sign(F.sub(x,y))==0 for x,y in zip(matvec(F,tuple(tuple(lift_K(F,x)for x in r)for r in G),row['c']),p,strict=True))
    assert F.sign(F.sub(dot(F,row['c'],p),F.rational(2)))==0
    assert F.sign(dot(F,row['c'],tuple(lift_K(F,x)for x in h)))>0
    block=record['block_names'];assert len(block)==len(set(block))==8
    values=[Z]*8
    for i,m in zip(record['block_support'],record['block_labels'],strict=True):values[i]=K.neg({3:O,5:P}[m])
    wants=values+[K.neg({2:Z,3:O,5:P}[record['anchor_label']])]
    for n,want in zip(block+[record['anchor_name']],wants,strict=True):
        common,er=embed(E,old[n]['field']);assert key(common)==key(E)
        c=tuple(er(x)for x in old[n]['c'])
        assert E.sign(E.sub(dot(E,c,u),lift_K(E,want)))==0 and E.sign(dot(E,c,v))==0
    checks={n:relation(row,r,pool)for n,r in old.items()}
    assert checks==record['pair_checks']and all(r['kind']in('finite','dotted')for r in checks.values())
    gram=[[parse(x)for x in r]for r in record['spherical_vertex_gram']]
    assert inertia_K(gram)==(9,0,0)and [r[-1]for r in gram[:-1]]==values
    for i,a in enumerate(block):
        for j,b in enumerate(block):
            ff,x=pair_profiles(old[a],old[b]);assert ff.sign(ff.sub(x,lift_K(ff,gram[i][j])))==0
    point_check(record['finite_vertex'],record,old|{row['name']:row},I,h)
    return row

def run():
    bank=read(PREFIX+'_caps.json');geometry=read(PREFIX+'_geometry.json');front=read(PREFIX+'_frontier.json')
    assert all(sha(p)==d for p,d in bank['source_sha256'].items())
    assert sha(PREFIX+'_caps.json')==geometry['source_sha256']
    G,I,h,old,_,oldcase,family=current();pool=angle_pool(ANGLE_SOURCE)
    assert all(sha(p)==d for p,d in family['source_sha256'].items())
    assert list(old)==bank['old_walls']and inertia_K(G)==(9,1,0)
    assert K.sign(pairing(K,h,I,h))<0
    for i,r in enumerate(G):assert matvec(K,I,r)==tuple(O if i==j else Z for j in range(10))
    rows=[cap_check(record,G,I,h,old,pool)for record in bank['caps']]
    selected=old|{r['name']:r for r in rows}
    assert list(selected)==geometry['walls']and len(selected)==geometry['facet_count']==79
    for r in selected.values():
        f=r['field'];gf=tuple(tuple(lift_K(f,x)for x in row)for row in G)
        assert all(f.sign(f.sub(x,y))==0 for x,y in zip(matvec(f,gf,r['c']),r['p'],strict=True))
        assert f.sign(f.sub(dot(f,r['c'],r['p']),f.rational(2)))==0
        assert f.sign(dot(f,r['c'],tuple(lift_K(f,x)for x in h)))>0
    pairs={}
    for a,b in combinations(selected.values(),2):
        pairs['|'.join(sorted((a['name'],b['name'])))]=relation(a,b,pool)
    assert pairs==geometry['pair_checks']and all(r['kind']in('finite','dotted')for r in pairs.values())
    census=census_check(list(selected),pairs,geometry['affine_census']);assert not geometry['affine_census']['affine_circuits']
    assert len(geometry['actual_finite_vertices'])==len(rows)
    for rec,proof in zip(bank['caps'],geometry['actual_finite_vertices'],strict=True):point_check(proof,rec,selected,I,h)
    assert geometry['geometry']['status']=='EXACT_STRICT_ESCAPES'
    for escape in geometry['geometry']['escapes']:
        cert=escape['certificate'];assert strict_escape(tuple(map(parse,cert['y'])),I,h,list(selected.values()))==cert
    starting=read(PREFIX+'_source.json')['strict_target']
    surviving=strict_escape(tuple(map(parse,starting['y'])),I,h,list(selected.values()))
    assert surviving['norm']==starting['norm']and surviving['center_pairing']==starting['center_pairing']
    # This fresh exact principal-Gram census certifies all16 reached vertices,
    # all19 finite edges, and all106 open ends; it does not trust floating times.
    assert frontier_run()==front
    bad=copy.deepcopy(bank['caps'][0]);F=rows[0]['field']
    bad['coefficients'][0]=F.encode(F.add(F.parse(bad['coefficients'][0]),F.one))
    try:cap_check(bad,G,I,h,old,pool)
    except AssertionError:pass
    else:raise AssertionError('Forged cap accepted')
    try:strict_escape(h,I,h,list(selected.values()))
    except AssertionError:pass
    else:raise AssertionError('Timelike escape accepted')
    forged=copy.deepcopy(geometry)
    e=next(r for r in forged['frontier']['edges']if r['classification']=='OPEN_TO_LIGHT_CONE_WITH_CURRENT_WALLS')
    e['classification']='CLOSED_AT_EXISTING_FINITE_VERTEX'
    try:frontier_check(forged)
    except AssertionError:pass
    else:raise AssertionError('Forged endpoint closure accepted')
    return {'status':'PASS_EXACT_SECOND_FRONTIER_CAPS_GEOMETRY_AND_EDGE_CENSUS',
            'source_sha256':{p:sha(p)for p in(PREFIX+'_caps.json',PREFIX+'_geometry.json',PREFIX+'_frontier.json')},
            'new_caps':len(rows),'facets':len(selected),'pair_checks':len(pairs),
            'affine_frontier_checks':census,'explicit_finite_vertex_points':2*len(rows),
            'exact_finite_vertices_in_reached_components':front['finite_vertices'],
            'exact_finite_finite_edges':front['finite_finite_edges'],'exact_finite_open_edges':front['finite_open_edges'],
            'exact_closed_old_open_edges':len(front['old_open_edges_now_closed']),
            'compact_two_faces':len(front['compact_two_faces']),
            'strict_spacelike_escapes':len(geometry['geometry']['escapes']),
            'starting74wall_strict_escape_also_survives':True,
            'forged_cap_rejected':True,'timelike_escape_rejected':True,'forged_edge_closure_rejected':True,
            'scope':'Literal79-wall chamber with exact caps, Coxeter angles, genuine facets,16 actual finite vertices, five bounded2-faces and exact local edge frontier. The strict escapes prove infinite volume. No first-root, residual-coset or unrestricted edge-search completeness claim.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();r=run();path=ROOT/(PREFIX+'_verify.json')
    if a.verify:assert json.loads(path.read_text())==r
    else:path.write_text(json.dumps(r,indent=2)+'\n')
    print(r)
