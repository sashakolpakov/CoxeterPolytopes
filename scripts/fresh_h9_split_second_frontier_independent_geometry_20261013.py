#!/usr/bin/env python3
"""Second exact arithmetic replay of the 79-wall split-seed chamber and its second33 finite vertex."""
from fresh_h9_replaced_core_independent_geometry_20261013 import (
    ROOT, Tower, K, Z, O, golden, field, inverse, sha)
import argparse
import json
import time

SOURCE='data/enumerations/fresh_h9_split_second_frontier_20261013_geometry.json'
INPUT=ROOT/'results/fresh_h9_split_second_frontier_independent_geometry_20261013_input.json'
OUT=ROOT/'results/fresh_h9_split_second_frontier_independent_geometry_20261013.json'


def merge(a,b):
    """Two structural tower embeddings, without assuming radical independence."""
    if a==b:return a,lambda x:x,lambda x:x
    if b.base is None:return a,lambda x:x,a.lift
    ancestor=a
    while ancestor.base is not None:
        ancestor=ancestor.base
        if ancestor==b:
            def emb(x):
                chain=[];p=a
                while p!=b:chain.append(p);p=p.base
                for p in reversed(chain):x=(x,p.base.lift(Z))
                return x
            return a,lambda x:x,emb
    middle,ea,eb=merge(a,b.base)
    target=Tower(middle,eb(b.radicand))
    return target,lambda x:(ea(x),middle.lift(Z)),lambda x:(eb(x[0]),eb(x[1]))


def extract():
    from scripts.fresh_h9_split_second_frontier_20261013_base import selected
    from scripts.fresh_h9_quartic_completion_20261013_arithmetic import K as SourceK
    def encode_field(f):
        return {'kind':'golden'}if f is SourceK else{
            'kind':'quadratic','base':encode_field(f.base),'radicand':f.base.encode(f.d)}
    gram,inv,h,ps,case=selected()
    roots=[]
    for name,row in ps.items():
        f=row['field'];roots.append({'name':name,'field':encode_field(f),
            'coefficients':list(map(f.encode,row['c'])),'seed_pairings':list(map(f.encode,row['p'])),
            'norm':SourceK.encode(row['q'])})
    sources={SOURCE,case['source'],
        'data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013_family.json',
        'data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013_family_geometry.json',
        'data/enumerations/fresh_h9_split_endpoint_walk_20261013_scout.json',
        'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'}
    scout=json.loads((ROOT/'data/enumerations/fresh_h9_split_endpoint_walk_20261013_scout.json').read_text())
    sources.update('data/enumerations/'+f for f,i in scout['source_snapshot']['refs'])
    cap=json.loads((ROOT/case['source']).read_text())['caps'][0]
    proof=case['actual_finite_vertices'][0]
    assert proof['point_field_from_cap']=='second33'
    data={'source':SOURCE,'source_sha256':{s:sha(s)for s in sorted(sources)},
          'gram':[[SourceK.encode(x)for x in row]for row in gram],
          'inverse':[[SourceK.encode(x)for x in row]for row in inv],
          'center_covector':list(map(SourceK.encode,h)),'roots':roots,
          'first_escape':case['geometry']['escapes'][0]['certificate'],'joint_finite_vertex':proof,
          'vertex_gram_order':cap['block_names']+['second33']}
    INPUT.write_text(json.dumps(data,indent=2)+'\n')


def replay():
    start=time.monotonic();data=json.loads(INPUT.read_text())
    assert all(sha(p)==s for p,s in data['source_sha256'].items())
    case=json.loads((ROOT/SOURCE).read_text())
    assert case['facet_count']==79 and case['walls']==[r['name']for r in data['roots']]
    assert data['first_escape']==case['geometry']['escapes'][0]['certificate']
    assert data['joint_finite_vertex']in case['actual_finite_vertices']
    gram=tuple(tuple(map(golden,row))for row in data['gram']);inv,signature=inverse(gram)
    assert signature==(9,1,0)and inv==tuple(tuple(map(golden,row))for row in data['inverse'])
    h=tuple(map(golden,data['center_covector']));y=tuple(map(golden,data['first_escape']['y']))
    center_norm=K.dot(h,K.mv(inv,h));escape_norm=K.dot(y,K.mv(inv,y));center_pair=K.dot(y,K.mv(inv,h))
    assert center_norm.sign()<0 and escape_norm.sign()>0 and center_pair.sign()<0
    assert escape_norm==golden(data['first_escape']['norm'])and center_pair==golden(data['first_escape']['center_pairing'])
    profiles={}
    for row in data['roots']:
        f=field(row['field']);c=tuple(map(f.parse,row['coefficients']));p=tuple(map(f.parse,row['seed_pairings']))
        actual=f.mv(tuple(tuple(f.lift(z)for z in r)for r in gram),c)
        assert all(f.sign(f.sub(x,z))==0 for x,z in zip(actual,p))
        assert f.sign(f.sub(f.dot(c,p),f.lift(golden(row['norm']))))==0
        assert golden(row['norm'])==2*O
        assert f.sign(f.dot(c,tuple(map(f.lift,h))))>0
        assert f.sign(f.dot(c,tuple(map(f.lift,y))))>0
        profiles[row['name']]=(f,c,p)
    assert data['first_escape']['all_inequality_signs']==[1]*79
    proof=data['joint_finite_vertex'];pf=profiles[proof['point_field_from_cap']][0];point=tuple(map(pf.parse,proof['y']))
    pnorm=pf.dot(point,pf.mv(tuple(tuple(pf.lift(z)for z in r)for r in inv),point))
    ppair=pf.dot(point,pf.mv(tuple(tuple(pf.lift(z)for z in r)for r in inv),tuple(map(pf.lift,h))))
    assert pf.sign(pnorm)<0 and pf.sign(ppair)<0
    assert pf.sign(pf.sub(pnorm,pf.parse(proof['norm'])))==0
    assert pf.sign(pf.sub(ppair,pf.parse(proof['center_pairing'])))==0
    signs={}
    for name,(f,c,p)in profiles.items():
        common,ep,er=merge(pf,f)
        signs[name]=common.sign(common.dot(tuple(map(ep,point)),tuple(map(er,c))))
    assert signs==proof['all_wall_signs']and all(s>=0 for s in signs.values())
    zero=sorted(name for name,s in signs.items()if s==0)
    assert zero==proof['zero_walls']and len(zero)==9
    vg=tuple(tuple(map(golden,row))for row in proof['gram'])
    assert inverse(vg)[1]==(9,0,0)
    order=data['vertex_gram_order'];assert sorted(order)==zero
    for i,name in enumerate(order):
        f,c,p=profiles[name]
        for j,name2 in enumerate(order):
            g,c2,p2=profiles[name2];common,ea,eb=merge(f,g)
            actual=common.dot(tuple(map(ea,p)),tuple(map(eb,c2)))
            assert common.sign(common.sub(actual,common.lift(vg[i][j])))==0
    assert 'second33' in order
    # Projection to each facet of a strict acute chamber preserves the other
    # strict inequalities and remains timelike. Acute hypotheses have their
    # own producer-free all-pair replay; this narrow audit does not repeat it.
    assert len(case['pair_checks'])==79*78//2
    assert all(r['kind']in('finite','dotted')for r in case['pair_checks'].values())
    assert not case['affine_census']['affine_circuits']
    assert Tower(K,O).sign((O,-O))==0 and Tower(K,Z).sign((-O,golden(['100','0'])))<0
    return {'status':'PASS_INDEPENDENT_SPLIT_79_WALL_GEOMETRY_AND_FINITE_VERTEX','source':SOURCE,
            'wall_count':79,'norm_gram_center_escape_checks':79,'gram_signature':list(signature),
            'center_norm':center_norm.to_json(),'escape_norm':escape_norm.to_json(),
            'center_escape_pairing':center_pair.to_json(),'joint_vertex_from':'second33',
            'joint_vertex_all_inequalities':79,'joint_vertex_rank':9,'joint_vertex_norm':proof['norm'],
            
            'seconds':round(time.monotonic()-start,3),
            'scope':'Separate recursive exact arithmetic for all 79 normal/Gram identities, common strict center, first strict spacelike escape, and one actual joint finite vertex with every inequality and its physical rank-nine positive definite Gram. Loading reused once for frozen input extraction; replay imports no construction producer. Cross-pair angle and affine completeness claims rely on the construction agent\'s separate producer-free audit. The displayed cone has infinite volume; arbitrary additions are not excluded.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--extract',action='store_true');p.add_argument('--verify',action='store_true');a=p.parse_args()
    if a.extract:extract()
    result=replay()
    if a.verify:
        saved=json.loads(OUT.read_text());saved.pop('seconds');result.pop('seconds');assert result==saved
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
