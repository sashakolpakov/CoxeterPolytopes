#!/usr/bin/env python3
"""Producer-free replay of the split-seed no-removal edge caps."""
import argparse
import copy
import json
from itertools import combinations
from scripts.fresh_h9_split_endpoint_walk_exact_20261013_base import *
from scripts.fresh_h9_vertex_walk_20261013_verify import point_check
from scripts.fresh_h9_replaced_core_20261013_verify import census_check
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K

def check_cap(record,G,I,h,old,target,pool):
    row=load_cap(record);F=row['field'];E=F.base
    assert K.sign(E.d)>0 and E.sign(F.d)>0
    u,v=[tuple(map(E.parse,record[k]))for k in ('line_constant','line_direction')]
    assert len(u)==len(v)==10
    IE=tuple(tuple(E.lift(x)for x in r)for r in I)
    A=pairing(E,v,IE,v);B=E.mul(E.rational(2),pairing(E,u,IE,v))
    C=E.sub(pairing(E,u,IE,u),E.rational(2))
    assert E.sign(A)
    assert list(map(E.encode,(A,B,C)))==record['norm_polynomial']
    D=E.sub(E.mul(B,B),E.mul(E.rational(4),E.mul(A,C)))
    assert E.sign(E.sub(D,F.d))==0
    branch=record['branch'];assert branch in (-1,1)
    t=(E.div(E.neg(B),E.mul(E.rational(2),A)),E.div(E.rational(branch),E.mul(E.rational(2),A)))
    assert F.encode(t)==record['parameter']
    p=tuple(F.add(F.lift(x),F.mul(t,F.lift(y)))for x,y in zip(u,v,strict=True))
    assert all(F.sign(F.sub(x,y))==0 for x,y in zip(p,row['p'],strict=True))
    gf=tuple(tuple(lift_K(F,x)for x in r)for r in G)
    assert all(F.sign(F.sub(x,y))==0 for x,y in zip(matvec(F,gf,row['c']),p,strict=True))
    assert F.sign(F.sub(dot(F,row['c'],p),F.rational(2)))==0
    assert F.sign(dot(F,row['c'],tuple(lift_K(F,x)for x in h)))>0
    for i in range(16):
        f,pair=pair_profiles(row,old[f'seed{i}']);assert f.sign(pair)<=0
    cut=dot(F,row['c'],tuple(lift_K(F,x)for x in target))
    assert F.sign(cut)<0 and F.encode(cut)==record['target_pairing']
    block=record['block_names'];assert len(block)==len(set(block))==8
    values=[Z]*8
    for i,m in zip(record['block_support'],record['block_labels'],strict=True):values[i]=K.neg({3:O,5:P}[m])
    wants=values+[K.neg({2:Z,3:O,5:P}[record['anchor_label']])]
    for n,want in zip(block+[record['anchor_name']],wants,strict=True):
        common,er=embed(E,old[n]['field']);assert key(common)==key(E)
        c=tuple(er(x)for x in old[n]['c'])
        assert E.sign(E.sub(dot(E,c,u),E.lift(want)))==0
        assert E.sign(dot(E,c,v))==0
    gram=[[parse(x)for x in r]for r in record['spherical_vertex_gram']]
    assert inertia_K(gram)==(9,0,0)
    assert [r[-1]for r in gram[:-1]]==values
    for i,a in enumerate(block):
        for j,b in enumerate(block):
            ff,x=pair_profiles(old[a],old[b]);assert ff.sign(ff.sub(x,lift_K(ff,gram[i][j])))==0
    checks={n:relation(row,r,pool)for n,r in old.items()}
    assert checks==record['pair_checks']
    assert all(r['kind']in ('finite','dotted')for r in checks.values())and record['dropped']==[]
    point_check(record['finite_vertex'],record,old|{row['name']:row},I,h)
    return row

def run(family=False):
    prefix=PREFIX+('_family'if family else '')
    bank=read(prefix+'.json');geometry=read(prefix+'_geometry.json')
    assert all(sha(p)==digest for p,digest in bank['source_sha256'].items())
    G,I,h,old,scout,banks=source();pool=angle_pool(ANGLE_SOURCE)
    assert list(BASIS)==bank['independent_seed_basis']and list(old)==bank['old_walls']
    assert [[enc(x)for x in r]for r in G]==bank['gram']
    assert [[enc(x)for x in r]for r in I]==bank['inverse']
    assert list(map(enc,h))==bank['strict_center']
    # Independently anchor the changed sixteen-wall form, not just the archive inverse.
    a=add(scale(O,Q(1,2)),P);b=add(scale(O,Q(-5,22)),scale(P,Q(1,11)))
    eps=(1,1,1,1,-1,-1,-1,-1)
    for i in range(16):
        for j in range(16):
            expect=scale(O,2)if i==j else Z
            if i//8!=j//8:expect=sub(sub(scale(O,-2 if i%8==j%8 else 0),a),scale(b,eps[i%8]*eps[j%8]))
            assert dot(K,old[f'seed{i}']['c'],old[f'seed{j}']['p'])==expect
    assert inertia_K(G)==(9,1,0)and K.sign(pairing(K,h,I,h))<0
    assert all(K.sign(dot(K,r['c'],h))>0 for n,r in old.items()if n.startswith('seed'))
    for i,r in enumerate(G):assert matvec(K,I,r)==tuple(O if i==j else Z for j in range(10))
    target=tuple(map(parse,bank['target']['y']))
    assert strict_escape(target,I,h,list(old.values()))==bank['target']
    rows=[check_cap(r,G,I,h,old,target,pool)for r in bank['caps']]
    if family:
        pairs={f'{i}|{j}':relation(rows[i],rows[j],pool)for i,j in combinations(range(len(rows)),2)}
        assert pairs==geometry['new_pair_relations']
        # Ten nodes are small enough for an independent exhaustive subset census.
        maximal=[]
        def compatible(i,j):return pairs[f'{min(i,j)}|{max(i,j)}']['kind']in('finite','dotted')
        for mask in range(1,1<<len(rows)):
            ids=[i for i in range(len(rows))if mask>>i&1]
            if not all(compatible(i,j)for i,j in combinations(ids,2)):continue
            if any(all(compatible(i,j)for i in ids)for j in range(len(rows))if j not in ids):continue
            maximal.append(ids)
        maximal.sort(key=lambda x:(-len(x),x))
        assert maximal==geometry['maximal_accepted_cliques']
        assert maximal==[c['new_caps']for c in geometry['cases']]
    else:assert relation(*rows,pool)==geometry['new_pair_relation']
    cache={};counts={'caps':len(rows),'configurations':0,'pair_checks':0,'distinct_pair_replays':0,
                    'actual_finite_vertices':len(rows),'strict_escapes':0,'affine_frontier_checks':0}
    for case in geometry['cases']:
        ids=case['new_caps'];selected=old|{rows[i]['name']:rows[i]for i in ids}
        assert list(selected)==case['walls']and len(selected)==case['facet_count']
        checks={}
        for a,b in combinations(selected.values(),2):
            kk='|'.join(sorted((a['name'],b['name'])))
            if kk not in cache:cache[kk]=relation(a,b,pool)
            checks[kk]=cache[kk]
        assert checks==case['pair_checks']and all(r['kind']in ('finite','dotted')for r in checks.values())
        counts['pair_checks']+=len(checks)
        counts['affine_frontier_checks']+=census_check(list(selected),checks,case['affine_census'])
        assert len(case['actual_finite_vertices'])==len(ids)
        for i,proof in zip(ids,case['actual_finite_vertices'],strict=True):
            point_check(proof,bank['caps'][i],selected,I,h);counts['actual_finite_vertices']+=1
        assert case['status']==case['geometry']['status']=='EXACT_STRICT_ESCAPES'
        for escape in case['geometry']['escapes']:
            cert=escape['certificate'];assert strict_escape(tuple(map(parse,cert['y'])),I,h,list(selected.values()))==cert
            counts['strict_escapes']+=1
        counts['configurations']+=1
    counts['distinct_pair_replays']=len(cache)
    if family:
        checkpoint=read(PREFIX+'_checkpoint.json')
        assert all(sha(p)==d for p,d in checkpoint['source_sha256'].items())
        all_selected=old|{r['name']:r for r in rows}
        assert list(all_selected)==checkpoint['all_walls']
        assert len(checkpoint['old_finite_vertices'])==4
        for item in checkpoint['old_finite_vertices']:
            point_check(item['certificate'],item['record'],all_selected,I,h)
            counts['actual_finite_vertices']+=1
        cert=checkpoint['common_strict_escape']
        assert strict_escape(tuple(map(parse,cert['y'])),I,h,list(all_selected.values()))==cert
        counts['strict_escapes']+=1
    bad=copy.deepcopy(bank['caps'][0]);F=rows[0]['field']
    bad['seed_coefficients'][0]=F.encode(F.add(F.parse(bad['seed_coefficients'][0]),F.one))
    try:check_cap(bad,G,I,h,old,target,pool)
    except AssertionError:pass
    else:raise AssertionError('Forged root accepted')
    try:strict_escape(h,I,h,list(old.values()))
    except AssertionError:pass
    else:raise AssertionError('Timelike escape accepted')
    artifacts=[prefix+'.json',prefix+'_geometry.json']
    if family:artifacts.append(PREFIX+'_checkpoint.json')
    return {'status':'PASS_PRODUCER_FREE_SPLIT_EDGE_CAPS_AND_GEOMETRY',**counts,
            'source_sha256':{p:sha(p)for p in artifacts},
            'forged_cap_rejected':True,'timelike_escape_rejected':True,
            'scope':'All literal new caps, retained angles, true finite vertices, complete affine censuses and strict escapes. No exhaustive continuous-family or finite-volume existence claim.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');parser.add_argument('--family',action='store_true');args=parser.parse_args()
    prefix=PREFIX+('_family'if args.family else '')
    result=run(args.family);path=ROOT/(prefix+'_verify.json')
    if args.verify:assert read(prefix+'_verify.json')==result
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result)
