#!/usr/bin/env python3
"""Complete raw357 relative-role inventory and bounded native marker tests.

Whole A2 roles and reflections are reduced before pair gates. No ordinary
pair outcome, capacity, column or global optimization is asserted here.
"""
import argparse,gc,json,time
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_raw357_extra_vertex_incidence_pilot_20261018 as G
import fresh_raw112_complete_marker_family_20261018 as L
C=G.C;A=G.A;R=G.R;K=G.K;ROOT=G.ROOT;T=L.T;H=L.H;U=L.U
TAG='fresh_raw357_relative_marker_pilot_20261018'
KINDS=(3574,3575,3576)
OLD=(56,57,76,1084,1085,1124,1125)
PINS={'scripts/'+G.TAG+'.py':'83242d583ac43fc0e0101dfa32c4a89fab8277ee78afe4bbb326bb048f75eeb6',
 'results/'+G.TAG+'.json':'468a7311c5e0295bcfbb8e582364f75ed363f539828d825a45204b5dfdb0c84b',
 'results/'+G.TAG+'.md':'90dbae64b2cc7cf6b68da672df06a2aa35bd3854e3d2a35c591b4032686fd083'}
normal=L.normal


def maps(row,g,extra):
    ids=extra['support'];om=extra['omitted'];kind=3570+om[1]
    target=R.P.forest_graph(row['target_forest'],dict(R.C.models()))
    J,ac,blocks,tens=T.projection_data(target)
    em=R.active_embeddings(C.active_graph(g,ids),target)
    assert len(em)==row['active_role_orientations']==C.expected_maps(extra['active_forest'],row['target_forest'])
    anchor=em[0];cov=[tuple(sorted((anchor[j],tuple(g[i][j]))for j in ids if g[i][j]!=A.RIGHT))for i in om]
    hs=[]
    for c in cov:
        f=L.field_for((c,),tens);z=L.projection(c,c,J,tens,f);assert not any(z[1:]);h=z[0]-2
        assert h.sign()>0;hs.append(h)
    f=L.field_for(cov,tens);z=L.projection(cov[0],cov[1],J,tens,f)
    prescribed=f.negative_contact(g[om[0]][om[1]])
    assert f.signs_difference(tuple(x-y for x,y in zip(z,prescribed)),hs[0]*hs[1])==0
    markers=[];seen=set()
    for mp in em:
        assert set(mp)==set(ids)and len(set(mp.values()))==8
        assert all(set(target[mp[i]])<=set(mp.values())for i in ids)
        assert all(target[mp[i]].get(mp[j],2)==g[i][j][0]for i,j in combinations(ids,2))
        cv=[tuple(sorted((mp[j],tuple(g[i][j]))for j in ids if g[i][j]!=A.RIGHT))for i in om]
        order=sorted(range(2),key=lambda i:cv[i])
        native=dict(contacts=[[[j,list(e)]for j,e in cv[i]]for i in order],
            squared_timelike_heights=[hs[i].out()for i in order],
            native_omitted_product=list(g[om[0]][om[1]]),occupied_pairs=[])
        key=H.frozen(native['contacts']);assert key not in seen;seen.add(key)
        markers.append(dict(index=len(markers),kind=kind,active_embedding=sorted(mp.items()),native=native,
            triangle_port_pairs=sorted(sorted((mp[i],mp[j]))for i,j in extra['forced_triangle_port_pairs'])))
    return target,markers


def components(target,marker):
    mp=dict(marker['active_embedding']);reverse={v:k for k,v in mp.items()}
    todo=set(mp.values());out=[]
    while todo:
        first=min(todo,key=reverse.get);part={first};stack=[first]
        while stack:
            i=stack.pop()
            for j in target[i]:
                assert j in reverse
                if j not in part:part.add(j);stack.append(j)
        todo-=part;out.append(tuple(sorted(part,key=reverse.get)))
    return sorted(out,key=lambda p:min(reverse[i]for i in p))


def signature(target,ka,a,kb,b):
    left=components(target,a);right=components(target,b)
    # The first marker fixes every node of its ordered active components.
    # New right components have no contacts to it and may be independently
    # permuted and reflected. Their right-source ordering fixes the payload.
    p=dict(a['active_embedding']);names={p[i]:j for j,i in enumerate(sorted(p))};pattern=[]
    for part in right:
        shared=next((j for j,c in enumerate(left)if set(c)==set(part)),None)
        if shared is None:
            assert not set(part)&set(names)
            pattern.append([-1,[]])
            for i in part:names[i]=len(names)
        else:pattern.append([shared,[left[shared].index(i)for i in part]])
    assert len(names)<=13
    ordered=sorted(names,key=names.get)
    gram=[[0 if i==j else target[i].get(j,2)for j in ordered]for i in ordered]
    cov=[];heights=[]
    for marker in(a,b):
        cs=[sorted([[names[i],e]for i,e in c])for c in marker['native']['contacts']]
        order=sorted(range(2),key=lambda i:cs[i]);cov.append([cs[i]for i in order])
        heights.append([marker['native']['squared_timelike_heights'][i]for i in order])
    payload=dict(contacted_native_gram=gram,full_contacted_covectors=cov,heights=heights,
        omitted_products=[a['native']['native_omitted_product'],b['native']['native_omitted_product']])
    return H.frozen([ka,kb,pattern]),pattern,payload


def old_inventory(vertices,paths):
    out=defaultdict(dict)
    for tag in G.PRIOR:
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):paths.add(f'{folder}/{tag}.{ext}')
        d=json.loads((ROOT/f'results/{tag}.json').read_text())
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
        for row in d['forest_records']:
            v=row['vertex']
            if v not in vertices:continue
            frame={int(i):{int(j):x for j,x in r.items()}for i,r in row['finite_frame'].items()}
            families=[dict(kind=76,markers=row['markers'])]if tag==G.PRIOR[0]else row['marker_families']
            for f in families:
                kind=f.get('kind',f.get('raw_family'));assert kind in OLD and kind not in out[v]
                out[v][kind]=dict(frame=frame,markers=f['markers'],source=tag)
        del d;gc.collect()
    return out


def run(preflight=False):
    started=time.monotonic();paths=set(PINS)
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    prior=json.loads((ROOT/f'results/{G.TAG}.json').read_text())
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert prior['counts']['source_classes']==47 and prior['counts']['marker_domains']==141
    geom=prior['geometry'];g=[[tuple(e)for e in row]for row in geom['native_gram']]
    grouped=defaultdict(list)
    for row in prior['extra_marker_inventory']:grouped[row['target_vertex']].append(row)
    assert len(grouped)==138
    old=old_inventory(set(grouped),paths)
    records=[];cases={};multiplicity=Counter();comparisons=0
    for v,rows in sorted(grouped.items()):
        fam={};target=None
        for row in rows:
            extra=next(x for x in geom['extra_vertices']if x['omitted']==row['omitted'])
            frame,markers=maps(row,g,extra)
            if target is not None:assert frame==target
            target=frame;kind=markers[0]['kind'];assert kind not in fam
            fam[kind]=markers
        oldfam=old.get(v,{})
        for kind,f in oldfam.items():assert f['frame']==target
        keys=[H.frozen(m['native']['contacts'])for ms in fam.values()for m in ms]
        assert len(keys)==len(set(keys))and not set(keys)&{H.frozen(m['native']['contacts'])for f in oldfam.values()for m in f['markers']}
        local=Counter()
        for ka,aa in sorted(fam.items()):
            anchor=aa[0]
            bs={kb:bb for kb,bb in fam.items()if kb>=ka};bs.update({kb:f['markers']for kb,f in oldfam.items()})
            for kb,bb in sorted(bs.items()):
                for b in bb:
                    if ka==kb and b['index']==anchor['index']:continue
                    key,pattern,payload=signature(target,ka,anchor,kb,b)
                    comparisons+=1;local[str((ka,kb))]+=1;multiplicity[key]+=1
                    if key in cases:assert cases[key]['payload']==normal(payload),(ka,kb,pattern)
                    else:cases[key]=dict(kinds=[ka,kb],pattern=pattern,payload=normal(payload),vertex=v,
                        forest=rows[0]['target_forest'],finite_frame=target,anchor=anchor,other=b)
        records.append(dict(vertex=v,forest=rows[0]['target_forest'],finite_frame=target,
            marker_families=[dict(kind=k,source_face=next(r['face']for r in rows if 3570+r['omitted'][1]==k),
                actual_extra_occurrences=1,markers=ms)for k,ms in sorted(fam.items())],
            old_marker_counts={str(k):len(f['markers'])for k,f in sorted(oldfam.items())},anchor_signature_counts=dict(local)))
        if len(records)%25==0:print('INVENTORY',len(records),'/',len(grouped),'cases',len(cases),'seconds',round(time.monotonic()-started,1),flush=True)
    assert sum(len(f['markers'])for r in records for f in r['marker_families'])==10624
    out=[];counts=Counter()
    for key,c in sorted(cases.items()):
        c['anchor_comparisons']=multiplicity[key];counts[str(tuple(c['kinds']))]+=1
        if not preflight:
            target={int(i):{int(j):m for j,m in r.items()}for i,r in c['finite_frame'].items()}
            J,ac,blocks,tens=T.projection_data(target)
            c['result']=L.pair_test(J,tens,c['anchor']['native'],c['other']['native'])
            print('PAIR',c['kinds'],c['pattern'],c['result']['allowed'],flush=True)
        out.append(c)
    for mod in(G,L,C,A,R,T,U,H):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    paths.add(str(Path(__file__).relative_to(ROOT)))
    return normal(dict(status='PASS_RAW357_COMPLETE_RELATIVE_PREFLIGHT'if preflight else'PASS_RAW357_COMPLETE_RELATIVE_MARKER_PILOT',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
        family=prior['complete_family_inventory'],pointed_family=prior['extra_marker_inventory'],
        geometry=geom,forest_records=records,relative_cases=out,
        counts=dict(source_classes=47,marker_domains=141,target_forests=138,role_orientations=10624,
            anchor_comparisons=comparisons,relative_cases=len(out),case_types=dict(counts),
            native_pair_tests=0 if preflight else len(out),retained_relative_cases=None if preflight else sum(c['result']['allowed']for c in out)),
        old_outcomes_changed=False,source_quad_demand=1,extra_occurrences_per_source=3,
        exported_columns=[],optimizer_invocations=0,
        scope='Complete relative marker signature inventory over all47 source classes and all141 target domains. New/new and new/old marker tests only; retained unknowns are compatible in the relaxation, not physically realized. No ordinary compatibility, capacity or one-budget column is asserted yet.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--preflight',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}{"_preflight"if args.preflight else""}.json'
    if not args.verify:assert not path.exists(),path
    out=run(args.preflight)
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
