#!/usr/bin/env python3
"""Current provenance and conservative native-domain merge for22 old classes.

No optimizer or new physical pair test. Frozen owner domains precede the
explicit unowned-pair successor, which precedes the added historical sources.
Only forests containing an added actual corner are rebuilt.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_current_historical_weighted_corners_20261017 as H

ROOT=H.ROOT;W=H.W
TAG='fresh_historical22_corner_domain_20261017'
BASE=H.TAG
PATCH='fresh_unowned_cover_pair_geometry_20261017'
EXTRA=('fresh_mixed_raw39_actual_corners_20261017','fresh_mixed_raw23_raw25_corners_20261017')
EXPECTED={26738,26755,26758,26779,26798,26933,28631,28640,28667,28676,28678,
          28681,28705,28712,28731,28740,28742,28760,28763,28798,28810,28813}
NATIVE=H.NATIVE_FIELDS

def normalized(r,q):
    r=json.loads(json.dumps(r));family=r['raw_family']
    assert family in(23,25,39)and family==q['raw_family']
    assert r['canonical_prism_core']==q['canonical_core']
    assert r['source_context']==q['source_context']and r['provenance']==q.get('provenance')
    assert {'existing_explicit_map':'explicit','unique_saved_context':'unique_saved_context'}[r['provenance_kind']]==q['provenance_kind']
    if r['provenance_kind']=='unique_saved_context':
        assert q['source_context_count']==1
        assert r['source_contexts']==[[q['source_context'],q['canonical_core']]]
    nc=len(r['corners']);assert nc%2==0 and 2<=nc//2<=5
    assert q['required']==(nc//2)*(nc//2-1)//2
    assert len({tuple(c['omitted'])for c in r['corners']})==nc
    assert all(c['face']==r['face']and c['provenance']==r['provenance']and c['raw_family']==family for c in r['corners'])
    sr=sum(2 if x.startswith('I2(')else int(x[1:])for x in r['spectator_forest'])
    assert sr==8 and r.get('raw_core_order',7)==7 and r.get('spectator_rank',sr)==sr
    r.update(raw_core_order=7,spectator_rank=sr)
    return r

def prepare():
    cert,restricted,alpha,ralpha,order,selected,owners,named,historical,bindings=H.prepare()
    base=W.read(BASE);patch=W.read(PATCH);support=W.read(H.SUPPORT)
    assert base['face_multipliers']==patch['face_multipliers']==cert['face_multipliers']
    assert base['restricted_face_multipliers']==patch['restricted_face_multipliers']==restricted['face_multipliers']
    assert base['source_precedence']==order and not base['open_forests']
    original=set(selected);selected=dict(selected);qrows={q['face']:q for q in support['quad_rows']}
    sources={};provenance={};presentations={}
    for tag in (BASE,PATCH,*EXTRA):
        doc=W.read(tag);assert doc['status'].startswith('PASS')
        for path,h in doc['sources'].items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for tag in EXTRA:
        for raw in W.read(tag)['selected_faces']:
            fi=raw['face']
            if fi not in alpha:continue
            r=normalized(raw,qrows[fi]);presentations.setdefault(fi,[]).append(tag)
            if fi in selected:
                assert r==selected[fi],(tag,fi,'historical/current record disagreement')
            else:selected[fi]=r;provenance[fi]=tag
    assert set(selected)-original==EXPECTED and len(selected)==270
    assert Counter(selected[fi]['raw_family']for fi in EXPECTED)=={25:6,39:16}
    assert all(fi not in named for fi in EXPECTED)
    assert sum(len(selected[fi]['corners'])for fi in EXPECTED)==112
    affected=sorted({c['vertex']for fi in EXPECTED for c in selected[fi]['corners']})
    assert len(affected)==75
    original_rows={r['vertex']:r for r in base['forest_records']}
    assert set(affected)-set(original_rows)=={207,374,433}
    added=[]
    for fi in sorted(EXPECTED):
        q=qrows[fi];r=selected[fi]
        added.append(dict(face=fi,mass=str(alpha[fi]),current_quad_record=q,
            primary_historical_source=provenance[fi],matching_historical_sources=presentations[fi],
            selected_face=r,actual_corner_count=len(r['corners']),
            summed_literal_orientations=sum(c['literal_orientations']for c in r['corners'])))
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    return dict(cert=cert,restricted=restricted,alpha=alpha,ralpha=ralpha,order=order,
        selected=selected,original=original,affected=affected,added=added,original_rows=original_rows,
        patch_rows={r['vertex']:r for r in patch['forest_records']},base=base,patch=patch,sources=sources)

def assemble_forest(v,p):
    order=p['order']+list(EXTRA)
    rows=[(tag,W.forests(tag)[v])for tag in order if v in W.forests(tag)]
    forest=rows[0][1]['forest'];assert all(r['forest']==forest for _,r in rows)
    refs=[];stars=[];keys={};groups=[];orbitkeys={};owner={};origins={}
    for tag,row in rows:
        for oldg in row['corner_orbits']:
            ff={fi for fi,m in oldg['positive_face_coefficients']};assert len(ff)==1;fi=next(iter(ff))
            if fi not in p['selected']:continue
            oo={H.origin_key(c):c for c in oldg['origins']}
            expected={H.origin_key(c):c for c in p['selected'][fi]['corners']if c['vertex']==v}
            assert len(oo)==len(oldg['origins'])and set(oo)<=set(expected)
            assert all(c==expected[k]for k,c in oo.items())
            assert oldg['positive_face_coefficients']==[[fi,len(oo)]]
            assert all(c['literal_orientations']==len(oldg['literal_stars'])for c in oo.values())
            nodes=[]
            for oldi in oldg['literal_stars']:
                star=row['literal_stars'][oldi];assert star['index']==oldi
                key=H.frozen(star['contacts'])
                if key not in keys:
                    i=len(stars);keys[key]=i;stars.append(star);refs.append([tag,oldi])
                else:
                    i=keys[key];assert all(star[k]==stars[i][k]for k in NATIVE),(v,tag,oldi)
                nodes.append(i)
            assert len(nodes)==len(set(nodes))and nodes
            ok=frozenset(nodes)
            if ok in orbitkeys:
                j=orbitkeys[ok];g=groups[j]
                assert g['face']==fi and g['origins']==oo
                g['source_orbits'].append([tag,oldg['index']])
            else:
                assert not any(i in owner for i in nodes),(v,tag,'partial orbit overlap')
                j=len(groups);orbitkeys[ok]=j
                for i in nodes:owner[i]=j
                groups.append(dict(index=j,face=fi,literal_stars=sorted(nodes),origins=oo,
                    positive_face_coefficients=oldg['positive_face_coefficients'],source_orbits=[[tag,oldg['index']]]))
                for k,c in oo.items():assert k not in origins;origins[k]=c
    expected={H.origin_key(c):c for r in p['selected'].values()for c in r['corners']if c['vertex']==v}
    assert origins==expected and set(owner)==set(range(len(stars)))
    domains=[]
    for tag,row in rows:
        mp={}
        for star in row['literal_stars']:
            key=H.frozen(star['contacts'])
            if key in keys:
                i=keys[key];assert i not in mp
                assert all(star[k]==stars[i][k]for k in NATIVE)
                mp[i]=star['index']
        if mp:domains.append(dict(source=tag,kind='complete_literal_domain',mapping=mp,bad=W.relation(tag,v)))
    olddomains=[d for d in domains if d['source']in p['order']]
    newdomains=[d for d in domains if d['source']in EXTRA]
    oldmap=[];oldnodes=set()
    if v in p['original_rows']:
        old=p['original_rows'][v]['union']
        assert old['forest']==forest
        for oi,(tag,si)in enumerate(old['literal_references']):
            star=W.forests(tag)[v]['literal_stars'][si];i=keys[H.frozen(star['contacts'])]
            assert all(star[k]==stars[i][k]for k in NATIVE)
            oldmap.append([oi,i]);oldnodes.add(i)
        assert len(oldmap)==len(oldnodes)
    patchmap=[];explicit={}
    if v in p['patch_rows']:
        rr=p['patch_rows'][v];assert rr['literal_references']==p['original_rows'][v]['union']['literal_references']
        tr=dict(oldmap);patchmap=oldmap
        for oi,oj,code in rr['newly_owned_pair_domain']:
            pair=tuple(sorted((tr[oi],tr[oj])))
            assert pair not in explicit
            assert not any(all(i in d['mapping']for i in pair)for d in olddomains)
            explicit[pair]=code
    # Membership chooses the first owner. Unknown is an owned outcome too.
    owned=Counter();excluded=Counter();bad=set();unowned=0
    graph=nx.Graph();graph.add_nodes_from(range(len(stars)))
    for i,j in combinations(range(len(stars)),2):
        chosen=None;reject=False
        for d in olddomains:
            mp=d['mapping']
            if i in mp and j in mp:
                chosen=d['source'];reject=tuple(sorted((mp[i],mp[j])))in d['bad'];break
        if chosen is None and(i,j)in explicit:chosen=PATCH;reject=explicit[i,j]is not None
        if chosen is None:
            for d in newdomains:
                mp=d['mapping']
                if i in mp and j in mp:
                    chosen=d['source'];reject=tuple(sorted((mp[i],mp[j])))in d['bad'];break
        if chosen is None:unowned+=1
        else:owned[chosen]+=1
        if reject:bad.add((i,j));excluded[chosen]+=1
        else:graph.add_edge(i,j)
    assert sum(owned.values())+unowned==len(stars)*(len(stars)-1)//2
    assert sum(excluded.values())==len(bad)
    for pair,code in explicit.items():assert(pair in bad)==(code is not None)
    for g in groups:g['origins']=[g['origins'][k]for k in sorted(g['origins'],key=str)]
    source_domains=[dict(source=d['source'],kind=d['kind'],new_to_old_star_indices=sorted(d['mapping'].items()),
        first_owned_pairs=owned[d['source']],first_owned_exclusions=excluded[d['source']])for d in domains]
    u=dict(vertex=v,forest=forest,literal_references=refs,corner_orbits=groups,
        source_domains=source_domains,baseline_to_new_star_indices=oldmap,
        explicit_patch=dict(source=PATCH,new_to_patch_star_indices=[[j,i]for i,j in patchmap],
            first_owned_pairs=owned[PATCH],first_owned_exclusions=excluded[PATCH],unknown_pairs=owned[PATCH]-excluded[PATCH]),
        counts=dict(literal_stars=len(stars),corner_orbits=len(groups),actual_corners=len(origins),
            added_literal_stars=len(stars)-len(oldnodes),added_orbits=sum(g['face']in EXPECTED for g in groups),
            added_actual_corners=sum(c['face']in EXPECTED for c in origins.values()),
            total_pairs=len(stars)*(len(stars)-1)//2,excluded_pairs=len(bad),
            inherited_owned_pairs=sum(owned[d['source']]for d in olddomains),
            explicit_patch_pairs=owned[PATCH],new_historical_owned_pairs=sum(owned[d['source']]for d in newdomains),
            owned_unknown_pairs=sum(owned.values())-len(bad),unowned_allowed_pairs=unowned),
        demands=[[str(sum((a.get(fi,F())*F(m)for fi,m in g['positive_face_coefficients']),F()))for g in groups]
            for a in(p['alpha'],p['ralpha'])])
    W.relation.cache_clear()
    return u,graph,owner

def main(verify=False):
    start=time.monotonic();p=prepare();records=[];counts=Counter()
    for k,v in enumerate(p['affected'],1):
        u,graph,owner=assemble_forest(v,p);records.append(u);counts.update(u['counts'])
        if k%15==0 or k==len(p['affected']):print('FORESTS',k,'/',len(p['affected']),'added stars',counts['added_literal_stars'],flush=True)
    assert counts['added_actual_corners']==112
    assert counts['added_literal_stars']==sum(r['summed_literal_orientations']for r in p['added'])
    out=dict(status='PASS_HISTORICAL22_CURRENT_PROVENANCE_AND_NATIVE_DOMAINS',sources=p['sources'],
        original_classes=248,added_classes=22,combined_positive_classes=270,added_faces=p['added'],
        affected_forests=p['affected'],affected_existing_forests=72,new_forests=[207,374,433],
        total_represented_forests=553,unchanged_forests=478,counts=dict(counts),forest_records=records,
        owner_precedence=[dict(source=t,kind='complete_literal_domain')for t in p['order']]+
            [dict(source=PATCH,kind='explicit_tested_pairs')]+[dict(source=t,kind='complete_literal_domain')for t in EXTRA],
        face_multipliers=p['base']['face_multipliers'],restricted_face_multipliers=p['base']['restricted_face_multipliers'],
        fixed_baseline_repairs=p['patch']['repairs'],
        scope='Current provenance and complete native orbit-domain merge only. Original owner outcomes and the explicit pair successor, including unknowns, precede added historical domains. Only75 changed forests reconstructed;478 old forests inherited unchanged. No new wall-pair tests, clique/cone optimization, repaired witness, full-model feasibility or dimension theorem.')
    out=json.loads(json.dumps(out,sort_keys=True));path=ROOT/f'results/{TAG}.json'
    if verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(out['status'],dict(counts),'seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');x=a.parse_args();main(x.verify)
