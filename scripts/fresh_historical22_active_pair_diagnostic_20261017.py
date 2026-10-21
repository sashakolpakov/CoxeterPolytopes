#!/usr/bin/env python3
"""Read-only priorities from the75 changed exact graph-cone certificates.

No new pair gate, clique census, local optimization or witness repair.
The complete unowned/owned-unknown orbit-pair domains are distinguished.
"""
import argparse,json,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_historical22_corner_domain_20261017 as D

ROOT=D.ROOT;W=D.W
TAG='fresh_historical22_active_pair_diagnostic_20261017'
CONES='fresh_historical22_corner_cones_20261017'

def normalize(x):return json.loads(json.dumps(x,sort_keys=True))

def run(verify=False):
    start=time.monotonic();p=D.prepare();cones=W.read(CONES);domain=W.read(D.TAG)
    assert cones['status']=='PASS_COMPLETE_HISTORICAL22_CHANGED_CORNER_CONES'and not cones['open_forests']
    assert cones['face_multipliers']==domain['face_multipliers']
    assert cones['restricted_face_multipliers']==domain['restricted_face_multipliers']
    sources={}
    for tag in(CONES,D.TAG,D.PATCH):
        doc=W.read(tag)
        for path,h in doc['sources'].items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for path in(Path(__file__),Path(D.__file__)):
        rel=str(path.relative_to(ROOT));sources[rel]=W.digest(rel)
    by={r['vertex']:r for r in cones['forest_records']}
    mus=[{v:F(x)for v,x in r['vertex_multipliers']}for r in cones['repairs']]
    rows=[];candidates=[];totals=Counter()
    for v in p['affected']:
        u,g,owner=D.assemble_forest(v,p);assert normalize(u)==next(r for r in domain['forest_records']if r['vertex']==v)
        olddom=[d for d in u['source_domains']if d['source']in p['order']]
        newdom=[d for d in u['source_domains']if d['source']in D.EXTRA]
        oldsets=[(d['source'],set(i for i,j in d['new_to_old_star_indices']))for d in olddom]
        newsets=[(d['source'],set(i for i,j in d['new_to_old_star_indices']))for d in newdom]
        explicit={}
        if v in p['patch_rows']:
            tr={j:i for i,j in u['explicit_patch']['new_to_patch_star_indices']}
            explicit={tuple(sorted((tr[i],tr[j]))):code for i,j,code in p['patch_rows'][v]['newly_owned_pair_domain']}
        def first_owner(i,j):
            for tag,ss in oldsets:
                if i in ss and j in ss:return tag
            if(i,j)in explicit:return D.PATCH
            for tag,ss in newsets:
                if i in ss and j in ss:return tag
            return None
        domains=defaultdict(list);pairkey={};owner_counts=Counter()
        for i,j in sorted(g.edges()):
            i,j=sorted((i,j));tag=first_owner(i,j);kind='unowned'if tag is None else'owned_unknown'
            key=(kind,*sorted((owner[i],owner[j])))
            domains[key].append([i,j]);pairkey[i,j]=key;owner_counts[kind]+=1
            if tag==D.PATCH:assert explicit[i,j]is None
        assert owner_counts['unowned']==u['counts']['unowned_allowed_pairs']
        assert owner_counts['owned_unknown']==u['counts']['owned_unknown_pairs']
        row=by[v];ccs=row['certificates'];masses=[defaultdict(F),defaultdict(F)]
        multiplicities=[defaultdict(F),defaultdict(F)];terms=[defaultdict(list),defaultdict(list)]
        cover_stats=[]
        for k,c in enumerate(ccs):
            subtotal=defaultdict(F);support=Counter()
            for term,(vec,amount)in enumerate(c['mixture']):
                lam=F(amount);clique=row['literal_representatives'][vec];keys=Counter()
                for ii,i in enumerate(clique):
                    for j in clique[ii+1:]:
                        pair=tuple(sorted((i,j)));assert g.has_edge(*pair)
                        keys[pairkey[pair]]+=1
                kinds={key[0]for key in keys}
                for kind in kinds:subtotal[kind]+=lam;support[kind]+=1
                for key,n in keys.items():
                    masses[k][key]+=lam;multiplicities[k][key]+=lam*n
                    terms[k][key].append(term)
            cover_stats.append(dict(optimal_mass=c['optimal_mass'],vertex_mass=str(mus[k].get(v,F())),
                slack=str(mus[k].get(v,F())-F(c['optimal_mass'])),
                active=mus[k].get(v,F())==F(c['optimal_mass'])and F(c['optimal_mass'])>0,
                terms=len(c['mixture']),terms_using_unowned=support['unowned'],
                mass_using_unowned=str(subtotal['unowned']),terms_using_owned_unknown=support['owned_unknown'],
                mass_using_owned_unknown=str(subtotal['owned_unknown'])))
        ids=[]
        for key,pairs in sorted(domains.items()):
            if not any(masses[k][key]for k in range(2)):continue
            kind,i,j=key;groups=[u['corner_orbits'][x]for x in(i,j)]
            idx=len(candidates);ids.append(idx)
            encoded=json.dumps(pairs,separators=(',',':')).encode()
            candidates.append(dict(index=idx,vertex=v,forest=u['forest'],kind=kind,corner_orbits=[i,j],
                face_coefficients=[x['positive_face_coefficients']for x in groups],
                endpoint_demands=[[u['demands'][k][i],u['demands'][k][j]]for k in range(2)],
                both_endpoint_demands_positive=[all(F(u['demands'][k][x])>0 for x in(i,j))for k in range(2)],
                complete_domain_pairs=len(pairs),pair_domain_sha256=sha256(encoded).hexdigest(),
                covering_mass=[str(masses[k][key])for k in range(2)],
                edge_multiplicity_mass=[str(multiplicities[k][key])for k in range(2)],
                affected_cover_terms=[terms[k][key]for k in range(2)],
                active_vertex=[x['active']for x in cover_stats],
                covered_mass_per_pair=[str(masses[k][key]/len(pairs))for k in range(2)]))
        rows.append(dict(vertex=v,forest=u['forest'],counts=dict(owner_counts),certificates=cover_stats,
            supported_candidate_domains=ids))
        totals.update(forests=1,unowned_pairs=owner_counts['unowned'],owned_unknown_pairs=owner_counts['owned_unknown'],
            supported_domains=len(ids),full_active=int(cover_stats[0]['active']),restricted_active=int(cover_stats[1]['active']))
    rankings={}
    for k,label in enumerate(('full','restricted')):
        for kind in('unowned','owned_unknown'):
            pool=[c for c in candidates if c['kind']==kind and c['active_vertex'][k]
                and c['both_endpoint_demands_positive'][k]and F(c['covering_mass'][k])>0]
            bymass=sorted(pool,key=lambda c:(-F(c['covering_mass'][k]),c['complete_domain_pairs'],c['vertex'],c['corner_orbits']))
            byeff=sorted(pool,key=lambda c:(-F(c['covered_mass_per_pair'][k]),-F(c['covering_mass'][k]),c['vertex'],c['corner_orbits']))
            rankings[label+'_'+kind]=dict(by_supported_mass=[c['index']for c in bymass],by_mass_per_pair=[c['index']for c in byeff])
    out=dict(status='PASS_HISTORICAL22_EXACT_COVER_PAIR_PRIORITY_DIAGNOSTIC',sources=sources,
        counts=dict(totals),forest_records=rows,candidate_domains=candidates,rankings=rankings,
        fixed_repairs=cones['repairs'],
        interpretation='Each candidate is the COMPLETE compatible literal-pair domain of the specified orbit pair and ownership kind. Scores count exact attaining-mixture mass using that domain, once per term; they are diagnostic dependencies, not predicted improvements. Rankings require positive demand at both endpoint orbits: a zero-demand orbit can be deleted from every covering clique, so its incident-edge exclusions cannot change the projected fixed-demand cone. Saved representatives need not be unique.',
        scope='Only75 expanded forests. No new physical pair tests, overwritten first-owner outcomes, clique enumeration, optimization, constraint, repaired witness or dimension theorem. Owned unknowns are reported separately and cannot be overwritten by a later inherited domain under the current protocol.')
    out=normalize(out);path=ROOT/f'results/{TAG}.json'
    if verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(out['status'],dict(totals),'seconds',round(time.monotonic()-start,2),flush=True)
    for label,r in rankings.items():
        print(label,'MASS',[(candidates[i]['vertex'],candidates[i]['corner_orbits'],candidates[i]['complete_domain_pairs'],candidates[i]['covering_mass'])for i in r['by_supported_mass'][:5]],flush=True)
        print(label,'CHEAP',[(candidates[i]['vertex'],candidates[i]['corner_orbits'],candidates[i]['complete_domain_pairs'],candidates[i]['covering_mass'])for i in r['by_mass_per_pair'][:5]],flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');x=a.parse_args();run(x.verify)
