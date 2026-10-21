#!/usr/bin/env python3
"""Finite exact attachment closure of new singular roots, with explicit frontier.

Roots are seven continuations of the51-class stage and the16 longer A-chain
stars. Tests inherit the frozen exact helpers. Every reduction to a known
sealed core carries a literal induced-subdiagram map.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as R
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_singular_extension_frontier_20261015 as F

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_singular_attachment_closure_20261015'
Q,C,T=F.Q,F.C,F.T

def known_catalogue(first,base,old):
    lookup={}
    roots=[('old',r['class_index'],base[r['union_index']])for r in old['records']if r['status']=='UNPAID_FINITE_CONTACT_DOMAIN'and r['singular']and not r['accepted']]
    roots +=[('new',r['index'],r['row'])for r in first['records']if not r['accepted']]
    for family,ident,row in roots:
        key,order=C.canonical(row);assert key not in lookup
        lookup[key]=dict(family=family,index=ident,row=row,order=order)
    assert len(lookup)==158
    return lookup

def sealed_match(row,lookup):
    g=Q.edges(row);sizes=sorted({key[0]-4 for key in lookup});matches=[];checked=0
    for size in sizes:
        for white in combinations(row['white'],size):
            ids=list(range(4))+list(white);sub=dict(nodes=len(ids),white=list(range(4,len(ids))),labels=[g[Q.pair(i,j)]for i,j in combinations(ids,2)])
            key,order=C.canonical(sub);checked+=1
            if key not in lookup:continue
            target=lookup[key];mp=[None]*len(ids)
            for a,b in zip(order,target['order']):mp[a]=b
            C.check_map(sub,target['row'],mp)
            outside=set(row['white'])-set(white)
            assert all(g[Q.pair(i,j)]==2 for i in outside for j in ids)
            matches.append(dict(family=target['family'],index=target['index'],nodes=ids,to_target=mp,outside_white=sorted(outside)))
    return matches,checked

def next_domain(row):
    n=row['nodes'];G=Q.tree_scaled(n,Q.edges(row));assert G is not None
    assert Q.inertia(G)==[n-2,1,1]and Q.inertia([[x.conj()for x in r]for r in G])==[n-1,0,1]
    assert F.corners(row)and not F.short(row)and not F.OLD.h4_payments(row)and not F.affine(row)
    domain,counts=F.profiles(row);trials=[];accepted=[]
    for contact in domain:
        ext=F.extension(row,contact);r=dict(contacts=contact)
        if(w:=F.OLD.h4_payments(ext)):r.update(status='PAID_H4',witness=w)
        elif(w:=F.affine(ext)):r.update(status='REJECT_AFFINE',witness=w)
        elif(M:=Q.tree_scaled(n+1,Q.edges(ext)))is None:r['status']='REJECT_FIELD'
        else:
            ker=Q.kernel([z[:-1]for z in M[:-1]]);assert len(ker)==1
            products=[sum((x*M[i][-1]for i,x in enumerate(v)),Q.K())for v in ker]
            if any(products):r.update(status='REJECT_RELATION',products=[x.enc()for x in products])
            else:
                pi=Q.inertia(M);co=Q.inertia([[x.conj()for x in z]for z in M]);r.update(physical=pi,conjugate=co)
                if pi[1]>1:r['status']='REJECT_TWO_NEGATIVE'
                elif co[1]:r['status']='REJECT_CONJUGATE'
                else:
                    assert pi==[n-1,1,1]and co==[n,0,1];r['status']='ACCEPT_UNPAID';accepted.append(contact)
        trials.append(r)
    return dict(domain_counts=counts,trials=trials,accepted=accepted)

def payment(row,original,sealed,extended):
    out=F.payment(row,original,sealed);upgrade={r['key']:list(map(R,r['new']))for r in extended['upgraded_orbits']}
    for c in out['corners']:
        if c['pair_orbit']in upgrade:
            v=upgrade[c['pair_orbit']];c['weights']=[str(v[int(c['core_has_a1'])]),str(v[1])]
    sums=[sum(R(c['weights'][i])for c in out['corners'])for i in(0,1)]
    out['sums_without_with_external_a1']=list(map(str,sums));out['paid']=min(sums)>=1
    return out

def run():
    load=lambda tag:json.loads((ROOT/f'results/{tag}.json').read_text())
    first=load(F.TAG);base=load('fresh_sqrt5_quad_exhaustion_20261015')['labelled_union_supports'];old=load('fresh_sqrt5_quad_contacts_20261015')
    original=load('fresh_sqrt5_forest_sensitive_weights_20261015');sealed=load('fresh_sqrt5_sealed_weights_20261015');extended=load('fresh_sqrt5_extended_forest_weights_20261015')
    chain=load('fresh_sqrt5_single_chain_completion_20261015');lookup=known_catalogue(first,base,old)
    assert extended['status']=='PASS_EXACT_EXTENDED_FOREST_RULE_WITH_RANK29_BUDGET84'
    assert extended['old_sealed_cores']==112 and extended['new_sealed_cores']==46 and extended['face_cases']==316
    initial_stops=[]
    for r in first['records']:
        if r['accepted']:
            pay=payment(r['row'],original,sealed,extended);assert pay['paid']
            initial_stops.append(dict(index=r['index'],sums=pay['sums_without_with_external_a1']))
    assert len(initial_stops)==5
    records=[];keys={};seedmaps=[]
    def add(row,origin):
        key,order=C.canonical(row)
        if key not in keys:
            keys[key]=len(records);records.append(dict(index=len(records),row=row,order=order,origins=[]))
        rec=records[keys[key]];mp=[None]*row['nodes']
        for a,b in zip(order,rec['order']):mp[a]=b
        C.check_map(row,rec['row'],mp);rec['origins'].append(dict(**origin,to_representative=mp))
        return rec['index'],mp
    for r in first['records']:
        for contact in r['accepted']:
            i,mp=add(F.extension(r['row'],contact),dict(kind='seven_continuations',source_index=r['index'],contacts=contact));seedmaps.append(dict(kind='seven_continuations',target=i,map=mp))
    for k,r in enumerate(chain['records']):
        if r.get('length',0)>1:
            assert r['status']=='UNIQUE_LENGTH_UNPAID_PRINCIPAL_STAR';i,mp=add(r['row'],dict(kind='longer_A_chain',source_record=k,length=r['length']));seedmaps.append(dict(kind='longer_A_chain',target=i,map=mp))
    assert Counter(r['kind']for r in seedmaps)=={'seven_continuations':7,'longer_A_chain':16}
    processed=0;allsubsets=0;trialcounts=Counter()
    while processed<len(records)and processed<100:
        rec=records[processed];row=rec['row'];matches,checks=sealed_match(row,lookup);allsubsets+=checks;rec['known_subcores']=matches
        if matches:rec['status']='REDUCES_TO_KNOWN_158_SEALED_CATALOGUE'
        else:
            out=next_domain(row);rec.update(out);trialcounts.update(t['status']for t in out['trials'])
            rec['payment']=payment(row,original,sealed,extended)
            if not out['accepted']:rec['status']='NEW_SEALED_PAID'if rec['payment']['paid']else'NEW_SEALED_UNPAID'
            else:
                rec['status']='EXPANDED_NONSEALED';rec['children']=[]
                before=Q.components(Q.graph(row['white'],Q.edges(row)))
                for contact in out['accepted']:
                    child=F.extension(row,contact);i,mp=add(child,dict(kind='recursive',source_index=processed,contacts=contact))
                    after=Q.components(Q.graph(child['white'],Q.edges(child)))
                    rec['children'].append(dict(index=i,contacts=contact,map=mp,white_components_before=len(before),white_components_after=len(after)))
        print('state',processed,rec['status'],'pending',len(records)-processed-1,flush=True);processed+=1
    counts=Counter(r.get('status','UNPROCESSED_FRONTIER')for r in records);complete=processed==len(records)
    assert complete and processed==17 and allsubsets==3054
    assert dict(counts)==dict(REDUCES_TO_KNOWN_158_SEALED_CATALOGUE=6,NEW_SEALED_PAID=11)
    assert dict(trialcounts)==dict(REJECT_RELATION=115,REJECT_FIELD=22,PAID_H4=1)
    paths=[f'scripts/{TAG}.py',f'scripts/{F.TAG}.py',f'results/{F.TAG}.json']
    for p in('fresh_sqrt5_quad_exhaustion_20261015','fresh_sqrt5_quad_contacts_20261015','fresh_sqrt5_forest_sensitive_weights_20261015','fresh_sqrt5_extended_forest_weights_20261015','fresh_sqrt5_sealed_weights_20261015','fresh_sqrt5_single_chain_completion_20261015'):paths.append(f'results/{p}.json')
    Q.manifest_check(paths)
    return dict(status='PASS_COMPLETE_SINGULAR_ATTACHMENT_FRONTIER'if complete else'PASS_BOUNDED_SINGULAR_ATTACHMENT_FRONTIER',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        initial_literal_seeds=23,seed_maps=seedmaps,initial_five_zero_further_contact_cases=initial_stops,processed=processed,marked_states=len(records),complete=complete,counts=dict(counts),trial_counts=dict(trialcounts),white_subsets_checked=allsubsets,records=records,
        scope='Exact recursive unpaid one-white closure from seven continuation profiles and16 longer single-A-chain stars. Every nonzero next contact is enumerated; orthogonal spectators are accounted for by the two externalA1 contexts at sealed leaves. Known sealed subcores carry literal maps. Complete means this explicit queue is exhausted, not that311 nonsingular first extensions or dotted families are closed. No global dimension bound; a processed-state cap would be reported as an open frontier, not a theorem.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text());print('PASS saved attachment-closure replay')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],data['trial_counts'])
