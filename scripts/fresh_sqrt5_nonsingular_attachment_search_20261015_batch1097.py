#!/usr/bin/env python3
"""Replay committed frontier proofs without repeating factorial subcore search."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_nonsingular_attachment_search_20261015 as S
import fresh_sqrt5_sealed_subcore_closure_20261015_review as J

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_nonsingular_attachment_search_20261015_batch1097'
ARCHIVE=f'data/enumerations/{TAG}.json'
F,Q,C=S.F,S.Q,S.C

def affine(row,w):
    ids=w['nodes'];g=Q.edges(row)
    sub={(i,j):g[Q.pair(ids[i],ids[j])]for i,j in combinations(range(len(ids)),2)}
    M=Q.tree_scaled(len(ids),sub);assert M is not None
    assert Q.inertia(M)==[len(ids)-1,0,1]
    k=Q.kernel(M);assert len(k)==1
    assert all(x.sign()>0 for x in k[0])or all(x.sign()<0 for x in k[0])

def run():
    a=json.loads((ROOT/ARCHIVE).read_text());assert a['state']['next_core']==1097 and len(a['records'])==1281
    assert a['sources']==S.manifest()
    base,lookup,weight=S.context();targets={(r['family'],r['index']):r for r in lookup.values()}
    raw=[r for r in S.load('fresh_sqrt5_sealed_subcore_closure_20261015')['records']if r['status']=='OPEN_NONSINGULAR']
    assert len(raw)==len(a['root_maps'])==311
    reached=set();maps=0
    for i,(r,z)in enumerate(zip(raw,a['root_maps'])):
        assert z['literal_index']==i and all(z[k]==r[k]for k in('source_union_index','source_class_index','trial_index','contacts'))
        J.check_map(J.matrix(F.extension(base[r['source_union_index']],r['contacts'])),J.matrix(a['records'][z['target']]['row']),z['map'])
        reached.add(z['target']);maps+=1
    assert reached==set(range(199))
    counts=Counter();tc=Counter();domains=0;affines=0;subcores=0;singular=[]
    for r in a['records']:
        i=r['index'];row=r['row'];old=r['result'];assert i in reached
        assert (old is not None)==(i<1097)
        if old is None:counts['PENDING']+=1;continue
        M,pi,co=S.classify(row);assert pi==old['physical']and co==old['conjugate']and bool(pi[2])==old['singular']
        assert F.corners(row)and not F.short(row)and not F.OLD.h4_payments(row)
        counts[old['status']]+=1
        if old['singular']:singular.append(i)
        if old['status']=='REJECT_CORE_AFFINE':affine(row,old['witness']);affines+=1;continue
        if old['status']=='PAID_BY_CLOSED_SUBCORE':
            w=old['closed_match'];target=targets[w['family'],w['index']];assert w['sealed']==target['sealed']
            ids=w['nodes'];g=Q.edges(row)
            sub=dict(nodes=len(ids),white=list(range(4,len(ids))),labels=[g[Q.pair(x,y)]for x,y in combinations(ids,2)])
            J.check_map(J.matrix(sub),J.matrix(target['row']),w['to_target'])
            if w['sealed']:assert all(g[Q.pair(x,y)]==2 for x in set(row['white'])-set(ids)for y in ids)
            subcores+=1;continue
        domain,dc=F.profiles(row);assert dc==old['domain_counts']and domain==[t['contacts']for t in old['trials']];domains+=len(domain)
        passed=[]
        for t in old['trials']:
            ext=F.extension(row,t['contacts']);status=t['status'];tc[status]+=1
            if status=='PAID_H4':assert F.OLD.h4_payments(ext)==t['witness'];continue
            if status=='REJECT_AFFINE':affine(ext,t['witness']);affines+=1;continue
            G=Q.tree_scaled(ext['nodes'],Q.edges(ext))
            if status=='REJECT_FIELD':assert G is None;continue
            assert G is not None
            ker=Q.kernel([z[:-1]for z in G[:-1]])if old['singular']else[]
            products=[sum((x*G[j][-1]for j,x in enumerate(v)),Q.K())for v in ker]
            if status=='REJECT_RELATION':assert any(products)and [x.enc()for x in products]==t['products'];continue
            assert not any(products)
            p=Q.inertia(G);c=Q.inertia([[x.conj()for x in z]for z in G]);assert p==t['physical']and c==t['conjugate']
            if status=='REJECT_TWO_NEGATIVE':assert p[1]>1
            elif status=='REJECT_CONJUGATE':assert p[1]<=1 and c[1]>0
            else:
                assert status=='ACCEPT_UNPAID'and p[1]==1 and p[2]in(0,1)and c[1]==0 and c[2]==p[2]
                passed.append(t['contacts'])
        assert passed==old['accepted'];children=old.get('children',[]);assert len(children)==len(passed)
        for c,z in zip(passed,children):
            assert z['contacts']==c
            J.check_map(J.matrix(F.extension(row,c)),J.matrix(a['records'][z['index']]['row']),z['map'])
            reached.add(z['index']);maps+=1
        if old['singular']:
            pay=S.A.payment(row,*weight);assert pay==old['payment_if_no_further_contact']
            expected='EXPANDED_SINGULAR'if passed else'NEW_SEALED_PAID'if pay['paid']else'UNPAID_SEALED_FRONTIER'
        else:expected='EXPANDED_NONSINGULAR'if passed else'REJECT_SEALED_NONSINGULAR_RANK'
        assert old['status']==expected
        if i and i%100==0:print('checked',i,flush=True)
    assert reached==set(range(1281))
    paths=[f'scripts/{TAG}.py',ARCHIVE,f'scripts/{J.TAG}.py',f'results/{J.TAG}.json']+list(a['sources'])
    Q.manifest_check(paths)
    return dict(status='PASS_COMMITTED_1097_CORE_WITNESS_AND_COMPLETE_DOMAIN_REPLAY',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in dict.fromkeys(paths)},
                processed=1097,marked_states=1281,pending=184,counts=dict(counts),trial_counts=dict(tc),domains_rebuilt=domains,
                literal_maps_checked=maps,affine_witnesses_checked=affines,closed_subcore_witnesses_checked=subcores,
                singular_states=singular,interruption=a['interruption'],
                scope='Immutable read-only snapshot after SIGINT during the uncommitted1097th index. All311 seed maps, every committed complete one-white domain, all stored exact exclusion witnesses, payments and child maps are checked. Closed-subcore inclusion is checked by its literal map, without repeating the expensive search for that witness; stored counts of failed subcore attempts are not replayed. Shared reviewed Fraction arithmetic, not an independent field engine or the full unchanged producer replay. All1281 states and184 pending cores are preserved. No all-size closure or global dimension theorem.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],data['trial_counts'])
