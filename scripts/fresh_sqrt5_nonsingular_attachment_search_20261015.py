#!/usr/bin/env python3
"""Resumable exact unpaid attachment frontier from311 nonsingular roots.

One whole marked core and all its outgoing profiles commit atomically.
No rank cutoff closes a branch. Nonempty pending queues and unpaid sealed
states are explicit. Complete would only refer to this seed family.
"""
import argparse,json,sqlite3,time
from collections import Counter
from fractions import Fraction as R
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_singular_attachment_closure_20261015 as A

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_nonsingular_attachment_search_20261015'
F,Q,C=A.F,A.Q,A.C
DB=ROOT/f'data/enumerations/{TAG}_running.sqlite';OUT=ROOT/f'results/{TAG}.json'
load=lambda p:json.loads((ROOT/f'results/{p}.json').read_text())

def manifest():
    names=['fresh_sqrt5_sealed_subcore_closure_20261015','fresh_sqrt5_quad_exhaustion_20261015',F.TAG,A.TAG,
           'fresh_sqrt5_quad_contacts_20261015','fresh_sqrt5_forest_sensitive_weights_20261015','fresh_sqrt5_extended_forest_weights_20261015','fresh_sqrt5_sealed_weights_20261015']
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py']+[f'results/{n}.json'for n in names]
    Q.manifest_check(paths)
    return{p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths}

def context():
    first=load(F.TAG);closure=load(A.TAG);old=load('fresh_sqrt5_quad_contacts_20261015');base=load('fresh_sqrt5_quad_exhaustion_20261015')['labelled_union_supports']
    lookup=A.known_catalogue(first,base,old)
    for r in lookup.values():r['sealed']=True
    for r in first['records']:
        if r['accepted']:
            key,order=C.canonical(r['row']);assert key not in lookup
            lookup[key]=dict(family='five_closed_starts',index=r['index'],row=r['row'],order=order,sealed=False)
    for r in closure['records']:
        if r['status']=='NEW_SEALED_PAID':
            key,order=C.canonical(r['row']);assert key not in lookup
            lookup[key]=dict(family='eleven_new_sealed',index=r['index'],row=r['row'],order=order,sealed=True)
    assert len(lookup)==174
    weight=[load(p)for p in('fresh_sqrt5_forest_sensitive_weights_20261015','fresh_sqrt5_sealed_weights_20261015','fresh_sqrt5_extended_forest_weights_20261015')]
    return base,lookup,weight

def classify(row):
    M=Q.tree_scaled(row['nodes'],Q.edges(row));assert M is not None
    pi=Q.inertia(M);co=Q.inertia([[x.conj()for x in r]for r in M]);n=row['nodes']
    assert pi[1]==1 and pi[2]in(0,1)and pi[0]+pi[2]==n-1
    assert co[1]==0 and co[2]==pi[2]
    return M,pi,co

def closed_match(row,lookup):
    g=Q.edges(row);sizes=sorted({k[0]-4 for k in lookup});checked=0
    for size in sizes:
        for white in combinations(row['white'],size):
            ids=list(range(4))+list(white);sub=dict(nodes=len(ids),white=list(range(4,len(ids))),labels=[g[Q.pair(i,j)]for i,j in combinations(ids,2)])
            key,order=C.canonical(sub);checked+=1
            if key not in lookup:continue
            r=lookup[key];mp=[None]*len(ids)
            for a,b in zip(order,r['order']):mp[a]=b
            C.check_map(sub,r['row'],mp)
            if r['sealed']:assert all(g[Q.pair(i,j)]==2 for i in set(row['white'])-set(white)for j in ids)
            return dict(family=r['family'],index=r['index'],nodes=ids,to_target=mp,sealed=r['sealed']),checked
    return None,checked

def evaluate(row,lookup,weight):
    n=row['nodes'];M,pi,co=classify(row);assert F.corners(row)and not F.short(row)and not F.OLD.h4_payments(row)
    out=dict(physical=pi,conjugate=co,singular=bool(pi[2]))
    if(w:=F.affine(row)):return dict(**out,status='REJECT_CORE_AFFINE',witness=w)
    if out['singular']:
        match,checked=closed_match(row,lookup);out.update(closed_subsets_checked=checked,closed_match=match)
        if match:return dict(**out,status='PAID_BY_CLOSED_SUBCORE')
    # A nonsingular Lorentz principal core cannot contain a singular
    # Lorentz core without introducing another negative direction. Its
    # complete closed-subcore search is therefore unnecessary.
    domain,dc=F.profiles(row);trials=[];accepted=[]
    for contact in domain:
        ext=F.extension(row,contact);r=dict(contacts=contact)
        if(w:=F.OLD.h4_payments(ext)):r.update(status='PAID_H4',witness=w)
        elif(w:=F.affine(ext)):r.update(status='REJECT_AFFINE',witness=w)
        elif(G:=Q.tree_scaled(n+1,Q.edges(ext)))is None:r['status']='REJECT_FIELD'
        else:
            ker=Q.kernel([z[:-1]for z in G[:-1]])if out['singular']else[]
            products=[sum((x*G[i][-1]for i,x in enumerate(v)),Q.K())for v in ker]
            if any(products):r.update(status='REJECT_RELATION',products=[x.enc()for x in products])
            else:
                physical=Q.inertia(G);conjugate=Q.inertia([[x.conj()for x in z]for z in G]);r.update(physical=physical,conjugate=conjugate)
                if physical[1]>1:r['status']='REJECT_TWO_NEGATIVE'
                elif conjugate[1]:r['status']='REJECT_CONJUGATE'
                else:
                    assert physical[1]==1 and physical[2]in(0,1)and physical[2]==conjugate[2]
                    r['status']='ACCEPT_UNPAID';accepted.append(contact)
        trials.append(r)
    out.update(domain_counts=dc,trials=trials,accepted=accepted)
    if out['singular']:
        out['payment_if_no_further_contact']=A.payment(row,*weight)
        if not accepted:out['status']='NEW_SEALED_PAID'if out['payment_if_no_further_contact']['paid']else'UNPAID_SEALED_FRONTIER'
        else:out['status']='EXPANDED_SINGULAR';out['zero_further_contact_unpaid']=not out['payment_if_no_further_contact']['paid']
    else:out['status']='EXPANDED_NONSINGULAR'if accepted else'REJECT_SEALED_NONSINGULAR_RANK'
    return out

def canon_key(row):
    key,order=C.canonical(row)
    return json.dumps(key,separators=(',',':')),list(order)
def add(con,row):
    key,order=canon_key(row);found=con.execute('SELECT id,row,ordering FROM cores WHERE key=?',(key,)).fetchone()
    if found is None:
        i=con.execute('SELECT COUNT(*) FROM cores').fetchone()[0]
        con.execute('INSERT INTO cores VALUES(?,?,?,?,NULL)',(i,key,json.dumps(row,separators=(',',':')),json.dumps(order)));target=row;target_order=order
    else:i,encoded,ordering=found;target=json.loads(encoded);target_order=json.loads(ordering)
    mp=[None]*row['nodes']
    for a,b in zip(order,target_order):mp[a]=b
    C.check_map(row,target,mp)
    return i,mp

def export(con,state):
    assert state['sources']==manifest()
    cores=[];counts=Counter();trials=Counter();unpaid=[]
    for i,row,payload in con.execute('SELECT id,row,result FROM cores ORDER BY id'):
        record=dict(index=i,row=json.loads(row),result=json.loads(payload)if payload else None);cores.append(record)
        if record['result']:
            r=record['result'];counts[r['status']]+=1;trials.update(t['status']for t in r.get('trials',[]))
            if r['status']=='UNPAID_SEALED_FRONTIER'or r.get('zero_further_contact_unpaid'):unpaid.append(i)
        else:counts['PENDING']+=1
    assert [r['index']for r in cores]==list(range(len(cores)))
    assert all((r['result']is not None)==(r['index']<state['next_core'])for r in cores)
    complete=state['next_core']==len(cores)
    out=dict(status='COMPLETE_SEED_FAMILY_FRONTIER'if complete else'PARTIAL_RESUMABLE_NONSINGULAR_ATTACHMENT_FRONTIER',sources=manifest(),state=state,
        database=str(DB.relative_to(ROOT)),database_sha256=sha256(DB.read_bytes()).hexdigest(),marked_states=len(cores),counts=dict(counts),trial_counts=dict(trials),
        complete_queue=complete,unpaid_stopping_cases=unpaid,root_maps=[json.loads(z[0])for z in con.execute('SELECT payload FROM roots ORDER BY literal_index')],records=cores,
        scope='Exact one-white necessary attachment search from311 literal nonsingular first extensions. Pending marked states and unpaid stopping cases are explicit. Whole cores and all outgoing profiles commit atomically. No rank cutoff is used to reject a branch; batch limits only pause. Known174 closed classes are inherited with literal maps. Complete queue alone would not establish payment if unpaid stopping cases remain. Dotted families and other seed families are outside this run; no global dimension bound.')
    OUT.write_text(json.dumps(out,indent=2)+'\n');return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');ap.add_argument('--verify',action='store_true');ap.add_argument('--max-cores',type=int,default=10);args=ap.parse_args()
    base,lookup,weight=context();sources=manifest()
    if args.verify:
        saved=json.loads(OUT.read_text());assert saved['sources']==sources
        con=sqlite3.connect(f'file:{DB}?mode=ro',uri=True);assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
        state=json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0]);assert state==saved['state']and sha256(DB.read_bytes()).hexdigest()==saved['database_sha256']
        assert state['sources']==sources
        initial=[r for r in load('fresh_sqrt5_sealed_subcore_closure_20261015')['records']if r['status']=='OPEN_NONSINGULAR'];assert len(initial)==len(saved['root_maps'])==311
        stored_roots=[json.loads(r[0])for r in con.execute('SELECT payload FROM roots ORDER BY literal_index')];assert stored_roots==saved['root_maps']
        root_targets=set()
        for i,(raw,mp)in enumerate(zip(initial,stored_roots)):
            assert mp['literal_index']==i and all(mp[k]==raw[k]for k in('source_union_index','source_class_index','trial_index','contacts'))
            C.check_map(F.extension(base[raw['source_union_index']],raw['contacts']),saved['records'][mp['target']]['row'],mp['map']);root_targets.add(mp['target'])
        assert root_targets==set(range(state['initial_marked_roots']))
        for row in saved['records']:
            stored=con.execute('SELECT row,result FROM cores WHERE id=?',(row['index'],)).fetchone();assert json.loads(stored[0])==row['row']
            if row['result']:
                result=evaluate(row['row'],lookup,weight);old=row['result'];children=old.get('children',[]);baseold={k:v for k,v in old.items()if k!='children'};assert result==baseold
                assert len(children)==len(result.get('accepted',[]))
                for z,c in zip(children,result.get('accepted',[])):
                    target=saved['records'][z['index']]['row'];C.check_map(F.extension(row['row'],c),target,z['map'])
                assert json.loads(stored[1])==old
            else:assert stored[1]is None
        con.close();print('PASS_SAVED_NONSINGULAR_ATTACHMENT_BATCH',state['next_core'],saved['marked_states']);return
    DB.parent.mkdir(parents=True,exist_ok=True)
    if not args.resume:assert not DB.exists(),'Use --resume; never overwrite a saved frontier.'
    con=sqlite3.connect(DB);con.execute('PRAGMA journal_mode=DELETE');con.execute('PRAGMA synchronous=FULL')
    if not args.resume:
        con.execute('BEGIN IMMEDIATE')
        con.execute('CREATE TABLE cores(id INTEGER PRIMARY KEY,key TEXT UNIQUE,row TEXT,ordering TEXT,result TEXT)')
        con.execute('CREATE TABLE roots(literal_index INTEGER PRIMARY KEY,payload TEXT)');con.execute('CREATE TABLE meta(id INTEGER PRIMARY KEY,value TEXT)')
        initial=[r for r in load('fresh_sqrt5_sealed_subcore_closure_20261015')['records']if r['status']=='OPEN_NONSINGULAR'];assert len(initial)==311
        for i,r in enumerate(initial):
            child=F.extension(base[r['source_union_index']],r['contacts']);j,mp=add(con,child)
            con.execute('INSERT INTO roots VALUES(?,?)',(i,json.dumps(dict(literal_index=i,source_union_index=r['source_union_index'],source_class_index=r['source_class_index'],trial_index=r['trial_index'],contacts=r['contacts'],target=j,map=mp))))
        state=dict(sources=sources,next_core=0,initial_literal_roots=311,initial_marked_roots=con.execute('SELECT COUNT(*) FROM cores').fetchone()[0],batches=[])
        con.execute('INSERT INTO meta VALUES(1,?)',(json.dumps(state),));con.commit()
    else:
        state=json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0]);assert state['sources']==sources
    start=state['next_core'];began=time.monotonic()
    for _ in range(args.max_cores):
        i=state['next_core'];record=con.execute('SELECT row FROM cores WHERE id=?',(i,)).fetchone()
        if record is None:break
        row=json.loads(record[0]);result=evaluate(row,lookup,weight);children=[]
        # No writes occur during expensive profile evaluation. This complete
        # result, all children, and the cursor share one SQLite transaction.
        con.execute('BEGIN IMMEDIATE')
        for contact in result.get('accepted',[]):
            j,mp=add(con,F.extension(row,contact));children.append(dict(index=j,contacts=contact,map=mp))
        if children:result['children']=children
        con.execute('UPDATE cores SET result=? WHERE id=?',(json.dumps(result,separators=(',',':')),i));state['next_core']=i+1
        con.execute('UPDATE meta SET value=? WHERE id=1',(json.dumps(state),));con.commit()
        print('core',i,result['status'],'children',len(children),'states',con.execute('SELECT COUNT(*) FROM cores').fetchone()[0],flush=True)
    state['batches'].append(dict(start=start,end=state['next_core']))
    con.execute('UPDATE meta SET value=? WHERE id=1',(json.dumps(state),));con.commit();out=export(con,state);con.close()
    print(out['status'],out['state']['next_core'],out['marked_states'],out['counts'],'seconds',round(time.monotonic()-began,3))

if __name__=='__main__':main()
