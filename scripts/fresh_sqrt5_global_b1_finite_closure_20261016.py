#!/usr/bin/env python3
"""All finite B1 white-attachment domains, conditional only on delta<=1/2.
No current-payment, B0 no-short, E6-orbit or small-rank pruning is used.
"""
import argparse,json,sqlite3,time
from collections import Counter
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_global_b1_finite_first_20261016 as F
import fresh_sqrt5_all_quad_contacts_20261015 as D
import fresh_sqrt5_adjacent_discount_frontier_20261016 as W
ROOT=F.ROOT;J=F.J
TAG='fresh_sqrt5_global_b1_finite_closure_20261016'
INPUT='results/fresh_sqrt5_global_b1_finite_first_20261016.json'
DB=ROOT/f'data/enumerations/{TAG}.sqlite';OUT=ROOT/f'results/{TAG}.json'


def manifest():
    paths=[f'scripts/{TAG}.py',INPUT,f'scripts/{F.TAG}.py',f'scripts/{D.TAG}.py',
           f'scripts/{W.TAG}.py',f'scripts/{J.TAG}.py']
    return {p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths}


def evaluate(row):
    physical=F.classify(row);assert physical['status']in('ACCEPT_SINGULAR_STOP','OPEN_NONSINGULAR')
    analysis=W.analysis(row);assert analysis['baseline']==1 and analysis['potentially_affected']
    singular=physical['status']=='ACCEPT_SINGULAR_STOP'
    out=dict(classification=physical,analysis=analysis,singular=singular,
        orthogonal_spectator_rank_budget=23-row['nodes'],trials=[],accepted=[])
    if row['nodes']==23:
        assert singular
        return dict(**out,status='RANK23_SINGULAR_STOP',domain_counts=dict(profiles=0,rank_boundary=True))
    domain,counters=D.domains(row);trials=[];accepted=[]
    for contact in domain:
        child=D.extension(row,contact);ca=W.analysis(child);W.transition_check(analysis,ca)
        if ca['baseline']>=2:
            trial=dict(status='HEREDITARY_SAFE_BASELINE_AT_LEAST_TWO',baseline=ca['baseline'])
        elif not ca['potentially_affected']:
            trial=dict(status='HEREDITARY_SAFE_NO_POSSIBLE_ADJACENT_DISCOUNT',baseline=ca['baseline'])
        else:
            assert ca['baseline']==1
            trial=F.classify(child)
            if trial['status']in('ACCEPT_SINGULAR_STOP','OPEN_NONSINGULAR'):accepted.append(contact)
        trials.append(dict(contacts=contact,**trial))
    out['trials']=trials;out['accepted']=accepted
    status=('EXPANDED_SINGULAR'if singular else'EXPANDED_NONSINGULAR')if accepted else('SEALED_SINGULAR_STOP'if singular else'REJECT_SEALED_NONSINGULAR_RANK')
    return dict(**out,status=status,domain_counts=dict(**counters,profiles=len(domain)))


class Registry(F.A.Registry):
    pass


def export(con,state):
    records=[];counts=Counter();tc=Counter();obligations=[]
    for i,s,payload in con.execute('SELECT id,row,result FROM cores ORDER BY id'):
        result=json.loads(payload)if payload else None;row=json.loads(s)
        assert(result is not None)==(i<state['next_core'])
        records.append(dict(index=i,row=row,result=result))
        if result:
            counts[result['status']]+=1;tc.update(t['status']for t in result['trials'])
            if result['singular']:obligations.append(i)
        else:counts['PENDING']+=1
    complete=state['next_core']==len(records)
    doc=dict(status='COMPLETE_FINITE_B1_ATTACHMENT_CLOSURE'if complete else'PARTIAL_ATOMIC_FINITE_B1_ATTACHMENT_CLOSURE',
        sources=state['sources'],state=state,dimension=21,maximum_support_nodes=23,
        conditional_weight_hypotheses=dict(distance2_discount='0',adjacent_discount='0 <= delta <= 1/2',bonuses='nonnegative'),
        database=str(DB.relative_to(ROOT)),database_sha256=sha256(DB.read_bytes()).hexdigest(),
        marked_states=len(records),complete_queue=complete,counts=dict(counts),trial_counts=dict(tc),
        root_maps=[json.loads(r[0])for r in con.execute('SELECT payload FROM roots ORDER BY literal_index')],
        singular_stopping_obligations=obligations,records=records,
        scope='Complete finite-corner attachment domains over labels2,3,4,5,6,10 from all1646 finite first-extension classes. Singular states are stopping obligations AND expanded states; they are not terminal merely because their current rank is correct. Every nonorthogonal attachment is tested without B0 no-short pruning. Baseline>=2 and permanently ineligible branches are safely omitted only under the displayed weight hypotheses; current E6 exemption or current payment is never a pruning reason. Exact physical/conjugate signatures, cyclic integrality, root-relation signs and explicit full marked isomorphism maps are checked. At rank23 every retained state is singular. Common orthogonal spectators consume23-n ranks. If partial, pending states are explicit. The dotted branch and actual weight evaluation remain separate; no dimension or compact-realizability claim.')
    OUT.write_text(json.dumps(doc,indent=2)+'\n');return doc


def verify():
    data=json.loads(OUT.read_text());assert data['sources']==manifest()
    assert sha256(DB.read_bytes()).hexdigest()==data['database_sha256']
    con=sqlite3.connect(f'file:{DB}?mode=ro',uri=True)
    assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    assert json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0])==data['state']
    dbrows=list(con.execute('SELECT id,row,result FROM cores ORDER BY id'));con.close()
    assert len(dbrows)==len(data['records'])
    prior=json.loads((ROOT/INPUT).read_text());reached=set()
    for m in data['root_maps']:
        source=prior['records'][m['literal_index']]['row'];target=data['records'][m['target']]['row']
        J.check_map(J.matrix(source),J.matrix(target),m['map']);reached.add(m['target'])
    assert len(data['root_maps'])==len(prior['records'])==1646
    counts=Counter();tc=Counter();obligations=[]
    for rec,stored in zip(data['records'],dbrows):
        i=rec['index'];assert stored[0]==i and json.loads(stored[1])==rec['row']and i in reached
        assert(json.loads(stored[2])if stored[2]else None)==rec['result']
        if rec['result']is None:counts['PENDING']+=1;continue
        old=rec['result'];new=evaluate(rec['row']);assert new=={k:v for k,v in old.items()if k!='children'}
        counts[new['status']]+=1;tc.update(t['status']for t in new['trials'])
        if new['singular']:obligations.append(i)
        children={tuple(c['contacts']):c for c in old['children']}
        assert set(children)==set(map(tuple,new['accepted']))and len(children)==len(old['children'])
        for contact in new['accepted']:
            ch=children[tuple(contact)];child=D.extension(rec['row'],contact)
            J.check_map(J.matrix(child),J.matrix(data['records'][ch['index']]['row']),ch['map']);reached.add(ch['index'])
        if i%100==0:print('FINITE_B1_REPLAY',i,flush=True)
    assert reached==set(range(len(data['records'])))and dict(counts)==data['counts']and dict(tc)==data['trial_counts']
    assert obligations==data['singular_stopping_obligations']
    print('PASS_SAVED_FINITE_B1_ATTACHMENT_CLOSURE',data['state']['next_core'],data['marked_states'],data['complete_queue'],flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');ap.add_argument('--verify',action='store_true')
    ap.add_argument('--max-cores',type=int,default=10000);ap.add_argument('--max-seconds',type=float,default=600);a=ap.parse_args()
    if a.verify:verify();return
    sources=manifest()
    if not a.resume:assert not DB.exists(),'Use --resume; do not overwrite an existing queue.'
    con=sqlite3.connect(DB);con.execute('PRAGMA journal_mode=DELETE');con.execute('PRAGMA synchronous=FULL')
    if not a.resume:
        con.execute('BEGIN IMMEDIATE');con.execute('CREATE TABLE cores(id INTEGER PRIMARY KEY,row TEXT,result TEXT)')
        con.execute('CREATE TABLE roots(literal_index INTEGER PRIMARY KEY,payload TEXT)');con.execute('CREATE TABLE meta(id INTEGER PRIMARY KEY,value TEXT)')
        registry=Registry(con);source=json.loads((ROOT/INPUT).read_text())
        for rec in source['records']:
            i=rec['index'];j,mp=registry.add(rec['row'])
            con.execute('INSERT INTO roots VALUES(?,?)',(i,json.dumps(dict(literal_index=i,target=j,map=mp))))
        assert registry.next_id==1646
        state=dict(sources=sources,next_core=0,initial_marked_roots=1646,batches=[])
        con.execute('INSERT INTO meta VALUES(1,?)',(json.dumps(state),));con.commit()
    else:
        state=json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0]);assert state['sources']==sources
        registry=Registry(con)
    start=state['next_core'];began=time.monotonic();stop='MAX_CORES'
    for k in range(a.max_cores):
        if k and time.monotonic()-began>=a.max_seconds:stop='TIME_GUARD_AFTER_COMPLETE_CORE';break
        i=state['next_core'];s=con.execute('SELECT row FROM cores WHERE id=?',(i,)).fetchone()
        if s is None:stop='QUEUE_EXHAUSTED';break
        row=json.loads(s[0]);result=evaluate(row);con.execute('BEGIN IMMEDIATE');children=[]
        for contact in result['accepted']:
            j,mp=registry.add(D.extension(row,contact));children.append(dict(index=j,contacts=contact,map=mp))
        result['children']=children;con.execute('UPDATE cores SET result=? WHERE id=?',(json.dumps(result,separators=(',',':')),i))
        state['next_core']=i+1;con.execute('UPDATE meta SET value=? WHERE id=1',(json.dumps(state),));con.commit()
        if i%25==0:print('FINITE_B1_CORE',i,'nodes',row['nodes'],'states',registry.next_id,'singular',result['singular'],flush=True)
    state['batches'].append(dict(start=start,end=state['next_core'],stop_reason=stop,elapsed_seconds=round(time.monotonic()-began,3)))
    con.execute('UPDATE meta SET value=? WHERE id=1',(json.dumps(state),));con.commit();out=export(con,state);con.close()
    print(out['status'],state['next_core'],out['marked_states'],len(out['singular_stopping_obligations']),flush=True)

if __name__=='__main__':main()
