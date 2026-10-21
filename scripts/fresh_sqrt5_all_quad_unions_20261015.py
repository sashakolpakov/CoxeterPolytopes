#!/usr/bin/env python3
"""Resumable ALL marked Qsqrt5 quad minimal-piece unions, including short paths.

The validated catalogue/corner arithmetic is inherited explicitly. New traversal
stores supports and the exact unfinished DFS stack in one SQLite transaction.
Only finite-edge cycle conflicts are pruned early. Complete finite-contact
leaves receive the separately reviewed physical/conjugate principal tests.
"""
import argparse,json,sqlite3,time
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_quad_exhaustion_20261015 as E
import fresh_sqrt5_quad_filters_20261015 as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_all_quad_unions_20261015'
ALPHABET=(2,3,4,5,6,10)
DB=ROOT/'data/enumerations'/f'{TAG}_running.sqlite';REPORT=ROOT/'results'/f'{TAG}.json'

@lru_cache(maxsize=30000)
def finite_word(n,word):return E.finite_word.__wrapped__(n,word)

def corners(n,g):
    return all(finite_word(len(ids),tuple(g.get(E.pair(i,j),2)for i,j in combinations(ids,2)))
               for a,b in product((0,2),(1,3))for ids in [tuple(sorted(list(range(4,n))+[a,b]))])

def cycle_field(n,g):
    # Dotted squareclasses are unknown. Ignoring those edges is a necessary
    # relaxation; every cycle consisting only of assigned finite edges is fixed.
    parent=list(range(n));off=[0]*n
    def find(i):
        k=0
        while parent[i]!=i:k^=off[i];i=parent[i]
        return i,k
    for (i,j),m in g.items():
        if m in(0,2):continue
        x,a=find(i);y,b=find(j);delta=a^b^A.BITS[m]
        if x==y:
            if delta:return False
        else:parent[x]=y;off[x]=delta
    return True

def prepare(bank,li,ri,mp):
    lp,rp=bank[li],bank[ri];lm=[0,2]+list(range(4,lp['order']+2));rm=[1,3]+list(mp);g={}
    for rec,mapping in ((lp,lm),(rp,rm)):
        for (i,j),m in rec['edges'].items():
            p=E.pair(mapping[i],mapping[j])
            if p in g and g[p]!=m:return None
            g[p]=m
    n=max(lm+rm)+1
    if not corners(n,g)or not cycle_field(n,g):return None
    # In contrast to the frozen no-short enumeration, adjacent black contacts
    # are NOT forced right. Every unspecified pair takes the full alphabet.
    free=[p for p in combinations(range(n),2)if p not in g]
    return n,g,free,lm,rm

def manifest():
    files=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_quad_exhaustion_20261015.py',
           'results/fresh_sqrt5_quad_exhaustion_20261015.json',
           'scripts/fresh_sqrt5_quad_filters_20261015.py','results/fresh_sqrt5_quad_filters_20261015.json',
           'results/fresh_sqrt5_quad_filters_20261015_review.json',
           'results/fresh_sqrt5_lanner_catalogue_20261015_review.json']
    return {p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files}

def initialize(resume):
    DB.parent.mkdir(parents=True,exist_ok=True)
    if not resume:assert not DB.exists(),'Database exists: use --resume; no saved frontier is overwritten.'
    con=sqlite3.connect(DB)
    con.execute('PRAGMA journal_mode=DELETE');con.execute('PRAGMA synchronous=FULL')
    if not resume:
        con.execute('CREATE TABLE meta (id INTEGER PRIMARY KEY CHECK(id=1), value TEXT NOT NULL)')
        con.execute('CREATE TABLE examined (key BLOB PRIMARY KEY,status TEXT NOT NULL,payload TEXT)')
        state=dict(sources=manifest(),left=0,right=0,mapping=0,stack=None,counts={},batches=[],complete=False)
        con.execute('INSERT INTO meta VALUES(1,?)',(json.dumps(state),));con.commit()
    else:
        state=json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0]);assert state['sources']==manifest()
    return con,state

def save(con,state,counts):
    state['counts']=dict(counts);con.execute('UPDATE meta SET value=? WHERE id=1',(json.dumps(state,sort_keys=True),));con.commit()

def current(bank,state):
    li,ri,mi=state['left'],state['right'],state['mapping']
    mappings=list(E.maps(list(range(4,bank[li]['order']+2)),bank[ri]['order']-2))
    return mappings,prepare(bank,li,ri,mappings[mi])if mi<len(mappings)else None

def leaf_status(n,g,lm,rm):
    if not E.connected_pieces(lm,rm,g):return 'REJECT_ORTHOGONAL_LANNER_PAIR'
    labels=[g[p]for p in combinations(range(n),2)]
    if 0 in labels:return 'PASS_DOTTED_PARAMETER_UNFILTERED'
    row=dict(nodes=n,white=list(range(4,n)),labels=labels)
    color=A.coloring(row);assert color['passes']
    checks=A.conditions(row,A.scaled_matrix(row,color['colors']))
    if not checks['physical']['passes']:return 'REJECT_TWO_NEGATIVE_DIRECTIONS'
    if not checks['conjugate']['passes']:return 'REJECT_CONJUGATE_NOT_PSD'
    return 'PASS_FINITE_PRINCIPAL_NECESSARY_TESTS'

def traversal(con,state,bank,args):
    counts=Counter(state['counts']);begin=counts['dfs_nodes'];startleft=state['left'];started=time.monotonic();last_save=begin
    accepted=con.execute("SELECT count(*) FROM examined WHERE status LIKE 'PASS_%'").fetchone()[0]
    reason=None
    while state['left']<len(bank):
        if counts['dfs_nodes']-begin>=args.max_nodes:reason='NODE_BUDGET';break
        if accepted>=args.max_supports:reason='SUPPORT_STORAGE_GUARD';break
        if state['left']-startleft>=args.max_left:reason='LEFT_PIECE_BATCH_LIMIT';break
        li,ri,mi=state['left'],state['right'],state['mapping']
        if ri>=len(bank):
            state.update(left=li+1,right=li+1,mapping=0,stack=None)
            counts['completed_left_pieces']+=1;save(con,state,counts)
            print('completed left',li,'visits',counts['dfs_nodes'],'accepted',accepted,flush=True);continue
        mappings,prep=current(bank,state)
        if mi>=len(mappings):state.update(right=ri+1,mapping=0,stack=None);continue
        if state['stack'] is None:
            counts['base_identifications']+=1
            if prep is None:state['mapping']+=1;continue
            counts['base_corner_and_field_pass']+=1;state['stack']=[[0,0]];counts['dfs_nodes']+=1
        assert prep is not None
        n,g,free,lm,rm=prep;stack=state['stack']
        # Parent next-label cursors encode the branch leading to each child.
        for depth in range(len(stack)-1):
            assert stack[depth][0]==depth and 1<=stack[depth][1]<=len(ALPHABET)
            g[free[depth]]=ALPHABET[stack[depth][1]-1]
        while stack:
            if counts['dfs_nodes']-begin>=args.max_nodes or accepted>=args.max_supports:break
            k,next_label=stack[-1];assert k==len(stack)-1
            if k==len(free):
                counts['labelled_completions']+=1;word=tuple(g[p]for p in combinations(range(n),2));key=bytes((n,)+word)
                previous=con.execute('SELECT status FROM examined WHERE key=?',(key,)).fetchone()
                if previous:counts['duplicate_complete_union']+=1
                else:
                    status=leaf_status(n,g,lm,rm);counts[status]+=1
                    payload=None
                    if status.startswith('PASS_'):
                        row=dict(nodes=n,white=list(range(4,n)),labels=list(word),left_piece=li,right_piece=ri,
                                 left_nodes=lm,right_nodes=rm,shared_white=sorted(set(lm)&set(rm)),
                                 added_contacts=[[a,b,g[a,b]]for a,b in free if g[a,b]!=2])
                        payload=json.dumps(row,separators=(',',':'));accepted+=1
                    con.execute('INSERT INTO examined VALUES(?,?,?)',(key,status,payload))
                stack.pop()
                if k:g.pop(free[k-1])
            elif next_label==len(ALPHABET):
                stack.pop()
                if k:g.pop(free[k-1])
            else:
                stack[-1][1]+=1;m=ALPHABET[next_label];p=free[k];g[p]=m;counts['label_trials']+=1
                if not corners(n,g):counts['corner_PD_prunes']+=1;g.pop(p)
                elif not cycle_field(n,g):counts['finite_cycle_field_prunes']+=1;g.pop(p)
                else:stack.append([k+1,0]);counts['dfs_nodes']+=1
            if counts['dfs_nodes']-last_save>=5000:
                save(con,state,counts);last_save=counts['dfs_nodes']
                print('cursor',li,ri,mi,'visits',counts['dfs_nodes'],'accepted',accepted,flush=True)
        if not stack:state.update(mapping=mi+1,stack=None)
    if state['left']==len(bank):state['complete']=True;reason='COMPLETE'
    state['batches'].append(dict(initial_nodes=begin,final_nodes=counts['dfs_nodes'],start_left=startleft,end_left=state['left'],reason=reason))
    save(con,state,counts)
    return reason,time.monotonic()-started

def validate(con,state,bank,full=False):
    assert state['sources']==manifest()
    assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    counts=Counter(dict(con.execute('SELECT status,count(*) FROM examined GROUP BY status')))
    for status,n in counts.items():assert n==state['counts'][status]
    accepted=[]
    for key,status,payload in con.execute("SELECT key,status,payload FROM examined WHERE status LIKE 'PASS_%' ORDER BY key"):
        row=json.loads(payload);n=row['nodes'];assert key==bytes((n,)+tuple(row['labels']))
        assert row['white']==list(range(4,n));g=dict(zip(combinations(range(n),2),row['labels']))
        assert corners(n,g)and cycle_field(n,g)and E.connected_pieces(row['left_nodes'],row['right_nodes'],g)
        for rec,mp in ((bank[row['left_piece']],row['left_nodes']),(bank[row['right_piece']],row['right_nodes'])):
            assert len(mp)==rec['order']and len(set(mp))==len(mp)
            for(i,j),m in rec['edges'].items():assert g[E.pair(mp[i],mp[j])]==m
        if full:assert leaf_status(n,g,row['left_nodes'],row['right_nodes'])==status
        accepted.append(row)
    if state['stack'] is not None:
        _,prep=current(bank,state);assert prep is not None
        n,g,free,lm,rm=prep
        for depth,frame in enumerate(state['stack']):
            assert frame[0]==depth and 0<=frame[1]<=len(ALPHABET)
            if depth+1<len(state['stack']):
                assert frame[1]>0;g[free[depth]]=ALPHABET[frame[1]-1]
                assert corners(n,g)and cycle_field(n,g)
    return dict(counts),accepted

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');ap.add_argument('--verify',action='store_true')
    ap.add_argument('--max-nodes',type=int,default=100000);ap.add_argument('--max-supports',type=int,default=50000);ap.add_argument('--max-left',type=int,default=1)
    args=ap.parse_args();bank=E.pieces();assert len(bank)==153
    if args.verify:
        con=sqlite3.connect(f'file:{DB}?mode=ro',uri=True);state=json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0])
        statuses,rows=validate(con,state,bank,True);saved=json.loads(REPORT.read_text())
        assert saved['state']==state and saved['distinct_status_counts']==statuses and saved['labelled_union_supports']==rows
        assert saved['database_sha256']==sha256(DB.read_bytes()).hexdigest()
        print('PASS_SAVED_FRONTIER_AND_ALL_RETAINED_UNION_ROWS',len(rows),state['counts']['dfs_nodes']);return
    con,state=initialize(args.resume);reason,elapsed=traversal(con,state,bank,args)
    statuses,rows=validate(con,state,bank);con.close()
    out=dict(status='COMPLETE_ALL_QUAD_UNIONS'if state['complete']else'PARTIAL_RESUMABLE_ALL_QUAD_UNIONS',sources=manifest(),
             database=str(DB.relative_to(ROOT)),database_sha256=sha256(DB.read_bytes()).hexdigest(),state=state,
             piece_catalogue=list(map(E.encode_piece,bank)),distinct_status_counts=statuses,labelled_union_supports=rows,
             scope='No no-short restriction. Full finite alphabet and shared-white identifications. Partial finite-cycle-field pruning, terminal exact principal filters; dotted parameters unfiltered. External white attachments and any dimension claim outside scope.',
             cursor_semantics='SQLite commits result rows and exact stack/counters atomically. Parent next-label values encode all unvisited siblings; resume checks source hashes. Only completed outer pieces count as completed.')
    REPORT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],reason,'elapsed_seconds',round(elapsed,3),'accepted',len(rows),'counts',state['counts'])
if __name__=='__main__':main()
