#!/usr/bin/env python3
"""Independent small-block census plus transactional resume protocol audit.

The active producer database is never opened. Producer code is invoked only on
fresh temporary databases for protocol tests; the comparison census uses a
separate recursive traversal, Dynkin classifier and Fraction field arithmetic.
Full153-left-block enumeration remains a separate unfinished obligation.
"""
import argparse,contextlib,importlib.util,io,json,sqlite3,sys,tempfile
from collections import Counter
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
from types import SimpleNamespace
ROOT=Path(__file__).resolve().parents[1];BASE='fresh_sqrt5_all_quad_unions_20261015';TAG=BASE+'_review'
sys.path.insert(0,str(ROOT/'scripts'))
import fresh_sqrt5_all_quad_unions_20261015 as P
import fresh_sqrt5_quad_contacts_20261015_review as Q
ALPHABET=(2,3,4,5,6,10)

def key(n,g):return bytes((n,)+tuple(g[ij]for ij in combinations(range(n),2)))
def piece_tuple(r):return (r['order'],tuple(r['edges'][p]for p in combinations(range(r['order']),2)))
def complete_pieces():
    data=json.loads((ROOT/'results/fresh_sqrt5_lanner_catalogue_20261015.json').read_text());found={(2,(0,))}
    for row in data['order_three']['accepted']:
        for labels in set(permutations(row['labels'])):found.add((3,labels))
    for block in data['higher_orders']:
        n=block['order'];pairs=list(combinations(range(n),2))
        for row in block['accepted_classes']:
            original=dict(zip(pairs,row['labels']))
            for a,b in permutations(range(n),2):
                whites=[i for i in range(n)if i not in(a,b)]
                words=[]
                for w in permutations(whites):
                    order=(a,b)+w;words.append(tuple(original[Q.pair(order[i],order[j])]for i,j in pairs))
                found.add((n,min(words)))
    bank=P.E.pieces();assert len(bank)==153 and {piece_tuple(z)for z in bank}==found
    # Producer order is only a common index convention; all literal pieces
    # were independently rebuilt from unordered black choices and white orbits.
    ours=[dict(order=n,edges=dict(zip(combinations(range(n),2),word)))for n,word in map(piece_tuple,bank)]
    return ours

def overlap_maps(left_count,right_count):
    left=list(range(4,4+left_count));out=[]
    for shared in range(min(left_count,right_count)+1):
        for slots in combinations(range(right_count),shared):
            for targets in permutations(left,shared):
                partial=dict(zip(slots,targets));fresh=iter(range(4+left_count,4+left_count+right_count-shared))
                out.append(tuple(partial[i]if i in partial else next(fresh)for i in range(right_count)))
    assert len(set(out))==len(out)
    return sorted(out)

def cycle_ok(n,g):
    adj=[[]for _ in range(n)];bits={3:0,4:1,5:0,6:2,10:4}
    for(a,b),label in g.items():
        if label in(0,2):continue
        adj[a].append((b,bits[label]));adj[b].append((a,bits[label]))
    color={}
    for start in range(n):
        if start in color:continue
        color[start]=0;todo=[start]
        while todo:
            a=todo.pop()
            for b,e in adj[a]:
                want=color[a]^e
                if b in color:
                    if color[b]!=want:return False
                else:color[b]=want;todo.append(b)
    return True

def corners(n,g):return all(Q.finite(list(range(4,n))+[a,b],g)for a,b in product((0,2),(1,3)))
def prepare(bank,li,ri,mp):
    left,right=bank[li],bank[ri];L=[0,2]+list(range(4,left['order']+2));R=[1,3]+list(mp);g={}
    for rec,images in((left,L),(right,R)):
        for ij,label in rec['edges'].items():
            p=Q.pair(*(images[i]for i in ij))
            if p in g and g[p]!=label:return None
            g[p]=label
    n=max(L+R)+1
    if not corners(n,g)or not cycle_ok(n,g):return None
    return n,g,[p for p in combinations(range(n),2)if p not in g],L,R

def status(n,g,L,R):
    if not set(L)&set(R)and all(g[Q.pair(a,b)]==2 for a,b in product(L,R)):return 'REJECT_ORTHOGONAL_LANNER_PAIR'
    if 0 in g.values():return 'PASS_DOTTED_PARAMETER_UNFILTERED'
    M=Q.tree_scaled(n,g);assert M is not None
    physical=Q.inertia(M);assert physical[1]>=1
    if physical[1]>1:return 'REJECT_TWO_NEGATIVE_DIRECTIONS'
    conjugate=Q.inertia([[a.conj()for a in row]for row in M])
    if conjugate[1]:return 'REJECT_CONJUGATE_NOT_PSD'
    return 'PASS_FINITE_PRINCIPAL_NECESSARY_TESTS'

def independent_two_blocks(bank):
    rows={};counts=Counter()
    for li in range(2):
        for ri in range(li,len(bank)):
            for mp in overlap_maps(bank[li]['order']-2,bank[ri]['order']-2):
                counts['base_identifications']+=1;prep=prepare(bank,li,ri,mp)
                if prep is None:continue
                counts['base_pass']+=1;n,g,free,L,R=prep
                def visit(depth):
                    counts['dfs_nodes']+=1
                    if depth==len(free):
                        counts['leaves']+=1;k=key(n,g)
                        if k not in rows:rows[k]=status(n,g,L,R)
                        return
                    ij=free[depth]
                    # Reverse label order differs from the producer protocol.
                    for label in reversed(ALPHABET):
                        g[ij]=label
                        if corners(n,g)and cycle_ok(n,g):visit(depth+1)
                    del g[ij]
                visit(0)
        print('Independent completed left',li,'nodes',counts['dfs_nodes'],'distinct',len(rows),flush=True)
    return rows,dict(counts)

def database_rows(con,include_payload=False):
    if include_payload:return list(con.execute('SELECT key,status,payload FROM examined ORDER BY key'))
    return dict(con.execute('SELECT key,status FROM examined ORDER BY key'))
def snapshot(con,state):
    statuses=Counter(dict(con.execute('SELECT status,count(*) FROM examined GROUP BY status')))
    assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    for k,n in statuses.items():assert state['counts'][k]==n
    return dict(cursor=[state['left'],state['right'],state['mapping']],stack=state['stack'],counts=state['counts'].copy(),complete=state['complete'],statuses=dict(statuses))

def run_protocol(bank,folder,mode):
    P.DB=Path(folder)/f'{mode}.sqlite';con,state=P.initialize(False);pauses=[];counter=0
    while state['left']<2:
        budgets=[1,2,7,37,4999,5000,17,1000000]
        nodes=1000000 if mode=='whole'else budgets[min(counter,len(budgets)-1)]
        maximum=1 if mode=='resumed'and counter==2 else 1000000
        args=SimpleNamespace(max_nodes=nodes,max_supports=maximum,max_left=2-state['left'])
        with contextlib.redirect_stdout(io.StringIO()):reason,_=P.traversal(con,state,bank,args)
        pauses.append(dict(reason=reason,cursor=[state['left'],state['right'],state['mapping']],stack_depth=len(state['stack']or[]),nodes=state['counts']['dfs_nodes']))
        # Actual close/reopen exercises JSON decoding and source checks.
        con.close();con,state=P.initialize(True);snapshot(con,state);counter+=1;assert counter<30
    result=database_rows(con,True);snap=snapshot(con,state);con.close()
    return result,snap,pauses

class InjectedCrash(Exception):pass
class InterruptAfterCommit:
    def __init__(self,connection):self.c=connection;self.commits=0;self.thrown=False
    def execute(self,sql,*a):
        ans=self.c.execute(sql,*a)
        if self.commits and sql.startswith('INSERT INTO examined'):
            self.thrown=True;raise InjectedCrash()
        return ans
    def commit(self):self.c.commit();self.commits+=1

def crash_protocol(bank,folder):
    P.DB=Path(folder)/'crash.sqlite';con,state=P.initialize(False);wrapper=InterruptAfterCommit(con)
    try:
        with contextlib.redirect_stdout(io.StringIO()):P.traversal(wrapper,state,bank,SimpleNamespace(max_nodes=1000000,max_supports=1000000,max_left=2))
        assert False,'Expected injected post-INSERT, pre-commit interruption'
    except InjectedCrash:pass
    assert wrapper.thrown and wrapper.commits>=1
    con.close() # rolls back row insertion after last atomic result+cursor commit
    con,state=P.initialize(True);saved=snapshot(con,state)
    assert saved['counts']['dfs_nodes']>=5000 and state['left']<2
    with contextlib.redirect_stdout(io.StringIO()):P.traversal(con,state,bank,SimpleNamespace(max_nodes=1000000,max_supports=1000000,max_left=2-state['left']))
    result=database_rows(con,True);snap=snapshot(con,state);con.close()
    return result,snap,saved

def old_subset_check(rows):
    old=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text());counts=Counter()
    for row in old['labelled_union_supports']:
        if row['left_piece']>=2:continue
        n=row['nodes'];g=Q.edges(row);assert corners(n,g)
        if not cycle_ok(n,g):counts['correctly_pruned_fixed_cycle']+=1;continue
        k=key(n,g);assert k in rows
        expected=status(n,g,row['left_nodes'],row['right_nodes']);assert rows[k]==expected
        counts[expected]+=1
    return dict(counts)

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    bank=complete_pieces();maps=[]
    for a,b in product(range(4),repeat=2):
        actual=overlap_maps(a,b);assert set(actual)==set(P.E.maps(list(range(4,4+a)),b));maps.append(dict(left_whites=a,right_whites=b,count=len(actual)))
    separate,counts=independent_two_blocks(bank)
    with tempfile.TemporaryDirectory(prefix='sqrt5_quad_protocol_')as folder:
        whole,snap,wp=run_protocol(P.E.pieces(),folder,'whole');resumed,rs,rp=run_protocol(P.E.pieces(),folder,'resumed');crashed,cs,rollback=crash_protocol(P.E.pieces(),folder)
        assert whole==resumed==crashed and snap==rs==cs
        assert dict((k,v)for k,v,_ in whole)==separate
        assert counts['dfs_nodes']==snap['counts']['dfs_nodes']and counts['base_identifications']==snap['counts']['base_identifications']
        assert counts['leaves']==snap['counts']['labelled_completions']
        assert any(z['reason']=='NODE_BUDGET'and z['stack_depth']>1 for z in rp)
        assert any(z['reason']=='SUPPORT_STORAGE_GUARD'for z in rp)
    retained=old_subset_check(separate)
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}_pause_test.json',
           'scripts/fresh_sqrt5_quad_contacts_20261015_review.py','results/fresh_sqrt5_quad_contacts_20261015_review.json',
           'scripts/fresh_sqrt5_lanner_catalogue_20261015.py','results/fresh_sqrt5_lanner_catalogue_20261015.json']+list(P.manifest())
    paths=list(dict.fromkeys(paths));Q.manifest_check(paths)
    digest=sha256(b''.join(k+v.encode()+b'\n'for k,v in sorted(separate.items()))).hexdigest()
    return dict(status='PASS_INDEPENDENT_TWO_LEFT_BLOCKS_AND_TRANSACTIONAL_RESUME_PROTOCOL_FULL_SEARCH_PENDING',sources={p:sha(p)for p in paths},
        piece_count=len(bank),all_overlap_map_domains=maps,independent_counts=counts,distinct_status_counts=dict(Counter(separate.values())),canonical_key_status_sha256=digest,
        protocol=dict(uninterrupted=snap,resumed_pauses=rp,interrupted_rollback_snapshot=rollback,comparison='All key/status/payload rows, counters and final cursor identical after uninterrupted, multipause and injected rollback runs.'),
        prior_no_short_completed_subset=retained,
        scope='Stable source audit plus independent complete left0/1 replay only. Active mutable producer DB untouched/unbound. Full153-block completion, full saved frontier and all old no-short-class retention remain separate final obligations. External white attachments/rank completion/global dimension not claimed.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text());print('PASS saved independent protocol review')
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['independent_counts'],out['distinct_status_counts'])
