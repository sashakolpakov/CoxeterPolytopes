#!/usr/bin/env python3
"""Bounded exact connected-white triangle extension BFS, with durable transactions.

This is a domain exploration, not complete dimension21 face coverage until the
total22-node frontier has actually been exhausted and independently reviewed.
"""
import argparse,json,sqlite3,time
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_sqrt5_triangle_weights_20261015 as T
import fresh_sqrt5_quad_filters_20261015 as Q

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_triangle_extension_frontier_20261016'
TRI='results/fresh_sqrt5_triangle_weights_20261015.json'
ALPHABET=(3,4,5,6,10)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def pack(x):return json.dumps(x,separators=(',',':'),sort_keys=True)
def edges(row):return [(i,j,m)for(i,j),m in zip(combinations(range(row['nodes']),2),row['labels'])if m!=2]
def from_edges(n,es):
    e={(min(i,j),max(i,j)):m for i,j,m in es}
    return dict(nodes=n,black=[0,1,2],white=list(range(3,n)),labels=[e.get(p,2)for p in combinations(range(n),2)])
def graph(row):
    g=nx.Graph();g.add_nodes_from((i,dict(color=int(i<3)))for i in range(row['nodes']));g.add_edges_from((i,j,dict(label=m))for i,j,m in edges(row));return g
def bucket(row):
    g=graph(row)
    return pack((row['nodes'],sorted((int(i<3),len(g[i]),tuple(sorted(g[i][j]['label']for j in g[i])))for i in g)))
def isomorphism(a,b):
    gm=nx.algorithms.isomorphism.GraphMatcher(graph(a),graph(b),node_match=lambda x,y:x['color']==y['color'],edge_match=lambda x,y:x['label']==y['label'])
    if not gm.is_isomorphic():return None
    mp=gm.mapping;assert set(mp)==set(range(a['nodes'])) and set(mp.values())==set(range(b['nodes']))
    ea={(i,j):m for i,j,m in edges(a)};eb={(i,j):m for i,j,m in edges(b)}
    assert all(ea.get((i,j),2)==eb.get(tuple(sorted((mp[i],mp[j]))),2)for i,j in combinations(range(a['nodes']),2))
    assert all((i<3)==(mp[i]<3)for i in mp)
    return [mp[i]for i in range(a['nodes'])]
def corner_info(row):
    es=edges(row);n=row['nodes'];out=[]
    for omitted in range(3):
        ids=[i for i in range(n)if i!=omitted];pair=[i for i in range(3)if i!=omitted]
        # T.graph expects edges already restricted to its literal node set.
        g=T.graph(ids,[(i,j,m)for i,j,m in es if i!=omitted and j!=omitted])
        comps=[];marked=None
        for comp in T.components(g):
            h=T.restrict(g,comp)
            try:name,extra=T.describe(h)
            except AssertionError:return None
            comps.append(dict(nodes=comp,type=name))
            if pair[0]in h and pair[1]in h:marked=h
        d=T.distance(g,*pair);base=int(d in (1,2));eligible=False;typ=None;rank=0
        if marked is not None:
            typ,extra=T.describe(marked);rank=len(marked)
            if typ=='H4':base+=extra.get(tuple(pair),0)
            eligible=rank>=6 and marked[pair[0]].get(pair[1])==3
        out.append(dict(omitted_black=omitted,vertices=ids,marked_pair=pair,components=comps,distance=d,baseline=base,eligible_adjacent3=eligible,marked_component_type=typ,marked_component_rank=rank))
    return dict(corners=out,baseline=sum(r['baseline']for r in out),eligible_count=sum(r['eligible_adjacent3']for r in out),black3_edges=sum(m==3 and j<3 for i,j,m in es))
def safe_reason(info):
    b,k=info['baseline'],info['black3_edges']
    if b>=4:return 'BASELINE_AT_LEAST4'
    if b==3 and k<2:return 'BASELINE3_FEWER_THAN_TWO_POTENTIAL_DISCOUNTS'
    if b==2 and k==0:return 'BASELINE2_NO_POTENTIAL_DISCOUNTS'
    return None
def arithmetic(row):
    col=Q.coloring(row)
    if not col['passes']:return dict(status='REJECT_CYCLE_FIELD',coloring=col)
    M=Q.scaled_matrix(row,col['colors']);corner=list(range(1,row['nodes']))
    try:S,ps=Q.schur(M,corner)
    except AssertionError:raise AssertionError('spherical classification disagrees with physical pivots')
    physical=S[0][0].sign();assert physical==-1
    C=[[x.conj()for x in r]for r in M]
    try:Sc,pc=Q.schur(C,corner)
    except AssertionError:return dict(status='REJECT_CONJUGATE_CORNER',colors=col['colors'],inertia=Q.inertia(C))
    other=Sc[0][0].sign();assert other!=0
    return dict(status='ACCEPT'if other>0 else'REJECT_CONJUGATE_NEGATIVE',colors=col['colors'],physical_last_pivot=S[0][0].enc(),conjugate_last_pivot=Sc[0][0].enc(),physical_inertia=[row['nodes']-1,1,0],conjugate_inertia=[row['nodes'],0,0]if other>0 else[row['nodes']-1,1,0])
def profiles(n):
    for u in range(n):
        for m in ALPHABET:yield [(u,m)]
    for us in combinations(range(n),2):
        for labs in [(3,3),(3,4),(4,3),(3,5),(5,3)]:yield list(zip(us,labs))
    for us in combinations(range(n),3):
        yield [(u,3)for u in us]
        if us==(0,1,2):
            for i in range(3):
                for m in [4,5]:yield [(u,m if j==i else 3)for j,u in enumerate(us)]
def append(row,p):return from_edges(row['nodes']+1,edges(row)+[(i,row['nodes'],m)for i,m in p])
def register(db,row,buckets,registry):
    key=bucket(row)
    for idx in buckets[key]:
        mp=isomorphism(row,registry[idx])
        if mp is not None:return idx,mp,False
    idx=len(registry);registry.append(row);buckets[key].append(idx)
    db.execute('INSERT INTO states(id,row,status) VALUES(?,?,?)',(idx,pack(row),'PENDING'))
    return idx,list(range(row['nodes'])),True
def initialize(db,maxnodes):
    db.executescript('CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT);CREATE TABLE IF NOT EXISTS states(id INTEGER PRIMARY KEY,row TEXT,status TEXT,payload TEXT);CREATE TABLE IF NOT EXISTS seeds(id INTEGER PRIMARY KEY,payload TEXT);')
    if db.execute('SELECT count(*) FROM meta').fetchone()[0]:
        assert json.loads(db.execute("SELECT value FROM meta WHERE key='max_nodes'").fetchone()[0])==maxnodes
        assert db.execute("SELECT value FROM meta WHERE key='source_hash'").fetchone()[0]==digest(f'scripts/{TAG}.py')
        return
    registry=[];buckets=defaultdict(list);cat=json.loads((ROOT/TRI).read_text())
    with db:
        db.execute('INSERT INTO meta VALUES(?,?)',('max_nodes',str(maxnodes)));db.execute('INSERT INTO meta VALUES(?,?)',('source_hash',digest(f'scripts/{TAG}.py')))
        for j,r in enumerate(cat['markings']):
            order=r['black']+r['white'];inv={v:i for i,v in enumerate(order)}
            es=[(inv[i],inv[k],m)for(i,k),m in zip(combinations(range(r['order']),2),r['labels'])if m!=2]
            row=from_edges(r['order'],es);info=corner_info(row);assert info and arithmetic(row)['status']=='ACCEPT'
            idx,mp,new=register(db,row,buckets,registry)
            db.execute('INSERT INTO seeds VALUES(?,?)',(j,pack(dict(marking_index=j,state=idx,original_order=order,map=mp))))
    assert db.execute('SELECT count(*) FROM seeds').fetchone()[0]==76
def export(db,maxnodes,output):
    rows=[];counts=Counter();byrank=Counter();trials=Counter()
    for idx,raw,status,payload in db.execute('SELECT id,row,status,payload FROM states ORDER BY id'):
        row=json.loads(raw);rec=dict(index=idx,row=row,status=status)
        if payload:rec.update(json.loads(payload));trials.update(rec.get('trial_counts',{}))
        rows.append(rec);counts[status]+=1;byrank[row['nodes']]+=1
    seeds=[json.loads(r[0])for r in db.execute('SELECT payload FROM seeds ORDER BY id')]
    src=[f'scripts/{TAG}.py',TRI,'scripts/fresh_sqrt5_triangle_weights_20261015.py','results/fresh_sqrt5_triangle_weights_20261015_review.json','scripts/fresh_sqrt5_quad_filters_20261015.py','results/fresh_sqrt5_quad_filters_20261015_review.json']
    d=dict(status='BOUNDED_EXACT_TRIANGLE_EXTENSION_FRONTIER',sources={p:digest(p)for p in src},maximum_nodes=maxnodes,target_full_triangle_nodes=22,connected_states=len(rows),pending=counts['PENDING'],complete_through_bound=not counts['PENDING'],counts=dict(counts),rank_counts=dict(sorted(byrank.items())),trial_counts=dict(trials),seed_maps=seeds,records=rows,scope='Connected white extensions of all76 minimal Lanner triangle markings, modulo exact black/white colored isomorphism. Orthogonal components remain common spectators for subsequent payment. Boundary states at max_nodes are OPEN for further attachments. No complete dimension21 claim or weight-vector validation.')
    output.write_text(json.dumps(d,indent=2)+'\n');return d
def run(dbpath,output,maxnodes,maxcores):
    db=sqlite3.connect(dbpath);db.execute('PRAGMA journal_mode=WAL');initialize(db,maxnodes)
    registry=[json.loads(r[0])for r in db.execute('SELECT row FROM states ORDER BY id')];buckets=defaultdict(list)
    for i,row in enumerate(registry):buckets[bucket(row)].append(i)
    done=0;start=time.monotonic()
    while done<maxcores:
        hit=db.execute("SELECT id,row FROM states WHERE status='PENDING' ORDER BY id LIMIT 1").fetchone()
        if hit is None:break
        idx,raw=hit;row=json.loads(raw);info=corner_info(row);assert info
        reason=safe_reason(info);payload=dict(**info)
        if reason:status='HEREDITARY_SAFE';payload['reason']=reason
        elif row['nodes']>=maxnodes:status='OPEN_RANK_BOUNDARY'
        else:
            status='EXPANDED';trialrows=[];counts=Counter()
            # One core and every newly registered child are committed atomically.
            with db:
                for p in profiles(row['nodes']):
                    child=append(row,p);ci=corner_info(child)
                    if ci is None:counts['REJECT_NONSPHERICAL_CORNER']+=1;continue
                    aa=arithmetic(child);item=dict(profile=p,**aa);counts[aa['status']]+=1
                    if aa['status']=='ACCEPT':
                        target,mp,new=register(db,child,buckets,registry)
                        item.update(target=target,map=mp,new_state=new)
                    trialrows.append(item)
                payload.update(trials=trialrows,trial_counts=dict(counts));db.execute('UPDATE states SET status=?,payload=? WHERE id=?',(status,pack(payload),idx))
        if status!='EXPANDED':
            with db:db.execute('UPDATE states SET status=?,payload=? WHERE id=?',(status,pack(payload),idx))
        done+=1
        if done%25==0:print('PROGRESS',done,'committed',idx+1,'registry',len(registry),'seconds',round(time.monotonic()-start,2),flush=True)
    result=export(db,maxnodes,output);db.execute('PRAGMA wal_checkpoint(TRUNCATE)');db.close()
    print(result['status'],result['counts'],result['rank_counts'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--max-nodes',type=int,default=9);ap.add_argument('--max-cores',type=int,default=100000);a=ap.parse_args();assert 5<=a.max_nodes<=22
    dbpath=ROOT/'data/enumerations'/f'{TAG}_rank{a.max_nodes}.sqlite';out=ROOT/'results'/f'{TAG}_rank{a.max_nodes}.json';dbpath.parent.mkdir(exist_ok=True,parents=True)
    run(dbpath,out,a.max_nodes,a.max_cores)
