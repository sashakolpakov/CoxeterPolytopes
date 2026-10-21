#!/usr/bin/env python3
"""Separate saved-case Bareiss/face replay plus all literal guidance witnesses."""
import argparse,json,os,time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from fresh_h9_rank5_gale_20261014_independent import replay

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_obstruction_guided_gale_20261014'
OLD='results/fresh_h9_rank5_gale_20261014.json'


@lru_cache(None)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def load(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d.get('sources',{}).items():assert sha(f)==h,f
    return d


@lru_cache(None)
def pattern(path,index):
    d=load(path)
    row=d['seven_node_obstruction'] if index=='seven_node_relation' else next(r for r in d['patterns'] if r['pattern']==index)
    assert row['verdict']=='EXACT_LOCAL_OBSTRUCTION'
    return row


def one(path,expected):
    r=replay(path,expected)
    d=load(path);mm=d['geometry']['minimal_nonfaces'];g=d['guidance'];supports=[]
    assert not g['large_missing_faces'] and not g['simplex_link_factor_count']
    for gate in g['obstructions']:
        u=gate['support'];inside=[m for m in mm if m&u==m]
        if gate['gate']=='two_lanner_blocks':
            a,b=gate['blocks']
            assert a|b==u and not a&b and set(inside)=={a,b}
            assert sorted((a.bit_count(),b.bit_count())) in ([3,4],[3,5],[4,4],[4,5],[5,5])
        else:
            assert gate['gate']=='literal_induced_pattern'
            p=pattern(gate['source'],gate['pattern']);order=gate['permutation']
            assert len(order)==len(set(order))==p['order'] and sum(1<<i for i in order)==u
            inverse={v:i for i,v in enumerate(order)}
            reduced=sorted(sum(1<<inverse[v] for v in order if m>>v&1) for m in inside)
            assert reduced==p['minimal_nonfaces']
        supports.append(u)
    assert len(supports)==len(set(supports))==g['recognized_support_count']==g['score']
    if supports:assert d['gate']==g['obstructions'][0]
    else:assert d['gate']['gate']=='complete_label_search'
    return {**r,'all_guidance_witnesses_rechecked':len(supports)}


def atomic(p,d):
    t=p.with_suffix('.tmp');t.write_text(json.dumps(d,separators=(',',':'))+'\n');t.replace(p)


def run(jobs,follow):
    primary=f'results/{TAG}.json';folder=ROOT/f'results/{TAG}_audit_cases';folder.mkdir(exist_ok=True)
    source_files=[f'scripts/{TAG}_independent.py','scripts/fresh_h9_rank5_gale_20261014_independent.py',
                  'scripts/fresh_h9_fifteen_facet_20261014_engine.py',
                  'scripts/fresh_h9_bugaenko_closure_20261014_large_first_audit.py',
                  'scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py',OLD]
    sources={p:sha(p) for p in source_files}
    old=load(OLD);known=set()
    for ref in old['cases']:
        assert sha(ref['path'])==ref['sha256']
        known.add(tuple(load(ref['path'])['geometry']['canonical_minimal_nonfaces']))
    records={}
    for p in folder.glob('case_*.json'):
        d=load(str(p.relative_to(ROOT)));assert d['sources']==sources
        records[d['result']['case']]=d['result']
    pending={};start=time.monotonic();initial=load(primary)
    def save(final=False):
        d=json.loads((ROOT/primary).read_text());refs={r['case']:r for r in d['cases']}
        for i,r in records.items():assert r['input_sha256']==refs[i]['sha256']
        paths=dict(sources)
        if final:paths[primary]=sha(primary)
        keys=set()
        for i in records:
            key=tuple(load(refs[i]['path'])['geometry']['canonical_minimal_nonfaces'])
            assert key not in known and key not in keys;keys.add(key)
        out=dict(status='PASS_ALL_SAVED_GUIDED_GALE_CASES' if final else 'RUNNING_SEPARATE_GUIDED_GALE_REPLAY',
                 sources=paths,cases_checked=len(records),primary_saved_cases=len(refs),
                 prior_keys_avoided=len(known),five_minors_rechecked=3003*len(records),
                 six_circuits_rechecked=5005*len(records),vertex_links_rechecked=15*len(records),
                 actual_vertices_rechecked=sum(r['vertices'] for r in records.values()),
                 guidance_witnesses_rechecked=sum(r['all_guidance_witnesses_rechecked'] for r in records.values()),
                 records=[records[i] for i in sorted(records)],
                 scope='Reuses the separate integer Bareiss/full-face/link-ideal checker, not the proposal geometry engine. Every saved induced gate and every archived guidance witness is literally rechecked. Old obstruction-library arithmetic remains inherited. Any complete label search is replayed with separate finite-tree and modular arithmetic; limits stay unresolved. Heuristic walk coverage is not a classification claim.')
        atomic(ROOT/f'results/{TAG}_independent.json',out);return out
    print('PID',os.getpid(),'resumed',len(records),flush=True)
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        while True:
            d=json.loads((ROOT/primary).read_text());assert d['sources']==initial['sources']
            terminal=not d['status'].startswith('RUNNING');waiting=set(pending.values())
            todo=[r for r in d['cases'] if r['case'] not in records and r['case'] not in waiting]
            for r in todo[:jobs-len(pending)]:pending[pool.submit(one,r['path'],r['sha256'])]=r['case']
            if not pending:
                if terminal or not follow:break
                time.sleep(2);continue
            done,_=wait(pending,timeout=3,return_when=FIRST_COMPLETED)
            for f in done:
                i=pending.pop(f);r=f.result();records[i]=r
                atomic(folder/f'case_{i:05d}.json',dict(sources=sources,result=r))
                if len(records)%20==0:print('checked',len(records),'seconds',round(time.monotonic()-start,1),flush=True)
            if done:save()
    return save(final=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--jobs',type=int,default=2);p.add_argument('--follow',action='store_true')
    a=p.parse_args();d=run(a.jobs,a.follow);print(d['status'],d['cases_checked'],flush=True)
