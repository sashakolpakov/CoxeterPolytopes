#!/usr/bin/env python3
"""Bounded resumable Gale mutation guided by literal induced obstructions.

Proposal ranking is heuristic. Saved realizations and exclusions are exact.
"""
import argparse,json,math,os,random,time
from collections import Counter
from pathlib import Path
import numpy as np
from fresh_h9_rank5_gale_20261014_core import ROOT,exact_minors,incidence,realize,simplex_link_factors,sha
from fresh_h9_obstruction_guided_gale_20261014_supports import SupportGates
from fresh_h9_fifteen_facet_20261014_engine import search

TAG='fresh_h9_obstruction_guided_gale_20261014'
OLD='results/fresh_h9_rank5_gale_20261014.json'
SEEDS=[908,638,748,661]


def atomic(path,data):
    p=ROOT/path;t=p.with_suffix('.tmp')
    t.write_text(json.dumps(data,separators=(',',':'))+'\n');t.replace(p)


def tuples(x):
    return tuple(map(tuples,x)) if isinstance(x,list) else x


def run(steps,limit):
    old=json.loads((ROOT/OLD).read_text())
    for p,h in old['sources'].items():assert sha(p)==h,p
    gates=SupportGates()
    files=list(old['sources'])+gates.paths+[OLD,f'scripts/{TAG}.py',f'scripts/{TAG}_supports.py']
    sources={p:sha(p) for p in dict.fromkeys(files)}
    known=set();oldrows={}
    for ref in old['cases']:
        assert sha(ref['path'])==ref['sha256']
        d=json.loads((ROOT/ref['path']).read_text())
        known.add(tuple(d['geometry']['canonical_minimal_nonfaces']))
        if ref['case'] in SEEDS:oldrows[ref['case']]=d
    assert len(known)==1063 and len(oldrows)==4
    options=dict(random_seed=2026101417,seed_cases=SEEDS,steps_per_walk=steps,
                 label_node_limit=limit,coordinate_bound=2048,
                 score='recognized_support_count+6*sum((large_missing_order-5)^2)+4*simplex_link_factor_count')
    folder=ROOT/f'results/{TAG}_cases';folder.mkdir(exist_ok=True)
    output=f'results/{TAG}.json';state=f'results/{TAG}_state.json'
    rng=random.Random(options['random_seed']);counts=Counter();rows=[];seen=set(known)
    walk=index=0;V=None;survivor=None;best=[];start=time.monotonic()
    if (ROOT/state).exists():
        previous=json.loads((ROOT/state).read_text())
        assert previous['sources']==sources and previous['options']==options
        rng.setstate(tuples(previous['random_state']));counts=Counter(previous['counts'])
        walk,index=previous['walk'],previous['index'];survivor=previous['survivor'];best=previous['best']
        V=np.array(previous['current'],dtype=np.int64) if previous['current'] is not None else None
        for p in sorted(folder.glob('case_*.json')):
            d=json.loads(p.read_text());assert d['sources']==sources
            rows.append(d);seen.add(tuple(d['geometry']['canonical_minimal_nonfaces']))

    def save(status):
        refs=[dict(case=r['case'],path=f'results/{TAG}_cases/case_{r["case"]:05d}.json',
                   sha256=sha(f'results/{TAG}_cases/case_{r["case"]:05d}.json'),
                   gate=r['gate']['gate'],verdict=r['gate'].get('search',{}).get('verdict'),
                   recognized_support_count=r['guidance']['recognized_support_count']) for r in rows]
        verdicts=Counter(r['gate'].get('search',{}).get('verdict','EXACT_INHERITED_OBSTRUCTION') for r in rows)
        d=dict(status=status,sources=sources,options=options,counts=dict(counts),cases=refs,
               new_saved_nerves=len(rows),known_keys_not_resaved=1063,first_full_category_survivor=survivor,
               verdicts=dict(verdicts),walk_index=walk,proposal_index=index,best_states=best,
               scope='Bounded heuristic mutation guided by exact literal obstruction supports. Every saved nerve is an exact realizable fifteen-facet candidate with missing faces<=5 and no simplex-boundary link factor. Previously saved1063 keys are not repeated. Inherited gates and complete DFS failures exclude only saved cases; node limits remain unresolved. No exhaustive Gale or fifteen-facet claim.')
        atomic(output,d)
        atomic(state,dict(status=status,sources=sources,options=options,counts=dict(counts),
                          random_state=rng.getstate(),walk=walk,index=index,
                          current=V.tolist() if V is not None else None,survivor=survivor,best=best))
        return d

    def analyze(matrix,origin,record=True):
        data,bad=incidence(matrix,require_small=False)
        if bad:
            if record:counts[bad]+=1
            return None
        mm=data['minimal_nonfaces'];large=[m for m in mm if m.bit_count()>5]
        links=simplex_link_factors(data['facets']) if not large else []
        obstructions=gates.all(mm)
        score=len(obstructions)+6*sum((m.bit_count()-5)**2 for m in large)+4*len(links)
        info=dict(score=score,recognized_support_count=len(obstructions),
                  large_missing_faces=len(large),simplex_link_factor_count=len(links))
        if not record:return data,obstructions,large,links,info
        if large:counts['large_minimal_nonface']+=1
        elif links:counts['simplex_link_factor']+=1
        else:
            counts['eligible_proposals']+=1
            geometry=realize(matrix,data);key=tuple(geometry['canonical_minimal_nonfaces'])
            if key in seen:
                counts['old_known_key' if key in known else 'new_batch_duplicate_key']+=1
            else:
                nonlocal survivor
                seen.add(key);case=len(rows);geometry['case']=case
                if obstructions:
                    gate=obstructions[0]
                else:
                    print('LABEL_START',case,origin,flush=True)
                    answer=search(geometry,limit)
                    gate=dict(gate='complete_label_search',search=answer)
                    print('LABEL_END',case,answer['verdict'],answer['nodes'],flush=True)
                    if answer['verdict']=='FULL_CATEGORY_METRIC_REQUIRED':survivor=case
                counts[gate['gate']]+=1
                row=dict(case=case,sources=sources,origin=origin,geometry=geometry,gate=gate,
                         guidance={**info,'obstructions':obstructions},label_node_limit=limit)
                atomic(f'results/{TAG}_cases/case_{case:05d}.json',row);rows.append(row)
                print('NEW',case,'walk',origin['walk'],'step',origin['index'],'score',score,
                      'gate',gate['gate'],'seconds',round(time.monotonic()-start,1),flush=True)
        if not best or score < best[-1]['score']:
            best.append(dict(walk=walk,index=index,seed=SEEDS[walk],**info,gale=matrix.tolist(),
                             minimal_nonfaces=mm,recognized_obstruction_supports=[r['support'] for r in obstructions]))
            print('BEST',info,'walk',walk,'step',index,flush=True)
        return data,obstructions,large,links,info

    def proposal(matrix,current):
        data,obs,large,links,info=current
        if large and (not obs or rng.random()<.25):
            support=rng.choice(large);target='large_missing_face'
        elif obs:
            support=rng.choice(obs)['support'];target='recognized_obstruction'
        elif links:
            link=rng.choice(links);support=link['factor']|(1<<link['vertex']);target='link_factor'
        else:
            support=32767;target='unconstrained'
        inside=[i for i in range(15) if support>>i&1]
        outside=[i for i in range(15) if not support>>i&1]
        i=rng.choice(inside)
        j=rng.choice(outside if outside and rng.random()<.75 else [j for j in range(15) if j!=i])
        col=rng.randrange(5)
        unit=matrix.copy();unit[i,col]+=1;unit[j,col]-=1
        direction=1
        if np.max(np.abs(unit))>2048:
            unit=matrix.copy();unit[i,col]-=1;unit[j,col]+=1;direction=-1
        if np.max(np.abs(unit))>2048:return None
        # Both row changes have the same column direction: this is a rank-one
        # matrix update, so every five-minor is exactly affine in delta.
        derivative=(exact_minors(unit)-data['det'])*direction
        candidates=set()
        for a,b in zip(data['det'],derivative):
            a,b=int(a),int(b)
            if b:
                s=-1 if a*b>0 else 1
                delta=s*(abs(a)//abs(b)+1)
                if abs(int(matrix[i,col])+delta)<=2048 and abs(int(matrix[j,col])-delta)<=2048:
                    candidates.add(delta)
        ordered=sorted(candidates,key=lambda x:(abs(x),x))[:16]
        if not ordered:return None
        if rng.random()<.1:
            delta=rng.choice([-1,1])*rng.choice([32,64,128,256,512]);kind='targeted_long_jump'
        else:
            rank=rng.choice([0,0,0,1,1,2,3,4,7,11,15]);delta=ordered[min(rank,len(ordered)-1)]
            kind='exact_five_minor_boundary_crossing'
        trial=matrix.copy();trial[i,col]+=delta;trial[j,col]-=delta
        if np.max(np.abs(trial))>2048:return None
        return trial,dict(target=target,target_support=support,rows=[i,j],column=col,
                          delta=delta,mutation=kind)

    print('PID',os.getpid(),'resume',walk,index,'options',options,flush=True)
    if survivor is not None:return save('STOPPED_AT_FULL_CATEGORY_SURVIVOR')
    while walk<len(SEEDS):
        if V is None:V=np.array(oldrows[SEEDS[walk]]['geometry']['gale'],dtype=np.int64);index=0
        current=analyze(V,dict(walk=walk,index=index),record=False);assert current
        counts['walks_started']+=int(index==0)
        if index==0:print('START_WALK',walk,'seed',SEEDS[walk],current[-1],flush=True)
        while index<steps:
            p=proposal(V,current);counts['proposals']+=1
            if p is None:
                counts['no_bounded_mutation']+=1
            else:
                trial,origin=p;counts[origin['mutation']]+=1
                result=analyze(trial,dict(walk=walk,index=index,**origin))
                if result is not None:
                    temperature=1.5+3*(1-(index%750)/750)
                    difference=result[-1]['score']-current[-1]['score']
                    if difference<=0 or rng.random()<math.exp(-difference/temperature):
                        V=trial;current=result;counts['accepted_moves']+=1
            index+=1
            if index%100==0 or survivor is not None:
                save('RUNNING_BOUNDED_OBSTRUCTION_GUIDED_GALE_WALK')
                print('PROGRESS',walk,index,'cases',len(rows),'current',current[-1],flush=True)
            if survivor is not None:return save('STOPPED_AT_FULL_CATEGORY_SURVIVOR')
        walk+=1;index=0;V=None;save('RUNNING_BOUNDED_OBSTRUCTION_GUIDED_GALE_WALK')
    return save('FINISHED_BOUNDED_OBSTRUCTION_GUIDED_GALE_WALK')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--steps',type=int,default=3000);p.add_argument('--nodes',type=int,default=10000)
    a=p.parse_args();out=run(a.steps,a.nodes)
    print(out['status'],out['new_saved_nerves'],out['verdicts'],flush=True)
