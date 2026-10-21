#!/usr/bin/env python3
"""Bounded positive-circuit proposals on targeted ten-vertex flip carriers.

Force five Gale rows to have a positive dependence, then perturb to either
side. Only balanced 5/5 complement partitions are proposed. This avoids most
minor crossings that cannot alter the actual facet triangulation.
"""
import argparse,json,os,random,time
from collections import Counter
from pathlib import Path
import numpy as np
from fresh_h9_rank5_gale_20261014_core import ROOT,incidence,realize,simplex_link_factors,sha
from fresh_h9_rank5_gale_20261014_independent import det
from fresh_h9_obstruction_guided_gale_20261014_supports import SupportGates
from fresh_h9_fifteen_facet_20261014_engine import search

TAG='fresh_h9_obstruction_guided_gale_20261014_circuits'
OLD='results/fresh_h9_rank5_gale_20261014.json'
GUIDED='results/fresh_h9_obstruction_guided_gale_20261014.json'


def atomic(path,d):
    p=ROOT/path;t=p.with_suffix('.tmp');t.write_text(json.dumps(d,separators=(',',':'))+'\n');t.replace(p)


def tuples(x):return tuple(map(tuples,x)) if isinstance(x,list) else x


def run(attempts,limit):
    gates=SupportGates();known=set();docs={};refs={};previous={}
    for p in [OLD,GUIDED]:
        d=json.loads((ROOT/p).read_text());previous[p]=d
        for f,h in d['sources'].items():assert sha(f)==h,f
        for ref in d['cases']:
            assert sha(ref['path'])==ref['sha256']
            row=json.loads((ROOT/ref['path']).read_text());known.add(tuple(row['geometry']['canonical_minimal_nonfaces']))
            docs[ref['path']]=row;refs[ref['path']]=ref
    eligible=sorted(previous[GUIDED]['cases'],key=lambda r:(r['recognized_support_count'],r['case']))
    seeds=[];features=set()
    for ref in eligible:
        g=docs[ref['path']]['geometry'];feature=tuple(sorted(Counter(m.bit_count() for m in g['minimal_nonfaces']).items()))
        if feature not in features:
            seeds.append(ref['path']);features.add(feature)
        if len(seeds)==8:break
    assert len(seeds)==8 and len(known)==1135
    files=list(previous[GUIDED]['sources'])+[OLD,GUIDED,f'scripts/{TAG}.py',
            'scripts/fresh_h9_rank5_gale_20261014_independent.py']+seeds
    sources={p:sha(p) for p in dict.fromkeys(files)}
    options=dict(random_seed=2026101429,boundary_attempts=attempts,seed_paths=seeds,
                 label_node_limit=limit,complement_partition=[5,5],perturbations=[-1,1],coordinate_bound=2048)
    outpath=f'results/{TAG}.json';statepath=f'results/{TAG}_state.json'
    folder=ROOT/f'results/{TAG}_cases';folder.mkdir(exist_ok=True)
    rng=random.Random(options['random_seed']);counts=Counter();rows=[];seen=set(known);index=0;survivor=None
    if (ROOT/statepath).exists():
        d=json.loads((ROOT/statepath).read_text());assert d['sources']==sources and d['options']==options
        rng.setstate(tuples(d['random_state']));counts=Counter(d['counts']);index=d['index'];survivor=d['survivor']
        for p in sorted(folder.glob('case_*.json')):
            d=json.loads(p.read_text());assert d['sources']==sources
            rows.append(d);seen.add(tuple(d['geometry']['canonical_minimal_nonfaces']))
    start=time.monotonic()

    def save(status):
        refrows=[dict(case=r['case'],path=f'results/{TAG}_cases/case_{r["case"]:05d}.json',
                      sha256=sha(f'results/{TAG}_cases/case_{r["case"]:05d}.json'),gate=r['gate']['gate'],
                      verdict=r['gate'].get('search',{}).get('verdict'),recognized_support_count=r['guidance']['recognized_support_count']) for r in rows]
        out=dict(status=status,sources=sources,options=options,counts=dict(counts),cases=refrows,
                 new_saved_nerves=len(rows),prior_keys_avoided=len(known),completed_boundary_attempts=index,
                 first_full_category_survivor=survivor,
                 verdicts=dict(Counter(r['gate'].get('search',{}).get('verdict','EXACT_INHERITED_OBSTRUCTION') for r in rows)),
                 scope='Bounded balanced positive-five-circuit proposals, both integer perturbation sides. Every retained nerve is realizable, has missing faces<=5 and lacks a simplex-boundary vertex-link factor. No claim that the boundary mechanism samples all oriented-matroid mutations. Inherited exact obstructions apply literally; node limits remain unresolved.')
        atomic(outpath,out);atomic(statepath,dict(status=status,sources=sources,options=options,counts=dict(counts),
                                                index=index,random_state=rng.getstate(),survivor=survivor))
        return out

    def record(matrix,origin):
        nonlocal survivor
        counts['perturbed_proposals']+=1
        if np.max(np.abs(matrix))>2048:counts['coordinate_bound']+=1;return
        data,bad=incidence(matrix)
        if bad:counts[bad]+=1;return
        links=simplex_link_factors(data['facets'])
        if links:counts['simplex_link_factor']+=1;return
        counts['eligible_proposals']+=1
        geometry=realize(matrix,data);key=tuple(geometry['canonical_minimal_nonfaces'])
        if key in seen:
            counts['old_known_key' if key in known else 'batch_duplicate_key']+=1;return
        seen.add(key);case=len(rows);geometry['case']=case;obstructions=gates.all(data['minimal_nonfaces'])
        if obstructions:gate=obstructions[0]
        else:
            print('LABEL_START',case,index,flush=True);answer=search(geometry,limit)
            print('LABEL_END',case,answer['verdict'],answer['nodes'],flush=True)
            gate=dict(gate='complete_label_search',search=answer)
            if answer['verdict']=='FULL_CATEGORY_METRIC_REQUIRED':survivor=case
        counts[gate['gate']]+=1
        row=dict(case=case,sources=sources,origin=origin,geometry=geometry,gate=gate,label_node_limit=limit,
                 guidance=dict(score=len(obstructions),recognized_support_count=len(obstructions),
                               large_missing_faces=0,simplex_link_factor_count=0,obstructions=obstructions))
        atomic(f'results/{TAG}_cases/case_{case:05d}.json',row);rows.append(row)
        print('NEW',case,'attempt',index,'supports',len(obstructions),'gate',gate['gate'],
              'vertices',len(data['facets']),'seconds',round(time.monotonic()-start,1),flush=True)

    print('PID',os.getpid(),'resume',index,'seeds',seeds,flush=True)
    if survivor is not None:return save('STOPPED_AT_FULL_CATEGORY_SURVIVOR')
    while index<attempts:
        seed=rng.choice(seeds);d=docs[seed];V=np.array(d['geometry']['gale'],dtype=np.int64)
        obstruction=rng.choice(d['guidance']['obstructions']);u=obstruction['support']
        target=rng.choice([m for m in d['geometry']['minimal_nonfaces'] if m&u==m])
        target_rows=[i for i in range(15) if target>>i&1]
        B=sorted(rng.sample([i for i in range(15) if not target>>i&1],5))
        X=[i for i in range(15) if i not in B]
        shift=V[B].sum(axis=0);possibilities=[]
        for i in B:
            four=[[int(x) for x in V[j]] for j in B if j!=i]
            normal=[(-1)**k*det([r[:k]+r[k+1:] for r in four]) for k in range(5)]
            assert any(normal)
            assert all(sum(a*b for a,b in zip(r,normal))==0 for r in four)
            products={k:sum(int(a)*b for a,b in zip(V[k],normal)) for k in X}
            adjustment=sum(int(a)*b for a,b in zip(shift,normal))
            for j in target_rows:
                signs={k:products[k]+(adjustment if k==j else 0) for k in X}
                if any(v==0 for v in signs.values()):continue
                positive=[k for k,v in signs.items() if v>0]
                negative=[k for k,v in signs.items() if v<0]
                if len(positive)!=5 or len(negative)!=5:continue
                if not (set(target_rows)<=set(positive) or set(target_rows)<=set(negative)):continue
                W=V.copy();W[i]-=shift;W[j]+=shift
                if np.max(np.abs(W))>2047:continue
                assert np.all(W[B].sum(axis=0)==0)
                assert all(sum(int(a)*b for a,b in zip(W[k],normal))==0 for k in B)
                possibilities.append((i,j,normal,positive,negative,W))
        if not possibilities:
            counts['no_balanced_target_carrier']+=1
        else:
            counts['balanced_positive_circuit_boundaries']+=1
            i,j,normal,positive,negative,W=rng.choice(possibilities)
            col=max(range(5),key=lambda k:abs(normal[k]))
            for side in [-1,1]:
                trial=W.copy();trial[i,col]+=side;trial[j,col]-=side
                record(trial,dict(attempt=index,seed=seed,target_obstruction_support=u,target_missing_face=target,
                                  positive_circuit=B,changed_rows=[i,j],hyperplane=normal,
                                  positive_complement=positive,negative_complement=negative,
                                  perturbation_column=col,perturbation_side=side))
                if survivor is not None:break
        index+=1
        if index%100==0 or survivor is not None:
            save('RUNNING_BOUNDED_BALANCED_POSITIVE_CIRCUIT_PROPOSALS')
            print('PROGRESS',index,'new',len(rows),'counts',dict(counts),flush=True)
        if survivor is not None:return save('STOPPED_AT_FULL_CATEGORY_SURVIVOR')
    return save('FINISHED_BOUNDED_BALANCED_POSITIVE_CIRCUIT_PROPOSALS')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--attempts',type=int,default=2000);p.add_argument('--nodes',type=int,default=10000)
    a=p.parse_args();d=run(a.attempts,a.nodes)
    print(d['status'],d['new_saved_nerves'],d['verdicts'],flush=True)
