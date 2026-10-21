#!/usr/bin/env python3
"""Physical retained outputs and exact count for the guarded D16 unit zero."""
import argparse
import fcntl
from fractions import Fraction as F
import json
from pathlib import Path
import signal
import fresh_d16_node3_guard_complete_finite_20261018 as J

N=J.N;S=N.S;V=N.V;M=N.M;U=N.U;ROOT=N.ROOT;P=S.P;B=S.B
TAG='fresh_d16_node3_rank16_large_d_guard_transfer_20261018'
value=J.value


def role(g,w,marks):
    ids=sorted(set(w)|set(marks))
    if S.Q.inertia(S.A.sub(g,ids))!=[len(ids),0,0]:return[]
    b=next((b for b in S.D.finite_blocks(g,ids)if b['type']=='D16'and set(marks)<=set(b['nodes'])),None)
    return[]if b is None else S.T.mark_pairs(g,b['nodes'],marks,16,'D')


def audit():
    later=[];sources=set();count=0
    for stem,keys in B.FINITE.items():
        doc=V.load(stem);sources.add(stem);rs=[]
        for key in keys:
            for path,p,g,w,bs in B.corners(doc[key],key):
                im=role(S.A.matrix(g),w,p['marks']);assert(3,13)not in im
                rs.append(dict(path=path,marks=p['marks'],D16_images=im));count+=1
        later.append(dict(source=stem,rows=rs))
    assert count==564
    docs={s:V.load(s)for s in V.O.LITERAL if s!=S.A.OVER};sources.update(docs)
    lp=V.load(S.TAG)['local']['selected_and_later_seals'];sources.add(S.TAG)
    g=next(p['gram']for p in lp if p['index']==1013)
    docs[S.A.OVER]={'isolation':{'cores':{'(25, 25, ((9, 9),))':{'gram':g}}}}
    literal=[]
    for stem,index,g,w,cs in S.T.literal_inputs(docs):
        rows=[]
        for marks,*_ in cs:
            im=role(g,w,marks);assert(3,13)not in im;rows.append(dict(marks=marks,D16_images=im))
        literal.append(dict(source=stem,index=index,corners=rows))
    assert len(literal)==216
    newer=[]
    collections=[(P.R.TAG,('local_proof','whole_proofs')),(P.H.TAG,('local_proof','whole_proofs')),
        (P.S.TAG,('local_proof','whole_proofs')),(P.N.TAG,('local_proof','proofs')),
        (P.R.R.TAG,('local_seals',)),(P.H.Z.U.TAG,('local_proof',)),
        (P.TAG,('local_proof','old_ideal17_seals')),(P.TAG,('local_proof','A16_seals'))]
    for stem,keys in collections:
        xs=V.load(stem);sources.add(stem)
        for k in keys:xs=xs[k]
        if isinstance(xs,dict):xs=[xs]
        for p in xs:
            pp=p.get('old_complete_seal',p);g=pp.get('gram',pp.get('old_seal',{}).get('gram'));w=pp.get('white',pp.get('old_seal',{}).get('white'))
            if g is None:continue
            rows=[]
            for c in p.get('corners',[]):
                im=role(S.A.matrix(g),w,c['marks']);assert(3,13)not in im;rows.append(dict(marks=c['marks'],D16_images=im))
            newer.append(dict(source=stem,index=p.get('index',p.get('name')),rows=rows))
    for stem,proofs in((S.TAG,lp),(V.TAG,V.load(V.TAG)['local_proof']['proofs'])):
        sources.add(stem)
        for p in proofs:
            rows=[]
            for c in p['corners']:
                im=role(S.A.matrix(p['gram']),p['white'],c['marks']);assert(3,13)not in im;rows.append(dict(marks=c['marks'],D16_images=im))
            newer.append(dict(source=stem,index=p['index'],rows=rows))
    chosen=[]
    for stem in(B.PART,B.REST):
        doc=V.load(stem);sources.add(stem);rs=[(path,c)for path,c in P.W.walk(doc)if F(c['payment'])>0]
        assert len(rs)==(280 if stem==B.PART else 52)
        assert all(not any(b['type']=='D16'for b in c['roles'])for _,c in rs)
        chosen.append(dict(source=stem,positive_outputs=len(rs),D16_outputs=0))
    return dict(later_packets=later,later_corners=count,literal=literal,literal_configurations=216,
        newer_seals=newer,chosen_D12=chosen,no_target_consumer=True),sources


def run():
    scout=V.load(N.TAG);assert V.O.enc(N.geometry())==scout['geometry'];assert V.O.enc(N.upper())==scout['upper']
    retained,sources=audit();finite=V.load(J.TAG)
    hstem='fresh_d16_a1_long_leaf_d17_zero_ambient58_frontier_20261018';h=V.load(hstem);sources.add(hstem)
    weights={tuple(r['pair']):F(r['weight'])for r in h['local_proof']['D16_ordinary_weighted_pairs']}
    assert weights[3,13]==1 and sum(weights.values())==F(217,2);weights[3,13]=F()
    assert sum(weights.values())==F(215,2)
    boundaries=[]
    for index in N.FULL+N.PROPER:
        c,core,g,w,left,right=N.get(index);boundaries.append(dict(index=index,left=left,right=right,selected_marks=c['preferred_corner']))
    paths={Path(__file__),Path(J.__file__),Path(N.__file__),Path(S.__file__),S.OVER,
        ROOT/f'results/{N.TAG}.json',ROOT/f'results/{J.TAG}.json'}
    paths.update(ROOT/f'results/{s}.json'for s in sources)
    return V.O.enc(dict(status='PASS_D16_NODE3_RANK16_LARGE_D_GUARD_TRANSFER',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(paths)},
        complete_local_geometry=scout['geometry'],native_boundaries=boundaries,current_retention=retained,
        ordinary_D16_cost='217/2',guarded_D16_cost='215/2',guarded_D16_cost_with_A1='217/2',
        final_guarded_D16_weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(weights.items())],
        changed_pair=[3,13],debit_guard='separate whole A16 orD16 and separate whole D_s>=18',
        coupon=dict(small='A2',large='A13',large_nodes=[3,9],whole_A1=True,increment='1',component_increment='4',same_external_witnesses=True),
        finite=dict(exact_global_maximum=finite['new_exact_global_formal_maximum'],attaining_forest=finite['attaining_forest'],
            increasing_family_upper=scout['upper']['increasing_family_upper']),
        ideal_weights_unchanged=True,triangle_minimum='2',quadrilateral_minimum='1',completed_dimension_bound=57,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);M.R.resources();U.memory_guard()
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();M.R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['exact_global_maximum'],out['current_retention']['later_corners'],flush=True)
