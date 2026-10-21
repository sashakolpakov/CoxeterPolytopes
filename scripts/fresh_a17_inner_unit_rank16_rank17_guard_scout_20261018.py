#!/usr/bin/env python3
"""Separate A17(1,15)/A15 endpoint guarded-exchange scout."""
import argparse
from collections import Counter
from fractions import Fraction as F
import fcntl
from itertools import combinations_with_replacement
import json
from pathlib import Path
import signal
import fresh_d16_rank17_witness_extension_scout_20261018 as R
import fresh_d16_node3_rank16_large_d_guard_transfer_20261018 as L

ROOT=R.ROOT;C=R.C;V=R.V;U=R.U;S=L.S;B=L.B
TAG='fresh_a17_inner_unit_rank16_rank17_guard_scout_20261018'
RANK15='fresh_rank15_sparse_payment_exchange_20261018'


def value(fs):
    co=Counter(fs);r16=co['A',16]+co['D',16];r17=co['A',17]+co['D',17]
    if not V.has_guard(fs):return R.value(fs)
    return R.value(fs)-(64*co['A',17]if r16 or r17>=2 else 0)+(64*co['A',15]if r16 or r17 else 0)


def role(g,w,marks):
    ids=sorted(set(w)|set(marks))
    if S.Q.inertia(S.A.sub(g,ids))!=[len(ids),0,0]:return[]
    b=next((b for b in S.D.finite_blocks(g,ids)if b['type']=='A17'and set(marks)<=set(b['nodes'])),None)
    return[]if b is None else S.T.mark_pairs(g,b['nodes'],marks,17,'A')


def geometry():
    p=S.native(1638);g=p['gram'];assert p['seal']['method']=='full_radical'
    old=next(p for p in V.load(RANK15)['selected_proofs']if p['index']==1638)
    assert old['gram']==S.V.O.enc(g)and p['seal']['allowed_contacts']==old['sealing']['allowed_contacts']
    cs=[]
    for c in p['corners']:
        marks=c['marks'];a=next((b for b in c['components']if set(marks)<=set(b['nodes'])),None)
        im=[]if a is None else S.T.mark_pairs(g,a['nodes'],marks,len(a['nodes']),a['type'][0])
        restore=a is not None and a['type']=='A15'and im==[(0,14)]
        debit=a is not None and a['type']=='A17'and im==[(1,15)]
        prev=next(z for z in old['corners']if z['marks']==marks);pay=F(prev['payment'])+int(restore)-int(debit)
        assert pay>=0;cs.append(dict(**c,canonical_images=im,old_payment=prev['payment'],restore=restore,debit=debit,payment=pay))
    assert sum(c['restore']for c in cs)==sum(c['debit']for c in cs)==1 and sum(c['payment']for c in cs)==1
    p['corners']=cs
    sup=V.load('fresh_a17_alternate_corner_pilot_20261018')['selected_support']
    assert sup['complete_eligible_patterns']==140 and sup['complete_maps']==48
    assert all([1,15]not in r['marked_images']for r in sup['rows'])
    base=V.load('fresh_finite_d17_long_pair_support_20261017');rec=next(r for r in base['records']if r['kind']=='A')
    assert [1,15]in rec['removed_pairs']and [1,15]not in rec['conservatively_required_pairs']
    return dict(consumer=p,complete_selected_patterns=140,complete_selected_maps=48,
        selected_domain_avoids_changed_pair=True,original_required_set_avoids_changed_pair=True)


def retention():
    hits=[];sources=set();counts=Counter()
    def one(stem,path,g,w,c):
        im=role(S.A.matrix(g),w,c['marks']);counts[stem]+=1
        if(1,15)in im:
            pay={k:c[k]for k in('payment','used_new_payment','retained_payment_lower_bound','hypothetical_payment','reason','selected')if k in c}
            hits.append(dict(source=stem,path=path,marks=c['marks'],images=im,saved=pay))
    for stem,keys in B.FINITE.items():
        doc=V.load(stem);sources.add(stem)
        for key in keys:
            for path,c,g,w,bs in B.corners(doc[key],key):one(stem,path,g,w,c)
    assert sum(counts.values())==564
    doc=V.load(RANK15);sources.add(RANK15)
    for path,c,g,w,bs in B.corners(doc['selected_proofs'],'selected_proofs'):one(RANK15,path,g,w,c)
    docs={s:V.load(s)for s in V.O.LITERAL if s!=S.A.OVER};sources.update(docs)
    lp=V.load(S.TAG)['local']['selected_and_later_seals'];sources.add(S.TAG)
    g=next(p['gram']for p in lp if p['index']==1013)
    docs[S.A.OVER]={'isolation':{'cores':{'(25, 25, ((9, 9),))':{'gram':g}}}}
    n=0
    for stem,index,g,w,cs in S.T.literal_inputs(docs):
        for j,(marks,*_)in enumerate(cs):one(stem,f'literal/{index}/{j}',g,w,dict(marks=marks))
        n+=1
    assert n==216
    P=S.P
    collections=[(P.R.TAG,('local_proof','whole_proofs')),(P.H.TAG,('local_proof','whole_proofs')),
        (P.S.TAG,('local_proof','whole_proofs')),(P.N.TAG,('local_proof','proofs')),
        (P.R.R.TAG,('local_seals',)),(P.H.Z.U.TAG,('local_proof',)),
        (P.TAG,('local_proof','old_ideal17_seals')),(P.TAG,('local_proof','A16_seals'))]
    for stem,keys in collections:
        xs=V.load(stem);sources.add(stem)
        for k in keys:xs=xs[k]
        if isinstance(xs,dict):xs=[xs]
        for i,p in enumerate(xs):
            pp=p.get('old_complete_seal',p);g=pp.get('gram',pp.get('old_seal',{}).get('gram'));w=pp.get('white',pp.get('old_seal',{}).get('white'))
            if g is None:continue
            for j,c in enumerate(p.get('corners',[])):one(stem,f'newer/{i}/{j}',g,w,c)
    for stem,proofs in((S.TAG,lp),(V.TAG,V.load(V.TAG)['local_proof']['proofs'])):
        sources.add(stem)
        for p in proofs:
            for j,c in enumerate(p['corners']):one(stem,f"newer/{p['index']}/{j}",p['gram'],p['white'],c)
    for group in V.load(L.N.TAG)['geometry'].values():
        if not isinstance(group,list):continue
        for p in group:
            if 'gram'not in p:continue
            for j,c in enumerate(p.get('corners',[])):one(L.N.TAG,f"D16/{p['index']}/{j}",p['gram'],p['white'],c)
            for mp in p.get('maps',[]):
                for j,c in enumerate(mp.get('corners',[])):one(L.N.TAG,f"D16completed/{p['index']}/{j}",mp['completed_gram'],mp['completed_white'],c)
    sources.add(L.N.TAG)
    chosen=[]
    for stem in(B.PART,B.REST):
        doc=V.load(stem);sources.add(stem);rs=[(path,c)for path,c in P.W.walk(doc)if F(c['payment'])>0]
        assert all(not any(b['type']=='A17'for b in c['roles'])for _,c in rs)
        chosen.append(dict(source=stem,positive_outputs=len(rs),A17_outputs=0))
    assert len(hits)==1 and hits[0]['source']==RANK15 and F(hits[0]['saved']['payment'])==1
    proof_number=int(hits[0]['path'].split('/')[1])
    assert V.load(RANK15)['selected_proofs'][proof_number]['index']==1638
    return dict(corner_counts=dict(counts),target_hits=hits,literal_configurations=n,chosen_D12=chosen),sources


def positive_family():
    small=[t for t in C.TYPES if t[1]<=8 and t!=('A',1)]
    assert all(C.actual_cross(a,b)==0 for a in small for b in small)
    rows=[];checks=preds=0
    for witness in(('A',16),('D',16),('A',17),('D',17)):
        for s in range(18,43-witness[1]):
            fixed=[('A',15),witness,('D',s)];left=42-witness[1]-s;tabs=[];opts=[]
            for a1 in range(left+1):
                anchor=fixed+[('A',1)]*a1;cap=left-a1;bc=value(anchor);ts=[t for t in small if t[1]<=cap];ws=[value(anchor+[t])-bc for t in ts]
                tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b;tabs.append(dict(a1=a1,table=tab))
                for i,a in enumerate(ts):
                    for j,b in enumerate(ts[i:],i):
                        if a[1]+b[1]<=cap:assert value(anchor+[a,b])==bc+ws[i]+ws[j]
                if tab['values'][cap]is not None:opts.append((bc+tab['values'][cap],a1))
            z,a1=max(opts);forest=fixed+[('A',1)]*a1;tab=tabs[a1]['table'];cap=left-a1
            while cap:t=tab['types'][tab['last'][cap]];forest.append(t);cap-=t[1]
            assert value(forest)==z and sum(t[1]for t in forest)==57
            rows.append(dict(witness=witness,reserved_D=s,tables=tabs,maximum=F(z,64),attaining_forest=forest))
    return dict(rows=rows,Bellman_inequalities=checks,attaining_predecessors=preds,complete_positive_family_upper=max(r['maximum']for r in rows))


def finite():
    old=V.load(R.TAG)['finite']['original_compact_class_maxima'];high=[r for r in old if F(r['maximum'])>371]
    new=[r for r in high if r['counts'][8]]
    assert{(r['counts'][7],r['counts'][8])for r in new}=={(0,1),(1,1),(0,2)}
    assert len(high)==6 and all(sum(r['counts'][:7])==r['k5']==r['a12']==0 for r in high)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS};res=[t for t in C.TYPES if C.representative(t)in reps and t!=('A',1)]
    assert ('A',15)not in res and ('A',17)not in res and ('D',17)in res
    assert C.representative(('A',17))==('A',17)and C.representative(('D',17))==('D',17)
    assert all(C.actual_cross(a,b)==0 for a in res+[('A',1)]for b in res+[('A',1)])
    fibers=[('A',16),('B',16),('D',16)];out=[];checks=preds=options=0
    for k16,k17 in((0,1),(1,1),(0,2)):
        rows=[]
        for chosen in combinations_with_replacement(fibers,k16):
            fixed=[('A',17)]*k17+list(chosen);rem=57-17*k17-16*k16;tabs=[];opts=[]
            for a1 in range(rem+1):
                for guard in(False,True):
                    if guard and rem-a1<18:continue
                    anchor=fixed+[('A',1)]*a1+([('D',18)]if guard else[]);cap=rem-a1-(18 if guard else 0);bc=value(anchor)
                    ts=[t for t in res if t[1]<=cap and(guard or not(t[0]=='D'and t[1]>=18))];ws=[value(anchor+[t])-bc for t in ts]
                    tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b;tabs.append(dict(a1=a1,guard=guard,table=tab))
                    for s in(range(18,rem-a1+1)if guard else[None]):
                        left=rem-a1-(s if guard else 0);z=tab['values'][left]
                        if z is None:continue
                        fs=fixed+[('A',1)]*a1+([('D',s)]if guard else[]);base=value(fs)
                        for t,w in zip(ts,ws):
                            if t[1]<=left:assert value(fs+[t])-base==w
                        opts.append((base+z,a1,guard,s,left,len(tabs)-1));options+=1
            z,a1,guard,s,left,ti=max(opts);forest=fixed+[('A',1)]*a1+([('D',s)]if guard else[]);tab=tabs[ti]['table']
            while left:t=tab['types'][tab['last'][left]];forest.append(t);left-=t[1]
            assert value(forest)==z and sum(t[1]for t in forest)==57
            rows.append(dict(fixed=fixed,tables=tabs,maximum=F(z,64),attaining_forest=forest))
        out.append(dict(fixed_A17_count=k17,fixed_rank16_count=k16,rows=rows,maximum=max(r['maximum']for r in rows)))
    prior=V.load(R.B.TAG)
    witness=next(r['attaining_forest']for c in prior['top_classes']['classes']for r in c['rows']if F(r['maximum'])==371 and value(list(map(tuple,r['attaining_forest'])))==64*371)
    return dict(old_above_371_classes=high,classes=out,inherited_three_rank16_class_upper='371',
        Bellman_inequalities=checks,attaining_predecessors=preds,explicit_options=options,
        attaining_forest=witness,conditional_complete_new_maximum='371')


def run():
    g=geometry();ret,ss=retention();upper=positive_family();top=finite()
    covers=list(map(F,V.load(R.TAG)['sequential_increasing_family_uppers']))+[upper['complete_positive_family_upper']]
    assert max(covers)<=371 and max(r['maximum']for r in top['classes'])<=371
    sources={Path(__file__),Path(R.__file__),Path(L.__file__),Path(S.__file__),S.OVER,
        ROOT/f'results/{R.TAG}.json',ROOT/f'results/{RANK15}.json'}
    sources.update(ROOT/f'results/{s}.json'for s in ss)
    return V.O.enc(dict(status='PASS_A17_INNER_UNIT_GUARD_SCOUT',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(sources)},
        rule='With a separate whole A16/D16/A17/D17 and D>=18, delete A17(1,15) unit and restore A15(0,14) unit under the same witnesses.',
        geometry=g,current_retention=ret,positive_family=upper,finite=top,sequential_increasing_family_uppers=covers,
        independent_review_complete=False,universal_composition_complete=False,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);R.M.R.resources();U.memory_guard();p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();R.M.R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'coupon',out['positive_family']['complete_positive_family_upper'],'maximum',out['finite']['conditional_complete_new_maximum'],
        'hits',len(out['current_retention']['target_hits']),[(r['fixed_A17_count'],r['fixed_rank16_count'],r['maximum'])for r in out['finite']['classes']],flush=True)
