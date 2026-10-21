#!/usr/bin/env python3
"""Retain A1/D5 long-end payment only with its sealed whole witnesses.

The unchanged complete count gives a monotone global upper. Fresh small
actual-type tables settle exactly the two affected D12 obstruction classes.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement,product
import fresh_quarter_center_rank13_ambient58_frontier_20261018 as C
import fresh_quarter_center_rank13_composed_ambient58_frontier_20261018 as G
import fresh_d5_d12_alternate_coupon_rebalance_20261018 as O

ROOT=C.ROOT;B=C.R.B;R=C.R;T=R.T;A=R.A
TAG='fresh_a1_d5_witness_context_retention_20261018'
STEMS=tuple(dict.fromkeys((C.TAG,G.TAG,O.TAG,R.TAG,C.Q.TAG,B.TAG,B.PART,B.REST,B.LOW,B.IDEAL,
    *O.STEMS,*B.STEMS,*R.STEMS)))
enc=O.enc

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def guard_types(types):
    co=Counter(types)
    return co['A11']>0 or(co['A9']>0 and co['A2']>0)
def guard(fs):
    co=Counter(fs)
    return co['A',11]>0 or(co['A',9]>0 and co['A',2]>0)
def ispair(a,b):return frozenset((a,b))==frozenset((('A',1),('D',5)))
def actual_base(t):return C.actual_base(t)
def actual_cross(a,b):return C.actual_cross(a,b)-(24 if ispair(a,b)else 0)
def context(fs):
    co=Counter(fs)
    return C.context(fs)+(24*co['A',1]*co['D',5]if guard(fs)else 0)
def value(fs):
    return sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)

def debit(blocks,marks,g=None):
    if not blocks:return F()
    aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
    if aa is bb or{aa['type'],bb['type']}!={'A1','D5'}:return F()
    if guard_types(b['type']for b in blocks):return F()
    d5,mark=(aa,marks[0])if aa['type']=='D5'else(bb,marks[1])
    if g is None:
        images=d5['images'];assert all(len(x)==1 for x in images)
        long=images==[[0]]
    else:long=T.mark_pairs(g,d5['nodes'],[mark],5,'D')==[(0,)]
    return -F(3,8)if long else F()

def walk(x,path=''):
    if isinstance(x,dict):
        if 'marks'in x and len(x['marks'])==2 and 'payment'in x and 'roles'in x:
            yield path,x
        for k,v in x.items():yield from walk(v,path+'/'+k)
    elif isinstance(x,list):
        for i,v in enumerate(x):yield from walk(v,path+'/'+str(i))

def local(docs):
    inherited=O.local_seals(docs)
    assert enc(inherited)==docs[O.TAG]['local_proofs']
    proofs=[]
    for p in inherited['proofs']:
        corners=[]
        for c in p['corners']:
            types=[b['type']for b in c['components']]
            if any(r['role']=='D5_long_A1_credit'for r in c['roles']):
                assert guard_types(types)and Counter(types)=={'A1':1,'A11':1,'D5':1}
            delta=debit(c['components'],c['marks'],A.matrix(p['old_seal']['gram']))
            assert delta==0
            corners.append(dict(marks=c['marks'],types=types,retained_lower=c['payment']))
        assert sum(c['retained_lower']for c in corners)==1
        proofs.append(dict(index=p['index'],corners=corners,total='1'))
    return proofs

def interfaces(docs):
    literals=[];changed=[]
    for prev,(stem,case,g,white,cs)in zip(docs[R.TAG]['literal_interfaces']['rows'],T.literal_inputs(docs),strict=True):
        assert(prev['source'],prev['case'])==(stem,enc(case))
        rows=[]
        for c,(marks,*_)in zip(prev['corners'],cs,strict=True):
            assert c['marks']==marks
            delta=debit(c['components'],marks,g);lower=F(c['parameter_values']['9/64'])+delta
            assert lower>=0
            rows.append(dict(marks=marks,before=c['parameter_values']['9/64'],debit=delta,lower=lower))
        total=sum(c['lower']for c in rows);assert total>=1
        row=dict(source=stem,case=enc(case),corners=rows,total=total);literals.append(row)
        if any(c['debit']for c in rows):changed.append(row)
    assert len(literals)==216 and len(changed)==2
    assert [r['case']for r in changed]==[[2372,0],[2540,0]]
    assert all(r['source']=='fresh_d16_partial_weight_transfer_20261017'and r['total']==1 for r in changed)
    later=[];hits=[];total=0
    for stem,keys in B.FINITE.items():
        count=0
        for section in keys:
            for path,c,gg,white,blocks in B.corners(docs[stem][section],section):
                g=A.matrix(gg);marks=c['marks'];assert sorted(A.components(g,sorted(set(white)|set(marks))))==sorted(b['nodes']for b in blocks)
                assert debit(blocks,marks,g)==0,(stem,path)
                aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
                if aa is not bb and{aa['type'],bb['type']}=={'A1','D5'}:
                    types=[b['type']for b in blocks];assert guard_types(types)
                    hits.append(dict(source=stem,path=path,marks=marks,types=types))
                count+=1;total+=1
        later.append(dict(source=stem,corners=count))
    assert len(later)==22 and total==564 and len(hits)==8
    chosen=[]
    for tag in(B.PART,B.REST):
        rows=[]
        for path,c in walk(docs[tag]):
            if F(c['payment'])<=0:continue
            assert not{'A1','D5'}<={b['type']for b in c['roles']},(tag,path)
            rows.append(dict(path=path,payment=c['payment'],reason=c['reason'],marked_types=[b['type']for b in c['roles']]))
        assert rows;chosen.append(dict(source=tag,positive_outputs=rows))
    for p in docs[C.Q.TAG]['whole_variants']:
        assert all(debit(c['components'],c['marks'])==0 for c in p['corners'])
    for p in docs[R.TAG]['local_seals']:
        assert all(debit(c['components'],c['marks'])==0 for c in p['corners'])
    return dict(literals=literals,changed_literals=changed,later_packets=later,later_corners=total,
        guarded_later_A1_D5_outputs=hits,D12_positive_outputs=chosen,quarter_variants=24,rank13_seals=4)

def refinements(docs,saved=None):
    old=docs[C.TAG];ref=old['actual_typed_refinement'];ff=old['finite']
    cs=(0,0,0,3,0,0,0,0,0);targets={(cs,3,0),(cs,4,0)}
    classes=[r for r in ref['classes']if(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count'])in targets]
    assert len(classes)==2
    fibers={}
    for t in C.TYPES:fibers.setdefault(C.representative(t),[]).append(t)
    savedrows={(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count']):r for r in saved}if saved else{}
    out=[];checks=attainers=0
    for row in classes:
        key=(tuple(row['counts']),row['fixed_A5_count'],row['fixed_A12_count'])
        upper=next(x for x in ff['above_curvature']if(tuple(x['counts']),x['fixed_A5_count'],x['fixed_A12_count'])==key)
        fixedrep=C.fixed_components(tuple(row['counts']),upper['key'])+[('A',5)]*row['fixed_A5_count']
        groups=sorted(Counter(fixedrep).items())
        fixedsets=sorted(tuple(sorted(t for group in choice for t in group))for choice in product(*(combinations_with_replacement(fibers[t],n)for t,n in groups)))
        assert fixedsets==sorted(tuple(map(tuple,x['fixed_actual']))for x in row['rows'])
        rem=58-sum(t[1]for t in fixedrep);assert rem in(2,7)
        types=[t for t in C.TYPES if t[1]<=rem and t[1]!=5]
        assert all(actual_cross(a,b)==0 for a in types for b in types)
        prior={tuple(map(tuple,x['fixed_actual'])):x for x in savedrows[key]['rows']}if saved else{}
        rows=[]
        for fixed in fixedsets:
            assert all(t[1]in(5,12)for t in fixed)
            # Residual rank at most7 cannot supply either A9 or A11; both
            # old contexts and the new A1/D5 context are identically zero.
            ws=[actual_base(t)+sum(actual_cross(t,b)for b in fixed)for t in types]
            if saved:
                v,last=prior[fixed]['values'],prior[fixed]['last']
            else:
                v=[0]+[-10**9]*rem
                for t,w in zip(types,ws):
                    for n in range(t[1],rem+1):v[n]=max(v[n],v[n-t[1]]+w)
                last=[None]+[next(i for i,(t,w)in enumerate(zip(types,ws))if t[1]<=n and v[n]==v[n-t[1]]+w)for n in range(1,rem+1)]
            assert len(v)==len(last)==rem+1 and v[0]==0 and last[0]is None
            for n in range(1,rem+1):
                for t,w in zip(types,ws):
                    if t[1]<=n:assert v[n]>=v[n-t[1]]+w;checks+=1
                i=last[n];assert 0<=i<len(types)and types[i][1]<=n and v[n]==v[n-types[i][1]]+ws[i];attainers+=1
            fs=list(fixed);left=rem
            while left:
                t=types[last[left]];fs.append(t);left-=t[1]
            base=sum(actual_base(t)for t in fixed)+sum(actual_cross(a,b)for i,a in enumerate(fixed)for b in fixed[i+1:])
            best=base+v[rem]
            assert not guard(fs)and context(fs)==0 and value(fs)==best and sum(n for k,n in fs)==58
            record=dict(fixed_actual=fixed,types=types,weights=ws,values=v,last=last,maximum=best,forest=fs)
            if saved:assert enc(record)==prior[fixed]
            rows.append(record)
        maximum=max(r['maximum']for r in rows)
        out.append(dict(counts=row['counts'],fixed_A5_count=row['fixed_A5_count'],fixed_A12_count=0,
            previous_maximum=row['actual_maximum'],new_maximum=F(maximum,64),rows=rows,
            maximizers=[r['forest']for r in rows if r['maximum']==maximum]))
    assert sorted(F(r['new_maximum'])for r in out)==[F(372),F(376)]
    others=[r for r in ref['classes']if r not in classes and r['actual_obstruction_retained']]
    assert len(others)==10
    for row in others:
        assert any(value(list(map(tuple,fs)))==64*F(row['actual_maximum'])for fs in row['maximizers'])
    witness=[('A',16)]*2+[('D',26)]
    assert value(witness)==C.value(witness)==64*F(761,2)
    return dict(classes=out,newly_below_curvature=2,retained_above_curvature=10,
        actual_fixed_assignments=sum(len(r['rows'])for r in out),Bellman_inequalities=checks,
        attaining_predecessors=attainers,global_exact_maximum='761/2',global_witness=witness,
        monotone_upper_from_complete_predecessor=True,other_above_curvature_class_maxima_unchanged=True,
        full_global_DP_repeated=False)

def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    local_proofs=local(docs);inter=interfaces(docs)
    for a in C.TYPES:
        for b in C.TYPES:assert actual_cross(a,b)==actual_cross(b,a)>=0
    refine=refinements(docs,saved['finite_refinement']['classes']if saved else None)
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_A1_D5_WHOLE_WITNESS_CONTEXT_RETENTION',
        sources={p:digest(p)for p in sorted(paths)},baseline=G.TAG,
        rule='A1/D5 long-end3/8 only with wholeA11 OR(wholeA9 ANDwholeA2)',
        once_per_marked_pair=True,local_original_seals=local_proofs,interfaces=inter,
        finite_refinement=refine,pointwise_total_cost_nonincreasing=True,
        ideal_weights_changed=False,same_component_weights_changed=False,
        triangle_minimum='2',quadrilateral_minimum='1',completed_dimension_bound=58,new_dimension_bound=None,
        scope='Chosen-output retention of current complete face proof, four original seals,216/564/24/4 saved interfaces, monotone global380.5 and fresh exact refinement of two D12 obstruction classes. Full arbitrary-support arguments are stated in the report; no new dimension bound or example.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists()
    saved=json.loads(p.read_text())if args.verify else None;out=run(saved)
    if args.verify:assert out==saved
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite_refinement']['global_exact_maximum'],out['finite_refinement']['retained_above_curvature'],flush=True)
