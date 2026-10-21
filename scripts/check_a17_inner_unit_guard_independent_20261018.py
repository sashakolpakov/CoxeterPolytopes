#!/usr/bin/env python3
"""Independent A17 guarded-unit geometry and arithmetic certificate check."""
from collections import Counter
import fcntl
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
import signal
import check_d16_rank17_witness_extension_independent_20261018 as X

ROOT=X.ROOT;W=X.W;C=X.C;OLD=X.OLD
TAG='fresh_a17_inner_unit_rank16_rank17_guard_scout_20261018'


def value(fs):
    fs=list(map(tuple,fs));co=Counter(fs);r16=co['A',16]+co['D',16];r17=co['A',17]+co['D',17]
    g=any(k=='D'and n>=18 for k,n in fs)
    return X.value(fs)+(64*co['A',15]if g and r16+r17 else 0)-(64*co['A',17]if g and(r16 or r17>=2)else 0)


def geometry(doc):
    gp=doc['geometry'];p=gp['consumer'];g=[list(map(F,r))for r in p['gram']];w=p['white'];z=list(map(F,p['seal']['radical']))
    assert len(g)==19 and W.inertia(g)==[17,1,1]and W.inertia(W.sub(g,w))==[15,0,0]
    assert all(sum(g[i][j]*z[j]for j in range(19))==0 for i in range(19))
    black=sorted({i for c in p['corners']for i in c['marks']});assert len(black)==4
    for b in black:assert W.inertia(W.sub(g,w+[b]))==[16,0,0]
    pairs={tuple(sorted(c['marks']))for c in p['corners']};assert len(pairs)==4
    u=black[0];left=[u]+[b for b in black[1:]if tuple(sorted((u,b)))not in pairs];right=[b for b in black if b not in left]
    assert len(left)==len(right)==2 and pairs=={tuple(sorted((a,b)))for a in left for b in right}
    ids=sorted(set(w)|set(p['selected_marks']));om=sorted(set(black)-set(p['selected_marks']))
    allowed=X.K.contacts(g,ids,om);assert allowed==p['seal']['allowed_contacts']and all(z[i]for i in allowed)
    du,dv=W.distances(g,left),W.distances(g,right)
    shorts={(i,j):min(du[i]+dv[j],du[j]+dv[i])+2 for i,j in combinations(allowed,2)if z[i]*z[j]<0}
    assert shorts and max(shorts.values())<=8
    assert shorts=={tuple(r['nodes']):r['length']for r in p['seal']['mixed_sign_shortcuts']}
    old=new=0;debits=credits=0
    for c in p['corners']:
        marks=c['marks'];ids=sorted(set(w)|set(marks));assert W.inertia(W.sub(g,ids))==[17,0,0]
        bs=W.blocks(g,ids);b=next((b for b in bs if set(marks)<=set(b['nodes'])),None);im=[]if b is None else W.images(b,marks)
        assert im==list(map(tuple,c['canonical_images']))
        debit=b is not None and b['type']=='A17'and im==[(1,15)]
        credit=b is not None and b['type']=='A15'and im==[(0,14)]
        if credit:assert sorted(b['type']for b in bs)==['A1','A1','A15']
        if debit:assert len(bs)==1
        assert debit==c['debit']and credit==c['restore']and F(c['old_payment'])==int(debit)
        assert F(c['payment'])==int(credit);old+=F(c['old_payment']);new+=F(c['payment']);debits+=debit;credits+=credit
    assert(old,new,debits,credits)==(1,1,1,1)
    sup=W.load('fresh_a17_alternate_corner_pilot_20261018')['selected_support']
    assert sup['complete_eligible_patterns']==140 and sup['complete_maps']==48
    assert all([1,15]not in r['marked_images']for r in sup['rows'])
    hist=next(r for r in W.load('fresh_finite_d17_long_pair_support_20261017')['records']if r['kind']=='A')
    assert [1,15]in hist['removed_pairs']and [1,15]not in hist['conservatively_required_pairs']
    hits=doc['current_retention']['target_hits'];assert len(hits)==1
    h=hits[0];assert h['source']=='fresh_rank15_sparse_payment_exchange_20261018'
    ri=int(h['path'].split('/')[1]);rp=W.load(h['source'])['selected_proofs'][ri];assert rp['index']==1638
    assert F(h['saved']['payment'])==1 and h['images']==[[1,15]]and h['marks']==next(c['marks']for c in p['corners']if c['debit'])
    return len(shorts),len(p['corners'])+5


def affine_check(anchor,ts,ws,cap):
    bc=value(anchor);count=0
    for i,a in enumerate(ts):
        for j in range(i,len(ts)):
            b=ts[j]
            if a[1]+b[1]<=cap:assert value(anchor+[a,b])==bc+ws[i]+ws[j];count+=1
    return count


def positive(doc):
    up=doc['positive_family'];small=[t for t in C.TYPES if t[1]<=8 and t!=('A',1)]
    assert all(C.actual_cross(a,b)==0 for a in small for b in small)
    expected={(t,s)for t in(('A',16),('D',16),('A',17),('D',17))for s in range(18,43-t[1])}
    assert{(tuple(r['witness']),r['reserved_D'])for r in up['rows']}==expected
    checks=preds=pairs=0
    for r in up['rows']:
        witness=tuple(r['witness']);s=r['reserved_D'];left=42-witness[1]-s;best=None
        assert{t['a1']for t in r['tables']}==set(range(left+1))
        for t in r['tables']:
            a1=t['a1'];anchor=[('A',15),witness,('D',s)]+[('A',1)]*a1;cap=left-a1
            ts=[t for t in small if t[1]<=cap];ws=[value(anchor+[t])-value(anchor)for t in ts];tab=t['table']
            assert list(map(tuple,tab['types']))==ts and tab['weights']==ws and len(tab['values'])==cap+1
            pairs+=affine_check(anchor,ts,ws,cap);a,b=X.K.bellman(tab);checks+=a;preds+=b
            if tab['values'][cap]is not None:z=value(anchor)+tab['values'][cap];best=z if best is None else max(best,z)
        assert F(best,64)==F(r['maximum'])and value(r['attaining_forest'])==best
    assert(checks,preds)==(up['Bellman_inequalities'],up['attaining_predecessors'])
    assert max(F(r['maximum'])for r in up['rows'])==F(1481,4)==F(up['complete_positive_family_upper'])
    return checks,preds,pairs


def finite(doc):
    top=doc['finite'];old=W.load(X.TAG)['finite']['original_compact_class_maxima'];high=[r for r in old if F(r['maximum'])>371]
    assert high==top['old_above_371_classes']and len(high)==6
    assert{(r['counts'][7],r['counts'][8])for r in high}=={(0,0),(1,0),(2,0),(0,1),(1,1),(0,2)}
    assert all(sum(r['counts'][:7])==r['k5']==r['a12']==0 for r in high)
    previous=W.load(X.K.J);assert max(F(r['maximum'])for r in previous['top_classes']['classes'])==371
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS};res=[t for t in C.TYPES if C.representative(t)in reps and t!=('A',1)]
    assert ('A',15)not in res and ('A',17)not in res and ('D',17)in res
    assert all(C.actual_cross(a,b)==0 for a in res+[('A',1)]for b in res+[('A',1)])
    fibers=[t for t in C.TYPES if C.representative(t)==('D',16)];assert fibers==[('A',16),('B',16),('D',16)]
    checks=preds=options=pairs=0
    for cl in top['classes']:
        k17,k16=cl['fixed_A17_count'],cl['fixed_rank16_count']
        assert{tuple(map(tuple,r['fixed']))for r in cl['rows']}=={(('A',17),)*k17+f for f in combinations_with_replacement(fibers,k16)}
        for r in cl['rows']:
            fixed=list(map(tuple,r['fixed']));rem=57-17*k17-16*k16;best=None
            assert{(t['a1'],t['guard'])for t in r['tables']}=={(a,g)for a in range(rem+1)for g in(False,True)if not g or rem-a>=18}
            for t in r['tables']:
                a1=t['a1'];g=t['guard'];anchor=fixed+[('A',1)]*a1+([('D',18)]if g else[]);cap=rem-a1-(18 if g else 0)
                ts=[x for x in res if x[1]<=cap and(g or not(x[0]=='D'and x[1]>=18))];ws=[value(anchor+[x])-value(anchor)for x in ts];tab=t['table']
                assert list(map(tuple,tab['types']))==ts and tab['weights']==ws and len(tab['values'])==cap+1
                pairs+=affine_check(anchor,ts,ws,cap);a,b=X.K.bellman(tab);checks+=a;preds+=b
                for s in(range(18,rem-a1+1)if g else[None]):
                    n=rem-a1-(s if g else 0);z=tab['values'][n]
                    if z is None:continue
                    fs=fixed+[('A',1)]*a1+([('D',s)]if g else[]);bc=value(fs)
                    for x,w in zip(ts,ws):
                        if x[1]<=n:assert value(fs+[x])-bc==w
                    z+=bc;best=z if best is None else max(best,z);options+=1
            assert F(best,64)==F(r['maximum'])and value(r['attaining_forest'])==best
        assert max(F(r['maximum'])for r in cl['rows'])==F(cl['maximum'])
    assert(checks,preds,options)==(top['Bellman_inequalities'],top['attaining_predecessors'],top['explicit_options'])
    assert list(map(lambda r:F(r['maximum']),top['classes']))==[F(1479,4),F(741,2),F(739,2)]
    assert value(top['attaining_forest'])==64*371 and top['conditional_complete_new_maximum']=='371'
    return checks,preds,options,pairs


def main():
    doc=W.load(TAG);hashes=0
    for p,h in doc['sources'].items():assert W.digest(ROOT/p)==h;hashes+=1
    a=geometry(doc);b=positive(doc);c=finite(doc)
    assert list(map(F,doc['sequential_increasing_family_uppers']))==[F(739,2),F(365),F(371),F(11513,32),F(22917,64),F(1481,4)]
    OLD.U.memory_guard();print('PASS_INDEPENDENT_A17_INNER_UNIT_GUARD',hashes,*a,*b,*c,'371',flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);OLD.R.resources();main();OLD.R.resources(start=False)
