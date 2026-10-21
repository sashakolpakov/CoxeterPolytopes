#!/usr/bin/env python3
"""Independent new-rule check; imports frozen A16 baseline, no D16 producer."""
from collections import Counter
import fcntl
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
import signal
import check_a16_double_large_d_independent_20261018 as W
import fresh_a16_double_large_d_upper_scout_20261018 as OLD

ROOT=W.ROOT;C=OLD.C
N='fresh_d16_node3_rank16_large_d_guard_scout_20261018'
J='fresh_d16_node3_guard_complete_finite_20261018'
P='fresh_d16_node3_rank16_large_d_guard_transfer_20261018'


def value(fs):
    co=Counter(map(tuple,fs));r=co['A',16]+co['D',16];g=any(k=='D'and n>=18 for k,n in fs)
    return OLD.value(list(map(tuple,fs)))-(64*co['D',16]if g and r>=2 else 0)+(256*co['A',2]*co['A',13]if g and r and co['A',1]else 0)


def bellman(tab):
    ts=list(map(tuple,tab['types']));ws=tab['weights'];v=tab['values'];last=tab['last'];checks=preds=0
    assert v[0]==0 and last[0]is None
    for n in range(1,len(v)):
        possible=[i for i,t in enumerate(ts)if t[1]<=n and v[n-t[1]]is not None]
        for i in possible:assert v[n]>=v[n-ts[i][1]]+ws[i];checks+=1
        if possible:i=last[n];assert i in possible and v[n]==v[n-ts[i][1]]+ws[i];preds+=1
        else:assert v[n]is None and last[n]is None
    return checks,preds


def contacts(g,ids,omitted):
    out=list(omitted)
    for b in W.blocks(g,ids):
        inv=W.inverse(W.sub(g,b['nodes']));assert all(x>=0 for row in inv for x in row)
        out.extend(i for k,i in enumerate(b['nodes'])if inv[k][k]<=2)
    return sorted(out)


def paid_corners(g,w,cs):
    checks=0
    for c in cs:
        marks=c['marks'];ids=sorted(set(w)|set(marks));assert W.inertia(W.sub(g,ids))==[len(ids),0,0];checks+=1
        bs=W.blocks(g,ids);a,b=[next(b for b in bs if m in b['nodes'])for m in marks]
        if c['selected_D16_zero']:assert a is b and a['type']=='D16'and W.images(a,marks)==[(3,13)]
        if c['coupon']:
            assert a is not b and any(b['type']=='A1'for b in bs)
            sm,la,m,n=(a,b,*marks)if a['type']=='A2'else(b,a,marks[1],marks[0])
            assert sm['type']=='A2'and la['type']=='A13'and W.images(la,[n])==[(3,),(9,)]
        assert c['payment']==int(c['coupon'])
    assert sum(c['payment']for c in cs)==1
    return checks


def geometry(scout,transfer):
    gp=scout['geometry'];assert gp==transfer['complete_local_geometry'];bounds={r['index']:r for r in transfer['native_boundaries']};checks=0
    sup=W.load('fresh_a16_d16_a4_a12_fractional_transfer_20261018')['local_proofs']['domains']['D']['complete_selected_support']
    rows=[r for r in sup['rows']if[3,13]in r['marked_images']]
    assert sorted(r['index']for r in rows)==sorted(p['index']for p in gp['full']+gp['proper'])and sum(len(r['maps'])for r in rows)==48==gp['context_maps']
    assert sum(m['marked_images']==[3,13]for r in rows for m in r['maps'])==24==gp['affected_maps']
    for p in gp['full']:
        g=[list(map(F,r))for r in p['gram']];w=p['white'];bd=bounds[p['index']];sel=bd['selected_marks'];black=bd['left']+bd['right'];om=sorted(set(black)-set(sel))
        ids=sorted(set(w)|set(sel));allowed=contacts(g,ids,om);du,dv=W.distances(g,bd['left']),W.distances(g,bd['right']);seal=p['seal']
        for ii in[w]+[w+[b]for b in black]:assert W.inertia(W.sub(g,ii))==[len(ii),0,0];checks+=1
        if seal['method']=='radical':
            assert W.inertia(g)==[len(g)-2,1,1];z=list(map(F,seal['radical']));assert all(sum(g[i][j]*z[j]for j in range(len(g)))==0 for i in range(len(g)))
            assert allowed==seal['allowed_contacts']and all(z[i]for i in allowed)
            expected={(i,j):min(du[i]+dv[j],du[j]+dv[i])+2 for i,j in combinations(allowed,2)if z[i]*z[j]<0}
            assert expected and max(expected.values())<=8 and expected=={tuple(r['nodes']):r['length']for r in seal['mixed_sign_shortcuts']}
        else:
            inv=W.inverse(g);assert inv==[list(map(F,r))for r in seal['inverse']]
            assert W.inertia(g)==[17,1,0]and sorted(inv[i][i]for i in om)==[4,8]
            a,b=seal['end_proof']['end_sets'];assert sorted(a+b)==allowed and not set(a)&set(b)
            assert all(inv[i][j]>0 for side in(a,b)for i in side for j in side)
        checks+=paid_corners(g,w,p['corners'])
    for p in gp['proper']:
        src=next(r for r in rows if r['index']==p['index'])
        assert sorted(tuple(m['images'])for m in p['maps'])==sorted(tuple(m['images'])for m in src['maps']if m['marked_images']==[3,13])
        g=[list(map(F,r))for r in p['gram']];iv=W.inverse(g);assert iv==[list(map(F,r))for r in p['inverse']]and W.inertia(g)==[len(g)-1,1,0]
        bd=bounds[p['index']];allowed=contacts(g,sorted(set(p['white'])|set(bd['selected_marks'])),p['omitted_black']);assert allowed==p['possible_contacts']
        du,dv=W.distances(g,bd['left']),W.distances(g,bd['right']);a,b=p['end_proof']['end_sets']
        assert sorted(a+b)==allowed and not set(a)&set(b)and all(iv[i][j]>0 for side in(a,b)for i in side for j in side)
        for mp in p['maps']:
            image=dict(zip(mp['images'],p['selected_nodes']));missing=sorted(set(range(16))-set(image));assert missing==mp['missing']
            expected=[]
            for new in missing:
                for old,node in image.items():
                    if abs(new-old)==1 and max(new,old)<=13 or set((new,old))in({13,14},{13,15}):expected.append((new,old,node,iv[node][node]))
            assert sorted(expected)==sorted((e['new'],e['old'],e['core'],F(e['diagonal']))for e in mp['forced_entries'])
            assert all(F(e['diagonal'])==iv[e['core']][e['core']]for e in mp['forced_entries'])
            outcome=mp['outcome']
            if outcome=='projection_above_two':assert all(F(e['diagonal'])==8 for e in mp['forced_entries'])
            elif outcome=='two_orthogonal_norm_saturated_leaves':
                leaves=[e for e in mp['forced_entries']if e['new']in(14,15)];assert len(leaves)==2 and len({e['core']for e in leaves})==1 and all(F(e['diagonal'])==2 for e in leaves)
            elif outcome=='norm_saturated_root_cannot_continue':
                k=mp['forced_entries'][0]['core'];assert iv[k][k]==2
                neg=[(j,min(du[k]+dv[j],dv[k]+du[j])+3)for j in allowed if iv[k][j]<0]
                assert neg and max(v for j,v in neg)<=8
            else:
                assert outcome=='forced_single_leaf_completion';k=mp['forced_entries'][0]['core'];n=len(g);assert iv[k][k]==2
                gg=[r+[F(-int(i==k))]for i,r in enumerate(g)]+[[F(-int(i==k))for i in range(n)]+[F(2)]]
                assert gg==[list(map(F,r))for r in mp['completed_gram']]and W.inertia(gg)==[n-1,1,1]
                z=list(map(F,mp['radical']));assert all(sum(gg[i][j]*z[j]for j in range(n+1))==0 for i in range(n+1))
                assert sorted(z[i]for i in mp['remaining_contacts'])==[-2,2]
                checks+=paid_corners(gg,mp['completed_white'],mp['corners'])
    return checks


def finite(doc):
    top=doc['top_classes'];old=top['inherited_compact_class_maxima'];assert len(old)==38
    high=[r for r in old if F(r['maximum'])>F(1487,4)];assert high==top['inherited_high_classes']and len(high)==3
    assert{tuple(r['counts'])for r in high}=={tuple([0]*7+[k,0])for k in(0,1,2)}
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS};res=[t for t in C.TYPES if C.representative(t)in reps and t!=('A',1)]
    fibers=[t for t in C.TYPES if C.representative(t)==C.representative(('D',16))];checks=preds=options=0
    for cl in top['classes']:
        k=cl['fixed_rank16_count'];assert{tuple(map(tuple,r['fixed']))for r in cl['rows']}==set(combinations_with_replacement(fibers,k))
        for row in cl['rows']:
            fixed=list(map(tuple,row['fixed']));rem=57-16*k;expected={(a1,g)for a1 in range(rem+1)for g in(False,True)if not g or rem-a1>=18};assert{(r['a1'],r['guard'])for r in row['tables']}==expected
            best=None
            for t in row['tables']:
                a1=t['a1'];g=t['guard'];anchor=fixed+[('A',1)]*a1+([('D',18)]if g else[]);cap=rem-a1-(18 if g else 0)
                ts=[x for x in res if x[1]<=cap and(g or not(x[0]=='D'and x[1]>=18))];ws=[value(anchor+[x])-value(anchor)for x in ts];tab=t['table']
                assert list(map(tuple,tab['types']))==ts and tab['weights']==ws and len(tab['values'])==cap+1
                a,b=bellman(tab);checks+=a;preds+=b
                for s in(range(18,rem-a1+1)if g else[None]):
                    n=rem-a1-(s if g else 0);v=tab['values'][n]
                    if v is None:continue
                    chosen=fixed+[('A',1)]*a1+([('D',s)]if g else[]);z=value(chosen)+v
                    for x,w in zip(ts,ws):
                        if x[1]<=n:assert value(chosen+[x])-value(chosen)==w
                    best=z if best is None else max(best,z);options+=1
            assert F(best,64)==F(row['maximum'])and value(row['attaining_forest'])==best
        assert max(F(r['maximum'])for r in cl['rows'])==F(cl['maximum'])
    assert(checks,preds,options)==(top['Bellman_inequalities'],top['attaining_predecessors'],top['explicit_options'])
    assert max(F(c['maximum'])for c in top['classes'])==371
    assert doc['new_exact_global_formal_maximum']=='1487/4'and value(doc['attaining_forest'])==64*F(1487,4)
    return checks,preds,options


def coupon_upper(scout):
    up=scout['upper'];cs=(0,0,0,0,1,0,0,1,0)
    key=next(key for counts,used,key in C.fixed_domain()if counts==cs and key[9]==0)
    observed={(r['k5'],r['reserved_D'])for r in up['rows']}
    assert observed=={(k,s)for k in(0,1)for s in range(18,26-5*k)}
    tables={r['k5']:r['table']for r in up['tables']};checks=preds=0
    for k,tab in tables.items():
        ws={r:C.item(r,key,k)+(128 if r==1 else 256 if r in(2,4)else 0)for r in C.ITEMS}
        assert dict(tab['weights'])==ws;v=tab['values'];last=tab['last'];assert v[0]==0 and last[0]is None
        for n in range(1,58):
            possible=[r for r in ws if r<=n and v[n-r]is not None]
            for r in possible:assert v[n]>=v[n-r]+ws[r];checks+=1
            if possible:r=last[n];assert r in possible and v[n]==v[n-r]+ws[r];preds+=1
            else:assert v[n]is None and last[n]is None
        fixed=C.fixed_components(cs,key)+[('A',5)]*k;fc=C.upper_value(fixed)
        for row in[r for r in up['rows']if r['k5']==k]:
            s=row['reserved_D'];n=28-5*k-s;assert n==row['residual']>=3
            opts=[]
            for a1 in range(1,n-1):
                for a2 in range(1,(n-a1)//2+1):
                    rest=n-a1-2*a2
                    if v[rest]is not None:opts.append(fc+C.actual_base(('D',s))+a1*ws[1]+a2*ws[2]+v[rest])
            assert F(max(opts),64)==F(row['upper'])
    assert checks==up['Bellman_inequalities']and preds==up['predecessors']
    assert max(F(r['upper'])for r in up['rows'])==F(11513,32)
    return checks,preds


def main():
    scout,count,transfer=W.load(N),W.load(J),W.load(P);hashes=0
    for doc in(scout,count,transfer):
        for p,h in doc['sources'].items():assert W.digest(ROOT/p)==h;hashes+=1
    geom=geometry(scout,transfer);dp=finite(count);cover=coupon_upper(scout)
    ret=transfer['current_retention'];assert ret['later_corners']==564 and ret['literal_configurations']==216
    for r in ret['later_packets']:
        assert all([3,13]not in c['D16_images']for c in r['rows'])
    for r in ret['literal']:
        assert all([3,13]not in c['D16_images']for c in r['corners'])
    for r in ret['newer_seals']:assert all([3,13]not in c['D16_images']for c in r['rows'])
    assert all(r['D16_outputs']==0 for r in ret['chosen_D12'])
    assert list(map(F,count['sequential_increasing_family_uppers']))==[F(739,2),F(365),F(371),F(11513,32)]
    OLD.U.memory_guard();print('PASS_INDEPENDENT_D16_NODE3_GUARD',hashes,geom,*dp,*cover,'1487/4',flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);OLD.R.resources();main();OLD.R.resources(start=False)
