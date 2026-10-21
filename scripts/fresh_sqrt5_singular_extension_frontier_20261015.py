#!/usr/bin/env python3
"""Complete next-white domains on the51 new singular unpaid Qsqrt5 cores.

Finite Dynkin and field/inertia arithmetic use the standalone independent
Fraction helper. Marked canonicalization and certified affine-witness search
are inherited read-only. No arbitrary-depth closure is inferred for survivors.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_quad_contacts_20261015_review as Q
import fresh_sqrt5_quad_contacts_20261015 as OLD
import fresh_sqrt5_all_quad_unions_20261015_classes as C
import fresh_sqrt5_sealed_weights_20261015 as T

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_singular_extension_frontier_20261015'

def extension(row,c):
    n=row['nodes'];g=Q.edges(row);g.update({(i,n):m for i,m in enumerate(c)})
    return dict(nodes=n+1,white=list(range(4,n+1)),labels=[g[p]for p in combinations(range(n+1),2)])
def corners(row):
    n=row['nodes'];g=Q.edges(row)
    return all(Q.finite(row['white']+[a,b],g)for a,b in product((0,2),(1,3)))
def short(row):
    g=Q.edges(row)
    return any((lambda d:d is not None and d<=2)(Q.distance(Q.graph(row['white']+[a,b],g),a,b))for a,b in product((0,2),(1,3)))
def profiles(row):
    n=row['nodes'];g=Q.edges(row);c=[2]*n;out=[];count=Counter()
    def visit(k):
        count['nodes']+=1
        if k==n:
            if any(x!=2 for x in c):out.append(c[:])
            return
        for label in(2,3,4,5,6,10):
            c[k]=label;g[k,n]=label;count['trials']+=1
            if label!=2:
                ids=list(range(4,n+1));bad=False
                for a,b in product((0,2),(1,3)):
                    H=Q.graph(ids+[a,b],g);d=Q.distance(H,a,b)
                    if d is not None and d<=2:count['short_prunes']+=1;bad=True;break
                    if not Q.finite(ids+[a,b],g):count['finite_prunes']+=1;bad=True;break
                if bad:continue
            visit(k+1)
        del g[k,n];c[k]=2
    assert corners(row)and not short(row);visit(0)
    assert len(out)==len(set(map(tuple,out)))
    return out,dict(count)

def affine(row):
    w=OLD.affine_witness(row)
    if w:
        g=Q.edges(row);ids=w['nodes'];sub={(i,j):g[Q.pair(ids[i],ids[j])]for i,j in combinations(range(len(ids)),2)}
        M=Q.tree_scaled(len(ids),sub);assert M is not None and Q.inertia(M)==[len(ids)-1,0,1]
        v=Q.kernel(M);assert len(v)==1 and (all(x.sign()>0 for x in v[0])or all(x.sign()<0 for x in v[0]))
    return w

def payment(row,weight_data,sealed):
    bonus={o['key']:F(o['weight'])for r in sealed['component_rows']for o in r['pair_orbits']}
    pure=set(weight_data['pure_d7_orbits']);half=set(weight_data['half_d7_orbits']);out=[];g=Q.edges(row)
    for a,b in product((0,2),(1,3)):
        H=T.graph(row['white']+[a,b],g);unseen=set(H);components=[]
        while unseen:
            ids=T.component(H,min(unseen));unseen.difference_update(ids)
            components.append({i:{j:m for j,m in H[i].items()if j in ids}for i in ids})
        a1=any(len(h)==1 for h in components);J=next(h for h in components if a in h);vals=[]
        if b not in J:vals=[F(),F()];typ=None;key=None
        else:
            typ=T.typename(J);key=T.keystr(T.tree_key(J,(a,b)))
            for outside_a1 in(False,True):
                extra=bonus.get(key,F())
                if typ=='D7'and not(a1 or outside_a1)and key not in pure:extra=F(1,2)if key in half else F()
                vals.append(F(T.baseline(J,a,b))+extra)
        out.append(dict(marked=[a,b],forest=[T.typename(h)for h in components],marked_component=typ,pair_orbit=key,core_has_a1=a1,weights=list(map(str,vals))))
    sums=[sum(F(r['weights'][i])for r in out)for i in(0,1)]
    return dict(corners=out,sums_without_with_external_a1=list(map(str,sums)),paid=min(sums)>=1)

def run():
    prior=json.loads((ROOT/'results/fresh_sqrt5_sealed_subcore_closure_20261015.json').read_text())
    base=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())['labelled_union_supports']
    forest=json.loads((ROOT/'results/fresh_sqrt5_forest_sensitive_weights_20261015.json').read_text())
    sealed=json.loads((ROOT/'results/fresh_sqrt5_sealed_weights_20261015.json').read_text())
    classes=[];lookup={}
    selected=[r for r in prior['records']if r['status']=='OPEN_SINGULAR'];assert len(selected)==61
    for ri,r in enumerate(selected):
        row=extension(base[r['source_union_index']],r['contacts']);key,order=C.canonical(row)
        if key not in lookup:
            lookup[key]=len(classes);classes.append(dict(row=row,canonical_order=order,origins=[]))
        rec=classes[lookup[key]];mp=[None]*row['nodes']
        for a,b in zip(order,rec['canonical_order']):mp[a]=b
        C.check_map(row,rec['row'],mp)
        rec['origins'].append(dict(open_index=ri,source_class_index=r['source_class_index'],source_union_index=r['source_union_index'],trial_index=r['trial_index'],to_representative=mp))
    assert len(classes)==51
    records=[];counts=Counter();trialcounts=Counter()
    for ci,rec in enumerate(classes):
        row=rec['row'];n=row['nodes'];g=Q.edges(row);M=Q.tree_scaled(n,g);assert M is not None
        assert Q.inertia(M)==[n-2,1,1]and Q.inertia([[x.conj()for x in r]for r in M])==[n-1,0,1]
        assert not affine(row)and not OLD.h4_payments(row)and not short(row)
        domain,dc=profiles(row);trials=[];accepted=[]
        for contact in domain:
            ext=extension(row,contact);paid=OLD.h4_payments(ext);trial=dict(contacts=contact)
            if paid:trial.update(status='BECOMES_PAID_ENTIRE_H4',witness=paid)
            elif(w:=affine(ext)):trial.update(status='REJECT_AFFINE',witness=w)
            else:
                G=Q.tree_scaled(n+1,Q.edges(ext))
                if G is None:trial['status']='REJECT_FINITE_CYCLE_FIELD'
                else:
                    K=Q.kernel([r[:-1]for r in G[:-1]]);assert len(K)==1
                    products=[sum((x*G[i][-1]for i,x in enumerate(v)),Q.K())for v in K]
                    if any(products):trial.update(status='REJECT_PHYSICAL_RELATION_AND_CONJUGATE_KERNEL',products=[z.enc()for z in products])
                    else:
                        pi=Q.inertia(G);co=Q.inertia([[z.conj()for z in r]for r in G]);trial.update(physical_inertia=pi,conjugate_inertia=co)
                        if pi[1]>1:trial['status']='REJECT_TWO_NEGATIVE'
                        elif co[1]:trial['status']='REJECT_CONJUGATE_NOT_PSD'
                        else:
                            assert pi==[n-1,1,1]and co==[n,0,1]
                            trial['status']='PASS_UNPAID_EXTRA_WHITE';accepted.append(contact)
            trials.append(trial);trialcounts[trial['status']]+=1
        weights=payment(row,forest,sealed)
        status='OPEN_SINGULAR_WITH_CONTACTS'if accepted else'SEALED_PAID_BY_FOREST_RULE'if weights['paid']else'SEALED_UNPAID_BY_CURRENT_FOREST_RULE'
        counts[status]+=1;records.append(dict(index=ci,**rec,profile_counts=dc,profiles=len(domain),trials=trials,accepted=accepted,payment=weights,status=status))
        print('class',ci,status,'profiles',len(domain),'accepted',len(accepted),'weights',weights['sums_without_with_external_a1'],flush=True)
    assert dict(counts)==dict(SEALED_PAID_BY_FOREST_RULE=32,OPEN_SINGULAR_WITH_CONTACTS=5,SEALED_UNPAID_BY_CURRENT_FOREST_RULE=14)
    assert sum(trialcounts.values())==502 and trialcounts['PASS_UNPAID_EXTRA_WHITE']==7
    assert [r['index']for r in records if r['status']=='SEALED_UNPAID_BY_CURRENT_FOREST_RULE']==[15,16,18,28,31,32,33,37,38,39,43,44,45,46]
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_quad_contacts_20261015_review.py','scripts/fresh_sqrt5_quad_contacts_20261015.py','scripts/fresh_sqrt5_all_quad_unions_20261015_classes.py','scripts/fresh_sqrt5_sealed_weights_20261015.py']
    for p in('fresh_sqrt5_sealed_subcore_closure_20261015','fresh_sqrt5_quad_contacts_20261015_review','fresh_sqrt5_all_quad_unions_20261015_review_complete','fresh_sqrt5_forest_sensitive_weights_20261015','fresh_sqrt5_sealed_weights_20261015','fresh_sqrt5_quad_exhaustion_20261015'):paths.append(f'results/{p}.json')
    Q.manifest_check(paths)
    return dict(status='PASS_COMPLETE_ONE_WHITE_DOMAINS_ON_51_NEW_SINGULAR_CORES',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        input_literal_extensions=61,marked_classes=51,counts=dict(counts),trial_counts=dict(trialcounts),records=records,
        scope='Complete finite-alphabet next-white necessary domains on the51 singular marked classes from61 first extensions. Short marked paths and entireH4 give persistent baseline payment. Every other profile is tested for a certified affine principal diagram, exact cyclic field, physical relation/conjugate kernel and full physical/conjugate inertia. Sealed cores have all external white nodes orthogonal unless already baseline-paid, so their entire corner components are fixed; both possible externalA1 flags are checked under the frozen forest-sensitive rule. Unpaid sealed cases are weight constraints, not claims of actual compact realization. Nonsealed cases, nonsingular first extensions and dotted parameters remain open; no global dimension bound.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text());print('PASS saved singular-frontier replay')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],data['trial_counts'])
