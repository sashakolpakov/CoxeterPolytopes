#!/usr/bin/env python3
"""Remove D17(2,12), paid by a whole-A9/A2 guarded A1/D5 unit.

Complete selected/output retention, a fresh conditional upper on every
potentially increasing class, and exact actual refinements of all39 old
above-curvature classes. No frozen module state is changed.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement,product
import fresh_rank15_a1_context_composed_ambient58_frontier_20261018 as G
import fresh_d17_existing_three_corner_transfer_20261018 as V
import fresh_order15_d16_d17_ambient59_frontier_20261018 as J

ROOT=G.ROOT;C=G.C;L=G.L;W=G.W;B=W.B;R=W.R;T=V.T;A=V.A;Q=V.Q;D=V.D
TAG='fresh_d17_a9_a2_context_zero_ambient58_frontier_20261018'
PAIR=(2,12);CASES=[1062,2525];N=58;SCALE=64
STEMS=tuple(dict.fromkeys((G.TAG,C.TAG,L.TAG,W.TAG,V.TAG,J.TAG,*W.STEMS,*L.STEMS,*B.STEMS)))
enc=W.enc

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return G.actual_base(t)-(40 if t==('D',17)else 0)
actual_cross=G.actual_cross
def context(fs):
    co=Counter(fs)
    return G.context(fs)+(40*co['A',1]*co['D',5]if co['A',9]and co['A',2]else 0)
def value(fs):
    return sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)

def changes(g,white,marks):
    ids=sorted(set(white)|set(marks))
    if Q.inertia(A.sub(g,ids))!=[len(ids),0,0]:return dict(delta=F(),roles=[],components=[])
    blocks=D.finite_blocks(g,ids);aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
    delta=F();roles=[]
    if aa is bb and aa['type']=='D17':
        im=T.mark_pairs(g,aa['nodes'],marks,17,'D')
        if PAIR in im:
            assert im==[PAIR];delta=-F(5,8);roles.append(dict(role='D17_node2_debit',images=im))
    elif aa is not bb and {'A1','D5'}=={aa['type'],bb['type']}:
        types={b['type']for b in blocks}
        small,large,m,lm=(aa,bb,*marks)if aa['type']=='A1'else(bb,aa,marks[1],marks[0])
        im=T.mark_pairs(g,large['nodes'],[lm],5,'D')
        if {'A9','A2'}<=types and im==[(0,)]:
            delta=F(5,8);roles.append(dict(role='guarded_whole_A1_D5_unit_increment',images=im))
    return dict(delta=delta,roles=roles,components=blocks)

def local(docs):
    over=docs[A.OVER];sup=T.selected(over,17)
    assert enc(sup)==docs[V.TAG]['complete_selected_support']
    assert(sup['complete_eligible_patterns'],sup['complete_maps'])==(140,168)
    affected=[r for r in sup['rows']if PAIR in map(tuple,r['marked_images'])]
    assert [r['index']for r in affected]==CASES
    assert all(len(r['base_nodes'])==17 and set(map(tuple,r['marked_images']))=={PAIR}for r in affected)
    proofs=[]
    for r in affected:
        c=over['isolation']['cases'][r['index']];old=V.seal(c,over['isolation']['cores'][c['core_key']])
        assert enc(old)==next(p for p in docs[V.TAG]['selected_proofs']if p['index']==r['index'])
        g=A.matrix(old['gram']);rows=[]
        for p in old['corners']:
            ch=changes(g,old['white'],p['marks']);new=p['payment']+ch['delta'];assert new>=0
            if ch['delta']>0:
                assert Counter(b['type']for b in ch['components'])=={'A1':1,'A2':1,'A9':1,'D5':1}
                assert p['reason']=='existing_D5_long_A1'and new==1
            if ch['delta']<0:assert p['reason']=='changed_selected_D17'and new==0
            rows.append(dict(marks=p['marks'],previous=p['payment'],**ch,payment=new))
        assert sorted(p['payment']for p in rows)==[0,0,0,1]
        assert sorted(p['delta']for p in rows)==[-F(5,8),0,0,F(5,8)]
        proofs.append(dict(index=r['index'],old_seal=old,gram=old['gram'],white=old['white'],corners=rows,total='1'))
    graph=T.target(17);adj,dist,_=V.P.distance_data(graph);autos=V.P.embeddings(graph,graph)
    assert dist[PAIR[0]][PAIR[1]]==10
    assert all(tuple(sorted(m[i]for i in PAIR))==PAIR for m in autos)
    weights={tuple(r['pair']):F(r['weight'])for r in docs[J.TAG]['combined_D17_weighted_pairs']}
    assert sum(weights.values())==F(G.actual_base(('D',17)),64)==F(919,8)
    assert weights[PAIR]==F(5,8);weights[PAIR]=F()
    assert sum(weights.values())==F(actual_base(('D',17)),64)==F(457,4)
    return dict(complete_selected_support=sup,affected_contexts=CASES,
        affected_maps=sum(len(r['maps'])for r in affected),proofs=proofs,
        old_D17_cost='919/8',new_D17_cost='457/4',new_weight='0',
        weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(weights.items())])

def interfaces(docs):
    literals=[]
    for prev,(stem,index,g,white,cs)in zip(docs[R.TAG]['literal_interfaces']['rows'],T.literal_inputs(docs),strict=True):
        assert(prev['source'],prev['case'])==(stem,enc(index));rows=[]
        for p,(marks,*_)in zip(prev['corners'],cs,strict=True):
            ch=changes(g,white,marks);assert ch['delta']>=0
            rows.append(dict(marks=marks,delta=ch['delta'],roles=ch['roles']))
        literals.append(dict(source=stem,case=index,corners=rows))
    assert len(literals)==216
    later=[];negative=[];total=0
    for stem,keys in B.FINITE.items():
        rows=[]
        for section in keys:
            for path,p,gg,white,blocks in B.corners(docs[stem][section],section):
                ch=changes(A.matrix(gg),white,p['marks'])
                assert sorted(b['nodes']for b in ch['components'])==sorted(b['nodes']for b in blocks)
                if ch['delta']<0:
                    assert stem==V.TAG and F(p['payment'])==F(5,8)
                    negative.append(dict(source=stem,path=path,marks=p['marks'],delta=ch['delta']))
                rows.append(dict(path=path,marks=p['marks'],delta=ch['delta'],roles=ch['roles']));total+=1
        later.append(dict(source=stem,corners=rows))
    assert len(later)==22 and total==564 and len(negative)==2
    additional=[]
    for stem,section,expected in((L.TAG,'local_proofs',24),(R.TAG,'local_seals',16)):
        # B.corners carries physical grams and white roots into every terminal corner.
        rows=[]
        for path,p,gg,white,blocks in B.corners(docs[stem][section],section):
            ch=changes(A.matrix(gg),white,p['marks']);assert ch['delta']>=0
            rows.append(dict(path=path,delta=ch['delta']))
        assert len(rows)>=expected
        additional.append(dict(source=stem,checked_corners=len(rows),rows=rows))
    variants=[]
    for p in docs[L.TAG]['current_interfaces']['retained_D14_variants']:
        assert all(F(c['delta'])>=0 for c in p['corners'])
    for p in docs[G.Q.TAG]['whole_variants']:
        # These records carry complete canonical images, but not a physical Gram.
        # Any D17 debit would require both marks in a D17 block.
        assert all(not any(b['type']=='D17'and set(c['marks'])<=set(b['nodes'])for b in c['components'])for c in p['corners'])
        variants.append(dict(index=p['index'],extra_contact=p['extra_contact'],negative_new_outputs=0))
    assert len(variants)==24
    chosen=[]
    for stem,expected in((B.PART,280),(B.REST,52)):
        rows=[(path,p)for path,p in W.walk(docs[stem])if F(p['payment'])>0]
        assert len(rows)==expected
        assert all(not any(b['type']=='D17'for b in p['roles'])for _,p in rows)
        chosen.append(dict(source=stem,positive_outputs=len(rows),D17_outputs=0))
    return dict(literal_configurations=216,literals=literals,later_packets=later,
        later_corners=total,replaced_own_negative_outputs=negative,
        additional_current_seals=additional,rank15_guard_variants=variants,
        chosen_D12_outputs=chosen,all_other_current_chosen_outputs_nondecreasing=True)

def conditional_upper(saved=None):
    """All forests with an A9 and D5, with an unconditional extra5/8*A1*D5.

    The old representative A1/rank5 allowance3/8 remains, so this relaxed
    polynomial dominates every potentially increasing actual forest.
    """
    domain=[(cs,u,key,k5)for cs,u,key in C.fixed_domain()if cs[0]
        for k5 in range(1,(N-u)//5+1)]
    assert len(domain)==759
    previous={(tuple(r['counts']),tuple(r['key']),r['fixed_A5_count']):r for r in saved['rows']}if saved else{}
    expected={(cs,key,k5)for cs,u,key,k5 in domain}
    if saved:assert set(previous)==expected and len(previous)==len(saved['rows'])
    rows=[];checks=preds=identities=options_count=0
    for cs,used,key,k5 in domain:
        fixed=C.fixed_components(cs,key)+[('A',5)]*k5;rem=N-used-5*k5
        weights={r:C.item(r,key,k5)+(40*k5 if r==1 else 0)for r in C.ITEMS}
        for r,w in weights.items():
            t=C.representative(('A'if r<=8 else'D',r))
            assert w==C.upper_base(t)+sum(C.upper_cross(t,b)for b in fixed)+(40*k5 if r==1 else 0)
            identities+=1
        old=previous.get((cs,key,k5));tables=[]
        for no2 in(False,True):
            ws={r:w for r,w in weights.items()if not no2 or r!=2}
            if saved:v,last=old['tables'][int(no2)]['values'],old['tables'][int(no2)]['last']
            else:
                v=[0]+[-10**9]*N
                for r,w in ws.items():
                    for n in range(r,N+1):v[n]=max(v[n],v[n-r]+w)
                last=[0]+[next(r for r,w in ws.items()if r<=n and v[n]==v[n-r]+w)for n in range(1,N+1)]
            assert len(v)==len(last)==N+1 and v[0]==last[0]==0
            for n in range(1,N+1):
                for r,w in ws.items():
                    if r<=n:assert v[n]>=v[n-r]+w;checks+=1
                r=last[n];assert r in ws and r<=n and v[n]==v[n-r]+ws[r];preds+=1
            tables.append(dict(exclude2=no2,values=v,last=last))
        opts=[dict(value=tables[0]['values'][rem],reserve=0,a2=0,no2=False)]
        for reserve in(1,2):
            oldguard=reserve==2 and bool(cs[1])
            if rem<reserve or(reserve==2 and not oldguard):continue
            for j in range((rem-reserve)//2+1):
                z=reserve*weights[1]+j*weights[2]+tables[1]['values'][rem-reserve-2*j]+64*j*k5
                if oldguard:z+=36*j*(j-1)
                opts.append(dict(value=z,reserve=reserve,a2=j,no2=True))
        options_count+=len(opts);win=max(opts,key=lambda z:(z['value'],z['reserve'],z['a2']))
        fs=fixed+[('A',1)]*win['reserve']+[('A',2)]*win['a2'];left=rem-win['reserve']-2*win['a2']
        tab=tables[int(win['no2'])]
        while left:
            r=tab['last'][left];fs.append(C.representative(('A'if r<=8 else'D',r)));left-=r
        fixedcost=sum(C.upper_base(t)for t in fixed)+sum(C.upper_cross(a,b)for i,a in enumerate(fixed)for b in fixed[i+1:])
        bound=fixedcost+win['value'];co=Counter(fs)
        assert sum(n for k,n in fs)==58
        assert C.value(fs,True)+40*co['A',1]*k5==bound and value(fs)<=bound
        row=dict(counts=cs,key=key,fixed_A5_count=k5,residual_rank=rem,
            weights=weights,tables=tables,context_options=opts,upper=bound,attaining_relaxed_forest=fs)
        if saved:assert enc(row)==old
        rows.append(row)
    maximum=max(r['upper']for r in rows)
    assert F(maximum,64)==F(2837,8)<377
    return dict(complete_domain_rows=len(domain),rows=rows,maximum=F(maximum,64),
        Bellman_inequalities=checks,attaining_predecessors=preds,item_identities=identities,context_options=options_count,
        every_potentially_increasing_forest_covered=True,extra_A1_rank5_weight='5/8',
        D17_decrease_ignored_in_relaxation=True,rank15_guard_relaxed_to_unconditional=True)

def refinements(docs,saved=None):
    old=docs[C.TAG];types=sorted(t for t in C.ENGINE.types()if t[1]<=58);fibers={}
    for t in types:fibers.setdefault(C.representative(t),[]).append(t)
    residual_reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS}
    residual=[t for t in types if C.representative(t)in residual_reps and t not in[('A',1),('A',2)]]
    assert not set(residual)&{('A',1),('A',2),('D',5),('A',9),('A',10),('A',11),('A',15),('D',15)}
    assert all(actual_cross(a,b)==0 for a in residual+[('A',1),('A',2)]for b in residual+[('A',1),('A',2)])
    oldrows={(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count']):r for r in saved['classes']}if saved else{}
    expected={(tuple(r['counts']),r['fixed_A5_count'],r['fixed_A12_count'])for r in old['finite']['above_curvature']}
    if saved:assert set(oldrows)==expected and len(oldrows)==len(saved['classes'])
    classes=[];checks=preds=options=0
    for u in old['finite']['above_curvature']:
        cs=tuple(u['counts']);k5=u['fixed_A5_count'];a12=u['fixed_A12_count'];key=(cs,k5,a12)
        groups=sorted(Counter(C.fixed_components(cs,u['key'])).items())
        if k5:groups.append((('A',5),k5))
        rem=58-sum(t[1]*n for t,n in groups)
        fixedsets=sorted(tuple(sorted(t for group in choice for t in group))for choice in product(*(combinations_with_replacement(fibers[t],n)for t,n in groups)))
        assert len(fixedsets)==len(set(fixedsets))
        prior={tuple(map(tuple,r['fixed_actual'])):r for r in oldrows[key]['rows']}if saved else{}
        if saved:assert set(prior)==set(fixedsets)and len(prior)==len(oldrows[key]['rows'])
        rows=[]
        for fixed in fixedsets:
            fixedcost=sum(actual_base(t)for t in fixed)+sum(actual_cross(a,b)for i,a in enumerate(fixed)for b in fixed[i+1:])
            ws=[actual_base(t)+sum(actual_cross(t,b)for b in fixed)for t in residual]
            if saved:v,last=prior[fixed]['values'],prior[fixed]['last']
            else:
                v=[0]+[None]*rem
                for t,w in zip(residual,ws):
                    for n in range(t[1],rem+1):
                        if v[n-t[1]]is not None:
                            z=v[n-t[1]]+w
                            if v[n]is None or z>v[n]:v[n]=z
                last=[None]+[next((i for i,(t,w)in enumerate(zip(residual,ws))if t[1]<=n and v[n-t[1]]is not None and v[n]==v[n-t[1]]+w),None)for n in range(1,rem+1)]
            assert len(v)==len(last)==rem+1 and v[0]==0 and last[0]is None
            for n in range(1,rem+1):
                possible=[]
                for i,(t,w)in enumerate(zip(residual,ws)):
                    if t[1]<=n and v[n-t[1]]is not None:
                        assert v[n]is not None and v[n]>=v[n-t[1]]+w;checks+=1;possible.append(i)
                if possible:
                    i=last[n];assert i in possible and v[n]==v[n-residual[i][1]]+ws[i];preds+=1
                else:assert v[n]is None and last[n]is None
            w1=actual_base(('A',1))+sum(actual_cross(('A',1),t)for t in fixed)
            w2=actual_base(('A',2))+sum(actual_cross(('A',2),t)for t in fixed)
            choices=[]
            for a1 in range(rem+1):
                for a2 in range((rem-a1)//2+1):
                    left=rem-a1-2*a2
                    if v[left]is None:continue
                    fs=list(fixed)+[('A',1)]*a1+[('A',2)]*a2
                    bonus=context(fs)
                    assert all(context(fs+[t])==bonus for t in residual)
                    choices.append((fixedcost+a1*w1+a2*w2+v[left]+bonus,a1,a2,left))
            options+=len(choices);best,a1,a2,left=max(choices)
            fs=list(fixed)+[('A',1)]*a1+[('A',2)]*a2
            while left:
                t=residual[last[left]];fs.append(t);left-=t[1]
            assert sum(n for k,n in fs)==58 and value(fs)==best
            row=dict(fixed_actual=fixed,residual_types=residual,weights=ws,values=v,last=last,
                sixtyfourths_actual_maximum=best,attaining_forest=sorted(fs),actual_A1_count=a1,actual_A2_count=a2)
            if saved:assert enc(row)==prior[fixed]
            rows.append(row)
        maximum=max(r['sixtyfourths_actual_maximum']for r in rows)
        classes.append(dict(counts=cs,fixed_A5_count=k5,fixed_A12_count=a12,residual_rank=rem,
            actual_maximum=F(maximum,64),rows=rows,actual_obstruction_retained=maximum>64*377,
            maximizers=[r['attaining_forest']for r in rows if r['sixtyfourths_actual_maximum']==maximum]))
    maximum=max(r['actual_maximum']for r in classes);assert maximum==379
    witness=[('A',16),('A',16),('D',26)];assert value(witness)==64*maximum
    return dict(classes=classes,complete_predecessor_upper_classes=len(classes),
        actual_fixed_assignments=sum(len(r['rows'])for r in classes),Bellman_inequalities=checks,
        attaining_predecessors=preds,explicit_A1_A2_options=options,
        exact_actual_maximum=maximum,actual_attaining_forest=witness,
        classes_still_above=sum(r['actual_obstruction_retained']for r in classes),
        complement_proof='Without an A9/A2 guarded increase, pointwise decrease from the complete predecessor; all potentially increasing forests are below377 by the separate complete conditional upper.',
        full_global_DP_repeated=False)

def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    geo=local(docs);inter=interfaces(docs)
    print('LOCAL',geo['affected_maps'],'maps',inter['later_corners'],'later corners',flush=True)
    upper=conditional_upper(saved['increasing_family_upper']if saved else None)
    print('INCREASING_FAMILY_UPPER',upper['maximum'],flush=True)
    finite=refinements(docs,saved['finite']if saved else None)
    ii=docs[C.TAG]['ideal'];assert ii['complete_profiles']==2552 and ii['all_positive']and ii['absolute_surplus']=='36'
    assert ii['weights_and_envelopes_unchanged']and ii['inherited_positive_low8_concentration_and_factor_fiber']
    excess=F(finite['exact_actual_maximum'])-377;assert excess==2
    minimum=-(-F(4205,excess).numerator//F(4205,excess).denominator);assert minimum==2103
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_D17_A9_A2_WITNESS_CONTEXT_ZERO_AMBIENT58_FRONTIER',
        sources={p:digest(p)for p in sorted(paths)},baseline=G.TAG,scale=64,
        local_proof=geo,interfaces=inter,increasing_family_upper=upper,finite=finite,
        rule='D17(2,12)=0; wholeA1/D5long=1 with wholeA9 and wholeA2; otherwise previous guarded3/8 rule',
        coupon_once_per_marked_pair=True,witness_multiplicity_not_counted=True,
        simultaneous_local_composition_proved=True,triangle_minimum='2',quadrilateral_minimum='1',
        ideal_weights_changed=False,positive_ideal_profiles=2552,actual_cusp_surplus='36',
        finite_curvature='377',finite_excess=str(excess),strict_global_constant='4205',
        necessary_ambient58_inequality='2*v_f >= 4205 + 36*v_i',minimum_finite_vertices=minimum,
        completed_dimension_bound=58,new_dimension_bound=None,new_finite_volume_example=False,
        scope='Complete selected domain and arbitrary-support radical seals, current chosen-output interfaces, fresh complete upper on759 potentially increasing representative classes, and fresh exact actual refinements of all39 predecessor above-curvature classes. The unchanged ideal theorem is bound, not replayed. No dimension57 theorem or new example.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists()
    saved=json.loads(p.read_text())if args.verify else None;out=run(saved)
    if args.verify:assert out==saved
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['exact_actual_maximum'],out['finite']['classes_still_above'],out['minimum_finite_vertices'],flush=True)
