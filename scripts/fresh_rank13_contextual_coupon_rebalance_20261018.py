#!/usr/bin/env python3
"""A contextual finite coupon repairs two rank13 reductions.

Exact local seals and complete rank59 forest certificates. The optional
A1/D13-neighbor increase is a counting preflight, not a geometric theorem.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_finite_rank14_existing_payment_transfer_20261018 as L
import fresh_ambient59_payment_frontier_20261018 as P
import fresh_order16_finite_cross_payments_20261017 as O
import fresh_cycle14_ideal5_transfer_20261017 as CY

ROOT=L.ROOT; T=L.T; A=L.A; Q=L.Q; D=L.D; C=P.C
TAG='fresh_rank13_contextual_coupon_rebalance_20261018'
DELTA=F(1,16); N=59; SCALE=16
STEMS=tuple(dict.fromkeys((*L.R.STEMS,L.TAG,L.L.TAG,P.TAG,P.COMPLETE,
    C.TAG,C.S.TAG,O.TAG,CY.TAG,*C.S.PROVENANCE,
    'fresh_order14_b_basic_distance9_transfer_20261018',
    'fresh_order14_c_basic_distance9_transfer_20261018')))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    return x


def key(stem,index):return stem,json.dumps(index,sort_keys=True)


def changes(g,white,marks):
    ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids))
    if sig!=[len(ids),0,0]:return dict(delta=F(),roles=[],components=[])
    bs=D.finite_blocks(g,ids);aa,bb=[next(b for b in bs if m in b['nodes']) for m in marks]
    roles=[];delta=F()
    if aa is not bb:
        for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
            if sm['type']=='A2' and la['type'] in ('A13','D13'):
                images=T.mark_pairs(g,la['nodes'],[lm],13,la['type'][0])
                wanted=[(6,)] if la['type']=='A13' else [(10,)]
                if images==wanted:
                    delta-=DELTA;roles.append('reduce_A2_'+la['type']+'_middle_or_center')
            if sm['type']=='A1' and la['type']=='A13':
                if T.mark_pairs(g,la['nodes'],[lm],13,'A')==[(6,)]:
                    delta+=DELTA;roles.append('increase_A1_A13_middle')
        types=Counter(b['type'] for b in bs)
        if aa['type']==bb['type']=='A2' and types['A10']>=1 and types['A1']>=2:
            delta+=2*DELTA;roles.append('context_A2_A2_with_A10_and_two_A1')
    return dict(delta=delta,roles=roles,components=bs)


def literal(docs):
    old={key(r['source'],r['case']):r for r in docs[L.TAG]['combined_literal_interface']['rows']}
    rows=[];hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=old[key(stem,index)];corners=[]
        assert len(prev['corners'])==len(cs)
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks']
            change=changes(g,white,marks);before=F(prev['corners'][j]['retained_payment_lower_bound'])
            after=before+change['delta'];assert after>=0
            hits.update(change['roles'])
            corners.append(dict(marks=marks,previous_lower=before,**change,retained_payment_lower_bound=after))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(rows)==216
    return dict(configurations=len(rows),changed_roles=dict(hits),rows=rows,
        minimum_total=min(r['total_lower_bound'] for r in rows),
        changed=[r for r in rows if any(c['delta'] for c in r['corners'])])


def local_seals(docs):
    # Fresh arithmetic for the full old order16 seal: positive affine kernel,
    # nonsingular Lorentz core and the inverse contact block (4,-4;-4,4).
    order16=O.exceptional();assert enc(order16)==docs[O.TAG]['exceptional_D15']
    g=A.matrix(order16['gram']);old=next(r for r in order16['cases'] if r['case']==384)
    white=sorted(set(range(len(g)))-set(old['black']))
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0] for b in old['black'])
    corners=[]
    for c in old['corners']:
        ch=changes(g,white,c['marks']);before=F(c['weight'])
        corners.append(dict(marks=c['marks'],previous=before,**ch,payment=before+ch['delta']))
    assert sorted(c['delta'] for c in corners)==[-DELTA,-DELTA,F(),2*DELTA]
    assert sum(c['payment'] for c in corners)==1
    context=next(c for c in corners if c['delta']==2*DELTA)
    assert Counter(b['type'] for b in context['components'])==dict(A1=2,A2=2,A10=1)
    out16=dict(gram=order16['gram'],white=white,black=old['black'],case=384,
        signature=order16['signature'],inverse_contact_block=order16['inverse_contact_block'],
        affine_nodes=order16['affine_nodes'],affine_positive_kernel=order16['affine_positive_kernel'],
        seal=order16['sealing'],corners=corners,total='1')

    # Reconstruct the cycle14 equality Gram and exact whole marked roles.
    g=[[F(2 if i==j else 0) for j in range(17)] for i in range(17)]
    for i,j in [(i,(i+1)%14) for i in range(14)]+[(0,14),(7,15),(14,16)]:g[i][j]=g[j][i]=F(-1)
    white=[i for i in range(14) if i not in (0,7)]+[16]
    assert Q.inertia(g)==[15,1,1] and white==docs[CY.TAG]['equality']['white']
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0] for b in (0,7,14,15))
    cycle=[]
    for old in docs[CY.TAG]['equality']['corners']:
        marks=old['marks'];ch=changes(g,white,marks);before=F(old['coupon'])
        cycle.append(dict(marks=marks,previous=before,**ch,payment=before+ch['delta']))
    assert sorted(c['delta'] for c in cycle)==[-DELTA,F(),F(),DELTA]
    assert sum(c['payment'] for c in cycle)==1

    cycle_gram=enc(g)
    # The two current D15 seals also have the contextual opposite corner.
    # Its two deltas more than restore the one lost center-coupon delta;
    # the existing selected same-component half stays unchanged.
    d15=[]
    for old in docs[L.L.TAG]['selected_proofs']:
        if old['index'] not in (1147,1269):continue
        g=A.matrix(old['gram']);white=old['white'];cs=[]
        assert Q.inertia(g)==[len(g)-1,1,0]
        inv=Q.inverse(g)
        assert sorted(inv[i][i] for i in old['remaining_contacts'])==[4,4,14]
        assert Q.inertia(A.sub(g,white))==[len(white),0,0]
        for c in old['corners']:
            ch=changes(g,white,c['marks']);before=F(c['payment'])
            cs.append(dict(marks=c['marks'],previous=before,**ch,payment=before+ch['delta']))
        assert sorted(c['delta'] for c in cs)==[-DELTA,F(),F(),2*DELTA]
        assert sum(c['payment'] for c in cs)==1+DELTA
        d15.append(dict(index=old['index'],old_seal=old,corners=cs,total=str(1+DELTA)))
    assert [r['index'] for r in d15]==[1147,1269]
    return dict(order16=out16,cycle14=dict(gram=cycle_gram,corners=cycle,total='1'),D15=d15)


def actual_base(t):
    value=2*C.actual_base(t)
    if t==('B',14):value-=64
    if t==('D',14):value-=8
    return value


def upper_base(t):return 2*C.upper_base(t)


def actual_cross(a,b,neighbor=False):
    value=2*C.actual_cross(a,b);pair={a,b}
    if pair in ({('A',2),('A',13)},{('A',2),('D',13)}):value-=2
    if pair=={('A',1),('A',13)}:value+=1
    if neighbor and pair=={('A',1),('D',13)}:value+=8
    return value


def upper_cross(a,b,neighbor=False):
    value=2*C.upper_cross(a,b)
    if {a,b}=={('A',2),('D',13)}:value-=2
    if {a,b}=={('A',1),('D',13)}:value+=8 if neighbor else 1
    return value


def context_cost(fs,upper=False):
    counts=Counter(fs);has10=any(r==10 for k,r in fs) if upper else counts['A',10]>0
    k=counts['A',2]
    return 4*k*(k-1) if has10 and counts['A',1]>=2 else 0


def value(fs,upper=False,neighbor=False):
    base,cross=(upper_base,upper_cross) if upper else (actual_base,actual_cross)
    return sum(base(t) for t in fs)+sum(cross(a,b,neighbor) for i,a in enumerate(fs) for b in fs[i+1:])+context_cost(fs,upper)


def item(r,key,neighbor=False):
    return 2*C.item(r,key)+((8 if neighbor else 1)*key[2] if r==1 else -2*key[2] if r==2 else 0)


def table(key,exclude2,neighbor):
    weights={r:item(r,key,neighbor) for r in P.ITEMS if not exclude2 or r!=2}
    v=[0]+[-10**9]*N
    for r,w in weights.items():
        for n in range(r,N+1):v[n]=max(v[n],v[n-r]+w)
    last=[0]+[next(r for r,w in weights.items() if r<=n and v[n]==v[n-r]+w) for n in range(1,N+1)]
    return dict(key=list(key),exclude2=exclude2,values=v,last=last)


def check_table(t,key,exclude2,neighbor):
    assert t['key']==list(key) and t['exclude2']==exclude2
    v,last=t['values'],t['last'];assert len(v)==len(last)==N+1 and v[0]==last[0]==0
    ws={r:item(r,key,neighbor) for r in P.ITEMS if not exclude2 or r!=2};checks=0
    for n in range(1,N+1):
        for r,w in ws.items():
            if r<=n:assert v[n]>=v[n-r]+w;checks+=1
        r=last[n];assert r in ws and r<=n and v[n]==v[n-r]+ws[r]
    return checks


def finite(neighbor=False,saved=None):
    rep=P.Z.B.B.representative;types=[t for t in P.Z.B.B.types() if t[1]<=N]
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0
        for b in types:assert upper_cross(rep(a),rep(b),neighbor)>=actual_cross(a,b,neighbor)==actual_cross(b,a,neighbor)>=0
    dom=P.fixed_domain();keys=sorted({key for cs,used,key in dom});assert len(dom)==961 and len(keys)==960
    tables={(tuple(t['key']),t['exclude2']):t for t in saved['tables']} if saved else {
        (k,e):table(k,e,neighbor) for k in keys for e in (False,True)}
    assert set(tables)=={(k,e) for k in keys for e in (False,True)} and len(tables)==1920
    if saved:assert len(saved['tables'])==1920
    checks=sum(check_table(tables[k,e],k,e,neighbor) for k in keys for e in (False,True))
    residual=[rep(('A' if r<=8 else 'D',r)) for r in P.ITEMS]
    assert all(upper_cross(a,b,neighbor)==0 for a in residual for b in residual)
    rows=[];items=choices=0
    for cs,used,k in dom:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for _,r in fixed)==used and k[2]==cs[4] and k[6]==cs[1]
        assert all(upper_cross(a,b,neighbor)==0 for a in fixed for b in fixed)
        for r,t in zip(P.ITEMS,residual):
            assert item(r,k,neighbor)==upper_base(t)+sum(upper_cross(t,b,neighbor) for b in fixed);items+=1
        rem=N-used;normal=tables[k,False];no2=tables[k,True]
        options=[(normal['values'][rem],0,0,False)]
        if cs[1] and rem>=2:
            options += [(2*item(1,k,neighbor)+j*item(2,k,neighbor)+no2['values'][rem-2-2*j]+4*j*(j-1),j,2,True)
                for j in range((rem-2)//2+1)]
        choices+=len(options);best,j,a,e=max(options);tab=tables[k,e]
        fs=fixed+[('A',2)]*j+[('A',1)]*a;rest=rem-2*j-a
        while rest:
            r=tab['last'][rest];fs.append(rep(('A' if r<=8 else 'D',r)));rest-=r
        upper=sum(upper_base(t) for t in fixed)+best
        # A no-bonus attaining forest can have the trigger; then the exact
        # positive-bonus option dominates it. Equality checks this explicitly.
        assert sum(r for _,r in fs)==N and value(fs,True,neighbor)==upper and value(fs,False,neighbor)<=upper
        rows.append(dict(counts=list(cs),key=list(k),forest=sorted(x+str(r) for x,r in fs),
            sixteenths_upper=upper,sixteenths_actual=value(fs,False,neighbor),triggered=e,
            fixed_A2_count=j,reserved_A1_count=a))
    maximum=max(r['sixteenths_upper'] for r in rows)
    assert maximum==max(r['sixteenths_actual'] for r in rows)==6360
    attainers=[['D16','D17','D26'],['A17','D16','D26'],['D16','D16','D27']]
    assert all(value([(s[0],int(s[1:])) for s in fs],False,neighbor)==maximum for fs in attainers)
    witnesses=[]
    for large in ('A','D'):
        fs=[('A',2)]*10+[(large,13)]*3
        assert value(fs,False,neighbor)==6244
        witnesses.append(dict(forest=sorted(k+str(r) for k,r in fs),old='394',new=str(F(6244,16))))
    assert max(r['sixteenths_upper'] for r in rows if r['counts'][4])<6264
    return dict(neighbor_half_increase_is_conditional_preflight=neighbor,exact_maximum=str(F(maximum,16)),
        curvature='783/2',fixed_vectors=961,interaction_keys=960,tables=[tables[k,e] for k in keys for e in (False,True)],
        simultaneous_type_pairs=len(types)**2,residual_item_checks=items,Bellman_upper_inequalities=checks,
        attaining_predecessors=2*len(keys)*N,context_options=choices,alternatives=rows,
        above_curvature=[r for r in rows if r['sixteenths_upper']>6264],actual_attainers=attainers,
        max_with_rank13=str(F(max(r['sixteenths_upper'] for r in rows if r['counts'][4]),16)),
        max_with_rank10=str(F(max(r['sixteenths_upper'] for r in rows if r['counts'][1]),16)),
        small_obstructions_removed=witnesses,actual_D15_cost=str(F(actual_base(('D',15)),16)),
        A2_A13_component_cost=str(F(actual_cross(('A',2),('A',13),neighbor),16)),
        A2_D13_component_cost=str(F(actual_cross(('A',2),('D',13),neighbor),16)),
        A1_A13_component_cost=str(F(actual_cross(('A',1),('A',13),neighbor),16)),
        A1_D13_component_cost=str(F(actual_cross(('A',1),('D',13),neighbor),16)),
        context_pair_payment='1/8',context_total='k*(k-1)/4 if whole A10 and at least two whole A1 are present; otherwise0',
        complete_context_count=True,geometric_realization_of_formal_forests_claimed=False)


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[P.COMPLETE]['dimension_bound']==59
    paths={Path(__file__),Path(L.__file__),Path(P.__file__),Path(C.__file__),Path(O.__file__),Path(CY.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    local=local_seals(docs);lit=literal(docs)
    return dict(status='PASS_RANK13_CONTEXTUAL_COUPON_REBALANCE',sources=sources,delta=str(DELTA),
        local_seals=local,combined_literal_interface=lit,
        finite=finite(False,saved['finite'] if saved else None),
        conditional_neighbor_preflight=finite(True,saved['conditional_neighbor_preflight'] if saved else None),
        triangle_minimum=2,quadrilateral_minimum=1,ideal_weights_unchanged=True,
        completed_unrestricted_bound=59,new_dimension_bound=None,new_finite_volume_example=False,
        scope='Exact reconstructed seals and all216 literal interfaces, plus complete rank59 conditional-component forest certificates. Whole-face composition and complete coupon-output inventory are written in the accompanying report. Optional A1/D13-neighbor half is only a finite-cost preflight; no geometry for it is claimed here.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    if not args.verify:assert not path.exists(),'Frozen packet exists; use --verify.'
    out=enc(run(saved))
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'rank13',out['finite']['max_with_rank13'],
        'literal',out['combined_literal_interface']['changed_roles'],flush=True)
