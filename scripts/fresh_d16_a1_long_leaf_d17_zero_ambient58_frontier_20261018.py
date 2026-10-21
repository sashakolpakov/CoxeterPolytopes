#!/usr/bin/env python3
"""Whole-A1 D16 long-leaf units replace the restored D17 eighths.

Retain an unconditional D16 half for older ideal obligations. Four complete
radical seals supply the A1 guard. One new residual table and inherited
complete conditional/fiber bounds prove the sharp ambient58 maximum378.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import fresh_final_rank14_rank16_composed_ambient58_frontier_20261018 as Z
import fresh_d16_d17_eighth_balance_20261018 as E

ROOT=Z.ROOT;R=Z.R;N=R.N;C=R.C;L=R.L;B=R.B;W=R.W
J=E.J;T=R.T;A=R.A;Q=R.Q;D=R.D
TAG='fresh_d16_a1_long_leaf_d17_zero_ambient58_frontier_20261018'
UP=E.UP;DOWN=E.DOWN;CASES=[1644,1664,1931,2324]
BALANCE='fresh_center_reversion_rank16_balance_ambient60_frontier_20261018'
STEMS=tuple(dict.fromkeys((Z.TAG,*Z.STEMS,*R.STEMS,E.TAG,J.TAG,E.SUPPORT,BALANCE)))
enc=R.enc

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return Z.actual_base(t)-{('D',16):48,('D',17):16}.get(t,0)
actual_cross=Z.actual_cross
def context(fs):return Z.context(fs)+(64*fs.count(('D',16)) if ('A',1)in fs else 0)
def value(fs):return sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)

def role(g,blocks,marks):
    b=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
    if b is None or b['type']not in('D16','D17'):return None
    typ=b['type'];images=T.mark_pairs(g,b['nodes'],marks,len(b['nodes']),'D')
    changes=UP if typ=='D16'else DOWN
    if not set(images)&changes:return None
    assert set(images)<=changes
    witnesses=[b for b in blocks if b['type']=='A1'and not set(marks)&set(b['nodes'])]
    old=F(7,8)if typ=='D16'else F(1,8)
    new=(F(1)if witnesses else F(1,2))if typ=='D16'else F()
    return dict(type=typ,images=images,whole_A1_witnesses=witnesses,old=old,new=new,delta=new-old)

def local(docs):
    over=docs[A.OVER];j=docs[J.TAG];sup=j['complete_selected_support']
    target=E.P.labelled(16,[(i,i+1,3)for i in range(13)]+[(13,14,3),(13,15,3)])
    assert J.V.validate_maps(sup,target,over)==578 and len(sup['rows'])==394
    affected=[r['index']for r in sup['rows']if set(map(tuple,r['marked_images']))&UP]
    assert affected==CASES
    d17=docs[E.SUPPORT];rows17=d17['complete_selected_support']['rows'];assert len(rows17)==140
    expected17=[c['index']for c in over['isolation']['cases']if(c.get('component_order_bound')or 0)>=17 and len(c['marked_component'])<=17]
    assert [r['index']for r in rows17]==expected17;maps17=0
    for row in rows17:
        c=over['isolation']['cases'][row['index']];g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        nodes=c['marked_component'];labels=T.T.labels(g,nodes);seen=set();pairs=set()
        for mp in row['maps']:
            images=mp['images'];assert len(images)==len(nodes)==len(set(images))and tuple(images)not in seen
            seen.add(tuple(images));assert all(0<=i<17 for i in images)
            assert all(labels[i][j]==d17['target'][images[i]][images[j]]for i in range(len(nodes))for j in range(len(nodes)))
            pair=tuple(sorted(images[nodes.index(m)]for m in c['preferred_corner']))
            assert list(pair)==mp['marked_images'];pairs.add(pair);maps17+=1
        assert sorted(pairs)==list(map(tuple,row['marked_images']))and not pairs&DOWN
    assert maps17==168
    old={p['index']:p for p in j['selected_proofs']};proofs=[]
    for index in CASES:
        c=over['isolation']['cases'][index];core=over['isolation']['cores'][c['core_key']]
        p=J.seal(c,core);assert J.enc(p)==old[index]and p['sealing']['method']=='radical'
        rows=[]
        for corner in p['corners']:
            rr=role(A.matrix(p['gram']),corner['components'],corner['marks'])
            if rr:
                if rr['type']=='D16':assert len(rr['whole_A1_witnesses'])==1 and rr['new']==1
                payment=rr['new']
            else:assert F(corner['payment'])==0;payment=F()
            rows.append(dict(marks=corner['marks'],components=corner['components'],role=rr,payment=payment))
        assert sorted(r['payment']for r in rows)==[0,0,0,1]
        proofs.append(dict(index=index,old_complete_seal=p,corners=rows,total='1'))
    old16={tuple(r['pair']):F(r['weight'])for r in docs[R.TAG]['local_proof']['ordinary_weighted_tables']['D16']}
    assert all(old16[p]==F(7,8)for p in UP)
    for p in UP:old16[p]=F(1,2)
    assert sum(old16.values())==F(actual_base(('D',16)),64)==F(217,2)
    old17={tuple(r['pair']):F(r['weight'])for r in docs[BALANCE]['local_proofs']['finite_weighted_pairs']['D17']}
    assert all(old17[p]==F(1,8)for p in DOWN)
    return dict(complete_D16_patterns=394,validated_D16_maps=578,complete_D17_patterns=140,validated_D17_maps=168,
        selected_D16_contexts=CASES,selected_D17_changed_maps=0,whole_proofs=proofs,
        D16_ordinary_weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(old16.items())],
        D16_unguarded_weight='1/2',D16_guarded_weight='1',D17_weight='0',
        ordinary_D16_cost='217/2',guarded_D16_cost='219/2',ordinary_D17_cost='114')

def interfaces(docs):
    literals=[];hits=[]
    previous=docs[R.TAG]['current_interfaces']['literals']
    for old,(stem,case,g,white,cs)in zip(previous,T.literal_inputs(docs),strict=True):
        assert(old['source'],old['case'])==(stem,enc(case));rows=[]
        for p,(marks,*_)in zip(old['corners'],cs,strict=True):
            assert p['marks']==marks
            rr=role(g,p['components'],marks);delta=rr['delta']if rr else F()
            payment=F(p['lower'])+delta;assert payment>=0
            if rr:
                assert rr['type']=='D16'and not rr['whole_A1_witnesses']and payment==F(1,2)
                hits.append(dict(source=stem,case=case,marks=marks,role=rr))
            rows.append(dict(marks=marks,role=rr,previous=p['lower'],lower=payment))
        total=sum(r['lower']for r in rows);assert total>=1
        literals.append(dict(source=stem,case=case,corners=rows,total=total))
    assert len(literals)==216 and len(hits)==4
    assert {(r['source'],r['case'])for r in hits}=={('fresh_order16_finite_cross_payments_20261017',i)for i in(663,676,829)}
    later=[];affected=[]
    for stem,keys in B.FINITE.items():
        rows=[]
        for key in keys:
            for path,p,gg,w,blocks in B.corners(docs[stem][key],key):
                rr=role(A.matrix(gg),blocks,p['marks'])
                if rr:
                    assert stem in(J.TAG,E.TAG)
                    assert rr['type']=='D17'or len(rr['whole_A1_witnesses'])==1
                    affected.append(dict(source=stem,path=path,role=rr))
                rows.append(dict(path=path,role=rr))
        later.append(dict(source=stem,corners=rows))
    assert len(later)==22 and sum(len(r['corners'])for r in later)==564 and len(affected)==24
    newer=[]
    for stem,key in((R.TAG,'local_proof'),(N.TAG,'local_proof'),(Z.U.TAG,'local_proof'),(R.R.TAG,'local_seals')):
        rows=[]
        for path,p,gg,w,blocks in B.corners(docs[stem][key],key):
            assert role(A.matrix(gg),blocks,p['marks'])is None
            rows.append(path)
        assert rows;newer.append(dict(source=stem,unchanged_corners=rows))
    variants=[]
    for p in docs[Z.U.S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            assert all(role(A.matrix(v['gram']),c['components'],c['marks'])is None for c in v['corners'])
            variants.append([p['index'],v['extra_contact']])
    assert len(variants)==24
    d12=[]
    for stem in(B.PART,B.REST):
        rows=[(path,p)for path,p in W.walk(docs[stem])if F(p['payment'])>0]
        assert all(not any(b['type']in('D16','D17')for b in p['roles'])for path,p in rows)
        d12.append(dict(source=stem,positive_outputs=len(rows)))
    assert sum(r['positive_outputs']for r in d12)==332
    return dict(literals=literals,changed_literal_corners=hits,later_packets=later,
        repeated_four_exchange_outputs=affected,newer_seals=newer,unchanged_D14_variants=variants,chosen_D12_outputs=d12)

def finite(docs,saved=None):
    # The old conditional relaxation already treats every rank16 as110.25
    # with A1. New D16 with A1 is only109.5; all cross costs are unchanged.
    conditional=docs[R.TAG]['increasing_family_upper'];rows=conditional['rows']
    assert len(rows)==539 and conditional['maximum']=='757/2'
    exceptional=[r for r in rows if F(r['upper'])>377];assert len(exceptional)==1
    e=exceptional[0];assert e['rank16_count']==2 and e['fixed_A5_count']==0
    gaps=[]
    for kind in'ABD':
        t=(kind,16);guarded=actual_base(t)+(128 if kind=='A'else 64 if kind=='D'else 0)
        gap=C.upper_base(('D',16))+48-guarded+C.upper_cross(('A',1),('D',16))-actual_cross(('A',1),t)
        assert gap>=32;gaps.append(dict(type=t,guarded_base=F(guarded,64),one_A1_gap=F(gap,64)))
    assert F(e['upper'])-2*min(r['one_A1_gap']for r in gaps)==F(755,2)
    # The old D17 increasing family stays below curvature even if every
    # possible D16 contributes the new maximum increase1/4.
    overlap=F(conditional['D17_family_upper_after_trade'])+F(1,4)*(58//16)
    assert overlap<377
    classes=[];count=0
    for cl in docs[N.TAG]['finite']['classes']:
        out=[]
        for r in cl['rows']:
            fs=list(map(tuple,r['fixed_actual']));co=Counter(fs)
            drop=80*co['A',16]+64*co['D',16]+(64 if cl['residual_rank']==0 else 32)*co['A',14]
            if cl['residual_rank']==0:drop+=64*co['D',14]
            upper=F(r['sixtyfourths_actual_maximum']-drop,64);count+=1
            out.append(dict(fixed_actual=fs,guaranteed_drop=F(drop,64),upper=upper))
        classes.append(dict(counts=cl['counts'],fixed_A5_count=cl['fixed_A5_count'],fixed_A12_count=cl['fixed_A12_count'],
            residual_rank=cl['residual_rank'],actual_fixed_fibers=out,upper=max(r['upper']for r in out)))
    assert len(classes)==39 and count==748
    costly=[r for r in classes if r['upper']>378]
    assert len(costly)==1 and not any(costly[0]['counts'])and costly[0]['fixed_A5_count']==0
    # Only this all-residual class needs fresh optimization. It has no fixed
    # context witness or cross interaction, including explicit A1/A2.
    inherited=next(cl for cl in docs[N.TAG]['finite']['classes']if not any(cl['counts'])and cl['fixed_A5_count']==0)
    assert inherited['residual_rank']==58 and len(inherited['rows'])==1
    prior=inherited['rows'][0];assert prior['fixed_actual']==[]
    types=[('A',1),('A',2)]+list(map(tuple,prior['residual_types']))
    assert len(types)==len(set(types))and all(actual_cross(a,b)==0 for a in types for b in types)
    assert not set(types)&{('A',9),('A',10),('A',11),('D',5),('A',14),('D',14),('A',15),('D',15),('A',16),('D',16)}
    assert all(context([a])==0 for a in types)
    weights=[actual_base(t)for t in types]
    if saved:v,last=saved['residual_table']['values'],saved['residual_table']['last']
    else:
        v=[0]+[None]*58
        for t,w in zip(types,weights):
            for n in range(t[1],59):
                if v[n-t[1]]is not None:v[n]=max(v[n]if v[n]is not None else-10**9,v[n-t[1]]+w)
        last=[None]+[next(i for i,(t,w)in enumerate(zip(types,weights))if t[1]<=n and v[n-t[1]]is not None and v[n]==v[n-t[1]]+w)for n in range(1,59)]
    checks=0
    assert len(v)==len(last)==59 and v[0]==0 and last[0]is None
    for n in range(1,59):
        for t,w in zip(types,weights):
            if t[1]<=n:assert v[n]>=v[n-t[1]]+w;checks+=1
        i=last[n];assert 0<=i<len(types)and types[i][1]<=n and v[n]==v[n-types[i][1]]+weights[i]
    assert F(v[58],64)==378
    witness=[('D',17),('D',17),('D',24)];plateau=[('A',13)]+[('A',15)]*3
    assert value(witness)==value(plateau)==64*378 and context(witness)==context(plateau)==0
    return dict(conditional_rank16_A1_classes=539,conditional_rank16_A1_upper='755/2',one_A1_type_gaps=gaps,
        inherited_D17_increasing_family_upper=overlap,complement_classes=classes,complete_actual_fixed_assignments=748,
        residual_table=dict(types=types,weights=weights,values=v,last=last,Bellman_inequalities=checks,attaining_predecessors=58),
        exact_actual_maximum='378',actual_attainers=[witness,plateau],new_obstruction_census=False,new_global_DP=False)

def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    proof=local(docs);inter=interfaces(docs);count=finite(docs,None if saved is None else saved['finite'])
    assert F(actual_base(('D',17)),64)==114
    for t in C.TYPES:
        assert actual_base(t)>=0
        for u in C.TYPES:assert actual_cross(t,u)==Z.actual_cross(t,u)
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_D16_A1_LONG_LEAF_D17_ZERO_AMBIENT58_FRONTIER',sources={p:digest(p)for p in sorted(paths)},
        baseline=Z.TAG,scale=64,local_proof=proof,current_interfaces=inter,finite=count,
        cross_costs_changed=False,new_cross_coupons=[],ideal_weights_changed=False,
        positive_ideal_profiles=2552,actual_cusp_surplus='36',finite_curvature='377',strict_global_constant='4205',
        triangle_minimum='2',quadrilateral_minimum='1',necessary_ambient58_inequality='v_f >= 4205 + 36*v_i',
        minimum_finite_vertices=4205,completed_dimension_bound=58,new_dimension_bound=None,
        scope='Four reconstructed full radical exchanges and complete current output retention, inherited complete conditional/fiber bounds, and one fresh all-residual table. Sharp formal378, not a dimension57 theorem or new obstruction census.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists()
    saved=json.loads(p.read_text())if args.verify else None;out=run(saved)
    if args.verify:assert out==saved
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['exact_actual_maximum'],flush=True)
