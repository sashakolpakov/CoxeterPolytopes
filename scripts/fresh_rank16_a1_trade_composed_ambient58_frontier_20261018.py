#!/usr/bin/env python3
"""Whole-A1 guarded A16 unit, D16 zero, and a stronger A4/A12 third coupon.

Composes the frozen D17 and D14 rules. Exact global finite maximum by complete
conditional and pointwise certificates; no new obstruction-class census.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import fresh_d17_a9_a2_context_zero_ambient58_frontier_20261018 as N
import fresh_d14_small_component_context_retention_20261018 as K

ROOT=N.ROOT;C=N.C;G=N.G;L=N.L;W=N.W;B=N.B;R=N.R;T=N.T;A=N.A;Q=N.Q;D=N.D
TAG='fresh_rank16_a1_trade_composed_ambient58_frontier_20261018'
STEMS=tuple(dict.fromkeys((N.TAG,K.TAG,*N.STEMS,*K.STEMS)))
PINS={N.TAG:('cf43d27af47f9f23216b89f6ef76fe35ccb5973e33a1f4216ff12fe60293cdac','f80edecac4af2021e39806a194befebcba0ab83c748138a6e40556a7dedc117c','d73be7ab2e27d029674cd236f3806116c400ecf2c9f478aa0dedd26fbcc0fcde'),
K.TAG:('0d887f3e8358522b024d13298cb30dc5ad10ce522e0fa41ab6c94231145a3de3','fb484aa601a397838fa63ac61c48c2a8d5f7fdc32c5d1de29505919126fee0e5','10fb43eb7646aeb36b1f14c776d3defc76a6f63923d8d96c2f29bbbe01170f02')}
enc=N.enc

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return N.actual_base(t)-{('A',16):80,('D',16):16,('D',14):64}.get(t,0)
def actual_cross(a,b):return N.actual_cross(a,b)-(32 if L.pair(a,b)else 0)
def context(fs):
    co=Counter(fs)
    return N.context(fs)+(128*co['A',16]if co['A',1]else 0)+(64*co['D',14]if any(co['A',r]for r in(1,2,3))else 0)
def value(fs):return sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)

def changes(g,white,marks):
    ids=sorted(set(white)|set(marks))
    if Q.inertia(A.sub(g,ids))!=[len(ids),0,0]:return dict(delta=F(),roles=[],components=[])
    blocks=D.finite_blocks(g,ids);aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
    delta=F();roles=[]
    if aa is bb and aa['type']in('A16','D16'):
        kind=aa['type'][0];im=T.mark_pairs(g,aa['nodes'],marks,16,kind)
        changed=L.ACHANGE if kind=='A'else L.DCHANGE
        if set(im)&changed:
            assert set(im)<=changed
            ws=[b for b in blocks if b['type']=='A1'and not set(marks)&set(b['nodes'])]
            new=F(1)if kind=='A'and ws else F()
            old=F(5,8)if kind=='A'else F(1,4);delta=new-old
            roles.append(dict(role=aa['type']+'_trade',images=im,whole_A1_witnesses=ws,old=old,new=new))
    elif aa is not bb:
        for sm,la,m,lm in((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
            if(sm['type'],la['type'])!=('A4','A12'):continue
            si=T.mark_pairs(g,sm['nodes'],[m],4,'A');li=T.mark_pairs(g,la['nodes'],[lm],12,'A')
            if si!=[(0,),(3,)]:continue
            if li==[(1,),(10,)]:old,new,role=F(3,8),F(),'A4end_A12second_removed'
            elif li==[(2,),(9,)]:old,new,role=F(3,4),F(1),'A4end_A12third_unit'
            else:continue
            delta=new-old;roles.append(dict(role=role,small_images=si,large_images=li,old=old,new=new))
    assert len(roles)<=1
    return dict(delta=delta,roles=roles,components=blocks)

def selected(docs):
    old=docs[L.TAG]['local_proofs'];assert len(old['whole_proofs'])==6
    assert [p['index']for p in old['proper_completion_proofs']]==L.DPROPER
    for kind,maps in(('A',238),('D',578)):
        sup=old['domains'][kind]['complete_selected_support']
        assert len(sup['rows'])==394 and sum(len(r['maps'])for r in sup['rows'])==maps
    rows=[]
    for p in old['whole_proofs']:
        g=A.matrix(p['gram']);cs=[]
        for c in p['corners']:
            ch=changes(g,p['white'],c['marks']);new=F(c['payment'])+ch['delta'];assert new>=0
            if p['index']in L.ACASES and ch['roles']and ch['roles'][0]['role']=='A16_trade':
                assert len(ch['roles'][0]['whole_A1_witnesses'])==1 and new==1
            cs.append(dict(marks=c['marks'],previous=c['payment'],**ch,payment=new))
        assert sorted(c['payment']for c in cs)==[0,0,0,1]
        assert sum(c['delta']for c in cs)==0
        rows.append(dict(index=p['index'],gram=p['gram'],white=p['white'],corners=cs,total='1',
            inherited_arbitrary_support_seal=L.TAG))
    assert [r['index']for r in rows]==L.ACASES+L.DWHOLE
    tables={}
    for typ,changed in(('A16',L.ACHANGE),('D16',L.DCHANGE)):
        before={tuple(r['pair']):F(r['weight'])for r in docs[L.TAG]['weighted_tables'][typ]}
        assert sum(before.values())==F(N.actual_base((typ[0],16)),64)==F(219,2)
        zero={p:F()if p in changed else w for p,w in before.items()}
        assert sum(zero.values())==F(actual_base((typ[0],16)),64)
        tables[typ]=[dict(pair=p,weight=w)for p,w in sorted(zero.items())]
    assert F(actual_base(('A',16)),64)==F(433,4)and F(actual_base(('D',16)),64)==F(437,4)
    return dict(complete_A_patterns=394,complete_A_maps=238,complete_D_patterns=394,complete_D_maps=578,
        whole_proofs=rows,inherited_proper_exclusions=L.DPROPER,ordinary_weighted_tables=tables,
        guarded_A16_pairs=sorted(L.ACHANGE),guarded_A16_weight='1',unguarded_A16_weight='0',D16_weight='0',
        A16_ordinary_cost='433/4',A16_cost_with_A1='441/4',D16_cost='437/4',
        A4_A12_second_weight='0',A4_A12_third_weight='1',A4_A12_total='8',
        guard_counts_whole_A1_presence_not_multiplicity=True)

def interfaces(docs):
    literals=[];changed=[];removed_unused=[]
    for prev,(stem,case,g,white,cs)in zip(docs[L.TAG]['current_interfaces']['literal_rows'],T.literal_inputs(docs),strict=True):
        assert(prev['source'],prev['case'])==(stem,enc(case));corners=[]
        for p,(marks,*_)in zip(prev['corners'],cs,strict=True):
            assert p['marks']==marks
            ch=changes(g,white,marks);lower=F(p['lower'])+ch['delta'];assert lower>=0
            # The later D17 and D14 packets separately prove these literal
            # outputs nondecreasing/retained. Their exact arrays are bound.
            corners.append(dict(marks=marks,previous=p['lower'],**ch,lower=lower))
        total=sum(c['lower']for c in corners);assert total>=1
        for oldc,newc in zip(prev['corners'],corners,strict=True):
            if newc['roles']and newc['roles'][0]['role']=='A4end_A12second_removed':
                assert F(oldc['before'])==newc['lower']==0
                removed_unused.append([stem,case,newc['marks']])
        if any(c['roles']and c['roles'][0]['role']=='A16_trade'and c['delta']<0 for c in corners):
            assert stem==L.LEAF and case in(6,7,12,13,18,19,24,25)
            u=prev['independent_retained_unit'];assert u['weight']=='1'
            assert next(c for c in corners if c['marks']==u['marks'])['lower']==1
            changed.append([stem,case])
        literals.append(dict(source=stem,case=case,corners=corners,total=total))
    assert len(literals)==216 and len(changed)==8 and len(removed_unused)==8
    assert len(docs[N.TAG]['interfaces']['literals'])==216
    for a,b in zip(literals,docs[N.TAG]['interfaces']['literals'],strict=True):
        assert(a['source'],enc(a['case']))==(b['source'],b['case'])
        assert [c['marks']for c in a['corners']]==[c['marks']for c in b['corners']]
        assert all(F(c['delta'])>=0 for c in b['corners'])
    later=[];debits=[];unused_later=[];total=0
    for stem,keys in B.FINITE.items():
        rows=[]
        for section in keys:
            for path,p,gg,white,blocks in B.corners(docs[stem][section],section):
                ch=changes(A.matrix(gg),white,p['marks'])
                if ch['delta']<0:
                    if ch['roles'][0]['role']=='A4end_A12second_removed':
                        assert F(p['payment'])==0
                        unused_later.append(dict(source=stem,path=path,old_required_payment='0',new_payment='0'))
                    else:
                        assert stem==L.H.TAG and ch['roles'][0]['role']=='D16_trade'
                        debits.append(dict(source=stem,path=path,delta=ch['delta']))
                rows.append(dict(path=path,marks=p['marks'],delta=ch['delta'],roles=ch['roles']));total+=1
        later.append(dict(source=stem,corners=rows))
    assert len(later)==22 and total==564 and len(debits)==4
    variants=[]
    current={(p['index'],p['extra_contact']):p for p in docs[G.TAG]['local_proof']['whole_variants']}
    for p in docs[L.C.R.S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            now=current[p['index'],v['extra_contact']];g=A.matrix(v['gram']);rows=[]
            for c,paid in zip(v['corners'],now['corners'],strict=True):
                assert c['marks']==paid['marks'];ch=changes(g,v['white'],c['marks'])
                assert ch['delta']==0 and F(paid['rank16_transfer_increment'])==0
                rows.append(dict(marks=c['marks'],delta=ch['delta'],lower=paid['joint_lower']))
            assert sum(F(c['lower'])for c in rows)>=1
            variants.append(dict(index=p['index'],extra_contact=v['extra_contact'],corners=rows))
    assert len(variants)==24
    seals=[]
    for p in docs[R.TAG]['local_seals']:
        g=A.matrix(p['gram']);rows=[]
        for c in p['corners']:
            ch=changes(g,p['white'],c['marks']);assert ch['delta']==0
            rows.append(dict(marks=c['marks'],delta=ch['delta']))
        seals.append(dict(name=p['name'],corners=rows))
    assert len(seals)==4
    d12=[]
    for stem,expected in((B.PART,280),(B.REST,52)):
        ps=[(path,p)for path,p in W.walk(docs[stem])if F(p['payment'])>0]
        assert len(ps)==expected
        for path,p in ps:
            types={b['type']for b in p['roles']}
            assert not types&{'A16','D16'}and not {'A4','A12'}<=types
        d12.append(dict(source=stem,positive_outputs=len(ps),no_new_trade_role=True))
    additional=[]
    for p in docs[N.TAG]['local_proof']['proofs']:
        rows=[]
        for c in p['corners']:
            ch=changes(A.matrix(p['gram']),p['white'],c['marks']);assert ch['delta']==0
            rows.append(dict(marks=c['marks'],delta=ch['delta']))
        additional.append(dict(index=p['index'],corners=rows))
    assert len(additional)==2
    assert len(docs[K.TAG]['selected_whole_variants'])==24
    for p in docs[K.TAG]['coupled_order11_boundary']['canonical_cases']:
        assert changes(A.matrix(p['gram']),p['white'],p['marks'])['delta']==0
    return dict(literal_configurations=216,literals=literals,changed_literal_rows=changed,
        later_packets=later,later_corners=total,replaced_own_D16_outputs=debits,unused_later_coupon_appearances=unused_later,
        retained_quarter_variants=variants,retained_rank13_seals=seals,chosen_D12_outputs=d12,
        unchanged_new_D17_seals=additional,unchanged_D14_order11_boundary_cases=4,
        removed_unused_literal_increments=removed_unused,no_other_required_second_coupon_consumer=True)

def conditional_upper(docs):
    """All rank16+A1 forests under C, globally adding3/4 per rank16.

    The separate D17 increasing family is bounded below curvature even after
    the largest possible new A16 increase, so its additional coupon need not
    be installed in these unchanged residual tables.
    """
    tabs={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in docs[C.TAG]['finite']['tables']}
    rows=[];options=0
    for cs,used,key in C.fixed_domain():
        large=C.fixed_components(cs,key);c16=large.count(('D',16))
        if not c16:continue
        for k5 in range((58-used)//5+1):
            fixed=large+[('A',5)]*k5;rem=58-used-5*k5
            if rem<1:continue
            normal=tabs[key,k5,False];no2=tabs[key,k5,True]
            opts=[(C.item(1,key,k5)+normal['values'][rem-1],1,0,False)]
            for reserve in(1,2):
                oldguard=reserve==2 and bool(cs[1]);newguard=bool(cs[0]and k5)
                if rem<reserve or not(oldguard or newguard)or(reserve==2 and not oldguard):continue
                for j in range((rem-reserve)//2+1):
                    z=reserve*C.item(1,key,k5)+j*C.item(2,key,k5)+no2['values'][rem-reserve-2*j]
                    z+=(36*j*(j-1)if oldguard else 0)+(64*j*k5 if newguard else 0)
                    opts.append((z,reserve,j,True))
            options+=len(opts);best,reserve,j,e=max(opts)
            fs=fixed+[('A',1)]*reserve+[('A',2)]*j;left=rem-reserve-2*j;tab=no2 if e else normal
            while left:
                r=tab['last'][left];fs.append(C.representative(('A'if r<=8 else'D',r)));left-=r
            fixedcost=sum(C.upper_base(t)for t in fixed)+sum(C.upper_cross(a,b)for i,a in enumerate(fixed)for b in fixed[i+1:])
            bound=fixedcost+48*c16+best
            assert sum(t[1]for t in fs)==58 and ('A',1)in fs
            assert C.value(fs,True)+48*c16==bound
            rows.append(dict(counts=cs,key=key,fixed_A5_count=k5,residual_rank=rem,rank16_count=c16,
                upper=F(bound,64),winner=dict(reserve=reserve,explicit_A2=j,exclude2=e),attaining_relaxed_forest=fs,
                complete_context_options=[dict(value=v,reserve=a1,explicit_A2=a2,exclude2=x)for v,a1,a2,x in opts]))
    assert len(rows)==539
    maximum=max(r['upper']for r in rows);assert maximum==F(757,2)
    above=[r for r in rows if r['upper']>377];assert len(above)==1
    assert above[0]['counts']==(0,0,0,0,0,0,0,2,0)
    d17old=F(docs[N.TAG]['increasing_family_upper']['maximum']);assert d17old==F(2837,8)
    d17new=d17old+F(3,4)*(58//16);assert d17new==F(2855,8)<377
    return dict(complete_rank16_A1_classes=539,rows=rows,complete_context_options=options,
        maximum=maximum,above_curvature_relaxed_rows=1,
        inherited_D17_guarded_family_classes=759,D17_family_upper_after_trade=d17new,
        residual_item_weights_unchanged=True,old_certified_tables_reused=True)

def finite(docs):
    old=docs[N.TAG]['finite'];assert len(old['classes'])==39
    classes=[];assignments=0
    for cl in old['classes']:
        rows=[]
        for r in cl['rows']:
            fs=list(map(tuple,r['fixed_actual']));co=Counter(fs)
            drop=80*co['A',16]+16*co['D',16]
            # Outside the A16+A1 family the two A16 roles are unconditionally
            # zero. D14 needs A3 as well as A1/A2, so subtract its drop only
            # when there are no residual roots at all.
            if cl['residual_rank']==0:
                assert not any(t in K.SMALL for t in fs);drop+=64*co['D',14]
            upper=F(r['sixtyfourths_actual_maximum']-drop,64)
            rows.append(dict(fixed_actual=fs,inherited_actual_maximum=F(r['sixtyfourths_actual_maximum'],64),
                guaranteed_drop_outside_A16_A1=F(drop,64),upper=upper));assignments+=1
        classes.append(dict(counts=cl['counts'],fixed_A5_count=cl['fixed_A5_count'],
            fixed_A12_count=cl['fixed_A12_count'],residual_rank=cl['residual_rank'],
            upper=max(r['upper']for r in rows),actual_fixed_fibers=rows))
    assert assignments==748
    assert max(r['upper']for r in classes)==F(757,2)
    witness=[('D',17),('D',17),('D',24)]
    assert value(witness)==64*F(757,2)and context(witness)==0
    return dict(complete_predecessor_upper_classes=39,complete_actual_fixed_assignments=748,
        complement_upper_classes=classes,complement_upper=F(757,2),
        exact_actual_maximum=F(757,2),actual_attaining_forest=witness,
        attained_costs=[F(actual_base(t),64)for t in witness],
        exact_obstruction_class_census_performed=False,new_residual_DP_performed=False,
        complement_argument='Outside the D17 increasing family and A16+A1 family, pointwise decrease from the complete frozen N count; all39 predecessor upper classes and all748 actual fixed fibers are bounded explicitly.')

def run():
    for stem,hs in PINS.items():
        for (folder,ext),h in zip((('scripts','py'),('results','json'),('results','md')),hs):assert digest(f'{folder}/{stem}.{ext}')==h
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    assert N.C is K.C and N.G is K.G and C.L is L
    local=selected(docs);inter=interfaces(docs);print('LOCAL',len(local['whole_proofs']),inter['later_corners'],flush=True)
    for t in C.TYPES:
        assert actual_base(t)>=0
        assert actual_base(t)-N.actual_base(t)==-{('A',16):80,('D',16):16,('D',14):64}.get(t,0)
        assert actual_base(t)+(128 if t==('A',16)else 0)<=C.upper_base(C.representative(t))+(48 if C.representative(t)==('D',16)else 0)
        for u in C.TYPES:
            assert 0<=actual_cross(t,u)==actual_cross(u,t)<=N.actual_cross(t,u)<=C.upper_cross(C.representative(t),C.representative(u))
    assert F(actual_cross(('A',4),('A',12)),64)==8
    upper=conditional_upper(docs);count=finite(docs)
    assert docs[C.TAG]['ideal']['complete_profiles']==2552 and docs[C.TAG]['ideal']['absolute_surplus']=='36'
    excess=count['exact_actual_maximum']-377;assert excess==F(3,2)
    vertex_ratio=F(4205)/excess;minimum=-(-vertex_ratio.numerator//vertex_ratio.denominator);assert minimum==2804
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_RANK16_A1_TRADE_COMPOSED_AMBIENT58_FRONTIER',sources={p:digest(p)for p in sorted(paths)},
        baseline=N.TAG,additional_frozen_D14_rule=K.TAG,scale=64,local_proof=local,current_interfaces=inter,
        type_count=len(C.TYPES),simultaneous_type_pairs=len(C.TYPES)**2,increasing_family_upper=upper,finite=count,
        ideal_weights_changed=False,inherited_positive_ideal_profiles=2552,actual_cusp_surplus='36',
        finite_curvature='377',strict_global_constant='4205',finite_excess='3/2',
        necessary_ambient58_inequality='(3/2)*v_f >= 4205 + 36*v_i',minimum_finite_vertices=minimum,
        completed_dimension_bound=58,new_dimension_bound=None,geometric_example_claimed=False,
        scope='Full current-output geometric trade and exact complete global maximum using two complete conditional families and all39 inherited upper classes/748 actual fixed fibers. No exact new obstruction census, residual DP, convolution or higher-dimensional theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['exact_actual_maximum'],out['minimum_finite_vertices'],flush=True)
