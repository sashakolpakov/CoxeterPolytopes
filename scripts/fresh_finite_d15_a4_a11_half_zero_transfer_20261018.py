#!/usr/bin/env python3
"""Zero one finite D15 half using a distinct whole A4/A11 endpoint half.

Complete selected support, inherited arbitrary-contact proofs reconstructed,
current literal interfaces and exact unchanged representative domination.
No forest optimization or dimension theorem is performed here.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_finite_d15_three_singleton_half_transfer_20261018 as L
import fresh_a14_cd16_ambient59_frontier_20261018 as C

ROOT=L.ROOT;T=L.T;A=L.A;Q=L.Q;D=L.D;P=L.P;X=L.X;enc=L.enc
TAG='fresh_finite_d15_a4_a11_half_zero_transfer_20261018'
CHANGE={(0,11)};CASES=[1030,1613,1641,1682,2501]
FINAL_MDS={L.TAG:'7dcfbd0056b95d4b53a04c34eddce21d15dc0895e89b89d724e3e0d918ebe70b',
 C.TAG:'257219a7c1e324f44d4c7ed1c39d846a0994da807d2328f07b2387bad9ddb6cf'}
STEMS=tuple(dict.fromkeys((*C.STEMS,C.TAG,*L.STEMS,L.TAG)))


def coupon(g,blocks,marks):
    aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks];out=[]
    if aa is bb:return out
    for sm,la,m,n in((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
        if sm['type']!='A4'or la['type']!='A11':continue
        si=T.mark_pairs(g,sm['nodes'],[m],4,'A');li=T.mark_pairs(g,la['nodes'],[n],11,'A')
        if si==[(0,),(3,)]and li==[(0,),(10,)]:
            out.append(dict(small=sm,large=la,small_mark=m,large_mark=n,
                small_images=si,large_images=li,old_weight='0',new_weight='1/2'))
    assert len(out)<=1
    return out


def support(rec,over):
    eligible=[c for c in over['isolation']['cases']
        if(c.get('component_order_bound')or 0)>=15 and len(c['marked_component'])<=15]
    assert len(eligible)==rec['complete_eligible_patterns']==536
    assert[c['index']for c in eligible]==[r['index']for r in rec['selected_patterns']]
    total=0;rows=[];affected=[]
    for c,old in zip(eligible,rec['selected_patterns']):
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        maps=T.T.mapped(g,c['marked_component'],c['preferred_corner'],rec['target'])
        assert maps==old['maps'];total+=len(maps)
        hit=[m for m in maps if tuple(m['marked_images'])in CHANGE]
        rows.append(dict(index=c['index'],maps=maps))
        if hit:affected.append(dict(index=c['index'],maps=hit))
    assert total==870 and[r['index']for r in affected]==CASES
    assert sum(len(r['maps'])for r in affected)==10
    return dict(complete_patterns=536,complete_maps=870,fresh_native_search=True,
        rows=rows,affected=affected,affected_maps=10)


def proofs(docs,sup):
    over=docs[A.OVER];target=next(r for r in docs[L.R.R.TAG]['records']if r['kind']=='D')['target']
    oldproof={p['index']:p for p in docs[L.TAG]['selected_proofs']}
    d16={tuple(r['pair']):F(r['weight'])for r in docs[L.B.TAG]['D16_weighted_pairs']}
    out=[]
    for hit in sup['affected']:
        i=hit['index'];c=over['isolation']['cases'][i];core=over['isolation']['cores'][c['core_key']]
        # Reconstruct all inherited maps for this context, including the
        # unaffected (1,12) images in1641, before restricting to the new role.
        oldhit=next(r for r in docs[L.TAG]['complete_selected_support']['affected']if r['index']==i)
        old=L.proof(c,core,oldhit['maps'],target,d16)
        assert enc(old)==oldproof[i]
        g=A.matrix(core['gram']);n=len(g)
        assert core['gauge']==[1]*n and all(g[j][j]==2 for j in range(n))
        if i not in(1682,2501):
            pr=dict(old);pr['maps']=[m for m in old['maps']if tuple(m['marked_images'])in CHANGE]
            assert len(pr['maps'])==len(hit['maps'])
            out.append(dict(**pr,affected_maps=hit['maps'],inherited_proof_reconstructed=True));continue
        assert n==18 and len(c['marked_component'])==15
        assert sorted(F(x[1])for x in old['remaining_inverse_diagonals'])==[4,4,14]
        corners=[]
        for r in old['corners']:
            matches=coupon(g,r['components'],r['marks'])
            selected=r['marks']==c['preferred_corner']
            oldhalf=F(r['payment'])if r['reason']=='existing_whole_D16_half'else F()
            if oldhalf:assert r['canonical_images']==[(1,12)]and oldhalf==F(1,2)
            if selected:assert r['canonical_images']==[(0,11)]and not matches
            payment=oldhalf+sum((F(m['new_weight'])for m in matches),F())
            corners.append(dict(marks=r['marks'],components=r['components'],selected=selected,
                retained_D16_images=r['canonical_images']if oldhalf else[],
                retained_D16_payment=oldhalf,new_A4_A11_halves=matches,payment=payment))
        assert sum(r['payment']for r in corners)==1
        assert sum(bool(r['new_A4_A11_halves'])for r in corners)==1
        assert len({tuple(r['marks'])for r in corners if r['payment']})==2
        out.append(dict(index=i,core_key=c['core_key'],gram=core['gram'],white=old['white'],
            inherited_whole_component_seal=old,affected_maps=hit['maps'],corners=corners,total='1'))
    return out


def literal(docs):
    old={C.C.R.key(r['source'],r['case']):r for r in docs[C.TAG]['combined_literal_interface']['rows']}
    rows=[];hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=old[C.C.R.key(stem,index)];assert len(cs)==len(prev['corners']);corners=[]
        for j,(marks,*_)in enumerate(cs):
            assert marks==prev['corners'][j]['marks'];ids=sorted(set(white)|set(marks))
            images=[];matches=[];debit=F()
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                blocks=D.finite_blocks(g,ids);shared=next((b for b in blocks if b['type']=='D15'and set(marks)<=set(b['nodes'])),None)
                if shared:
                    images=T.mark_pairs(g,shared['nodes'],marks,15,'D')
                    if set(images)&CHANGE:assert set(images)<=CHANGE;debit=F(1,2)
                matches=coupon(g,blocks,marks)
            inc=sum((F(m['new_weight'])for m in matches),F())
            before=F(prev['corners'][j]['retained_payment_lower_bound']);after=before-debit+inc
            assert after>=0;hits['debits']+=bool(debit);hits['credits']+=bool(inc)
            corners.append(dict(marks=marks,previous_lower=before,D15_images=images,debit=debit,
                new_A4_A11_coupons=matches,increment=inc,retained_payment_lower_bound=after))
        total=sum(r['retained_payment_lower_bound']for r in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(rows)==216 and hits['debits']==0
    return dict(configurations=216,changed_roles=dict(hits),rows=rows,native_Grams_reconstructed=True,
        minimum_total=min(r['total_lower_bound']for r in rows))


def actual_base(t):return C.actual_base(t)-(8 if t==('D',15)else 0)
def upper_base(t):return C.upper_base(t)
def actual_cross(a,b):return C.actual_cross(a,b)+(32 if {a,b}=={('A',4),('A',11)}else 0)
def upper_cross(a,b):return C.upper_cross(a,b)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    for s,h in FINAL_MDS.items():assert sha256((ROOT/f'results/{s}.md').read_bytes()).hexdigest()==h,s
    assert docs[C.TAG]['finite']['exact_maximum']=='787/2'
    weights={tuple(r['pair']):F(r['weight'])for r in docs[L.TAG]['weighted_pairs']}
    assert sum(weights.values())==F(201,2)and weights[0,11]==F(1,2)
    weights[0,11]=F();assert sum(weights.values())==100
    rec=next(r for r in docs[L.R.R.TAG]['records']if r['kind']=='D');over=docs[A.OVER]
    sup=support(rec,over);prs=proofs(docs,sup)
    adj,dist,_=P.distance_data(rec['target']);assert dist[0][11]==11
    assert all(weights[tuple(p)]==1 for p in rec['retained_basic_pairs'])
    assert all({tuple(sorted(mp[i]for i in p))for p in CHANGE}==CHANGE for mp in rec['automorphisms'])
    # Explicit current positive obligations retain their distinct marked roles.
    d16now=docs['fresh_finite_d16_alternate_half_zero_transfer_20261018']
    d16weights={tuple(r['pair']):F(r['weight'])for r in d16now['D16_weighted_pairs']}
    assert d16weights[1,12]==F(1,2)
    d14=docs['fresh_finite_rank14_existing_payment_transfer_20261018']
    uses=[]
    for p in d14['selected_proofs']:
        for branch in[p,*p.get('extra_white_branches',[])]:
            for c in branch.get('corners',[]):
                if c.get('reason')=='existing_whole_D15_half':
                    assert c['images']==[[1,12]]and weights[1,12]==F(1,2)
                    uses.append(dict(index=p['index'],pair=[1,12],retained='1/2'))
    assert uses
    ideal15=[]
    for corner in docs['fresh_order15_bd_basic_fractional_transfer_20261018']['D_whole_finite_corners']:
        for role in corner['roles']:
            if role['type']=='D15':
                assert role['pair_images']==[[0,10]]and weights[0,10]==1
                ideal15.append(dict(case=corner['case'],pair=[0,10],retained='1'))
    assert ideal15
    types=[t for t in C.C.P.Z.B.B.types()if t[1]<=59];rep=C.C.P.Z.B.B.representative
    assert C.actual_cross(('A',4),('A',11))==0
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0
        for b in types:assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    assert actual_base(('D',15))==1600 and upper_base(('D',15))==1616
    assert actual_cross(('A',4),('A',11))==32 and upper_cross(('A',4),('D',11))==128
    cases=over['isolation']['cases']
    partition=[sum(0<(c.get('component_order_bound')or 0)<=14 for c in cases),
        sum(c.get('component_order_bound')is None for c in cases),
        sum((c.get('component_order_bound')or 0)>=15 and len(c['marked_component'])>15 for c in cases),536]
    assert len(cases)==2614 and partition==[1691,83,304,536]
    paths={Path(__file__),Path(L.__file__),Path(C.__file__),Path(T.__file__),Path(T.T.__file__),
        Path(A.__file__),Path(P.__file__),Path(Q.__file__),Path(D.__file__),Path(X.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_FINITE_D15_A4_A11_HALF_ZERO_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        target=rec['target'],changed_pairs=[[0,11]],old_weight='1/2',new_weight='0',old_cost='201/2',new_cost='100',
        weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(weights.items())],
        complete_selected_support=sup,selected_proofs=prs,whole_contexts=2,proper_contexts=3,
        complete_high_overlap_partition=partition,current_D14_D15_interfaces=uses,
        current_ideal15_D15_interfaces=ideal15,
        retained_D16_half=[1,12],new_coupon=dict(small='A4',small_marks=[0,3],large='A11',large_marks=[0,10],
            old_weight='0',new_weight='1/2',component_pair_increment='2',whole_components_required=True),
        typed_accounting=dict(simultaneous_type_pairs=len(types)**2,scale=16,new_actual_D15='100',
            unchanged_rank15_representative_upper='101',new_actual_A4_A11='2',unchanged_A4_D11_upper='8',
            all_upper_bases_and_cross_costs_unchanged=True,finite_forest_recomputed=False),
        current_literal_interface=literal(docs),triangle_minimum=2,quadrilateral_minimum=1,
        all_other_weights_unchanged=True,ideal_weights_unchanged=True,new_dimension_bound=None,
        scope='Complete536-pattern/870-map native reconstruction; five selected contexts with inherited arbitrary-support seals/proper exclusions reconstructed, distinct whole A4/A11 half plus retained whole D16 half, current216 literal and32,400 typed comparisons. Written complete routing/current-priority scope remains explicit; no new forest optimization or dimension claim.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    result=run()
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'D15',result['new_cost'],'maps',result['complete_selected_support']['affected_maps'],flush=True)
