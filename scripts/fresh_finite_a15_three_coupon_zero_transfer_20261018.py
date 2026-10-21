#!/usr/bin/env python3
"""Two A15 halves become zero using three explicitly charged finite coupons.

Complete selected maps, fresh whole-core seals, exact current literal
interfaces and simultaneous typed costs. No D15 reduction is assumed.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_rank15_sparse_payment_exchange_20261018 as L
import fresh_a14_cd16_ambient59_frontier_20261018 as C
ROOT=L.ROOT;T=L.T;A=L.A;P=L.P;Q=L.Q;D=L.D
TAG='fresh_finite_a15_three_coupon_zero_transfer_20261018'
SCOUT=L.R.R.TAG;CHANGES={(0,12),(2,14)}
STEMS=tuple(dict.fromkeys((*C.STEMS,C.TAG,*L.STEMS,L.TAG,SCOUT)))
enc=C.enc
COUPONS=[
 dict(small='A1',small_marks=[0],large='A12',large_marks=[0,11],increment='1/4',component_pair_increment='1/2'),
 dict(small='A1',small_marks=[0],large='A13',large_marks=[1,11],increment='1/4',component_pair_increment='1/2'),
 dict(small='A4',small_marks=[1,2],large='A13',large_marks=[1,11],increment='1/2',component_pair_increment='2')]


def coupons(g,blocks,marks):
    aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
    if aa is bb:return []
    out=[]
    for sm,la,m,n in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
        for cp in COUPONS:
            if sm['type']!=cp['small']or la['type']!=cp['large']:continue
            si=T.mark_pairs(g,sm['nodes'],[m],int(sm['type'][1:]),'A')
            li=T.mark_pairs(g,la['nodes'],[n],int(la['type'][1:]),'A')
            if {p[0]for p in si}==set(cp['small_marks'])and {p[0]for p in li}==set(cp['large_marks']):
                out.append(dict(small=sm,large=la,marks=[m,n],small_images=si,large_images=li,increment=cp['increment']))
    assert len(out)<=1;return out


def actual_cross(a,b):
    pair={a,b};inc=8 if pair in ({('A',1),('A',12)},{('A',1),('A',13)})else 32 if pair=={('A',4),('A',13)}else 0
    return C.actual_cross(a,b)+inc


def upper_cross(a,b):
    pair={a,b};return C.upper_cross(a,b)+(1 if pair=={('A',1),('D',13)}else 16 if pair=={('A',4),('D',13)}else 0)


def actual_base(t):return C.actual_base(t)-(16 if t==('A',15)else 0)
def upper_base(t):return C.upper_base(t)-(8 if t==('D',15)else 0)


def support(rec,over):
    patterns=[c for c in over['isolation']['cases']if(c.get('component_order_bound')or 0)>=15 and len(c['marked_component'])<=15]
    assert len(patterns)==rec['complete_eligible_patterns']==536
    assert [c['index']for c in patterns]==[r['index']for r in rec['selected_patterns']]
    rows=[];hits=[];total=0
    for c,old in zip(patterns,rec['selected_patterns']):
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        maps=T.T.mapped(g,c['marked_component'],c['preferred_corner'],rec['target'])
        assert maps==old['maps'];total+=len(maps);rows.append(dict(index=c['index'],maps=maps))
        changed=[m for m in maps if tuple(m['marked_images'])in CHANGES]
        if changed:
            assert len(changed)==len(maps)==2;hits.append(dict(index=c['index'],maps=changed))
    assert total==482 and [r['index']for r in hits]==[1029,1612,1671,2416]
    return dict(complete_patterns=536,complete_maps=482,affected_maps=8,rows=rows,affected=hits,fresh_native_search=True)


def whole(c,core,a16,prior):
    old=L.a_proof(c,core,a16);assert enc(old)==prior
    assert c['index']in(1029,1612,1671,2416)and len(c['marked_component'])==15
    g=A.matrix(core['gram']);rows=[]
    for corner in old['corners']:
        marks=corner['marks'];cp=coupons(g,corner['components'],marks)
        debit=F(1,2)if marks==c['preferred_corner']else F()
        before=F(corner['payment']);inc=sum((F(r['increment'])for r in cp),F())
        pay=before-debit+inc;assert pay>=0
        rows.append(dict(marks=marks,whole_components=corner['components'],
            previous_payment=before,A15_debit=debit,new_coupons=cp,increment=inc,payment=pay,
            retained_A16_images=corner['canonical_images']if corner['reason']=='retained_whole_A16_half'else []))
    assert sum(r['payment']for r in rows)==1
    assert sum(r['A15_debit']for r in rows)==sum(r['increment']for r in rows)==F(1,2)
    assert sum(bool(r['new_coupons'])for r in rows)==(2 if c['index']in(1029,1612)else 1)
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],white=old['white'],
        inherited_seal_freshly_rechecked=old['sealing'],corners=rows,total='1')


def literal(docs):
    prior={C.C.R.key(r['source'],r['case']):r for r in docs[C.TAG]['combined_literal_interface']['rows']}
    rows=[];hits=Counter();seen=set()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        key=C.C.R.key(stem,index);assert key not in seen;seen.add(key)
        before=prior[key];corners=[];assert len(cs)==len(before['corners'])
        for j,(marks,*_)in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks));debit=F();cp=[];images=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                # The native literal catalogue is simply laced in every
                # corner where a new coupon can apply; recognize the two
                # marked components separately before finite_blocks.
                components=A.components(g,ids)
                shared=next((b for b in components if set(marks)<=set(b)),None)
                if shared and len(shared)==15:
                    images=T.mark_pairs(g,shared,marks,15,'A')
                    if set(images)&CHANGES:assert set(images)<=CHANGES;debit=F(1,2)
                aa,bb=[next(b for b in components if m in b)for m in marks]
                if aa is not bb and sorted([len(aa),len(bb)])in([1,12],[1,13],[4,13]):
                    # Use only the marked blocks; a disconnected exceptional
                    # spectator need not be classified by this helper.
                    tryblocks=[]
                    for b in(aa,bb):
                        images1=T.mark_pairs(g,b,[next(m for m in marks if m in b)],len(b),'A')
                        if images1:tryblocks.append(dict(nodes=b,type='A'+str(len(b))))
                    if len(tryblocks)==2:cp=coupons(g,tryblocks,marks)
            inc=sum((F(r['increment'])for r in cp),F())
            pay=F(before['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert pay>=0
            hits['A15_debits']+=bool(debit);hits['finite_credits']+=bool(inc)
            corners.append(dict(marks=marks,A15_images=images,debit=debit,increment=inc,new_coupons=cp,retained_payment_lower_bound=pay))
        total=sum(r['retained_payment_lower_bound']for r in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert seen==set(prior)and len(rows)==216
    return dict(configurations=216,changes=dict(hits),rows=rows,minimum_total=min(r['total_lower_bound']for r in rows),native_Grams_reconstructed=True)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    rec=next(r for r in docs[SCOUT]['records']if r['kind']=='A');over=docs[A.OVER]
    oldrule=next(r for r in docs[L.TAG]['rules']if r['kind']=='A')
    weights={tuple(r['pair']):F(r['weight'])for r in oldrule['weighted_pairs']}
    assert sum(weights.values())==101
    for p in CHANGES:assert weights[p]==F(1,2);weights[p]=F()
    assert sum(weights.values())==100 and weights[1,13]==1
    assert all(weights[tuple(p)]==1 for p in rec['retained_basic_pairs'])
    assert all({tuple(sorted(mp[i]for i in p))for p in CHANGES}==CHANGES for mp in rec['automorphisms'])
    _,dist,features=P.distance_data(rec['target']);assert all(dist[i][j]==12 for i,j in CHANGES)
    assert all(weights[p]==1 for p in weights if dist[p[0]][p[1]]<=9)
    a16={tuple(r['pair']):F(r['weight'])for r in docs[L.A16]['weighted_pairs']}
    assert a16[1,13]==a16[2,14]==F(1,2)
    sup=support(rec,over);prior={p['index']:p for p in docs[L.TAG]['selected_proofs']}
    proofs=[]
    for hit in sup['affected']:
        c=over['isolation']['cases'][hit['index']];core=over['isolation']['cores'][c['core_key']]
        p=whole(c,core,a16,prior[c['index']]);p['affected_maps']=hit['maps'];proofs.append(p)
    types=[t for t in C.C.P.Z.B.B.types()if t[1]<=59];rep=C.C.P.Z.B.B.representative
    for a in types:
        assert actual_base(a)<=upper_base(rep(a))
        for b in types:assert actual_cross(a,b)==actual_cross(b,a)<=upper_cross(rep(a),rep(b))
    assert actual_base(('A',15))==1600 and upper_base(('D',15))==1608
    costs=[]
    for cp in COUPONS:
        a=('A',int(cp['small'][1:]));b=('A',int(cp['large'][1:]))
        assert F(cp['component_pair_increment'])==len(cp['small_marks'])*len(cp['large_marks'])*F(cp['increment'])
        costs.append(dict(coupon=cp,old_actual=F(C.actual_cross(a,b),16),new_actual=F(actual_cross(a,b),16),
            old_representative=F(C.upper_cross(rep(a),rep(b)),16),new_representative=F(upper_cross(rep(a),rep(b)),16)))
    paths={Path(__file__),Path(L.__file__),Path(C.__file__),Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(Q.__file__),Path(P.__file__),Path(D.__file__),Path(L.R.X.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_FINITE_A15_THREE_COUPON_ZERO_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        target=rec['target'],automorphisms=rec['automorphisms'],changed_pairs=sorted(CHANGES),
        old_A15_cost='101',new_A15_cost='100',current_D15_cost='201/2',rank15_representative_upper='201/2',
        weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(weights.items())],retained_A15_unit=[1,13],retained_A16_halves=[[1,13],[2,14]],
        complete_selected_support=sup,selected_proofs=proofs,new_coupons=COUPONS,component_pair_accounting=costs,
        current_literal_interface=literal(docs),simultaneous_type_pairs=len(types)**2,
        inherited_high_overlap_routing=docs[L.TAG]['complete_high_overlap_routing'],
        all_ideal_weights_unchanged=True,other_finite_bases_unchanged=True,triangle_minimum=2,quadrilateral_minimum=1,
        D15_new_reduction_assumed=False,finite_forest_recomputed=False,new_dimension_bound=None,
        scope='Complete536-pattern/482-map A15 support, four freshly reconstructed arbitrary-white seals, retained whole A16 halves and three additive coupons. Current216 native literals and full type-pair domination. No uncounted coupon absorption or ambient theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_A15_cost'],out['current_literal_interface']['changes'],flush=True)
