#!/usr/bin/env python3
"""Balance a sealed D16/D17 payment as5/8+3/8 instead of1/2+1/2."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_finite_d16_second_pair_zero_transfer_20261018 as N

C=N.C;J=N.J;T=N.T;A=N.A;Q=N.Q;D=N.D;P=N.P;ROOT=N.ROOT
TAG='fresh_d16_d17_eighth_balance_20261018'
UP={(0,14),(0,15)};DOWN={(1,15),(1,16)}
SUPPORT='fresh_finite_d17_five_unit_reduction_20261018'
STEMS=(N.TAG,C.TAG,J.TAG,SUPPORT,C.L.TAG,T.TAG,A.OVER,C.Z.DEPTH)


def literal(current):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in C.L.STEMS}
    old={(r['source'],C.L.N.normalize(r['case'])):r for r in current['literal_interface']['rows']}
    rows=[];changed=[];debits=credits=0
    for stem,index,g,white,cs in T.literal_inputs(docs):
        before=old[stem,C.L.N.normalize(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks));delta=F();images=[];typ=None
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                bs=D.finite_blocks(g,ids);b=next((b for b in bs if set(marks)<=set(b['nodes'])),None)
                if b and b['type'] in ('D16','D17'):
                    typ=b['type'];images=T.mark_pairs(g,b['nodes'],marks,len(b['nodes']),'D');change=UP if typ=='D16' else DOWN
                    if set(images)&change:
                        assert set(images)<=change;delta=F(1,8) if typ=='D16' else F(-1,8)
                        credits+=delta>0;debits+=delta<0
            value=F(before['corners'][j]['retained_payment_lower_bound'])+delta;assert value>=0
            corners.append(dict(marks=marks,type=typ,canonical_images=images,delta=delta,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['delta'] for c in corners):changed.append(row)
    assert len(rows)==216 and debits==0 and credits==4 and len(changed)==3
    return dict(configurations=216,debited_corners=debits,credited_corners=credits,changed=changed,rows=rows)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    weights16={tuple(r['pair']):F(r['weight']) for r in docs[N.TAG]['weighted_pairs']}
    weights17={tuple(r['pair']):F(r['weight']) for r in docs[C.L.TAG]['weighted_pairs']}
    assert sum(weights16.values())==110 and sum(weights17.values())==F(235,2)
    for p in UP:assert weights16[p]==F(1,2);weights16[p]=F(5,8)
    for p in DOWN:assert weights17[p]==F(1,2);weights17[p]=F(3,8)
    assert sum(weights16.values())==F(441,4) and sum(weights17.values())==F(469,4)
    assert weights16[1,11]==0 and weights16[1,14]==weights16[1,15]==F(1,2)
    target=docs[SUPPORT]['target'];sup=docs[SUPPORT]['complete_selected_support'];over=docs[A.OVER]
    assert sup['complete_eligible_patterns']==140 and sup['complete_maps']==168
    maps=0;hits=[]
    for row in sup['rows']:
        c=over['isolation']['cases'][row['index']];core=over['isolation']['cores'][c['core_key']]
        g=A.matrix(core['gram']);base=c['marked_component'];graph=T.T.labels(g,base);seen=set();pairs=set()
        for saved in row['maps']:
            mp=saved['images'];assert len(mp)==len(base)==len(set(mp))
            assert all(0<=i<17 for i in mp) and tuple(mp) not in seen;seen.add(tuple(mp));maps+=1
            assert all(graph[i][j]==target[mp[i]][mp[j]] for i in range(len(mp)) for j in range(len(mp)))
            pair=tuple(sorted(mp[base.index(m)] for m in c['preferred_corner']))
            assert list(pair)==saved['marked_images'];pairs.add(pair)
            if pair in DOWN:hits.append(dict(index=c['index'],map=mp))
        assert sorted(pairs)==list(map(tuple,row['marked_images']))
    assert maps==168 and not hits
    assert all({tuple(sorted(m[i] for i in p)) for p in DOWN}==DOWN for m in P.embeddings(target,target))
    old={r['index']:r for r in docs[J.TAG]['selected_proofs']};proofs=[]
    for i in (1644,1664,1931,2324):
        c=over['isolation']['cases'][i];core=over['isolation']['cores'][c['core_key']]
        p=J.seal(c,core);assert J.enc(p)==old[i]
        assert p['sealing']['method']=='radical' and set(p['selected_images'])==UP
        corners=[]
        for r in p['corners']:
            value=F()
            if r['reason']=='selected_D16_half':value=F(5,8);assert set(r['canonical_images'])==UP
            elif r['reason']=='restored_whole_D17_half':value=F(3,8);assert set(r['canonical_images'])==DOWN
            else:assert r['payment']==0
            corners.append(dict(marks=r['marks'],components=r['components'],canonical_images=r['canonical_images'],payment=value,reason=r['reason']))
        assert sum(r['payment'] for r in corners)==1
        proofs.append(dict(index=i,old_seal=p,new_corners=corners,total=1))
    assert docs[C.Z.DEPTH]['feature_depth']==10
    paths={Path(__file__),Path(N.__file__),Path(C.__file__),Path(J.__file__),Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return J.enc(dict(status='PASS_D16_D17_EIGHTH_BALANCE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        old_costs=dict(D16='110',D17='235/2'),new_costs=dict(D16='441/4',D17='469/4'),
        D16_upgraded_pairs=sorted(UP),D16_pair_weight='5/8',D17_reduced_pairs=sorted(DOWN),D17_pair_weight='3/8',
        D16_weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights16.items())],
        D17_weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights17.items())],
        complete_selected_D17_patterns=140,validated_saved_D17_maps=168,changed_selected_D17_maps=hits,
        affected_later_D16_contexts=[1644,1664,1931,2324],reweighted_proofs=proofs,
        literal_interface=literal(docs[N.TAG]),retained_rank15_D16_halves=[[1,14],[1,15]],
        all_ideal_weights_unchanged=True,all_cross_component_coupons_unchanged=True,
        triangle_minimum=2,quadrilateral_minimum=1,compatible_finite7_feature_depth=10,
        finite_forest_bound_proved=False,new_dimension_bound=None,
        scope='Complete saved140/168 D17 map validation gives no changed selected image. Reconstruct four later D16 radical proofs with distinct5/8+3/8 payments;216 current literals have no debit and four credits. All other finite/ideal interfaces remain inherited; no forest or ambient theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_costs'],'literal debits',out['literal_interface']['debited_corners'],flush=True)
