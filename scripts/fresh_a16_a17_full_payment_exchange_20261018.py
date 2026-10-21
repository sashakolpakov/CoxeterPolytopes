#!/usr/bin/env python3
"""A16 cost110, A17 cost117: finish a sealed two-corner half exchange."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_d16_second_pair_zero_transfer_20261018 as N

C=N.C;J=N.J;V=J.V;T=J.T;A=J.A;Q=J.Q;D=J.D;P=J.P;ROOT=N.ROOT
TAG='fresh_a16_a17_full_payment_exchange_20261018'
A16='fresh_finite_a16_endpoint_sparse_half_transfer_20261018'
RANK15='fresh_rank15_sparse_payment_exchange_20261018'
ZERO={(0,11),(4,15)};RESTORE={(1,12),(4,15)}
STEMS=(C.TAG,V.TAG,T.TAG,T.T.TAG,A.OVER,A16,RANK15,C.Z.DEPTH)


def literal(current):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in C.L.STEMS}
    before={(r['source'],C.L.N.normalize(r['case'])):r for r in current['combined_literal_interface']['rows']}
    rows=[];changed=[]
    for stem,index,g,white,cs in T.literal_inputs(docs):
        old=before[stem,C.L.N.normalize(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==old['corners'][j]['marks'];ids=sorted(set(white)|set(marks));delta=F();images=[];typ=None
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                bs=D.finite_blocks(g,ids);shared=next((b for b in bs if set(marks)<=set(b['nodes'])),None)
                if shared and shared['type'] in ('A16','A17'):
                    typ=shared['type'];images=T.mark_pairs(g,shared['nodes'],marks,int(typ[1:]),'A')
                    change=ZERO if typ=='A16' else RESTORE
                    if set(images)&change:
                        assert set(images)<=change;delta=F(-1,2) if typ=='A16' else F(1,2)
            value=F(old['corners'][j]['retained_payment_lower_bound'])+delta;assert value>=0
            corners.append(dict(marks=marks,type=typ,canonical_images=images,delta=delta,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['delta'] for c in corners):changed.append(row)
    assert len(rows)==216
    return dict(configurations=216,changed=changed,rows=rows)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    current=docs[C.TAG];assert current['finite']['D17_cost']=='235/2'
    a16={tuple(r['pair']):F(r['weight']) for r in docs[A16]['weighted_pairs']}
    assert sum(a16.values())==111 and all(a16[p]==F(1,2) for p in ZERO)
    baseline=next(r for r in docs[T.T.TAG]['records'] if r['kind']=='A')
    a17={p:F(int(list(p) in baseline['retained_pairs'])) for p in combinations(range(17),2)}
    assert sum(a17.values())==114
    for p in V.ARESTORE:assert a17[p]==0;a17[p]=F(1,2)
    r15=tuple(docs[RANK15]['A17_restored_unit_pair']);assert r15==(1,15) and a17[r15]==0;a17[r15]=F(1)
    assert sum(a17.values())==116 and all(a17[p]==F(1,2) for p in RESTORE)
    for p in ZERO:a16[p]=F()
    for p in RESTORE:a17[p]=F(1)
    assert sum(a16.values())==110 and sum(a17.values())==117
    assert a16[1,13]==a16[2,14]==F(1,2) and a17[1,15]==1
    assert ZERO==V.AHALF and RESTORE==V.ARESTORE
    sup=V.support('A',docs[T.TAG]);target=docs[A16]['target'];over=docs[A.OVER]
    assert V.validate_maps(sup,target,over)==238 and sup['complete_eligible_patterns']==394
    _,dist,_=P.distance_data(target);assert all(dist[i][j]==11 for i,j in ZERO)
    assert all({tuple(sorted(m[i] for i in p)) for p in ZERO}==ZERO for m in P.embeddings(target,target))
    oldproof={r['index']:r for r in docs[V.TAG]['selected_proofs'] if r['kind']=='A'}
    hits=[];proofs=[]
    for row in sup['rows']:
        maps=[m for m in row['maps'] if tuple(m['marked_images']) in ZERO]
        if not maps:continue
        c=over['isolation']['cases'][row['index']];core=over['isolation']['cores'][c['core_key']]
        p=V.core_proof(c,core,'A',ZERO)
        assert V.enc(p)==oldproof[c['index']]
        assert p['sealing']['method']=='radical' and len(c['marked_component'])==16
        corners=[]
        for r in p['corners']:
            value=F()
            if r['reason']=='retained_selected_half':assert set(r['canonical_images'])<=ZERO
            elif r['reason']=='restored_whole_A17_half':assert set(r['canonical_images'])<=RESTORE;value=F(1)
            else:assert r['payment']==0
            corners.append(dict(marks=r['marks'],components=r['components'],canonical_images=r['canonical_images'],payment=value,reason=r['reason']))
        assert sum(r['payment'] for r in corners)==1
        proofs.append(dict(index=c['index'],old_seal=p,new_corners=corners,total=1))
        hits.append(dict(index=c['index'],maps=maps))
    assert [r['index'] for r in hits]==[1657,2067] and sum(len(r['maps']) for r in hits)==4
    assert docs[C.Z.DEPTH]['feature_depth']==10
    paths={Path(__file__),Path(N.__file__),Path(C.__file__),Path(J.__file__),Path(V.__file__),Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return J.enc(dict(status='PASS_A16_A17_FULL_PAYMENT_EXCHANGE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        old_costs=dict(A16='111',A17='116'),new_costs=dict(A16='110',A17='117'),
        A16_zero_pairs=sorted(ZERO),A17_restored_unit_pairs=sorted(RESTORE),
        A16_weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(a16.items())],
        A17_weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(a17.items())],
        complete_selected_patterns=394,validated_saved_maps=238,affected_contexts=hits,selected_proofs=proofs,
        literal_interface=literal(current),retained_rank15_A16_half_pairs=[[1,13],[2,14]],
        retained_rank15_A17_unit_pair=[1,15],all_ideal_weights_unchanged=True,all_cross_component_coupons_unchanged=True,
        finite_D17_cost_unchanged='235/2',triangle_minimum=2,quadrilateral_minimum=1,
        compatible_finite7_feature_depth=10,finite_forest_bound_proved=False,new_dimension_bound=None,
        scope='Complete saved394/238 A16 map validation, two freshly reconstructed arbitrary-support radical seals and216 current literal interfaces. Finish selected-half/alternate-half as selected-zero/alternate-unit. A17 increase1 must be counted; no forest solve, new ideal cost or dimension claim.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_costs'],'literal changes',len(out['literal_interface']['changed']),flush=True)
