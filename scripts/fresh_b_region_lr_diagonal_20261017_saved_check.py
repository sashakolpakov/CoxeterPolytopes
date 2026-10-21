#!/usr/bin/env python3
"""Saved incidence/count checks only; no producer or numerical imports."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_b_region_lr_diagonal_20261017'


def main():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    old=json.loads((ROOT/'results/fresh_b_region_dotted7_four_cap_20261017.json').read_text())
    assert d['exact_deletion_debit']==old['chosen_debit']
    p=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vertices=p['retained_old_finite_supports']+[r['support']for r in p['new_finite_vertices']]
    assert len(vertices)==199
    lost=[v for v in vertices if 47 in v]
    assert lost==d['exact_deletion_debit']['lost_finite_supports'] and len(lost)==2
    assert all(8 not in v for v in lost)
    fs=d['disk']['triangles'];es=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));bc=d['disk']['boundary']
    assert (len(set(sum(fs,[]))),len(es),len(fs))==(8,17,10)
    assert (68,69)not in es and es[67,70]==2
    assert {e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    assert all(n in(1,2)for n in es.values())
    for rec in d['disk']['links']:
        v=rec['vertex'];link=[tuple(i for i in f if i!=v)for f in fs if v in f]
        assert sorted(link)==sorted(map(tuple,rec['edges']))
        degree=Counter(i for e in link for i in e);seen={next(iter(degree))}
        while True:
            new=seen|{i for e in link if set(e)&seen for i in e}
            if new==seen:break
            seen=new
        assert seen==set(degree)
        assert sorted(degree.values())==([1,1]+[2]*(len(degree)-2)if v in bc else[2]*len(degree))
    assert [sum(i in f for f in fs)for i in[67,68,69,70]]==[5,4,4,5]
    for ledger in d['boundary_ledgers']:
        assert not ledger['warnings']
        assert dict(Counter(r['kind']for r in ledger['rays']))==ledger['counts']
        finite={tuple(sorted(r['active']))for r in ledger['rays']if r['kind']=='finite'}
        assert finite==set(map(tuple,ledger['PD_finite_supports']))
        assert ledger['open_ends']==ledger['PD_open_ends']
        assert ledger['open_ends']==sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in ledger['edges'])
    assert len(d['brackets'])==len(d['events'])==22 and len(d['samples'])==124
    assert d['counts']==dict(L_samples=420,native_L_samples=124,P_evaluations=9292,native_P_evaluations=2974,
                            R_screen_evaluations=19648,native_R_without_L_evaluations=6266,LR_event_roots=22)
    assert d['bank_counts']==dict(norm_equations=1012,norm_branches=2024)
    assert sum(x['count']for x in d['bank_blockers'])==2024
    assert not d['completions'] and not d['improvements'] and not d['unresolved'] and not d['bank_unresolved']
    assert Counter(e['LR_label']for e in d['events'])=={4:1,5:9,6:3,8:3,12:3,24:3}
    exterior={(30,33,66),(30,52,70),(30,59,66),(52,67,70),(58,60,66),(59,60,66)}
    for e in d['events']:
        assert e['L_profile']==0 and e['R'][0]==0 and e['R'][1][:2]==[2,2]
        l=d['boundary_ledgers'][e['score']]
        assert l['counts']=={'finite':30,'exterior':6}and l['open_ends']==10
        assert {tuple(r['active'])for r in l['rays']if r['kind']=='exterior'}==exterior
    near=d['near_native_rejections']
    assert Counter(tuple(i for i,v in r['blockers'])for r in near)=={(70,):440,(68,70):27,(7,70):106}
    assert all(0<=r['event']<22 and r['center_product']<0 and dict(r['blockers'])[70]>0 for r in near)
    only=[r for r in near if len(r['blockers'])==1]
    print('PASS 14 bindings, same two-vertex deletion debit, flipped disk, five ledgers, 22 events, 2024 rejected Q branches')
    print('440 R-only saved QR products:',min(r['blockers'][0][1]for r in only),max(r['blockers'][0][1]for r in only))


if __name__=='__main__':main()
