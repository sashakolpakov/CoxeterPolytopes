#!/usr/bin/env python3
"""Saved P-duty/event accounting only; no project imports or geometry replay."""
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_p12_right_duty_20261018'

def ledger(l,finite,ends):
    assert l['counts']=={'finite':finite,'exterior':6 if ends==10 else 7}
    assert l['open_ends']==l['PD_open_ends']==ends and not l['warnings']
    assert{tuple(r['active'])for r in l['rays']if r['kind']=='finite'}==set(map(tuple,l['PD_finite_supports']))

def run():
    path=ROOT/f'results/{TAG}.json';d=json.loads(path.read_text());c=d['closure'];counts=Counter()
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h;counts['bindings']+=1
    assert d['P_old_recipe']['labels']==[3,3,2]and d['P_new_recipe']['labels']==[3,2,2]
    assert d['P_new_recipe']['frame']==[7,12,67]and d['P_new_recipe']['profile']==0
    assert len(d['u_grid'])==131 and len(d['records'])==262 and len(d['accepted'])==131
    assert len(set((r['u'],r['P_sheet'])for r in d['records']))==262
    assert{(r['u'],r['P_sheet'])for r in d['records']}=={(u,s)for u in d['u_grid']for s in(-1,1)}
    assert d['deleted_walls']==[47]and d['retained_old_wall_count']==62
    assert d['exact_deletion_debit']['lost_finite_count']==2 and d['exact_deletion_debit']['lost_supports_containing8']==0
    for i in d['accepted']:
        r=d['records'][i];assert r['index']==i and r['P_sheet']==1 and not r['L_blockers']and not r['P_blockers']
        assert r['L_center']<0 and r['P_center']<0 and r['P_norm_error']<2e-5 and r['frame_error']<2e-5
        ledger(d['boundary_ledgers'][r['score']],27,11);counts['accepted_prefixes']+=1
    assert len(c['domain']['u_grid'])==122 and c['domain']['LQ_labels']==[3,4,12,24]
    assert len(c['branches'])==144 and len({tuple(r['key'])for r in c['branches']})==144
    assert sum(r['real_samples']for r in c['branches'])==17568
    assert sum(len(r['brackets'])for r in c['branches'])==len(c['brackets'])==33
    assert len(c['rejected_events'])+len(c['events'])==33 and len(c['events'])==8
    prior=json.loads((ROOT/'results/fresh_qj_connector_20261018.json').read_text())
    kr=[]
    for e in c['events']:
        assert e['index']==counts['events'];counts['events']+=1
        assert abs(e['QR'])<2e-7
        for r in e['native']:
            assert not r['blockers']and r['center_product']<0 and r['norm_error']<2e-5
            counts['event_normals']+=1
        ledger(c['boundary_ledgers'][e['score']],32,10)
        assert len(e['K_checks'])==len(e['J_bank']['caps'])
        for j,k in zip(e['J_bank']['caps'],e['K_checks']):
            assert k['J']==j['index']and j['profile']==0 and j['labels'][:2]==[2,2]
            assert j['labels'][2]in(3,4)and j['center_product']<0
            assert dict(j['native_products'])[69]<-2
            ledger(c['boundary_ledgers'][j['score']],34,10);counts['J']+=1
            domain=prior['joint_finite_domains'][','.join(map(str,k['domain_key']))]
            assert len(k['records'])==len(domain)and k['condition']<1e10 and not k['accepted']
            for i,(r,rec)in enumerate(zip(k['records'],domain)):
                assert r['domain_record']==i and rec['high_position']is None
                assert r['profile']==rec['profile']and r['labels']==rec['labels']
                assert abs(F(str(r['norm_residual'])))>F(3,500)
                kr.append((e['index'],j['index'],r));counts['K_recipes']+=1
    assert counts['J']==5 and counts['K_recipes']==255
    assert min(kr,key=lambda r:abs(r[2]['norm_residual']))[:2]==(1,0)
    assert not d['unresolved']and not c['unresolved']and not c['improvements']
    return dict(status='PASS_SAVED_P12_RIGHT_DUTY',counts=dict(counts),source_hash=sha256(Path(__file__).read_bytes()).hexdigest(),
        certificate_hash=sha256(path.read_bytes()).hexdigest(),scope='Checks saved hashes, changed duty, sample/event/recipe coverage, native acceptance records, residual gaps, deletion debit and saved cone/PD incidence agreement. No producer/numerical imports, native-normal reconstruction, root/cone replay or all-parameter certification.')
if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
