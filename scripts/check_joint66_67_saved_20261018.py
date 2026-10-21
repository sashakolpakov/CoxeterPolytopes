#!/usr/bin/env python3
"""Independent saved-value/accounting check; no project imports or cone replay."""
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_joint66_67_native_prefixes_20261018'
M=np.diag([-1.,1.,1.,1.]);J=np.array([[2.,1.],[1.,2.]])/3


def native(x):
    if x < -2-2e-7:return 'dotted'
    if abs(x+2)<=2e-7:return 'parallel_or_unresolved'
    if abs(x)<=2e-7:return 2
    if x>0 or x < -2:return None
    m=int(round(np.pi/np.arccos(-x/2)))
    if m<3:return None
    theta=lambda n:2*np.cos(np.pi/n)
    gap=min(abs(theta(m)-theta(m-1)),abs(theta(m+1)-theta(m)))
    return m if abs(x+theta(m))<=min(2e-7,gap/100)else None


def ledger(l):
    assert not l['warnings'] and l['PD_status']=='numerically_separated_acute_incidence'
    assert l['open_ends']==l['PD_open_ends']==sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in l['edges'])
    assert Counter(r['kind']for r in l['rays'])==Counter(l['counts'])
    assert {tuple(r['active'])for r in l['rays']if r['kind']=='finite'}=={tuple(r)for r in l['PD_finite_supports']}


def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    assert d['deleted_walls']==[59,60,61,66] and d['retained_old_walls']==59 and not d['unresolved']
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vs=[tuple(v)for v in pri['retained_old_finite_supports']]+[tuple(v['support'])for v in pri['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    lost={i for i,v in enumerate(vs)if set(v)&{59,60,61,66}}
    assert len(lost)==9 and all(8 not in vs[i]for i in lost)
    edges=Counter(e for v in vs for e in combinations(v,8));faces=defaultdict(list)
    for i,v in enumerate(vs):
        for k in range(8):
            for f in combinations(v,k):faces[f].append(i)
    touched=[]
    for f,indices in faces.items():
        if not(lost&set(indices))or len(indices)<10-len(f):continue
        if all(edges[tuple(j for j in vs[i]if j!=k)]==2 for i in indices for k in vs[i]if k not in f):
            touched.append(dict(support=list(f),dimension=9-len(f),finite_supports=[list(vs[i])for i in indices]))
    debit=d['exact_deletion_debit']
    assert debit['lost_finite_count']==9 and debit['lost_finite_supports']==[list(vs[i])for i in sorted(lost)]
    assert touched==debit['compact_faces_touched'] and Counter(r['dimension']for r in touched)=={2:6,3:1}
    assert all(8 not in r['support']for r in touched)
    for l in d['boundary_ledgers']:ledger(l)
    assert len(d['boundary_ledgers'])==6
    assert d['boundary_ledgers'][d['control_score']]['counts']==dict(finite=26,exterior=6)
    assert d['boundary_ledgers'][d['control_score']]['open_ends']==10
    assert d['boundary_ledgers'][d['deletion_score']]['counts']==dict(finite=21,exterior=5)
    assert d['boundary_ledgers'][d['deletion_score']]['open_ends']==9
    contacts=[np.array([0.,0.])]
    for side in(0,1):
        for m in(3,4,5):
            c=np.zeros(2);c[side]=-2*np.cos(np.pi/m);contacts.append(c)
    squares=[2-float(c@J@c)for c in contacts]
    assert d['low_labels']==[2,3,4,5] and len(d['banks'])==len(d['bank_specs'])==2
    eqs=branches=caps=0;oldkeys=None
    for bi,b in enumerate(d['banks']):
        spec=d['bank_specs'][bi]
        assert all(b[k]==spec[k]for k in('new_wall','frame','duties'))
        assert b['frame']==([55,30,58]if bi==0 else[7,33,58]) and b['duties']==[[0,1],[0,2]]
        assert b['complete_low_recipe_count']==len(b['recipes'])==[96,13][bi]
        keys={(r['profile'],*r['labels'])for r in b['recipes']}
        assert len(keys)==len(b['recipes'])
        refs=[]
        for ri,r in enumerate(b['recipes']):
            eqs+=1
            assert r['index']==ri and 0<=r['profile']<7 and len(r['labels'])==3 and all(2<=v<=5 for v in r['labels'])
            assert r['norm']['status']=='two_real_norm_branches' and sorted(br['branch']for br in r['branches'])==[-1,1]
            for br in r['branches']:
                branches+=1
                eligible=not br['blockers']and br['center_product']<-1e-8 and br['norm_residual']<2e-5
                assert ('cap'in br)==eligible
                assert all(native(val)in(None,'parallel_or_unresolved')for wall,val in br['blockers'])
                if 'cap'in br:
                    cap=b['caps'][br['cap']];refs.append(br['cap'])
                    assert cap['recipe']==ri and cap['branch']==br['branch'] and cap['profile']==r['profile'] and cap['labels']==r['labels']
                    assert all(cap[k]==br[k]for k in('center_product','norm_residual'))
        assert sorted(refs)==list(range(len(b['caps'])))and len(b['caps'])==[70,5][bi]
        for ci,c in enumerate(b['caps']):
            caps+=1;assert c['index']==ci and c['center_product']<-1e-8 and c['norm_residual']<2e-5
            vals=dict(c['native_products']);assert len(vals)==len(c['native_products'])==53 and not(set(vals)&{59,60,61,66})
            if oldkeys is None:oldkeys=set(vals)
            assert set(vals)==oldkeys and all(native(v)not in(None,'parallel_or_unresolved')for v in vals.values())
            assert all(native(vals[w])==m for w,m in zip(b['frame'],c['labels']))
            v=np.array(c['normal']);assert abs(float(v@M@v)-squares[c['profile']])<2e-5
            l=d['boundary_ledgers'][c['score']];assert l['counts']==dict(finite=23,exterior=5)and l['open_ends']==9
            assert c['restores_old66_numerically']==(bi==1 and ci==2)
    assert(eqs,branches,caps)==(109,218,75)
    assert d['counts']==dict(norm_equations=109,norm_branches=218)
    assert [c['labels']for c in d['banks'][1]['caps']]==[[2,3,3],[2,3,4],[2,4,2],[2,4,3],[2,4,4]]
    digest=sha256();hist=Counter();accepted={r['pair']:r for r in d['accepted']};seen=set();restores=0
    assert len(accepted)==len(d['accepted'])==299
    exterior={(7,30,47),(30,33,68),(30,47,52),(30,67,68),(58,67,68)}
    cycle=[30,67,58,68,33];openedges={tuple(sorted(e))for e in zip(cycle,cycle[1:]+cycle[:1])}|{(7,47),(7,30),(47,52),(30,52)}
    assert len(d['pairs'])==d['complete_pair_count']==350
    for ix,x in enumerate(d['banks'][0]['caps']):
        for iy,y in enumerate(d['banks'][1]['caps']):
            pi=ix*5+iy;p=d['pairs'][pi]
            val=float(np.array(x['normal'])@M@np.array(y['normal']))+float(contacts[x['profile']]@J@contacts[y['profile']])
            rel=native(val);hist[str(rel)]+=1
            assert(p['index'],p['x'],p['y'],p['native_product'],p['relation'])==(pi,ix,iy,val,rel)
            digest.update(json.dumps([ix,iy,val,rel],separators=(',',':')).encode())
            if rel is None:assert p['score']is None and pi not in accepted;continue
            assert rel=='dotted' and pi in accepted
            a=accepted[pi];seen.add(pi);l=d['boundary_ledgers'][p['score']]
            assert a['score']==p['score'] and a['old66_restored']==y['restores_old66_numerically'];restores+=a['old66_restored']
            assert l['counts']==dict(finite=25,exterior=5)and a['finite_rays']==25 and l['open_ends']==a['open_ends']==9
            assert a['exterior_supports']==[r['active']for r in l['rays']if r['kind']=='exterior']
            assert {tuple(r)for r in a['exterior_supports']}==exterior
            assert a['open_edges']==[e for e in l['edges']if e['kind']in('one_open','two_open')]
            assert all(e['kind']=='one_open'for e in a['open_edges'])and{tuple(sorted(e['sides']))for e in a['open_edges']}==openedges
            assert len(a['whole_face_records'])==3
            for w,f in zip((55,67,68),a['whole_face_records']):
                assert f['wall']==w and f['rays']==[r for r in l['rays']if w in r['active']]
                assert f['all_rays_finite']==all(r['kind']=='finite'for r in f['rays'])==(w==55)
    assert seen==set(accepted)and restores==48
    assert hist==Counter({'dotted':299,'None':51})and dict(hist)==d['pair_histogram']and digest.hexdigest()==d['complete_pair_digest']
    return dict(status='PASS_JOINT66_67_SAVED',recipes=eqs,branches=branches,caps=caps,pairs=350,native_pairs=299,
        genuinely_changed66=251,old66_restorations=48,complete_pair_hashes=1,exact_debit=[9,6,1],
        scope='No project imports. Recomputed199-support/compact-face debit; saved low-domain recipe/branch/cap indices and native values; all350 mutual products and hash; all299 saved open-boundary/whole-face records and stored cone-PD equality. No physical normal/root, classifier, cone or interval reconstruction.')


if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
