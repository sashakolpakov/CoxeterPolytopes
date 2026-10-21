#!/usr/bin/env python3
"""Independent stdlib accounting of saved changed-bond and companion records."""
import json,math
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
A='fresh_wall59_bond_replacement_20261018';B='fresh_wall59_low_bond_companions_20261018'

def native(x):
    tol=2e-7
    if x<-2-tol:return 'dotted'
    if abs(x+2)<=tol:return 'parallel_or_unresolved'
    if abs(x)<=tol:return 2
    if x>0 or x<-2:return None
    m=round(math.pi/math.acos(-x/2))
    if m<3:return None
    th=lambda k:2*math.cos(math.pi/k)
    spacing=min(abs(th(m)-th(m-1)),abs(th(m+1)-th(m)))
    return m if abs(x+th(m))<=min(tol,spacing/100)else None

def ledger(l):
    assert not l['warnings'] and l['open_ends']==l['PD_open_ends']
    assert l['PD_status']=='numerically_separated_acute_incidence'
    assert Counter(r['kind']for r in l['rays'])==Counter(l['counts'])
    assert {tuple(sorted(r['active']))for r in l['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in l['PD_finite_supports']}
    assert sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in l['edges'])==l['open_ends']

def topology(t):
    fs=t['triangles'];ee=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));vs=set(sum(fs,[]));bc=t['boundary']
    assert (len(vs),len(ee),len(fs))==(t['vertices'],t['edges'],t['faces'])
    assert len(vs)-len(ee)+len(fs)==1 and set(ee.values())=={1,2}
    assert {e for e,n in ee.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    for v in vs:
        le=[tuple(x for x in f if x!=v)for f in fs if v in f];deg=Counter(x for e in le for x in e);seen={next(iter(deg))}
        while True:
            new=seen|{b for a,b in le if a in seen}|{a for a,b in le if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(deg) and sorted(deg.values())==([1,1]+[2]*(len(deg)-2)if v in bc else[2]*len(deg))

def run():
    pa=ROOT/f'results/{A}.json';pb=ROOT/f'results/{B}.json';a=json.loads(pa.read_text());b=json.loads(pb.read_text());count=Counter()
    for d in(a,b):
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h;count['source_bindings']+=1
        for l in d['boundary_ledgers']:ledger(l);count['saved_ledgers']+=1
        assert not d['unresolved']
    old=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vs=[tuple(v)for v in old['retained_old_finite_supports']]+[tuple(r['support'])for r in old['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199 and all(len(v)==9 for v in vs)
    edge=Counter(e for v in vs for e in combinations(v,8));faces=defaultdict(list)
    for v in vs:
        for k in range(8):
            for w in combinations(v,k):faces[w].append(v)
    for debit in a['exact_deletion_debits']+[b['exact_deletion_debit']]:
        release=set(debit['deleted_walls']);lost=[v for v in vs if release.intersection(v)]
        assert list(map(list,lost))==debit['lost_finite_supports'] and len(lost)==debit['lost_finite_count']==7
        assert all(8 not in v for v in lost)
        expected={}
        for w,vertices in faces.items():
            if not any(release.intersection(v)for v in vertices)or len(vertices)<10-len(w):continue
            if all(edge[tuple(x for x in v if x!=q)]==2 for v in vertices for q in v if q not in w):expected[w]=vertices
        actual={tuple(f['support']):list(map(tuple,f['finite_supports']))for f in debit['compact_faces_touched']}
        assert expected==actual and Counter(9-len(w)for w in expected)=={2:6,3:1};count['exact_deletion_scopes']+=1
    assert len(a['records'])==91 and a['frame']==[55,30,58] and a['profile']==0
    assert [r['labels']for r in a['records']]==[[2,m,n]for m,n in product(a['label30_samples'],a['label58_samples'])]
    forward=Counter();bybranch={}
    for i,r in enumerate(a['records']):
        assert r['index']==i and r['norm']['status']=='two_real_norm_branches'
        assert sorted(x['branch']for x in r['branches'])==[-1,1]
        for x in r['branches']:
            assert x['norm_error']<2e-5 and len(x['normal'])==4
            vals=dict(x['native_products']);assert len(vals)==56
            assert x['blockers']==[j for j in vals if native(vals[j])in(None,'parallel_or_unresolved')]
            v=x['normal'];assert abs(-v[0]**2+sum(z*z for z in v[1:])-2)<2e-5
            bybranch[i,x['branch']]=x
            if x['center_product']<0:forward[tuple(x['blockers'])]+=1
            count['replacement_branches']+=1
    assert forward=={(60,61):88,(60,61,66):2,():1}
    for i,r in enumerate(a['accepted']):
        assert r['index']==i and r['normal']==bybranch[r['record'],r['branch']]['normal']
        x=bybranch[r['record'],r['branch']];vals=dict(x['native_products']);release=set(r['release'])
        assert not(release-set([59,60,61])) and x['center_product']<0
        assert all(native(v)not in(None,'parallel_or_unresolved')for j,v in vals.items()if j not in release)
        assert r['finite_old_contacts']==[[j,native(v)]for j,v in sorted(vals.items())if j not in release and isinstance(native(v),int)]
        l=a['boundary_ledgers'][r['score']];rs=[x for x in l['rays']if 55 in x['active']]
        assert r['whole55_quadrilateral'] and all(x['kind']=='finite'for x in rs)
        assert {tuple(sorted(x['active']))for x in rs}=={(7,30,55),(7,55,58),(30,55,67),(55,58,67)}
        count['accepted_replacements']+=1
    assert Counter(tuple(r['release'])for r in a['accepted'])=={(59,):1,(59,60):1,(59,60,61):89}
    assert len(a['low58_label_accepted_indices'])==47
    for t in b['topologies']:topology(t);count['disk_topologies']+=1
    assert len(b['topologies'])==10 and b['seed_labels']==[[2,2,3],[2,4,2]]
    expectedpairs=set();bankcounts=[]
    for si,s in enumerate(b['seeds']):
        seed=a['accepted'][s['source_index']];assert seed['release']==[59,60,61] and seed['labels']==s['labels'] and seed['normal']==s['normal']
        l=b['boundary_ledgers'][s['score']];edges={tuple(sorted(e['sides'])):e['kind']for e in l['edges']};bc=b['cycle']
        ce={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
        assert all(edges[e]=='one_open'for e in ce)
        assert {e for e,k in edges.items()if k in['one_open','two_open']}-ce=={(7,47),(7,30),(30,52),(47,52)}
        bankcounts.append([len(v['caps'])for v in s['banks']])
        for rot,bank in enumerate(s['banks']):
            assert bank['frame']==[bc[(rot+j)%5]for j in range(3)] and bank['duties']==[[0,1],[1,2]]
            for ci,c in enumerate(bank['caps']):
                assert c['index']==ci and c['center_product']<0 and c['norm_residual']<2e-5
                assert len(c['native_products'])==55 and all(native(x)not in(None,'parallel_or_unresolved')for j,x in c['native_products'])
                assert b['boundary_ledgers'][c['score']]['open_ends']==9;count['companion_caps']+=1
            other=s['banks'][(rot+2)%5]
            expectedpairs|={(si,rot,i,j)for i in range(len(bank['caps']))for j in range(len(other['caps']))}
    assert bankcounts==[[0,62,0,36,0],[0,138,0,36,0]]
    actualpairs=set();ph=Counter()
    for r in b['compatible_pairs']:
        key=(r['seed'],r['rotation'],r['x'],r['y']);assert key not in actualpairs;actualpairs.add(key)
        assert native(r['native_product'])==r['relation']=='dotted'
        l=b['boundary_ledgers'][r['score']];assert l['open_ends']==9 and l['counts']==dict(finite=29,exterior=5)
        assert {(7,30,47),(30,47,52)}<={tuple(sorted(z['active']))for z in l['rays']if z['kind']=='exterior'}
        ph[r['seed'],r['rotation'],r['relation']]+=1;count['compatible_pairs']+=1
    assert actualpairs==expectedpairs and len(actualpairs)==7200
    assert ph==Counter({(r['seed'],r['rotation'],r['relation']):r['count']for r in b['pair_histograms']})
    assert not b['finite_junctions']and not b['dynamic_z_banks']and not b['completions']
    assert count['companion_caps']==272 and count['replacement_branches']==182
    return dict(status='PASS_SAVED_WALL59_REPLACEMENT_AND_COMPANIONS',counts=dict(count),
        source_hash=sha256(Path(__file__).read_bytes()).hexdigest(),certificate_hashes={A:sha256(pa.read_bytes()).hexdigest(),B:sha256(pb.read_bytes()).hexdigest()},
        scope='Independent stdlib saved-record verification: reconstruct old199 incidence deletion/compact-face debit, check native values and profile0 norm residuals, exact disk links, complete saved branch and pair index domains, boundary/open-end accounting and saved cone-versus-PD equality. No physical normal/root/cone or inherited spherical-classifier replay; floating saved values are not interval certificates or an ambient finite-volume proof.')

if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
