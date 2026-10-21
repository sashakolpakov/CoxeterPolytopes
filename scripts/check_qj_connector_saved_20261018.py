#!/usr/bin/env python3
"""Independent saved-record checks; no numerical normal or cone reconstruction."""
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_qj_connector_20261018'

def run():
    p=ROOT/f'results/{TAG}.json';d=json.loads(p.read_text());checks=Counter()
    for name,h in d['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h;checks['bindings']+=1
    t=d['disk'];fs=t['triangles'];vs=set(sum(fs,[]));es=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2))
    assert len(fs)==14 and len(vs)==10 and len(es)==23 and len(vs)-len(es)+len(fs)==1
    assert (69,71)not in es and set(es.values())=={1,2}
    b=t['boundary'];assert{e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(b,b[1:]+b[:1])}
    for v in vs:
        le=[tuple(i for i in f if i!=v)for f in fs if v in f];degree=Counter(x for e in le for x in e);seen={next(iter(degree))}
        while True:
            nxt=seen|{i for e in le if set(e)&seen for i in e}
            if nxt==seen:break
            seen=nxt
        assert seen==set(degree)
        assert sorted(degree.values())==([1,1]+[2]*(len(degree)-2)if v in b else[2]*len(degree))
        checks['vertex_links']+=1
    assert d['deleted_walls']==[47]and d['retained_old_wall_count']==62
    debit=d['exact_deletion_debit'];assert debit['lost_finite_count']==len(debit['lost_finite_supports'])==2
    assert all(47 in s and 8 not in s for s in debit['lost_finite_supports'])and debit['lost_supports_containing8']==0
    assert len(debit['compact_faces_touched'])==1 and debit['compact_faces_by_dimension']=={'2':1}
    assert len(d['seeds'])==15 and len(d['records'])==855
    counts=Counter();seen=set()
    domains=d['joint_finite_domains'];assert{key:len(v)for key,v in domains.items()}=={'3,4':43,'3,3':63,'2,4':88}
    for seed in d['seeds']:
        assert seed['index']==checks['seeds'];checks['seeds']+=1
        assert seed['condition']<1e10
        ledger=d['boundary_ledgers'][seed['score']]
        assert ledger['counts']=={'finite':34,'exterior':6}and ledger['open_ends']==10 and not ledger['warnings']
        assert ledger['PD_open_ends']==10
        assert{tuple(r['active'])for r in ledger['rays']if r['kind']=='finite'}==set(map(tuple,ledger['PD_finite_supports']))
    for r in d['records']:
        pair=(r['seed'],r['domain_record']);assert pair not in seen;seen.add(pair)
        recipe=domains[d['seeds'][r['seed']]['domain_key']][r['domain_record']]
        assert all(r[k]==recipe[k]for k in('profile','labels','high_position'))
        counts['recipes']+=1
        if r['high_position']is None:
            counts['low_recipes']+=1
            assert abs(F(str(r['norm_residual'])))>F(1,125)
        else:
            counts['high_families']+=1;assert r['high_position']in(0,1)and r['profile']==0
            a,b,c=map(lambda x:F(str(x)),r['quadratic'])
            # Exact arithmetic on SAVED floating coefficients only. This does
            # not enclose the physical coefficients or establish exact geometry.
            lo=F(433,250);assert lo*lo<3 and a>0 and b>=0
            assert a*lo*lo+b*lo+c>0
            assert not r['roots_in_high_domain']and all(x<lo or x>=2 for x in r['real_roots'])
            checks['saved_quadratic_tail_signs']+=1
    assert counts==d['counts']=={'recipes':855,'low_recipes':847,'high_families':8}
    expected={(s['index'],i)for s in d['seeds']for i in range(len(domains[s['domain_key']]))};assert seen==expected
    assert not d['candidates']and not d['improvements']and not d['unresolved']
    closest=min((r for r in d['records']if r['high_position']is None),key=lambda r:abs(r['norm_residual']))
    assert closest==d['diagnosis']['closest_low_record']and closest['seed']==14 and closest['labels']==[4,2,2,2]
    checks.update(counts)
    return dict(status='PASS_SAVED_QJ_CONNECTOR_RECORDS',checks=dict(checks),source_hash=sha256(Path(__file__).read_bytes()).hexdigest(),
        certificate_hash=sha256(p.read_bytes()).hexdigest(),scope='Saved topology, source hashes, complete saved recipe index domains, residual separations, high-quadratic signs for rounded saved coefficients, and recorded cone/PD incidence agreement. No producer imports, physical reconstruction, numerical root/cone replay, inherited classifier proof or exact geometric exclusion.')
if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
