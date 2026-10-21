#!/usr/bin/env python3
"""Bounded extra-blocker/deletion diagnostic after the five-66 endpoint scout."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_alternative66_endcaps_20261017 as E

H,B,M,ROOT=E.H,E.B,E.M,E.ROOT
TAG='fresh_nonarithmetic_alternative66_release59_20261017'

def topology(old,new):
    eo={tuple(x['sides']):x['kind']for x in old['noncompact_edges']}
    en={tuple(x['sides']):x['kind']for x in new['noncompact_edges']}
    return dict(old_noncompact_removed=[list(k)for k in eo if k not in en],
                new_noncompact_added=[dict(sides=list(k),kind=v)for k,v in en.items()if k not in eo])

def run():
    path=ROOT/'results'/f'{E.TAG}.json';saved=json.loads(path.read_text())
    z,labels,ns,h,profiles,primary=H.setup();old=H.A.setup()[0]
    gram=np.array([[B.mid(x)for x in row]for row in old['residual_gram']]);vals,vecs=np.linalg.eigh((gram+gram.T)/2);tr=vecs@np.diag(1/np.sqrt(abs(vals)))
    proofs=json.loads((ROOT/'results'/f'{E.G.TAG}.json').read_text())['records'][1]['exact_recipes']
    squares={r['wall']:B.mid(r['square'])for r in z['projections']};actual=z['actual_sides'];after=[i for i in actual if i!=59]
    cost=H.deletion(primary,z,(59,));stages=[];improvements=[]
    for s in saved['stages']:
        choice=s['choice'];proof=next(r for r in proofs if r['profile']==0 and r['labels']==choice)
        branch=next(r for r in proof['branches']if r['sign']==-1)
        ns[66]=np.linalg.solve(tr,np.array([B.mid(v)for v in branch['coordinates']]))
        control=s['control'];deleted=B.ledger(np.array([ns[i]for i in after]),after)
        left=[]
        for a in s['accepted_low_caps']:
            assert a['group']==0
            left.append(dict(index=a['index'],profile=a['profile'],labels=a['assigned_labels'],
                incidence=a['acute_incidence'],boundary=a['projective_boundary'],topology=topology(control,a['projective_boundary'])))
        right=[]
        for a in s['rejected_low_caps']:
            if a['group']!=1 or a['center_product']>=-1e-8 or {r['wall']for r in a['blockers']}!={59}:continue
            y=np.array(a['projection']);pi=a['profile']
            inc=H.incidence_score([ns[i]for i in after]+[y],after+[68],[squares[i]for i in after]+[profiles[pi]['square']])
            boundary=B.ledger(np.array([ns[i]for i in after]+[y]),after+[68])
            a=dict(a,index=len(right),incidence=inc,boundary=boundary,topology=topology(control,boundary));right.append(a)
            if inc.get('open_ends',99)<10:improvements.append(dict(choice=choice,kind='single_right_after_delete59',record=a))
        pairs=[];counts=Counter()
        for a,b in product(s['accepted_low_caps'],right):
            x,y=np.array(a['projection']),np.array(b['projection'])
            native=H.A.native(x,a['profile'],y,b['profile'],profiles);rel=B.relation(native);counts[str(rel)]+=1
            if not isinstance(rel,int):continue
            ty=E.shared_corner(a,b,profiles,rel)
            if ty is None:counts['finite_but_nonfinite_shared_corner']+=1;continue
            inc=H.incidence_score([ns[i]for i in after]+[x,y],after+[67,68],
                    [squares[i]for i in after]+[profiles[a['profile']]['square'],profiles[b['profile']]['square']])
            boundary=B.ledger(np.array([ns[i]for i in after]+[x,y]),after+[67,68])
            rec=dict(left_index=a['index'],right_index=b['index'],mutual_product=native,mutual_label=rel,
                     shared_corner_types=ty,incidence=inc,boundary=boundary,topology=topology(control,boundary))
            pairs.append(rec)
            if inc.get('open_ends',99)<10:improvements.append(dict(choice=choice,kind='paired_after_delete59',record=rec,roots=[a,b]))
        stages.append(dict(choice=choice,delete59_only=deleted,left_caps=left,right_caps_after_delete59=right,pair_counts=dict(counts),finite_pairs=pairs))
        print('CHOICE',choice,'DELETE59',deleted['open_ends'],'RIGHT',len(right),Counter(x['incidence'].get('open_ends',x['incidence']['status'])for x in right),'PAIRS',dict(counts),'VALID',len(pairs),flush=True)
    paths=[Path(__file__),Path(E.__file__),path,H.PARENT,H.BOUNDARY,ROOT/'results'/f'{E.G.TAG}.json']
    return dict(status='BOUNDED_NUMERICAL_EXTRA_BLOCKER59_AND_LEFT_TOPOLOGY_DIAGNOSTIC',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        exact_baseline_deletion59_cost=cost,stages=stages,improvements=improvements,
        scope='All left low caps and every future right low cap blocked solely by59 from the predecessor scout. Numerical full projective ledger plus independent PD-support incidence cross-check. Deletion cost is exact only for the frozen baseline66 control. Other four66 alternatives have only numerical boundary checks. No all-label tail exclusion, exact new geometry, or finite-volume claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/(TAG+'.json')
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'improvements',len(out['improvements']),flush=True)
