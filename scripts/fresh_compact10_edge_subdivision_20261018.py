#!/usr/bin/env python3
"""Replace the 10--30 hinge by a finite four-corner connector; bounded scout."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import numpy as np
import fresh_joint30_33_compact10_first_20261018 as K

H,T,ROOT=K.H,K.T,K.ROOT
TAG='fresh_compact10_edge_subdivision_20261018'


def topology():
    old=[[10,67,33],[10,67,44]]
    triangles=[[10,68,33],[67,68,33],[10,68,44],[67,68,44]]
    boundary=[10,33,67,44]
    edges=Counter(tuple(sorted(e)) for f in triangles for e in combinations(f,2))
    assert (len(set(sum(triangles,[]))),len(edges),len(triangles))==(5,8,4)
    assert set(edges.values())=={1,2} and (10,67) not in edges
    assert {e for e,n in edges.items() if n==1}=={tuple(sorted(e)) for e in zip(boundary,boundary[1:]+boundary[:1])}
    links=[]
    for v in sorted(set(sum(triangles,[]))):
        ee=[tuple(x for x in f if x!=v) for f in triangles if v in f]
        degree=Counter(x for e in ee for x in e);seen={next(iter(degree))}
        while True:
            new=seen|{b for a,b in ee if a in seen}|{a for a,b in ee if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(degree)
        assert sorted(degree.values())==([2]*4 if v==68 else [1,1]+[2]*(len(degree)-2))
        links.append(dict(vertex=v,edges=ee))
    return dict(old_triangles=old,triangles=triangles,boundary=boundary,removed_edge=[10,67],connector=68,links=links)


def frame_checks(s):
    p=s.primary;dec=lambda x:tuple(map(F,x))
    def mul(a,b):
        v=[x*y for x in a for y in b];return min(v),max(v)
    h=list(map(dec,p['common_center']));g=[[dec(q) for q in row] for row in p['computational_gram']]
    normals={}
    for r in p['facet_witnesses']:
        hp=dec(p['center_products'][str(r['wall'])]);assert hp[1]<0
        normals[r['wall']]=[mul((2*(a[0]-b[1]),2*(a[1]-b[0])),(1/hp[1],1/hp[0])) for a,b in zip(h,map(dec,r['point']))]
    pairs=sorted({tuple(sorted(v)) for frame in ([7,55,33],[10,44,33]) for v in combinations(frame,2)})
    out=[]
    for i,j in pairs:
        q=(F(),F())
        for a,x in enumerate(normals[i]):
            for b,y in enumerate(normals[j]):
                z=mul(mul(x,g[a][b]),y);q=(q[0]+z[0],q[1]+z[1])
        val=s.dot(s.old[i],s.ps(i,{}),s.old[j],s.ps(j,{}))
        assert float(q[0])-2e-5<=val<=float(q[1])+2e-5
        out.append(dict(walls=[i,j],profiles=[s.ps(i,{}),s.ps(j,{})],native_product=val,
            relation=T.G.relation(val),exact_native_enclosure=list(map(str,q))))
    return out


def run():
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s);checks=frame_checks(s)
    debit=T.H.deletion(s.primary,s.cur,[30])
    assert debit['lost_finite_count']==13 and debit['compact_faces_by_dimension']=={2:10,3:1}
    assert all(8 not in v for v in debit['lost_finite_supports'])
    base={i:v for i,v in s.old.items() if i!=30};ids=[i for i in s.cur['actual_sides'] if i!=30]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63')
    deleted=s.score(base,{},ids,'delete30');assert control is not None
    rbank=K.bank(s,base,ids,[7,55,33],[(0,1)],[],67,'R_bank',30)
    print('SUBDIVISION_R',len(rbank['recipes']),len(rbank['caps']),flush=True)
    jbank=K.bank(s,base,ids,[10,44,33],[(0,1),(0,2)],[],68,'J_bank',30)
    print('SUBDIVISION_J',len(jbank['recipes']),len(jbank['caps']),flush=True)
    selected=[]
    for r in rbank['caps']:
        vals=dict(r['native_products']);relations={i:T.G.relation(vals[i]) for i in (10,44)}
        r['changed_incidence_relations']=relations
        dot=float(np.array(r['normal'])@T.M@s.old[10])
        determinant=s.profiles[r['profile']]['square']*s.profiles[s.ps(10,{})]['square']-dot*dot
        r['projected10_pair']=dict(product=dot,determinant=determinant)
        if abs(determinant)<=1e-8:s.unresolved.append(dict(stage='R10_pair',cap=r['index'],determinant=determinant))
        r['eligible_for_subdivision']=dot<0 and determinant < -1e-8 and isinstance(relations[44],int)
        if r['eligible_for_subdivision']:selected.append(r['index'])
    pairs=[];hist=Counter();digest=sha256();coherent=[];completed=[]
    duties=topology()['triangles']+[[7,55,67]]
    for ri,ji in product(selected,range(len(jbank['caps']))):
        r,j=rbank['caps'][ri],jbank['caps'][ji]
        value=s.dot(np.array(r['normal']),r['profile'],np.array(j['normal']),j['profile']);rel=T.G.relation(value)
        row=dict(index=len(pairs),R=ri,J=ji,native_product=value,relation=rel,score=None);pairs.append(row)
        hist[str(rel)]+=1;digest.update(json.dumps([ri,ji,value,rel],separators=(',',':')).encode())
        if rel in (None,'parallel_or_unresolved'):continue
        nn={**base,67:np.array(r['normal']),68:np.array(j['normal'])};types={67:r['profile'],68:j['profile']}
        for new in (67,68):
            vals,bad,cp=s.native(nn[new],types[new],{i:v for i,v in nn.items() if i!=new},types)
            assert not bad and cp<-1e-8
        sid=s.score(nn,types,ids+[67,68],'joint_subdivision');row['score']=sid
        if sid is None:continue
        ledger=s.ledgers[sid];failed=[f for f in duties if s.corner(f,nn,types) is None]
        faces=[H.face(ledger,i) for i in (8,10,44,55,63,64,67,68)]
        expected={tuple(sorted(f)) for f in topology()['triangles']}
        actual={tuple(sorted(x['active'])) for x in faces[-1]['rays']}
        full_connector=faces[-1]['all_rays_finite'] and actual==expected
        rec=dict(pair=row['index'],score=sid,finite_rays=ledger['counts'].get('finite',0),open_ends=ledger['open_ends'],
            failed_corners=failed,whole_face_records=faces,connector_complete=full_connector,
            compact8_preserved=faces[0]['all_rays_finite'],compact10_restored=faces[1]['all_rays_finite'],
            exterior_supports=[x['active'] for x in ledger['rays'] if x['kind']=='exterior'],
            open_edges=[e for e in ledger['edges'] if e['kind'] in ('one_open','two_open')])
        coherent.append(rec)
        if not failed and full_connector:completed.append(row['index'])
    paths=[Path(__file__),Path(K.__file__),Path(H.__file__),Path(T.__file__),Path(T.H.__file__),Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),
        Path(T.H.A.__file__),Path(T.H.A.C.__file__),T.H.A.INPUT,T.H.A.C.PARENT,T.H.PARENT,T.H.BOUNDARY,Path(T.S.__file__),Path(T.S.C.__file__),
        ROOT/'results/fresh_nonarithmetic_release30_66_hinge_20261017.md',ROOT/'results/fresh_h9_fixed_core_cusp_obstruction_20261013.md']
    return dict(status='BOUNDED_FINITE_COMPACT10_EDGE_SUBDIVISION',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        decision_record='An affine10/44/new30 proposal was rejected before any bank execution: the frozen permanent-cusp obstruction excludes finite-volume completion of an affine null direction orthogonal to retainedH4. This packet instead replaces the finite10--30 incidence by a connector, without deleting8 or forcing10/new30 finite.',
        deleted_walls=[30],retained_old_walls=62,exact_deletion_debit=debit,old_frame_interface=interface,frame_native_checks=checks,
        topology=topology(),control_score=control,deletion_score=deleted,low_labels=[2,3,4,5],banks=[rbank,jbank],selected_R=selected,
        selection_rule='All R bank candidates whose intrinsic10/R Gram has negative determinant and negative product, and whose native44/R relation is finite. The removed10/R target-face edge must be absent; a finite native10/R label alone would not decide this after white projection. Near-zero determinants stay unresolved. All J bank candidates enter. No ranking or truncation.',
        required_corners=duties,pairs=pairs,complete_pair_count=len(pairs),complete_pair_digest=digest.hexdigest(),pair_histogram=dict(hist),coherent_pairs=coherent,
        completed_connector_pairs=completed,boundary_ledgers=s.ledgers,counts=dict(s.counts),
        histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n) for (a,b,c,d),n in sorted(s.hist.items())],unresolved=s.unresolved,
        scope='Changed finite incidence subdivides10--new30; all seven inherited profiles and2..5 frame labels, both norm branches and all62 retained old walls plus mutual native/common-center gates. Complete selected Cartesian pair domain; every coherent final pair has full numericalS6 cone and independentPD support/open-end agreement and eight whole-face records. Exact old13/10/1 certificate debit; no high/continuous completeness, exact new normals, ambient finite census or finite-volume/newdimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'banks',[len(b['caps']) for b in d['banks']],'selected',len(d['selected_R']),'pairs',len(d['pairs']),'completed',len(d['completed_connector_pairs']),flush=True)
