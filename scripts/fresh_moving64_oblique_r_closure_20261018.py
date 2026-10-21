#!/usr/bin/env python3
"""Change replacement30's actual duties after moving64.

All63 frozen J/L/X/M prefixes are reconstructed. R70 uses10/J/X with
only10/J/R andJ/X/R duties; old33 remains a tested wall but33/R may be
dotted. Seven-profile low cubes, numerical full-cone/PD and whole faces.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_joint64_native_compact_closure_20261018 as C

ROOT=C.ROOT;T=C.T;Q=C.Q;H=C.H
TAG='fresh_moving64_oblique_r_closure_20261018'
DELETED=C.DELETED


def run():
    old=json.loads((ROOT/f'results/{C.TAG}.json').read_text())
    flip=json.loads((ROOT/f'results/{Q.TAG}.json').read_text())
    seed=json.loads((ROOT/f'results/{Q.P.TAG}.json').read_text())
    for doc in (old,flip,seed):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert old['summary']['M_native']==63 and old['summary']['native_five_caps']==0
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s);assert interface==old['old_frame_interface']
    assert all(i not in s.index for i in (67,68,69,70,71))
    debit=T.H.deletion(s.primary,s.cur,DELETED)
    assert json.loads(json.dumps(debit))==old['exact_deletion_debit']
    base={i:v for i,v in s.old.items() if i not in DELETED}
    ids=[i for i in s.cur['actual_sides'] if i not in DELETED]
    score=lambda nd,tt,ii,stage:Q.P.B.B.score(s,nd,tt,ii,stage,[i for i in C.WALLS if i in nd])
    control=score(s.old,{},s.cur['actual_sides'],'original63_control')
    assert control['score'] is not None and control['open_ends']==10 and control['finite_rays']==26
    deletion=score(base,{},ids,'delete30_44_63_64_65')
    top=C.topology();records=[];accepted=[];reverse_candidates={};duplicates8=[]
    def completion(r,reference):
        full,ft,rec=reference
        cs=[s.corner(f,full,ft) for f in top['triangles']]
        exterior=[[7,55,70],[55,59,70],[10,33,70]]
        ec=[s.corner(f,full,ft) for f in exterior]
        faces={f['wall']:f for f in r['whole_faces']}
        compact={str(w):faces[w]['all_rays_finite'] for w in (8,10,11,33,52,55,66,67,68,69,70,71)}
        capsets={str(w):faces[w]['all_rays_finite'] and sorted(sorted(q['active']) for q in faces[w]['rays'])==sorted(sorted(f) for f in top['triangles'] if w in f) for w in (67,68,69,71)}
        return dict(**rec,corner_types=cs,all_disk_corners=all(c is not None for c in cs),
            diagnostic_external_duties=exterior,diagnostic_external_corner_types=ec,
            old33_product=s.dot(full[70],ft[70],full[33],s.ps(33,ft)),
            whole_compact_faces=compact,exact_whole_new_caps=capsets,whole_J_closed=capsets['67'],
            clean=r['score'] is not None,
            **{k:r[k] for k in ('score','finite_rays','open_ends','exterior_supports','open_edges','whole_faces','cone_diagnostic')})
    for prefix in old['records']:
        sel=prefix['selection'];fr=flip['records'][sel['prefix']]
        j=seed['J_caps'][fr['J']];l=seed['L_bank']['caps'][fr['L']];x=fr['bank']['caps'][sel['X']]
        jv,jcheck=C.solve_saved(s,[8,10,63],0,[2,3,2],j['branch'],s.old,{},j['normal'])
        lv,lcheck=C.solve_saved(s,[8,52,64],l['profile'],l['labels'],l['branch'],s.old,{},l['normal'])
        nd={**base,67:jv,68:lv};tt={67:0,68:l['profile']}
        xv,xcheck=C.solve_saved(s,[8,67,68],x['profile'],x['labels'],x['branch'],nd,tt,x['normal'])
        nd[69]=xv;tt[69]=x['profile']
        for m in prefix['bank']['caps']:
            mv,mcheck=C.solve_saved(s,[52,68,69],m['profile'],m['labels'],m['branch'],nd,tt,m['normal'])
            nn={**nd,71:mv};types={**tt,71:m['profile']}
            for i in types:
                _,bad,cp=s.native(nn[i],types[i],{k:v for k,v in nn.items() if k!=i},types)
                assert not bad and cp<-1e-8
            before=score(nn,types,ids+[67,68,69,71],'frozen_moving64_prefix')
            assert before['score'] is not None and before['finite_rays']==19 and before['open_ends']==13
            rb=C.bank(s,[10,67,69],[(0,1),(1,2)],nn,types,ids+[67,68,69,71],70,'oblique_new30')
            row=dict(index=len(records),prefix=prefix['index'],M=m['index'],selection=sel,
                reconstruction=[jcheck,lcheck,xcheck,mcheck],before=before,bank=rb)
            records.append(row)
            for q in rb['recipes']:
                for b in q['branches']:
                    bad=[i for i,_ in b['blockers']]
                    if bad==[8]:
                        difference=float(np.max(abs(np.array(b['normal'])-s.old[8])))
                        duplicates8.append(dict(record=row['index'],recipe=q['index'],branch=b['branch'],
                            profile=q['profile'],old8_profile=s.ps(8,{}),max_projected_difference=difference,
                            same_full_profile=q['profile']==s.ps(8,{}),
                            duplicates_old8=q['profile']==s.ps(8,{}) and difference<2e-6))
                    if bad==[71] and b['center_product']<-1e-8:
                        key=(prefix['index'],q['profile'],*np.round(b['normal'],8))
                        if key not in reverse_candidates:
                            reverse_candidates[key]=dict(record=row['index'],recipe=q['index'],branch=b['branch'],
                                prefix=prefix['index'],profile=q['profile'],labels=q['labels'],normal=b['normal'],
                                removed_M_product=dict(b['native_products'])[71],
                                nd={**nd,70:np.array(b['normal'])},types={**tt,70:q['profile']})
            for r in rb['caps']:
                full={**nn,70:np.array(r['normal'])};ft={**types,70:r['profile']}
                accepted.append(completion(r,(full,ft,dict(record=row['index'],R=r['index']))))
        print('PREFIX',prefix['index'],'M',len(prefix['bank']['caps']),'cumulative_R',len(accepted),flush=True)
    assert len(records)==63
    reverse=[];reverse_accepted=[]
    for candidate in reverse_candidates.values():
        nd=candidate['nd'];tt=candidate['types']
        for i in tt:
            _,bad,cp=s.native(nd[i],tt[i],{k:v for k,v in nd.items() if k!=i},tt)
            assert not bad and cp<-1e-8
        before=score(nd,tt,ids+[67,68,69,70],'native_R_without_previous_M')
        bank=C.bank(s,[68,69,70],[(0,1),(1,2)],nd,tt,ids+[67,68,69,70],71,'reverse_M_closure')
        row=dict(index=len(reverse),selection={k:v for k,v in candidate.items() if k not in ('nd','types')},before=before,bank=bank)
        reverse.append(row)
        for m in bank['caps']:
            full={**nd,71:np.array(m['normal'])};ft={**tt,71:m['profile']}
            reverse_accepted.append(completion(m,(full,ft,dict(reverse=row['index'],M=m['index']))))
        print('REVERSE',row['index'],'recipes',len(bank['recipes']),'nativeM',len(bank['caps']),flush=True)
    common=set(base)
    for r in records:
        for q in r['bank']['recipes']:
            for b in q['branches']:common&={i for i,_ in b['blockers']}
    summary=dict(frozen_selected_JLX_prefixes=9,reconstructed_M_prefixes=len(records),
        R_recipes=sum(len(r['bank']['recipes']) for r in records),
        R_branches=sum(sum(len(q['branches']) for q in r['bank']['recipes']) for r in records),
        native_R=len(accepted),clean_R=sum(a['clean'] for a in accepted),
        whole_J_closed=sum(a['clean'] and a['whole_J_closed'] for a in accepted),
        completed_disks=sum(a['clean'] and a['all_disk_corners'] and all(a['exact_whole_new_caps'].values()) for a in accepted),
        fewer_than_ten_ends=sum(a['clean'] and a['open_ends']<10 for a in accepted),
        common_retained_old_blockers=sorted(common),sole_old8_branches=len(duplicates8),old8_duplicates=sum(r['duplicates_old8'] for r in duplicates8),
        reverse_prefixes=len(reverse),reverse_M_recipes=sum(len(r['bank']['recipes']) for r in reverse),
        reverse_M_branches=sum(sum(len(q['branches']) for q in r['bank']['recipes']) for r in reverse),
        reverse_native_M=len(reverse_accepted),reverse_clean_M=sum(a['clean'] for a in reverse_accepted),
        reverse_completed_disks=sum(a['clean'] and a['all_disk_corners'] and all(a['exact_whole_new_caps'].values()) for a in reverse_accepted),
        reverse_fewer_than_ten_ends=sum(a['clean'] and a['open_ends']<10 for a in reverse_accepted),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(C.__file__),Path(Q.__file__),Path(H.__file__),Path(T.__file__),
        Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),Path(T.H.__file__)}
    paths.update(ROOT/p for p in old['sources'])
    paths.update(ROOT/f'{d}/{tag}.{e}' for tag in (C.TAG,Q.TAG,Q.P.TAG) for d,e in (('scripts','py'),('results','json'),('results','md')))
    return dict(status='BOUNDED_MOVING64_OBLIQUE_R_CLOSURE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        deleted_walls=list(DELETED),retained_old_walls=len(base)+6,replacement64_id=71,replacement30_id=70,
        exact_deletion_debit=debit,old_frame_interface=interface,topology=top,control=control,deletion=deletion,
        prior_selected=old['selected'],records=records,accepted=accepted,old8_duplicate_checks=duplicates8,
        reverse_records=reverse,reverse_accepted=reverse_accepted,summary=summary,
        boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        numerical_snapshot=dict(metric=T.M.tolist(),center=s.center.tolist(),old_normals={str(i):v.tolist() for i,v in s.old.items()},labels=s.labels,profiles=s.profiles,actual_sides=s.cur['actual_sides']),
        scope='All63 frozen native moving64 prefixes over nine selected J/L/X geometries; changedR70 frame10/J/X requires only10/J/R andJ/X/R. ReverseMstage takes every distinctprefix/Rgeometry blockedsolelybyM71, removesoldM andusesL/X/R withL/X/M,X/R/Mduties, thenchecksallremainingdiskcorners. Complete seven-profile labels2..5 banks with both real norm branches; all58 retainedoldwalls andmutualnative/strictcenter before fullnumericalcone/PD/wholeface records. Old33/R maybedotted; no10/XorL/Redge imposed. Exact originalcertificate debit19/12/1. No high-label/allprefix/continuous completeness, exactnewnormals or ambientH9 finitevolume certificate.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists(),'Certificate exists; use --verify.'
    result=json.loads(json.dumps(run()))
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['summary'],flush=True)
