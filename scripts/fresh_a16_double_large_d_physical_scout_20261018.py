#!/usr/bin/env python3
"""Bounded physical/current-interface reconstruction for the proposed A16 rule.

No geometric theorem is inferred from this finite interface certificate.
The accompanying report supplies arbitrary-support and routing arguments.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import fcntl
from itertools import combinations
import json
from pathlib import Path
import signal

import fresh_a16_double_large_d_upper_scout_20261018 as M
import fresh_a16_d16_a4_a12_fractional_transfer_20261018 as L
import fresh_finite_a16_endpoint_sparse_half_transfer_20261018 as E
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

V=M.V;U=M.U;R=M.R;ROOT=M.ROOT;A=V.A;Q=V.Q;D=V.D;T=V.T;P=V.P;B=V.B
TAG='fresh_a16_double_large_d_physical_scout_20261018'
OVER=ROOT/f'results/{A.OVER}.json'
PAIRS={(0,13),(2,15),(1,13),(2,14),(1,12),(3,14),(2,13)}
SELECTED=[1010,1013,1028,1053,1138,1167,1611,2015,2403]
DC=[1035,1905,1964,2421];AC=[1029,1612,1671,2416]
RESTORE={(0,12),(2,14)}


def extract(marker,back):
    """One minified JSON object from a streamed64KiB marker search."""
    with OVER.open('rb')as stream:
        buffer=b''
        while True:
            chunk=stream.read(65536)
            if not chunk:raise ValueError(('marker absent',marker))
            buffer+=chunk;at=buffer.find(marker)
            if at<0:
                buffer=buffer[-len(marker):];continue
            buffer=buffer[at+len(marker)-back:]
            while True:
                try:return json.JSONDecoder().raw_decode(buffer.decode())[0]
                except json.JSONDecodeError:
                    chunk=stream.read(65536)
                    if not chunk:raise
                    buffer+=chunk
                    assert len(buffer)<4*1024**2,'Selected object unexpectedly large'


def native(index):
    marker=('{"index":'+str(index)+',').encode();c=extract(marker,len(marker))
    key=c['core_key'];core=extract((json.dumps(key,separators=(',',':'))+':{').encode(),1)
    g=A.matrix(core['gram']);n=len(g);assert core['gauge']==[1]*n
    assert all(g[i][i]==2 for i in range(n))
    left=c['left'][1];right=[c['right_copy_map'][i]for i in c['right'][1]]
    white=sorted(set(range(n))-set(left+right))
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,white+[b]))==[len(white)+1,0,0]for b in left+right)
    actual=Q.inspect_corner(g,[1]*n,left+right,*c['preferred_corner'])
    for k in('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
    if c['method']=='radical':
        assert n==19 and Q.inertia(g)==[17,1,1]
        z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j]for j in range(n))==0 for i in range(n))
        allowed=c['possible_contacts'];assert all(z[i]for i in allowed)
        du=Q.graph_dist(g,left);dv=Q.graph_dist(g,right);short=[]
        for i,j in combinations(allowed,2):
            if z[i]*z[j]<0:
                length=min(du[i]+dv[j],du[j]+dv[i])+2;assert length<=8
                short.append(dict(nodes=[i,j],length=length))
        assert short
        seal=dict(method='full_radical',radical=z,allowed_contacts=allowed,mixed_sign_shortcuts=short)
    else:
        gg,J,du,dv,end=X.constraints(c,core);assert gg==g and Q.inertia(g)==[17,1,0]
        remaining=sorted(set(c['possible_contacts'])-set(c['marked_component']))
        assert remaining and all(J[i][i]>2 for i in remaining)
        if index==1010:
            assert remaining==[1,11]and[J[i][i]for i in remaining]==[8,8]
        seal=dict(method='whole_selected_component_then_all_contact_projection',
            whole_selected_nodes=c['marked_component'],remaining_contacts=remaining,
            inverse_diagonals=[(i,J[i][i])for i in remaining],inverse=J,end_proof=end)
    corners=[]
    for p in c['all_corners']:
        marks=p['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        corners.append(dict(marks=marks,components=D.finite_blocks(g,ids)))
    return dict(index=index,core_key=key,gram=g,white=white,selected_marks=c['preferred_corner'],
        selected_nodes=c['marked_component'],seal=seal,corners=corners)


def images(g,block,mark):
    return T.mark_pairs(g,block['nodes'],[mark],len(block['nodes']),block['type'][0])


def change(g,white,marks,external_a1=False):
    ids=sorted(set(white)|set(marks))
    if Q.inertia(A.sub(g,ids))!=[len(ids),0,0]:return dict(delta=F(),roles=[],components=[])
    bs=D.finite_blocks(g,ids);co=Counter(b['type']for b in bs)
    a,b=[next(b for b in bs if m in b['nodes'])for m in marks]
    has1=bool(co['A1']or external_a1);roles=[];delta=F();exact_before=None
    if a is b and a['type']=='A16':
        im=T.mark_pairs(g,a['nodes'],marks,16,'A')
        if set(im)&PAIRS:
            assert set(im)<=PAIRS
            old=F(5,8)if set(im)&{(1,12),(3,14)}and has1 else F()if set(im)&{(1,12),(3,14)}else F(1)if im==[(2,13)]else F(1,2)
            delta=-old;exact_before=old
            roles.append(dict(role='two_A16_large_D_debit',images=im,old=old))
    elif a is b and a['type']=='A15':
        im=T.mark_pairs(g,a['nodes'],marks,15,'A')
        if set(im)&RESTORE:
            assert set(im)<=RESTORE
            delta=F(1,2);exact_before=F()
            roles.append(dict(role='A15_half_restore_for_A16_consumer',images=im))
    elif a is not b:
        for sm,la,m,n in((a,b,*marks),(b,a,marks[1],marks[0])):
            st,lt=sm['type'],la['type']
            if(st,lt)not in(('A1','A13'),('A4','A12'),('A2','A14'),('A4','A13'),('A1','D13'),('A4','D13')):continue
            si=images(g,sm,m);li=images(g,la,n);inc=F();role=None
            if(st,lt)==('A1','A13'):
                if li==[(0,),(12,)]:inc=F(1,2);role='A1_A13_end_increment'
                elif li==[(2,),(10,)]and co['A2']:inc=F(1,2);role='A1_A13_third_with_A2'
            elif(st,lt)==('A4','A12')and si==[(0,),(3,)]and li==[(1,),(10,)]and has1:
                inc=F(5,8);role='A4end_A12_second_with_A1'
            elif(st,lt)==('A2','A14')and li==[(2,),(11,)]and has1:
                inc=F(1,2);role='A2_A14_third_with_A1'
            elif(st,lt)==('A4','A13')and si==[(1,),(2,)]and li==[(2,),(10,)]:
                inc=F(1,2);role='A4inner_A13_third'
            elif(st,lt)==('A1','D13')and li==[(11,),(12,)]and co['A2']:
                inc=F(1,2);role='A1_D13_short_with_A2'
            elif(st,lt)==('A4','D13')and si==[(1,),(2,)]and li==[(11,),(12,)]:
                inc=F(1,2);role='A4inner_D13_short'
            if inc:
                delta+=inc;roles.append(dict(role=role,small_images=si,large_images=li))
    assert len(roles)<=1
    return dict(delta=delta,roles=roles,components=bs,exact_current_before=exact_before)


def amended(g,white,marks,previous,external_a1=False):
    ch=change(g,white,marks,external_a1)
    before=F(previous)if ch.get('exact_current_before')is None else ch['exact_current_before']
    pay=before+ch['delta'];assert pay>=0
    return dict(marks=marks,saved_previous=previous,current_before=before,**ch,payment=pay)


def local():
    old=V.load(L.TAG);prior=V.load(P.TAG);endpoint=V.load(E.TAG)
    support=old['local_proofs']['domains']['A']['complete_selected_support']
    hits=[p for p in support['rows']if set(map(tuple,p['marked_images']))&PAIRS]
    assert [p['index']for p in hits]==SELECTED and sum(len(p['maps'])for p in hits)==18
    assert all(len(p['base_nodes'])==16 for p in hits)
    baselines={p['index']:p for p in endpoint['selected_proofs']}
    baselines.update({p['index']:p for p in prior['local_proof']['A16_seals']})
    d16=V.load('fresh_finite_d16_alternate_half_zero_transfer_20261018')
    dproofs={p['index']:p for p in d16['selected_proofs']['proofs']}
    a15=V.load('fresh_finite_a15_three_coupon_zero_transfer_20261018')
    aproofs={p['index']:p for p in a15['selected_proofs']}
    rows=[]
    for index in SELECTED+DC+AC:
        p=native(index);g=p['gram'];white=p['white']
        old=baselines.get(index,dproofs.get(index,aproofs.get(index)))
        previous={tuple(c['marks']):c['payment']for c in old['corners']}if old else{tuple(c['marks']):int(c['marks']==p['selected_marks'])for c in p['corners']}
        variants=[]
        for ext in(False,True):
            cs=[amended(g,white,c['marks'],previous[tuple(c['marks'])],ext)for c in p['corners']]
            total=sum(c['payment']for c in cs);assert total>=1,(index,ext,total)
            variants.append(dict(external_A1=ext,corners=cs,total=total))
        p['variants']=variants;rows.append(p)
    return dict(complete_A16_patterns=394,complete_A16_maps=238,affected_contexts=SELECTED,
        affected_maps=18,selected_and_later_seals=rows,selected_seals=9,D16_consumers=DC,A15_consumers=AC)


def literals(local_proof):
    current=V.load(V.TAG);gram=local_proof['selected_and_later_seals'][1]['gram']
    docs={s:V.load(s)for s in V.O.LITERAL if s!=A.OVER}
    docs[A.OVER]={'isolation':{'cores':{'(25, 25, ((9, 9),))':{'gram':gram}}}}
    out=[];counts=Counter();negative=[]
    for old,(stem,index,g,white,cs)in zip(current['current_literal_interface']['rows'],T.literal_inputs(docs),strict=True):
        assert(old['source'],old['case'])==(stem,V.O.enc(index))
        modes=[]
        for ext in(False,True):
            corners=[]
            for prior,(marks,*_)in zip(old['corners'],cs,strict=True):
                assert prior['marks']==marks
                row=amended(g,white,marks,prior['lower'],ext);corners.append(row)
                if row['delta']<0:
                    assert stem==T.T.NEW_D
                    negative.append(dict(source=stem,case=index,external_A1=ext,marks=marks))
            total=sum(c['payment']for c in corners);assert total>=1,(stem,index,ext,total)
            modes.append(dict(external_A1=ext,corners=corners,total=total))
        out.append(dict(source=stem,case=index,modes=modes))
    assert len(out)==216
    return dict(configurations=216,guarded_variants=432,rows=out,negative_A16_outputs=negative)


def later():
    hits=[];packets=[];count=0
    for stem,keys in B.FINITE.items():
        doc=V.load(stem);n=0
        for sec in keys:
            for path,p,g,w,bs in B.corners(doc[sec],sec):
                b=next((b for b in bs if b['type']=='A16'and set(p['marks'])<=set(b['nodes'])),None)
                if b:
                    im=T.mark_pairs(A.matrix(g),b['nodes'],p['marks'],16,'A')
                    if set(im)&PAIRS:hits.append(dict(source=stem,path=path,images=im))
                n+=1;count+=1
        packets.append(dict(source=stem,corners=n))
    assert count==564 and len(hits)==20
    assert set(p['source']for p in hits)=={E.TAG,'fresh_rank16_existing_payment_exchange_20261018','fresh_finite_d16_alternate_half_zero_transfer_20261018','fresh_finite_a15_three_coupon_zero_transfer_20261018'}
    # Newer current seals are checked in both external-A1 states. Only the
    # two own1167/2015 contexts may contain a changed positive A16 role.
    newer=[];own=[]
    collections=[(P.R.TAG,('local_proof','whole_proofs')),(P.H.TAG,('local_proof','whole_proofs')),
        (P.S.TAG,('local_proof','whole_proofs')),(P.N.TAG,('local_proof','proofs')),
        (P.R.R.TAG,('local_seals',)),(P.H.Z.U.TAG,('local_proof',)),
        (P.TAG,('local_proof','old_ideal17_seals')),(P.TAG,('local_proof','A16_seals'))]
    for stem,keys in collections:
        doc=V.load(stem);xs=doc
        for k in keys:xs=xs[k]
        if isinstance(xs,dict):xs=[xs]
        rows=[]
        for p in xs:
            physical=p.get('old_complete_seal',p)
            gg=physical.get('gram',physical.get('old_seal',{}).get('gram'))
            ww=physical.get('white',physical.get('old_seal',{}).get('white'))
            if gg is None:continue
            g=A.matrix(gg)
            for c in p.get('corners',[]):
                ch=change(g,ww,c['marks'],True)
                if ch['delta']<0:
                    assert p.get('index')in(1167,2015)
                    own.append(dict(source=stem,index=p['index'],marks=c['marks']))
                rows.append(dict(index=p.get('index',p.get('name')),marks=c['marks'],roles=ch['roles']))
        newer.append(dict(source=stem,path=keys,rows=rows))
    chosen=[]
    for stem in(B.PART,B.REST):
        doc=V.load(stem);rs=[(path,c)for path,c in P.W.walk(doc)if F(c['payment'])>0]
        assert len(rs)==(280 if stem==B.PART else 52)
        assert all(not any(b['type']=='A16'for b in c['roles'])for _,c in rs)
        chosen.append(dict(source=stem,positive_outputs=len(rs),A16_outputs=0))
    return dict(later_packets=packets,later_corners=count,target_hits=hits,newer_seals=newer,
        replaced_own_newer_outputs=own,chosen_D12=chosen)


def run():
    lp=local();literal=literals(lp);interfaces=later();upper=V.load(M.TAG)
    assert upper['complete_classes']==16 and F(upper['increasing_family_upper'])==371
    paths={Path(__file__),Path(M.__file__),Path(V.__file__),Path(L.__file__),Path(E.__file__),Path(X.__file__),OVER}
    paths.update(ROOT/f'results/{s}.json'for s in{M.TAG,V.TAG,L.TAG,E.TAG,P.TAG,*V.O.LITERAL,*B.FINITE,B.PART,B.REST})
    return V.O.enc(dict(status='PASS_A16_DOUBLE_LARGE_D_PHYSICAL_SCOUT',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(paths)},
        local=lp,current_literal_interface=literal,current_later_interfaces=interfaces,
        arithmetic=upper['witnesses'],complete_increasing_family_upper='371',
        universal_geometry_in_written_report=True,independent_review_complete=False,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);R.resources();U.memory_guard()
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=run();U.memory_guard();R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['local']['selected_and_later_seals']),len(out['current_literal_interface']['negative_A16_outputs']))
