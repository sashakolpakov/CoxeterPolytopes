#!/usr/bin/env python3
"""C16 two free zeros and four D16 sixteenth savings with one small coupon.

Complete changed native support, arbitrary-contact inherited completions,
whole-factor roles and simultaneous current literal/ideal interfaces.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_order16_bd_sparse_unit_transfer_20261018 as K
import fresh_order17_center_contextual_zero_transfer_20261018 as Z

ROOT=K.ROOT;D=K.D;X=K.X;P=K.P;Q=K.Q;L=K.L;C=Z.C;T=C.T;A=C.A
TAG='fresh_order16_cd_sparse_sixteenth_transfer_20261018'
CHANGES={'C':{(2,12),(3,13)},'D':{(1,14),(1,15),(10,12),(10,13)}}
DEXPECT=[1036,1050,1647,1667,1906,1934,1975,2299,2327,2368,2505,2512]
SPARSE='fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018'
STEMS=tuple(dict.fromkeys((*Z.STEMS,Z.TAG,K.TAG,*K.STEMS,SPARSE)))
FROZEN_Z_MD='7d4e96b2d7be94925208016f00f16c3d8f00d3b2f2c918ea168cf32f5a2bf021'
enc=Z.enc


def edge_checks(c,g):
    black=set(c['left'][1])|{c['right_copy_map'][i] for i in c['right'][1]};white=sorted(set(range(len(g)))-black)
    assert len(black)==4 and Q.inertia(L.sub(g,white))==[len(white),0,0]
    for b in black:assert Q.inertia(L.sub(g,white+[b]))==[len(white)+1,0,0]
    return dict(white=white,black=sorted(black),white_positive=True,four_edge_positive=True)


def c_case(c,core,target,maps):
    g,J,du,dv,ends=X.constraints(c,core);assert c['method']=='lorentz_residual'
    edges=edge_checks(c,g);rows=[]
    for mp in maps:
        image={t:c['marked_component'][i] for i,t in enumerate(mp)};entries=[]
        for new in sorted(set(range(16))-set(image)):
            for old,i in sorted(image.items()):
                label=target[new][old]
                if not label:continue
                assert label in (3,4) and i in c['possible_contacts']
                entries.append(dict(new_target=new,old_target=old,old_core_node=i,label=label,
                    physical_contact_square=1 if label==3 else 2,projection_lower_bound=(1 if label==3 else 2)*J[i][i]))
        witnesses=[r for r in entries if r['projection_lower_bound']>2]
        assert len(witnesses)==1 and witnesses[0]['label']==4 and witnesses[0]['projection_lower_bound']==F(8,3)
        rows.append(dict(embedding=mp,required_entries=entries,witness=witnesses[0],outcome='forced_doubled_terminal_exceeds_norm2'))
    return dict(case=c['index'],core_key=c['core_key'],gram=core['gram'],inverse=core['inverse'],
        edge_checks=edges,arbitrary_contact_end_proof=ends,rows=rows)


def new_coupon(g,blocks,marks):
    a,b=[next(b for b in blocks if m in b['nodes']) for m in marks]
    if a is b:return []
    out=[]
    for sm,la,m,n in ((a,b,*marks),(b,a,marks[1],marks[0])):
        if sm['type']=='A1' and la['type']=='D15':
            images=T.mark_pairs(g,la['nodes'],[n],15,'D')
            if images==[(11,)]:out.append(dict(small=sm,large=la,marks=[m,n],images=images,weight='1/16'))
    assert len(out)<=1;return out


def d_payments(index,proof,g):
    assert index in (1975,2368,2505,2512) and proof['outcome']=='sealed_completion'
    assert not proof['new_roots'] and proof['signature']==[17,1,0]
    corners=[]
    for row in proof['corners']:
        if row['kind']=='affineD15_order16':
            corners.append(dict(marks=row['marks'],kind=row['kind'],payment=F(7,16)));continue
        marks=row['marks'];blocks=D.finite_blocks(g,sorted(set(proof['white'])|set(marks)))
        assert blocks==row['components'];new=new_coupon(g,blocks,marks);old=[]
        a,b=[next(b for b in blocks if m in b['nodes']) for m in marks]
        if a is not b:
            for sm,la,m,n in ((a,b,*marks),(b,a,marks[1],marks[0])):
                if sm['type']=='A4' and la['type'] in ('A12','D12'):
                    si=T.mark_pairs(g,sm['nodes'],[m],4,'A');li=T.mark_pairs(g,la['nodes'],[n],12,la['type'][0])
                    if si==[(0,),(3,)] and li==([(0,),(11,)] if la['type']=='A12' else [(10,),(11,)]):
                        old.append(dict(small=sm,large=la,marks=[m,n],images=[si,li],payment='1' if la['type']=='A12' else '1/2'))
        assert len(old)<=1 and not(old and new)
        pay=sum((F(p['payment']) for p in old),F())+F(len(new),16)
        corners.append(dict(marks=marks,kind='finite',whole_components=blocks,existing_coupon=old,new_coupon=new,
            payment=pay,uses_finite_same_component_payment=False))
    total=sum(r['payment'] for r in corners);assert total>=1
    if index in (2368,2512):assert total==1 and sum(bool(r.get('new_coupon')) for r in corners)==1
    else:assert total==F(31,16) and any(r['payment']==1 for r in corners)
    return dict(case=index,corners=corners,total=total)


def d_case(c,core,target,maps,old):
    g,con=D.core_constraints(c,core);assert c['method']=='lorentz_residual' and con==old['constraints']
    edges=edge_checks(c,g);rows=[];payments=[]
    for mp in maps:
        proof,full=D.complete(c,core,g,con,target,mp)
        previous=next(r for r in old['completions'] if any(e['map']==mp for e in r['embeddings']))
        # The predecessor groups maps with identical completed Grams; its
        # representative target ordering can permute the same fork leaves.
        ordered={'target_to_literal','affine_nodes','positive_affine_kernel'}
        assert {k:v for k,v in proof.items() if k not in ordered}=={k:v for k,v in previous.items() if k not in ordered|{'embeddings'}}
        if 'affine_nodes' in proof:
            assert dict(zip(proof['affine_nodes'],proof['positive_affine_kernel']))==dict(zip(previous['affine_nodes'],previous['positive_affine_kernel']))
        if proof['outcome']=='sealed_completion':payments.append(d_payments(c['index'],proof,full))
        rows.append(dict(embedding=mp,proof=proof))
    return dict(case=c['index'],core_key=c['core_key'],gram=core['gram'],constraints=con,edge_checks=edges,
        completions=rows,current_payments=payments)


def local(docs):
    small=docs[D.SMALL];over=docs[D.OVER]['isolation'];patterns=[p for p in small['patterns'] if p['order_bound']>=16 and len(p['nodes'])<=16]
    assert len(patterns)==394;oldcases={r['index']:r for r in docs[D.TAG]['cases']};rules=[]
    for kind in ('C','D'):
        rec=next(r for r in small['records'] if (r['kind'],r['order'])==(kind,16));target=rec['target']
        if kind=='C':ws={p:F(p in set(map(tuple,rec['retained_paid_pairs']))) for p in combinations(range(16),2)}
        else:ws={tuple(r['pair']):F(r['weight']) for r in next(r for r in docs[K.TAG]['records'] if r['kind']=='D')['weighted_pairs']}
        assert sum(ws.values())==103;before=dict(ws)
        for p in CHANGES[kind]:assert ws[p]==(1 if kind=='C' else F(1,2));ws[p]=F() if kind=='C' else F(7,16)
        assert sum(ws.values())==(101 if kind=='C' else F(411,4))
        autos=P.embeddings(target,target);_,dist,features=P.distance_data(target);basic=L.S.basic_pairs(target)
        assert not basic&CHANGES[kind] and all(ws[p]==before[p] for p in basic)
        assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGES[kind]}==CHANGES[kind] for mp in autos)
        assert all(dist[i][j]==(10 if kind=='C' else 11) for i,j in CHANGES[kind])
        assert all(ws[p]==1 for p in ws if dist[p[0]][p[1]]<=8)
        triangles={p for p in ws if dist[p[0]][p[1]]<=5 or (dist[p[0]][p[1]]<=8 and any(dist[i][f]<=5 for i in p for f in features))}
        assert all(ws[p]==1 for p in triangles)
        saved={r['case']:r for r in rec['case_embeddings']};support=[];cases=[];outcomes=Counter()
        for p in patterns:
            maps=P.embeddings(p['labels'],target);images=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
            if maps:assert maps==saved[p['index']]['embeddings'] and images==list(map(tuple,saved[p['index']]['marked_images']))
            else:assert p['index'] not in saved
            support.append(dict(case=p['index'],embeddings=maps,marked_images=images))
            changed=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks'])) in CHANGES[kind]]
            if not changed:continue
            c=over['cases'][p['index']];core=over['cores'][c['core_key']]
            proof=c_case(c,core,target,changed) if kind=='C' else d_case(c,core,target,changed,oldcases[c['index']])
            cases.append(proof)
            if kind=='C':outcomes.update(r['outcome'] for r in proof['rows'])
            else:outcomes.update(r['proof']['outcome'] for r in proof['completions'])
        assert [r['case'] for r in cases]==([1031,1614] if kind=='C' else DEXPECT)
        assert outcomes==(dict(forced_doubled_terminal_exceeds_norm2=4) if kind=='C' else
            dict(excluded_projection_above_norm=32,excluded_lorentz_index=16,sealed_completion=32))
        rules.append(dict(kind=kind,target=target,changed_pairs=sorted(CHANGES[kind]),old_cost='103',new_cost=str(sum(ws.values())),
            weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(ws.items())],automorphisms=autos,
            complete_patterns=394,successful_patterns=len(saved),complete_maps=sum(len(r['embeddings']) for r in support),
            support=support,cases=cases,affected_maps=sum(outcomes.values()),outcomes=dict(outcomes),
            retained_basic_pairs=sorted(basic),retained_triangle_pairs=sorted(triangles)))
    return rules


def literal(docs,rules):
    old={C.R.key(r['source'],r['case']):r for r in docs[Z.TAG]['current_literal_interface']['rows']};rows=[];hits=Counter()
    targets={r['kind']:r['target'] for r in rules}
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prior=old[C.R.key(stem,index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==prior['corners'][j]['marks'];ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids));debit=F();cp=[];images={}
            if sig==[len(ids),0,0]:cp=new_coupon(g,D.finite_blocks(g,ids),marks)
            elif sig[1]==0:
                b=next((b for b in A.components(g,ids) if set(marks)<=set(b)),None)
                if b and len(b)==16:
                    labels=T.T.labels(g,b)
                    if labels is not None:
                        for kind,target in targets.items():
                            ps=sorted({tuple(sorted(mp[b.index(m)] for m in marks)) for mp in P.embeddings(labels,target)})
                            if ps:
                                images[kind]=ps
                                if set(ps)&CHANGES[kind]:assert set(ps)<=CHANGES[kind];debit=F(1) if kind=='C' else F(1,16)
            inc=F(len(cp),16);pay=F(prior['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert pay>=0
            if debit:hits['ideal_debits']+=1
            if inc:hits['finite_credits']+=1
            corners.append(dict(marks=marks,images=images,debit=debit,increment=inc,new_coupon=cp,retained_payment_lower_bound=pay))
        total=sum(r['retained_payment_lower_bound'] for r in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(rows)==216
    return dict(configurations=216,changes=dict(hits),rows=rows,minimum_total=min(r['total_lower_bound'] for r in rows))


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert sha256((ROOT/f'results/{Z.TAG}.md').read_bytes()).hexdigest()==FROZEN_Z_MD
    assert docs[K.TAG]['new_order16_family_costs']==dict(A=80,B=102,C=103,D=103)
    assert any(r['name']=='A4end_A12end' and r['weight']=='1' for r in docs[SPARSE]['coupons'])
    rules=local(docs);types=[t for t in C.P.Z.B.B.types() if t[1]<=59];rep=C.P.Z.B.B.representative
    for a in types:
        for b in types:
            old13=8 if {a,b}=={('A',1),('B',12)} else 0
            inc=int({a,b}=={('A',1),('D',15)})
            upperinc=int({rep(a),rep(b)}=={('A',1),('D',15)})
            assert C.actual_cross(a,b)+old13+inc<=C.upper_cross(rep(a),rep(b))+upperinc
    assert C.actual_cross(('A',1),('D',15))==C.upper_cross(('A',1),('D',15))==24
    rows=[]
    for r in docs[Z.TAG]['ideal']['rows']:
        rs=tuple(r['ranks']);saving=F(rs.count(15),64);m=F(r['normalized_margin'])+saving
        rows.append(dict(ranks=rs,link_vertices=r['link_vertices'],previous_normalized_margin=r['normalized_margin'],
            normalized_saving=saving,normalized_margin=m,absolute_margin=m*r['link_vertices']))
    summary=C.U.summary(rows);assert summary['count']==1890 and summary['negative_count']==summary['zero_count']==0
    assert summary['normalized_minimum']==F(1,340) and summary['absolute_minimum']==F(91,2)
    paths={Path(__file__),Path(K.__file__),Path(Z.__file__),Path(D.__file__),Path(X.__file__),Path(P.__file__),Path(Q.__file__),Path(L.__file__),Path(T.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_ORDER16_C101_D411_OVER4_SPARSE_SIXTEENTH_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        rules=rules,new_family_costs=dict(A='80',B='102',C='101',D='411/4'),new_order16_maximum='411/4',new_rank15_envelope='411/64',
        new_coupon=dict(small='whole_A1',large='whole_D15',large_marks=[11],large_role='long_arm_neighbor_of_fork_center',
            old_weight='0',new_weight='1/16',component_pair_increment='1/16',old_actual_component_cost='3/2',new_actual_component_cost='25/16'),
        finite_accounting=dict(simultaneous_type_pairs=len(types)**2,representative_A1_D15_allowance='25/16',previous_allowance='3/2',
            scale16_item1_increment='fixed_rank15_count',forest_count_proved=False,current_contextual17_cost_unchanged=True),
        current_literal_interface=literal(docs,rules),ideal_frontier_update=dict(rows=rows,summary=summary,complete_saved_domain=1890,new_convolution=False),
        inherited_routing=docs[K.TAG]['complete_high_overlap_routing'],triangle_minimum=2,quadrilateral_minimum=1,
        pending_finite6_rule_assumed=False,completed_unrestricted_bound=59,new_dimension_bound=None,
        scope='Four forced C terminal exclusions,80 complete D completions and four whole-factor payments; one ordinary A1/D15 sixteenth. Full current literal and terminal ideal update; new finite forest budget remains separate.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_family_costs'],out['current_literal_interface']['changes'],flush=True)
