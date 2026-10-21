#!/usr/bin/env python3
"""Whole A16/D16 fractional transfers with distinct A4/A12 coupons.

Complete saved selected domains, arbitrary-support local seals, and current
output interfaces. No forest optimization or new dimension bound.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import fresh_quarter_center_rank13_ambient58_frontier_20261018 as C
import fresh_a1_d5_witness_context_retention_20261018 as W
import fresh_finite_d16_fork_exchange_20261018 as H

ROOT=C.ROOT;V=H.V;T=H.T;A=H.A;Q=H.Q;D=H.D;P=H.P;B=W.B
TAG='fresh_a16_d16_a4_a12_fractional_transfer_20261018'
ACHANGE={(1,12),(3,14)};DCHANGE={(2,12)}
ACASES=[1167,2015];DWHOLE=[1060,2436,2489,2515];DPROPER=[1011,1031,1614]
DW='fresh_finite_d16_alternate_half_zero_transfer_20261018'
LEAF='fresh_affine_d11_leaf_pair_transfer_20261017'
STEMS=tuple(dict.fromkeys((C.TAG,W.TAG,H.TAG,DW,LEAF,C.Q.TAG,C.R.TAG,
    *W.STEMS,*H.STEMS,*B.STEMS)))
enc=W.enc

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def pair(a,b):return frozenset((a,b))==frozenset((('A',4),('A',12)))
def actual_base(t):return C.actual_base(t)-({('A',16):48,('D',16):16}.get(t,0))
def upper_base(t):return C.upper_base(t)-(48 if t==('D',16)else 0)
def actual_cross(a,b):return W.actual_cross(a,b)+(160 if pair(a,b)else 0)
def upper_cross(a,b):return C.upper_cross(a,b)+(160 if pair(a,b)else 0)
context=W.context

def changes(g,white,marks):
    ids=sorted(set(white)|set(marks))
    if Q.inertia(A.sub(g,ids))!=[len(ids),0,0]:return dict(delta=F(),roles=[],components=[])
    blocks=D.finite_blocks(g,ids);aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
    delta=F();roles=[]
    if aa is bb and aa['type']in('A16','D16'):
        kind=aa['type'][0];images=T.mark_pairs(g,aa['nodes'],marks,16,kind)
        changed=ACHANGE if kind=='A'else DCHANGE
        if set(images)&changed:
            assert set(images)<=changed
            delta=-F(3,8)if kind=='A'else-F(1,4)
            roles.append(dict(role=aa['type']+'_debit',images=images,delta=delta))
    elif aa is not bb:
        for sm,la,m,lm in((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
            if(sm['type'],la['type'])!=('A4','A12'):continue
            si=T.mark_pairs(g,sm['nodes'],[m],4,'A');li=T.mark_pairs(g,la['nodes'],[lm],12,'A')
            if si!=[(0,),(3,)]:continue
            if li==[(1,),(10,)]:delta=F(3,8);role='A4end_A12second'
            elif li==[(2,),(9,)]:delta=F(1,4);role='A4end_A12third_upgrade'
            else:continue
            roles.append(dict(role=role,small_images=si,large_images=li,delta=delta))
    assert len(roles)<=1
    return dict(delta=delta,roles=roles,components=blocks)

def a_seal(c,core):
    g=A.matrix(core['gram']);n=len(g)
    assert n==19 and core['gauge']==[1]*n and all(g[i][i]==2 for i in range(n))
    assert c['method']=='radical'and Q.inertia(g)==[17,1,1]
    left=c['left'][1];right=[c['right_copy_map'][i]for i in c['right'][1]]
    white=sorted(set(range(n))-set(left+right));assert len(white)==15
    assert Q.inertia(A.sub(g,white))==[15,0,0]
    assert all(Q.inertia(A.sub(g,white+[b]))==[16,0,0]for b in left+right)
    checked=Q.inspect_corner(g,[1]*n,left+right,*c['preferred_corner'])
    for k in('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert checked[k]==c[k]
    images=T.mark_pairs(g,c['marked_component'],c['preferred_corner'],16,'A')
    assert len(c['marked_component'])==16 and set(images)==ACHANGE
    z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j]for j in range(n))==0 for i in range(n))
    allowed=c['possible_contacts'];assert all(z[i]for i in allowed)
    du,dv=Q.graph_dist(g,left),Q.graph_dist(g,right);short=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],dv[i]+du[j])+2;assert length<=6
            short.append(dict(nodes=[i,j],length=length))
    assert short
    rows=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ch=changes(g,white,marks)
        assert ch['components']
        val=F();reason='unused'
        if marks==c['preferred_corner']:
            assert ch['delta']==-F(3,8);val=F(5,8);reason='selected_A16_five_eighths'
        elif ch['delta']==F(3,8):val=F(3,8);reason='distinct_whole_A4end_A12second'
        rows.append(dict(marks=marks,**ch,payment=val,reason=reason))
    assert sorted(r['payment']for r in rows)==[F(),F(),F(3,8),F(5,8)]
    return dict(index=c['index'],core_key=c['core_key'],gram=g,white=white,radical=z,
        allowed_contacts=allowed,mixed_sign_shortcuts=short,corners=rows,total='1',
        arbitrary_support_sealed_unless_independent_short6=True)

def selected(docs):
    over=docs[A.OVER];rules={r['kind']:r for r in docs[H.S.TAG]['rules']if r['rank']==16}
    partition=Counter();eligible=[]
    for c in over['isolation']['cases']:
        bound=c.get('component_order_bound')
        if bound is None:partition[c['method']]+=1
        elif bound<16:partition['bound_below16']+=1
        elif len(c['marked_component'])>16:partition['display_above16']+=1
        else:partition['eligible']+=1;eligible.append(c['index'])
    assert partition==dict(bound_below16=2045,display_above16=92,eligible=394,
        finite_leaf9_or_ideal_feature5=74,nine_context_adaptive_finite_payment=9)
    out={};proofs=[];proper=[]
    for kind,changed,expected,total in(('A',ACHANGE,ACASES,238),('D',DCHANGE,sorted(DWHOLE+DPROPER),578)):
        sup=V.support(kind,docs[T.TAG]);target=rules[kind]['target']
        assert eligible==[r['index']for r in sup['rows']]
        assert V.validate_maps(sup,target,over)==total
        affected=[r for r in sup['rows']if set(map(tuple,r['marked_images']))&changed]
        assert [r['index']for r in affected]==expected
        for r in affected:
            c=over['isolation']['cases'][r['index']];core=over['isolation']['cores'][c['core_key']]
            if kind=='A':proofs.append(a_seal(c,core))
            elif c['index']in DPROPER:
                p=H.proper(c,core,r['maps'],target)
                assert enc(p)==next(x for x in docs[H.TAG]['proper_completion_proofs']if x['index']==c['index'])
                proper.append(p)
            else:
                old=H.seal(c,core)
                assert enc(old)==next(x for x in docs[H.TAG]['selected_proofs']if x['index']==c['index'])
                g=A.matrix(old['gram']);rows=[]
                for q in old['corners']:
                    ch=changes(g,old['white'],q['marks']);new=q['payment']+ch['delta'];assert new>=0
                    rows.append(dict(marks=q['marks'],previous=q['payment'],**ch,payment=new,reason=q['reason']))
                assert sorted(r['payment']for r in rows)==[F(),F(),F(1,4),F(3,4)]
                proofs.append(dict(index=c['index'],old_seal=old,gram=g,white=old['white'],corners=rows,total='1'))
        _,dist,_=P.distance_data(target)
        assert all(dist[i][j]==(11 if kind=='A'else 10)for i,j in changed)
        assert all({tuple(sorted(m[i]for i in p))for p in changed}==changed for m in rules[kind]['automorphisms'])
        assert not changed&set(map(tuple,rules[kind]['retained_basic_pairs']))
        out[kind]=dict(complete_selected_support=sup,affected_contexts=expected,
            affected_maps=sum(sum(tuple(m['marked_images'])in changed for m in r['maps'])for r in affected))
    return dict(complete_high_partition=dict(partition),domains=out,whole_proofs=proofs,proper_completion_proofs=proper)

def interfaces(docs):
    rows=[];changed=[]
    old=docs[C.R.TAG]['literal_interfaces']['rows']
    for prev,(stem,case,g,white,cs)in zip(old,T.literal_inputs(docs),strict=True):
        assert(prev['source'],prev['case'])==(stem,enc(case));corners=[]
        for prior,(marks,*_)in zip(prev['corners'],cs,strict=True):
            assert marks==prior['marks'];ch=changes(g,white,marks)
            wd=W.debit(prior['components'],marks,g)
            before=F(prior['parameter_values']['9/64'])+wd
            lower=before+ch['delta'];assert lower>=0
            corners.append(dict(marks=marks,before=before,**ch,lower=lower))
        total=sum(c['lower']for c in corners);assert total>=1
        row=dict(source=stem,case=enc(case),corners=corners,total=total)
        rows.append(row)
        if any(c['delta']<0 for c in corners):
            assert stem==LEAF and case in(6,7,12,13,18,19,24,25)
            retained=[]
            for c in corners:
                shared=next((b for b in c['components']if set(c['marks'])<=set(b['nodes'])),None)
                if shared and shared['type']=='D13':
                    im=T.mark_pairs(g,shared['nodes'],c['marks'],13,'D')
                    assert set(im)=={(1,11),(1,12)}and c['lower']==1 and c['delta']==0
                    retained.append(dict(marks=c['marks'],whole_type='D13',images=im,weight='1'))
            assert len(retained)==1
            row['independent_retained_unit']=retained[0];changed.append(row)
    assert len(rows)==216 and len(changed)==8
    assert all(r['total']>=1 for r in rows)
    later=[];debits=[];total=0
    for stem,keys in B.FINITE.items():
        cs=[]
        for section in keys:
            for path,c,gg,white,blocks in B.corners(docs[stem][section],section):
                g=A.matrix(gg);marks=c['marks'];ch=changes(g,white,marks)
                assert sorted(A.components(g,sorted(set(white)|set(marks))))==sorted(b['nodes']for b in blocks)
                if ch['delta']<0:
                    assert stem==H.TAG and ch['roles'][0]['role']=='D16_debit'
                    assert F(c['payment'])==F(1,2)
                    debits.append(dict(source=stem,path=path,marks=marks,delta=ch['delta']))
                cs.append(dict(path=path,marks=marks,delta=ch['delta'],roles=ch['roles']));total+=1
        later.append(dict(source=stem,corners=cs))
    assert len(later)==22 and total==564 and len(debits)==4
    variants=[]
    for p in docs[C.R.S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            g=A.matrix(v['gram']);cs=[]
            for c in v['corners']:
                ch=changes(g,v['white'],c['marks']);assert ch['delta']>=0
                cs.append(dict(marks=c['marks'],delta=ch['delta']))
            variants.append(dict(index=p['index'],extra_contact=v['extra_contact'],corners=cs))
    assert len(variants)==24
    for p in docs[C.R.TAG]['local_seals']:
        g=A.matrix(p['gram'])
        assert all(changes(g,p['white'],c['marks'])['delta']>=0 for c in p['corners'])
    d12=[]
    for tag,expected in((B.PART,280),(B.REST,52)):
        positive=[(path,c)for path,c in W.walk(docs[tag])if F(c['payment'])>0]
        assert len(positive)==expected
        assert all(b['type']not in('A16','D16')for _,c in positive for b in c['roles'])
        d12.append(dict(source=tag,positive_outputs=len(positive),
            marked_types=sorted({b['type']for _,c in positive for b in c['roles']})))
    return dict(literal_rows=rows,changed_literal_rows=changed,later_finite=later,
        later_corners=total,replaced_own_D16_outputs=debits,retained_D14_variants=variants,
        retained_rank13_seals=4,unchanged_chosen_D12_outputs=d12)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    local=selected(docs);inter=interfaces(docs)
    tables={}
    for typ,source,field,changed,new in(('A16',C.Q.TAG,'weighted_tables',ACHANGE,F(5,8)),
                                      ('D16',DW,'D16_weighted_pairs',DCHANGE,F(1,4))):
        data=docs[source][field][typ]if field=='weighted_tables'else docs[source][field]
        old={tuple(r['pair']):F(r['weight'])for r in data}
        assert sum(old.values())==F(C.actual_base((typ[0],16)),64)
        assert all(old[p]==(1 if typ=='A16'else F(1,2))for p in changed)
        weights={p:new if p in changed else w for p,w in old.items()}
        assert sum(weights.values())==F(actual_base((typ[0],16)),64)==F(219,2)
        tables[typ]=[dict(pair=p,weight=w)for p,w in sorted(weights.items())]
    assert C.actual_cross(('A',4),('A',12))==64*6
    assert actual_cross(('A',4),('A',12))==upper_cross(('A',4),('A',12))==64*F(17,2)
    for a in C.TYPES:
        assert 0<=actual_base(a)<=upper_base(C.representative(a))
        for b in C.TYPES:
            assert 0<=actual_cross(a,b)==actual_cross(b,a)<=upper_cross(C.representative(a),C.representative(b))
    assert actual_base(('D',13))==64*78
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in
        (('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_A16_D16_A4_A12_FRACTIONAL_FACE_TRANSFER',
        sources={p:digest(p)for p in sorted(paths)},baseline=W.TAG,scale=64,
        changed_pairs=dict(A16=sorted(ACHANGE),D16=sorted(DCHANGE)),weighted_tables=tables,
        local_proofs=local,current_interfaces=inter,
        new_same_component_costs=dict(A16='219/2',D16='219/2'),rank16_representative='219/2',
        coupons=[dict(small='whole A4 endpoint',small_nodes=[0,3],large='whole A12 second node',large_nodes=[1,10],
                    old='0',new='3/8',component_increment='3/2'),
                 dict(small='whole A4 endpoint',small_nodes=[0,3],large='whole A12 third node',large_nodes=[2,9],
                    old='1/2',new='3/4',component_increment='1')],
        actual_A4_A12_total='17/2',previous_A4_A12_total='6',type_pairs=len(C.TYPES)**2,
        ideal_weights_changed=False,context_rules_changed=False,triangle_minimum='2',quadrilateral_minimum='1',
        complete_forest_count_performed=False,new_dimension_bound=None,
        scope='Complete saved394/238 and394/578 selected domains, six whole seals and three proper exclusions; eight old affine-D11 literal uses repaired by a distinct whole D13 unit; complete current output inventory and arbitrary-support proof in report. New finite coupons must be counted; no forest maximum or dimension theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_same_component_costs'],len(out['current_interfaces']['changed_literal_rows']),flush=True)
