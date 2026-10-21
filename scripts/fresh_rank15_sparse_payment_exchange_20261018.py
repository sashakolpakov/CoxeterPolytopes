#!/usr/bin/env python3
"""Rank15 local exchange with its A17 restoration and typed coupon counted.

Consumes complete saved selected domains and finite Bellman certificates.
No historical census, native-map enumeration or new optimizer is executed.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_b15_d15_existing_payment_transfer_20261018 as U
import fresh_finite_d16_fork_exchange_20261018 as J
import fresh_ambient60_payment_frontier_20261018 as Z

R=U.R; T=R.T; A=R.A; Q=R.Q; D=R.D; P=R.P; ROOT=R.ROOT
TAG='fresh_rank15_sparse_payment_exchange_20261018'
A15='fresh_finite_a15_endpoint_coupon_zero_20261018'
A16='fresh_finite_a16_endpoint_sparse_half_transfer_20261018'
COMPLETE='fresh_order11_order16_dimension60_20261018'
DEPTH='fresh_depth10_finite7_zero_transfer_20261018'
AHALF={(0,12),(2,14)}
AZERO={(0,14)}
BZ={(0,12),(1,13),(1,14),(2,14),(3,14),(4,14)}
DZ={(1,13),(1,14)}
CHANGES=dict(A=AHALF|AZERO,B=BZ,D=DZ)
AWHOLE={1029,1612}; ARADICAL={1638,1671,2416}
STEMS=tuple(dict.fromkeys((*J.STEMS,J.TAG,R.TAG,R.R.TAG,U.TAG,A15,A16,
                         COMPLETE,DEPTH,Z.TAG,R.T.T.TAG)))


def enc(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [enc(v) for v in x]
    return x


def selected_support(rec,over):
    """Validate every saved map, including rows with no induced embedding."""
    assert rec['complete_eligible_patterns']==536
    target=rec['target']; maps=0; affected=[]
    for row in rec['selected_patterns']:
        c=over['isolation']['cases'][row['index']]
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        base=c['marked_component']; graph=T.T.labels(g,base); seen=set(); pairs=set(); hit=[]
        for m in row['maps']:
            mp=m['images']; assert len(mp)==len(base)==len(set(mp))
            assert all(0<=i<15 for i in mp) and tuple(mp) not in seen
            seen.add(tuple(mp)); maps+=1
            assert all(graph[i][j]==target[mp[i]][mp[j]] for i in range(len(base)) for j in range(len(base)))
            pair=tuple(sorted(mp[base.index(i)] for i in c['preferred_corner']))
            assert list(pair)==m['marked_images']; pairs.add(pair)
            if pair in CHANGES[rec['kind']]: hit.append(mp)
        assert sorted(pairs)==list(map(tuple,row['marked_images']))
        if hit:
            assert pairs<=CHANGES[rec['kind']]
            affected.append(dict(index=c['index'],maps=hit))
    assert maps==dict(A=482,B=240,D=870)[rec['kind']]
    expected=dict(A=sorted(AWHOLE|ARADICAL),B=[],D=sorted(U.PROPER+U.WHOLE))
    assert [r['index'] for r in affected]==expected[rec['kind']]
    return dict(kind=rec['kind'],patterns=536,maps=maps,affected=affected)


def a_proof(c,core,a16):
    g=A.matrix(core['gram']); n=len(g)
    assert core['gauge']==[1]*n and all(g[i][i]==2 for i in range(n))
    u=c['left'][1]; v=[c['right_copy_map'][i] for i in c['right'][1]]
    white=sorted(set(range(n))-set(u+v))
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0] for b in u+v)
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for key in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[key]==c[key]
    assert len(c['marked_component'])==15
    images=T.mark_pairs(g,c['marked_component'],c['preferred_corner'],15,'A')
    assert images and set(images)<=CHANGES['A']
    if c['index'] in ARADICAL:
        assert c['method']=='radical' and n==19 and Q.inertia(g)==[17,1,1]
        z=list(map(F,core['radical'][0]))
        assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
        allowed=c['possible_contacts']; assert all(z[i] for i in allowed)
        du=Q.graph_dist(g,u); dv=Q.graph_dist(g,v); shorts=[]
        for i,j in combinations(allowed,2):
            if z[i]*z[j]<0:
                length=min(du[i]+dv[j],du[j]+dv[i])+2
                assert length<=8; shorts.append(dict(nodes=[i,j],length=length))
        assert shorts
        seal=dict(method='actual_radical',radical=z,allowed_contacts=allowed,mixed_sign_shortcuts=shorts)
    else:
        assert c['index'] in AWHOLE and n==18 and Q.inertia(g)==[17,1,0]
        gg,inv,du,dv,end=R.X.constraints(c,core); assert gg==g
        remain=sorted(set(c['possible_contacts'])-set(c['marked_component']))
        assert sorted(inv[i][i] for i in remain)==[4,8,14]
        seal=dict(method='whole_A15_positive_end_projection',whole_selected_nodes=c['marked_component'],
                  end_projection=end,remaining_contacts=remain,remaining_diagonals=[[i,inv[i][i]] for i in remain])
    corners=[]
    for old in c['all_corners']:
        marks=old['preferred_corner']; ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        blocks=D.finite_blocks(g,ids)
        shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
        payment=F(); reason=None; pairs=[]
        if marks==c['preferred_corner']:
            pairs=images; payment=F() if c['index']==1638 else F(1,2)
            reason='new_A15_zero' if not payment else 'new_A15_half'
        elif shared and shared['type']=='A16' and c['index']!=1638:
            pairs=T.mark_pairs(g,shared['nodes'],marks,16,'A')
            assert pairs==[(1,13),(2,14)] and all(a16[p]==F(1,2) for p in pairs)
            payment=F(1,2); reason='retained_whole_A16_half'
        elif shared and shared['type']=='A17' and c['index']==1638:
            pairs=T.mark_pairs(g,shared['nodes'],marks,17,'A')
            assert pairs==[(1,15)]; payment=F(1); reason='restored_whole_A17_unit'
        corners.append(dict(marks=marks,components=blocks,payment=payment,reason=reason,canonical_images=pairs))
    assert sum(x['payment'] for x in corners)==1
    assert sum(bool(x['payment']) for x in corners)==(1 if c['index']==1638 else 2)
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],white=white,
                selected_images=images,sealing=seal,corners=corners,total='1')


def finite_count(doc,newcost):
    """Same residual Bellman tables; only fixed D15/A17 base costs change."""
    def base(t): return int(2*F(newcost[t[0]+str(t[1])])) if t[0]+str(t[1]) in newcost else Z.B.base(t)
    def cross(a,b): return Z.current_cross(a,b)+(4 if {a,b}=={('A',2),('D',13)} else 0)
    types=[t for t in Z.B.B.types() if t[1]<=60]
    for a in types:
        assert base(Z.B.B.representative(a))>=base(a)
        for b in types:
            assert cross(a,b)==cross(b,a)
            assert cross(a,b)<=Z.B.upper_cross(Z.B.B.representative(a),Z.B.B.representative(b))
    dom=Z.large_domain(); assert len(dom)==1044
    old={tuple(r['counts']):r for r in doc['finite']['alternatives']}
    tables={tuple(r['key']):r for r in doc['finite']['tables']}
    assert set(old)=={cs for cs,used,key in dom} and len(tables)==792
    assert {key for cs,used,key in dom}==set(tables)
    checks=sum(Z.check_table(tables[key],key) for key in sorted(tables))
    residual=[Z.B.B.representative(('A' if r<=8 else 'D',r)) for r in Z.ITEMS]
    assert all(base(t)==Z.B.base(t) for t in residual)
    assert all(Z.B.upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[]
    for cs,used,key in dom:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert all(Z.B.upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(Z.ITEMS,residual):
            assert Z.B.item(r,key)==base(t)+sum(Z.B.upper_cross(t,b) for b in fixed)
        value=sum(base(t) for t in fixed)+tables[key]['values'][60-used]
        assert value==old[cs]['doubled_upper']+2*(cs[8]-cs[6])
        rows.append(dict(counts=list(cs),key=list(key),doubled_upper=value))
    assert max(r['doubled_upper'] for r in rows)==826
    witness=[('D',17),('D',17),('D',26)]
    assert sum(base(t) for t in witness)+sum(cross(a,b) for i,a in enumerate(witness) for b in witness[i+1:])==826
    assert cross(('A',2),('D',13))==Z.B.upper_cross(('A',2),('D',13))==10
    return dict(exact_maximum='413',curvature='405',complete_fixed_vectors=1044,
                reused_residual_tables=792,Bellman_inequalities_checked=checks,
                simultaneous_type_pairs=len(types)**2,alternatives=rows,
                above_curvature=sum(r['doubled_upper']>810 for r in rows),
                unchanged_attainer=['D17','D17','D26'],new_optimization=False,
                scope='Exact upper and formal-forest attainment; depth10 deletions only lower costs and do not touch the attainer. No geometric realization claimed.')


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items(): assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[COMPLETE]['dimension_bound']==60 and docs[DEPTH]['feature_depth']==10
    assert docs[J.TAG]['new_costs']==dict(D16='111',D17='119')
    over=docs[A.OVER]; scouts={r['kind']:r for r in docs[R.R.TAG]['records']}
    oldrules={r['kind']:r for r in docs[U.TAG]['rules']}; oldrules['A']=docs[A15]
    rules=[]; supports=[]; proofs=[]
    a16={tuple(x['pair']):F(x['weight']) for x in docs[A16]['weighted_pairs']}
    for kind in 'ABD':
        old={tuple(r['pair']):F(r['weight']) for r in oldrules[kind]['weighted_pairs']}
        assert sum(old.values())==103
        new=dict(old)
        for p in CHANGES[kind]:
            assert old[p]==(F(1,2) if kind=='D' else F(1))
            new[p]=F(1,2) if kind=='A' and p in AHALF else F()
        assert sum(new.values())==dict(A=101,B=97,D=102)[kind]
        rec=scouts[kind]; _,dist,_=P.distance_data(rec['target'])
        assert all(dist[i][j]>=10 for i,j in CHANGES[kind])
        assert all(new[tuple(p)]==1 for p in rec['retained_basic_pairs'])
        assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGES[kind]}==CHANGES[kind] for mp in rec['automorphisms'])
        sup=selected_support(rec,over); supports.append(sup)
        for hit in sup['affected']:
            c=over['isolation']['cases'][hit['index']]; core=over['isolation']['cores'][c['core_key']]
            if kind=='A': proof=a_proof(c,core,a16)
            else:
                proof=U.proof(c,core,hit['maps'],rec['target'])
                previous=next(p for p in docs[U.TAG]['selected_proofs'] if p['index']==c['index'])
                assert enc(proof)==previous
                if 'corners' in proof:
                    for corner in proof['corners']:
                        if corner['reason']=='selected_D15_half': corner['payment']=F();corner['reason']='new_D15_zero'
                        elif corner['reason']=='existing_whole_A2_D13_short_fork_leaf_half':
                            corner['payment']=F(1);corner['reason']='upgraded_whole_A2_D13_short_fork_leaf_unit'
                    for coupon in proof['coupons']:coupon['weight']='1'
                    assert sum(x['payment'] for x in proof['corners'])==1
            proofs.append(proof)
        rules.append(dict(kind=kind,target=rec['target'],automorphisms=rec['automorphisms'],
                          changed_pairs=sorted(CHANGES[kind]),weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(new.items())],
                          new_cost=str(sum(new.values())),retained_basic_pairs=rec['retained_basic_pairs']))
    olda17=next(r for r in docs[R.T.T.TAG]['records'] if r['kind']=='A')
    assert [1,15] not in olda17['retained_pairs'] and (1,15) not in R.E.ARESTORE
    assert docs[R.E.TAG]['new_costs']['A17']=='115'
    previous={(r['source'],U.E.N.normalize(r['case'])):r for r in docs[J.TAG]['literal_rows']}
    literal=[]; hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=previous[stem,U.E.N.normalize(index)]; corners=[]
        assert len(cs)==len(prev['corners'])
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks']; ids=sorted(set(white)|set(marks)); debit=F();images=[];typ=None
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                b=next((b for b in D.finite_blocks(g,ids) if set(marks)<=set(b['nodes'])),None)
                if b and b['type'] in ('A15','B15','D15'):
                    typ=b['type']; images=T.mark_pairs(g,b['nodes'],marks,15,typ[0])
                    if set(images)&CHANGES[typ[0]]:
                        assert set(images)<=CHANGES[typ[0]]
                        debit=F(1,2) if typ=='D15' or typ=='A15' and set(images)<=AHALF else F(1)
                        hits[typ]+=1
            value=F(prev['corners'][j]['retained_payment'])-debit;assert value>=0
            corners.append(dict(marks=marks,type=typ,images=images,debit=debit,retained_payment_lower_bound=value))
        total=sum(r['retained_payment_lower_bound'] for r in corners);assert total>=1
        literal.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(literal)==216 and not hits
    # Earlier positive replacements are retained; changed U selected halves
    # are explicitly replaced above, with their proper exclusions unchanged.
    interface_records=[]
    for stem in (R.TAG,R.E.TAG,U.E.TAG,J.TAG,A16,A15):
        for proof in docs[stem].get('selected_proofs',[]):
            if 'corners' not in proof: continue
            for corner in proof['corners']:
                if not F(corner.get('payment','0')):continue
                for b in corner['components']:
                    if b['type'] in ('A15','B15','D15') and set(corner['marks'])<=set(b['nodes']):
                        im=T.mark_pairs(A.matrix(proof['gram']),b['nodes'],corner['marks'],15,b['type'][0])
                        assert not set(im)&CHANGES[b['type'][0]]
                        interface_records.append(dict(source=stem,case=proof['index'],type=b['type'],images=im))
    routing=docs[R.TAG]['complete_high_overlap_routing']
    assert routing==dict(total=2614,bound_at_most14=1691,protected_or_adaptive=83,bound_at_least15=840,displayed_larger_than15=304,eligible_rank15=536)
    newcost=dict(A15='101',B15='97',D15='102',A17='116')
    count=finite_count(docs[Z.TAG],newcost)
    paths={Path(__file__),Path(U.__file__),Path(J.__file__),Path(R.__file__),Path(R.X.__file__),Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__),Path(Z.__file__),Path(Z.B.__file__),Path(Z.B.B.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_RANK15_SPARSE_PAYMENT_EXCHANGE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        new_costs=newcost,rules=rules,complete_selected_support=supports,selected_proofs=proofs,
        A17_restored_unit_pair=[1,15],A17_cost_increase='1',
        upgraded_coupon=dict(small='whole A2',small_marks=[0,1],large='whole D13',large_marks=[11,12],old_weight='1/2',new_weight='1',component_pair_cost_increase='2',new_total_component_pair_cost='5'),
        literal_configurations=216,literal_debited_corners=dict(hits),literal_rows=literal,
        retained_positive_interfaces=interface_records,complete_high_overlap_routing=routing,
        finite_rank60_count=count,triangle_minimum=2,quadrilateral_minimum=1,ideal_weights_unchanged=True,
        compatible_deep_finite7_feature_depth=10,new_dimension_bound=None,new_finite_volume_example=False,
        scope='Complete saved rank15 map validation, five new A15 arbitrary-support seals, reconstructed eight D15 proper/whole proofs, simultaneous216 current literal interfaces and positive finite alternatives. All A17/coupon increases counted in the complete1044-vector finite comparison using checked saved792 Bellman tables. Universal ideal/short-path and nonselected reductions remain written inherited theorem inputs. No new historical census or global dimension claim.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:
        assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['new_costs'],'finite',result['finite_rank60_count']['exact_maximum'],flush=True)
