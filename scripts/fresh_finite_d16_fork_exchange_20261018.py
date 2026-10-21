#!/usr/bin/env python3
"""D16 cost111, D17 cost119, and a whole A4-end/A12-third half coupon.

Consumes saved complete selected maps and audits simultaneous literal payments.
No historical census or finite-forest optimization is repeated.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_finite_d17_a4_d13_half_transfer_20261018 as E
import fresh_rank16_existing_payment_exchange_20261018 as V

ROOT=E.ROOT;T=E.T;A=E.A;Q=E.Q;D=E.D;P=E.P;S=E.S;X=V.X
TAG='fresh_finite_d16_fork_exchange_20261018'
HALF={(0,14),(0,15),(2,12)};RESTORE={(1,15),(1,16)}
PROPER=[1011,1031,1614]
WHOLE=[1060,2489]
RADICAL=[1644,1664,1931,2324,2436,2515]
A15='fresh_finite_a15_endpoint_coupon_zero_20261018'
A16='fresh_finite_a16_endpoint_sparse_half_transfer_20261018'
BD15='fresh_finite_b15_d15_existing_payment_transfer_20261018'
DEPTH='fresh_depth8_finite8_zero_transfer_20261018'
STEMS=tuple(dict.fromkeys((*E.STEMS,E.TAG,V.TAG,A15,A16,BD15,DEPTH)))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    return x


def coupon(g,ids,marks):
    blocks=D.finite_blocks(g,ids);aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks];out=[]
    if aa is bb:return out
    for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
        if sm['type']!='A4' or la['type']!='A12':continue
        si=T.mark_pairs(g,sm['nodes'],[m],4,'A');li=T.mark_pairs(g,la['nodes'],[lm],12,'A')
        if si==[(0,),(3,)] and li==[(2,),(9,)]:
            out.append(dict(small_type='whole_A4',small_nodes=sm['nodes'],small_mark=m,small_images=si,
                large_type='whole_A12',large_nodes=la['nodes'],large_mark=lm,large_images=li,weight='1/2'))
    assert len(out)<=1
    return out


def proper(c,core,maps,target):
    g,J,du,dv,end=X.constraints(c,core);base=c['marked_component'];rows=[]
    for mp in maps:
        if tuple(mp['marked_images'])!=(2,12):continue
        image=dict(zip(mp['images'],base));missing=sorted(set(range(16))-set(image));entries=[]
        for new in missing:
            for old,node in image.items():
                if target[new][old]:
                    assert target[new][old]==3 and node in c['possible_contacts']
                    entries.append(dict(new_node=new,old_node=old,core_node=node,inverse_diagonal=J[node][node]))
        assert entries
        row=dict(map=mp['images'],missing_target_nodes=missing,forced_entries=entries)
        above=[e for e in entries if e['inverse_diagonal']>2]
        if above:
            assert above[0]['inverse_diagonal'] in (F(6),F(10,3))
            row.update(outcome='first_entry_projection_above_norm2',witness=above[0])
        else:
            assert c['index']in (1031,1614) and len(entries)==1
            e=entries[0];k=e['core_node'];assert e['new_node']==13 and e['inverse_diagonal']==2
            assert missing==[13,14,15] and target[13][14]==target[13][15]==3 and target[14][15]==0
            same=next(z for z in c['end_sets']if k in z)
            assert all(J[k][j]>0 for j in same)
            # Norm saturation forces the center's contact vector to e_k and
            # its residual to zero. Every negative projection coefficient
            # of a possible new leaf gives an independently paid short path.
            tests=[]
            for j in c['possible_contacts']:
                length=min(du[k]+dv[j],dv[k]+du[j])+3
                if J[k][j]<0:assert length<=8
                tests.append(dict(node=j,coefficient=J[k][j],negative_coefficient_shortcut=length if J[k][j]<0 else None))
            black=set(c['left'][1]+[c['right_copy_map'][i]for i in c['right'][1]])
            old_other=sorted((set(range(len(g)))-black|set(c['preferred_corner']))-set(base))
            row.update(outcome='zero_residual_fork_center',center_contact=k,projection_square='2',residual_square='0',
                sole_center_contact=True,old_other_corner_nodes=old_other,leaf_projection_tests=tests,
                conclusion='An old other-corner root is orthogonal to the center. Any new leaf with negative center product uses a negative coefficient and gives short8; otherwise the required simple bond is impossible.')
        rows.append(row)
    assert rows
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],inverse=core['inverse'],
        base_nodes=base,preferred_corner=c['preferred_corner'],end_projection_proof=end,maps=rows)


def seal(c,core):
    g=A.matrix(core['gram']);n=len(g);assert core['gauge']==[1]*n and all(g[i][i]==2 for i in range(n))
    u=c['left'][1];v=[c['right_copy_map'][i]for i in c['right'][1]];black=set(u+v);white=sorted(set(range(n))-black)
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0]for b in black)
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
    selected=T.mark_pairs(g,c['marked_component'],c['preferred_corner'],16,'D')
    assert len(c['marked_component'])==16 and selected and set(selected)<=HALF
    if c['method']=='radical':
        assert c['index']in RADICAL and n==19 and Q.inertia(g)==[17,1,1]
        z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j]for j in range(n))==0 for i in range(n))
        allowed=c['possible_contacts'];assert all(z[i]for i in allowed)
        du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);short=[]
        for i,j in combinations(allowed,2):
            if z[i]*z[j]<0:
                length=min(du[i]+dv[j],dv[i]+du[j])+2;assert length<=8
                short.append(dict(nodes=[i,j],length=length))
        assert short;sealing=dict(method='radical',radical=z,allowed_contacts=allowed,mixed_sign_shortcuts=short)
    else:
        assert c['index']in WHOLE and n==18 and Q.inertia(g)==[17,1,0]
        gg,J,du,dv,end=X.constraints(c,core);assert gg==g
        assert sorted(set(white)|set(c['preferred_corner']))==c['marked_component']
        remain=sorted(set(c['possible_contacts'])-set(c['marked_component']))
        assert remain==c['omitted_black'] and sorted(J[i][i]for i in remain)==[4,8]
        assert all(sum(i in ee for i in remain)==1 for ee in c['end_sets'])
        sealing=dict(method='whole_selected_D16_projection',whole_nodes=c['marked_component'],end_projection=end,
            remaining_contacts=remain,remaining_diagonals=[[i,J[i][i]]for i in remain])
    corners=[]
    for cc in c['all_corners']:
        marks=cc['preferred_corner'];ids=sorted(set(white)|set(marks));assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        blocks=D.finite_blocks(g,ids);shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
        val=F();reason=None;images=[];coupons=[]
        if marks==c['preferred_corner']:val=F(1,2);reason='selected_D16_half';images=selected
        elif c['index']in (1644,1664,1931,2324):
            if shared and shared['type']=='D17':
                images=T.mark_pairs(g,shared['nodes'],marks,17,'D');assert images and set(images)<=RESTORE
                val=F(1,2);reason='restored_whole_D17_half'
        else:
            coupons=coupon(g,ids,marks)
            if coupons:val=F(1,2);reason='whole_A4_end_A12_third_half'
        corners.append(dict(marks=marks,components=blocks,payment=val,reason=reason,canonical_images=images,new_coupons=coupons))
    assert sum(r['payment']for r in corners)==1 and sum(bool(r['payment'])for r in corners)==2
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],white=white,selected_images=selected,
        sealing=sealing,corners=corners,total='1',whole_core_sealed_unless_short8=True)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[DEPTH]['feature_depth']==8 and docs[E.TAG]['new_cost']=='118'
    assert docs[A15]['new_cost']=='103' and docs[A16]['new_cost']=='111' and docs[BD15]['new_costs']==dict(B15='103',D15='103')
    rec=next(r for r in docs[S.TAG]['rules']if r['kind']=='D'and r['rank']==16)
    old={tuple(r['pair']):F(r['weight'])for r in rec['weighted_pairs']};assert sum(old.values())==114
    old.update({p:F(1,2)for p in V.DHALF});assert sum(old.values())==F(225,2)
    assert all(old[p]==1 for p in HALF);weights={p:(F(1,2)if p in HALF else w)for p,w in old.items()}
    assert sum(weights.values())==111
    old17={tuple(r['pair']):F(r['weight'])for r in docs[E.TAG]['weighted_pairs']}
    assert sum(old17.values())==118 and all(old17[p]==0 for p in RESTORE)
    new17={p:(F(1,2)if p in RESTORE else w)for p,w in old17.items()};assert sum(new17.values())==119
    target=rec['target'];adj,dist,_=P.distance_data(target)
    assert all(dist[i][j]>=10 for i,j in HALF)
    assert all({tuple(sorted(m[i]for i in p))for p in HALF}==HALF for m in rec['automorphisms'])
    assert all({tuple(sorted(m[i]for i in p))for p in RESTORE}==RESTORE for m in docs[E.TAG]['automorphisms'])
    assert all(weights[tuple(p)]==1 for p in rec['retained_basic_pairs'])
    assert weights[1,14]==weights[1,15]==F(1,2)
    half=docs[T.TAG];sup=V.support('D',half);over=docs[A.OVER]
    assert V.validate_maps(sup,target,over)==578
    affected=[r for r in sup['rows']if set(map(tuple,r['marked_images']))&HALF]
    assert sorted(r['index']for r in affected)==sorted(PROPER+WHOLE+RADICAL)
    proofs=[];prop=[]
    for r in affected:
        c=over['isolation']['cases'][r['index']];core=over['isolation']['cores'][c['core_key']]
        if c['index']in PROPER:prop.append(proper(c,core,r['maps'],target))
        else:proofs.append(seal(c,core))
    previous={E.key(r):r for r in docs[E.TAG]['literal_rows']};assert len(previous)==216
    for stem in (A15,A16,BD15):
        assert len(docs[stem]['literal_rows'])==216
        for r in docs[stem]['literal_rows']:
            oldrow=previous[E.key(r)];assert len(r['corners'])==len(oldrow['corners'])
            for x,y in zip(r['corners'],oldrow['corners']):
                assert x['marks']==y['marks'] and F(x['debit'])==0
                assert F(x.get('retained_payment',x.get('retained_payment_lower_bound')))==F(y['retained_payment'])
    literal=[];hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=previous[stem,E.N.normalize(index)];assert len(cs)==len(prev['corners']);corners=[]
        for j,(marks,*_)in enumerate(cs):
            assert marks==prev['corners'][j]['marks'];ids=sorted(set(white)|set(marks));delta=F();images=[];typ=None;coupons=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                blocks=D.finite_blocks(g,ids);shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
                if shared and shared['type']in('D16','D17'):
                    typ=shared['type'];images=T.mark_pairs(g,shared['nodes'],marks,len(shared['nodes']),'D')
                    changes=HALF if typ=='D16'else RESTORE
                    if set(images)&changes:
                        assert set(images)<=changes;delta=F(-1,2)if typ=='D16'else F(1,2);hits[typ]+=1
                coupons=coupon(g,ids,marks)
                if coupons:delta+=F(1,2);hits['new_coupon']+=1
            val=F(prev['corners'][j]['retained_payment'])+delta;assert val>=0
            corners.append(dict(marks=marks,type=typ,images=images,old_payment=prev['corners'][j]['retained_payment'],delta=delta,new_coupons=coupons,retained_payment=val))
        total=sum(r['retained_payment']for r in corners);assert total>=1,(stem,index,total)
        literal.append(dict(source=stem,case=index,corners=corners,total=total))
    assert len(literal)==216
    # New local alternatives of the intervening transfers never consume a
    # changed D16 unit; restorations only increase every old D17 payment.
    for stem in (V.TAG,E.R.TAG,E.TAG,A15,A16,BD15):
        for proof in docs[stem].get('selected_proofs',[]):
            if 'corners'not in proof:
                assert stem==BD15 and 'maps'in proof  # excluded proper completions spend no alternate payment
                continue
            g=A.matrix(proof['gram'])
            for c in proof['corners']:
                if not F(c.get('payment','0')):continue
                for b in c['components']:
                    if b['type']=='D16'and set(c['marks'])<=set(b['nodes']):
                        images=T.mark_pairs(g,b['nodes'],c['marks'],16,'D');assert not set(images)&HALF,(stem,proof['index'])
    paths=[Path(__file__),Path(E.__file__),Path(V.__file__),Path(X.__file__),Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__),Path(S.__file__)]
    paths +=[ROOT/f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return enc(dict(status='PASS_FINITE_D16_FORK_EXCHANGE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        D16_half_pairs=sorted(HALF),D17_restored_half_pairs=sorted(RESTORE),old_costs=dict(D16='225/2',D17='118'),new_costs=dict(D16='111',D17='119'),
        D16_weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(weights.items())],D17_weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(new17.items())],
        complete_selected_support=sup,affected_contexts=[r['index']for r in affected],proper_completion_proofs=prop,selected_proofs=proofs,
        proper_outcomes=dict(Counter(m['outcome']for r in prop for m in r['maps'])),new_cross_component_coupons=[dict(name='A4end_A12third',small='whole A4',small_nodes=[0,3],large='whole A12',large_nodes=[2,9],weight='1/2',component_pair_cost='2')],
        literal_configurations=216,literal_changed_corners=dict(hits),literal_rows=literal,
        retained_rank15_alternatives=dict(D16_pairs=[[1,14],[1,15]],weight='1/2'),compatible_feature_depth=8,
        triangle_minimum=2,quadrilateral_minimum=1,finite_forest_recomputed=False,new_dimension_bound=None,
        scope='Saved complete394-pattern/578-map D16 domain, three proper completion exclusions, two whole-selected Lorentz seals and six radical seals; all216 simultaneous latest literal interfaces. Existing D17 zeros are restored to halves and a new whole A4endpoint/A12third half costs2 per component pair. Written current-priority composition is required; no historical census, finite-forest optimization or ambient bound.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_costs'],out['proper_outcomes'],'literal',out['literal_changed_corners'],flush=True)
