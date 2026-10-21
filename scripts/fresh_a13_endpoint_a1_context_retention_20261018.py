#!/usr/bin/env python3
"""Retain a whole A13 endpoint unit only with a whole unmarked A1.

Complete native selected support, four radical seals and one saturated
two-witness seal; current chosen-output retention and monotone accounting.
No new dynamic program, obstruction census or dimension bound.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import fresh_d16_a1_long_leaf_d17_zero_ambient58_frontier_20261018 as H
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

ROOT=H.ROOT;R=H.R;T=R.T;A=R.A;Q=R.Q;D=R.D;B=R.B;W=R.W
TAG='fresh_a13_endpoint_a1_context_retention_20261018'
PAIR=(0,12);CASES=[623,1349,1554,1627,1640]
STEMS=tuple(dict.fromkeys((H.TAG,*H.STEMS,*W.STEMS,X.TAG)))
enc=R.enc

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return H.actual_base(t)-(64 if t==('A',13)else 0)
actual_cross=H.actual_cross
def context(fs):return H.context(fs)+(64*fs.count(('A',13))if('A',1)in fs else 0)
def value(fs):
    z=sum(actual_base(t)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)
    assert z==H.value(fs)-(0 if('A',1)in fs else 64)*fs.count(('A',13))
    return z

def endpoint(g,blocks,marks):
    shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
    if shared is None or shared['type']!='A13':return False
    images=T.mark_pairs(g,shared['nodes'],marks,13,'A')
    if PAIR in images:assert images==[PAIR];return True
    return False

def basic(c,core):
    g=A.matrix(core['gram']);n=len(g)
    assert core['gauge']==[1]*n and all(g[i][i]==2 for i in range(n))
    u=c['left'][1];v=[c['right_copy_map'][i]for i in c['right'][1]]
    white=sorted(set(range(n))-set(u+v))
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0]for b in u+v)
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for k in('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
    assert len(c['marked_component'])==13
    blocks=D.finite_blocks(g,sorted(set(white)|set(c['preferred_corner'])))
    assert endpoint(g,blocks,c['preferred_corner'])
    witnesses=[b for b in blocks if b['type']=='A1'and not set(c['preferred_corner'])&set(b['nodes'])]
    assert witnesses
    return g,white,u,v,blocks,witnesses

def radical(c,core):
    g,white,u,v,blocks,witnesses=basic(c,core);n=len(g)
    assert c['index']in CASES[:-1]and c['method']=='radical'and Q.inertia(g)==[n-2,1,1]
    assert len(core['radical'])==1
    z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j]for j in range(n))==0 for i in range(n))
    allowed=c['possible_contacts'];assert all(z[i]for i in allowed)
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);short=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],dv[i]+du[j])+2;assert length<=5
            short.append(dict(nodes=[i,j],length=length))
    assert short and max(r['length']for r in short)==5
    return dict(index=c['index'],core_key=c['core_key'],gram=g,white=white,signature=Q.inertia(g),
        selected_marks=c['preferred_corner'],selected_components=blocks,whole_A1_witnesses=witnesses,
        radical=z,allowed_contacts=allowed,mixed_sign_shortcuts=short,
        all_further_whites_core_orthogonal_unless_short5=True,selected_payment='1',total='1')

def saturated(c,core):
    g,white,u,v,blocks,witnesses=basic(c,core)
    gg,J,du,dv,end=X.constraints(c,core)
    assert gg==g and c['index']==1640 and len(g)==17 and Q.inertia(g)==[16,1,0]
    remaining=sorted(set(c['possible_contacts'])-set(c['marked_component']))
    assert remaining==[2,3,12,13]and[J[i][i]for i in remaining]==[6,2,6,2]
    assert all(J[i][j]>0 for e in c['end_sets']for i in e for j in e)
    assert end['both_ends_shortcut_bound']<=4
    isolated=sorted(b['nodes'][0]for b in witnesses);assert isolated==[3,13]
    assert J[3][13]==-2
    # A nonzero additional contact has exactly one simple bond at3 or13.
    # Its projection has norm2, so positive-definite residual space forces
    # residual zero. These are the two forced full-core-span roots.
    roots=[]
    for node in isolated:
        coefficients=[-J[i][node]for i in range(len(g))]
        assert all(sum(g[i][j]*coefficients[j]for j in range(len(g)))==-int(i==node)for i in range(len(g)))
        assert sum(coefficients[i]*g[i][j]*coefficients[j]for i in range(len(g))for j in range(len(g)))==2
        roots.append(dict(contact=node,coefficients=coefficients,contact_square='1',projection_square='2',residual_square='0'))
    # Repeated same-position roots have product2; opposite candidates have
    # product-2 and cannot both belong to a positive-definite white Gram.
    products=[[J[i][j]for j in isolated]for i in isolated]
    assert products==[[2,-2],[-2,2]]and Q.inertia(products)==[1,0,1]
    variants=[]
    for node in[None,*isolated]:
        full=[r[:]for r in g];ww=white[:]
        if node is not None:
            for i,row in enumerate(full):row.append(F(-int(i==node)))
            full.append([F(-int(i==node))for i in range(len(g))]+[F(2)]);ww.append(len(g))
        ids=sorted(set(ww)|set(c['preferred_corner']))
        assert Q.inertia(A.sub(full,ww))==[len(ww),0,0]
        assert Q.inertia(A.sub(full,ids))==[len(ids),0,0]
        current=D.finite_blocks(full,ids);assert endpoint(full,current,c['preferred_corner'])
        a1=[b for b in current if b['type']=='A1'and not set(c['preferred_corner'])&set(b['nodes'])]
        assert len(a1)==(2 if node is None else 1)
        assert all(set(b['nodes'])<={3,13}for b in a1)
        variants.append(dict(extra_contact=node,gram=full,white=ww,components=current,whole_A1_witnesses=a1,payment='1'))
    return dict(index=1640,core_key=c['core_key'],gram=g,white=white,signature=[16,1,0],
        selected_marks=c['preferred_corner'],end_projection=end,remaining_contacts=remaining,
        remaining_inverse_diagonals=[[i,J[i][i]]for i in remaining],forced_roots=roots,
        forced_root_Gram=products,maximum_further_core_contact_whites=1,
        core_orthogonal_whites_orthogonal_to_forced_roots=True,variants=variants,
        at_least_one_original_A1_remains_whole=True,total='1')

def local(docs):
    over=docs[A.OVER];target=T.target(13,'A')
    patterns=[c for c in over['isolation']['cases']if(c.get('component_order_bound')or 0)>=13 and len(c['marked_component'])<=13]
    assert len(patterns)==801
    rows=[];hits=[];maps_count=0
    for c in patterns:
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        maps=T.T.mapped(g,c['marked_component'],c['preferred_corner'],target);maps_count+=len(maps)
        affected=[m for m in maps if tuple(m['marked_images'])==PAIR]
        rows.append(dict(index=c['index'],maps=maps,affected_maps=affected))
        if affected:
            assert len(c['marked_component'])==13 and len(affected)==len(maps)==2
            hits.append(c)
    assert maps_count==522 and[c['index']for c in hits]==CASES
    radical_proofs=[radical(c,over['isolation']['cores'][c['core_key']])for c in hits[:-1]]
    last=hits[-1];sat=saturated(last,over['isolation']['cores'][last['core_key']])
    return dict(complete_patterns=801,complete_native_maps=522,target=target,selected_patterns=rows,
        affected_contexts=CASES,affected_maps=10,all_affected_components_already_whole_A13=True,
        radical_proofs=radical_proofs,saturated_two_witness_proof=sat,
        arbitrary_additional_whites_and_contact_magnitudes_covered=True,guarded_weight='1',unguarded_weight='0')

def interfaces(docs):
    literal=[]
    for old,(stem,case,g,white,cs)in zip(docs[H.TAG]['current_interfaces']['literals'],T.literal_inputs(docs),strict=True):
        assert(old['source'],old['case'])==(stem,enc(case));rows=[]
        for prev,(marks,*_)in zip(old['corners'],cs,strict=True):
            assert prev['marks']==marks
            ids=sorted(set(white)|set(marks));finite=Q.inertia(A.sub(g,ids))==[len(ids),0,0]
            blocks=D.finite_blocks(g,ids)if finite else[]
            assert not endpoint(g,blocks,marks),(stem,case,marks)
            rows.append(dict(marks=marks,finite=finite,components=blocks,retained_lower=prev['lower']))
        assert sum(F(r['retained_lower'])for r in rows)==F(old['total'])>=1
        literal.append(dict(source=stem,case=case,corners=rows,total=old['total']))
    assert len(literal)==216
    later=[]
    for stem,keys in B.FINITE.items():
        rows=[]
        for key in keys:
            for path,p,gg,w,blocks in B.corners(docs[stem][key],key):
                g=A.matrix(gg);actual=D.finite_blocks(g,sorted(set(w)|set(p['marks'])))
                assert actual==blocks and not endpoint(g,blocks,p['marks']),(stem,path)
                rows.append(dict(path=path,marks=p['marks'],A13_endpoint=False))
        later.append(dict(source=stem,corners=rows))
    assert len(later)==22 and sum(len(r['corners'])for r in later)==564
    newer=[]
    for stem,key in((R.L.TAG,'local_proofs'),(R.R.TAG,'local_seals'),(R.TAG,'local_proof'),
                    (R.N.TAG,'local_proof'),(H.Z.U.TAG,'local_proof')):
        rows=[]
        for path,p,gg,w,blocks in B.corners(docs[stem][key],key):
            g=A.matrix(gg);ids=sorted(set(w)|set(p['marks']))
            finite=Q.inertia(A.sub(g,ids))==[len(ids),0,0]
            assert(D.finite_blocks(g,ids)if finite else[])==blocks
            assert not endpoint(g,blocks,p['marks']),(stem,path)
            rows.append(path)
        assert rows;newer.append(dict(source=stem,unchanged_corners=rows))
    exchanges=[]
    for p in docs[H.TAG]['local_proof']['whole_proofs']:
        g=A.matrix(p['old_complete_seal']['gram']);white=p['old_complete_seal']['white']
        for c in p['corners']:
            blocks=D.finite_blocks(g,sorted(set(white)|set(c['marks'])))
            assert blocks==c['components']and not endpoint(g,blocks,c['marks'])
        exchanges.append(p['index'])
    assert exchanges==H.CASES
    variants=[]
    for p in docs[H.Z.U.S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            g=A.matrix(v['gram'])
            for c in v['corners']:
                blocks=D.finite_blocks(g,sorted(set(v['white'])|set(c['marks'])))
                assert blocks==c['components']and not endpoint(g,blocks,c['marks'])
            variants.append([p['index'],v['extra_contact']])
    assert len(variants)==24
    d12=[];near=0
    for stem in(B.PART,B.REST):
        outputs=[];n=0
        for path,p in W.walk(docs[stem]):
            if F(p['payment'])<=0:continue
            n+=1
            if any(b['type']=='A13'for b in p['roles']):
                assert all(tuple(pair) in((0,11),(1,12))for pair in p['images'])
                assert PAIR not in list(map(tuple,p['images']));near+=1
                outputs.append(dict(path=path,images=p['images'],payment=p['payment']))
        d12.append(dict(source=stem,positive_outputs=n,retained_A13_near_endpoints=outputs))
    assert sum(r['positive_outputs']for r in d12)==332 and near==10
    return dict(historical_literals=literal,later_finite=later,newer_seals=newer,
        unchanged_D16_exchanges=exchanges,unchanged_D14_variants=variants,chosen_D12_outputs=d12,
        chosen_A13_near_endpoint_outputs=near,changed_current_chosen_outputs=0)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    proof=local(docs);inter=interfaces(docs)
    assert docs[H.TAG]['finite']['exact_actual_maximum']=='378'
    for t in H.C.TYPES:
        assert actual_base(t)==H.actual_base(t)-(64 if t==('A',13)else 0)and actual_base(t)>=0
        for u in H.C.TYPES:assert actual_cross(t,u)==H.actual_cross(t,u)
    assert F(actual_base(('A',13)),64)==77
    witness=[('D',17),('D',17),('D',24)];plateau=[('A',13)]+[('A',15)]*3
    assert value(witness)==H.value(witness)==378*64
    assert H.value(plateau)==378*64 and value(plateau)==377*64
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return enc(dict(status='PASS_A13_ENDPOINT_A1_CONTEXT_RETENTION_AND_MONOTONE_FRONTIER',
        sources={p:digest(p)for p in sorted(paths)},baseline=H.TAG,scale=64,
        changed_pair=PAIR,guard='A whole unmarked A1 in the same finite corner; presence only',
        ordinary_A13_cost='77',guarded_A13_cost='78',local_proof=proof,current_interfaces=inter,
        ordinary_base_identities=len(H.C.TYPES),cross_identities=len(H.C.TYPES)**2,
        cost_identity='new=old-count(A13) if no A1; new=old otherwise',pointwise_nonincreasing=True,
        exact_actual_maximum='378',actual_attaining_forest=witness,
        removed_formal_plateau=dict(forest=plateau,previous_cost='378',new_cost='377'),
        complete_finite_upper_inherited=H.TAG,new_DP_or_obstruction_census=False,
        cross_costs_changed=False,new_cross_coupons=[],ideal_weights_changed=False,
        triangle_minimum='2',quadrilateral_minimum='1',positive_ideal_profiles=2552,actual_cusp_surplus='36',
        finite_curvature='377',strict_global_constant='4205',necessary_ambient58_inequality='v_f >= 4205 + 36*v_i',
        minimum_finite_vertices=4205,completed_dimension_bound=58,new_dimension_bound=None,
        scope='Fresh complete801/522 native selected search; four arbitrary-contact radical seals and one saturated two-A1 witness proof; current chosen-output retention and monotone sharp378 inheritance. No global DP, exact obstruction census, or dimension57 theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['exact_actual_maximum'],flush=True)
