#!/usr/bin/env python3
"""Independent Fraction replay of the three-free two-anchor cap dictionary.

No cap producer or floating-point arithmetic is imported. This checks saved
accepted caps only, not completeness of the numerical proposal family.
"""
import argparse
import json
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,dot,mv,pair,matrix,inverse,rref)
from scripts.fresh_h8_alternate_first_cap_review_20261013 import ea,es,em,sgn,pe,ve,ee
from scripts.fresh_h8_alternate_core_angles_independent_20261013 import load,bsign
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results/fresh_h8_threefree_independent_audit_20261013.json'
PERMITTED={Z:2,O:3,scale(O,2):4,scale(O,3):6,add(O,P):5,add(scale(O,2),P):10}
FOUR=scale(O,4)

def classify(value,d):
    c=sgn(es(value,(FOUR,Z)),d)
    if c>=0:return 'parallel' if c==0 else 'dotted'
    for target,label in PERMITTED.items():
        if sgn(es(value,(target,Z)),d)==0:return label
    raise AssertionError(('forbidden angle',value,d))

def run(verify=False):
    raw=json.loads((ROOT/'data/enumerations/fresh_h8_threefree_caps_20261013.json').read_text())
    base=json.loads((ROOT/raw['base_source']).read_text());Q=matrix(base['host_gram']);R=matrix(base['seed_normals'])
    G=matrix(base['seed_gram']);I=inverse(G);center=ve(base['strict_center_y'])
    assert G==tuple(tuple(pair(r,Q,s) for s in R) for r in R)
    assert I==matrix(base['seed_gram_inverse']) and len(rref(R)[1])==10
    assert all(G[i][i]==scale(O,2) for i in range(10))
    assert Q[0][0]==scale(P,-2) and all(Q[0][i]==Q[i][0]==Z for i in range(1,10))
    positive=[list(row[1:]) for row in Q[1:]]
    for i in range(9):
        t=positive[i][i];assert sign(t)>0
        for j in range(i+1,9):
            for k in range(i+1,9):positive[j][k]=sub(positive[j][k],div(mul(positive[j][i],positive[i][k]),t))
    assert all(sign(x)>0 for x in center) and sign(pair(center,I,center))<0
    graph,_,nodes=load();proof=json.loads((ROOT/raw['witness_source']).read_text())
    witnesses=[ve(r['escape_y']) for r in proof['five_complete_escape_clauses']]
    assert all(sign(pair(y,I,y))>0 and sign(pair(y,I,center))<0 and all(sign(t)>0 for t in y) for y in witnesses)
    old_by_dotted={}
    for i,n in enumerate(nodes):
        actual_q=(add(dot(n['a'],n['pa']),mul(n['d'],dot(n['b'],n['pb']))),scale(dot(n['a'],n['pb']),2))
        assert sgn(es(actual_q,(n['q'],Z)),n['d'])==0 and sign(n['q'])>0
        support=[]
        for j in range(10):
            s=em((n['pa'][j],n['pb'][j]),(n['pa'][j],n['pb'][j]),n['d'])
            s=tuple(mul(div(scale(O,2),n['q']),t) for t in s)
            if sgn(es(s,(FOUR,Z)),n['d'])>0:support.append(j)
        old_by_dotted.setdefault(tuple(support),[]).append(i)
    checks=[]
    for index,row in enumerate(raw['caps']):
        d=parse(row['discriminant']);assert sign(d)>0
        w0,w1=ve(row['line_constant']),ve(row['line_direction']);A,B,C=ve(row['norm_polynomial_descending'])
        assert (A,B,C)==(pair(w1,I,w1),scale(pair(w0,I,w1),2),sub(pair(w0,I,w0),scale(O,2)))
        assert d==sub(mul(B,B),scale(mul(A,C),4))
        x=pe(row['parameter']);assert sgn(ea(ea(em((A,Z),em(x,x,d),d),em((B,Z),x,d)),(C,Z)),d)==0
        ps=tuple((neg(add(a,mul(b,x[0]))),neg(mul(b,x[1]))) for a,b in zip(w0,w1))
        pa,pb=tuple(t[0] for t in ps),tuple(t[1] for t in ps)
        assert pa==ve(row['pairings_rational']) and pb==ve(row['pairings_radical'])
        ca,cb=mv(I,pa),mv(I,pb)
        assert ca==ve(row['coefficients_rational']) and cb==ve(row['coefficients_radical'])
        q=(add(dot(pa,ca),mul(d,dot(pb,cb))),scale(dot(pa,cb),2))
        assert sgn(es(q,(scale(O,2),Z)),d)==0
        assert sgn((dot(ca,center),dot(cb,center)),d)>0
        seed_labels=[]
        for p in ps:
            assert sgn(p,d)<=0;seed_labels.append(classify(em(p,p,d),d))
        dotted=[i for i,c in enumerate(seed_labels) if c=='dotted']
        assert sorted(row['free_seed_indices'])==dotted and len(dotted)==3
        assert row['free_original_walls']==[base['retained_original_walls'][i] for i in row['free_seed_indices']]
        fixed=ve(row['fixed_pairing_magnitudes'])
        assert sum(v!=Z for v in fixed)<=raw['fixed_pairing_maximum_support']
        for j,p in enumerate(ps):
            if j not in dotted:assert p==(neg(fixed[j]),Z)
        anchor_labels=[]
        for anchor,target in zip(row['anchor_node_indices'],ve(row['anchor_pairing_magnitudes'])):
            n=nodes[anchor];assert n['d']==Z and all(b==Z for b in n['b'])
            p=(dot(pa,n['a']),dot(pb,n['a']))
            assert sgn(ea(p,(target,Z)),d)==0 and sgn(p,d)<=0
            value=tuple(mul(div(scale(O,2),n['q']),t) for t in em(p,p,d))
            anchor_labels.append(classify(value,d))
        assert row['anchor_node_indices'][0]==348
        other=row['anchor_node_indices'][1];assert other in raw['second_anchor_node_indices']
        assert f'{other},348' in graph['edges']
        na,nb=nodes[other],nodes[348]
        pp=dot(na['pa'],nb['a']);assert sign(pp)<=0
        anchor_anchor=classify((div(scale(mul(pp,pp),4),mul(na['q'],nb['q'])),Z),d)
        vals=[(dot(ca,y),dot(cb,y)) for y in witnesses]
        assert vals==[pe(v) for v in row['core_constraint_values']]
        cuts=[j for j,v in enumerate(vals) if sgn(v,d)<0]
        assert cuts==row['cut_core_witness_indices'] and set(cuts)&set(raw['target_core_witnesses'])
        # Equal mirrors must have identical normalized squared seed pairings.
        # First discard old entries with a different strictly dotted support.
        old_matches=old_by_dotted.get(tuple(dotted),[]);separators=[]
        for oi in old_matches:
            n=nodes[oi];factor=div(scale(O,2),n['q']);different=None
            for j,p in enumerate(ps):
                a=em(p,p,d);b=em((n['pa'][j],n['pb'][j]),(n['pa'][j],n['pb'][j]),n['d'])
                diff=(sub(a[0],mul(factor,b[0])),a[1],neg(mul(factor,b[1])),Z)
                if bsign(diff,d,n['d'])!=0:different=j;break
            assert different is not None
            separators.append({'old_node':oi,'different_seed_index':different})
        checks.append({'cap':index,'norm':ee(q),'seed_labels':seed_labels,'anchor_nodes':row['anchor_node_indices'],
            'anchor_labels':anchor_labels,'anchor_anchor_label':anchor_anchor,'cut_core_witness_indices':cuts,
            'old_same_dotted_support_count':len(old_matches),'normalized_pairing_separators':separators})
    result={'status':'PASS','source':'data/enumerations/fresh_h8_threefree_caps_20261013.json',
        'arithmetic':'Independent Fraction arithmetic in Q(phi) and real quadratic extensions; exact biquadratic comparisons for old-mirror exclusion.',
        'caps_checked':len(checks),'old_dictionary_size':len(nodes),'seed_rank':10,'seed_inertia':[9,1],
        'all_caps_norm_two':True,'strict_common_timelike_center_y':base['strict_center_y'],
        'common_center_norm':enc(pair(center,I,center)),'checks':checks,
        'scope':'Every accepted cap is an exact new mirror outside the saved493 bank, compatible with the ten seed walls and its two stated compatible anchors, and cuts a recorded strict spacelike witness. This is an accepted-root audit, not a complete family search or a finite-volume result.'}
    if verify:assert json.loads(OUTPUT.read_text())==result
    else:OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','caps_checked','old_dictionary_size','scope')}));return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
