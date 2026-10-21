#!/usr/bin/env python3
"""Standalone exact replay of the deleted-wall two-height construction.

No producer, NumPy, optimizer, or SAT solver is imported. The arithmetic
comes from earlier independent Fraction verifiers, with fresh replay logic.
"""
import argparse,hashlib,json
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,total,dot,mv,pair,matrix,inverse,determinant,
)
from scripts.fresh_h8_release_h9_free_caps_20261013 import esign,emul,esub,bsign,bmul

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/enumerations/fresh_h8_twoheight_lift_20261013.json'
FREE=ROOT/'data/enumerations/fresh_h8_twoheight_free_caps_20261013.json'
OUTPUT=ROOT/'results/fresh_h8_twoheight_independent_audit_20261013.json'

def bsub(x,c):return (sub(x[0],c),*x[1:])
def relation(a,b,label,polys):
    d,e=a['d'],b['d'];v=(dot(a['pa'],b['a']),dot(a['pb'],b['a']),dot(a['pa'],b['b']),dot(a['pb'],b['b']))
    assert bsign(v,d,e)<=0
    ratio=tuple(mul(div(scale(O,4),mul(a['q'],b['q'])),x)for x in bmul(v,v,d,e))
    if label=='dotted':assert bsign(bsub(ratio,scale(O,4)),d,e)>0;return
    assert label!='parallel'
    label=int(label);poly,lo,hi=polys[label];value=(Z,Z,Z,Z)
    for coef in poly:
        value=bmul(value,ratio,d,e);value=(add(value[0],(coef,F(0))),*value[1:])
    assert bsign(value,d,e)==0
    assert bsign(bsub(ratio,(lo,F(0))),d,e)>=0 and bsign(bsub(ratio,(hi,F(0))),d,e)<=0

def escape_check(row,I,center,nodes):
    assert row['status']=='EXACT_STRICT_SPACELIKE_ESCAPE';y=tuple(map(parse,row['y']))
    assert all(sign(x)>0 for x in y)
    q=pair(y,I,y);cp=pair(y,I,center)
    assert sign(q)>0 and sign(cp)<0
    if 'norm'in row:assert parse(row['norm'])==q
    if 'center_pairing'in row:assert parse(row['center_pairing'])==cp
    for n in nodes:assert esign((dot(n['a'],y),dot(n['b'],y)),n['d'])>0
    return y

def check(base,free):
    old=matrix(json.loads((ROOT/base['source']).read_text())['wall_gram']);ids=base['retained_original_walls'];assert ids==[0,1,2,4,5,6,7,8,9,10]
    Q=add(scale(O,3),scale(P,2));left=(0,1,2,8);expected=[[old[i][j]for j in ids]for i in ids]
    for i in left:
        for j in left:expected[i][j]=mul(Q,expected[i][j])
    for i,j in ((0,1),(0,8),(8,9)):expected[i][j]=expected[j][i]=neg(Q)
    G=matrix(base['seed_gram']);I=matrix(base['seed_inverse']);assert G==tuple(map(tuple,expected)) and I==inverse(G)
    assert determinant(G)==parse(base['seed_determinant'])
    x0,x9=map(parse,base['height_squares']);assert x0==scale(P,2) and x9==div(add(scale(O,2),scale(P,8)),scale(O,11))
    n0,n9=add(scale(O,2),x0),add(scale(O,2),x9)
    assert div(scale(mul(P,P),4),scale(n0,2))==O
    assert div(scale(mul(x0,x9),4),mul(n0,n9))==O
    assert div(scale(mul(scale(P,2),scale(P,2)),4),scale(n9,2))==Q
    for conj in (False,True):
        A=[list(row)for row in G];signs=[]
        for k in range(10):
            p=A[k][k];assert p!=Z;signs.append(sign((p[0]+p[1],-p[1])if conj else p))
            for i in range(k+1,10):
                for j in range(k+1,10):A[i][j]=sub(A[i][j],div(mul(A[i][k],A[k][j]),p))
        assert signs.count(1)==9 and signs.count(-1)==1
    center=tuple(map(parse,base['strict_center_y']));assert all(sign(x)>0 for x in center)and sign(pair(center,I,center))<0
    source=json.loads((ROOT/'results/fresh_h8_alternate_core_angles_independent_20261013.json').read_text());polys={r['label']:(tuple(map(F,r['coefficients_descending'])),F(r['largest_root_interval'][0]),F(r['largest_root_interval'][1]))for r in source['real_cyclotomic_polynomials']}
    seednodes=[]
    for i in range(10):
        a=tuple(O if i==j else Z for j in range(10));seednodes.append({'q':G[i][i],'d':Z,'a':a,'b':(Z,)*10,'pa':G[i],'pb':(Z,)*10})
    labels={Z:2,O:3,scale(O,2):4,add(O,P):5,scale(O,3):6,add(scale(O,2),P):10}
    for i,j in combinations(range(10),2):
        value=div(scale(mul(G[i][j],G[i][j]),4),mul(G[i][i],G[j][j]));label='dotted'if sign(sub(value,scale(O,4)))>0 else labels[value];relation(seednodes[i],seednodes[j],label,polys)
    nodes=[]
    for row in base['caps']:
        n={'q':parse(row['norm']),'d':Z,'a':tuple(map(parse,row['coefficients'])),'b':(Z,)*10,'pa':tuple(map(parse,row['pairings'])),'pb':(Z,)*10}
        assert mv(G,n['a'])==n['pa'] and dot(n['a'],n['pa'])==n['q'] and sign(n['q'])>0
        assert sign(dot(n['a'],center))>0
        for i,label in enumerate(row['seed_angles']):relation(n,seednodes[i],label,polys)
        nodes.append(n)
    assert len({n['a']for n in nodes})==len(nodes)
    assert all(n['q']in(scale(O,2),scale(Q,2))for n in nodes)
    assert mul(Q,(Q[0]+Q[1],-Q[1]))==scale(O,11)
    # Norm(Q)=11 is nonsquare, so different norm types cannot be rescalings
    # by an element of Q(phi); same-type equality was checked above.
    initial=len(nodes);y=escape_check(base['all_caps_escape'],I,center,nodes);assert tuple(map(parse,free['target_y']))==y
    base_pairs={(i,j):label for i,j,label in base['compatible_pairs']}
    for row in base['selections']:
        chosen=row['caps'];assert row['facets']==10+len(chosen)
        for i,j in combinations(chosen,2):relation(nodes[i],nodes[j],base_pairs[i,j],polys)
        escape_check(row['volume_gate'],I,center,[nodes[i]for i in chosen])
    for row in free['caps']:
        n={'q':parse(row['norm']),'d':parse(row['discriminant']),'a':tuple(map(parse,row['coefficients_rational'])),'b':tuple(map(parse,row['coefficients_radical'])),'pa':tuple(map(parse,row['pairings_rational'])),'pb':tuple(map(parse,row['pairings_radical']))}
        assert sign(n['d'])>0 and sign(n['q'])>0
        assert mv(G,n['a'])==n['pa']and mv(G,n['b'])==n['pb']
        q=(add(dot(n['a'],n['pa']),mul(n['d'],dot(n['b'],n['pb']))),add(dot(n['a'],n['pb']),dot(n['b'],n['pa'])))
        assert esign(esub(q,(n['q'],Z)),n['d'])==0
        assert esign((dot(n['a'],center),dot(n['b'],center)),n['d'])>0
        assert esign((dot(n['a'],y),dot(n['b'],y)),n['d'])<0
        for i,label in enumerate(row['seed_angles']):relation(n,seednodes[i],label,polys)
        nodes.append(n)
    for row in free['cases']:
        chosen=row['chosen_indices'];assert row['facets']==10+len(chosen)
        for i,j in combinations(chosen,2):
            rec=free['exact_pair_checks'][f'{i},{j}'];assert rec['compatible'];relation(nodes[i],nodes[j],rec['label']if rec['kind']=='finite'else rec['kind'],polys)
        escape_check(row['volume_gate'],I,center,[nodes[i]for i in chosen])
    escape_check(free['all_dictionary_escape'],I,center,nodes)
    # Different target signs exclude every new mirror from the old dictionary.
    aliases=[]
    for i,j in combinations(range(initial,len(nodes)),2):
        a,b=nodes[i],nodes[j];equal=True
        for k in range(10):
            aa=emul((a['pa'][k],a['pb'][k]),(a['pa'][k],a['pb'][k]),a['d']);bb=emul((b['pa'][k],b['pb'][k]),(b['pa'][k],b['pb'][k]),b['d'])
            diff=(sub(div(aa[0],a['q']),div(bb[0],b['q'])),div(aa[1],a['q']),neg(div(bb[1],b['q'])),Z)
            if bsign(diff,a['d'],b['d'])!=0:equal=False;break
        if equal:aliases.append([i,j])
    return {'base_caps':initial,'new_cap_entries':len(nodes)-initial,'new_mirror_alias_pairs':aliases,'base_configurations':len(base['selections']),'new_configurations':len(free['cases']),'all_dictionary_escape_verified':True,'facet_range':[min(r['facets']for r in base['selections']+free['cases']),max(r['facets']for r in base['selections']+free['cases'])]}

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();base=json.loads(BASE.read_text());free=json.loads(FREE.read_text());counts=check(base,free)
    mutations=[]
    for name,change in [('wrong_seed_pairing',lambda b,f:b['seed_gram'][0].__setitem__(8,enc(Z))),('timelike_escape',lambda b,f:b['all_caps_escape'].__setitem__('y',b['strict_center_y']))]:
        b,f=deepcopy(base),deepcopy(free);change(b,f)
        try:check(b,f)
        except (AssertionError,KeyError):mutations.append(name)
        else:raise AssertionError(f'mutation accepted: {name}')
    result={'status':'PASS_INDEPENDENT_EXACT_TWOHEIGHT_REPLAY',**counts,'rejected_mutations':mutations,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(BASE,FREE)}}
    if a.verify:assert result==json.loads(OUTPUT.read_text())
    else:OUTPUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
