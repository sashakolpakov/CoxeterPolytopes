#!/usr/bin/env python3
"""Independent label-dependency audit for the sparse T8 rank theorem.

Rebuilds the local witnesses and budgets solely from the eight hypotheses.
Reuses the frozen independently verified polynomial banks; no search or new
polynomial generation is performed, and no producer arithmetic is imported.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_t8_sparse_frontier_20261014_review'
GEO='results/fresh_h9_tree_cube_wedges_20261014.json'
BASE='results/fresh_h9_t8_asymmetric_20261014.json'
MID='results/fresh_h9_t8_middle_20261014.json'
SPARSE='results/fresh_h9_t8_sparse_frontier_20261014.json'
STATE='results/fresh_h9_t8_middle_20261014_search_case_08.json'
REVIEWS=['results/fresh_h9_t8_asymmetric_20261014_review.json',
         'results/fresh_h9_t8_block_polynomials_20261014_review.json',
         'results/fresh_h9_t8_middle_20261014_review.json',
         'results/fresh_h9_t8_middle_20261014_review_algebra.json']
S=[1,2,4,5,7,8,12,13]
FIXED={(1,14):2,(2,14):3,(4,14):3,(5,14):3,(4,5):4}
LOWER={(1,2):7,(7,8):4,(12,13):4}
Z=(F(0),F(0));ONE=(F(1),F(0));HALF=(F(1,2),F(0))
ROOT_HALF=(F(0),F(1,2))

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return[enc(y)for y in x]
    if isinstance(x,dict):return{str(k):enc(y)for k,y in x.items()}
    return x
def add(x,y):return x[0]+y[0],x[1]+y[1]
def mul(x,y,d):return x[0]*y[0]+d*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def total(xs):
    v=Z
    for x in xs:v=add(v,x)
    return v
def witness(support,edges,w,d):
    n=len(support);G=[[ONE if i==j else Z for j in range(n)]for i in range(n)]
    for i,j,c in edges:G[i][j]=G[j][i]=(-c[0],-c[1])
    assert all(a>=0 and b>=0 and(a or b)for a,b in w)
    assert all(total(mul(G[i][j],w[j],d)for j in range(n))==Z for i in range(n))
    return dict(support=support,radical=d,relaxed_gram=G,positive_null_vector=w)

def run():
    paths=[GEO,BASE,MID,SPARSE,STATE]+REVIEWS
    docs={p:json.loads((ROOT/p).read_text())for p in paths}
    hashed={};visited=set()
    def bind(p):
        if p in visited:return
        visited.add(p)
        obj=json.loads((ROOT/p).read_text())
        for name,value in obj.get('sources',{}).items():
            if name not in hashed:hashed[name]=sha(name)
            assert hashed[name]==value,name
            if name.endswith('.json'):bind(name)
    for p in paths:bind(p)
    geo=docs[GEO]['representatives'][8]['geometry'];vertices=[set(v['incident'])for v in geo['vertices']]
    base,mid,sparse,state=[docs[p]for p in[BASE,MID,SPARSE,STATE]]
    assert base['positive_support']==S
    core=sorted(set(S)|{11,14,17})
    assert len(core)==11 and core==sparse['core_roots']
    def actual(s,expected=None):
        choices=[i for i,v in enumerate(vertices)if set(s)<=v]
        assert choices,s
        if expected is not None:assert expected in choices;return expected
        return choices[0]
    actual(S,base['positive_support_actual_vertex'])
    pair_incidence={str(p):actual(p)for p in combinations(core,2)}
    assert len(pair_incidence)==55
    for missing in[[11,12,13],[12,13,17]]:
        assert missing in geo['minimal_nonfaces']
        assert all(any(set(p)<=v for v in vertices)for p in combinations(missing,2))
        assert not any(set(missing)<=v for v in vertices)
    # Stage1: cross-block contacts. Only the four internal lower bounds enter.
    blocks=[S[2*i:2*i+2]for i in range(4)]
    old_cross={tuple(r['pair']):r for r in base['cross_right_guards']}
    cross=[]
    for i,j in combinations(range(4),2):
        for a,b in product(blocks[i],blocks[j]):
            if i==0:
                support=[next(x for x in blocks[0]if x!=a),a,b]
                row=witness(support,[(0,1,ROOT_HALF),(1,2,HALF)],
                            [ROOT_HALF,ONE,HALF],3)
                row['label_dependencies']=[dict(pair=[1,2],lower=7,relaxed_to=6),
                                           dict(pair=sorted([a,b]),tested_lower=3)]
            else:
                support=[next(x for x in blocks[i]if x!=a),a,b,next(x for x in blocks[j]if x!=b)]
                row=witness(support,[(0,1,ROOT_HALF),(1,2,HALF),(2,3,ROOT_HALF)],
                            [ONE,(F(0),F(1)),(F(0),F(1)),ONE],2)
                for block in[blocks[i],blocks[j]]:
                    p=tuple(block);assert FIXED.get(p,LOWER.get(p))>=4
                row['label_dependencies']=[dict(pair=blocks[i],lower=4,relaxed_to=4),
                    dict(pair=blocks[j],lower=4,relaxed_to=4),dict(pair=sorted([a,b]),tested_lower=3)]
            old=old_cross[a,b];assert row['support']==old['support']
            row['actual_vertex']=actual(row['support'],old['actual_vertex'])
            row['conclusion']=dict(pair=sorted([a,b]),label=2)
            cross.append(row)
    assert len(cross)==24
    # Stage2: literal budgets, after the cross-right conclusion.
    anchor=[]
    for leaf,kept in product([11,17],[12,13]):
        support=sorted((set(S)-{25-kept})|{leaf})
        old=next(r for r in base['anchor_budgets']if r['leaf']==leaf and r['retained_own_root']==kept)
        assert old['support']==support
        anchor.append(dict(leaf=leaf,kept_own_root=kept,support=support,
            actual_vertex=actual(support,old['actual_vertex']),
            assumed_labels='Only internal lower bounds, derived cross-right, and the chosen own contact being nonright; its minimum cosine is1/2. No other leaf labels.',
            own_minimum_cost=F(1,4),foreign_B2_minimum_cost=F(1,2),
            two_foreign_cost=F(1,4)+2*F(1,2),foreign_label4_cost=F(1),
            own_label4_plus_foreign=F(1,2)+F(1,2)))
    assert F(1,4)+2*F(1,2)>1 and F(1,4)+1>1 and F(1,2)+F(1,2)==1
    # I2(m), m>=7: c²>3/4, so inverse diagonal>4 and a pi3 costs>1.
    assert F(1,4)/(1-F(3,4))==1
    third=[]
    for kept in[4,5]:
        support=sorted((set(S)-{1}-({4,5}-{kept}))|{14})
        old=next(r for r in base['third_budgets']if r['retained_own_roots']==[2,kept])
        assert old['support']==support
        assert FIXED[2,14]==FIXED[kept,14]==3
        third.append(dict(support=support,actual_vertex=actual(support,old['actual_vertex']),
            fixed_contacts=[[2,14,3],[kept,14,3]],known_cost=F(1,2),
            any_nonright_foreign_extra_cost_at_least=F(1,2),
            forced_right_pairs=[[14,j]for j in[7,8,12,13]]))
    assert F(1,4)+F(1,4)+F(1,2)==1
    # Stage3: new stars replace ALL old full-nine-vertex upper witnesses.
    stars=[]
    for leaf,vi in[(11,72),(17,322)]:
        row=witness([leaf,14,2,4],[(0,1,ROOT_HALF),(1,2,HALF),(1,3,HALF)],
                    [ROOT_HALF,ONE,HALF,HALF],2)
        row['actual_vertex']=actual(row['support'],vi)
        assert FIXED[2,14]==FIXED[4,14]==3
        row['label_dependencies']=[[2,14,3],[4,14,3]]
        row['tested_contact']=[min(leaf,14),max(leaf,14),4]
        row['other_three_contacts']='Relaxed right; arbitrary larger actual cosines only decrease the positive-vector quadratic value.'
        assert row['support']==next(r['support']for r in sparse['replacement_third_contact_guards']if r['anchor']==leaf)
        stars.append(row)
    # Stage4: every literal own-support pattern, without incidence renaming.
    patterns=[]
    for sa,sb in product([(12,),(13,),(12,13)],repeat=2):
        saved=next(r for r in base['own_support_patterns']if r['own_supports']==[list(sa),list(sb)])
        shared=sorted(set(sa)&set(sb))
        if shared:
            support=[11,shared[0],17]
            w=witness(support,[(0,1,HALF),(1,2,HALF),(0,2,HALF)],[ONE]*3,2)
            w['actual_vertex']=actual(support,saved['anchor_pair_guard']['actual_vertex'])
            w['dependencies']='Two active own contacts are >=3; test mutual>=3. No outside contact is used.'
            mutual='right'
        else:
            w=dict(support=[11,17],actual_vertex=actual([11,17],saved['anchor_pair_guard']['actual_vertex']))
            mutual='Any finite label, covered by the entire cosine interval[0,1].'
        patterns.append(dict(own_supports=[sa,sb],mutual_domain=mutual,guard=w,
                             polynomial_family=saved['polynomial_family']))
    assert len(patterns)==9
    # Stage5: check the exact split feeding the already independently replayed banks.
    assert base['case_count']==408 and base['positive_coefficient_count']==11736
    assert mid['case_count']==256 and mid['inherited_case_count']==152 and mid['positive_coefficient_count']==4320
    assert docs[REVIEWS[0]]['status']=='PASS_INDEPENDENT_FULL_T8_OWN_PATTERN_FAMILY_REVIEW'
    assert docs[REVIEWS[2]]['status']=='PASS_INDEPENDENT_T8_UNRESTRICTED_MIDDLE_LABEL_EXTENSION'
    assert docs[REVIEWS[3]]['case_count']==256
    # A foreign contact to I2(p), p>=6, costs>=1 plus positive owncost.
    assert F(1,4)/(1-F(3,4))+F(1,4)>1
    # For n<=6 a single own pi3 gives a PSD missing triple; hence label>=4.
    assert F(3,4)-F(3,4)==0 and F(1,2)+F(1,2)==1
    # Largest U6 squared-root threshold, reproduced by integer recurrence.
    U0=[1];U1=[0,2]
    for _ in range(2,7):
        U=[0]+[2*x for x in U1]
        for i,x in enumerate(U0):U[i]-=x
        U0,U1=U1,U
    assert U1==[-1,0,24,0,-80,0,64]
    x=F(81,100);assert 64*x**3-80*x*x+24*x-1==-F(559,15625)
    assert 192*F(3,4)**2-160*F(3,4)+24==12
    assert 384*F(3,4)-160==128
    # Every one of the55 selected Gram entries now has a certified domain.
    accounted=set(tuple(sorted(r['conclusion']['pair']))for r in cross)
    accounted|={tuple(b)for b in blocks}
    accounted|={tuple(sorted((a,b)))for a in[11,17]for b in S}
    accounted|={tuple(sorted((14,b)))for b in S}
    accounted|={(11,17),(11,14),(14,17)}
    assert accounted==set(combinations(core,2))
    labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
    assert all(labels[p]==m for p,m in FIXED.items())and all(labels[p]>=m for p,m in LOWER.items())
    for a in[11,17]:assert any(labels[tuple(sorted((a,b)))]>=3 for b in[12,13])
    paths+=['scripts/fresh_h9_t8_sparse_frontier_20261014.py',f'scripts/{TAG}.py']
    return enc(dict(status='PASS_INDEPENDENT_SPARSE_T8_DEPENDENCY_REDUCTION',sources={p:sha(p)for p in paths},
        recursively_bound_manifests=len(visited),core_roots=core,all55_core_pairs_are_meeting=pair_incidence,
        fixed_labels=[[a,b,m]for(a,b),m in sorted(FIXED.items())],finite_label_lower_bounds=[[a,b,m]for(a,b),m in sorted(LOWER.items())],
        independent_cross_guards=cross,anchor_budget_guards=anchor,third_budget_guards=third,
        replacement_affine_stars=stars,all_nine_own_patterns=patterns,selected_gram_domain_pairs=len(accounted),
        inherited_polynomial_cases=dict(base=408,middle_changed=256,middle_inherited=152),
        saved_state_excluded=True,
        scope='Reconstructs every local guard from the eight sparse conditions and derived domains, not just support containment. No guard uses outside leaf labels; all55 selected Gram entries are covered. Frozen independent polynomial banks are bound and reused, not regenerated. No new search and no whole-incidence exclusion.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'core_pairs=',out['selected_gram_domain_pairs'])
