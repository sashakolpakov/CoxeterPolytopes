#!/usr/bin/env python3
"""Independent SymPy review of the four remaining cases and coefficient10 binder.

No producer arithmetic imports. Literal graphs and recorded results are inputs;
all critical inverses and residual matrices are independently reconstructed.
"""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_finite_edge_cusp_tradeoff10_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
EASY='fresh_finite_edge_seven_simplex_20261017'
HARD='fresh_seven_simplex_affine_and_d8_20261017'
PREV='fresh_finite_edge_cusp_tradeoff12_20261017'

def load(t):return json.loads((ROOT/f'results/{t}.json').read_text())
def gram(L):
    G=2*s.eye(L['order'])
    for a,b,m in L['edges']:G[a,b]=G[b,a]=-1 if m==3 else -s.sqrt(2)
    return G
def rational(x):
    x=s.simplify(x);assert x.is_Rational;return str(x)
def encoded(x):
    x=s.expand(x);a=x.coeff(s.sqrt(2),0);b=x.coeff(s.sqrt(2))
    assert s.simplify(x-a-b*s.sqrt(2))==0 and a.is_Rational and b.is_Rational
    return [str(a),str(b)]
def data(M):return [[encoded(M[i,j]) for j in range(M.cols)] for i in range(M.rows)]
def psd2(M):return M[0,0]>=0 and M[1,1]>=0 and s.simplify(M.det())>=0
def costs(B,profiles):
    out=[]
    for b in range(B.rows):
        ids=[i for i in range(B.rows) if i!=b];H=B.extract(ids,ids)
        assert all(H[:k,:k].det()>0 for k in range(1,H.rows+1))
        J=H.inv();assert H*J==s.eye(H.rows)
        C=s.Matrix([[p[i] for p in profiles] for i in ids])
        out.append(s.simplify(C.T*J*C))
    return out

def independent_review(census,hard):
    D=census['diagrams'];aff=[]
    expected=[(20,'B7',[(0,),(1,)]),(21,'D7',[(0,),(1,),(6,),(7,)]),
              (22,'E7',[(4,),(7,),(4,7)])]
    for index,name,supports in expected:
        record=next(r for r in hard['affine_black_cases'] if r['diagram']==index)
        B=gram(D[index])[:8,:8]
        assert B.det()==0 and B.rank()==7
        profiles=[[int(i in p) for i in range(8)] for p in supports]
        provenance=[p for p in census['attachment_provenances'] if p['affine_source']==name]
        assert sorted(p['contacts'] for p in provenance)==sorted(
            [3 if i in p else 2 for i in range(8)] for p in supports)
        for p in provenance:
            mp=dict(p['isomorphism']);G=gram(D[p['diagram']])
            assert G.extract([mp[i] for i in range(8)],[mp[i] for i in range(8)])==B
        K=costs(B,profiles);assert [data(k) for k in K]==record['cost_matrices']
        singles=[]
        for i in range(len(profiles)):
            c=[k[i,i] for k in K];hi=max(c);ties=[b for b,v in enumerate(c) if v==hi]
            assert hi>1 and hi<=2
            rec=record['single_profiles'][i]
            assert rec['costs']==[rational(v) for v in c]
            assert rec['outside_diagonal_maximum']==rational(2-hi)
            assert rec['possible_ideal_blacks']==ties
            b=c.index(hi);duplicate=2*s.eye(2)-hi*s.ones(2)
            assert duplicate.det()<0
            assert encoded(duplicate.det())==record['repeated_profile_exclusions'][i]['determinant']
            singles.append(dict(profile=i,maximum_cost=rational(hi),maximizing_black=ties))
        pair_records=[]
        for i,j in combinations(range(len(profiles)),2):
            rec=next(p for p in record['pair_profiles'] if p['profiles']==[i,j])
            H0=[2*s.eye(2)-k.extract([i,j],[i,j]) for k in K]
            bad=[b for b,H in enumerate(H0) if not psd2(H)]
            if bad:
                assert rec['allowed_native_magnitudes']==[]
                b=rec['exclusion_black'];assert H0[b].det()<0
                assert encoded(H0[b].det())==rec['determinant_at_right']
                pair_records.append(dict(profiles=[i,j],excluded_at_black=b));continue
            H1=[H-s.Matrix([[0,1],[1,0]]) for H in H0]
            allowed=[0,1] if all(psd2(H) for H in H1) else [0]
            assert allowed==rec['allowed_native_magnitudes']
            b=rec['all_label_cap_black'];cap=s.Rational(rec['all_label_cap'])
            H=H0[b]-s.Matrix([[0,cap],[cap,0]])
            assert H.det()==0 and K[b][i,j]>0
            statuses=[]
            for a in allowed:
                HH=H0 if a==0 else H1;ideal=[b for b,H in enumerate(HH) if H.det()==0]
                rr=next(r for r in rec['statuses'] if r['native']==a)
                assert ideal==rr['original_ideal_blacks']
                if ideal:
                    b=rr['positive_null_at'];ker=HH[b].nullspace();assert len(ker)==1
                    v=ker[0]
                    assert all(x>0 for x in v) or all(x<0 for x in v)
                    assert rr['maximum_ideals']==len(ideal)
                else:assert index==22 and [i,j]==[0,1] and a==0
                statuses.append(dict(native=a,ideal_black=ideal))
            pair_records.append(dict(profiles=[i,j],allowed_native_magnitudes=allowed,statuses=statuses))
        aff.append(dict(diagram=index,single_costs=singles,pairs=pair_records))

    B=gram(D[22])[:8,:8];K=costs(B,[[int(i==j) for i in range(8)] for j in (4,7)])
    H=[2*s.eye(2)-k for k in K];assert [data(h) for h in H]==hard['E7_right_pair_outside_white_certificate']['residuals']
    half=s.Rational(1,2);bounds=[]
    assert H[4]==s.diag(2,half) and H[7]==s.diag(half,2)
    for b in (0,2,3,5,6):
        assert H[b][0,1]==0
        lower=s.factor((H[b][0,0]-half)*(H[b][1,1]-half)-half**2);assert lower>0
        bounds.append([b,str(lower)])
    assert bounds==[[0,'5/16'],[2,'11/36'],[3,'1/4'],[5,'11/36'],[6,'1/4']]
    assert H[1]==s.Matrix([[s.Rational(9,8),-s.Rational(1,8)],[-s.Rational(1,8),s.Rational(9,8)]])
    assert (H[4]-s.ones(2)/2).det()==-s.Rational(1,4)

    L=gram(D[23]);mp=[2,0,1,9,3,4,5,6,7];E=L.extract(mp,mp)
    assert E==gram(D[25])[:9,:9]
    assert [L[8,i] for i in mp]==[0]*8+[-s.sqrt(2)]
    profiles=[[int(i==8) for i in range(9)],
              [s.sqrt(2)*int(i==8) for i in range(9)],[int(i==7) for i in range(9)]]
    K=costs(E,profiles);dets=[s.factor((2*s.eye(2)-K[1].extract([1,p],[1,p])).det()) for p in range(3)]
    assert dets==[-s.Rational(4,3),-s.Rational(28,9),-s.Rational(10,9)]
    assert data(K[1])==hard['finite_D8_case']['delete1_cost_gram']
    c=[[j,str(K[j][1,1])] for j in range(9) if j!=3]
    assert c==hard['finite_D8_case']['existing_profile_costs']
    assert max(s.Rational(v) for j,v in c)==s.Rational(16,9)
    assert [j for j,v in c if s.Rational(v)==s.Rational(16,9)]==[1]

    rank9=[]
    for rec in hard['finite_black_rank9_alternatives']:
        L=D[rec['diagram']];w=rec['internal_white'];G=gram(L)
        B=[i for i in range(9) if i!=w];assert G.extract(B,B).det()>0
        I=L['ideal_deletions'];assert w not in I
        if len(I)>=2:mx=len(I)
        else:
            a=I[0];ts=[]
            for j in B:
                if j==a:continue
                ids=[i for i in range(9) if i!=j];J=G.extract(ids,ids).inv()
                ts.append([j,rational(1/J[ids.index(a),ids.index(a)])])
            assert ts==rec['thresholds'];low=min(s.Rational(v) for _,v in ts)
            ties=[j for j,v in ts if s.Rational(v)==low]
            assert ties==rec['ties'];mx=1+len(ties)
        assert mx==rec['maximum_ideals']<=3;rank9.append([rec['diagram'],w,mx])
    assert len(rank9)==31
    return dict(affine_profile_and_pair_reconstructions=aff,
        E7_right_pair_residuals=[data(h) for h in H],E7_strict_determinant_bounds=bounds,
        E7_excluded_equality_black4_determinant='-1/4',
        D8_standard_to_L=mp,D8_delete1_cost_gram=data(K[1]),
        D8_pair_determinants=list(map(str,dets)),D8_costs=c,
        finite_black_rank9_alternatives=rank9)

def recursive_sources(paths):
    checked={};todo=list(paths);seen_json=set()
    while todo:
        p=todo.pop();full=ROOT/p;assert full.exists(),p
        digest=hashlib.sha256(full.read_bytes()).hexdigest()
        if p in checked:assert checked[p]==digest
        checked[p]=digest
        if p.endswith('.json') and p not in seen_json:
            seen_json.add(p);d=json.loads(full.read_text())
            for child,h in d.get('sources',{}).items():
                assert isinstance(h,str) and hashlib.sha256((ROOT/child).read_bytes()).hexdigest()==h,child
                todo.append(child)
    return checked

def run():
    census=load(CENSUS);easy=load(EASY);hard=load(HARD);prev=load(PREV)
    assert easy['status']=='PASS_EXACT_SEVEN_SIMPLEX_162_EASY_AND_NINE_FINITE_BLACK_REDUCTIONS'
    assert hard['status']=='PASS_FOUR_HARD_SEVEN_SIMPLEX_CASES_AT_MOST_THREE_IDEAL'
    assert prev['status']=='PASS_EXACT_EIGHT_SIMPLEX_AT_MOST_THREE_IDEAL_AND_FINITE_EDGE_COEFFICIENT12'
    review=independent_review(census,hard);rows=[]
    special={(20,(8,)):1,(21,(8,)):2,(22,(8,)):3,(23,(8,9)):3}
    for case in easy['cases']:
        key=(case['diagram'],tuple(case['internal_whites']))
        if 'possible_ideal_counts' in case:counts=case['possible_ideal_counts'];method=case['method']
        elif key in special:counts=list(range(special[key]+1));method='four_remaining_hard_cases'
        else:
            assert key[0]==25 and key[1][-1]==9
            counts=[0,1,2,3];method='finite_black_E8_reduction'
        assert all(i<=3 or i==8 for i in counts)
        finite_edge_counts=[i for i in counts if 8-i>=2]
        assert all(i<=3 for i in finite_edge_counts)
        rows.append(dict(diagram=key[0],internal_whites=list(key[1]),
                         method=method,allowed_ideal_count_envelope=counts,
                         allowed_for_finite_edge_blocks=finite_edge_counts))
    assert len(rows)==175 and len({(r['diagram'],tuple(r['internal_whites'])) for r in rows})==175
    coefficient=[]
    for vertices in range(2,11):
        for f in range(2,vertices+1):
            i=vertices-f
            if vertices in (8,9) and i>3:continue
            if vertices==10 and (f,i) not in [(7,3),(8,2),(9,1)]:continue
            c=2*(F(i)+F(f-2,3));assert c<=10
            coefficient.append(dict(dimension=vertices-1,f=f,i=i,coefficient=str(c)))
    assert [r for r in coefficient if r['coefficient']=='10']==[dict(dimension=6,f=2,i=5,coefficient='10')]
    paths=[f'scripts/{TAG}.py']
    for t in (EASY,HARD,PREV):paths += [f'scripts/{t}.py',f'results/{t}.json',f'results/{t}.md']
    paths+=['results/fresh_finite_ended_simplex_ports_20261017.md',
            'results/fresh_finite_volume_simplex_truncation_20261017.md']
    sources=recursive_sources(paths)
    return dict(status='PASS_COMPLETE_SEVEN_SIMPLEX_IDEAL_DICHOTOMY_AND_FINITE_EDGE_COEFFICIENT10',
        scope='Complete internally reviewed lemma for actual simplex faces and finite-edge triangle accounting; no ambient dimension bound.',
        sources=sources,independent_four_case_review=review,complete_seven_simplex_cases=rows,
        seven_simplex_case_count=175,seven_simplex_ideal_count_envelope=[0,1,2,3,8],
        all_ideal_boundary_is_retained=True,coefficient_cases=coefficient,
        uniform_finite_edge_coefficient='10',formal_unresolved_equality=dict(dimension=6,f=2,i=5),
        inequality='2(T0+T1) <= 10 E_ff = 5 n v_f - 5 E_fi <= 5 n v_f',
        residual_budget='2(T2+T3)+Q >= (b_n-5n)v_f+5E_fi+alpha_n sum E(L_v)',
        independent_review_scope='Separate SymPy inverses/costs/pair determinants/E7 endpoint inequality/D8 map and31rank9 alternatives, no producer arithmetic imports; not a new replay of the146262-attachment census.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(out.read_text())==d
    else:out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],'saved_replay=',args.verify)
