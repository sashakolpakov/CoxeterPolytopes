#!/usr/bin/env python3
"""Separate exact domains and polynomial replay for two high T8 labels."""
import argparse,json
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
from fresh_h9_star_cube_metric_20261014_cut_profiles_review import (
    F,Z,O,CO,add,sub,neg,mul,scale,total,inverse,enc)
from fresh_h9_star_cube_metric_20261014_all_cut_profiles_review import sign,det_gauss
from fresh_h9_t8_metric_20261014_review import q2,padd,pmul,pdet
from fresh_h9_t8_metric_20261014_two_labels_review import mpdet,mpeval
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_metric_20261014_two_high'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    src=f'results/{TAG}.json';data=json.loads((ROOT/src).read_text())
    assert data['status']=='PASS_T8_TWO_HIGH_LABELS_M_N_GE_7_EXCLUDED'
    bound=set()
    def bind(p):
        if p in bound:return
        bound.add(p)
        for s,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():
            assert sha(s)==h,s
            if s.endswith('.json'):bind(s)
    bind(src)
    gp='results/fresh_h9_tree_cube_wedges_20261014.json'
    sp='results/fresh_h9_t8_metric_20261014_family_search_case_08.json'
    g=json.loads((ROOT/gp).read_text())['representatives'][8]['geometry']
    state=json.loads((ROOT/sp).read_text());assert state['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
    labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
    assert {p for p,m in labels.items()if m==7}=={(1,2),(12,13)}
    assert labels[11,12]==labels[12,17]==2 and labels[11,13]==labels[13,17]==3
    V=[v['incident']for v in g['vertices']]
    S=[1,2,4,5,7,8,12,13];blocks=[S[2*i:2*i+2]for i in range(4)]
    assert data['positive_support']==S and set(S)<=set(V[data['positive_support_actual_vertex']])
    assert [labels[tuple(p)]for p in blocks]==[7,4,4,7]
    cross={(a,b)for i,j in combinations(range(4),2)for a,b in product(blocks[i],blocks[j])}
    seen=set();n_high=n_relaxed=0
    zero,one,x=(Z,),(O,),(Z,O)
    for row in data['forced_cross_right_witnesses']:
        a,b=row['pair'];assert(a,b)in cross and(a,b)not in seen;seen.add((a,b))
        s=row['support'];assert set(s)<=set(V[row['actual_vertex']])
        i=S.index(a)//2;j=S.index(b)//2
        if i==0:
            assert s==[next(v for v in blocks[0]if v!=a),a,b]
            H=[[one,(Z,neg(O)),zero],[(Z,neg(O)),one,(neg(CO[3]),)],
               [zero,(neg(CO[3]),),one]]
            assert pdet(H)==(scale(O,F(3,4)),Z,neg(O));n_high+=1
        else:
            assert s==[next(v for v in blocks[i]if v!=a),a,b,next(v for v in blocks[j]if v!=b)]
            assert labels[tuple(blocks[i])]>=4 and labels[tuple(blocks[j])]>=4
            H=[[O if a==b else Z for b in range(4)]for a in range(4)]
            for k,c in enumerate([CO[4],CO[3],CO[4]]):H[k][k+1]=H[k+1][k]=neg(c)
            null=[O,scale(CO[4],2),scale(CO[4],2),O]
            assert all(total(mul(a,b)for a,b in zip(r,null))==Z for r in H)
            assert H==[[q2(v)for v in r]for r in row['gram']]
            assert null==list(map(q2,row['null_vector']));n_relaxed+=1
    assert seen==cross and n_high==n_relaxed==12
    own=data['anchor_pair_right_witness'];assert own['support']==[11,13,17]
    assert own['pair']==[11,17]and own['test_label']==3
    assert set(own['support'])<=set(V[own['actual_vertex']])
    H=[[O if i==j else neg(CO[3])for j in range(3)]for i in range(3)]
    assert H==[[q2(v)for v in r]for r in own['gram']]and all(total(r)==Z for r in H)
    def gram(nodes,lab):
        assert not {1,2}<=set(nodes)and not {12,13}<=set(nodes)
        return [[O if a==b else neg(CO[lab.get(tuple(sorted((a,b))),2)])for b in nodes]for a in nodes]
    b2labels={(4,5):4,(7,8):4}
    budgets=[]
    for row in data['anchor_budget_witnesses']:
        leaf=row['leaf'];assert leaf in[11,17]and row['retained_own_root']==13
        assert set(row['support'])==(set(S)-{12})|{leaf}<=set(V[row['actual_vertex']])
        assert labels[tuple(sorted((leaf,13)))]==3
        J=inverse(gram([4,5,7,8,13],b2labels))
        assert J[-1][-1]==O and all(J[k][k]==scale(O,2)for k in range(4))
        assert all(sign(v)>=0 for r in J for v in r)
        budgets.append(row)
    assert len(budgets)==2 and {r['leaf']for r in budgets}=={11,17}
    # First I2 inverse diagonal>4 forbids all anchor contacts there; the
    # retained root13 costs1/4, leaving at most one foreign B2 label3.
    assert F(1,4)/(1-F(3,4))==1
    assert F(1,4)+2*F(1,4)*2>1 and F(1,4)+2*F(1,2)>1
    third=[]
    for row in data['third_leaf_budget_witnesses']:
        a,kept=row['retained_own_roots'];assert a==2 and kept in[4,5]and row['leaf']==14
        assert set(row['support'])=={2,kept,7,8,12,13,14}<=set(V[row['actual_vertex']])
        assert labels[2,14]==labels[kept,14]==3
        third.append(row)
    assert len(third)==2 and {r['retained_own_roots'][1]for r in third}=={4,5}
    assert F(1,4)+F(1,4)+2*F(1,4)==1
    # For the varied last block the inverse diagonal>5, so the old >=2
    # bound used in these third-leaf supports is a valid relaxation.
    upper=[]
    for row in data['anchor_contact_upper_witnesses']:
        pair=tuple(row['pair']);assert pair in[(11,14),(14,17)]and row['test_label']==4
        v=V[row['actual_vertex']];assert set(pair)<=set(v)
        A=gram(v,labels|{pair:4});k=row['leading_order']
        d=det_gauss([r[:k]for r in A[:k]])
        p=det_gauss([r[:k-1]for r in A[:k-1]])if k>1 else O
        assert sign(p)>0 and sign(d)<=0 and mul(q2(row['last_pivot']),p)==d
        upper.append(dict(pair=list(pair),actual_vertex=row['actual_vertex'],prefix=k,determinant=enc(d)))
    assert len(upper)==2 and {tuple(r['pair'])for r in upper}=={(11,14),(14,17)}

    # Independent threshold proof using c=cos(2pi/7), rather than the
    # primary U6 recurrence: f(c)=8c^3+4c^2-4c-1 follows from the seven
    # roots of unity. c>1/2, f is increasing there, and f(3/5)<0.
    f=lambda c:8*c**3+4*c*c-4*c-1
    assert f(F(3,5))==F(-29,125)
    assert 24*F(1,2)**2+8*F(1,2)-4==6
    assert 48*F(1,2)+8==32
    assert (1+F(3,5))/2==F(4,5)and 1/(4*(1-F(4,5)))==F(5,4)
    # Polynomial substitution c=2u-1 also independently verifies the
    # primary cosine-square cubic and its recorded rational checks.
    c=(neg(O),scale(O,2));c2=pmul(c,c);c3=pmul(c2,c)
    coeff=padd(padd(tuple(scale(v,8)for v in c3),tuple(scale(v,4)for v in c2)),
               padd(tuple(scale(v,-4)for v in c),(neg(O),)))
    assert coeff==tuple(scale(O,v)for v in[-1,24,-80,64])
    assert data['cosine_square_cubic']==[-1,24,-80,64]
    assert data['chebyshev_second_kind_degree6']==[-1,0,24,0,-80,0,64]
    assert F(data['cubic_at_four_fifths'])==F(-29,125)
    assert data['derivative_at_three_quarters']==12 and data['second_derivative_at_three_quarters']==128

    B=[[O if i==j else neg(CO[4])if i//2==j//2 else Z for j in range(4)]for i in range(4)]
    J=inverse(B)
    def dot(a,b):return total(mul(a[i],mul(J[i][j],b[j]))for i in range(4)for j in range(4))
    patterns={None:[Z]*4}
    for i in range(4):
        p=[Z]*4;p[i]=neg(CO[3]);patterns[i+2]=p
    c=[neg(CO[3]),neg(CO[3]),Z,Z]
    saved={(tuple(r['foreign_contacts']),tuple(r['contact_labels'])):r for r in data['cases']}
    sa={tuple(r['foreign_contacts']):r for r in data['anchors']}
    assert len(saved)==100==len(data['cases'])and len(sa)==25==len(data['anchors'])
    kmin=scale(O,F(5,4));qmin=neg(add(O,CO[4]));records=[];anchor_records=[];values=[]
    for ia,ib in product(patterns,repeat=2):
        a,b=patterns[ia],patterns[ib];u,v,p=dot(a,a),dot(b,b),dot(a,b)
        qa={(0,0):sub(O,u),(1,0):neg(O)};qb={(0,0):sub(O,v),(1,0):neg(O)}
        ab={(0,0):neg(p),(1,0):neg(O)}
        A=mpdet([[qa,ab],[ab,qb]])
        assert all(j==0 and i<=1 for i,j in A)
        a1=A[(1,0)];av=mpeval(A,kmin,Z)
        assert sign(a1)<0 and sign(av)<0
        assert a1==q2(sa[ia,ib]['negative_k_slope'])and av==q2(sa[ia,ib]['determinant_at_five_quarters'])
        anchor_records.append(dict(foreign_contacts=[ia,ib],negative_k_slope=enc(a1),negative_boundary=enc(av)))
        for ma,mb in product([2,3],repeat=2):
            ac=sub(neg(CO[ma]),dot(a,c));bc=sub(neg(CO[mb]),dot(b,c))
            D=mpdet([[qa,ab,{(0,0):ac}],[ab,qb,{(0,0):bc}],
                     [{(0,0):ac},{(0,0):bc},{(0,1):O}]])
            sq=mul(sub(ac,bc),sub(ac,bc))
            assert all(i<=1 and j<=1 for i,j in D)
            assert D.get((1,1),Z)==a1 and D.get((1,0),Z)==sq
            assert D.get((0,1),Z)==A.get((0,0),Z)
            value=mpeval(D,kmin,qmin);slope=add(mul(a1,qmin),sq)
            assert sign(value)>0 and sign(slope)>0 and sign(sq)>=0
            row=saved[(ia,ib),(ma,mb)]
            for key,val in [('negative_anchor_slope',a1),('nonnegative_square',sq),
                            ('positive_k_slope_at_boundary',slope),('positive_minimum',value)]:assert q2(row[key])==val
            values.append(value)
            records.append(dict(foreign_contacts=[ia,ib],contact_labels=[ma,mb],
                                polynomial_in_k_q=[dict(k_degree=i,q_degree=j,coefficient=enc(v))for(i,j),v in sorted(D.items())],
                                positive_boundary=enc(value)))
    minimum=add(scale(O,F(1,8)),scale(CO[4],F(1,4)))
    assert minimum in values and all(sign(sub(v,minimum))>=0 for v in values)
    assert minimum==q2(data['minimum_boundary'])and len(records)==100
    paths=[src,f'scripts/{TAG}.py',gp,sp,
           'scripts/fresh_h9_star_cube_metric_20261014_cut_profiles_review.py',
           'scripts/fresh_h9_star_cube_metric_20261014_all_cut_profiles_review.py',
           'scripts/fresh_h9_t8_metric_20261014_review.py',
           'scripts/fresh_h9_t8_metric_20261014_two_labels_review.py',f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_T8_TWO_HIGH_LABEL_DOMAIN_AND_RANK_REVIEW',
                sources={p:sha(p)for p in paths},recursively_checked_manifests=sorted(bound),
                high_label_cross_guards=n_high,relaxed_affine_cross_guards=n_relaxed,
                mutual_anchor_root=13,changed_anchor_budgets=budgets,third_leaf_budgets=third,
                fresh_upper_failure_minors=upper,
                independent_cosine_threshold={'variable':'c=cos(2pi/7)>1/2','cubic':[-1,-4,4,8],
                    'cubic_at_three_fifths':'-29/125','positive_derivative_bound':'6','conclusion':'cos²(pi/n)>4/5 and k>5/4 for every finite n>=7'},
                anchors=anchor_records,exact_bivariate_cases=records,minimum_boundary=enc(minimum),
                scope='The new saved T8 assignment with both category7 labels finite m,n>=7 and anchors own(2,3). Fresh literal domains plus100 independently rebuilt bivariate determinant identities. k>5/4,t>3/4; determinants strictly increase in both and are positive at the boundary. No numerical or finite-angle cutoff. Other assignments and the full incidence remain outside.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}_review.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['exact_bivariate_cases']),flush=True)
