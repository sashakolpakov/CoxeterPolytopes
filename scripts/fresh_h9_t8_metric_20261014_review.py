#!/usr/bin/env python3
"""Separate Fraction/polynomial verification of the continuous T8 gate.

No imports from the primary producer or its arithmetic implementation.
Uses the frozen independent square-root-basis engine, with only Q(sqrt2)
coefficients needed here. The high-label parameter is kept symbolic.
"""
import argparse,json
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
from fresh_h9_star_cube_metric_20261014_cut_profiles_review import (
    F,Z,O,CO,add,sub,neg,mul,scale,total,inverse,determinant,enc)
from fresh_h9_star_cube_metric_20261014_all_cut_profiles_review import sign,det_gauss
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_metric_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def q2(x):assert len(x)==2;return (F(x[0]),F(x[1]),F(0),F(0))
def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==Z:p.pop()
    return tuple(p)
def padd(a,b):return trim([add(a[i] if i<len(a) else Z,b[i] if i<len(b) else Z)for i in range(max(len(a),len(b)))])
def pneg(a):return tuple(neg(x)for x in a)
def pmul(a,b):
    out=[Z]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=add(out[i+j],mul(x,y))
    return trim(out)
def pdet(A):
    n=len(A);out=(Z,)
    for p in permutations(range(n)):
        val=(O,)
        for i,j in enumerate(p):val=pmul(val,A[i][j])
        if sum(p[i]>p[j]for i in range(n)for j in range(i+1,n))%2:val=pneg(val)
        out=padd(out,val)
    return out

def run():
    src=f'results/{TAG}.json';data=json.loads((ROOT/src).read_text())
    assert data['status']=='PASS_T8_SAVED_ASSIGNMENT_ALL_M_GE_7_RANK_ELEVEN'
    manifests=set()
    def bind(p):
        if p in manifests:return
        manifests.add(p)
        for s,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():
            assert sha(s)==h,s
            if s.endswith('.json'):bind(s)
    bind(src)
    gp='results/fresh_h9_tree_cube_wedges_20261014.json'
    sp='results/fresh_h9_tree_cube_wedges_20261014_case_08.json'
    g=json.loads((ROOT/gp).read_text())['representatives'][8]['geometry']
    state=json.loads((ROOT/sp).read_text())
    assert state['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
    labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
    assert labels[1,2]==7 and [p for p,m in labels.items()if m==7]==[(1,2)]
    V=[v['incident']for v in g['vertices']]
    S=[1,2,4,5,7,8,12,13];blocks=[S[2*i:2*i+2]for i in range(4)]
    assert data['positive_support']==S and set(S)<=set(V[data['positive_support_actual_vertex']])
    assert [labels[tuple(p)]for p in blocks]==[7,4,4,4]
    # Polynomial variable u=cos(pi/m); t=u^2 is introduced only afterward.
    zero,one,u=(Z,),(O,),(Z,O)
    I2=[[one,pneg(u)],[pneg(u),one]]
    detI2=pdet(I2);assert detI2==(O,Z,neg(O))
    adj=[[one,u],[u,one]]
    for i in range(2):
        for j in range(2):
            val=zero
            for k in range(2):val=padd(val,pmul(I2[i][k],adj[k][j]))
            assert val==(detI2 if i==j else zero)
    b2=[[O,neg(CO[4])],[neg(CO[4]),O]]
    assert determinant(b2)==scale(O,F(1,2))
    # Three B2 factors give det B=(1-u^2)/8. All blocks are PD for 0<u<1.
    detB=tuple(scale(v,F(1,8))for v in detI2)
    assert detB==(scale(O,F(1,8)),Z,scale(O,F(-1,8)))
    cross_expected={tuple(sorted((a,b)))for i,j in combinations(range(4),2)for a,b in product(blocks[i],blocks[j])}
    seen=set();high=0;affine=0
    for row in data['forced_cross_right_witnesses']:
        a,b=row['pair'];assert(a,b)in cross_expected and(a,b)not in seen;seen.add((a,b))
        support=row['support'];assert set(support)<=set(V[row['actual_vertex']])
        i=S.index(a)//2;j=S.index(b)//2
        if i==0:
            assert support==[next(x for x in blocks[0]if x!=a),a,b]
            H=[[one,pneg(u),zero],[pneg(u),one,(neg(CO[3]),)],[zero,(neg(CO[3]),),one]]
            assert pdet(H)==(scale(O,F(3,4)),Z,neg(O))
            high+=1
        else:
            assert support==[next(x for x in blocks[i]if x!=a),a,b,next(x for x in blocks[j]if x!=b)]
            H=[[O if x==y else Z for y in range(4)]for x in range(4)]
            for k,c in enumerate([CO[4],CO[3],CO[4]]):H[k][k+1]=H[k+1][k]=neg(c)
            null=[O,scale(CO[4],2),scale(CO[4],2),O]
            assert all(total(mul(x,y)for x,y in zip(r,null))==Z for r in H)
            assert H==[[q2(x)for x in r]for r in row['gram']]
            assert null==list(map(q2,row['null_vector']))
            affine+=1
    assert seen==cross_expected and high==affine==12
    own=data['anchor_pair_right_witness'];assert own['pair']==[11,17]and own['test_label']==3
    assert set(own['support'])=={11,12,17}<=set(V[own['actual_vertex']])
    assert labels[11,12]==labels[12,17]==3
    A=[[O if i==j else neg(CO[3])for j in range(3)]for i in range(3)]
    assert A==[[q2(x)for x in r]for r in own['gram']] and all(total(r)==Z for r in A)

    def gram(nodes,lab):
        assert not {1,2}<=set(nodes)
        return [[O if a==b else neg(CO[lab.get(tuple(sorted((a,b))),2)])for b in nodes]for a in nodes]
    base={tuple(p):4 for p in blocks[1:]}
    seen_budget=set()
    for row in data['anchor_budget_witnesses']:
        leaf,kept=row['leaf'],row['retained_own_root']
        assert leaf in[11,17]and kept in[12,13]
        assert labels[tuple(sorted((leaf,kept)))]==3
        assert set(row['support'])==(set(S)-({12,13}-{kept}))|{leaf}
        assert set(row['support'])<=set(V[row['actual_vertex']])
        # I2 inverse is entrywise positive and has diagonal 1/(1-t)>4.
        # A contact of magnitude >=1/2 therefore already costs >1.
        nodes=[4,5,7,8,kept];J=inverse(gram(nodes,base))
        assert J[-1][-1]==O and all(sign(x)>=0 for r in J for x in r)
        assert all(J[k][k]==scale(O,2)for k in range(4))
        seen_budget.add((leaf,kept))
    assert seen_budget==set(product([11,17],[12,13]))and len(data['anchor_budget_witnesses'])==4
    assert F(1,4)/(1-F(3,4))==1
    assert F(1,4)+2*F(1,4)*2>1 and F(1,4)+2*F(1,2)>1
    seen_third=set()
    for row in data['third_leaf_budget_witnesses']:
        assert row['leaf']==14
        a,kept=row['retained_own_roots'];assert a==2 and kept in[4,5]
        assert set(row['support'])=={2,kept,7,8,12,13,14}<=set(V[row['actual_vertex']])
        nodes=[2,kept,7,8,12,13];J=inverse(gram(nodes,base))
        assert J[0][0]==J[1][1]==O and all(sign(x)>=0 for r in J for x in r)
        assert all(J[k][k]==scale(O,2)for k in range(2,6))
        assert labels[2,14]==labels[tuple(sorted((kept,14)))]==3
        seen_third.add(kept)
    assert seen_third=={4,5}and len(data['third_leaf_budget_witnesses'])==2
    assert F(1,4)+F(1,4)+2*F(1,4)==1
    assert labels[1,14]==2
    upper=[]
    for row in data['anchor_contact_upper_witnesses']:
        pair=tuple(row['pair']);assert pair in[(11,14),(14,17)]and row['test_label']==4
        v=V[row['actual_vertex']];assert set(pair)<=set(v)
        H=gram(v,labels|{pair:4});k=row['leading_order']
        d=det_gauss([r[:k]for r in H[:k]])
        p=det_gauss([r[:k-1]for r in H[:k-1]])if k>1 else O
        assert sign(p)>0 and sign(d)<=0 and mul(q2(row['last_pivot']),p)==d
        upper.append(dict(pair=list(pair),actual_vertex=row['actual_vertex'],prefix=k,determinant=enc(d)))
    assert {tuple(r['pair'])for r in upper}=={(11,14),(14,17)}and len(upper)==2

    # Independent 6x6 inverse; the I2 block is handled by its symbolic adjugate.
    nodes=[4,5,7,8,12,13];J=inverse(gram(nodes,base))
    def dot(x,y):return total(mul(x[i],mul(J[i][j],y[j]))for i in range(6)for j in range(6))
    patterns={}
    for extra in[None,2,3,4,5]:
        x=[Z]*6;x[4]=x[5]=neg(CO[3])
        if extra is not None:x[extra-2]=neg(CO[3])
        patterns[extra]=x
    third=[neg(CO[3]),neg(CO[3]),Z,Z,Z,Z]
    qthird=sub(O,dot(third,third));assert qthird==neg(CO[4])
    primary_anchors={tuple(r['foreign_contacts']):r for r in data['anchors']}
    primary_cases={(tuple(r['foreign_contacts']),tuple(r['contact_labels'])):r for r in data['residual_cases']}
    assert len(primary_anchors)==25==len(data['anchors']) and len(primary_cases)==100==len(data['residual_cases'])
    records=[];anchor_values=set();boundaries=set()
    for ia,ib in product(patterns,repeat=2):
        x,y=patterns[ia],patterns[ib]
        qa=sub(O,dot(x,x));qb=sub(O,dot(y,y));ab=neg(dot(x,y))
        small=[[qa,ab],[ab,qb]];A=determinant(small)
        saved=primary_anchors[ia,ib]
        assert small==[[q2(v)for v in r]for r in saved['residual_gram']]
        assert A==q2(saved['determinant'])and sign(A)<0;anchor_values.add(A)
        for ma,mb in product([2,3],repeat=2):
            ac=sub(neg(CO[ma]),dot(x,third));bc=sub(neg(CO[mb]),dot(y,third))
            # Polynomial variable h=1/[4(1-t)].  This is an independent
            # permutation determinant, so linearity/coefficient signs are checked.
            H=[[(qa,),(ab,),(ac,)],[(ab,),(qb,),(bc,)],[(ac,),(bc,),(qthird,neg(O))]]
            coeff=pdet(H);assert len(coeff)==2 and coeff[1]==neg(A)
            D0=coeff[0];positive=scale(neg(A),F(1,4));boundary=sub(D0,A)
            assert sign(positive)>0 and sign(boundary)>0
            saved=primary_cases[(ia,ib),(ma,mb)]
            assert D0==q2(saved['constant'])
            assert positive==q2(saved['positive_inverse_denominator_coefficient'])
            assert boundary==q2(saved['boundary_at_three_quarters'])
            # D(t)=D0+positive/(1-t), t in(3/4,1): denominator positive,
            # derivative positive/(1-t)^2>0 and D(t)>D(3/4)>0.
            assert add(D0,scale(positive,4))==boundary
            boundaries.add(boundary)
            records.append(dict(foreign_contacts=[ia,ib],contact_labels=[ma,mb],
                                anchor_determinant=enc(A),determinant_in_h=list(map(enc,coeff)),
                                boundary_at_three_quarters=enc(boundary)))
    minimum=add(scale(O,F(3,4)),scale(CO[4],F(5,4)))
    assert minimum in boundaries and all(sign(sub(v,minimum))>=0 for v in boundaries)
    assert minimum==q2(data['minimum_boundary'])
    paths=[src,f'scripts/{TAG}.py',gp,sp,
           'scripts/fresh_h9_star_cube_metric_20261014_cut_profiles_review.py',
           'scripts/fresh_h9_star_cube_metric_20261014_all_cut_profiles_review.py',f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_T8_CONTINUOUS_DOMAIN_AND_SCHUR_REVIEW',
                sources={p:sha(p)for p in paths},recursively_checked_manifests=sorted(manifests),
                high_label_cross_guards=high,affine_cross_guards=affine,anchor_mutual_right_guards=1,
                anchor_budget_guards=4,third_budget_guards=2,upper_failure_minors=upper,
                parameter='u=cos(pi/m), t=u^2; 3/4<t<1 for every finite integer m>=7',
                symbolic_high_path_determinant_in_u=list(map(enc,[scale(O,F(3,4)),Z,neg(O)])),
                symbolic_positive_support_determinant_in_u=list(map(enc,detB)),
                exact_anchor_patterns=25,distinct_anchor_determinants=len(anchor_values),
                exact_residual_cases=records,distinct_boundary_values=len(boundaries),minimum_boundary=enc(minimum),
                conclusion='All100 residual determinants are strictly positive throughout3/4<t<1, hence the PD8 Schur decomposition forces rank11. No finite m cutoff and no rank-zero cases.',
                scope='Only the saved mandatory assignment of T-shaped incidence case8. All omitted contacts used in the proof have literal actual-vertex and Schur-budget guards. Other label assignments and the entire T incidence are not excluded. Reuses independent Fraction arithmetic, with no primary arithmetic imports or numerical acceptance predicates.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}_review.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['exact_residual_cases']),flush=True)
