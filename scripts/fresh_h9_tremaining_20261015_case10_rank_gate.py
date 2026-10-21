#!/usr/bin/env python3
"""Exact metric exclusion of the first T10 mandatory-label family, all m>=7.

The large-label bound is a separately reviewed theorem. This checker verifies
the remaining domains and 100 affine-in-d residual determinants. Chebyshev
identities are used only as necessary equations; no numerical root isolation
or cyclotomic-degree assumption is needed.
"""
import argparse
import json
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import fresh_h20_joint_route_20261015_exact_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_tremaining_20261015_case10_rank_gate'
PRE='fresh_h9_tremaining_20261015_case10_high_bound'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
STATE='results/fresh_h9_tremaining_20261015_search_case_10_batch_001.json'
DOM='results/fresh_h9_tremaining_20261015_case10_domains.json'
HELPER='scripts/fresh_h20_joint_route_20261015_exact_review.py'
Z=R.ZERO;O=R.ONE;H=R.dec('1/2')
SQ=R.dec([0,1,0,0]);s=R.mul(H,SQ);u=R.add(O,s)
COS={2:Z,3:H,4:s,6:R.dec([0,0,'1/2',0])}
S=[3,4,6,7,9,10,12,13]


def enc(x):return [str(v)for v in x]


def determinant(M):
    ans=Z;n=len(M)
    for p in permutations(range(n)):
        t=O
        for i in range(n):t=R.mul(t,M[i][p[i]])
        if sum(p[i]>p[j]for i in range(n)for j in range(i+1,n))%2:t=R.neg(t)
        ans=R.add(ans,t)
    return ans


def chebyshev(n,z):
    a,b=O,z
    for _ in range(n-1):a,b=b,R.sub(R.mul(R.dec(2),R.mul(z,b)),a)
    return b if n else a


def run():
    def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
    seen=set()
    def read(p):
        d=json.loads((ROOT/p).read_text())
        if p not in seen:
            seen.add(p)
            for q,h in d.get('sources',{}).items():
                assert sha(q)==h,q
                if q.endswith('.json'):read(q)
        return d
    geo,state,domains,high,review=map(read,(GEOM,STATE,DOM,
        f'results/{PRE}.json',f'results/{PRE}_review.json'))
    assert high['status']=='PASS_T10_SAVED_FAMILY_HIGH_LABEL_AT_MOST_ELEVEN'
    assert high['remaining_integer_high_labels']==[7,8,9,10,11]
    assert review['status'].startswith('PASS')
    g=geo['representatives'][10]['geometry']
    V=[set(v['incident'])for v in g['vertices']]
    labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
    assert labels=={tuple(r[:2]):r[2]for r in high['fixed_mandatory_family']}
    assert [labels[tuple(b)]for b in [[3,4],[6,7],[9,10],[12,13]]]==[4,4,7,4]
    for leaf,own in [(2,[3,4]),(11,[12,13])]:
        assert all(labels[tuple(sorted((leaf,r)))]==3 for r in own)
    core=S+[11,16,2]
    assert len(set(core))==11
    assert all(any(set(p)<=v for v in V)for p in combinations(core,2))
    rows={tuple(r['pair']):r for r in domains['all_domains']}
    extra=[]
    for pair in [(2,11),(2,16)]:
        row=rows[pair];assert row['domain']==[2,3]
        w=next(r for r in row['failures']if r['label']==4)
        supp=w['support'];assert set(supp)==V[w['actual_vertex']]
        M=[]
        for a in supp:
            line=[]
            for b in supp:
                if a==b:line.append(O);continue
                p=tuple(sorted((a,b)));m=4 if p==pair else labels.get(p,2)
                if m==7:m=6
                line.append(R.neg(COS[m]))
            M.append(line)
        A=[r[:]for r in M];pivots=[]
        while A:
            p=A[0][0];pivots.append(p)
            if R.sign(p)<=0:break
            q=R.inv(p)
            A=[[R.sub(A[i][j],R.mul(R.mul(A[i][0],q),A[0][j]))
                for j in range(1,len(A))]for i in range(1,len(A))]
        assert R.sign(pivots[-1])<=0
        extra.append(dict(pair=list(pair),trial_label=4,actual_vertex=w['actual_vertex'],
                          support=supp,relaxed_high_label=6,pivots=list(map(enc,pivots))))
    budgets=[]
    for leaf,own,foreign in [(11,[12,13],[3,4,6,7]),(2,[3,4],[6,7,12,13])]:
        for omit in own:
            support=sorted(set(S)-{omit}|{leaf})
            vi=next(i for i,v in enumerate(V)if set(support)<=v)
            budgets.append(dict(leaf=leaf,removed_own_root=omit,support=support,actual_vertex=vi,
                                retained_own_cost='1/4',foreign_roots=foreign))
    # Each foreign low-B2 root costs at least 2*(1/2)^2=1/2;
    # nonnegative inverse entries make these additive lower bounds.
    # Two such contacts, or one label>=4, contradict the actual budget.
    # Each high-C contact costs >=4*(1/2)^2=1 already at m=6,
    # again contradicting the extra retained-own cost 1/4.
    assert R.sign(R.sub(R.add(R.dec('1/4'),R.mul(R.dec(2),R.dec('1/2'))),O))>0
    assert R.sign(R.sub(R.add(R.dec('1/4'),O),O))>0
    J=[[Z for _ in S]for _ in S]
    for a in [0,2,6]:
        J[a][a]=J[a+1][a+1]=R.dec(2)
        J[a][a+1]=J[a+1][a]=SQ
        G=[[O,R.neg(s)],[R.neg(s),O]]
        assert [[R.dot(row,col)for col in zip(*[[R.dec(2),SQ],[SQ,R.dec(2)]])]
                for row in G]==[[O,Z],[Z,O]]
    def pats(own,foreign):
        out=[]
        for hit in [None]+foreign:
            p=[Z]*8
            for r in own:p[S.index(r)]=H
            if hit is not None:p[S.index(hit)]=H
            out.append((hit,p))
        return out
    P=pats([12,13],[3,4,6,7]);T=pats([3,4],[6,7,12,13])
    v=pats([12,13],[])[0][1]
    cost=lambda x,y:R.dot(x,R.mv(J,y))
    assert cost(v,v)==u
    census=[];survivors=[]
    for (ah,a),(th,t),l1,l2 in product(P,T,[2,3],[2,3]):
        aa=R.sub(O,cost(a,a));tt=R.sub(O,cost(t,t))
        at=R.sub(R.neg(COS[l1]),cost(a,t))
        av=R.neg(cost(a,v));vt=R.sub(R.neg(COS[l2]),cost(v,t))
        def matrix(d):return [[aa,av,at],[av,R.sub(R.sub(O,u),d),vt],[at,vt,tt]]
        B=determinant(matrix(Z));A=R.sub(determinant(matrix(O)),B)
        assert determinant(matrix(R.dec(2)))==R.add(B,R.mul(R.dec(2),A))
        assert A!=Z
        d=R.mul(R.neg(B),R.inv(A));assert d!=Z
        z=R.sub(O,R.inv(R.mul(R.dec(2),d)))
        tests={str(m):chebyshev(m,z)for m in [7,8,9,10,11]}
        eligible=[int(m)for m,value in tests.items()if value==O]
        assert not eligible or eligible==[8]
        record=dict(anchor11_foreign=ah,third2_foreign=th,pair_11_2=l1,pair_16_2=l2,
                    coefficient_A=enc(A),coefficient_B=enc(B),unique_d_root=enc(d),
                    cosine2theta=enc(z),chebyshev_values={m:enc(x)for m,x in tests.items()},
                    necessary_integer_labels=eligible)
        if eligible:
            assert d==u and determinant(matrix(u))==Z
            minor=R.sub(R.mul(aa,tt),R.mul(at,at))
            assert R.sign(aa)<0 and R.sign(tt)<0 and R.sign(minor)>0
            record['negative_plane_11_2']=dict(norm11=enc(aa),norm2=enc(tt),inner=enc(at),determinant=enc(minor))
            survivors.append([ah,th,l1,l2])
        census.append(record)
    assert len(census)==100 and len(survivors)==4
    assert survivors==[[6,None,2,2],[6,7,2,3],[7,None,2,2],[7,6,2,3]]
    paths=[GEOM,STATE,DOM,HELPER,f'scripts/{PRE}.py',f'results/{PRE}.json',
           f'scripts/{PRE}_review.py',f'results/{PRE}_review.json',f'scripts/{TAG}.py']
    return dict(status='PASS_T10_FIRST_MANDATORY_FAMILY_ALL_HIGH_LABELS_EXCLUDED',
                sources={p:sha(p)for p in paths},case=10,fixed_mandatory_family=state['solution']['pair_labels'],
                positive_support=S,residual_normals=[11,16,2],core_meeting_pairs=55,
                inherited_forced_right_guards=31,new_upper_guards=extra,budget_supports=budgets,
                projection_patterns=dict(anchor11=5,third2=5,anchor16=1),
                high_parameter='d=1/(4sin²(pi/m)), m>=7',
                determinant_identity='det(residual)=A*d+B, A!=0 in all100 cases',
                necessary_angle_identity='T_m(1-1/(2*d))=1',cases=census,
                counts=dict(cases=100,nonzero_slopes=100,chebyshev_tests=500,
                            angle_survivors=4,negative_plane_closures=4),
                scope='All omitted meeting labels and all dotted distances remain free. '
                      'The separately reviewed negative-plane theorem excludes every m>=12. '
                      'For m7..11, positive-S8 Schur rank at most2 forces the displayed affine '
                      'determinant to vanish. The necessary Chebyshev equation leaves four m8 '
                      'cases; each has a negative definite residual2-plane. This excludes the '
                      'saved mandatory-label family in compact H9, not the entire T10 incidence. '
                      'The rank step is dimension9-specific; the final negative-plane step is not.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'])
