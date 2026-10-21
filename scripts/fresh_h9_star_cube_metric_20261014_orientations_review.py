#!/usr/bin/env python3
"""Separate Fraction review of the sixteen fixed-B2 leaf orientations.

No producer or SymPy imports.  Reuses frozen independent square-root-basis
arithmetic, but rebuilds the geometric domains and every distinct residual.
"""
import argparse
import json
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
from fresh_h9_star_cube_metric_20261014_cut_profiles_review import (
    F,Z,O,CO,add,sub,neg,mul,scale,total,inverse,determinant,anp,enc)
from fresh_h9_star_cube_metric_20261014_all_cut_profiles_review import sign,det_gauss

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_cube_metric_20261014_orientations'
def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    src=f'results/{TAG}.json'
    data=json.loads((ROOT/src).read_text())
    assert data['status']=='PASS_ALL_SIXTEEN_ORIENTATIONS_EXCLUDED'
    manifests=set()
    def bind_manifest(path):
        if path in manifests: return
        manifests.add(path)
        obj=json.loads((ROOT/path).read_text())
        for p,h in obj.get('sources',{}).items():
            assert sha(p)==h,p
            if p.endswith('.json'): bind_manifest(p)
    bind_manifest(src)
    gp='results/fresh_h9_tree_cube_wedges_20261014.json'
    g=json.loads((ROOT/gp).read_text())['representatives'][0]['geometry']
    V=[v['incident'] for v in g['vertices']]
    M={tuple(m) for m in g['minimal_nonfaces']}
    S=[3,4,6,7,9,10,12,13]; P=[2,5,8,11]; C=[14,15,16,17]
    first=json.loads((ROOT/'results/fresh_h9_star_cube_metric_20261014.json').read_text())
    assert (S,P,C)==(first['fixed_positive_support'],first['positive_leaves'],first['cut_leaves'])
    assert all(tuple(sorted((leaf,*S[2*i:2*i+2]))) in M
               for i in range(4) for leaf in [P[i],C[i]])
    base={(S[2*i],S[2*i+1]):4 for i in range(4)}
    def gram(nodes,labels):
        return [[O if a==b else neg(CO[labels.get(tuple(sorted((a,b))),2)])
                 for b in nodes] for a in nodes]
    @lru_cache(None)
    def prefix_dets(A,k):
        d=det_gauss([r[:k] for r in A[:k]])
        p=det_gauss([r[:k-1] for r in A[:k-1]]) if k>1 else O
        return d,p
    failures=[]
    def failure(row,labels,pair_gate=True):
        v=V[row['actual_vertex']]
        lab=labels.copy()
        if pair_gate:
            pair=tuple(row['pair']); assert set(pair)<=set(v)
            lab[pair]=row['test_label']
        A=tuple(map(tuple,gram(v,lab))); k=row['leading_order']
        d,p=prefix_dets(A,k)
        assert sign(p)>0 and sign(d)<=0
        assert mul(anp(row['last_pivot']),p)==d
        failures.append(dict(actual_vertex=row['actual_vertex'],prefix=k,determinant=enc(d)))

    # This preliminary gate has no leaf-contact premise at all.
    cross={tuple(sorted((a,b))) for i,j in combinations(range(4),2)
           for a in S[2*i:2*i+2] for b in S[2*j:2*j+2]}
    assert {tuple(r['pair']) for r in data['forced_cross_right_witnesses']}==cross
    assert len(data['forced_cross_right_witnesses'])==24
    for row in data['forced_cross_right_witnesses']:
        assert row['test_label']==3
        failure(row,base)
        a,b=row['pair']; i=S.index(a)//2; j=S.index(b)//2
        aa=next(x for x in S[2*i:2*i+2] if x!=a)
        bb=next(x for x in S[2*j:2*j+2] if x!=b)
        chain=[aa,a,b,bb]
        assert any(set(chain)<=set(v) for v in V)
        A=gram(chain,base|{tuple(sorted((a,b))):3})
        root2=scale(CO[4],2); null=[O,root2,root2,O]
        assert all(total(mul(x,y) for x,y in zip(r,null))==Z for r in A)

    profiles=[pr for pr in product(range(2,6),repeat=2)
              if F(1,4)+sum(F(1,m) for m in pr)<1]
    assert len(profiles)==11
    assert set(profiles)==set(product([3,4,5],repeat=2))|{(2,5),(5,2)}
    B=gram(S,base); J=inverse(B)
    assert J==[[(scale(O,2) if i==j else scale(CO[4],2)) if i//2==j//2 else Z
                for j in range(8)] for i in range(8)]
    budgets=[]; triples={}
    for i in range(4):
        for leaf in [P[i],C[i]]:
            for j in range(2):
                keep=S[2*i+j]; nodes=[x for x in S if x!=S[2*i+1-j]]
                support=sorted(nodes+[leaf]); vi=next(k for k,v in enumerate(V) if set(support)<=set(v))
                K=inverse(gram(nodes,base)); own=nodes.index(keep)
                assert K[own][own]==O and all(sign(x)>=0 for r in K for x in r)
                assert all(K[k][k]==scale(O,2) for k,x in enumerate(nodes) if x not in S[2*i:2*i+2])
                budgets.append(dict(leaf=leaf,retained_own_root=keep,support=support,actual_vertex=vi))
        for j in range(2):
            tri=sorted([P[i],C[i],S[2*i+j]])
            triples[i,j]=next(k for k,v in enumerate(V) if set(tri)<=set(v))
    assert budgets==data['all_leaf_both_orientation_budget_guards']
    assert len(budgets)==16
    for row in data['base_own_contact_guards']:
        i=C.index(row['cut']); j=S.index(row['increased_own_root'])%2
        assert row['path_actual_vertex']==triples[i,j]
        assert set(row['label6_path'])=={P[i],C[i],S[2*i+j]}
        assert set(row['eight_support'])==(set(S)-{S[2*i+1-j]})|{C[i]}
        assert set(row['eight_support'])<=set(V[row['actual_vertex']])
    assert len(data['base_own_contact_guards'])==8
    # Generic own-pair gate precedes this bound: 3-2*(1/2+1/2+c)<=0
    # whenever the free missing triangle has an own c>=1/2.
    assert 3-6*F(1,2)==0
    assert 1-F(1,4)-F(3,4)==0  # actual 3--6 path
    assert F(1,4)+2*F(1,4)*2>1  # two foreign contacts
    assert F(1,4)+2*F(1,2)>1    # one foreign label>=4
    assert F(1,2)+2*F(1,4)==1   # non33 profile: no foreign contact

    @lru_cache(None)
    def pdot(x,y):
        return total(mul(x[i],mul(J[i][j],y[j])) for i in range(8) for j in range(8)
                     if x[i]!=Z and J[i][j]!=Z and y[j]!=Z)
    @lru_cache(None)
    def patterns(group,pr):
        out=[]
        for extra in ([None]+[k for k in range(8) if k//2!=group] if pr==(3,3) else [None]):
            p=[Z]*8; p[2*group]=neg(CO[pr[0]]); p[2*group+1]=neg(CO[pr[1]])
            if extra is not None: p[extra]=neg(CO[3])
            out.append((extra,tuple(p)))
        return tuple(out)
    dets={}; anchor_dets={}; occurrences=0; zeros=[]; uppers=0; summaries=[]; own_checked=0
    assert [o['orientation'] for o in data['orientations']]==list(range(16))
    for ori in data['orientations']:
        o=ori['orientation']; fixed=[C[i] if o&(1<<i) else P[i] for i in range(4)]
        assert fixed==ori['fixed_33_leaves']
        assert ori['status']=='EXACT_ORIENTATION_EXCLUDED' and not ori['missing_label4_guards']
        lab0=base.copy()
        for i,leaf in enumerate(fixed):
            for root in S[2*i:2*i+2]: lab0[tuple(sorted((leaf,root)))]=3
            free=C[i] if leaf==P[i] else P[i]
            for pr in profiles:
                j=next(j for j,m in enumerate(pr) if m>=3)
                tri=[P[i],C[i],S[2*i+j]]
                assert set(tri)<=set(V[triples[i,j]])
                lab=lab0|{tuple(sorted((P[i],C[i]))):3,
                          tuple(sorted((free,S[2*i+j]))):pr[j]}
                assert sign(total(x for r in gram(tri,lab) for x in r))<=0
                own_checked+=1
        thirds=[(3,3)] if o&1 else profiles
        expected_rows={}
        for third,free in product(thirds,profiles):
            ap,bp=(free,(3,3)) if o&8 else ((3,3),free)
            expected_rows[third,ap,bp]=None
        got_rows={}; local_n=0; local_zero=[]
        for row in ori['records']:
            third,ap,bp=map(tuple,[row['C0_profile'],row['P3_profile'],row['C3_profile']])
            keyrow=third,ap,bp
            assert keyrow in expected_rows and keyrow not in got_rows
            got_rows[keyrow]=True
            labels=lab0.copy()
            for leaf,group,pr in [(C[0],0,third),(P[3],3,ap),(C[3],3,bp)]:
                for j,m in enumerate(pr): labels[tuple(sorted((leaf,S[2*group+j])))]=m
            assert {tuple(r['pair']) for r in row['upper_witnesses']}=={(11,14),(14,17)}
            for r in row['upper_witnesses']:
                assert r['test_label']==4; failure(r,labels); uppers+=1
            expected={}
            for ia,x in patterns(3,ap):
                for ib,y in patterns(3,bp):
                    qx=sub(O,pdot(x,x)); qy=sub(O,pdot(y,y)); xy=neg(pdot(x,y))
                    ak=ap,bp,ia,ib
                    if ak not in anchor_dets:
                        anchor_dets[ak]=determinant([[qx,xy],[xy,qy]])
                    assert anchor_dets[ak]!=Z
                    for ic,z in patterns(0,third):
                        qz=sub(O,pdot(z,z))
                        for ma,mb in product([2,3],repeat=2):
                            dk=(ap,bp,third,ia,ib,ic,ma,mb)
                            if dk not in dets:
                                xz=sub(neg(CO[ma]),pdot(x,z)); yz=sub(neg(CO[mb]),pdot(y,z))
                                dets[dk]=determinant([[qx,xy,xz],[xy,qy,yz],[xz,yz,qz]])
                            expected[(ia,ib,ic),(ma,mb)]=dets[dk]
            saved={}
            for case in row['cases']:
                k=tuple(case['foreign_contacts']),tuple(case['anchor_labels'])
                assert k not in saved; saved[k]=anp(case['residual_determinant'])
            assert saved==expected
            local_n+=len(expected)
            for k,d in expected.items():
                if d==Z: local_zero.append((third,ap,bp,*k))
        assert set(got_rows)==set(expected_rows)
        assert local_n==ori['residual_cases'] and len(local_zero)==4
        saved_zero=[]
        for row in ori['zero_rank_cases']:
            key=tuple(row['C0_profile']),tuple(row['P3_profile']),tuple(row['C3_profile']),tuple(row['foreign_contacts']),tuple(row['anchor_labels'])
            assert key in local_zero and key not in saved_zero; saved_zero.append(key)
            third,ap,bp,(ia,ib,ic),(ma,mb)=key
            labels=lab0.copy()
            for leaf,group,pr,extra in [(C[0],0,third,ic),(P[3],3,ap,ia),(C[3],3,bp,ib)]:
                for j,m in enumerate(pr): labels[tuple(sorted((leaf,S[2*group+j])))]=m
                if extra is not None: labels[tuple(sorted((leaf,S[extra])))]=3
            labels[11,14]=ma; labels[14,17]=mb
            witness=row['vertex_exclusion']; assert witness is not None
            failure(witness,labels,False)
            # Independently locate a positive affine null minor where one is
            # present.  The full-vertex determinant gate covers every zero.
            v=V[witness['actual_vertex']]; cycles=[]
            for nodes in combinations(v,4):
                H=gram(nodes,labels)
                if all(total(r)==Z for r in H): cycles.append(list(nodes))
            zeros.append(dict(orientation=o,C0_profile=list(third),P3_profile=list(ap),C3_profile=list(bp),
                              foreign_contacts=list(key[3]),anchor_labels=list(key[4]),
                              actual_vertex=witness['actual_vertex'],affine_four_cycles=cycles))
        assert set(saved_zero)==set(local_zero)
        occurrences+=local_n
        summaries.append(dict(orientation=o,records=len(got_rows),residual_occurrences=local_n,zeros=len(local_zero)))
        print('review orientation',o,'occurrences',local_n,'zeros',len(local_zero),flush=True)
    assert occurrences==91392 and len(dets)==12852==data['unique_determinants_computed']
    assert uppers==2112 and len(zeros)==64 and own_checked==704
    paths=[src,f'scripts/{TAG}.py',gp,
           'scripts/fresh_h9_star_cube_metric_20261014_cut_profiles_review.py',
           'scripts/fresh_h9_star_cube_metric_20261014_all_cut_profiles_review.py',
           f'scripts/{TAG}_review.py']
    return dict(status='PASS_SEPARATE_SIXTEEN_ORIENTATION_DOMAIN_AND_RANK_REVIEW',
                sources={p:sha(p) for p in paths},recursively_checked_manifests=sorted(manifests),
                exact_field='Fraction basis 1,sqrt2,sqrt5,sqrt10; rational interval signs',
                cross_block_guards=24,all_leaf_budget_guards=16,own_pair_checks=own_checked,
                label6_paths=8,fresh_anchor_upper_guards=uppers,
                orientation_summaries=summaries,residual_occurrences=occurrences,
                independent_distinct_residual_cases=len(dets),distinct_nondegenerate_anchor_cases=len(anchor_dets),
                rank_nonzero_occurrences=occurrences-len(zeros),zero_cases_closed=zeros,
                failure_minor_checks=len(failures),unique_prefix_minor_calculations=prefix_dets.cache_info().currsize,
                scope='Fixed four B2 blocks; in each original positive/cut leaf pair at least one leaf has own contacts (3,3). The partner has any allowable Coxeter contacts. Rebuilds every literal orientation domain, all12852 distinct residual determinants and compares every91392 saved occurrence; checks64 actual-vertex failure minors. No incidence automorphism is assumed. Compact H9 realizations in this family are excluded, not arbitrary changed B2/leaf profiles or the full incidence.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}_review.json'
    if a.verify: assert d==json.loads(out.read_text())
    else: out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['residual_occurrences'],len(d['zero_cases_closed']),flush=True)
