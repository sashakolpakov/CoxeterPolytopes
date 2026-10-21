#!/usr/bin/env python3
"""Saved two-duty checker; imports only the frozen stand-alone interval checker.

Reconstructs all physical projections, both sheets and original positive
blockers. It imports no producer module and does not rerun a cone census.
Original normals and their Coxeter-native premise remain inherited inputs.
"""
import json
from hashlib import sha256
from itertools import combinations
import check_exact66_quadrilateral_cap_saved_20261018 as Q

ROOT=Q.ROOT;TAG='fresh_exact66_old7_boundary_cap_obstruction_20261018'
V,add,sub,mul,div,neg,sqrt,dot=Q.iv,Q.add,Q.sub,Q.mul,Q.div,Q.neg,Q.sqrt,Q.dot

def refine(raw,label):
    if label=='norm2':
        assert raw[0]<=2<=raw[1];return V(2)
    if label==2:
        assert -1<raw[0]<=0<=raw[1]<1;return V(0)
    target=neg(Q.theta(label))
    assert raw[0]>-Q.theta(label+1)[0]and raw[1]<-Q.theta(label-1)[1]
    out=max(raw[0],target[0]),min(raw[1],target[1]);assert out[0]<=out[1]
    return out

def main():
    assert sha256((ROOT/'scripts/check_exact66_quadrilateral_cap_saved_20261018.py').read_bytes()).hexdigest()=='6d84aea1130bdd8afcca3f325204cbdc6ca321ba9cb3015b87d20fa632809b75'
    d=Q.read(TAG);seed=Q.read(Q.SEED);current=Q.read(Q.CURRENT)
    assert d['status']=='PASS_EXACT_OLD7_TWO_FINITE_DUTY_CAP_OBSTRUCTION'
    assert len(d['sources'])==7
    for doc in(d,seed,current):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert current['physical']['summary']['total_walls']==66
    assert current['physical']['summary']['complete_finite_vertices']==205
    p=seed['physical_input'];metric=[list(map(V,r))for r in p['metric']]
    roots={int(i):list(map(V,v))for i,v in p['roots'].items()}
    hc={int(i):V(x)for i,x in p['center_products'].items()}
    white=d['fixed_white_walls'];assert white==[0,1,2,3,5,6]
    frame=d['frame'];assert frame==[7,30,52]
    labels={}
    for r in seed['native_isolation']+d['native_isolation']:
        key=tuple(sorted(r['pair']));old=labels.setdefault(key,r['identity']);assert old==r['identity']
    def native(i,j):
        raw=dot(roots[i],metric,roots[j]);label='norm2'if i==j else labels.get(tuple(sorted((i,j))))
        return raw if label is None else refine(raw,label)
    first=white+[7,30];expected={(0,1):5,(1,2):3,(2,3):3,(5,6):3,(6,7):3,(7,30):5}
    for i,j in combinations(first,2):
        label=expected.get((i,j),2);assert labels[i,j]==label
        refine(dot(roots[i],metric,roots[j]),label)
    assert labels[30,52]==4
    # The two connected components are exactly H4; the remaining fixed
    # label4 three-root path is finite only for K52=2 or3.
    assert d['forced_H4_factors']==[[0,1,2,3],[5,6,7,30]]
    assert d['forced_profile']==0 and d['forced_contacts']=={'7':2,'30':2}
    assert d['complete_remaining_labels52']==[2,3]
    wi,wp=Q.gauss([[native(i,j)for j in white]for i in white]);assert all(x[0]>0 for x in wp)
    ids=frame+[4,47];wc={i:[native(w,i)for w in white]for i in ids}
    hw=[hc[w]for w in white];hp={i:sub(hc[i],dot(hw,wi,wc[i]))for i in ids}
    hp2=sub(V(p['center_square']),dot(hw,wi,hw));assert hp2[1]<0
    g=[[sub(native(i,j),dot(wc[i],wi,wc[j]))for j in frame]for i in frame]
    inv,piv=Q.gauss(g);assert [1 if z[0]>0 else -1 for z in piv]==[1,1,-1]
    adj,det=Q.inverse(g);assert det[1]<0
    for i in range(3):
        for j in range(3):
            inv[i][j]=max(inv[i][j][0],adj[i][j][0]),min(inv[i][j][1],adj[i][j][1])
            assert inv[i][j][0]<=inv[i][j][1]
    hf=[hp[i]for i in frame];q=sub(hp2,dot(hf,inv,hf));assert q[0]>0
    cov={j:[sub(native(i,j),dot(wc[i],wi,wc[j]))for i in frame]for j in(4,47)}
    rows=d['branches'];assert len(rows)==4
    assert {(r['label52'],r['sheet_relative_to_center_residual'])for r in rows}=={(m,s)for m in(2,3)for s in(-1,1)}
    for r in rows:
        m,s=r['label52'],r['sheet_relative_to_center_residual'];rhs=[V(0),V(0),neg(Q.theta(m))]
        rem=sub(V(2),dot(rhs,inv,rhs));assert rem[0]>0
        factor=sqrt(div(rem,q));cp=add(dot(rhs,inv,hf),mul(V(s),mul(factor,q)))
        positive=[]
        for j in(4,47):
            z=add(dot(rhs,inv,cov[j]),mul(V(s),mul(factor,sub(hp[j],dot(hf,inv,cov[j])))))
            saved=V(r['products'][str(j)]);assert max(z[0],saved[0])<=min(z[1],saved[1])
            if z[0]>0:positive.append(j)
        assert positive==r['positive_original_blockers']and positive
        saved=V(r['center_product']);assert max(cp[0],saved[0])<=min(cp[1],saved[1])
    assert d['deletion_debit']=={'deleted_original_walls':[],'old_finite_vertices_lost':0}
    assert d['scope']['common_center_assumption_needed']is False
    assert d['scope']['added_walls75_76_77_needed']is False
    assert d['scope']['all_norm_two_roots_not_only_seven_profiles']is True
    assert d['scope']['finite_volume_completion_excluded']is False
    print('PASS',TAG,'4 exact branches, both sheets excluded by original positive products')

if __name__=='__main__':main()
