#!/usr/bin/env python3
"""All-label two-duty obstruction using only the original physical walls."""
import argparse,json
from hashlib import sha256
from itertools import combinations
import fresh_remote60_exact_original_extension_20261018 as A

ROOT=A.ROOT;TAG='fresh_exact66_old7_boundary_cap_obstruction_20261018'
CURRENT='fresh_exact65_escape33_cap_20261018'
V,add,sub,mul,div,neg,sqrt,dot=A.V,A.add,A.sub,A.mul,A.div,A.neg,A.sqrt,A.dot
FRAME=[7,30,52];BLOCKERS=[4,47]

def run():
    seed=json.loads((ROOT/f'results/{A.TAG}.json').read_text())
    current=json.loads((ROOT/f'results/{CURRENT}.json').read_text())
    for d in(seed,current):
        for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert current['physical']['summary']['total_walls']==66
    assert current['physical']['summary']['complete_finite_vertices']==205
    data=seed['physical_input'];metric=[list(map(V,r))for r in data['metric']]
    roots={int(i):list(map(V,v))for i,v in data['roots'].items()}
    h=list(map(V,data['center']));hc={int(i):V(x)for i,x in data['center_products'].items()}
    white=list(seed['fixed_white_walls']);assert white==[0,1,2,3,5,6]
    wi=[list(map(V,r))for r in seed['white_inverse']]
    wc={int(i):list(map(V,r))for i,r in seed['white_covectors'].items()}
    hw=[hc[i]for i in white];hp={i:sub(hc[i],dot(hw,wi,wc[i]))for i in roots}
    hp2=sub(V(data['center_square']),dot(hw,wi,hw));assert hp2[1]<0
    ids=white+[7,30];identities=[];native={}
    expected={(0,1):5,(1,2):3,(2,3):3,(5,6):3,(6,7):3,(7,30):5,(30,52):4}
    for i,j in combinations(ids+[52],2):
        z,proof=A.refine(dot(roots[i],metric,roots[j]),i,j)
        if proof:identities.append(proof)
        native[i,j]=native[j,i]=z
        if i in ids and j in ids or {i,j}=={30,52}:
            assert proof and proof['identity']==expected.get(tuple(sorted((i,j))),2)
    # W+7+30 is H4+H4. Finite adjoining K must be orthogonal to all
    # eight roots. The other finite duty has a fixed bond4, leaving K52=2/3.
    gram=[[sub(V(2)if i==j else native[i,j],dot(wc[i],wi,wc[j]))for j in FRAME]for i in FRAME]
    inv,piv=A.I.inverse(gram);assert [1 if z[0]>0 else -1 for z in piv]==[1,1,-1]
    hf=[hp[i]for i in FRAME]
    q=sub(hp2,dot(hf,inv,hf));assert q[0]>0
    projected_roots={i:A.vector_sub(roots[i],A.linear(A.row_times(wc[i],wi),[roots[w]for w in white]))for i in FRAME}
    ph=A.vector_sub(h,A.linear(A.row_times(hw,wi),[roots[w]for w in white]))
    residual=A.vector_sub(ph,A.linear(A.row_times(hf,inv),[projected_roots[i]for i in FRAME]))
    cs={}
    for j in BLOCKERS:
        cs[j]=[]
        for i in FRAME:
            z,proof=A.refine(dot(roots[i],metric,roots[j]),i,j)
            if proof:identities.append(proof)
            cs[j].append(sub(z,dot(wc[i],wi,wc[j])))
    rows=[]
    for label in(2,3):
        rhs=[V(0),V(0),neg(A.E.theta(label))]
        rem=sub(V(2),dot(rhs,inv,rhs));assert rem[0]>0
        factor=sqrt(div(rem,q));part=A.linear(A.row_times(rhs,inv),[projected_roots[i]for i in FRAME])
        for sign in(-1,1):
            k=[add(x,mul(V(sign),mul(factor,y)))for x,y in zip(part,residual)]
            cp=add(dot(rhs,inv,hf),mul(V(sign),mul(factor,q)))
            products={};raw={};direct={}
            for j in BLOCKERS:
                raw[j]=add(dot(rhs,inv,cs[j]),mul(V(sign),mul(factor,sub(hp[j],dot(hf,inv,cs[j])))))
                direct[j]=dot(k,metric,roots[j])
                products[j]=(max(raw[j][0],direct[j][0]),min(raw[j][1],direct[j][1]))
                assert products[j][0]<=products[j][1]
            positive=[j for j in BLOCKERS if products[j][0]>0]
            assert positive
            rows.append(dict(label52=label,sheet_relative_to_center_residual=sign,rhs=rhs,remainder=rem,
                factor=factor,normal=k,center_product=cp,schur_products=raw,direct_products=direct,
                products=products,positive_original_blockers=positive))
    assert len(rows)==4
    paths=[f'scripts/{TAG}.py']+[f'{folder}/{s}.{ext}'for s in(A.TAG,CURRENT)
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_EXACT_OLD7_TWO_FINITE_DUTY_CAP_OBSTRUCTION',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        fixed_white_walls=white,first_finite_duty=[7,30,'K'],second_finite_duty=[30,52,'K'],
        forced_H4_factors=[[0,1,2,3],[5,6,7,30]],forced_profile=0,
        forced_contacts={7:2,30:2},complete_remaining_labels52=[2,3],native_isolation=identities,
        frame=FRAME,frame_Gram=gram,frame_inverse=inv,frame_pivots=piv,
        projected_center_square=hp2,projected_center_contacts=hf,positive_center_complement_square=q,
        blocker_frame_covectors=cs,branches=rows,
        deletion_debit=dict(deleted_original_walls=[],old_finite_vertices_lost=0),
        scope=dict(all_finite_labels_for_two_duties=True,all_norm_two_roots_not_only_seven_profiles=True,
            original_wall_positive_products_exclude_both_sheets=True,common_center_assumption_needed=False,
            added_walls75_76_77_needed=False,other_incidence_excluded=False,
            finite_volume_completion_excluded=False,new_geometric_example=False))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert a.verify or not path.exists();out=A.enc(run())
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['branches']),flush=True)
