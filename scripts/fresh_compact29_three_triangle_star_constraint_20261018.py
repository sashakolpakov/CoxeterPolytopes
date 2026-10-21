#!/usr/bin/env python3
"""Complete three-angle triangle domains at an E7^3+E8 compact vertex."""
import argparse,json
from hashlib import sha256
from itertools import combinations,product
import fresh_compact29_small_faces_20261015 as P
import fresh_compact29_small_faces_20261015_scout as S
import fresh_compact29_coupled_first_roots_20261018 as C

ROOT=P.ROOT;TAG='fresh_compact29_three_triangle_star_constraint_20261018'
MARKS=[(2,4),(2,5),(3,5)]

def run():
    prior=json.loads((ROOT/f'results/{C.TAG}.json').read_text())
    original=json.loads((ROOT/f'results/{P.TAG}.json').read_text())
    for d in(prior,original):
        for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    n,e=P.diagram('E7');g=P.base_gram(n,e,n);J=P.inverse(g);rows=[]
    degree={i:sum(i in edge for edge in e)for i in range(n)}
    height=P.K((P.F(1,4),0,P.F(3,4),0))
    for a,b in MARKS:
        assert (a,b)not in e
        # Any high contact at a white root would force that root isolated
        # after both distinct marked deletions. At a mark it requires an
        # adjacent leaf. These nonadjacent pairs have neither possibility.
        assert all(not(degree[x]==1 and tuple(sorted((x,y)))in e)for x,y in((a,b),(b,a)))
        eligible=[];rejected=0;ledger=sha256()
        for ms in product(range(2,6),repeat=7):
            ee=e|{(i,7):m for i,m in enumerate(ms)if m>2}
            corners=[S.finite([i for i in range(8)if i!=omit],ee)for omit in(a,b)]
            ledger.update(json.dumps([ms,corners],separators=(',',':')).encode()+b'\n')
            if not all(corners):rejected+=1;continue
            contacts=[-P.COS[m]for m in ms];h=P.dot(contacts,contacts,J)-2
            row=dict(labels=ms,corner_finite=corners,height_squared=h)
            if h.sign()<=0:row['status']='nonpositive_height';eligible.append(row);continue
            gg=P.base_gram(7,ee,8);parabolic=None;tests=0
            for size in range(1,9):
                for ix in combinations(range(8),size):
                    tests+=1
                    if P.psd_status(P.sub(gg,ix))==0:parabolic=list(ix);break
                if parabolic:break
            row['principal_tests']=tests
            if parabolic:row.update(status='parabolic_obstruction',parabolic_subset=parabolic)
            else:
                assert ms==(2,2,2,2,2,5,2)and h==height
                row.update(status='unique_necessary_triangle_root',compact_principals_all_nonsingular=True)
            eligible.append(row)
        counts={s:sum(r['status']==s for r in eligible)for s in('nonpositive_height','parabolic_obstruction','unique_necessary_triangle_root')}
        assert counts['nonpositive_height']==3 and counts['unique_necessary_triangle_root']==1
        assert(len(eligible),counts['parabolic_obstruction'])==((12,8)if(a,b)==(2,4)else(7,3))
        assert len(eligible)+rejected==4**7
        frozen=next(r for r in original['cases']if r['recipe']['type']=='E7'and r['recipe']['marked']==[a,b])
        assert frozen['recipe']['kind']=='triangle'and frozen['recipe']['contacts']==[[5,5]]
        rows.append(dict(marked=[a,b],all_low_vectors=4**7,finite_rejections=rejected,
            finite_filter_ledger_sha256=ledger.hexdigest(),finite_eligible=eligible,status_counts=counts,
            no_high_contact_possible=True,frozen_local_witness=frozen['recipe']))
    squared=height*height;assert squared==P.K((P.F(23,8),0,P.F(3,8),0))
    lo11,hi11=C.cosine_square_bounds(11);lo12,hi12=C.cosine_square_bounds(12)
    assert hi11<P.F(37,10)and(squared-P.F(37,10)).sign()>0
    assert(P.F(93,25)-squared).sign()>0 and P.F(93,25)<lo12
    paths=[f'scripts/{TAG}.py',f'scripts/{S.TAG}.py']+[f'{folder}/{s}.{ext}'for s in(P.TAG,C.TAG)
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return P.enc(dict(status='PASS_COMPACT29_UNIVERSAL_THREE_TRIANGLE_STAR_CONSTRAINT',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        forest=['E7','E7','E7','E8'],ambient_dimension=29,marked_pairs=MARKS,E7_Gram=g,E7_inverse=J,
        complete_triangle_domains=rows,all_low_vectors=3*4**7,total_finite_eligible=26,
        total_nonpositive_heights=9,total_parabolic_obstructions=14,total_unique_roots=3,
        forced_contact_node=5,forced_contact_label=5,forced_height_squared=height,
        cross_factor_squared_product=squared,native_gap=dict(labels=[11,12],
            label11_squared=[lo11,hi11],label12_squared=[lo12,hi12],lower='37/10',upper='93/25'),
        triangle_factor_capacity=1,designated_angles=9,designated_triangle_upper=3,
        designated_nonsimplex_faces_lower=6,
        geometric_scope=dict(other_component_contacts_excluded_by_finite_E8_maximality=True,
            all_labels_for_these_three_triangle_angles=True,compactness_used_to_exclude_parabolics=True,
            no_arithmeticity_assumption=True,quadrilateral_stars_classified=False,
            new_face_charge_saving_proved=False,new_compact_dimension_bound=None)))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert a.verify or not p.exists();out=run()
    if a.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['all_low_vectors'],out['designated_triangle_upper'],flush=True)
