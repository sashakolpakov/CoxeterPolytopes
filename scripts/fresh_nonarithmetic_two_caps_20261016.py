#!/usr/bin/env python3
"""Exact first-stage obstruction for the eight frozen two-cap recipes.

Normals are kept as sqrt(s)*p + sign*sqrt(d)*k.  Every comparison is reduced
to signs in the frozen octic field; no finite-angle alphabet is assumed.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_h9_joint_corner_20261015_bridge_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_two_caps_20261016'
TARGET='results/fresh_nonarithmetic_escape_target_20261016.json'
BASE='results/fresh_h9_joint_pair_20261015_two_extension.json'
BRIDGE='results/fresh_h9_joint_corner_20261015_bridge.json'
C,K,Z,O=R.C,R.K,R.Z,R.O
T=[0,1,2,5,6,10,30]

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d['sources'].items():assert sha(f)==h,f
    return d
def radical_sign(a,b,s,d,certificate=False):
    """Sign of a*sqrt(s)+b*sqrt(d), for nonnegative s,d in the real field."""
    assert R.sg(s)>=0 and R.sg(d)>=0
    sa=R.sg(a) if s else 0;sb=R.sg(b) if d else 0
    difference=None;difference_sign=None
    if not sa:answer=sb;reason='first term zero'
    elif not sb:answer=sa;reason='second term zero'
    elif sa==sb:answer=sa;reason='same signs'
    else:
        difference=s*a*a-d*b*b;difference_sign=R.sg(difference)
        answer=sa*difference_sign;reason='opposite signs; compare squared magnitudes'
    if not certificate:return answer
    return dict(a=R.enc(a),b=R.enc(b),s=R.enc(s),d=R.enc(d),
        first_term_sign=sa,second_term_sign=sb,
        squared_difference=None if difference is None else R.enc(difference),
        squared_difference_sign=difference_sign,sign=answer,reason=reason)

def run():
    target,base,bridge=map(read,(TARGET,BASE,BRIDGE))
    assert target['status']=='PASS_BOUNDED_NONARITHMETIC_ESCAPE_TARGET_GEOMETRY'
    proposed=target['proposed_bounded_experiment']
    assert proposed['support']==T and proposed['first_cap_extra_anchor']==[3,7]
    assert proposed['first_anchor_labels']==[2,3,4,5]
    assert proposed['algebraic_recipe_count']==128 and proposed['maximum_isolated_branch_pairs']==512
    Q,N,H=R.m4(base['form']),R.m4(base['normals']),R.v4(base['center'])
    N.append(R.v8(bridge['cap']['normal']))
    E=R.v4(base['strict_escape']['point'])
    assert len(N)==33
    pair=lambda x,y:C.pair(x,y,Q)
    cov=[C.mv(Q,n) for n in N]
    norms=[pair(n,n) for n in N]
    gram=[[C.dot(n,c) for c in cov] for n in N]
    assert all(R.sg(q)>0 for q in norms)
    assert all(R.sg(gram[i][j])<=0 for i in range(33) for j in range(i))
    assert R.sg(pair(H,H))<0 and R.sg(pair(E,E))>0
    assert all(R.sg(pair(v,n))<0 for v in (H,E) for n in N)
    phi=(O+R.BASIS[1])/2
    assert norms[3]==norms[7]==norms[31]==2*O
    svalues={2:Z,3:O,4:2*O,5:phi*phi}
    finite_angles={2:Z,3:O,4:2*O,5:phi*phi,6:3*O}
    # T+31 is spherical. Every norm2 normal orthogonal to it yields a PD9
    # endpoint Gram, regardless of the separate old-wall compatibility gates.
    endpoint_ids=T+[31]
    endpoint_gram=[[gram[i][j] for j in endpoint_ids] for i in endpoint_ids]
    positive,endpoint_pivots=R.pd(endpoint_gram);assert positive
    JT=C.inverse([[gram[i][j] for j in T] for i in T])
    projected31=norms[31]-C.dot([gram[31][i] for i in T],C.mv(JT,[gram[i][31] for i in T]))
    assert projected31==O/2
    frames=[];recipes=[];counts=Counter()
    # Each expected entry is an exact interval between CONSECUTIVE labels.
    gaps={(3,3):(2,3),(3,4):(4,5),(3,5):(7,8),
          (7,2):(2,3),(7,3):(3,4),(7,4):(3,4),(7,5):(4,5)}
    for anchor in (3,7):
        ids=T+[31,anchor];rows=[cov[i] for i in ids]
        _,pivots=C.eliminate(rows);assert len(pivots)==9
        k=C.null_one(rows);kk=pair(k,k);assert R.sg(kk)>0
        p=C.mv(C.inverse(rows+[C.mv(Q,k)]),[Z]*8+[-O,Z])
        pp=pair(p,p);assert R.sg(pp)<0 and not pair(p,k)
        assert C.mv(rows,k)==[Z]*9 and C.mv(rows,p)==[Z]*8+[-O]
        assert any(k)
        old_a=[pair(p,n) for n in N];old_b=[pair(k,n) for n in N]
        frame_index=len(frames)
        frames.append(dict(anchor=anchor,linear_wall_ids=ids,linear_rank=9,
            pivot_columns=pivots,particular=R.enc(p),kernel=R.enc(k),
            particular_norm=R.enc(pp),kernel_norm=R.enc(kk),cross_pairing=R.enc(Z),
            old_pairing_a=R.enc(old_a),old_pairing_b=R.enc(old_b),
            center_coefficients=R.enc([pair(p,H),pair(k,H)]),
            escape_coefficients=R.enc([pair(p,E),pair(k,E)])))
        for m,s in svalues.items():
            d=(2*O-s*pp)/kk
            assert R.sg(d)>0 and s*pp+d*kk==2*O
            branches=[]
            for sign in (-1,1):
                signs=[radical_sign(a,sign*b,s,d) for a,b in zip(old_a,old_b)]
                assert all(signs[i]==0 for i in T+[31])
                assert signs[anchor]==(0 if m==2 else -1)
                center_sign=radical_sign(pair(p,H),sign*pair(k,H),s,d)
                escape_sign=radical_sign(pair(p,E),sign*pair(k,E),s,d)
                obtuse=[i for i,value in enumerate(signs) if value>0]
                row=dict(sign=sign,old_wall_pairing_signs=signs,center_sign=center_sign,
                    escape_sign=escape_sign,norm=R.enc(2*O),positive_old_pairings=obtuse)
                if obtuse:
                    blocker=obtuse[0]
                    cert=radical_sign(old_a[blocker],sign*old_b[blocker],s,d,True)
                    assert cert['sign']==1
                    row.update(status='REJECT_OBTUSE',blocker_wall=blocker,blocker_sign_certificate=cert)
                    counts['obtuse_branches']+=1
                else:
                    assert sign==-1 and center_sign==-1
                    blocker=7 if anchor==3 else 3
                    a,b=old_a[blocker],sign*old_b[blocker]
                    assert signs[blocker]==-1
                    rad=s*d
                    ca=2*(s*a*a+d*b*b)/norms[blocker]
                    cb=4*a*b/norms[blocker]
                    # c=4(u,n)^2/(u^2*n^2)=ca+cb*sqrt(rad).
                    assert radical_sign(ca,cb,O,rad)>0
                    assert radical_sign(ca-4*O,cb,O,rad)<0
                    lower_label,upper_label=gaps[anchor,m]
                    if (lower_label,upper_label)==(7,8):
                        lower,upper=K.convert(F(33,10)),K.convert(F(10,3))
                        kind='rational bracket strictly inside the label7/label8 gap'
                    else:
                        lower,upper=finite_angles[lower_label],finite_angles[upper_label]
                        kind='exact consecutive finite-angle Cartan cycles'
                    lower_cert=radical_sign(ca-lower,cb,O,rad,True)
                    upper_cert=radical_sign(ca-upper,cb,O,rad,True)
                    assert lower_cert['sign']==1 and upper_cert['sign']==-1
                    row.update(status='REJECT_ALL_LABEL_ANGLE_GAP',blocker_wall=blocker,
                        cartan_cycle=dict(a=R.enc(ca),b=R.enc(cb),radicand=R.enc(rad)),
                        consecutive_labels=[lower_label,upper_label],gap_kind=kind,
                        lower_bound=R.enc(lower),upper_bound=R.enc(upper),
                        lower_comparison=lower_cert,upper_comparison=upper_cert)
                    counts['acute_angle_gap_branches']+=1
                branches.append(row);counts['isolated_branches']+=1
            recipes.append(dict(anchor=anchor,label=m,frame=frame_index,s=R.enc(s),d=R.enc(d),
                normal_recipe='u=sqrt(s)*particular + sign*sqrt(d)*kernel; both square roots nonnegative',
                real_branches=2,degenerate=False,branches=branches))
    assert len(recipes)==8 and counts==dict(obtuse_branches=9,acute_angle_gap_branches=7,isolated_branches=16)
    # For c7=2+2cos(2pi/7), P(c7)=0 and 3<c7<4. P is increasing on [3,infinity):
    # P'(3)=3 and P''(c)=6c-10>=8. Thus P(33/10)>0 proves c7<33/10.
    polynomial=lambda c:c**3-5*c*c+6*c-1
    assert polynomial(F(3))==-1 and polynomial(F(33,10))==F(287,1000)
    assert 3*F(3)**2-10*F(3)+6==3
    assert 6*F(3)-10==8
    # c8=2+sqrt2 >10/3 because 2>(4/3)^2.
    assert F(2)-F(4,3)**2==F(2,9)>0
    paths=[f'scripts/{TAG}.py',TARGET,BASE,BRIDGE,
        'scripts/fresh_h9_joint_corner_20261015_bridge_review.py',
        'scripts/fresh_h9_h8_all_high_20261014_review.py',
        'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_COMPLETE_BOUNDED_TWO_CAP_FIRST_STAGE_OBSTRUCTION',
        sources={p:sha(p) for p in paths},
        arithmetic='Formal two-radical vectors over the frozen exact octic field; opposite-sign comparisons reduce to exact squared-magnitude comparisons with 384-bit outward rational sign bounds.',
        old_wall_count=33,target_support=T,first_endpoint_wall=31,
        endpoint_spherical_gram=R.enc(endpoint_gram),endpoint_positive_pivots=R.enc(endpoint_pivots),
        projected_wall31_norm=R.enc(projected31),linear_frames=frames,recipes=recipes,
        label7_label8_gap=dict(polynomial_coefficients=[1,-5,6,-1],
            polynomial_at3='-1',polynomial_at33_over10='287/1000',
            first_derivative_at3='3',second_derivative_lower_bound_on_ge3='8',
            lower_rational='33/10',upper_rational='10/3',
            label8_square_margin='2/9',
            proof='c7 satisfies the displayed cubic and is >3; the cubic is strictly increasing on [3,infinity), so c7<33/10. Also c8=2+sqrt2>10/3. Every last-branch cycle lies strictly between these rationals.'),
        counts=dict(first_stage_recipes=8,**counts,degenerate_linear_systems=0,
            continuous_solution_sets=0,old_wall_sign_checks=528,accepted_first_caps=0,
            unresolved_first_caps=0,second_stage_recipes_attempted=0,
            ordered_two_cap_recipes_excluded=128),
        conclusion='Every isolated branch of every proposed first-cap recipe is incompatible with the retained33 walls. Therefore none of the128 frozen ordered two-cap recipes can succeed.',
        scope='Complete only for the literal T-perpendicular norm2 first-cap recipes, right to31 and label2..5 to anchor3or7. Exclusions use obtuseness or gaps between consecutive integer Coxeter labels, not an assumed field alphabet. No claim for other anchor labels, nonright support contacts, other supports, changed retained walls, larger clusters, or arbitrary finite-volume completions. No second-cap search or new wall was constructed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text());print('PASS saved first-stage obstruction replay')
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
