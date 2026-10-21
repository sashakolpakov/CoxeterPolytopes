#!/usr/bin/env python3
"""Exact interval exclusion of the fixed-scaffold finite bridge family.

No new normal is approximated. Original physical certificate intervals are
projected off all six white roots before an invariant one-line norm formula.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path
import fresh_compact10_cross_bridge_20261018 as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact10_cross_bridge_high_interval_20261018'
PRIMARY='fresh_nonarithmetic_release49_62_replacement_20261017'
WHITE=[0,1,2,3,5,6]
FRAME=[8,10,33]
IDS=WHITE+FRAME+[7]
PREC=10**35


def iv(x):
    if isinstance(x,(tuple,list)):return tuple(map(F,x))
    return F(x),F(x)


def rounded(a):
    lo,hi=a;assert lo<=hi
    return F((lo.numerator*PREC)//lo.denominator,PREC),F(-((-hi.numerator*PREC)//hi.denominator),PREC)


def add(a,b):return rounded((a[0]+b[0],a[1]+b[1]))
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    z=[x*y for x in a for y in b]
    return rounded((min(z),max(z)))
def div(a,b):
    assert b[0]*b[1]>0
    return mul(a,(1/b[1],1/b[0]))
def total(xs):
    out=iv(0)
    for x in xs:out=add(out,x)
    return out
def sqrt(a):
    assert a[0]>=0
    lows=isqrt((a[0].numerator*PREC**2)//a[0].denominator)
    highs=isqrt((a[1].numerator*PREC**2)//a[1].denominator)
    if F(highs,PREC)**2<a[1]:highs+=1
    out=F(lows,PREC),F(highs,PREC)
    assert out[0]**2<=a[0] and out[1]**2>=a[1]
    return out


def inverse(g):
    n=len(g);a=[row[:]+[iv(int(i==j)) for j in range(n)] for i,row in enumerate(g)];pivots=[]
    for i in range(n):
        pivot=a[i][i];assert pivot[0]*pivot[1]>0;pivots.append(pivot)
        a[i]=[div(x,pivot) for x in a[i]]
        for j in range(n):
            if j==i:continue
            q=a[j][i];a[j]=[sub(x,mul(q,y)) for x,y in zip(a[j],a[i])]
    return [row[n:] for row in a],pivots


def dot(x,g,y):return total(mul(mul(a,g[i][j]),b) for i,a in enumerate(x) for j,b in enumerate(y))
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def original_projection(primary):
    h=list(map(iv,primary['common_center']));g=[list(map(iv,row)) for row in primary['computational_gram']]
    normals={}
    for w in primary['facet_witnesses']:
        i=w['wall']
        if i not in IDS:continue
        hp=iv(primary['center_products'][str(i)]);assert hp[1]<0
        normals[i]=[div(mul(iv(2),sub(x,y)),hp) for x,y in zip(h,map(iv,w['point']))]
    assert set(normals)==set(IDS) and len(g)==10
    native=[[dot(normals[i],g,normals[j]) for j in IDS] for i in IDS]
    for i in range(10):assert native[i][i][0]<=2<=native[i][i][1]
    wg=[r[:6] for r in native[:6]];wi,wp=inverse(wg)
    assert all(q[0]>0 for q in wp)
    residual=[[sub(native[i][j],dot(native[i][:6],wi,[native[k][j] for k in range(6)])) for j in range(6,10)] for i in range(6,10)]
    frame=[r[:3] for r in residual[:3]];fi,fp=inverse(frame)
    assert fp[0][0]>0 and fp[1][0]>0 and fp[2][1]<0
    contact=[residual[i][3] for i in range(3)];projection=[total(mul(fi[i][j],contact[j]) for j in range(3)) for i in range(3)]
    rem7=sub(residual[3][3],dot(contact,fi,contact));assert rem7[0]>0
    # The inherited exact native control isolates these side profiles. These
    # comparisons check the physical/projection interface, not new equalities
    # inferred merely because a floating-point value is small.
    expected={8:[0]*6,10:[0,0,0,0,-1,0],33:[0]*6,7:[0,0,0,0,0,-1]}
    for i,w in enumerate(FRAME+[7],6):
        for j,x in enumerate(expected[w]):assert native[i][j][0]<=x<=native[i][j][1]
    assert residual[1][1][0]<=F(4,3)<=residual[1][1][1]
    assert residual[3][3][0]<=F(4,3)<=residual[3][3][1]
    assert residual[0][2][1]<-2
    return dict(walls=IDS,physical_gram=native,white_inverse=wi,white_pivots=wp,
        residual_walls=FRAME+[7],residual_gram=residual,frame_inverse=fi,frame_pivots=fp,
        seven_projection_coefficients=projection,seven_residual_square=rem7,
        old_profiles={str(i):v for i,v in expected.items()})


def run():
    prior=json.loads((ROOT/f'results/{B.TAG}.json').read_text())
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert prior['status']=='BOUNDED_COMPACT10_CHANGED_CROSS_BRIDGE'
    assert primary['status']=='PASS_EXACT_JOINT49_62_REPLACEMENT_63_WALLS_199_FINITE'
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    symbolic=B.T.E.recipes((0,1,0),tuple([2]*9),((0,1),(1,2)))
    high=[(pi,list(ms),list(h)) for pi,ms,h in symbolic if h]
    expected=[(0,[a,2,b],[i for i,m in enumerate((a,2,b)) if m==6]) for a,b in product(range(2,7),repeat=2) if max(a,b)==6]
    assert high==expected and len(high)==prior['bank']['symbolic_high_families_not_tested']==9
    assert len(symbolic)-len(high)==len(prior['bank']['recipes'])==42
    geom=original_projection(primary);fi=geom['frame_inverse'];v=geom['seven_projection_coefficients'];s7=geom['seven_residual_square']
    rt2,rt3,rt5=(sqrt(iv(n)) for n in (2,3,5));phi=div(add(iv(1),rt5),iv(2))
    magnitudes={2:iv(0),3:iv(1),4:rt2,5:phi,6:(rt3[0],F(2))}
    rows=[]
    for pi,labels,slots in high:
        a,b=magnitudes[labels[0]],magnitudes[labels[2]]
        rhs=[neg(a),iv(0),neg(b)]
        linear=total(mul(rhs[i],v[i]) for i in range(3))
        remainder=sub(iv(2),dot(rhs,fi,rhs));assert remainder[0]>0
        radical=sqrt(mul(remainder,s7))
        minus,plus=sub(linear,radical),add(linear,radical)
        assert plus[0]>0
        gaplabels=(4,5) if labels[2]==2 else (5,6)
        uppernegative=neg(magnitudes[gaplabels[0]])
        lowernegative=neg(rt3 if gaplabels[1]==6 else phi)
        # Strict separation from both adjacent native values. Higher labels,
        # parallel and dotted products lie still farther left than this gap.
        margin_left=minus[0]-lowernegative[1]
        margin_right=uppernegative[0]-minus[1]
        assert margin_left>0 and margin_right>0,(labels,enc(minus),gaplabels)
        rows.append(dict(index=len(rows),profile=pi,symbolic_labels=labels,high_slots=slots,
            parameter_a=a,parameter_b=b,linear=linear,norm_remainder=remainder,radical=radical,
            plus_native7=plus,minus_native7=minus,plus_positive_margin=plus[0],
            excluded_between_native_labels=list(gaplabels),minus_left_margin=margin_left,minus_right_margin=margin_right))
    # Certify the already-enumerated low recipes with the same physical
    # metric. The seven prescribed white covectors are finite point domains;
    # these checks do not enumerate another numerical bank.
    low=[];native=geom['physical_gram'];wi=geom['white_inverse']
    profiles=[[iv(0)]*6]
    for pos in (4,5):
        for m in (3,4,5):
            p=[iv(0)]*6;p[pos]=neg(magnitudes[m]);profiles.append(p)
    lows=[(pi,list(ms)) for pi,ms,h in symbolic if not h]
    assert lows==[(r['profile'],r['labels']) for r in prior['bank']['recipes']]
    decisions={'positive':0,'native_gap':0}
    for ri,(pi,labels) in enumerate(lows):
        wc=profiles[pi];square=sub(iv(2),dot(wc,wi,wc));assert square[0]>0
        offsets=[dot(wc,wi,[native[k][6+i] for k in range(6)]) for i in range(4)]
        rhs=[sub(neg(magnitudes[m]),offsets[i]) for i,m in enumerate(labels)]
        linear=add(total(mul(rhs[i],v[i]) for i in range(3)),offsets[3])
        remainder=sub(square,dot(rhs,fi,rhs));assert remainder[0]>0
        radical=sqrt(mul(remainder,s7));branches=[]
        for sign,value in ((-1,sub(linear,radical)),(1,add(linear,radical))):
            if value[0]>0:
                proof=dict(kind='positive',margin=value[0])
            else:
                gaps=[]
                for m in (2,3,4,5):
                    lower=neg(magnitudes[m+1] if m+1<6 else rt3);upper=neg(magnitudes[m])
                    left=value[0]-lower[1];right=upper[0]-value[1]
                    if left>0 and right>0:gaps.append(dict(kind='native_gap',labels=[m,m+1],left_margin=left,right_margin=right))
                assert len(gaps)==1,(pi,labels,enc(value));proof=gaps[0]
            decisions[proof['kind']]+=1;branches.append(dict(sign=sign,native7=value,exclusion=proof))
        # Comparison to old numerics is an interface diagnostic only. The
        # exclusion itself uses the exact interval inequalities above.
        old=sorted(dict(b['native_products'])[7] for b in prior['bank']['recipes'][ri]['branches'])
        now=sorted(float(sum(b['native7'])/2) for b in branches)
        assert max(abs(x-y) for x,y in zip(old,now))<2e-7
        low.append(dict(index=ri,profile=pi,labels=labels,white_covector=wc,projected_square=square,
            white_offsets=offsets,rhs=rhs,linear=linear,norm_remainder=remainder,radical=radical,branches=branches))
    assert len(low)==42 and decisions==dict(positive=45,native_gap=39)
    debit=prior['exact_deletion_debit'];assert debit['lost_finite_count']==16 and debit['compact_faces_by_dimension']=={'2':12,'3':1}
    assert sum(8 in v for v in debit['lost_finite_supports'])==2
    paths={Path(__file__),Path(B.__file__),Path(B.T.E.__file__),Path(B.T.P.__file__),
           ROOT/f'results/{B.TAG}.json',ROOT/f'results/{B.TAG}.md',ROOT/f'results/{PRIMARY}.json',ROOT/f'results/{PRIMARY}.md'}
    return enc(dict(status='PASS_EXACT_INTERVAL_HIGH_BRIDGE_EXCLUSION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        precision_denominator=PREC,frame=FRAME,finite_corner_duties=[[8,10,67],[10,33,67]],nonedge=[8,33],
        old_control_walls=63,retained_old_walls=61,deleted_walls=[30,44],exact_deletion_debit=debit,
        projection=geom,root_intervals=dict(sqrt2=rt2,sqrt3=rt3,sqrt5=rt5,phi=phi),
        symbolic_high_families=high,continuous_boxes=rows,covered_high_families=9,covered_high_norm_signs=18,
        low_point_domains=low,low_branch_decisions=decisions,covered_low_recipes=42,covered_low_norm_signs=84,
        total_symbolic_domains=51,total_norm_sign_domains=102,
        two_parameter_boxes=1,minimum_rational_separation=min(min(r['plus_positive_margin'],r['minus_left_margin'],r['minus_right_margin']) for r in rows),
        new_normals=0,new_cone_checks=0,prior_control_open_ends=prior['control']['open_ends'],
        prior_deletion_provisional=prior['deletion']['score'] is None,
        low_label_scope='The predecessor supplied the complete42 low recipe domain; this successor independently encloses both native7 products using original physical intervals and certifies their nonnative gaps or positive signs.',
        scope='Exact rational interval necessary obstruction for all42 low recipes and all9 high-label families on this fixed scaffold, including independent simultaneous labels and closed infinite-label limits. Inherited seven-profile/finite-diagram classification and original exact physical control are explicit inputs. Both algebraic norm signs fail retained wall7 before native-center/cone acceptance. No new normal, replacement chamber, arbitrary-incidence nonexistence, compactness or H9 finite-volume claim.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'families',out['covered_high_families'],'low_recipes',out['covered_low_recipes'],'signs',out['total_norm_sign_domains'],
          'minimum_separation',out['minimum_rational_separation'],flush=True)
