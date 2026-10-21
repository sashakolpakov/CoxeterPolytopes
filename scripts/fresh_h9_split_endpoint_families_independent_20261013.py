#!/usr/bin/env python3
"""Independent Fraction audit of the 728 endpoint-cap entries.

No producer imports. Uses the earlier independent Q(phi) arithmetic and
integer-square-root interval helper, with every used enclosure checked by
exact field signs. Pairwise interval separation plus exact compositum
equality proves distinctness and literal overlap with the older bank.
"""
import argparse
from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import fresh_h9_split_leaf_3040_independent_20261013 as K
import fresh_h9_split_leaf_3034_safe_graph_20261013 as I

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h9_split_endpoint_families_20261013_exact.json'
OLD=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
OUT=ROOT/'results/fresh_h9_split_endpoint_families_independent_20261013.json'
Z,O,P=K.Z,K.O,K.PHI
WITNESS=[8,8,8,1,-3,-3,-3,-3,-3,-3,-3,5,5,5,5,5]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def parse_vector(v):
    assert len(v)==16 and all(len(x)==2 for x in v)
    return tuple(tuple(map(K.parse,x))for x in v)


def seed():
    a=(F(1,2),F(1));b=K.neg(K.mul(a,K.inv(K.add(O,K.sc(a,12)))))
    assert b==(-F(5,22),F(1,11))
    eps=[1]*4+[-1]*4
    g=tuple(tuple(K.sc(O,2)if i==j else Z if i//8==j//8 else
        K.neg(K.add(K.add(a,K.sc(b,eps[i%8]*eps[j%8])),K.sc(O,2 if i%8==j%8 else 0)))
        for j in range(16))for i in range(16))
    assert all(K.sign(K.add(g[i][j],K.sc(O,2)))<0 for i in range(8)for j in range(8,16))
    center=(O,)*16
    assert K.mv(g,center)==(K.sc(a,-8),)*16
    square=K.dot(center,K.mv(g,center));assert square==K.sc(a,-128)and K.sign(square)<0
    # An explicit congruence basis proves inertia (9,1,6), hence H9 after
    # quotienting the six-dimensional coefficient kernel.
    s=[1]*16;e=eps+eps;d=[1]*8+[-1]*8;f=eps+[-x for x in eps]
    positive=[];kernel=[]
    for start in(0,4):
        for j in range(start+1,start+4):
            p=[0]*16;k=[0]*16
            p[j]=p[start+8]=1;p[start]=p[j+8]=-1
            k[j]=k[j+8]=1;k[start]=k[start+8]=-1
            positive.append(p);kernel.append(k)
    basis=[s,e,d,f]+positive+kernel
    bg=[[K.dot([K.sc(O,x)for x in v],K.mv(g,[K.sc(O,x)for x in w]))for w in basis]for v in basis]
    diag=[K.sc(a,-128),K.sc(b,-128),K.add(K.sc(O,64),K.sc(a,128)),K.add(K.sc(O,64),K.sc(b,128))]
    assert [K.sign(x)for x in diag]==[-1,1,1,1]
    for i in range(16):
        for j in range(16):
            expected=diag[i]if i==j and i<4 else K.sc(O,16 if i==j else 8)if 4<=i<10 and 4<=j<10 and(i-4)//3==(j-4)//3 else Z
            assert bg[i][j]==expected
    # Each three-difference block plus its sum spans its four coordinates;
    # sum/difference of U and V is invertible, proving this is a basis.
    return g,square


def check_cap(raw,scout,g):
    D=K.parse(raw['radicand']);field=K.Quadratic(D)
    normal=parse_vector(raw['normal']);cov=parse_vector(raw['covector'])
    computed=tuple((K.dot(row,[x[0]for x in normal]),K.dot(row,[x[1]for x in normal]))for row in g)
    assert computed==cov
    assert field.sign(field.sub(field.dot(normal,cov),(K.sc(O,2),Z)))==0
    assert field.sign(field.total(cov))<0
    side=raw['seed_side'];assert type(side)is int and side in(0,1)
    assert side==scout['seed_side']and raw['finite_profile']==scout['finite_profile']
    labels=[];counts=Counter()
    for i,v in enumerate(cov):
        label=next((m for m,c in((2,Z),(3,K.neg(O)),(5,K.neg(P)))if field.sign(field.sub(v,(c,Z)))==0),None)
        if label is None:
            assert field.sign(field.add(v,(K.sc(O,2),Z)))<0
            counts['strictly_dotted']+=1
        else:counts['pi_over_'+str(label)]+=1
        labels.append(label)
        # Independent positive-root choice: exact comparison with rational
        # neighborhoods of the original decimal proposal, not float math.
        target=F(str(scout['seed_covector'][i]));tol=F(1,10**6)
        assert field.sign(field.sub(v,((target-tol,F(0)),Z)))>0
        assert field.sign(field.sub(v,((target+tol,F(0)),Z)))<0
    assert labels==raw['finite_seed_labels']
    # Orthogonal eight-wall block + cap: Schur complement q is positive,
    # so the resulting nine-wall diagram is spherical, with full rank9.
    local=cov[side*8:side*8+8]
    assert all(labels[i]is not None for i in range(side*8,side*8+8))
    q=field.sub((K.sc(O,2),Z),field.mul(((F(1,2),F(0)),Z),field.dot(local,local)))
    assert field.sign(q)>0
    assert field.sign(field.sub(q,(K.parse(raw['endpoint_projection_norm']),Z)))==0
    # Every single-cap affine subdiagram would be a star within one block:
    # the two seed blocks are strictly dotted to each other. A semidefinite
    # finite-label star requires k3+k5*phi^2=4, thus k5=0 and k3=4.
    assert all(sum(m==3 for m in labels[start:start+8])<4 for start in(0,8))
    escape=field.total(field.mul(((F(v),F(0)),Z),c)for v,c in zip(WITNESS,cov))
    assert escape==tuple(map(K.parse,raw['old_universal_escape_pairing']))
    escape_sign=field.sign(escape);assert escape_sign==raw['old_universal_escape_sign']
    return D,normal,cov,labels,counts,escape_sign


def equal_element(x,D,y,E):
    """Sign zero in K(sqrt(D),sqrt(E)), including dependent radicals."""
    field=K.Quadratic(D)
    a=(K.sub(x[0],y[0]),x[1]);b=(K.neg(y[1]),Z)
    sa,sb=field.sign(a),field.sign(b)
    if E==Z:return sa==0
    if not sa:return sb==0
    if not sb or sa==sb:return False
    gap=field.sub(field.mul(a,a),field.mul((E,Z),field.mul(b,b)))
    return field.sign(gap)==0


def bounds(rows,enc,checked):
    result=[]
    for raw in rows:
        D=K.parse(raw['radicand']);field=K.Quadratic(D);row=[]
        for v in parse_vector(raw['covector']):
            lo,hi=enc.element(v,D);key=(v,D,lo,hi)
            if key not in checked:
                assert field.sign(field.sub(v,((F(lo,I.SCALE),F(0)),Z)))>=0
                assert field.sign(field.sub(((F(hi,I.SCALE),F(0)),Z),v))>=0
                checked.add(key)
            row.append((lo,hi))
        result.append(tuple(row))
    return result


def same_possible(a,b):return all(al<=bh and bl<=ah for(al,ah),(bl,bh)in zip(a,b))
def exact_same(a,b):
    D=K.parse(a['radicand']);E=K.parse(b['radicand'])
    return all(equal_element(x,D,y,E)for x,y in zip(parse_vector(a['covector']),parse_vector(b['covector'])))


def run():
    data=json.loads(SOURCE.read_text());scout_path=ROOT/data['source'];scout=json.loads(scout_path.read_text())
    assert data['source_sha256']==sha(scout_path)
    rows=data['caps'];assert len(rows)==728
    g,center_square=seed();assert g==tuple(tuple(map(K.parse,row))for row in data['gram'])
    totals=Counter();escape_signs=Counter();cutters=[];source_ids=[]
    for i,raw in enumerate(rows):
        assert raw['index']==i and type(raw['scout_index'])is int
        sid=raw['scout_index'];assert 0<=sid<len(scout['caps']);source_ids.append(sid)
        _,_,_,_,counts,sgn=check_cap(raw,scout['caps'][sid],g)
        totals.update(counts);escape_signs[sgn]+=1
        if sgn>0:cutters.append(i)
        if i%150==0:print('Exact endpoint-cap replay',i,flush=True)
    assert len(set(source_ids))==728
    assert data['positive_old_escape_cutters']==len(cutters)==10
    witness=tuple(K.sc(O,v)for v in WITNESS)
    witness_square=K.dot(witness,K.mv(g,witness));witness_center=K.total(K.mv(g,witness))
    assert witness_square==(F(7242,11),F(-2800,11))and K.sign(witness_square)>0
    assert witness_center==(-F(116),-F(232))and K.sign(witness_center)<0
    assert all(K.sign(x)<0 for x in K.mv(g,witness))
    enc=I.Enclosures();checked=set();new_bounds=bounds(rows,enc,checked)
    old=json.loads(OLD.read_text());old_rows=old['caps'];assert len(old_rows)==3040
    assert g==tuple(tuple(map(K.parse,row))for row in old['gram'])
    assert sha(OLD)=='13012086e9203f5dbe23e06143addeb4b0858e851bfcdf7d800d5b4d779ccf44'
    old_bounds=bounds(old_rows,enc,checked)
    duplicates=[];new_unresolved=[]
    for i in range(728):
        for j in range(i):
            if same_possible(new_bounds[i],new_bounds[j]):
                if exact_same(rows[i],rows[j]):duplicates.append([j,i])
                else:new_unresolved.append([j,i])
    overlap=[];old_unresolved=[]
    for i in range(728):
        matches=[]
        for j in range(3040):
            if same_possible(new_bounds[i],old_bounds[j]):
                if exact_same(rows[i],old_rows[j]):
                    # Verify that the old covector here also comes from its
                    # archived normal, rather than relying only on an ID.
                    n=parse_vector(old_rows[j]['normal']);v=parse_vector(old_rows[j]['covector'])
                    assert tuple((K.dot(r,[z[0]for z in n]),K.dot(r,[z[1]for z in n]))for r in g)==v
                    matches.append(j)
                else:old_unresolved.append([i,j])
        if matches:overlap.append({'new_index':i,'old_indices':matches})
    assert not duplicates
    # False-positive overlapping intervals are allowed and then resolved by
    # exact arithmetic; these lists should happen to be empty for this bank.
    assert len(overlap)==32
    forged=deepcopy(rows[0]);forged['normal'][0][0]['a']=str(F(forged['normal'][0][0]['a'])+1)
    try:check_cap(forged,scout['caps'][forged['scout_index']],g)
    except AssertionError:pass
    else:raise AssertionError('Forged normal accepted')
    forged=deepcopy(rows[0]);forged['finite_seed_labels'][0]=3
    try:check_cap(forged,scout['caps'][forged['scout_index']],g)
    except AssertionError:pass
    else:raise AssertionError('Forged seed label accepted')
    forged=deepcopy(rows[0]);forged['old_universal_escape_sign']=1
    try:check_cap(forged,scout['caps'][forged['scout_index']],g)
    except AssertionError:pass
    else:raise AssertionError('Forged escape sign accepted')
    return {'status':'PASS_INDEPENDENT_EXACT_728_ENDPOINT_CAPS','standard_library_only':True,
        'source':str(SOURCE.relative_to(ROOT)),'source_sha256':sha(SOURCE),
        'scout_source':str(scout_path.relative_to(ROOT)),'scout_sha256':sha(scout_path),
        'seed_signature':[9,1,6],'strict_center_square':K.encode(center_square),
        'norm_two_caps':728,'rank_nine_spherical_extensions':728,
        'cap_seed_pair_counts':dict(sorted(totals.items())),
        'caps_without_single_cap_affine_seed_subdiagrams':728,
        'old_escape_square':K.encode(witness_square),'old_escape_center_pairing':K.encode(witness_center),
        'old_escape_sign_counts':dict(sorted(escape_signs.items())),
        'positive_old_escape_cutters':cutters,
        'distinct_endpoint_cap_mirrors':728,'duplicate_pairs':duplicates,
        'old_bank_source':str(OLD.relative_to(ROOT)),'old_bank_sha256':sha(OLD),
        'literal_old_bank_overlap_count':len(overlap),'literal_old_bank_overlap':overlap,
        'distinct_mirrors_outside_old_bank':728-len(overlap),
        'pairwise_interval_separation_checks':728*727//2+728*3040,
        'independently_certified_covector_intervals':len(checked),
        'coefficient_interval_scale':str(I.SCALE),
        'interval_overlaps_resolved_as_unequal':{'within_new_bank':new_unresolved,'with_old_bank':old_unresolved},
        'forged_normal_rejected':True,'forged_seed_label_rejected':True,'forged_escape_sign_rejected':True,
        'scope':'Exact individual cap normals, seed relations, common strict center, spherical nine-wall extensions, absence of single-cap affine seed subdiagrams, and literal mirror-bank comparison only. No cap-cap compatibility, completeness, finite-volume, compactness or global novelty claim.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();result=run()
    if a.verify:assert json.loads(OUT.read_text())==json.loads(json.dumps(result))
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if k not in('literal_old_bank_overlap',)},indent=2))
