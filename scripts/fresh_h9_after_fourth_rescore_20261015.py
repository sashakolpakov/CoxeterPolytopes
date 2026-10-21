#!/usr/bin/env python3
"""Bounded numerical rescore after the fourth cap, plus two exact field tests.

No temporary input or project arithmetic import. The old 68111-row numerical
bank is inherited; this does not regenerate its bounded proposal search.
"""
import argparse,gzip,json,math
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import numpy as np
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_after_fourth_rescore_20261015'
BASE='results/fresh_h9_fourth_escape_extension_20261015.json'
BTAG='fresh_h9_after_third_vertex_free_scout_20261015'
BANK=f'data/enumerations/{BTAG}_bank.json.gz'
RAD=6+s.Rational(16,5)*s.sqrt(5)
K=s.QQ.algebraic_field(s.sqrt(5),s.sqrt(RAD));Z,O=K.zero,K.one
BASIS=[O,K.from_sympy(s.sqrt(5)),K.from_sympy(s.sqrt(RAD)),K.from_sympy(s.sqrt(5*RAD))]

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def dec(x):return sum((K.convert(F(v))*b for v,b in zip(x[0]+x[1],BASIS)),Z)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Z)
def mv(M,v):return [dot(r,v)for r in M]
def dm(M):return DomainMatrix(M,(len(M),len(M[0])),K)
def coeff(x):return [F(str(v))for v in [s.S.Zero]*(4-len(x.to_list()))+x.to_list()]
EB=s.Matrix.hstack(*[s.Matrix(coeff(b))for b in BASIS]).inv()
def enc(x):
    if isinstance(x,list):return [enc(a)for a in x]
    c=list(EB*s.Matrix(coeff(x)));return[[str(c[0]),str(c[1])],[str(c[2]),str(c[3])]]
def scalar(x):
    sq=mp.sqrt(5);r=6+16*sq/5
    return mp.mpf(x[0][0])+mp.mpf(x[0][1])*sq+(mp.mpf(x[1][0])+mp.mpf(x[1][1])*sq)*mp.sqrt(r)
def smooth(x):return float(format(float(x),'.14g'))
def rational_square(x):
    if x<0:return False
    a,b=math.isqrt(x.numerator),math.isqrt(x.denominator)
    return a*a==x.numerator and b*b==x.denominator

def recipe_check(d,bank,cutter_ids):
    Q=[[dec(x)for x in r]for r in d['form']]
    N=[[dec(x)for x in r]for r in d['normals']]
    cov=[mv(Q,n)for n in N]
    pair=lambda v,w:dot(v,mv(Q,w))
    phi=(O+BASIS[1])/2
    records=[]
    for label,a,b in [(3,F(50),F(22)),(5,F(245,4),F(111,4))]:
        matches=[i for i in cutter_ids if bank[i]['contacts']==[[27,label]]
                 and bank[i].get('anchor')==2 and bank[i].get('anchor_label')==3]
        assert len(matches)==1
        idx=matches[0];row=bank[idx];support=row['support'];ids=support+[2]
        assert support==[0,1,3,5,6,7,10,27]and len(set(ids))==9
        M=[cov[i]for i in ids];rr,pivots=dm(M).rref();rr=rr.to_list()
        assert len(pivots)==9
        free=next(i for i in range(10)if i not in pivots)
        v=[Z]*10;v[free]=O
        for i,j in enumerate(pivots):v[j]=-rr[i][free]
        assert all(not dot(r,v)for r in M)
        rhs=[-phi if i==27 and label==5 else -O if i in(2,27)else Z for i in ids]
        p=mv(dm(M+[mv(Q,v)]).inv().to_list(),rhs+[Z])
        assert all(dot(cov[i],p)==x for i,x in zip(ids,rhs))and not pair(p,v)
        vq=pair(v,v);pq=pair(p,p);delta=(2*O-pq)/vq
        expected=K.convert(a)+K.convert(b)*BASIS[1]
        assert delta==expected and scalar(enc(delta))>0
        # For delta in F, a square in F(sqrt(D)) is either a square in F,
        # or delta/D is a square in F. Both are ruled out by their F-norms.
        normdelta=a*a-5*b*b;normrad=F(-76,5);normratio=normdelta/normrad
        assert not rational_square(normdelta)and not rational_square(normratio)
        assert normdelta==(F(80)if label==3 else F(-395,4))
        assert normratio==(F(-100,19)if label==3 else F(1975,304))
        # Numerical branch matching is recorded separately from the exact identities.
        vals=[]
        for sign in[-1,1]:
            x=[scalar(enc(u))+sign*mp.sqrt(scalar(enc(delta)))*scalar(enc(z))for u,z in zip(p,v)]
            gg=[sum(x[j]*scalar(enc(cov[i][j]))for j in range(10))*mp.sqrt(2/scalar(d['norms'][i]))for i in range(28)]
            vals.append((max(abs(float(g)-t)for g,t in zip(gg,row['normalized_wall_pairings'])),sign))
        error,branch=min(vals);assert error<1e-7
        records.append(dict(source_bank_index=idx,source_face=row['face_id'],support=support,
            anchor=2,anchor_label=3,support_contact_label=label,
            nine_finite_contacts=[[i,label if i==27 else 3 if i==2 else 2]for i in ids],
            free_coordinate=free,kernel_v=enc(v),particular_p=enc(p),
            kernel_norm=enc(vq),particular_norm=enc(pq),orthogonality='p.v=0',
            delta=enc(delta),F_norm_delta=str(normdelta),F_norm_delta_over_RAD=str(normratio),
            square_in_K=False,numerically_matching_branch=branch,
            branch_pairing_max_error=smooth(error),
            normal_recipe='p + branch*sqrt(delta)*v; norm2 exactly, extension-field Coxeter/endpoint checks NOT certified'))
    return records

def run():
    mp.mp.dps=100
    d=json.loads((ROOT/BASE).read_text());meta=json.loads((ROOT/f'results/{BTAG}.json').read_text())
    for doc in[d,meta]:
        for p,h in doc.get('sources',{}).items():assert sha(p)==h,p
    assert sha(BANK)==meta['candidate_bank']['sha256']
    bank=json.loads(gzip.decompress((ROOT/BANK).read_bytes()))
    assert len(bank)==68111 and len(d['normals'])==29
    exactQ=[[dec(x)for x in r]for r in d['form']]
    exactN=[[dec(x)for x in r]for r in d['normals']]
    cov=[mv(exactQ,n)for n in exactN[:25]]
    _,basis=dm(cov).transpose().rref();basis=list(basis);assert len(basis)==10
    Q=mp.matrix([[scalar(x)for x in r]for r in d['form']])
    N=[mp.matrix([scalar(x)for x in r])*mp.sqrt(2/scalar(d['norms'][i]))for i,r in enumerate(d['normals'])]
    M=mp.matrix([[x for x in Q*N[i]]for i in basis]);J=M**-1
    E=mp.matrix([scalar(x)for x in d['strict_escape']['point']]);H=mp.matrix([scalar(x)for x in d['center']])
    transform=np.array([[float(x)for x in v.T*Q*J]for v in[N[28],E,H]])
    V=np.array([r['normalized_wall_pairings']for r in bank]);assert V.shape==(68111,28)
    new=V[:,basis]@transform.T
    finite=-2*np.cos(np.pi/np.arange(2,257))
    errors=np.min(abs(new[:,0,None]-finite[None,:]),axis=1)
    angle_ok=(new[:,0]<-2-1e-7)|(errors<1e-7)
    ok=angle_ok&(new[:,2]<-1e-8)
    cutters=ok&(new[:,1]>1e-7)
    ids=list(map(int,np.where(cutters)[0]));assert int(ok.sum())==67056 and len(ids)==152
    out=[]
    for i in ids:
        r=bank[i]
        keep={k:v for k,v in r.items()if k not in('normalized_wall_pairings','endpoint_other_pairings')}
        out.append(dict(source_bank_index=i,**keep,new_pairings=list(map(smooth,new[i]))))
    out.sort(key=lambda r:(-r['new_pairings'][1],r['source_bank_index']))
    recipes=recipe_check(d,bank,ids)
    paths=[BASE,f'scripts/fresh_h9_fourth_escape_extension_20261015.py',BANK,
           f'results/{BTAG}.json',f'scripts/{BTAG}.py',f'scripts/{TAG}.py']
    return dict(status='PASS_BOUNDED_NUMERICAL_29_WALL_RESCORE_WITH_EXACT_NONSQUARE_RECIPES',
        sources={p:sha(p)for p in paths},field=d['field'],
        numerical_scope='Replays only the saved bounded68111proposal bank; inherited old28 gates are numerical. No all-label/all-parameter cap search and no exact new cap.',
        arithmetic='Standalone100-digit mpmath frame reconstruction then binary64 NumPy rescore; separate SymPy algebraic-field nine-contact identities, no project arithmetic imports.',
        independent_pairing_basis=basis,transform_outputs=['norm2-normalized wall28','current escape','current strict center'],
        tolerances=dict(finite_labels=[2,256],finite_angle_fit=1e-7,dotted_beyond_minus2=1e-7,
                        center_negative=1e-8,escape_positive=1e-7),
        inherited_search=dict(anchor_labels=meta['anchor_labels'],sampled_dihedral_tail_labels=meta['tail_labels']),
        bank_size=len(bank),numerically_compatible_count=int(ok.sum()),numerical_escape_cutters=len(ids),
        accepted_mask_sha256=sha256(np.packbits(ok).tobytes()).hexdigest(),
        cutter_mask_sha256=sha256(np.packbits(cutters).tobytes()).hexdigest(),
        cutters=out,exact_nine_contact_recipes=recipes,
        exact_scope='Nine linear contacts, norm2 quadratic recipes, and nonsquareness in the inherited quartic field are exact. The required quadratic extension is not an obstruction to geometric existence. Full new-wall pair gates, actual endpoints and a coherent finite-volume completion remain unproved.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved bounded29rescore and exact recipe tests')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['numerically_compatible_count'],d['numerical_escape_cutters'])
