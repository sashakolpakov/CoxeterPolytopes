#!/usr/bin/env python3
"""Four-wall noncentral section proposals with two auxiliary orthogonalities.

Numerical scouts select finitely many rational normal pencils. Every real
root of the declared exact angle equations is then isolated. Geometric
candidates require a strict timelike center; three angle repairs are claimed
only when all three actual ridges are certified. Rejection of a geometric
candidate requires an actual forbidden ridge.
No completeness over four-wall supports or coefficient choices is claimed.
"""
from __future__ import annotations
import argparse
from collections import Counter
from itertools import combinations, permutations
import json
from types import SimpleNamespace
import numpy as np
import sympy as sp
import fresh_h20_triple_central_20261013 as source
import fresh_h20_triple_angle_targets_20261013 as old
import fresh_h20_bivariate_nonelliptic_20261013 as ridge
from fresh_h20_global_section_parameters_20261013 import Field, coeffs

ROOT=source.ROOT
OUTPUT=ROOT/'data/enumerations/fresh_h20_fourwall_sections_20261013.json'
t=old.t; P=old.P
SAMPLES=tuple(sp.Rational(x) for x in ('-8','-4','-2','-1','-1/2','-1/4','1/4','1/2','1','2','4','8'))


def supports(G):
    """Literal first representatives by weighted Gram, not claimed orbits."""
    seen={}
    for orbit in json.loads(source.ORBIT_OUTPUT.read_text())['orbits']:
        ids=[v-1 for v in orbit['representative']]
        if np.linalg.eigvalsh(np.asarray(G[np.ix_(ids,ids)],float)).min()>1e-7:continue
        for fourth in range(210):
            if fourth in ids:continue
            js=tuple(sorted((*ids,fourth)));M=G[np.ix_(js,js)]
            key=min(tuple(int(M[a,b]) for a in p for b in p) for p in permutations(range(4)))
            seen.setdefault(key,js)
    return list(seen.values())


def nullspace(rows):
    M=sp.Matrix(rows); ns=M.nullspace()
    if len(ns)!=2:return None
    if any(ns[0][i]==ns[1][i]==0 for i in range(4)):return None
    return tuple(tuple(v) for v in ns)


def lower_support_plane(W,ids,basis):
    E=sp.Matrix([[sum(c*W[i,d] for c,i in zip(v,ids)) for d in range(22)] for v in basis])
    reduced,pivots=E.rref();assert len(pivots)==2
    groups={};zeros=[]
    for i,w in enumerate(W):
        v=sp.Matrix([list(w)])
        for k,pivot in enumerate(pivots):v-=v[pivot]*reduced[k,:]
        first=next((x for x in v if x),None)
        if first is None:zeros.append(i)
        else:groups.setdefault(tuple(x/first for x in v),[]).append(i)
    for group in [zeros,*[js+zeros for js in groups.values()]]:
        for size in (2,3):
            for js in combinations(group,size):
                M=sp.Matrix([list(W[i]) for i in js])
                if M.col_join(E).rank()==M.rank():return tuple(i+1 for i in js)
    return None


def section(Q,W,G,h,ids,basis):
    vectors=[sum((c*W[i] for c,i in zip(v,ids)),np.zeros(22,dtype=object)) for v in basis]
    e0,e1=vectors;e=[P(e0[d]+t*e1[d]) for d in range(22)]
    q=P(e0@Q@e0+2*t*(e0@Q@e1)+t*t*(e1@Q@e1))
    vals=[P(w@Q@e0+t*(w@Q@e1)) for w in W]
    he=P(h@Q@e0+t*(h@Q@e1))
    def gr(a,b):return q*int(G[a,b])-vals[a]*vals[b]
    C=[P(q*int(h[d])-he*e[d]) for d in range(22)]
    return (e,q,vals,gr),C,he


def numerical_ridges(H,cp,pairs):
    for a,b in pairs:
        aa,bb,ab=H[a,a],H[b,b],H[a,b];D=aa*bb-ab*ab
        if min(aa,bb,D)<=1e-8:continue
        A=bb*cp[a]-ab*cp[b];B=aa*cp[b]-ab*cp[a]
        vals=D*cp-A*H[:,a]-B*H[:,b]
        if max(vals)<1e-8*max(1,max(abs(vals))):yield a,b


def repaired_root_scout(sd,target,ids,basis,aux,gf,hp):
    """Find numerical three-angle proposals, never a validity certificate."""
    a,b=target;gr=sd[3];N=gr(a,b)**2;den=gr(a,a)*gr(b,b)
    equations={3:4*N-den,4:2*N-den,5:16*N*N-12*N*den+den*den,6:4*N-3*den}
    bf=np.asarray(basis,float);M=gf[np.ix_(ids,ids)];g=gf[:,ids];hs=hp[list(ids)];best=None
    for m,eq in equations.items():
        if eq.is_zero:continue
        for z in np.roots([float(x) for x in eq.all_coeffs()]):
            if abs(z.imag)>1e-7:continue
            value=float(z.real);c=bf[0]+value*bf[1]
            if min(abs(c))<1e-7:continue
            q=float(c@M@c);he=float(hs@c)
            if q<=1e-7 or abs(he)<1e-7:continue
            v=g@c;cp=q*hp-he*v
            if max(cp)>=-1e-7:continue
            H=gf-np.outer(v,v)/q
            aa,bb,ab=H[a,a],H[b,b],H[a,b]
            if min(aa,bb)<=1e-7 or ab>=0:continue
            if abs(ab*ab/(aa*bb)-np.cos(np.pi/m)**2)>1e-6:continue
            if next(numerical_ridges(H,cp,[target]),None) is None:continue
            rights=[]
            for l in aux:
                pairs=((min(l,o),max(l,o)) for o in range(210) if o!=l and gf[l,o]==0)
                pair=next(numerical_ridges(H,cp,pairs),None)
                if pair is not None:rights.append(pair)
            if len(rights)!=2 or rights[0]==rights[1]:continue
            diag=np.diag(H);denom=np.outer(diag,diag);sq=H*H
            elliptic=(diag[:,None]>1e-8)&(diag[None,:]>1e-8)&(denom-sq>1e-8)
            close=np.zeros((210,210),bool)
            for cos2 in (0,.25,.5,(3+5**.5)/8,.75):close|=abs(sq-cos2*denom)<1e-6
            score=int(np.triu(elliptic&((H>1e-8)|(~close)),1).sum())
            key=(score,m,abs(value),value)
            if best is None or key<best:best=key
    return best


def scout(Q,W,G,h,limit):
    gf=np.asarray(G,float);hp=np.asarray(W@Q@h,float);hh=float(h@Q@h)
    allsupports=supports(G);scouts=[];census=Counter()
    for ids in allsupports:
        patterns={}
        for l in range(210):
            if l in ids:continue
            row=tuple(int(G[l,i]) for i in ids)
            if any(row):patterns.setdefault(row,[]).append(l)
        bestpatterns=sorted(patterns,key=lambda r:(-len(patterns[r]),r))[:8]
        plane_seen=set();support_scouts=[]
        for ra,rb in combinations(bestpatterns,2):
            basis=nullspace((ra,rb))
            if basis is None or basis in plane_seen:continue
            plane_seen.add(basis);census['full_four_coordinate_planes']+=1
            bf=np.asarray(basis,float);M=gf[np.ix_(ids,ids)];g=gf[:,ids];hs=hp[list(ids)]
            if np.max(abs(bf@hs))<1e-10:continue
            best=None
            for value in SAMPLES:
                c=bf[0]+float(value)*bf[1]
                if min(abs(c))<1e-10:continue
                q=float(c@M@c);he=float(hs@c)
                if q<=1e-8 or abs(he)<1e-9:continue
                v=g@c;cp=q*hp-he*v
                if max(cp)>=-1e-8:continue
                H=gf-np.outer(v,v)/q
                diag=np.diag(H);den=np.outer(diag,diag);N=H*H
                elliptic=(diag[:,None]>1e-8)&(diag[None,:]>1e-8)&(den-N>1e-8)
                close=np.zeros((210,210),bool)
                for cos2 in (0,.25,.5,(3+5**.5)/8,.75):close|=abs(N-cos2*den)<1e-7
                bad=np.triu(elliptic&((H>1e-8)|(~close)),1)
                score=int(bad.sum())
                proposal=(score,-int(np.count_nonzero(abs(v)<1e-9)),abs(he)/q,value,H,cp,bad)
                if best is None or proposal[:4]<best[:4]:best=proposal
            if best is not None:
                aux=(patterns[ra][0],patterns[rb][0]);support_scouts.append((best,ids,basis,aux))
        # Cap the expensive actual-ridge scouts, retaining two pencils per support.
        for best,ids,basis,aux in sorted(support_scouts,key=lambda p:p[0][:4])[:2]:
            score,minuszero,offset,value,H,cp,bad=best
            right=[]
            for l in aux:
                pairs=((min(l,o),max(l,o)) for o in range(210) if o!=l and G[l,o]==0)
                pair=next(numerical_ridges(H,cp,pairs),None)
                if pair is None:break
                right.append(pair)
            if len(right)!=2 or right[0]==right[1]:continue
            sd,C,he=section(Q,W,G,h,ids,basis)
            target=None;target_scout=None
            badpairs=[tuple(map(int,p)) for p in np.argwhere(bad)]
            badpairs.sort(key=lambda p:(H[p]>0,abs(H[p]**2/(H[p[0],p[0]]*H[p[1],p[1]])-.25),p))
            tries=0
            for a,b in numerical_ridges(H,cp,badpairs):
                if H[a,b]>=-1e-9:continue
                ratio=sp.cancel(sd[3](a,b).as_expr()**2/(sd[3](a,a)*sd[3](b,b)).as_expr())
                if t not in ratio.free_symbols:continue
                tries+=1;candidate=repaired_root_scout(sd,(a,b),ids,basis,aux,gf,hp)
                if candidate is not None and (target_scout is None or candidate<target_scout):
                    target=(a,b);target_scout=candidate
                if tries==8:break
            if target is None:continue
            census['pencils_with_numerical_nonright_three_angle_repair']+=1
            scouts.append((target_scout[0],minuszero,offset,ids,basis,aux,value,target,right))
        print('support',tuple(i+1 for i in ids),'pencils',len(scouts),flush=True)
    # One pencil per literal support, ranked by complete-pair heuristic only.
    chosen=[];used=set();reduced_planes=[]
    for row in sorted(scouts,key=lambda p:(p[:4],p[4])):
        if row[3] in used:continue
        reduction=lower_support_plane(W,row[3],row[4])
        if reduction is not None:
            reduced_planes.append({'support':[i+1 for i in row[3]],'spanning_walls':reduction});continue
        chosen.append(row);used.add(row[3])
        if len(chosen)==limit:break
    return chosen,{'support_gram_representatives':len(allsupports),**census,
                   'scout_pencils_with_two_right_ridges':len(scouts),'selected':len(chosen),
                   'discarded_planes_in_three_wall_spans':reduced_planes,
                   'sample_parameters':[str(x) for x in SAMPLES],
                   'note':'Numerical scouts and Gram-type deduplication select proposals only, not geometric exclusions or a symmetry-orbit classification.'}


def alternate_center(Q,W,h,sd,F):
    """LP proposes an interior; rational projection and exact signs certify it."""
    from scipy.optimize import linprog
    value=float((F.lo+F.hi)/2);e,q,v,gr=sd
    ef=np.asarray([float(p.eval(value)) for p in e]);Qf=np.asarray(Q,float);Wf=np.asarray(W,float);hf=np.asarray(h,float)
    lp=linprog(np.r_[np.zeros(22),-1.],A_ub=np.c_[Wf@Qf,np.ones(210)],b_ub=np.zeros(210),
               A_eq=np.c_[np.asarray([hf@Qf,ef@Qf]),np.zeros(2)],b_eq=[-1.,0.],
               bounds=[(None,None)]*22+[(0,None)],method='highs')
    if not lp.success or lp.x[-1]<1e-7:return None
    for denominator in (100,10000,1000000):
        raw=[P(sp.Rational(str(float(x))).limit_denominator(denominator)) for x in lp.x[:-1]]
        pairing=old.pdot(Q,raw,e,F);C=[F.rem(q*x-pairing*y) for x,y in zip(raw,e)]
        if F.sign(old.pdot(Q,C,C,F))>=0 or F.sign(old.pdot(Q,C,[P(int(x)) for x in h],F))>=0:continue
        if all(F.sign(old.pdot(Q,C,[P(int(x)) for x in w],F))<0 for w in W):return C
    return None


def run_pencil(Q,W,G,h,proposal):
    score,minuszero,offset,ids,basis,aux,value,target,scout_right=proposal
    sd,rawC,he=section(Q,W,G,h,ids,basis);e,q,v,gr=sd;a,b=target
    N=gr(a,b)**2;den=gr(a,a)*gr(b,b)
    equations={2:N,3:4*N-den,4:2*N-den,5:16*N*N-12*N*den+den*den,6:4*N-3*den}
    factors={}
    for eq in equations.values():
        assert not eq.is_zero
        for f,power in sp.factor_list(eq)[1]:factors[tuple(coeffs(f.monic()))]=f.monic()
    rows=[]
    for key,poly in sorted(factors.items()):
        for (lo,hi),mult in poly.intervals(eps=sp.Rational(1,10**18)):
            if lo==hi:lo-=sp.Rational(1,10**20);hi+=sp.Rational(1,10**20)
            F=Field(poly,lo,hi);row={'minimal_polynomial':list(key),'isolating_interval':[str(lo),str(hi)],'approximate_parameter':float((lo+hi)/2)}
            if any(F.sign(P(basis[0][i]+t*basis[1][i]))==0 for i in range(4)):
                row['status']='LOWER_SUPPORT_PARAMETER';rows.append(row);continue
            if F.sign(he)==0:
                row['status']='CENTRAL_PARAMETER';rows.append(row);continue
            if F.sign(q)<=0:
                row['status']='NO_SPACELIKE_NORMAL';rows.append(row);continue
            if min(F.sign(gr(a,a)),F.sign(gr(b,b)))<=0 or F.sign(gr(a,b))>0:
                row['status']='TARGET_NOT_ACUTE_ELLIPTIC_PAIR';rows.append(row);continue
            labels=[m for m,eq in equations.items() if F.sign(eq)==0 and (m!=5 or F.sign(2*N-den)>0)]
            if not labels:
                row['status']='CONJUGATE_ANGLE_BRANCH';rows.append(row);continue
            C=[F.rem(p) for p in rawC];cp=[old.pdot(Q,C,[P(int(z)) for z in w],F) for w in W];method='projected_original_center'
            if any(F.sign(p)>=0 for p in cp):
                C=alternate_center(Q,W,h,sd,F)
                if C is None:
                    row.update(status='UNRESOLVED_STRICT_CENTER_SCOUT',target_labels=labels);rows.append(row);continue
                cp=[old.pdot(Q,C,[P(int(z)) for z in w],F) for w in W];method='rational_LP_point_projected_exactly'
            target_proof=ridge.actual_ridge(Q,W,h,sd,C,F,a,b,cp)
            rights=[]
            for l in aux:
                found=None
                for other in range(210):
                    if other==l or G[l,other]!=0:continue
                    found=ridge.actual_ridge(Q,W,h,sd,C,F,l,other,cp)
                    if found is not None:break
                if found is not None:rights.append(found)
            bad=None
            for x,y in old.candidates(Q,W,h,sd,C,F):
                bad=old.witness(Q,W,h,sd,C,F,x,y,cp)
                if bad is not None:break
            repaired=target_proof is not None and len(rights)==2
            status=('NON_COXETER_ACTUAL_RIDGE' if repaired else 'NON_COXETER_WITH_INCOMPLETE_REPAIR_CERTIFICATES') if bad else 'OPEN_COXETER_CANDIDATE'
            row.update(isolating_interval=[str(F.lo),str(F.hi)],target_labels=labels,strict_interior=[coeffs(p) for p in C],interior_method=method,
                       repaired_target_ridge=target_proof,repaired_right_ridges=rights,actual_forbidden_ridge=bad,
                       status=status)
            rows.append(row)
            print('candidate',tuple(i+1 for i in ids),labels,row['approximate_parameter'],row['status'],flush=True)
    return {'support_walls_one_based':[i+1 for i in ids],
            'support_gram':[[int(G[i,j]) for j in ids] for i in ids],
            'normal_coefficient_basis':[[str(x) for x in v] for v in basis],
            'auxiliary_orthogonal_walls_one_based':[i+1 for i in aux],
            'center_pairing_polynomial':coeffs(he),'scout_parameter':str(value),'scout_raw_pair_conflicts':score,
            'scout_score_note':'The seed parameter selects a target pair; the score counts raw pair conflicts at the best numerical repaired-angle root, not necessarily at that seed. Neither count is an actual-ridge census.',
            'identically_orthogonal_wall_count':sum(p.is_zero for p in v),
            'target_pair_one_based':[a+1,b+1],'target_labels':[2,3,4,5,6],
            'target_equations':{str(m):coeffs(p) for m,p in equations.items()},'candidates':rows}


def build(limit,replay=False):
    Q,W,G,h=source.base.load()
    if replay:
        prior=json.loads(OUTPUT.read_text());assert prior['complete'];scouting=prior['scouting'];selected=[]
        for r in prior['pencils']:
            selected.append((r['scout_raw_pair_conflicts'],0,0,tuple(i-1 for i in r['support_walls_one_based']),
                             tuple(tuple(map(sp.Rational,v)) for v in r['normal_coefficient_basis']),
                             tuple(i-1 for i in r['auxiliary_orthogonal_walls_one_based']),sp.Rational(r['scout_parameter']),
                             tuple(i-1 for i in r['target_pair_one_based']),[]))
    else:selected,scouting=scout(Q,W,G,h,limit)
    rows=[]
    for proposal in selected:
        rows.append(run_pencil(Q,W,G,h,proposal))
        OUTPUT.write_text(json.dumps({'complete':False,'scouting':scouting,'pencils':rows},indent=2)+'\n')
    out={'schema':'noncentral-four-wall-three-angle-sections-v1','complete':True,
         'source':str(source.base.SOURCE.relative_to(ROOT)),'scouting':scouting,
         'pencil_count':len(rows),'pencils':rows,'statuses':dict(Counter(c['status'] for p in rows for c in p['candidates'])),
         'scope':'Finite selection of full four-coordinate noncentral pencils, two independent auxiliary orthogonalities and one exact angle equation. Complete real-root census of the declared equations only; no global four-wall or support-orbit exclusion. Every geometric candidate has a strict timelike center and an actual forbidden ridge. Three repaired actual ridges are asserted only when all three certificates are present.'}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print('DONE',json.dumps(out['statuses']),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--limit',type=int,default=12)
    ap.add_argument('--replay-pencils',action='store_true',help='Keep the archived literal pencils and repeat their exact root/geometry calculations.')
    a=ap.parse_args();build(a.limit,a.replay_pencils)
