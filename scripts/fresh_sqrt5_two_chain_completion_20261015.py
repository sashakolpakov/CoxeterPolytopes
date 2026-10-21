#!/usr/bin/env python3
"""Exact all-length necessity for two disjoint simply-laced white chains.

Only each chain's first node touches a nonsingular unpaid finite core.
Infinite Diophantine families, if present, are retained explicitly, not sampled.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import gcd, isqrt, lcm
from pathlib import Path
import fresh_sqrt5_quad_contacts_20261015 as C
import fresh_sqrt5_quad_filters_20261015 as A
from fresh_sqrt5_single_chain_completion_20261015 import inverse

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_two_chain_completion_20261015'
K = A.K


def chains(core, contacts, lengths):
    n = core['nodes']; total = n + sum(lengths)
    old = C.edges(core)
    g = {p: old.get(p, 2) for p in combinations(range(total), 2)}
    start = n
    for contact, length in zip(contacts, lengths):
        assert length >= 1
        for i, m in enumerate(contact):
            g[i, start] = m
        for j in range(start, start + length - 1):
            g[j, j + 1] = 3
        start += length
    return dict(nodes=total, white=list(range(4, total)),
                labels=[g[p] for p in combinations(range(total), 2)])


def polynomial(coeff, ell, emm):
    a, b, c = coeff
    return a*ell*emm + b*ell + c*emm + 1


def positive_integer(x):
    return not x.b and x.a.denominator == 1 and x.a > 0


def integer_roots(q):
    """All positive integer roots of a nonzero rational quadratic."""
    aa, bb, cc = q
    assert aa or bb or cc
    if not aa:
        roots = [] if not bb else [-cc/bb]
    else:
        delta = bb*bb - 4*aa*cc
        if delta < 0:
            return []
        u, v = isqrt(delta.numerator), isqrt(delta.denominator)
        if u*u != delta.numerator or v*v != delta.denominator:
            return []
        s = F(u, v)
        roots = [(-bb-s)/(2*aa), (-bb+s)/(2*aa)]
    return sorted({int(x) for x in roots if x > 0 and x.denominator == 1})


def extended_gcd(a, b):
    oldr, r, olds, s, oldt, t = abs(a), abs(b), 1, 0, 0, 1
    while r:
        q = oldr//r
        oldr, r = r, oldr-q*r
        olds, s = s, olds-q*s
        oldt, t = t, oldt-q*t
    return oldr, olds*(1 if a >= 0 else -1), oldt*(1 if b >= 0 else -1)


def rational_bilinear(coeff):
    assert all(not v.b for v in coeff)
    values = [v.a for v in coeff] + [F(1)]
    den = lcm(*(v.denominator for v in values))
    aa, bb, cc, dd = [int(v*den) for v in values]
    gg = gcd(gcd(abs(aa), abs(bb)), gcd(abs(cc), abs(dd)))
    aa, bb, cc, dd = [v//gg for v in (aa, bb, cc, dd)]
    out = dict(method='RATIONAL_BILINEAR', integer_coefficients=[aa, bb, cc, dd],
               finite_pairs=[], families=[])
    if aa:
        # (a*l+c)(a*m+b)=bc-ad. Both divisor signs are needed.
        nn = bb*cc-aa*dd
        out['factor_product'] = nn
        if not nn:
            if -cc % aa == 0 and -cc//aa > 0:
                out['families'].append(dict(kind='FIXED_FIRST', first=-cc//aa, second_min=1))
            if -bb % aa == 0 and -bb//aa > 0:
                out['families'].append(dict(kind='FIXED_SECOND', second=-bb//aa, first_min=1))
        else:
            # Exact factorization, with the resulting divisors checked below.
            from sympy import divisors
            for xx0 in divisors(abs(nn)):
                for xx in (int(xx0), -int(xx0)):
                    assert nn % xx == 0
                    yy = nn//xx
                    if (xx-cc) % aa == 0 and (yy-bb) % aa == 0:
                        ell, emm = (xx-cc)//aa, (yy-bb)//aa
                        if ell > 0 and emm > 0:
                            out['finite_pairs'].append([ell, emm])
    elif not bb and not cc:
        assert dd != 0
    elif not bb or not cc:
        variable, value = (cc, 'FIXED_SECOND') if not bb else (bb, 'FIXED_FIRST')
        if -dd % variable == 0 and -dd//variable > 0:
            out['families'].append(dict(kind=value, **(
                dict(second=-dd//variable, first_min=1) if not bb else
                dict(first=-dd//variable, second_min=1))))
    else:
        gg, u, v = extended_gcd(bb, cc)
        assert bb*u+cc*v == gg
        if (-dd) % gg == 0:
            ell0, emm0 = u*(-dd//gg), v*(-dd//gg)
            sl, sm = cc//gg, -bb//gg
            lo = hi = None
            for x0, step in ((ell0, sl), (emm0, sm)):
                if step > 0:
                    bound = -((x0-1)//step)
                    lo = bound if lo is None else max(lo, bound)
                else:
                    bound = (1-x0)//step
                    hi = bound if hi is None else min(hi, bound)
            if lo is None or hi is None:
                out['families'].append(dict(kind='AFFINE_PARAMETER', first_start=ell0,
                    second_start=emm0, first_step=sl, second_step=sm,
                    parameter_min=lo, parameter_max=hi))
            elif lo <= hi:
                out['finite_pairs'] = [[ell0+sl*t, emm0+sm*t] for t in range(lo, hi+1)]
    out['finite_pairs'] = sorted(set(map(tuple, out['finite_pairs'])))
    out['finite_pairs'] = [list(p) for p in out['finite_pairs']]
    return out


def solve_lengths(coeff):
    aa, bb, cc = coeff
    if not any(v.b for v in coeff):
        out = rational_bilinear(coeff)
    else:
        q = (aa.a*bb.b-aa.b*bb.a,
             cc.a*bb.b-aa.b-cc.b*bb.a, -cc.b)
        out = dict(method='NONRATIONAL_ELIMINATION', elimination=[str(v) for v in q],
                   finite_pairs=[], families=[])
        if any(q):
            roots = integer_roots(q)
            out['positive_integer_elimination_roots'] = roots
            for ell in roots:
                denominator, numerator = aa*ell+cc, bb*ell+1
                if denominator:
                    emm = -numerator/denominator
                    if positive_integer(emm):
                        out['finite_pairs'].append([ell, int(emm.a)])
                elif not numerator:
                    out['families'].append(dict(kind='FIXED_FIRST', first=ell, second_min=1))
        else:
            # With constant term1, the only nonrational degeneration is
            # (b*l+1)(c0*m+1). This includes c0=0 (no integer solution).
            assert not cc.b and bb.b and aa == cc.a*bb
            out['method'] = 'NONRATIONAL_FACTORIZATION'
            if cc.a:
                emm = -1/cc.a
                if emm.denominator == 1 and emm > 0:
                    out['families'].append(dict(kind='FIXED_SECOND', second=int(emm), first_min=1))
    assert all(not polynomial(coeff, *p) for p in out['finite_pairs'])
    for fam in out['families']:
        # The exact coefficient identities prove all lengths, not just samples.
        if fam['kind'] == 'FIXED_FIRST':
            ell = fam['first']; assert not aa*ell+cc and not bb*ell+1
        elif fam['kind'] == 'FIXED_SECOND':
            emm = fam['second']; assert not aa*emm+bb and not cc*emm+1
        else:
            ell, emm = fam['first_start'], fam['second_start']
            sl, sm = fam['first_step'], fam['second_step']
            assert not aa and not polynomial(coeff, ell, emm) and not bb*sl+cc*sm
    return out


def solver_checks():
    """Exercise each exact Diophantine branch against a bounded direct oracle."""
    examples = [
        [K(-1), K(), K()], [K(1), K(-2), K(-2)],
        [K(1), K(-1), K(-1)], [K(), K(-1), K(1)],
        [K(), K(F(-1,8)), K(F(-1,12))], [K(), K(), K(-1)],
        [K(), K(), K()], [K(), K(-1), K()],
        [K(0,1), K(0,-1), K(-1)], [K(), K(0,1), K()],
        [K(1,1), K(-2), K(0,-1)], [K(1,1), K(0,-1), K(-1)],
    ]
    def includes(out, ell, emm):
        if [ell, emm] in out['finite_pairs']:return True
        for r in out['families']:
            if r['kind']=='FIXED_FIRST' and ell==r['first']:return True
            if r['kind']=='FIXED_SECOND' and emm==r['second']:return True
            if r['kind']=='AFFINE_PARAMETER':
                if (ell-r['first_start']) % r['first_step']:continue
                t=(ell-r['first_start'])//r['first_step']
                if (r['parameter_min'] is None or t>=r['parameter_min']) and (
                    r['parameter_max'] is None or t<=r['parameter_max']) and emm==r['second_start']+r['second_step']*t:return True
        return False
    for co in examples:
        out=solve_lengths(co)
        for ell in range(1,31):
            for emm in range(1,31):
                assert includes(out,ell,emm)==(not polynomial(co,ell,emm)),(co,out,ell,emm)
    return len(examples)*900


def completed_star(core, contact, lengths):
    row=chains(core,contact,lengths);out=dict(lengths=lengths, nodes=row['nodes'])
    g=C.edges(row)
    if not C.E.corners_finite(row['white'],g):
        out['status']='REJECT_NONFINITE_CORNER';return out
    color=A.coloring(row);assert color['passes']
    G=A.scaled_matrix(row,color['colors'])
    physical=A.inertia(G)
    conjugate=A.inertia([[x.conj() for x in r] for r in G])
    out.update(physical_inertia=physical,conjugate_inertia=conjugate)
    if physical != [row['nodes']-2,1,1] or conjugate != [row['nodes']-1,0,1]:
        out['status']='REJECT_FULL_STAR_INERTIA';return out
    aff=C.affine_witness(row)
    if aff:out.update(status='REJECT_AFFINE_WITNESS',affine_witness=aff)
    elif not C.E.no_short(row['white'],g):out['status']='PAID_SHORT_PATH'
    elif C.h4_payments(row):out['status']='PAID_H4'
    else:out.update(status='UNPAID_PRINCIPAL_STAR',row=row)
    return out


def bound_family_by_corner(core, contact, family):
    """A length-two non-PD principal corner excludes every longer length.

    This is a certified lower-bound failure, not a length cutoff. Any family
    not admitting this witness is retained open without inference.
    """
    kind=family['kind']
    if kind not in ('FIXED_FIRST','FIXED_SECOND'):
        return dict(status='OPEN',family=family)
    lens=[family['first'],2] if kind=='FIXED_FIRST' else [2,family['second']]
    row=chains(core,contact,lens)
    cc=A.coloring(row);assert cc['passes']
    M=A.scaled_matrix(row,cc['colors'])
    for a,b in product((0,2),(1,3)):
        ids=row['white']+[a,b];G=A.sub(M,ids);ii=A.inertia(G)
        if ii!=[len(ids),0,0]:
            below=[family['first'],1] if kind=='FIXED_FIRST' else [1,family['second']]
            out=dict(status='ALL_FREE_LENGTHS_AT_LEAST_TWO_REJECTED',family=family,
                     failed_lengths=lens,marked_corner=[a,b],corner_nodes=ids,
                     corner_gram=A.encode(G),corner_inertia=ii,
                     remaining_positive_pairs=[below])
            kernels=A.nullspace(G)
            if kernels:
                out['corner_kernel']=[[x.enc() for x in v] for v in kernels]
                assert all(not A.dot(r,v) for v in kernels for r in G)
            return out
    return dict(status='OPEN',family=family)


def run():
    checks=solver_checks()
    contacts=json.loads((ROOT/'results/fresh_sqrt5_quad_contacts_20261015.json').read_text())
    unions=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    selected=[r for r in contacts['records'] if r['status']=='UNPAID_FINITE_CONTACT_DOMAIN'
              and not r['singular'] and r['accepted']]
    records=[];counts=Counter();completed=Counter();lengths=Counter();methods=Counter()
    profile_count=0;family_count=0;open_families=0;bounded_families=0
    for ri,rec in enumerate(selected):
        core=unions['labelled_union_supports'][rec['union_index']]
        color=A.coloring(core);assert color['passes']
        H=A.scaled_matrix(core,color['colors']);Hi=inverse(H);n=core['nodes']
        assert A.inertia(H)==[n-1,1,0]
        assert A.inertia([[x.conj() for x in r] for r in H])==[n,0,0]
        profiles=[(i,t['contacts']) for i,t in enumerate(rec['trials'])
                  if t['status']=='PASS_EXTRA_WHITE_NECESSARY_TESTS']
        assert [v for _,v in profiles]==rec['accepted']
        profile_count+=len(profiles)
        for (ti,v),(tj,w) in combinations_with_replacement(profiles,2):
            out=dict(source_class_index=rec['class_index'],union_index=rec['union_index'],
                     trial_indices=[ti,tj],contacts=[v,w])
            head=chains(core,[v,w],[1,1]);g=C.edges(head)
            if not C.E.corners_finite(head['white'],g):out['status']='HEAD_REJECT_NONFINITE_CORNER'
            elif not C.E.no_short(head['white'],g):out['status']='HEAD_PAID_SHORT_PATH'
            elif C.h4_payments(head):out['status']='HEAD_PAID_H4'
            else:
                cc=A.coloring(head)
                if not cc['passes']:out['status']='HEAD_REJECT_CYCLE_FIELD'
                else:
                    assert cc['colors'][:n]==color['colors']
                    G=A.scaled_matrix(head,cc['colors'])
                    tests=A.conditions(head,G);out['head_conditions']=tests
                    if not tests['physical']['passes']:out['status']='HEAD_REJECT_TWO_NEGATIVE'
                    elif not tests['conjugate']['passes']:out['status']='HEAD_REJECT_CONJUGATE_PSD'
                    else:
                        aff=C.affine_witness(head)
                        if aff:out.update(status='HEAD_REJECT_AFFINE',affine_witness=aff)
                        else:
                            vv=[G[i][n] for i in range(n)];ww=[G[i][n+1] for i in range(n)]
                            d,e=G[n][n],G[n+1][n+1]
                            hv=[A.dot(r,vv) for r in Hi];hw=[A.dot(r,ww) for r in Hi]
                            c11,c22,c12=A.dot(vv,hv),A.dot(ww,hw),A.dot(vv,hw)
                            assert c12==A.dot(ww,hv)
                            aa,bb,ff=2*c11/d,2*c22/e,4*(c11*c22-c12*c12)/(d*e)
                            coeff=[1-aa-bb+ff,1-aa,1-bb]
                            for ell,emm in product((1,2),repeat=2):
                                residual_det=(d*F(ell+1,2*ell)-c11)*(e*F(emm+1,2*emm)-c22)-c12*c12
                                assert 4*ell*emm*residual_det/(d*e)==polynomial(coeff,ell,emm)
                            out['schur_data']=dict(first_norms=[d.enc(),e.enc()],
                                costs=[c11.enc(),c22.enc(),c12.enc()],a=aa.enc(),b=bb.enc(),f=ff.enc())
                            out['polynomial_coefficients']=[x.enc() for x in coeff]
                            solution=solve_lengths(coeff);out['length_solutions']=solution
                            methods[solution['method']]+=1
                            out['family_corner_bounds']=[bound_family_by_corner(core,[v,w],fam)
                                                        for fam in solution['families']]
                            finite_pairs=[list(p) for p in solution['finite_pairs']]
                            for bound in out['family_corner_bounds']:
                                if bound['status']=='OPEN':open_families+=1
                                else:
                                    bounded_families+=1
                                    finite_pairs.extend(bound['remaining_positive_pairs'])
                            finite_pairs=[list(p) for p in sorted(set(map(tuple,finite_pairs)))]
                            out['completed_stars']=[]
                            for pair in finite_pairs:
                                star=completed_star(core,[v,w],pair)
                                out['completed_stars'].append(star)
                                completed[star['status']]+=1
                                if star['status']=='UNPAID_PRINCIPAL_STAR':lengths[tuple(pair)]+=1
                            family_count+=len(solution['families'])
                            if any(b['status']=='OPEN' for b in out['family_corner_bounds']):out['status']='OPEN_INFINITE_RANK_FAMILY'
                            elif any(s['status']=='UNPAID_PRINCIPAL_STAR' for s in out['completed_stars']):out['status']='FINITE_UNPAID_RANK_COMPLETIONS'
                            elif finite_pairs:out['status']='ALL_FINITE_LENGTHS_FILTERED'
                            else:out['status']='NO_POSITIVE_INTEGER_LENGTH_PAIR'
            counts[out['status']]+=1;records.append(out)
        if (ri+1)%10==0:print('cores',ri+1,'pairs',len(records),dict(counts),flush=True)
    assert len(selected)==130 and profile_count==458 and len(records)==1455
    assert sum(methods.values())==239 and family_count==bounded_families==5 and open_families==0
    assert completed==Counter(UNPAID_PRINCIPAL_STAR=112)
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_quad_contacts_20261015.py',
           'scripts/fresh_sqrt5_quad_filters_20261015.py','scripts/fresh_sqrt5_quad_exhaustion_20261015.py',
           'scripts/fresh_sqrt5_single_chain_completion_20261015.py',
           'results/fresh_sqrt5_quad_contacts_20261015.json',
           'results/fresh_sqrt5_quad_contacts_20261015_review.json',
           'results/fresh_sqrt5_quad_exhaustion_20261015.json',
           'results/fresh_sqrt5_single_chain_completion_20261015.json',
           'results/fresh_sqrt5_single_chain_completion_20261015_review.json']
    return dict(status='PASS_EXACT_ALL_LENGTH_TWO_WHITE_CHAIN_NECESSITY',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        source_cores=len(selected),source_profiles=profile_count,unordered_profile_pairs=len(records),
        solver_branch_checks=checks,counts=dict(sorted(counts.items())),
        methods=dict(sorted(methods.items())),completed_star_counts=dict(sorted(completed.items())),
        surviving_length_pairs={','.join(map(str,k)):v for k,v in sorted(lengths.items())},
        determinant_identity_checks=4*sum(methods.values()),
        infinite_rank_families=family_count,corner_bounded_rank_families=bounded_families,
        unresolved_infinite_rank_families=open_families,records=records,
        determinant_identity='For chain blocks(d/2)Cartan(A_l),(e/2)Cartan(A_m), only first nodes touch H. With c11=v^T H^-1v,c22=w^T H^-1w,c12=v^T H^-1w, a=2c11/d,b=2c22/e,f=4(c11c22-c12^2)/(de), det(full)/(det(H)det(chain1)det(chain2))=P/((l+1)(m+1)), P=(1-a-b+f)lm+(1-a)l+(1-b)m+1. Equality is in Q(phi), hence both embeddings simultaneously.',
        scope='Complete all-length rank necessity for pairs with replacement of the458 accepted first-contact profiles on the130 old nonsingular unpaid finite cores, with two disjoint white simply-laced chains, only each first node touching the core and no contacts between chains. Orthogonal spectator components do not affect rank. Explicit infinite rank families remain open if present; branching, multiple touching nodes in one component, extra contacts, dotted cores, paid alternatives and arbitrary completions are outside this classification. Surviving finite principal stars are not claimed to be realized compact quadrilateral faces or polytopes. Affine witness absence is not an exhaustive no-affine certificate.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--solver-check',action='store_true');args=ap.parse_args()
    if args.solver_check:print('PASS_SOLVER_BRANCH_CHECKS',solver_checks())
    else:
        out=run();path=ROOT/'results'/f'{TAG}.json'
        if args.verify:assert out==json.loads(path.read_text())
        else:path.write_text(json.dumps(out,indent=2)+'\n')
        print(out['status'],out['counts'],out['completed_star_counts'],
              'rank_families',out['infinite_rank_families'],
              'unresolved',out['unresolved_infinite_rank_families'])
