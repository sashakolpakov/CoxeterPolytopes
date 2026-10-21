#!/usr/bin/env python3
"""Exact continuous incompatibility of two cavity endpoint families.

This is a scoped family gate, not a chamber noncompletion theorem.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_s6_release34_exact_20261017 as X

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_release34_cavity_gate_20261017'
INPUT = ROOT/'results/fresh_nonarithmetic_release34_s6_boundary_20261017.json'
P, U, I = X.P, X.U, X.I


def dec(xs):
    return [I(F(a), F(b)) for a, b in xs]


def dot(a, b):
    return U.rnd(sum((x*y for x, y in zip(a, b)), I()))


def mv(g, x):
    return [dot(row, x) for row in g]


def combine(c, rows):
    return [dot(c, col) for col in zip(*rows)]


def coarse(x):
    d = 10**10
    lo = (x.l*d).__floor__()
    hi = (x.h*d).__ceil__()
    return I(F(lo, d), F(hi, d))


def reduced():
    z = json.loads(INPUT.read_text())['after65']
    g = [dec(row) for row in z['residual_gram']]
    inv, pivots = P.C53.inverse(g)
    cov = {r['wall']: dec(r['covector']) for r in z['projections']}
    ns = {i:mv(inv, row) for i, row in cov.items()}
    h = dec(z['projected_common_center'])
    hp = mv(g, h)
    old_labels=X.D.E.labels()
    assert [old_labels[7][i]for i in X.S]==[2,2,2,2,2,3]
    assert all(old_labels[i][j]==2 for i in[8,12,47]for j in X.S)
    assert old_labels[5][6]==3
    assert all(old_labels[i][j]==2 for i,j in[(7,12),(7,47),(8,12)])
    # A2 inverse diagonal2/3 gives projected old7 square2-2/3=4/3.
    # Replacement52 is right to S6+8 by its bound parent recipe.
    def pair(a, b):
        return dot(a, mv(g, b))
    frames, data, selfs, blockers = [], [], [], []
    for ids, blocker in [([7, 12, 47], 8), ([8, 12, 52], 7)]:
        # The central old side is orthogonal to both remaining sides.
        q = pair(ns[ids[0]], ns[ids[0]])
        diagonal = I(F(4, 3)) if ids[0] == 7 else I(2)
        assert q.l <= diagonal.l <= q.h
        d = pair(ns[ids[1]], ns[ids[2]])
        assert d.h < -2
        gram = [[diagonal,I(),I()],[I(),I(2),d],[I(),d,I(2)]]
        gi, gp = P.C53.inverse(gram)
        ps = [combine([-v for v in gi[j]], [ns[i] for i in ids]) for j in [1,2]]
        # Fixed coordinate gives an independently checkable kernel vector;
        # orthogonality and the unit square are algebraic identities.
        matrix = [cov[i] for i in ids]+[[I(1),I(),I(),I()]]
        ki, kp = P.C53.inverse(matrix)
        k = mv(ki, [I(),I(),I(),I(1)])
        kq = pair(k,k)
        assert kq.l > 0
        k = [U.rnd(v/U.sqrt(kq)) for v in k]
        kh = dot(k,hp)
        assert kh.l*kh.h > 0
        if kh.l > 0:
            k = [-v for v in k]
        frame = ps+[k]
        frames.append(frame)
        selfs.append([[coarse(gi[i][j]) for j in [1,2]] for i in [1,2]])
        blockers.append([coarse(dot(v,cov[blocker])) for v in frame])
        data.append(dict(old_sides=ids, blocker=blocker, gram=[P.enc(r) for r in gram],
            inverse_pivots=gp, kernel_inverse_pivots=kp, kernel_raw_square=kq.enc(),
            reduced_vectors=[P.enc(r) for r in frame], blocker_products=P.enc(blockers[-1]),
            particular_gram=[P.enc(r) for r in selfs[-1]]))
    cross = [[coarse(pair(a,b)) for b in frames[1]] for a in frames[0]]
    return selfs, cross, blockers, dict(input_inverse_pivots=pivots, frames=data, cross_gram=[P.enc(r) for r in cross])


def radical(q, x):
    r = U.sqrt(2-dot(x,mv(q,x)))
    return r, [-v/r for v in mv(q,x)]


def derivative_bounds(q, c, box):
    a, b = box[:2], box[2:]
    ra, da = radical(q[0],a)
    rb, db = radical(q[1],b)
    aa, bb = a+[ra], b+[rb]
    cb = mv(c,bb)
    ac = [dot(aa,col) for col in zip(*c)]
    return [U.rnd(cb[i]+da[i]*cb[2]) for i in range(2)]+[U.rnd(ac[i]+db[i]*ac[2]) for i in range(2)]


def prove(q, c, blockers):
    start = [I(F(7,5),2),I(F(7,5),2),I(F(7,5),2),I(1,2)]
    stack, leaves = [start], []
    while stack:
        box = stack.pop()
        bounds = derivative_bounds(q,c,box)
        if all(v.h < 0 for v in bounds):
            leaves.append(dict(box=P.enc(box), derivatives=P.enc(bounds)))
            continue
        j = max(range(4),key=lambda i:box[i].h-box[i].l)
        assert box[j].h-box[j].l > F(1,2**16), ('unproved derivative',P.enc(box),P.enc(bounds))
        mid = (box[j].l+box[j].h)/2
        for interval in [I(box[j].l,mid),I(mid,box[j].h)]:
            nxt=box.copy();nxt[j]=interval;stack.append(nxt)
    endpoint=[I(2),I(2)]
    ra,_=radical(q[0],endpoint);rb,_=radical(q[1],endpoint)
    value=dot(endpoint+[ra],mv(c,endpoint+[rb]))
    assert value.l > F(13,100)
    branch_proofs=[]
    for j in range(2):
        domains=[start[:2],start[2:]][j]
        queue=[domains];out=[]
        while queue:
            box=queue.pop();r,_=radical(q[j],box)
            # Opposite sign to the coherent branch is obtuse to old8/7.
            product=dot(box+[-r],blockers[j])
            if product.l > 0:
                out.append(dict(box=P.enc(box),product=product.enc()));continue
            k=max(range(2),key=lambda i:box[i].h-box[i].l)
            assert box[k].h-box[k].l > F(1,2**16)
            mid=(box[k].l+box[k].h)/2
            for interval in [I(box[k].l,mid),I(mid,box[k].h)]:
                nxt=box.copy();nxt[k]=interval;queue.append(nxt)
        branch_proofs.append(out)
    return dict(domain=P.enc(start), derivative_partition=leaves,
                minimum_at_upper_corner=value.enc(), rejected_opposite_branches=branch_proofs)


def run():
    q,c,b,data=reduced();proof=prove(q,c,b)
    return dict(status='PASS_EXACT_CONTINUOUS_HIGH_CONTACT_CAVITY_FAMILY_INCOMPATIBILITY',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),INPUT,Path(X.__file__)]},
        reductions=data,proof=proof,
        theorem='Any norm2 right-S6 cap A with A.7=0, -A.12 and -A.47 in[7/5,2], and any norm2 right-S6 cap B with B.8=0,-B.12 in[7/5,2],-B.52 in[1,2], cannot both be acute to the retained chamber and each other. Acuteness selects their single branches, and A.B>13/100.',
        coxeter_scope='In particular A labels(2,m,n),m,n>=4, and B labels(2,p,q),p>=4,q>=3, including parallel limits, are incompatible. Other profiles and other cluster topologies are not excluded.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=run();path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['proof']['derivative_partition']),[len(x)for x in out['proof']['rejected_opposite_branches']],out['proof']['minimum_at_upper_corner'],flush=True)


if __name__=='__main__':main()
