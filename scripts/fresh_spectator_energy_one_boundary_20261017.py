#!/usr/bin/env python3
"""Exact identities for the all-rank strict spectator boundary D=1."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import fresh_spectator_energy_below_one_20261017 as B

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_spectator_energy_one_boundary_20261017'


def energy(g, x, inverse_x):
    for i, row in enumerate(g):
        assert sum((a*b for a,b in zip(row,inverse_x)),F()) == x[i]
    return sum((a*b for a,b in zip(x,inverse_x)),F())


def run():
    prior = json.loads((ROOT/f'results/{B.TAG}.json').read_text())
    assert prior['status'] == 'PASS_ALL_RANK_SPECTATOR_ENERGY_BELOW_ONE'
    for name, h in prior['sources'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest() == h, name
    one = []
    for kind, start in (('B',2),('D',3)):
        for n in range(start,65):
            if kind == 'B':
                edges = [(i,i+1,2 if i==n-2 else 1) for i in range(n-1)]
                g = B.matrix(n,edges,True)
                v = [F(1)]*(n-1)+[F(1,2)]
            else:
                edges = [(i,i+1,1) for i in range(n-3)] + [(n-3,n-2,1),(n-3,n-1,1)]
                g = B.matrix(n,edges)
                v = [F(1)]*(n-2)+[F(1,2)]*2
            x = [F(int(i==0)) for i in range(n)]
            assert energy(g,x,v)==1 and all(t>0 for t in v)
            pivots = B.tree_pivots(g)
            one.append(dict(kind='A3_middle' if kind=='D' and n==3 else kind,
                rank=n,contact=0,energy='1',strict_extension_schur='1',
                positive_inverse_contact_column=list(map(str,v)),positive_pivots=pivots))

    twins=[]
    for n in range(3,65):
        edges=[(i,i+1,1) for i in range(n-3)]+[(n-3,n-2,1),(n-3,n-1,1)]
        g=B.matrix(n,edges)
        x=[F(int(i==n-2)) for i in range(n)]
        y=[F(int(i==n-1)) for i in range(n)]
        ix=[F(i+1,2) for i in range(n-2)]+[F(n,4),F(n-2,4)]
        iy=ix[:-2]+list(reversed(ix[-2:]))
        xx=energy(g,x,ix);yy=energy(g,y,iy)
        z=[a-b for a,b in zip(x,y)];iz=[a-b for a,b in zip(ix,iy)]
        assert xx==yy==F(n,4) and energy(g,z,iz)==1
        assert iz==[t/2 for t in z]
        allowed=xx<2 and yy<2
        assert allowed==(n<=7)
        twins.append(dict(kind='A3_endpoints' if n==3 else 'D_short_fork_leaves',rank=n,
            separate_projection_norm=str(xx),separate_extension_schur=str(2-xx),
            energy='1',strictly_allowed=allowed,inverse_difference=list(map(str,iz))))

    isolated=[]
    for name,diagonal,x,y in (
        ('one_sqrt2_contact_on_A1_scaled',[4],[2],[0]),
        ('two_unit_contacts_same_side_on_A1_plus_A1',[2,2],[1,1],[0,0]),
        ('two_unit_contacts_opposite_sides_on_A1_plus_A1',[2,2],[1,0],[0,1]),
    ):
        g=[[F(diagonal[i] if i==j else 0) for j in range(len(diagonal))] for i in range(len(diagonal))]
        x=list(map(F,x));y=list(map(F,y))
        ix=[a/b for a,b in zip(x,diagonal)];iy=[a/b for a,b in zip(y,diagonal)]
        xx=energy(g,x,ix);yy=energy(g,y,iy)
        assert xx<2 and yy<2
        assert energy(g,[a-b for a,b in zip(x,y)],[a-b for a,b in zip(ix,iy)])==1
        isolated.append(dict(kind=name,gram=list(map(str,diagonal)),x=list(map(str,x)),y=list(map(str,y)),
            separate_projection_norms=list(map(str,(xx,yy))),energy='1'))

    # The elementary Cauchy equality used to classify nonadjacent contacts:
    # G^{-1}(x-y)=(e_i-e_j)/2. Every other coordinate of G(e_i-e_j)/2
    # is in(-1,1), hence cannot be a nonzero native contact.
    assert F(1)>2*F(2,3)**2  # adjacent nonunit energy is strictly greater than1
    assert len(one)==125 and len(twins)==62
    paths=[Path(__file__)]+[ROOT/f'{folder}/{B.TAG}.{ext}' for folder,ext in
        (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ALL_RANK_SPECTATOR_ENERGY_ONE_BOUNDARY',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        hypotheses=prior['hypotheses']|dict(energy='D=1'),
        one_nonempty_support=one,opposite_twin_contacts=twins,isolated_cases=isolated,
        counts=dict(one_contact_BD_identities=125,twin_identities=62,isolated_identities=3),
        theorem='Up to exchanging x,y and diagram automorphism: one side only has one sqrt2 contact on whole A1, two unit contacts on distinct whole A1 components, or one unit contact at whole A3 middle, whole B_r remote leaf(r>=2), or whole D_r long leaf(r>=4). Two nonempty sides have unit contacts on distinct whole A1 components, on whole A3 endpoints, or the two short leaves of whole D_r with4<=r<=7. No other contacts occur.',
        arbitrary_rank_proof='Written Cauchy equality, native gap, strict positive inverse-column sealing, fork/path and exact all-rank inverse-column identities. Finite identity samples corroborate the formulas.',
        scope='Universal local boundary classification. No completed corner transfer, dimension improvement, finite-volume example or independent full arithmetic reconstruction.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],flush=True)
