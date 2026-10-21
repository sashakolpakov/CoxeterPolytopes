#!/usr/bin/env python3
"""Complete intrinsic S6 ledger after the exact release54/replacement12."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_release54_replace12_20261017 as C
import fresh_nonarithmetic_release34_s6_boundary_20261017 as B

ROOT=C.ROOT
TAG='fresh_nonarithmetic_release54_replace12_boundary_20261017'


def run():
    primary=json.loads((ROOT/'results'/(C.TAG+'.json')).read_text())
    assert primary['status']=='PASS_EXACT_64_WALL_RELEASE54_REPLACEMENT12'
    oldpath=ROOT/'results'/(B.TAG+'.json')
    before=json.loads(oldpath.read_text())['after65']
    z=C.geometry()
    assert primary['new_normal']==C.P.enc(z['normals'][12])
    assert primary['common_center']==C.P.enc(z['common_center'])
    after=B.ledger(z['normals'],z['pair'],z['common_center'],True)
    assert before['open_ends']==17 and after['open_ends']==16
    assert after['ray_counts']==dict(finite_vertex=24,exterior_ray=8)
    assert after['edge_counts']==dict(compact_geodesic_edge=30,one_open_end=12,two_open_ends=2,outside_hyperbolic_space=4)
    assert len(after['actual_sides'])==18 and len(after['triples'])==816
    assert not after['singular_triples'] and not after['nonsimple_rays'] and not after['ideal_rays']
    a={tuple(r['sides']):r for r in after['edges']}
    b={tuple(r['sides']):r for r in before['edges']}
    assert a.keys()==b.keys()
    changes=[dict(sides=list(k),before=b[k]['kind'],after=a[k]['kind'])for k in sorted(a)if a[k]['kind']!=b[k]['kind']]
    assert changes==[
        dict(sides=[8,12],before='one_open_end',after='compact_geodesic_edge'),
        dict(sides=[8,52],before='one_open_end',after='compact_geodesic_edge'),
        dict(sides=[12,52],before='outside_hyperbolic_space',after='one_open_end')]
    compact_before={r['wall'] for r in before['facets']if r['compact']}
    compact_after={r['wall'] for r in after['facets']if r['compact']}
    assert compact_before<=compact_after
    paths=[Path(__file__),Path(C.__file__),Path(B.__file__),oldpath,ROOT/'results'/(C.TAG+'.json')]
    return dict(status='PASS_COMPLETE_EXACT_S6_BOUNDARY_17_TO_16_OPEN_ENDS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        before65=before,after64=after,changed_edges=changes,
        compact_side_faces_lost=sorted(compact_before-compact_after),compact_side_faces_gained=sorted(compact_after-compact_before),
        completeness='Frozen exact residual4 chart/triple enumeration is reused with the new exact normal. All816 triple systems are nonsingular;32simple rays and48projective edges are exhausted. All causal projected walls are strictly future timelike and redundant on the future H3. Exact positive-cone Gram tests distinguish compact, one-open, two-open and outside edges.',
        scope='Complete intrinsic S6 comparison only. Old54 is not an S6 side; label12 is a replacement. Two old open edges close, one previously outside edge becomes one-open. Complete ambient finite census197 and a strict escape are supplied by the parent geometry; no full ambient exterior cone or finite volume is claimed.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)


if __name__=='__main__':main()
