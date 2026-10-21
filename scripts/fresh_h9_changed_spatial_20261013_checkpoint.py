#!/usr/bin/env python3
"""Final D9 geometry, live third edge and scoped corner-transport checks."""
from fresh_h9_changed_spatial_20261013_steps import *

SOURCE=ROOT/'results'/f'{PREFIX}_diagonal.json'
CHECKPOINT=ROOT/'results'/f'{PREFIX}_checkpoint.json'


def run(verify=False):
    d=json.loads(SOURCE.read_text());fd=json.loads(OUT.read_text())
    walls=[tuple(map(K.from_json,r))for r in d['all_walls']]
    inside=tuple(map(K.from_json,fd['strict_interior']))
    escape=(O,*(K(F(x,100))for x in(80,65,50,35,20,10,5,3,1)))
    assert pair(inside,G,inside)<0 and all(pair(inside,G,r)<0 for r in walls)
    assert pair(escape,G,escape)==K(F(4977,2000))-T>0
    assert pair(escape,G,inside)<0 and all(pair(escape,G,r)<0 for r in walls)
    for i,r in enumerate(walls):
        assert root_reason(r,pair(r,G,r))is None
        assert all(compatible(r,s)for s in walls[:i])
    # The original third edge points along (1,1,1,0,...,0). Its simple
    # spacelike continuation obeys every current wall, so no existing cap
    # closes that ray before the light cone.
    edge_escape=(O,K(F(27,50)),K(F(27,50)),K(F(27,50)),*(Z for _ in range(6)))
    assert pair(edge_escape,G,edge_escape)>0
    assert all(pair(edge_escape,G,r)<=0 for r in walls)
    assert [i for i in range(9)if pair(edge_escape,G,walls[i])==0]==[0,1,3,4,5,6,7,8]
    _,_,_,_,_,eps=read_interface();A=eye()
    for j in range(1,10):A=mm(boost(j,eps),A)
    P=[list(row)for row in eye()]
    for i in range(1,10):
        P[i]=[Z]*10;P[i][10-i]=-O
    P=tuple(map(tuple,P));H=mm(inverse(A),P)
    assert mm(H,H)==eye() and mm(mm(tuple(zip(*H)),G),H)==G
    assert mv(H,(O,*(Z for _ in range(9))))[0]>0
    initial_image=[mv(H,r)for r in INITIAL]
    expected=[mv(inverse(A),r)for r in corner(9)]
    assert set(initial_image)==set(expected)
    failures=[{'basis_index':i,'image':enc(mv(H,b)),'spatial_sum_over_two':(sum(mv(H,b)[1:],Z)/2).to_json()}
              for i,b in enumerate(DBASIS)if not lattice(mv(H,b))]
    assert failures
    swaps=[]
    for k in(2,7):
        B=eye()
        for j in range(1,k+1):B=mm(boost(j,eps),B)
        R=[list(row)for row in eye()]
        ids=(1,2)if k==2 else(8,9)
        i,j=ids
        R[i][i]=(-1 if k==2 else 1)/SQRT2
        R[i][j]=1/SQRT2;R[j][i]=1/SQRT2
        R[j][j]=(1 if k==2 else -1)/SQRT2
        R=tuple(map(tuple,R));J=mm(mm(inverse(B),R),B)
        local=corner(k);short,long=(local[0],local[1])if k==2 else(local[-1],local[-2])
        assert mv(R,short)==tuple(SQRT2*x for x in long)
        assert mv(R,long)==tuple(x/SQRT2 for x in short)
        assert mm(J,J)==eye() and mm(mm(tuple(zip(*J)),G),J)==G
        bad=[i for i,b in enumerate(DBASIS)if not lattice(mv(J,b))]
        assert bad
        swaps.append({'corner_stage':k,'matrix':menc(J),'nonintegral_basis_images':bad})
    result={'status':'PASS_EXACT_D9_19_FACET_CHECKPOINT_WITH_STRICT_ESCAPE',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,OUT)},
        'walls':19,'all_pairs_checked':171,'strict_interior':enc(inside),
        'strict_escape':enc(escape),'strict_escape_norm':pair(escape,G,escape).to_json(),
        'strict_escape_pairings':enc(tuple(pair(escape,G,r)for r in walls)),
        'known_complete_first_edges':10,'known_full_finite_corners':11,
        'next_initial_edge':{'omitted_initial_wall':2,'support':[0,1,3,4,5,6,7,8],
            'type':'A2+B6','initial_corner':enc((O,*(Z for _ in range(9)))),
            'spacelike_continuation':enc(edge_escape),'norm':pair(edge_escape,G,edge_escape).to_json(),
            'all_wall_pairings':enc(tuple(pair(edge_escape,G,r)for r in walls)),
            'candidate_relative_field':'K(sqrt(3*(T-1)))',
            'scope':'This current-cone edge reaches the light cone; its first mirror in the full arithmetic group remains uncomputed.'},
        'opposite_B9_transport':{'matrix':menc(H),'square_is_identity':True,
            'form_and_oriented_sector_map':True,'D9_lattice_failures':failures},
        'B2_diagram_swap_tests':swaps,
        'symmetry_scope':'The tested opposite-B9 frame transport is an involution and does not preserve the D9 lattice. The two natural B2 swaps also fail lattice preservation. These tests do not establish global nonreflectivity or exclude other chamber symmetries.',
        'volume_scope':'All19 normals are primitive reflective, all171 pairs pass exact Coxeter-angle gates, and a common strict timelike point exists. Acute projection makes every inequality a facet. The strict future spacelike point proves infinite volume of this prefix; full-host reflectivity remains open.'}
    if verify:assert result==json.loads(CHECKPOINT.read_text())
    else:CHECKPOINT.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();run(args.verify)
