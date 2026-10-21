#!/usr/bin/env python3
"""Independently reconstruct two full-chamber symmetries of the new H9 host."""
import argparse
import json
import fresh_h9_beta131_h300_symmetry_20261016_review as I

K, ROOT = I.K, I.ROOT
TAG = 'fresh_h9_witt_twist_symmetry_20261016_review'
PREFIX = 'results/fresh_h9_witt_twist_pilot_20261016.json'
MAPS = 'results/fresh_h9_witt_twist_pilot_20261016_maps.json'
PREFIX_REVIEW = 'results/fresh_h9_witt_twist_pilot_20261016_review.json'
SETUP_REVIEW = 'results/fresh_h9_witt_twist_20261016_review.json'
FRAMES = [([0,1,2,4,6,7,8,9,11], [4,3,2,0,6,7,8,11,9]),
          ([0,1,2,3,4,5,6,8,10], [10,14,2,3,4,5,6,8,0])]


def run():
    prefix, maps, census, setup = map(I.read, (PREFIX, MAPS, PREFIX_REVIEW, SETUP_REVIEW))
    assert census['status'] == 'PASS_INDEPENDENT_COMPLETE_WITT_TWIST_PREFIX'
    assert census['controls'] == 13 and census['walls'] == 16 and census['height_bound'] == '20'
    assert setup['status'] == 'PASS_INDEPENDENT_DISTINCT_WITT_HOST_NORMS_AND_ZERO_SECTOR'
    phi, beta = K(0,1), K(4,-1)
    diagonal = [-phi]+[K(1)]*7+[beta]*2
    form = [[diagonal[i] if i == j else K() for j in range(10)] for i in range(10)]
    assert form == [I.dec(row) for row in prefix['form']] == [I.dec(row) for row in setup['form']]
    roots = [I.dec(row) for row in prefix['roots']]
    h = [K(1)]+[K()]*9
    assert I.dec(prefix['basepoint']) == h

    def pair(a, b):
        return sum((d*x*y for d,x,y in zip(diagonal,a,b)),K())

    assert all(d.conjugate() > 0 for d in diagonal)
    assert pair(h,h) < 0 and all(pair(h,n) <= 0 for n in roots)
    assert all(pair(a,b) <= 0 for i,a in enumerate(roots) for b in roots[:i])
    center = I.dec(prefix['strict_center'])
    assert pair(center,center) < 0 and pair(center,h) < 0
    assert all(pair(center,n) < 0 for n in roots)
    records = []
    for which,(left,right) in zip((14,3),FRAMES):
        saved = maps['integral_maps'][which]
        assert saved['source_ids'] == left and saved['target_ids'] == right
        corners = []
        spherical = None
        for ids in (left,right):
            ns = [roots[i] for i in ids]
            gram = [[pair(a,b) for b in ns] for a in ns]
            pivots = I.pivots(gram)
            inverse,_ = I.inverse(gram)
            assert all(x >= 0 for row in inverse for x in row)
            b = I.mv(inverse,[pair(h,n) for n in ns])
            assert all(x <= 0 for x in b)
            v = [h[k]-sum((b[j]*ns[j][k] for j in range(9)),K()) for k in range(10)]
            norm = pair(v,v)
            assert norm < 0 and pair(v,h) == norm
            vals = [pair(v,n) for n in roots]
            assert all(x <= 0 for x in vals)
            assert [i for i,x in enumerate(vals) if not x] == sorted(ids)
            if spherical is None:
                spherical = gram
            else:
                assert spherical == gram
            corners.append(dict(ids=ids,apex=v,norm=norm,coefficients=b,pairings=vals,pivots=pivots))
        v,w = [c['apex'] for c in corners]
        a = I.trans([roots[i] for i in left]+[v])
        b = I.trans([roots[i] for i in right]+[w])
        ai,da = I.inverse(a)
        _,db = I.inverse(b)
        scale = da/db
        if scale < 0:
            scale = -scale
        assert scale > 0 and scale*scale*pair(w,w) == pair(v,v)
        for row in b:
            row[-1] *= scale
        t = I.mm(b,ai)
        ti,det = I.inverse(t)
        assert all(x.is_integral for row in t+ti for x in row)
        assert I.mm(I.trans(t),I.mm(form,t)) == form
        assert I.mm(t,ti) == I.mm(ti,t) == I.identity(10)
        for i,j in zip(left,right):
            assert I.mv(t,roots[i]) == roots[j]
        assert I.mv(t,v) == [scale*x for x in w] and pair(v,w) < 0
        assert t == [I.dec(row) for row in saved['matrix']]
        assert ti == [I.dec(row) for row in saved['inverse']]
        records.append(dict(map_index=which,corners=corners,matrix=t,inverse=ti,determinant=det,
                            positive_apex_scale=scale,trace=sum((t[i][i] for i in range(10)),K())))
    product = I.mm(records[0]['matrix'],records[1]['matrix'])
    inverse,det = I.inverse(product)
    assert all(x.is_integral for row in product+inverse for x in row)
    assert det == -1 and I.mm(I.trans(product),I.mm(form,product)) == form
    trace = sum((product[i][i] for i in range(10)),K())
    assert trace == 50+76*phi and trace > 10
    assert product == [I.dec(row) for row in maps['hit']['matrix']]
    assert maps['hit']['kind'] == 'product' and maps['hit']['maps'] == [14,3]
    paths = [f'scripts/{TAG}.py','scripts/fresh_h9_beta131_h300_symmetry_20261016_review.py',
             'src/coxeter_polytopes/exact/golden.py',PREFIX,MAPS,PREFIX_REVIEW,SETUP_REVIEW]
    return I.enc(dict(status='PASS_INDEPENDENT_WITT_TWIST_FULL_CHAMBER_NONREFLECTIVITY',
                      sources={p:I.sha(p) for p in paths},form=form,basepoint=h,generators=records,
                      product=product,inverse_product=inverse,product_determinant=det,trace=trace,
                      checked_corners=4,checked_wall_signs=64,
                      conclusion='This locally distinct maximal-integral H9 lattice is nonreflective; its full chamber stabilizer contains the displayed infinite-order product.',
                      scope='Four literal finite sectors, complete independently reviewed height20 prefix, all-height primitive norm alphabet. No dependency on full cone/corner enumeration or exhaustive map search. No conclusion for other lattices or arbitrary geometric walls.'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    args = ap.parse_args()
    result = run()
    path = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['trace'])
