#!/usr/bin/env python3
"""Separate exact check of the literal E8/F4-squared host isometry."""
import argparse
import json
from hashlib import sha256
from pathlib import Path
import fresh_silver_f4_independent_verify as I
import fresh_silver_e8_full_chamber_20261016_review as E

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = 'fresh_silver_e8_f4_identification_20261016'
TAG = PRIMARY + '_review'


def run():
    data = I.read(f'results/{PRIMARY}.json')
    source = I.host(I.UNIT,2,True);target = E.G
    assert source == I.matrix(data['source_gram']) and target == I.matrix(data['target_gram'])
    matrix = I.matrix(data['matrix']);inverse = I.inverse(matrix)
    assert inverse == I.matrix(data['inverse'])
    assert all(I.integral(x) for row in matrix+inverse for x in row)
    assert I.mm(matrix,inverse) == I.mm(inverse,matrix) == I.eye(10)
    assert I.mm(I.tr(matrix),I.mm(target,matrix)) == source
    assert I.determinant(matrix) == I.O
    base = (I.O,)+(I.Z,)*9;image = I.mv(matrix,base)
    assert I.sign(I.pair(image,target,base))<0
    assert I.pair(image,target,image) == I.pair(base,source,base) == I.neg(I.UNIT)
    assert image == tuple(map(I.parse,data['target_basepoint']))
    paths = [f'scripts/{TAG}.py',f'results/{PRIMARY}.json',
             'scripts/fresh_silver_f4_independent_verify.py',f'scripts/{E.TAG}.py']
    for p,h in data['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return dict(status='PASS_SEPARATE_FRACTION_E8_F4_INTEGRAL_ISOMETRY_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                rank=10,integral_inverse=True,gram_congruence=True,determinant=I.enc(I.O),
                future_orientation=True,basepoint_image=list(map(I.enc,image)),
                conclusion='The full Lorentz lattices, including time coordinate, are integrally isometric. Their integral reflection groups are conjugate; the E8 host duplicates the previously excluded F4-squared host.',
                scope='Independent matrix/Gram arithmetic using the older standalone field engine and separately reconstructed E8 Gram. Does not repeat how the map was discovered or the finite-corner census; neither is needed to verify this explicit lattice isomorphism.')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
