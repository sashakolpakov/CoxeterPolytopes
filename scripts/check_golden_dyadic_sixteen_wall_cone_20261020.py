#!/usr/bin/env python3
"""Fixed all-sixteen spacelike witness, not a cone or height enumeration.

Run only through the parent-lock bounded wrapper. --write refuses replacement;
--verify recomputes every fixed identity in a separate process. No H9 example
or full-lattice obstruction follows from this finite-prefix certificate.
"""
import argparse
import hashlib
import json
import os
import signal
from fractions import Fraction as F
from pathlib import Path
from check_arithmetic_dn_opposed_fork_models_20261019 import field, inverse

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/golden_dyadic_sixteen_wall_cone_20261020.json'
PINS = {
    'results/golden_dyadic_sixteen_wall_caps_20261020.md':
        'c288b5de619471772c8381d98afa3996ac8851b4d2d146ad520d6474ee8ba92b',
    'results/golden_dyadic_sixteen_wall_caps_20261020.json':
        '80f18b60003375876f5bae9c122c9466ff1eb1c69fed582d87c7facf50fd5745',
    'results/golden_dyadic_sixteen_wall_cone_preparation_20261020.md':
        '6cdeff97990a342b3f5e7afc5026c89cf17bc0f747b853a563fd80e618d9233d',
    'results/golden_dyadic_d4_fourth_vertex_and_sealed_h3_face_20261020.md':
        'd6b148954b09e259d0e0598cc4f754904ba1977b373a3f491196ef87f93f9501',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
        '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def check():
    sources = {p: digest(ROOT/p) for p in PINS}
    assert sources == PINS, 'pinned scientific/control source changed'
    own_path = Path(__file__).relative_to(ROOT).as_posix()
    sources[own_path] = digest(ROOT/own_path)
    data = json.loads((ROOT/'results/golden_dyadic_sixteen_wall_caps_20261020.json').read_text())
    Q = field(5)
    zero, phi = Q(0), Q(F(1, 2), F(1, 2))
    p = [Q(1)]
    for _ in range(8):
        p.append(p[-1]*phi)
    assert phi*phi == phi+1 and phi.conjugate().sign() < 0

    def dec(v):
        return Q(F(v[0]))+Q(F(v[1]))*phi

    def dv(v):
        return [dec(x) for x in v]

    def enc(v):
        return [x.encoded() for x in v]

    def add(v, w):
        return [a+b for a, b in zip(v, w)]

    def scale(a, v):
        return [a*x for x in v]

    def coordinate(i):
        return [Q(int(i == j)) for j in range(10)]

    g3 = [[Q(4), Q(-2), zero],
          [Q(-2), Q(4), -2*phi], [zero, -2*phi, Q(4)]]

    def dot(v, w):
        return (-phi*v[0]*w[0]
                +sum((v[i]*w[i] for i in range(1, 7)), zero)
                +sum((v[7+i]*g3[i][j]*w[7+j]
                      for i in range(3) for j in range(3)), zero))

    def original_dot(v, w):
        return -phi*v[0]*w[0]+sum((v[i]*w[i] for i in range(1, 10)), zero)

    # Columns of the exact original R5 embedding, not just a K-similarity.
    embedding = [
        [p[3]]+[zero]*5+[-p[2], -phi, zero, -p[3]],
        [-phi]+[zero]*5+[phi, zero, zero, phi],
    ]
    embedding += [coordinate(i) for i in range(1, 6)]
    embedding += [scale(Q(2), coordinate(8)),
                  [zero]*6+[Q(-1), Q(-1), Q(-1), Q(1)],
                  [Q(2)]+[zero]*8+[-2*phi]]
    assert embedding == [dv(v) for v in data['basis_in_original_coordinates']]
    assert dec(data['metric']['time_square']) == -phi
    assert data['metric']['spatial_identity_rank'] == 6
    assert [dv(row) for row in data['metric']['last_block']] == g3
    assert data['basis_order'] == [
        'V0', 'e1=-c', 'e2=I1', 'e3=I2', 'e4=I3', 'e5=I4', 'e6=I5', 'a3', 's', 'r']
    assert all(original_dot(embedding[i], embedding[j]) == dot(coordinate(i), coordinate(j))
               for i in range(10) for j in range(10))

    def original(v):
        return [sum((v[i]*embedding[i][j] for i in range(10)), zero)
                for j in range(10)]

    def in_lattice(v):
        return all(x.integral() for x in v) and (sum(v[6:], zero)/2).integral()

    # This is a complete O-basis of <-phi>+I5+D4, in original coordinates.
    lattice_basis = [coordinate(i) for i in range(6)]
    for i, j, sign in ((6, 7, -1), (7, 8, -1), (8, 9, -1), (8, 9, 1)):
        lattice_basis.append(add(coordinate(i), scale(Q(sign), coordinate(j))))
    assert len(lattice_basis) == 10 and all(in_lattice(v) for v in lattice_basis)

    beta = [add(scale(Q(-1), coordinate(i)), coordinate(i+1))
            for i in range(1, 6)]+[scale(Q(-1), coordinate(6))]
    wp = [-(phi+2)/2, -p[2], -p[3]/2]
    normals = {**{'beta'+str(i+1): v for i, v in enumerate(beta)},
               'a3': coordinate(7), 's': coordinate(8), 'r': coordinate(9),
               'u': [phi-1, phi]+[zero]*8,
               'w': [p[2]]+[phi]*5+[zero]*4,
               'h': [2*p[3], 2*p[2]]+[p[3]]*5+[zero]*3,
               'a2': [Q(1)]+[zero]*6+wp,
               'A_a2': [12*phi+7]+[6*phi+4]*6+wp,
               'W': [2*p[3]]+[p[3]]*6+[-p[2], -2*p[2], -p[3]],
               'Z': [2*p[3]]+[p[3]]*6+scale(phi, wp)}
    names = list(normals)
    old_names, new_names = names[:9], names[9:]
    assert len(names) == 16 and data['initial_corner'] == old_names
    assert data['added_order'] == new_names
    assert [row['id'] for row in data['normals']] == names
    original_normals, norms = {}, {}
    for row in data['normals']:
        name, v = row['id'], normals[row['id']]
        assert v == dv(row['coordinates']) and len(v) == 10
        q = dot(v, v)
        assert q == dec(row['norm']) and q.sign() > 0
        ov = original(v)
        assert in_lattice(ov) and original_dot(ov, ov) == q
        for b in lattice_basis:
            reflected = add(b, scale(-2*original_dot(b, ov)/q, ov))
            assert in_lattice(reflected), 'reflection does not preserve R5: '+name
        original_normals[name], norms[name] = ov, q
    assert [[dot(normals[n], normals[b]) for b in old_names] for n in new_names] == [
        dv(row) for row in data['new_to_corner_products']]
    assert [[dot(normals[n], normals[m]) for m in new_names[i+1:]]
            for i, n in enumerate(new_names)] == [
        dv(row) for row in data['new_to_new_upper_triangle']]
    pairs = []
    for i, n in enumerate(names):
        for m in names[i+1:]:
            value = dot(normals[n], normals[m])
            assert value == original_dot(original_normals[n], original_normals[m])
            assert value.sign() <= 0
            product = norms[n]*norms[m]
            if value == 0:
                kind = 'right'
            elif (value*value-product).sign() > 0:
                kind = 'ultraparallel'
            else:
                x = 4*value*value/product
                kind = next((str(label) for label, target in ((3, Q(1)), (4, Q(2)), (5, p[2]))
                             if x == target), None)
                assert kind is not None, 'unexpected native angle or tangency: '+n+'/'+m
            pairs.append({'walls': [n, m], 'product': value.encoded(), 'type': kind})
    assert len(pairs) == 120

    # PD old Gram. Positive pivots justify the definitions-only helper inverse.
    gram = [[dot(normals[n], normals[m]) for m in old_names] for n in old_names]
    rows, pivots = [list(row) for row in gram], []
    for i in range(9):
        pivot = rows[i][i]
        assert pivot.sign() > 0
        pivots.append(pivot)
        for j in range(i+1, 9):
            ratio = rows[j][i]/pivot
            for k in range(i+1, 9):
                rows[j][k] = rows[j][k]-ratio*rows[i][k]
    c = inverse(gram, Q)
    assert all(x.sign() >= 0 for row in c for x in row)
    dual = [[sum((c[i][k]*normals[old_names[k]][j] for k in range(9)), zero)
             for j in range(10)] for i in range(9)]
    assert all(dot(dual[i], normals[old_names[j]]) == int(i == j)
               for i in range(9) for j in range(9))

    def corner(v):
        return [v[0]]+[-dot(v, normals[n]) for n in old_names]

    def from_corner(y):
        return [y[0]*int(j == 0)-sum((y[i+1]*dual[i][j] for i in range(9)), zero)
                for j in range(10)]

    def value(row, v):
        return sum((a*b for a, b in zip(row, v)), zero)

    native_rows = [[dot(coordinate(0), normals[n])]+[-dot(v, normals[n]) for v in dual]
                   for n in names]
    all_rows = [[Q(-1)]+[zero]*9]+native_rows
    eps = Q(F(1, 100))
    strict_corner = [Q(1)]+[eps]*9
    strict = from_corner(strict_corner)
    assert strict == dv(data['strict_point']) and corner(strict) == strict_corner
    strict_values = [dot(strict, normals[n]) for n in names]
    assert all(x.sign() < 0 for x in strict_values)
    strict_norm = -phi+eps*eps*(383+24*phi)/4
    assert dot(strict, strict) == strict_norm and strict_norm.sign() < 0
    assert original(strict)[0].sign() > 0

    escape = [Q(1)]+[Q(1)/phi]*3+[Q(1)/p[2]]*3+[
        (1-phi)/2, Q(F(-1, 2)), (4-3*phi)/2]
    esc = data['escape']
    assert escape == dv(esc['coordinates'])
    assert escape[7:] == add(scale(-Q(1)/p[3], inverse(g3, Q)[0]),
                               scale(-Q(1)/p[5], inverse(g3, Q)[2]))
    assert dot(escape[:7]+[zero]*3, escape[:7]+[zero]*3) == 21-13*phi
    assert dot([zero]*7+escape[7:], [zero]*7+escape[7:]) == 26-16*phi
    square = dot(escape, escape)
    assert square == 47-29*phi == dec(esc['norm']) and square.sign() > 0
    twice_original = scale(Q(2), original(escape))
    assert twice_original == dv(esc['original_integral_coordinates'])
    assert in_lattice(twice_original) and twice_original[0].sign() > 0
    assert original_dot(twice_original, twice_original) == 4*square
    y = corner(escape)
    assert y == [Q(1), zero, zero, 2*phi-3, zero, zero, 2-phi, 2*phi-3, zero, 5*phi-8]
    assert all(x.sign() >= 0 for x in y) and from_corner(y) == escape
    values = [dot(escape, normals[n]) for n in names]
    assert values == dv(esc['wall_products_in_normal_order'])
    assert all(x.sign() <= 0 for x in values)
    assert [value(row, y) for row in native_rows] == values
    assert [value(row, strict_corner) for row in native_rows] == strict_values
    active = [i for i, v in enumerate(values) if v == 0]
    assert [names[i] for i in active] == ['beta1', 'beta2', 'beta4', 'beta5', 's', 'u', 'w', 'W', 'Z']
    # One 9x10 exact rank calculation, not a matrix or extreme-ray farm.
    reduced, pivot_columns, top = [list(native_rows[i]) for i in active], [], 0
    for column in range(10):
        pivot = next((i for i in range(top, 9) if reduced[i][column] != 0), None)
        if pivot is None:
            continue
        reduced[top], reduced[pivot] = reduced[pivot], reduced[top]
        divisor = reduced[top][column]
        reduced[top] = [x/divisor for x in reduced[top]]
        for i in range(top+1, 9):
            ratio = reduced[i][column]
            reduced[i] = [a-ratio*b for a, b in zip(reduced[i], reduced[top])]
        pivot_columns.append(column)
        top += 1
        if top == 9:
            break
    assert top == 9
    mask = sum(1 << (i+1) for i in active)  # bit0 is the strict future-time row.
    assert [i for i, row in enumerate(all_rows) if value(row, y) == 0] == [i+1 for i in active]
    assert all(digest(ROOT/path) == pin for path, pin in sources.items()), 'source changed'
    return {
        'status': 'EXACT_SIXTEEN_WALL_SPACELIKE_EXTREME_RAY_PASS',
        'owns_lock': False, 'sources': sources,
        'field_encoding': '[a,b] means a+b*sqrt5; rational strings',
        'input_encoding': '[a,b] means a+b*phi, independently decoded',
        'scope': 'Infinite-volume finite16 prefix; full integral R5 chamber remains open.',
        'metric_signature_positive_negative_zero': [9, 1, 0],
        'conjugate_original_metric_positive_diagonal': enc([-phi.conjugate()]+[Q(1)]*9),
        'normal_names': names, 'norms': enc([norms[n] for n in names]),
        'original_normals': {n: enc(v) for n, v in original_normals.items()},
        'reflections_on_full_integral_basis_checked': 160,
        'native_pair_count': len(pairs), 'native_pairs': pairs,
        'old_corner_gram': [enc(row) for row in gram],
        'old_corner_gram_inverse': [enc(row) for row in c], 'pd_pivots': enc(pivots),
        'corner_constraint_names': ['future_time']+names,
        'corner_constraint_rows': [enc(row) for row in all_rows],
        'strict_point': enc(strict), 'strict_point_corner': enc(strict_corner),
        'strict_point_original': enc(original(strict)), 'strict_point_norm': strict_norm.encoded(),
        'strict_point_wall_products': enc(strict_values),
        'escape': enc(escape), 'escape_corner': enc(y),
        'escape_original_integral_2E': enc(twice_original), 'escape_norm': square.encoded(),
        'escape_wall_products': enc(values), 'escape_full_active_mask': mask,
        'escape_active_native_walls': [names[i] for i in active],
        'escape_active_rank': top, 'escape_active_rank_pivot_columns': pivot_columns,
        'written_implication': 'strict future interior plus strict spacelike perturbation gives an open ideal-boundary patch; no full-arrangement conclusion',
        'searches_or_DD_stages': 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.write and CERT.exists():
        parser.error('refusing to overwrite existing certificate')
    assert int(os.environ.get('COXETER_CHECK_OWNER', '-1')) == os.getppid(), 'bounded parent wrapper required'
    signal.alarm(45)
    data = check()
    payload = json.dumps(data, sort_keys=True, indent=2)+'\n'
    assert len(payload.encode()) <= 131072, 'certificate exceeds128KiB cap'
    if args.write:
        with CERT.open('x') as stream:
            stream.write(payload)
    else:
        assert CERT.stat().st_size <= 131072
        assert data == json.loads(CERT.read_text()), 'saved certificate differs'
    print(json.dumps({'status': 'PASS', 'mode': 'write' if args.write else 'verify',
                      'native_pairs': 120, 'active_rank': 9, 'result': 'SPACELIKE_ESCAPE',
                      'certificate_sha256': hashlib.sha256(payload.encode()).hexdigest(),
                      'certificate_bytes': len(payload.encode())}))


if __name__ == '__main__':
    main()
