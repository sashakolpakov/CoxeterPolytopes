#!/usr/bin/env python3
"""Exact small local models; run through the bounded research wrapper.

No polytope completion or exhaustive classification is asserted. Arithmetic
uses only Fraction in Q(sqrt2) or Q(sqrt5). Generation refuses replacement;
--verify recomputes the same certificate in a separate process.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'results/arithmetic_dn_opposed_fork_models_20261019.json'


def field(discriminant):
    class Q:
        __slots__ = ('a', 'b')

        def __init__(self, a=0, b=0):
            self.a, self.b = F(a), F(b)

        def __add__(self, other):
            other = other if isinstance(other, Q) else Q(other)
            return Q(self.a + other.a, self.b + other.b)

        __radd__ = __add__

        def __neg__(self):
            return Q(-self.a, -self.b)

        def __sub__(self, other):
            return self + (-other if isinstance(other, Q) else Q(-other))

        def __rsub__(self, other):
            return -self + other

        def __mul__(self, other):
            other = other if isinstance(other, Q) else Q(other)
            return Q(self.a * other.a + discriminant * self.b * other.b,
                     self.a * other.b + self.b * other.a)

        __rmul__ = __mul__

        def conjugate(self):
            return Q(self.a, -self.b)

        def __truediv__(self, other):
            other = other if isinstance(other, Q) else Q(other)
            norm = other.a * other.a - discriminant * other.b * other.b
            assert norm
            product = self * other.conjugate()
            return Q(product.a / norm, product.b / norm)

        def __eq__(self, other):
            other = other if isinstance(other, Q) else Q(other)
            return self.a == other.a and self.b == other.b

        def sign(self):
            sign = lambda a: (a > 0) - (a < 0)
            sa, sb = sign(self.a), sign(self.b)
            if not sb:
                return sa
            if not sa or sa == sb:
                return sb
            difference = self.a * self.a - discriminant * self.b * self.b
            return sa * sign(difference)

        def encoded(self):
            return [str(self.a), str(self.b)]

        def integral(self):
            if discriminant == 2:
                return self.a.denominator == self.b.denominator == 1
            aa, bb = 2 * self.a, 2 * self.b
            return (aa.denominator == bb.denominator == 1
                    and (int(aa) - int(bb)) % 2 == 0)

    return Q


def inertia(matrix, conjugate=False):
    a = [[x.conjugate() if conjugate else x for x in row] for row in matrix]
    counts = [0, 0, 0]  # positive, negative, zero
    for i in range(len(a)):
        pivot = a[i][i]
        sign = pivot.sign()
        counts[{1: 0, -1: 1, 0: 2}[sign]] += 1
        if not sign:
            assert all(a[i][j] == 0 for j in range(i, len(a)))
            continue
        for j in range(i + 1, len(a)):
            for k in range(j, len(a)):
                a[j][k] = a[j][k] - a[j][i] * a[i][k] / pivot
                a[k][j] = a[j][k]
    return counts


def inverse(matrix, Q):
    n = len(matrix)
    rows = [list(row) + [Q(int(i == j)) for j in range(n)]
            for i, row in enumerate(matrix)]
    for i in range(n):
        pivot = rows[i][i]
        assert pivot != 0
        rows[i] = [entry / pivot for entry in rows[i]]
        for j in range(n):
            if j == i:
                continue
            ratio = rows[j][i]
            rows[j] = [x - ratio * y for x, y in zip(rows[j], rows[i])]
    assert all(rows[i][j] == int(i == j) for i in range(n) for j in range(n))
    return [row[n:] for row in rows]


def bilinear(left, matrix, right):
    return sum(left[i] * sum(matrix[i][j] * right[j] for j in range(len(right)))
               for i in range(len(left)))


def submatrix(g, indices):
    return [[g[i][j] for j in indices] for i in indices]


def finite_types(g, indices):
    adjacency = {i: {j for j in indices if i != j and g[i][j] != 0}
                 for i in indices}
    assert all(g[i][j] == -F(1, 2) for i in indices for j in adjacency[i])
    remaining, types = set(indices), []
    while remaining:
        stack, component = [min(remaining)], set()
        while stack:
            i = stack.pop()
            if i not in component:
                component.add(i)
                stack.extend(adjacency[i] - component)
        remaining -= component
        degrees = [len(adjacency[i]) for i in component]
        n = len(component)
        assert sum(degrees) == 2 * (n - 1)
        if max(degrees) <= 2:
            types.append('A' + str(n))
        else:
            assert degrees.count(3) == 1 and max(degrees) == 3
            fork = next(i for i in component if len(adjacency[i]) == 3)
            arms = []
            for first in adjacency[fork]:
                previous, current, length = fork, first, 1
                while len(adjacency[current]) == 2:
                    nxt = next(j for j in adjacency[current] if j != previous)
                    previous, current, length = current, nxt, length + 1
                arms.append(length)
            assert sorted(arms) == [1, 1, n - 3]
            types.append('D' + str(n))
    return sorted(types)


def model(n, discriminant):
    Q = field(discriminant)
    x = Q(F(1, 4), F(1, 4)) if discriminant == 5 else Q(F(1, 2), F(1, 2))
    beta = x + F(1, 2) if discriminant == 5 else x
    assert x * (beta - 1) == F(1, 4)
    a, b, e, c, d = 0, 1, 2, n, n + 1
    chain = list(range(3, n))
    g = [[Q(int(i == j)) for j in range(n + 2)] for i in range(n + 2)]
    for i, j in [(a, b), (a, e), (a, chain[0]), *zip(chain, chain[1:]),
                 (c, e), (c, chain[-1]), (c, d)]:
        g[i][j] = g[j][i] = Q(-F(1, 2))
    g[a][c] = g[c][a] = -x
    g[b][d] = g[d][b] = -beta
    assert all((2 * entry).integral() for row in g for entry in row)
    assert all(g[i][j].sign() <= 0 for i in range(n + 2) for j in range(i))
    physical, conjugate = inertia(g), inertia(g, True)
    assert physical == [n, 1, 1] and conjugate == [n + 1, 0, 1]
    kernel = [Q(1) for _ in range(n + 2)]
    kernel[b] = kernel[d] = -2 * x
    assert all(sum(entry * coefficient for entry, coefficient in zip(row, kernel)) == 0
               for row in g)
    assert all(coefficient != 0 and coefficient.conjugate() != 0 for coefficient in kernel)
    white, black = [e, *chain], [a, b, c, d]
    wi = inverse(submatrix(g, white), Q)
    projected = [[g[i][j] - bilinear([g[i][k] for k in white], wi,
                                    [g[k][j] for k in white]) for j in black] for i in black]
    s = Q(F(n, 4 * (n - 2)))
    assert projected == [[s, Q(-F(1, 2)), -x - s, Q()],
                         [Q(-F(1, 2)), Q(1), Q(), -beta],
                         [-x - s, Q(), s, Q(-F(1, 2))],
                         [Q(), -beta, Q(-F(1, 2)), Q(1)]]
    corner_types = []
    for i in range(4):
        corner = white + [black[i], black[(i + 1) % 4]]
        assert inertia(submatrix(g, corner)) == [n, 0, 0]
        corner_types.append(finite_types(g, corner))
    assert corner_types == [['D' + str(n)], ['A1', 'A' + str(n - 1)],
                            ['D' + str(n)], ['A1', 'A' + str(n - 1)]]
    assert (2 * (1 - beta)).sign() < 0
    assert all((g[i][b] + g[i][d]).sign() < 0 for i in black)
    old_inverse = inverse(submatrix(g, list(range(n))), Q)
    contact = g[c][:n]
    swapped = list(contact)
    swapped[b], swapped[e] = swapped[e], swapped[b]
    q = bilinear(contact, old_inverse, contact)
    qq = bilinear(swapped, old_inverse, swapped)
    cross = bilinear(contact, old_inverse, swapped)
    assert q == qq and (q - 1).sign() > 0
    assert q - cross == F(1, 4)
    assert cross - (q - 1) == F(3, 4)
    mixed_gap = None
    if discriminant == 5:
        triangle = [Q() for _ in range(n)]
        triangle[e] = -x
        triangle_cost = bilinear(triangle, old_inverse, triangle)
        mixed_cross = bilinear(contact, old_inverse, triangle)
        assert (triangle_cost - 1).sign() > 0 and mixed_cross.sign() > 0
        mixed_gap = mixed_cross * mixed_cross - (q - 1) * (triangle_cost - 1)
        assert mixed_gap == (22 * (n - 1) * x + 7 * (n - 2) + 3) / 16
        assert mixed_gap.sign() > 0
    return dict(n=n, field_discriminant=discriminant, x=x.encoded(), beta=beta.encoded(),
                physical_inertia=physical, conjugate_inertia=conjugate,
                finite_corner_types=corner_types, white_projection_checked=True,
                integral_norm_two_gram=True, full_support_kernel_checked=True,
                strict_timelike_quad_center_checked=True,
                twin_wall_product=(cross - (q - 1)).encoded(),
                golden_mixed_triangle_quad_gap=mixed_gap.encoded() if mixed_gap is not None else None,
                baseline_corner_weights=[1, 0, 1, 0])


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(65536), b''):
            result.update(chunk)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    reports = ['results/arithmetic_dn_golden_opposed_fork_quad_survivor_20261019.md',
               'results/arithmetic_dn_silver_opposed_fork_quad_survivor_20261019.md',
               'results/arithmetic_dn_golden_triangle_quad_capacity_20261019.md']
    result = dict(status='PASS_EXACT_LOCAL_MODELS_NOT_GLOBAL_POLYTOPES',
                  checker_sha256=digest(Path(__file__)),
                  report_sha256={name: digest(ROOT / name) for name in reports},
                  models=[model(n, d) for d in (5, 2) for n in (19, 20)],
                  scope='Four local face stars and their twin-wall capacity only; no complete census, global compact completion or dimension exclusion.')
    if args.verify:
        assert json.loads(CERT.read_text()) == result, 'Saved certificate differs'
    else:
        with CERT.open('x') as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write('\n')
    print(json.dumps(dict(status=result['status'], models=len(result['models']),
                          replay=args.verify, certificate=str(CERT.relative_to(ROOT)))))


if __name__ == '__main__':
    main()
