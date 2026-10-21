#!/usr/bin/env python3
"""Fraction arithmetic in real quadratic towers; no search or numeric imports."""
from fractions import Fraction
from scripts.fresh_h4_blocks_independent_verify import (
    Z, O, P, add, sub, neg, mul, div, sign, scale, parse, enc)


class GoldenField:
    zero, one = Z, O
    add = staticmethod(add)
    sub = staticmethod(sub)
    neg = staticmethod(neg)
    mul = staticmethod(mul)
    div = staticmethod(div)
    sign = staticmethod(sign)
    encode = staticmethod(enc)
    parse = staticmethod(parse)

    @staticmethod
    def rational(value):
        return scale(O, Fraction(value))


K = GoldenField()


class Quadratic:
    """The distinguished real embedding with positive square root.

    Sign and polynomial evaluation remain valid when a radicand is a square.
    Division is used only for genuine field presentations or nonzero norm.
    """
    def __init__(self, base, radicand):
        assert base.sign(radicand) >= 0
        self.base, self.d = base, radicand
        self.zero = base.zero, base.zero
        self.one = base.one, base.zero

    def lift(self, value):
        return value, self.base.zero

    def rational(self, value):
        return self.lift(self.base.rational(value))

    def add(self, x, y):
        return self.base.add(x[0], y[0]), self.base.add(x[1], y[1])

    def neg(self, x):
        return self.base.neg(x[0]), self.base.neg(x[1])

    def sub(self, x, y):
        return self.add(x, self.neg(y))

    def mul(self, x, y):
        b = self.base
        return (b.add(b.mul(x[0], y[0]), b.mul(self.d, b.mul(x[1], y[1]))),
                b.add(b.mul(x[0], y[1]), b.mul(x[1], y[0])))

    def div(self, x, y):
        b = self.base
        norm = b.sub(b.mul(y[0], y[0]), b.mul(self.d, b.mul(y[1], y[1])))
        assert b.sign(norm)
        numerator = self.mul(x, (y[0], b.neg(y[1])))
        return b.div(numerator[0], norm), b.div(numerator[1], norm)

    def sign(self, x):
        b = self.base
        if not b.sign(self.d):
            return b.sign(x[0])
        sa, sb = b.sign(x[0]), b.sign(x[1])
        if not sb:
            return sa
        if not sa or sa == sb:
            return sb
        difference = b.sub(b.mul(x[0], x[0]), b.mul(self.d, b.mul(x[1], x[1])))
        return sa * b.sign(difference)

    def encode(self, x):
        return [self.base.encode(x[0]), self.base.encode(x[1])]

    def parse(self, x):
        assert len(x) == 2
        return self.base.parse(x[0]), self.base.parse(x[1])


def total(field, values):
    out = field.zero
    for value in values:
        out = field.add(out, value)
    return out


def dot(field, a, b):
    return total(field, (field.mul(x, y) for x, y in zip(a, b, strict=True)))


def matvec(field, matrix, vector):
    return tuple(dot(field, row, vector) for row in matrix)


def pairing(field, a, matrix, b):
    return dot(field, a, matvec(field, matrix, b))


def finite_label(field, value, polynomials):
    """Identify 4*cos(pi/m)^2 by its exact largest real root."""
    for label, coefficients, lower, upper in polynomials:
        if field.sign(field.sub(value, field.rational(lower))) < 0:
            continue
        if field.sign(field.sub(value, field.rational(upper))) > 0:
            continue
        result = field.zero
        for coefficient in coefficients:
            result = field.add(field.mul(result, value), field.rational(coefficient))
        if field.sign(result) == 0:
            return label
    return None


def old_pair_gate(field, pa, old, polynomials, allow_parallel=False):
    """New square-2 normal in field F versus one old quadratic-over-K normal."""
    # Lift K through E to F, then adjoin the old radical. The resulting
    # extension has degree at most16 over Q; dependence causes no sign error.
    lift = lambda x: field.lift(field.base.lift(x))
    extension = Quadratic(field, lift(old['d']))
    value = (
        dot(field, pa, tuple(lift(x) for x in old['a'])),
        dot(field, pa, tuple(lift(x) for x in old['b'])))
    if extension.sign(value) > 0:
        return {'compatible': False, 'kind': 'positive', 'label': None}
    ratio = extension.mul(extension.mul(value, value),
                          extension.lift(lift(div(scale(O, 2), old['q']))))
    comparison = extension.sign(extension.sub(ratio, extension.rational(4)))
    if comparison > 0:
        return {'compatible': True, 'kind': 'dotted', 'label': None}
    if comparison == 0:
        return {'compatible': allow_parallel, 'kind': 'parallel', 'label': None}
    label = finite_label(extension, ratio, polynomials)
    return {'compatible': label is not None,
            'kind': 'finite' if label is not None else 'noncoxeter',
            'label': label}


def new_pair_gate(field, pa, other_radicand, other_coefficients, polynomials,
                  allow_parallel=False):
    """Two square-2 new normals sharing the same quartic base E."""
    extension = Quadratic(field, field.lift(other_radicand))
    value = (dot(field, pa, tuple(field.lift(x[0]) for x in other_coefficients)),
             dot(field, pa, tuple(field.lift(x[1]) for x in other_coefficients)))
    if extension.sign(value) > 0:
        return {'compatible': False, 'kind': 'positive', 'label': None}
    ratio = extension.mul(value, value)  # 4p²/(2*2) = p².
    comparison = extension.sign(extension.sub(ratio, extension.rational(4)))
    if comparison > 0:
        return {'compatible': True, 'kind': 'dotted', 'label': None}
    if comparison == 0:
        return {'compatible': allow_parallel, 'kind': 'parallel', 'label': None}
    label = finite_label(extension, ratio, polynomials)
    return {'compatible': label is not None,
            'kind': 'finite' if label is not None else 'noncoxeter',
            'label': label}
