"""Exact certificate coordinates and sparse elimination without svec radicals.

Matrix variables are actual upper-triangle entries. Off-diagonal trace terms
carry a factor of two. All equation arithmetic therefore stays over Q.
"""
from fractions import Fraction as F
from functools import reduce
import heapq
from math import gcd, lcm


def primitive(row):
    """Clear denominators locally and remove integer content, including sign."""
    row = {int(k): F(v) for k, v in row.items() if v}
    if not row:
        return {}
    denominator = lcm(*(v.denominator for v in row.values()))
    answer = {k: v.numerator * (denominator // v.denominator) for k, v in row.items()}
    content = reduce(gcd, map(abs, answer.values()))
    if answer[min(answer)] < 0:
        content = -content
    return {k: v // content for k, v in answer.items()}


def heights(rows):
    result = dict(nonzero_coefficients=0, maximum_numerator_bits=0,
                  maximum_denominator_bits=0, maximum_row_terms=0)
    for row in rows:
        result['maximum_row_terms'] = max(result['maximum_row_terms'], len(row))
        for value in row.values():
            value = F(value)
            result['nonzero_coefficients'] += 1
            result['maximum_numerator_bits'] = max(
                result['maximum_numerator_bits'], abs(value.numerator).bit_length())
            result['maximum_denominator_bits'] = max(
                result['maximum_denominator_bits'], value.denominator.bit_length())
    return result


def source_equations(problem, layout, embedding):
    """Return <F_j,Z>=0 and 1+<F_0,Z>=0 on the proved dual face."""
    import numpy as np
    count = layout['scalar_cones'] + sum(n*(n+1)//2 for n in layout['PSD_sizes'])
    inverse = np.full(count, -1, dtype=np.int64)
    inverse[embedding] = np.arange(len(embedding))
    position = {v: i for i, v in enumerate(problem['variables'])}
    equations = [{} for _ in range(len(position)+1)]
    equations[-1][-1] = F(1)
    scalar, full, base = 0, 0, layout['scalar_cones']
    for block in problem['blocks']:
        entries = block['entries']
        n = len(entries)
        if n == 1:
            coordinates = [(scalar, 0, 0)]
            scalar += 1
        else:
            if layout['PSD_sizes'][full] != n:
                raise ValueError('Different exact source and cone layout')
            coordinates = []
            for i in range(n):
                for j in range(i, n):
                    coordinates.append((base, i, j))
                    base += 1
            full += 1
        for old, i, j in coordinates:
            col = int(inverse[old])
            if col < 0:
                continue
            weight = 1 if i == j else 2
            for name, value in entries[i][j].items():
                row = position[name] if name else len(position)
                value = weight * F(value)
                if value:
                    equations[row][col] = value
    if scalar != layout['scalar_cones'] or full != len(layout['PSD_sizes']) or base != count:
        raise ValueError('Incomplete source cone traversal')
    return equations


def peel(rows, scalar_coordinates, progress=None, cascade=False):
    """Eliminate scalar columns of degree one, with exact triangular witnesses.

    The default uses only original degree-one columns: substitution cannot fill
    another equation, and its total size is bounded by the source coefficient
    count. Cascading is an explicit experiment; it can create substantial fill.
    """
    incidence = [set() for _ in range(scalar_coordinates)]
    remaining = {r for r, row in enumerate(rows) if row}
    for r in remaining:
        for col in rows[r]:
            if 0 <= col < scalar_coordinates:
                incidence[col].add(r)
    queue = []
    def offer(col):
        if len(incidence[col]) == 1:
            row = next(iter(incidence[col]))
            coefficient = abs(rows[row][col])
            heapq.heappush(queue, (len(rows[row]), coefficient.bit_length(), row, col))
    for col in range(scalar_coordinates):
        offer(col)
    steps = []
    while queue:
        _, _, row, col = heapq.heappop(queue)
        if row not in remaining or incidence[col] != {row}:
            continue
        steps.append(dict(source_row=row, pivot=col))
        remaining.remove(row)
        for other in rows[row]:
            if 0 <= other < scalar_coordinates:
                incidence[other].remove(row)
                if cascade:
                    offer(other)
        if progress and len(steps) % 10000 == 0:
            progress(stage='leaf_elimination', pivots=len(steps), remaining_rows=len(remaining))
    return steps, sorted(remaining)


def verify_peeling(rows, scalar_coordinates, steps, remaining):
    incidence = [set() for _ in range(scalar_coordinates)]
    active = {r for r, row in enumerate(rows) if row}
    for r in active:
        for col in rows[r]:
            if 0 <= col < scalar_coordinates:
                incidence[col].add(r)
    pivots = set()
    for step in steps:
        row, col = step['source_row'], step['pivot']
        if (not 0 <= col < scalar_coordinates or col in pivots or row not in active
                or incidence[col] != {row} or not rows[row].get(col)):
            raise ValueError('Invalid degree-one rational elimination witness')
        active.remove(row)
        pivots.add(col)
        for other in rows[row]:
            if 0 <= other < scalar_coordinates:
                incidence[other].remove(row)
    if active != set(remaining):
        raise ValueError('A surviving exact equation was lost')
    return pivots


def substitute(row, substitutions):
    result = {}
    for col, coefficient in row.items():
        replacement = substitutions.get(col)
        if replacement is None:
            result[col] = result.get(col, F(0)) + coefficient
            if not result[col]:
                del result[col]
            continue
        for other, value in replacement.items():
            result[other] = result.get(other, F(0)) + coefficient * value
            if not result[other]:
                del result[other]
    return result


def back_substitution(rows, steps, maximum_terms=15000000, maximum_bits=4096, progress=None):
    substitutions = {}
    terms = 0
    for step in reversed(steps):
        row, pivot = rows[step['source_row']], step['pivot']
        coefficient = F(row[pivot])
        expression = substitute({k: F(-v)/coefficient for k, v in row.items() if k != pivot},
                                substitutions)
        if any(max(abs(v.numerator).bit_length(), v.denominator.bit_length()) > maximum_bits
               for v in expression.values()):
            raise ValueError('Explicit rational coefficient growth limit reached')
        substitutions[pivot] = expression
        terms += len(expression)
        if terms > maximum_terms:
            raise ValueError('Explicit rational substitution fill limit reached')
        if progress and len(substitutions) % 10000 == 0:
            progress(stage='back_substitution', pivots=len(substitutions), terms=terms)
    return substitutions


def exact_remainders(rows, steps, remaining, substitutions):
    """Check all eliminated rows and preserve every residual equation over Q."""
    pivots = {s['source_row'] for s in steps}
    if len(pivots) != len(steps):
        raise ValueError('An exact source equation was used twice')
    residual = {}
    remaining = set(remaining)
    for index, row in enumerate(rows):
        value = substitute(row, substitutions)
        if index in pivots or not row:
            if value:
                raise ValueError('The exact substitution violates a removed equation')
        elif index in remaining:
            value = primitive(value)
            if value:
                residual[index] = value
        else:
            raise ValueError('An exact source equation is unaccounted for')
    return residual


def encode(substitutions, residual, variables, layout):
    """Convert an exactly reduced certificate to floating solver input once."""
    import numpy as np
    import scipy.sparse as sp
    free = np.asarray([i for i in range(variables) if i not in substitutions], dtype=np.int64)
    inverse = np.full(variables, -1, dtype=np.int64)
    inverse[free] = np.arange(len(free))
    z = len(residual)
    factors = np.ones(variables)
    base = layout['scalar_cones']
    for n in layout['PSD_sizes']:
        i, j = np.triu_indices(n)
        factors[base:base+len(i)][i != j] = np.sqrt(2.)
        base += len(i)
    if base != variables:
        raise ValueError('The reduced certificate cone has a different dimension')
    rr, cc, vv = [], [], []
    rhs = np.zeros(z+variables)
    for r, row in enumerate(residual.values()):
        rhs[r] = -float(row.get(-1, 0))
        for col, value in row.items():
            if col >= 0:
                if inverse[col] < 0:
                    raise ValueError('A residual equation contains an eliminated coordinate')
                rr.append(r); cc.append(int(inverse[col])); vv.append(float(value))
    for col in range(variables):
        row = substitutions.get(col, {col: F(1)})
        r = z + col
        rhs[r] = factors[col] * float(row.get(-1, 0))
        for other, value in row.items():
            if other >= 0:
                if inverse[other] < 0:
                    raise ValueError('A reconstruction still contains an eliminated coordinate')
                rr.append(r); cc.append(int(inverse[other])); vv.append(-factors[col]*float(value))
    A = sp.coo_matrix((vv, (rr, cc)), shape=(z+variables, len(free))).tocsc()
    return dict(A=A, b=rhs, c=np.zeros(len(free))), dict(
        zero_cones=z, scalar_cones=layout['scalar_cones'], PSD_sizes=layout['PSD_sizes'],
        mode='rationally_reduced_normalized_certificate_feasibility',
        exact_coordinates='Actual matrix entries; trace factor two; svec sqrt(2) only at export',
        no_claim_of_original_primal_feasibility_from_reduced_dual_infeasibility=True), free
