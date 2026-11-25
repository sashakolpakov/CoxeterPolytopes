"""Exact affine conic compiler and certificate checking for the global model.

Models use L_0 + sum x_i L_i >= 0. Equalities are eliminated over Q.
Native SDPA output is only a proposal; check_point checks original data.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re

ZERO = F(0)
ONE = F(1)


def form(constant=0, **terms):
    result = {key: F(value) for key, value in terms.items() if value}
    if constant:
        result[''] = F(constant)
    return result


def plus(*forms):
    result = {}
    for expression in forms:
        for key, value in expression.items():
            result[key] = result.get(key, ZERO) + value
            if not result[key]:
                del result[key]
    return result


def times(scalar, expression):
    scalar = F(scalar)
    return {key: scalar * value for key, value in expression.items() if scalar * value}


def value(expression, point):
    return expression.get('', ZERO) + sum(
        (coefficient * point[name] for name, coefficient in expression.items() if name), ZERO)


def digest_bytes(data):
    return hashlib.sha256(data).hexdigest()


def plain(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, dict):
        return {str(key): plain(item) for key, item in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [plain(item) for item in obj]
    return obj


def canonical(obj):
    return json.dumps(plain(obj), sort_keys=True, separators=(',', ':')).encode()


def write_json(path, obj):
    with Path(path).open('x') as stream:
        json.dump(plain(obj), stream, indent=2, sort_keys=True)
        stream.write('\n')


def model(variables, blocks, equations=(), metadata=None):
    names = list(variables)
    if len(set(names)) != len(names) or any(not name for name in names):
        raise ValueError('Repeated or empty variable')
    allowed = set(names) | {''}
    for block in blocks:
        entries = block['entries']
        if not entries or any(len(row) != len(entries) for row in entries):
            raise ValueError('Nonsquare PSD block')
        for i, row in enumerate(entries):
            for j, entry in enumerate(row):
                if entry != entries[j][i] or not set(entry) <= allowed:
                    raise ValueError('Nonsymmetric or unregistered affine entry')
    if any(not set(row) <= allowed for row in equations):
        raise ValueError('Unregistered equality coefficient')
    return dict(schema='coxeter-global-affine-sdp-v1', variables=names,
        blocks=list(blocks), equations=list(equations), metadata=metadata or {})


def read_model(path):
    raw = json.loads(Path(path).read_text())
    for block in raw['blocks']:
        block['entries'] = [[{key: F(v) for key, v in entry.items()} for entry in row]
                            for row in block['entries']]
    raw['equations'] = [{key: F(v) for key, v in row.items()} for row in raw['equations']]
    return model(raw['variables'], raw['blocks'], raw['equations'], raw['metadata'])


def eliminate(problem, extra_equations=()):
    """Sparse rational RREF with an exact original-coordinate substitution."""
    variables = problem['variables']
    rank = {name: i for i, name in enumerate(variables)}
    pivots = {}
    rows = problem['equations'] + list(extra_equations)
    for original in rows:
        row = dict(original)
        for pivot in sorted(pivots, key=rank.get):
            if pivot in row:
                row = plus(row, times(-row[pivot], pivots[pivot]))
        keys = [name for name in row if name]
        if not keys:
            if row.get('', ZERO):
                raise ValueError('Inconsistent exact linear equations')
            continue
        pivot = min(keys, key=rank.get)
        row = times(ONE / row[pivot], row)
        for old in list(pivots):
            if pivot in pivots[old]:
                pivots[old] = plus(pivots[old], times(-pivots[old][pivot], row))
        pivots[pivot] = row
    free = [name for name in variables if name not in pivots]
    substitutions = {}
    for name in variables:
        if name in pivots:
            substitutions[name] = times(-ONE, {k: v for k, v in pivots[name].items() if k != name})
        else:
            substitutions[name] = form(**{name: ONE})

    def substitute(expression):
        answer = form(expression.get('', ZERO))
        for name, coefficient in expression.items():
            if name:
                answer = plus(answer, times(coefficient, substitutions[name]))
        return answer

    blocks = [dict(name=block['name'], entries=[[substitute(entry) for entry in row]
                    for row in block['entries']]) for block in problem['blocks']]
    used = {name for block in blocks for row in block['entries'] for entry in row for name in entry if name}
    null = [name for name in free if name not in used]
    free = [name for name in free if name in used]
    if null:
        # Redundant original equalities give free dual gauges that do not
        # change any cone entry. Fix them to zero, and preserve the full
        # original-coordinate reconstruction with that same gauge choice.
        removed = set(null)
        substitutions = {name: {k: v for k, v in expression.items() if k not in removed}
            for name, expression in substitutions.items()}
    compiled = model(free, blocks, metadata=dict(original_variables=variables,
        exact_equality_rank=len(pivots), original_model_sha256=digest_bytes(canonical(problem))))
    compiled['metadata']['free_null_directions_fixed_to_zero'] = null
    return compiled, substitutions


def decimal_text(number, precision=90):
    with localcontext() as context:
        context.prec = precision
        return str(Decimal(number.numerator) / Decimal(number.denominator))


def write_sdpa(path, problem):
    if problem['equations']:
        raise ValueError('Eliminate equalities before exporting SDPA')
    variables = problem['variables']
    if not variables:
        raise ValueError('SDPA needs a nonempty free variable vector')
    # SDPA's negative block size represents a diagonal nonnegative cone.
    # Pack scalar cones together, preserving every original constraint.
    full = [b for b in problem['blocks'] if len(b['entries']) > 1]
    scalar = [b for b in problem['blocks'] if len(b['entries']) == 1]
    sizes = [len(b['entries']) for b in full]+([-len(scalar)] if scalar else [])
    indices = {name: i+1 for i, name in enumerate(variables)}
    indices[''] = 0
    nonzeros = 0
    with Path(path).open('x') as stream:
        stream.write('"Global Coxeter SDP; numerical decimals, exact model is separate\n')
        stream.write(f'{len(variables)}\n{len(sizes)}\n')
        stream.write(' '.join(map(str, sizes)) + '\n')
        stream.write(' '.join('0' for _ in variables) + '\n')
        def entry(block, i, j, expression):
            nonlocal nonzeros
            for variable, coefficient in expression.items():
                if coefficient:
                    matrix = indices[variable]
                    coefficient = coefficient if matrix else -coefficient
                    stream.write(f'{matrix} {block} {i} {j} {decimal_text(coefficient)}\n')
                    nonzeros += 1
        for block_index, block in enumerate(full, 1):
            for i, row in enumerate(block['entries']):
                for j in range(i, len(row)):
                    entry(block_index, i+1, j+1, row[j])
        for i, block in enumerate(scalar, 1):
            entry(len(full)+1, i, i, block['entries'][0][0])
    return dict(free_variables=len(variables), full_PSD_blocks=len(full),
        packed_scalar_constraints=len(scalar), sparse_nonzeros=nonzeros)


def farkas(problem):
    """Exact alternative for ALL primal PSD blocks and affine equalities.

    For F_b(x)>=0, E(x)=0, seek Z_b>=0 and free y such that
    every x coefficient in sum <F_b,Z_b>+y E vanishes and its
    constant is -1. Off-diagonal trace entries carry multiplicity two.
    """
    variables, blocks = [], []
    equations = {name: {} for name in problem['variables']}
    normalization = form(1)
    for i, row in enumerate(problem['equations']):
        name = f"potential_{i}"
        variables.append(name)
        normalization = plus(normalization, form(**{name: row.get('', ZERO)}))
        for old, coefficient in row.items():
            if old:
                equations[old] = plus(equations[old], form(**{name: coefficient}))
    for b, block in enumerate(problem['blocks']):
        size = len(block['entries'])
        entries = [[{} for _ in range(size)] for _ in range(size)]
        for i in range(size):
            for j in range(i, size):
                name = f"dual_b{b}_{i}_{j}"
                variables.append(name)
                entries[i][j] = entries[j][i] = form(**{name: ONE})
                weight = ONE if i == j else F(2)
                for old, coefficient in block['entries'][i][j].items():
                    if old:
                        equations[old] = plus(equations[old], form(**{name: weight*coefficient}))
                    else:
                        normalization = plus(normalization, form(**{name: weight*coefficient}))
        blocks.append(dict(name='dual_'+block['name'], entries=entries))
    return model(variables, blocks, list(equations.values())+[normalization], metadata=dict(
        problem_role='infeasibility_certificate_for_global_necessary_data',
        necessary_model_sha256=digest_bytes(canonical(problem)),
        primal_equalities_included=True,
        corner_weights_are_free_equality_potentials=True,
        interpretation='Feasibility proves the necessary-data primal infeasible; failure to find a certificate is inconclusive.'))


def psd_pivots(matrix):
    """Exact Schur/LDL test, including singular matrices without jitter."""
    working = [list(row) for row in matrix]
    size = len(working)
    diagonal = []
    for i in range(size):
        if any(working[i][j] != working[j][i] for j in range(size)):
            raise ValueError('Nonsymmetric rational matrix')
        pivot = working[i][i]
        if pivot < 0:
            raise ValueError(f'Negative exact PSD pivot at {i}')
        diagonal.append(pivot)
        if pivot == 0:
            if any(working[i][j] for j in range(i+1, size)):
                raise ValueError(f'Nonzero row at zero PSD pivot {i}')
            continue
        for j in range(i+1, size):
            for k in range(j, size):
                entry = working[j][k] - working[j][i] * working[i][k] / pivot
                working[j][k] = working[k][j] = entry
    return diagonal


def check_point(problem, point):
    if set(point) != set(problem['variables']):
        raise ValueError('Wrong complete candidate variable set')
    for index, row in enumerate(problem['equations']):
        if value(row, point):
            raise ValueError(f'Nonzero exact equality {index}')
    checks = []
    for block in problem['blocks']:
        evaluated = [[value(entry, point) for entry in row] for row in block['entries']]
        try:
            pivots = psd_pivots(evaluated)
        except ValueError as error:
            raise ValueError(block['name'] + ': ' + str(error)) from error
        checks.append(dict(name=block['name'], order=len(evaluated), exact_ldl_pivots=pivots))
    return dict(equalities_checked=len(problem['equations']), blocks=checks,
        positive_semidefinite=True, candidate_sha256=digest_bytes(canonical(point)))


NUMBER = re.compile(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?')


def native_vector(path):
    text = Path(path).read_text()
    if len(text) > 128 * 1024**2:
        raise ValueError('Oversized SDPA output')
    phase = re.search(r'phase\.value\s*=\s*(\S+)', text)
    match = re.search(r'xVec\s*=\s*\{([^}]*)\}', text, re.S)
    if not match:
        raise ValueError('SDPA output has no xVec')
    values = [F(token) for token in NUMBER.findall(match.group(1))]
    return values, None if phase is None else phase.group(1)


def reconstruct(problem, compiled, substitutions, raw, extra_equations=()):
    if len(raw) != len(compiled['variables']):
        raise ValueError('Wrong native vector length')
    attempts = []
    for digits in (6, 10, 16, 24, 36, 48, 60):
        grid = 10**digits
        free = {name: F(round(number * grid), grid)
                for name, number in zip(compiled['variables'], raw, strict=True)}
        point = {name: value(expression, free) for name, expression in substitutions.items()}
        try:
            checked = check_point(problem, point)
            if any(value(row, point) for row in extra_equations):
                raise ValueError('Polishing equation failed')
            return point, dict(method='common-grid-rounding-and-exact-equality-projection',
                decimal_grid_digits=digits, attempts=attempts, check=checked)
        except ValueError as error:
            attempts.append(dict(decimal_grid_digits=digits, reason=str(error)))
    raise ValueError('No exact PSD reconstruction: ' + json.dumps(attempts))


def reconstruct_on_face(problem, compiled, substitutions, raw):
    """Propose rational kernels from the numerical point; verify original SDP.

    Kernel restrictions only reduce the candidate search. No PSD jitter,
    matrix perturbation or relaxed equality is accepted in the final check.
    """
    import numpy as np
    original_free = dict(zip(compiled['variables'], raw, strict=True))
    approximate = {name: value(expression, original_free)
                   for name, expression in substitutions.items()}
    failures = []
    for tolerance in (1e-10, 1e-16, 1e-24):
        extra = []
        kernels = []
        for block in problem['blocks']:
            entries = block['entries']
            size = len(entries)
            matrix = np.array([[float(value(entry, approximate)) for entry in row]
                               for row in entries], dtype=float)
            bound = tolerance * max(1.0, float(np.max(np.abs(matrix))))
            if size == 1:
                if abs(matrix[0, 0]) <= bound:
                    extra.append(entries[0][0])
                continue
            work = matrix.copy()
            pivot_columns = []
            row_index = 0
            for column in range(size):
                if row_index == size:
                    break
                pivot_row = row_index + int(np.argmax(np.abs(work[row_index:, column])))
                if abs(work[pivot_row, column]) <= bound:
                    continue
                work[[row_index, pivot_row]] = work[[pivot_row, row_index]]
                work[row_index] /= work[row_index, column]
                for other in range(size):
                    if other != row_index:
                        work[other] -= work[other, column]*work[row_index]
                pivot_columns.append(column)
                row_index += 1
            for column in range(size):
                if column in pivot_columns:
                    continue
                vector = [ZERO] * size
                vector[column] = ONE
                for r, pivot in enumerate(pivot_columns):
                    number = float(-work[r, column])
                    vector[pivot] = ZERO if abs(number) <= tolerance else F(number).limit_denominator(10**8)
                floating = np.array([float(x) for x in vector])
                if float(np.max(np.abs(matrix @ floating))) > 100*bound*max(1.0, float(np.max(np.abs(floating)))):
                    continue
                kernels.append(dict(block=block['name'], vector=vector))
                for row in entries:
                    extra.append(plus(*(times(coefficient, entry)
                                        for coefficient, entry in zip(vector, row, strict=True))))
        extra = [row for row in extra if row]
        try:
            refined, refined_substitutions = eliminate(problem, extra)
            refined_raw = [approximate[name] for name in refined['variables']]
            point, check = reconstruct(problem, refined, refined_substitutions, refined_raw, extra)
            check.update(method='numerical-face-proposal-followed-by-exact-equality-projection',
                proposed_kernel_tolerance=str(tolerance), proposed_kernels=kernels,
                additional_reconstruction_equalities=extra, earlier_face_attempts=failures)
            return point, check
        except ValueError as error:
            failures.append(dict(tolerance=str(tolerance), reason=str(error)))
    raise ValueError('Face reconstruction failed: ' + json.dumps(failures))


def source_snapshot(paths):
    answer = {}
    for path in paths:
        location = Path(path)
        if not location.is_file() or location.stat().st_size > 4 * 1024**2:
            raise ValueError('Missing or oversized source input')
        answer[str(path)] = dict(sha256=digest_bytes(location.read_bytes()),
                                 bytes=location.stat().st_size)
    return answer
