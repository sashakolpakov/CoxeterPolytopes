"""Direct sparse encodings of the same gated SDP for independent solvers.

The canonical form is min c'x, A x + s = b, s in K. All rows are
retained. PSD rows use lower-column svec with sqrt(2) off-diagonals;
Clarabel uses a separately checked upper-column permutation. Conversion of
the floating native input is not an exact certificate for the rational model.
"""
import json
from pathlib import Path


def read_native(path):
    import numpy as np
    import scipy.sparse as sp
    with Path(path).open('rb') as stream:
        comment = stream.readline()
        n = int(stream.readline())
        number_blocks = int(stream.readline())
        sizes = np.fromstring(stream.readline(), dtype=np.int64, sep=' ')
        c = np.fromstring(stream.readline(), dtype=float, sep=' ')
        entries = np.fromstring(stream.read(), dtype=float, sep=' ')
    if not comment.startswith(b'"Global Coxeter SDP'):
        raise ValueError('Unexpected native header')
    if len(sizes) != number_blocks or len(c) != n or np.any(sizes == 0) or entries.size % 5:
        raise ValueError('Incomplete native input')
    entries = entries.reshape((-1, 5))
    if not np.isfinite(entries).all() or not np.isfinite(c).all():
        raise ValueError('Nonfinite native data')
    coordinates = entries[:, :4].astype(np.int64)
    if not np.array_equal(coordinates, entries[:, :4]):
        raise ValueError('Noninteger native indices')
    matrix, block, i, j = coordinates.T
    block = block-1
    i = i-1
    j = j-1
    if np.any(matrix < 0) or np.any(matrix > n) or np.any(block < 0) or np.any(block >= len(sizes)):
        raise ValueError('Native index out of range')
    dimension = np.abs(sizes[block])
    if np.any(i < 0) or np.any(j < i) or np.any(j >= dimension):
        raise ValueError('Expected an upper triangular native entry')
    full = sizes[block] > 0
    if np.any(i[~full] != j[~full]):
        raise ValueError('Off-diagonal entry in a diagonal cone')
    scalar = int(-sizes[sizes < 0].sum())
    offsets = np.zeros(len(sizes), dtype=np.int64)
    scalar_position, psd_position = 0, scalar
    for k, size in enumerate(sizes):
        if size < 0:
            offsets[k] = scalar_position
            scalar_position -= int(size)
        else:
            offsets[k] = psd_position
            psd_position += int(size*(size+1)//2)
    rows = offsets[block]+i
    rows[full] = (offsets[block[full]] + i[full]*dimension[full]
                  - i[full]*(i[full]+1)//2 + j[full])
    values = entries[:, 4].copy()
    values[full & (i != j)] *= np.sqrt(2.)
    b = np.zeros(psd_position)
    constant = matrix == 0
    np.add.at(b, rows[constant], -values[constant])
    A = sp.coo_matrix((-values[~constant], (rows[~constant], matrix[~constant]-1)),
                      shape=(len(b), n)).tocsc()
    A.sum_duplicates()
    A.eliminate_zeros()  # Only literal floating zeros, never small coefficients.
    if not np.isfinite(A.data).all():
        raise ValueError('Nonfinite canonical data')
    return {'A': A, 'b': b, 'c': c}, {
        'native_block_sizes': sizes.tolist(), 'scalar_cones': scalar,
        'PSD_sizes': sizes[sizes > 0].tolist(), 'native_sparse_entries': len(entries),
        'canonical_rows': A.shape[0], 'canonical_columns': A.shape[1],
        'canonical_nonzeros': A.nnz,
        'canonical_form': 'min c^T x; A x+s=b; s in nonnegative x PSD cones',
        'source_F0_is_negative_affine_constant': True,
        'PSD_encoding': 'lower column-major svec, sqrt(2) off diagonal',
        'no_small_coefficient_threshold': True}


def save_encoding(folder, data, layout):
    import numpy as np
    import scipy.sparse as sp
    from . import conic, solver_gate
    folder = Path(folder)
    sp.save_npz(folder/'A.npz', data['A'])
    np.savez_compressed(folder/'vectors.npz', b=data['b'], c=data['c'])
    conic.write_json(folder/'layout.json', layout)
    return {name: solver_gate.file_digest(folder/name)
            for name in ['A.npz', 'vectors.npz', 'layout.json']}


def load_encoding(folder):
    import numpy as np
    import scipy.sparse as sp
    folder = Path(folder)
    vectors = np.load(folder/'vectors.npz', allow_pickle=False)
    return {'A': sp.load_npz(folder/'A.npz'), 'b': vectors['b'], 'c': vectors['c']}, json.loads((folder/'layout.json').read_text())


def original_feasibility(data, layout):
    """Set margin=0 and remove its redundant upper bound 1>=0 exactly."""
    import numpy as np
    scalar = layout['scalar_cones']
    row = scalar-1
    A, b, c = data['A'], data['b'], data['c']
    if row < 0 or c[-1] != -1 or np.count_nonzero(c) != 1:
        raise ValueError('Not the expected common-margin problem')
    bound = A.getrow(row)
    if b[row] != 1 or bound.nnz != 1 or bound.indices[0] != A.shape[1]-1 or bound.data[0] != 1:
        raise ValueError('The removed row is not precisely 1-margin>=0')
    keep = np.ones(A.shape[0], dtype=bool)
    keep[row] = False
    changed = dict(layout, scalar_cones=scalar-1, mode='original_feasibility',
                   exact_relation_to_auxiliary='Set last coordinate to zero, delete only 1>=0, set objective zero')
    return {'A': A[keep, :-1].tocsc(), 'b': b[keep], 'c': np.zeros(A.shape[1]-1)}, changed


def certificate_feasibility(data, layout):
    """A^T y=0, b^T y=-1, y in K is a normalized Farkas search."""
    import numpy as np
    import scipy.sparse as sp
    original, layout = original_feasibility(data, layout)
    A, b = original['A'], original['b']
    equalities = A.shape[1]+1
    matrix = sp.vstack([A.T, sp.csc_matrix(b.reshape((1, -1))),
                        -sp.eye(A.shape[0], format='csc')], format='csc')
    rhs = np.zeros(matrix.shape[0])
    rhs[equalities-1] = -1
    return {'A': matrix, 'b': rhs, 'c': np.zeros(A.shape[0])}, dict(
        layout, zero_cones=equalities, mode='normalized_certificate_feasibility',
        exact_relation_to_auxiliary='Original feasibility alternative: A^T y=0, b^T y=-1, y in K')


def upper_permutation(layout):
    import numpy as np
    start = layout.get('zero_cones', 0)+layout['scalar_cones']
    total = start+sum(n*(n+1)//2 for n in layout['PSD_sizes'])
    permutation = np.arange(total)
    for n in layout['PSD_sizes']:
        i, j = np.triu_indices(n)
        target = j*(j+1)//2+i
        permutation[start+target] = start+np.arange(len(i))
        start += len(i)
    return permutation


def native_arrays(path, layout):
    import numpy as np
    from . import conic, feasibility
    numbers, _ = conic.native_vector(path)
    x = np.array([float(v) for v in numbers])
    arrays = {}
    for field, name in [('xMat', 's'), ('yMat', 'y')]:
        text = feasibility.native_array(path, field)
        raw = np.fromstring(text.translate(str.maketrans({'{':' ', '}':' ', ',':' '})), sep=' ')
        packed = np.empty(layout['canonical_rows'])
        scalar_position, psd_position = 0, layout['scalar_cones']
        position = 0
        for n in layout['native_block_sizes']:
            if n < 0:
                packed[scalar_position:scalar_position-n] = raw[position:position-n]
                scalar_position -= n
                position -= n
            else:
                M = raw[position:position+n*n].reshape((n,n))
                if not np.allclose(M, M.T, rtol=0, atol=1e-12):
                    raise ValueError('Nonsymmetric native cone array')
                i, j = np.triu_indices(n)
                values = M[i,j].copy()
                values[i != j] *= np.sqrt(2.)
                packed[psd_position:psd_position+len(values)] = values
                psd_position += len(values)
                position += n*n
        if position != len(raw) or not np.isfinite(packed).all():
            raise ValueError('Incomplete or nonfinite native array')
        arrays[name] = packed
    return {'x': x, **arrays}


def cone_minimum(vector, layout):
    import numpy as np
    zero = layout.get('zero_cones', 0)
    scalar = layout['scalar_cones']
    minimum = float(np.min(vector[zero:zero+scalar])) if scalar else float('inf')
    position = zero+scalar
    for n in layout['PSD_sizes']:
        i,j = np.triu_indices(n)
        values = vector[position:position+len(i)].copy()
        values[i != j] /= np.sqrt(2.)
        M = np.zeros((n,n))
        M[i,j] = values
        M[j,i] = values
        minimum = min(minimum, float(np.linalg.eigvalsh(M)[0]))
        position += len(i)
    if position != len(vector):
        raise ValueError('Wrong complete cone dimension')
    return minimum


def diagnostics(data, layout, solution):
    import numpy as np
    A,b,c = data['A'], data['b'], data['c']
    x,y,s = [np.asarray(solution[k]) for k in ['x','y','s']]
    if not all(np.isfinite(v).all() for v in [x,y,s]):
        return {'finite_candidate':False, 'exact_certificate':False}
    affine = b-A@x
    adjoint = A.T@y+c
    return {'finite_candidate':True, 'primal_affine_residual_max':float(np.max(np.abs(affine-s))),
            'dual_adjoint_residual_max':float(np.max(np.abs(adjoint))),
            'primal_objective':float(c@x), 'dual_objective':float(-b@y),
            'gap':float(c@x+b@y), 'complementarity':float(s@y),
            'saved_primal_cone_minimum_eigenvalue':cone_minimum(s,layout),
            'actual_affine_cone_minimum_eigenvalue':cone_minimum(affine,layout),
            'saved_dual_cone_minimum_eigenvalue':cone_minimum(y,layout),
            'exact_certificate':False}
