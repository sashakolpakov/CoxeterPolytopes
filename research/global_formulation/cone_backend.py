"""Direct sparse conic solver adapters for an unchanged affine SDP.

Clarabel uses upper-triangular column-major svec; SCS uses lower-triangular
column-major svec. Both multiply off-diagonal entries by sqrt(2).
Every order-two PSD cone may equivalently use the Lorentz vector
(a+c, 2b, a-c). Dual matrices are reconstructed by the adjoint map.
"""
from contextlib import contextmanager
from fractions import Fraction as F
import os
import sys
import time

from . import conic


@contextmanager
def native_output(path):
    sys.stdout.flush(); sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    try:
        with path.open("x") as stream:
            os.dup2(stream.fileno(), 1); os.dup2(stream.fileno(), 2)
            yield
            sys.stdout.flush(); sys.stderr.flush()
    finally:
        for fd, old in zip((1, 2), saved, strict=True):
            os.dup2(old, fd); os.close(old)


def pack(problem, solver, soc2=True):
    import numpy as np
    from scipy.sparse import coo_matrix
    columns = {n: i for i, n in enumerate(problem["variables"])}
    scalars, soc, psd = [], [], []
    for b, block in enumerate(problem["blocks"]):
        r = len(block["entries"])
        (scalars if r == 1 else soc if r == 2 and soc2 else psd).append((b, block))
    ii, jj, vv, rhs, layout = [], [], [], [], []
    def row(q, factor=1.0):
        i = len(rhs)
        rhs.append(factor*float(q.get("", F(0))))
        for v, c in q.items():
            if v: ii.append(i); jj.append(columns[v]); vv.append(-factor*float(c))
    for b, block in scalars:
        start = len(rhs); row(block["entries"][0][0]); layout.append(("scalar", b, start, 1))
    for b, block in soc:
        M = block["entries"]; start = len(rhs)
        row(conic.plus(M[0][0], M[1][1])); row(conic.times(2, M[0][1]))
        row(conic.plus(M[0][0], conic.times(-1, M[1][1])))
        layout.append(("soc", b, start, 2))
    for b, block in psd:
        M = block["entries"]; r = len(M); start = len(rhs)
        for j in range(r):
            for i in (range(j+1) if solver == "clarabel" else range(j, r)):
                row(M[i][j], 1.0 if i == j else np.sqrt(2.0))
        layout.append(("psd", b, start, r))
    A = coo_matrix((vv, (ii, jj)), shape=(len(rhs), len(columns))).tocsc()
    return A, np.asarray(rhs), layout, dict(nonnegative=len(scalars), soc=[3]*len(soc),
        psd=[len(b["entries"]) for _, b in psd])


def unpack_dual(problem, layout, vector, solver):
    import numpy as np
    if vector is None or not np.all(np.isfinite(vector)): return None
    matrices = [None]*len(problem["blocks"])
    for kind, b, start, r in layout:
        if kind == "scalar": M = [[float(vector[start])]]
        elif kind == "soc":
            u, v, w = map(float, vector[start:start+3])
            M = [[u+w, v], [v, u-w]]
        else:
            M = [[0.0]*r for _ in range(r)]; k = start
            for j in range(r):
                for i in (range(j+1) if solver == "clarabel" else range(j, r)):
                    M[i][j] = M[j][i] = float(vector[k])/(1.0 if i == j else np.sqrt(2.0)); k += 1
        matrices[b] = [[F(str(v)) for v in row] for row in M]
    return matrices


def solve(problem, out, seconds, progress, solver="clarabel"):
    import numpy as np
    from scipy.sparse import csc_matrix
    A, b, layout, sizes = pack(problem, solver)
    conic.write_json(out/"numerical_cone_layout.json", layout)
    progress("sparse_conic_solver_start", dict(solver=solver, variables=A.shape[1], rows=A.shape[0],
        sparse_nonzeros=A.nnz, cones=sizes, order_two_PSD_cones_converted_to_SOC=True))
    started = time.monotonic()
    if solver == "clarabel":
        import clarabel
        cones = ([clarabel.NonnegativeConeT(sizes["nonnegative"])] if sizes["nonnegative"] else [])
        cones += [clarabel.SecondOrderConeT(3) for _ in sizes["soc"]]
        cones += [clarabel.PSDTriangleConeT(r) for r in sizes["psd"]]
        settings = clarabel.DefaultSettings()
        settings.verbose = True; settings.max_iter = 150; settings.time_limit = seconds
        settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-9
        settings.direct_solve_method = "qdldl"
        if hasattr(settings, "max_threads"): settings.max_threads = 1
        with native_output(out/"solver.stdout"):
            instance = clarabel.DefaultSolver(csc_matrix((A.shape[1], A.shape[1])),
                np.zeros(A.shape[1]), A, b, cones, settings)
            result = instance.solve()
        info = dict(solver="Clarabel", version=getattr(clarabel, "__version__", None), status=str(result.status),
            iterations=result.iterations, solve_time=result.solve_time,
            primal_residual=result.r_prim, dual_residual=result.r_dual,
            objective=result.obj_val, dual_objective=result.obj_val_dual)
        x, z = result.x, result.z
    elif solver == "scs":
        import scs
        cone = dict(l=sizes["nonnegative"], q=sizes["soc"], s=sizes["psd"])
        with native_output(out/"solver.stdout"):
            result = scs.solve(dict(A=A, b=b, c=np.zeros(A.shape[1])), cone,
                verbose=True, time_limit_secs=seconds, max_iters=200000,
                eps_abs=1e-8, eps_rel=1e-8, normalize=True)
        info = dict(solver="SCS", version=scs.__version__, **result["info"])
        x, z = result["x"], result["y"]
    else: raise ValueError("Unknown conic backend")
    info["wall_seconds"] = time.monotonic()-started
    info["cone_sizes"] = sizes
    point = None if x is None or not np.all(np.isfinite(x)) else {
        v: F(str(float(value))) for v, value in zip(problem["variables"], x, strict=True)}
    return point, unpack_dual(problem, layout, z, solver), info
