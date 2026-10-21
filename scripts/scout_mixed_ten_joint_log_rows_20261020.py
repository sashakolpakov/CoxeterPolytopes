#!/usr/bin/env python3
"""Prepared bounded numerical relaxation, NOT an exact nonexistence proof.

Root execution only, through the shared resource/lock wrapper. Importing this
file does not import scipy or run an optimizer. --plan is solver-free.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import comb, isqrt
from pathlib import Path
import json
import signal

ROOT = Path(__file__).resolve().parents[1]
INPUTS = [
    "results/right_angled_ten_joint_log_transport_20261020.md",
    "results/right_angled_r7_four_direction_transport_research_20261020.md",
    "results/right_angled_three_face_upper_d_two_excess_20261020.md",
    "results/right_angled_finite_vertex_seed_and_face_cover_20261020.md",
    "results/ideal_right_angled_dimension4_proof_20261020.md",
    "manuscripts/right_angled_dimension10.tex",
    "results/right_angled_mixed_incidence_feasibility_research_20261018.md",
]


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def sqrt_upper(x):
    scale = 10**12
    n = isqrt(x.numerator * scale * scale // x.denominator)
    assert F(n*n, scale*scale) <= x < F((n+1)**2, scale*scale)
    return F(n+1, scale)


def log_bounds(x):
    assert x > 0
    y = (x-1)/(x+1)
    partial = 2*sum((y**(2*j+1)/F(2*j+1) for j in range(48)), F(0))
    tail = 2*abs(y)**97/(97*(1-y*y))
    return partial-tail, partial+tail


def log_coefficients():
    left, right = [F(0)]*11, [F(0)]*11
    for k in range(3, 10):
        mu = F(11, 10)-F(3, 10)*sqrt_upper(F(10-k, k))
        left[k] = -choose(k, 3)*log_bounds(mu)[0]
    left[9], left[10] = F(0), -F(5815, 616)
    for m in range(4, 10):
        mu = F(7, 5)-F(4, 5)*sqrt_upper(F(10-m, m))
        right[m] = -choose(m, 4)*log_bounds(mu)[0]
    right[8], right[10] = F(0), -F(1070, 21)
    return left, right, log_bounds(F(25, 16))[1], log_bounds(F(12, 5))[1]


class Rows:
    def __init__(self):
        self.names, self.bounds, self.eq, self.ub = [], [], [], []

    def var(self, name, lower=0, upper=None):
        self.names.append(name)
        self.bounds.append((lower, upper))
        return name

    def add(self, target, name, terms, rhs=0):
        row = {}
        for key, value in terms:
            row[key] = row.get(key, F(0))+F(value)
        target.append((name, row, F(rhs)))

    def equal(self, name, terms, rhs=0):
        self.add(self.eq, name, terms, rhs)

    def le(self, name, terms, rhs=0):
        self.add(self.ub, name, terms, rhs)


def model():
    r = Rows()
    fk = [r.var(f"f{k}") for k in range(11)]
    nk = [r.var(f"N{k}") for k in range(10)]
    zk = [r.var(f"z{k}") for k in range(11)]
    bk = [r.var(f"b4_{k}") for k in range(11)]
    for name in ("E", "D", "H", "Z", "c", "B", "R", "X", "C", "Q", "g", "Pbig", "Sbig", "Fbig", "u"):
        r.var(name, None if name in ("D", "R") else 0)
    for j in range(1, 10):
        r.var(f"a{j}")
    r.equal("finite_normalization", [(v, 1) for v in fk], 1)
    r.equal("D_vertex", [("D", 1)]+[(fk[k], -choose(k, 2)*(10-k)+3*choose(10-k, 3)) for k in range(11)])
    r.equal("H_missing", [("H", 1), ("B", 2)]+[(fk[k], -choose(k, 3)) for k in range(11)])
    r.equal("Z_sum", [("Z", 1)]+[(v, -1) for v in zk])
    r.equal("R7_centers", [("C", 1), ("Z", 1)]+[(fk[k], -choose(10-k, 3)) for k in range(11)])
    for k in range(11):
        m = 10-k
        r.le(f"z_capacity_{k}", [(zk[k], 1), (fk[k], -choose(m, 3))])
        r.le(f"b4_capacity_{k}", [(bk[k], 1), (fk[k], -choose(m, 4))])
        r.le(f"b4_union_{k}", [(bk[k], 1), (zk[k], -max(m-3, 0))])
    for k in range(10):
        r.le(f"irregular_capacity_{k}", [(nk[k], 1), (fk[k], -(10-k))])
    for k in range(5):
        j = 9-k
        r.equal(f"FF_pair_{k}_{j}", [(fk[k], 10-k), (nk[k], -1), (fk[j], -(k+1)), (nk[j], 1)])
    r.le("D_lower", [("D", -1), ("E", -30)])
    r.le("D_upper", [("D", 1), ("E", -2)])
    r.le("H_budget", [("H", 1), ("E", -2)])
    r.le("Z_budget", [("Z", 4), ("E", -10), ("D", 1)])
    r.le("FF_union_coarse", [(v, 1) for v in nk]+[("E", -30)])
    r.le("FF_union_refined", [(v, 1) for v in nk]+[("Pbig", -2), ("Z", 2), ("D", 1), ("H", 2)])
    r.le("big_A2_nonnegative", [("D", -1), ("Z", -3)])
    r.le("big_A1_nonnegative", [("Z", 4), ("D", 1), ("H", 1), ("Pbig", -1)])
    r.le("big_FI_capacity", [("Pbig", 1), ("Z", 2), ("D", 1), ("H", 2), ("Sbig", -4)])
    L, M, penalty3, penalty4 = log_coefficients()
    r.le("P3_log", [(fk[k], -L[k]) for k in range(11)]+[("H", -penalty3)])
    r.le("R7_log", [(fk[k], -M[10-k]) for k in range(11)]+[(v, -penalty4) for v in bk])
    r.equal("three_face_partition", [("a3", 1), ("B", -1), ("Q", -1), ("C", -1), ("g", -1)])
    r.equal("finite_three_incidence", [("B", 2), ("Q", 8), ("C", 4), ("Pbig", 1)], 120)
    r.equal("cusp_three_incidence", [("B", 3), ("Q", 2), ("C", 3), ("Sbig", 1), ("c", -4608)])
    r.equal("facet_three_incidence", [("B", 6), ("Q", 8), ("C", 7), ("Fbig", 1), ("a2", -8)])
    r.equal("big_excess", [("E", 1), ("Fbig", -1), ("Sbig", -1), ("g", 10)])
    r.equal("big_Euler", [("Pbig", 1), ("Fbig", -2), ("Sbig", 2), ("g", 4)])
    r.equal("signed_R", [("R", 1), ("C", -1), ("g", -8), ("Fbig", 1)])
    r.equal("exceptional_corners", [("X", 1), ("C", -4), ("Pbig", -1)])
    r.le("X_budget", [("X", 1), ("E", -18), ("R", -4)])
    r.le("big_count", [("g", 1), ("E", -1)])
    r.le("big_facets", [("g", 7), ("Fbig", -1)])
    r.le("big_finite_vertices", [("Pbig", 1), ("E", -10)])
    # Truncated simple ten-polytope: normalize 1/f as the nonnegative u.
    poly = [{} for _ in range(11)]
    for degree in range(11):
        poly[degree]["u"] = F(choose(10, degree)*(-1)**(10-degree))
    for j in range(10):
        base = [(f"a{j}", 1)] if j else [(v, 1) for v in fk]
        base += [("c", 2**(9-j)*choose(9, j))]
        for degree in range(j+1):
            mult = choose(j, degree)*(-1)**(j-degree)
            for key, value in base:
                poly[degree][key] = poly[degree].get(key, F(0))+mult*value
    for j in range(5):
        r.equal(f"DS_{j}", list(poly[j].items())+[(v, -c) for v, c in poly[10-j].items()])
    for j in range(11):
        r.le(f"h_nonnegative_{j}", [(v, -c) for v, c in poly[j].items()])
    for j in range(1, 6):
        r.le(f"g_nonnegative_{j}", list(poly[j-1].items())+[(v, -c) for v, c in poly[j].items()])
    minima = {3: 6, 4: 10, 5: 16, 6: 27, 7: 46, 8: 82, 9: 218, 10: 1923}
    for j, minimum in minima.items():
        high = f"a{j}" if j < 10 else "u"
        r.le(f"facet_min_{j}", [(high, minimum), (f"a{j-1}", -(11-j))])
    nu_minima = {3: 9, 4: 15, 5: 26, 6: 45, 7: 81, 8: 217, 9: 1922, 10: 183967}
    for j, minimum in nu_minima.items():
        high = f"a{j}" if j < 10 else "u"
        r.le(f"facet_cusp_min_{j}", [(high, minimum), (f"a{j-1}", -(11-j)), ("c", -2**(10-j)*choose(9, j-1))])
    finite_minima = {5: 6, 6: 27, 7: 178, 8: 1825, 9: 44206, 10: 8500814}
    for j, minimum in finite_minima.items():
        high = f"a{j}" if j < 10 else "u"
        r.le(f"finite_face_cover_{j}", [(high, minimum)], choose(10, j))
    r.equal("edge_incidence", [("a1", 2), ("c", -512)], 10)
    return r


def digest(path):
    h = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", action="store_true")
    parser.add_argument("--write", type=Path, help="exclusive numerical scout ledger; never a proof certificate")
    args = parser.parse_args()
    r = model()
    pins = {p: digest(ROOT/p) for p in INPUTS+[str(Path(__file__).relative_to(ROOT))]}
    if args.plan:
        print(json.dumps({"status": "SOURCE_PLAN_ONLY", "variables": len(r.names), "equalities": len(r.eq), "inequalities": len(r.ub), "owns_lock": False, "inputs": pins}, sort_keys=True))
        return
    if args.write is None:
        parser.error("use --plan or an exclusive --write path; no implicit optimizer execution")
    if args.write.exists():
        parser.error("refusing to overwrite an existing ledger")
    signal.alarm(45)
    from scipy.optimize import linprog
    index = {v: j for j, v in enumerate(r.names)}
    def arrays(rows):
        return ([[float(row.get(v, 0)) for v in r.names] for _, row, _ in rows], [float(rhs) for _, _, rhs in rows])
    aub, bub = arrays(r.ub)
    aeq, beq = arrays(r.eq)
    objective = [0.0]*len(r.names)
    objective[index["E"]] = 1.0
    result = linprog(objective, A_ub=aub, b_ub=bub, A_eq=aeq, b_eq=beq, bounds=r.bounds, method="highs", options={"time_limit": 25, "threads": 1})
    out = {"status": "NUMERICAL_RELAXATION_ONLY", "solver_status": int(result.status), "message": result.message, "objective": "min E / finite_vertices", "variables": len(r.names), "equalities": len(r.eq), "inequalities": len(r.ub), "pins": pins}
    if result.x is not None:
        out["candidate"] = {v: float(result.x[j]) for j, v in enumerate(r.names)}
        out["maximum_float_equality_residual"] = max(abs(sum(float(c)*result.x[index[v]] for v, c in row.items())-float(rhs)) for _, row, rhs in r.eq)
        out["maximum_float_inequality_violation"] = max(0.0, max(sum(float(c)*result.x[index[v]] for v, c in row.items())-float(rhs) for _, row, rhs in r.ub))
    assert pins == {p: digest(ROOT/p) for p in pins}, "source changed during scout"
    payload = json.dumps(out, indent=2, sort_keys=True)+"\n"
    assert len(payload.encode()) < 1_048_576
    with args.write.open("x") as stream:
        stream.write(payload)
    print(json.dumps({key: out[key] for key in ("status", "solver_status", "message", "variables", "equalities", "inequalities")}))


if __name__ == "__main__":
    main()
