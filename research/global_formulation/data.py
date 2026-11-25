"""The complete intrinsic component schemas, with symbolic unbounded tails.

The fixed corner tables come from the manuscript's retained literal source.
Small connected components are fixed input bins, not selected vertex/face
diagrams. No list of forests or marked stars is generated.
"""
from collections import deque
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LITERAL = "scripts/finite57_literal_weights_audit_20261024.py"
LITERAL_SHA256 = "f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0"
FINITE_TAIL = {"A": (21, 23), "B": (21, 26), "D": (22, 30)}
AFFINE_MINIMUM = 25
SAME_TYPES = (("A", 13), ("A", 14), ("A", 15), ("A", 15), ("D", 15), ("D", 15),
    ("A", 16), ("A", 16), ("D", 16), ("D", 16), ("D", 14), ("D", 14),
    ("D", 13), ("D", 13), ("D", 28))
SAME_BASELINES = tuple(F(1, 2) if j in (12, 13) else F(1) if j in (16, 17) else F(0) for j in range(4, 19))
SAME_DISTANCE = (12, 13, 13, 13, 12, 12, 11, 11, 14, 14, 12, 12, 11, 11, 12)
CROSS = ((19, ("A", 4), ("A", 12), 4, "has_a1"),
    (20, ("A", 2), ("A", 4), 4, "a1_a10"),
    (21, ("A", 2), ("D", 5), 2, "a1_a9"),
    (22, ("A", 1), ("D", 5), 1, "a11_or_a9_a2"),
    (23, ("A", 1), ("D", 5), 1, "a9_a2"),
    (24, ("A", 12), ("D", 16), 2, "two15"))


def finite_types():
    result = [(kind, r) for kind, lower in (("A", 1), ("B", 2), ("D", 4))
        for r in range(lower, FINITE_TAIL[kind][0])]
    return tuple(sorted((*result, ("D", 28), ("E", 6), ("E", 7), ("E", 8),
        ("F", 4), ("H", 3), ("H", 4), ("I", 2))))


def affine_types():
    return tuple(sorted([(k, r) for k, lower in (("A", 1), ("B", 3), ("C", 2), ("D", 4))
        for r in range(lower, AFFINE_MINIMUM)]+[("E", 6), ("E", 7), ("E", 8), ("F", 4), ("G", 2)]))


def literal_module():
    path = ROOT/LITERAL
    if hashlib.sha256(path.read_bytes()).hexdigest() != LITERAL_SHA256:
        raise ValueError("The manuscript's retained intrinsic table source changed")
    spec = importlib.util.spec_from_file_location("_global_coxeter_intrinsic_tables", path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def graph(kind, rank, affine=False):
    order = rank+int(affine)
    adjacency = [set() for _ in range(order)]
    def edge(a, b):
        adjacency[a].add(b); adjacency[b].add(a)
    def arms(lengths):
        following = 1
        for length in lengths:
            previous = 0
            for _ in range(length):
                edge(previous, following); previous, following = following, following+1
        if following != order:
            raise ValueError("Wrong exceptional connected-template order")
    if kind == "E":
        arms(({6: (2, 2, 2), 7: (3, 3, 1), 8: (5, 2, 1)} if affine else
              {6: (2, 2, 1), 7: (3, 2, 1), 8: (4, 2, 1)})[rank])
    elif not affine and kind == "D":
        for i in range(order-3): edge(i, i+1)
        edge(order-3, order-2); edge(order-3, order-1)
    elif affine and kind == "A" and rank >= 2:
        for i in range(order): edge(i, (i+1) % order)
    elif affine and kind == "B":
        for i in range(order-3): edge(i, i+1)
        edge(0, order-2); edge(0, order-1)
    elif affine and kind == "D":
        for i in range(order-5): edge(i, i+1)
        for leaf in (order-4, order-3): edge(0, leaf)
        for leaf in (order-2, order-1): edge(order-5, leaf)
    else:
        for i in range(order-1): edge(i, i+1)
    return adjacency


def prokhorov_statistics(adjacency, deleted=None):
    W, J = F(0), F(0)
    for start in range(len(adjacency)):
        if start == deleted: continue
        distances, queue = {start: 0}, deque([start])
        while queue:
            a = queue.popleft()
            for b in adjacency[a]:
                if b != deleted and b not in distances:
                    distances[b] = distances[a]+1; queue.append(b)
        for end, distance in distances.items():
            if start < end:
                weight = F(1) if distance <= 7 else F(1, 3) if distance <= 15 else F(0)
                W += weight; J += (distance+1)*weight
    return W, J


def intrinsic_tables():
    source = literal_module()
    finite, affine, cross = {}, {}, {}
    for t in finite_types():
        S = sum((source.finite_same_weight(t, i, j) for i, j in combinations(range(t[1]), 2)), F(0))
        W, J = prokhorov_statistics(graph(*t))
        finite[t] = dict(ordinary=S, W=W, J=J)
    for t in affine_types():
        order = t[1]+1
        S = sum((source.affine_weight(t[0], order, i, j) for i, j in combinations(range(order), 2)), F(0))
        adjacency = graph(*t, affine=True)
        deletion = sum((prokhorov_statistics(adjacency, deleted=i)[0] for i in range(order)), F(0))
        affine[t] = dict(ordinary_per_order=S/order, deletion_per_order=deletion/order)
    small = [t for t in finite_types() if t[1] <= 17]
    for a, b in combinations(small, 2):
        if min(a[1], b[1]) > 8: continue
        value = sum((source.finite_cross_weight(a, i, b, j) for i in range(a[1]) for j in range(b[1])), F(0))
        if value: cross[a, b] = value
    return dict(finite=finite, affine=affine, ordinary_cross=cross,
        source=LITERAL, source_sha256=LITERAL_SHA256)


def atoms():
    result = []
    for j in range(26):
        v = [F(0)]*26; v[0] = 1
        if j: v[j] = 1
        if 4 <= j <= 18: v[1] = SAME_BASELINES[j-4]
        if j == 23: v[22] = 1
        result.append(tuple(v))
    return tuple(result)
