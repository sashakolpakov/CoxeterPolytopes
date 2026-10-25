"""Independent transcription of finite_volume57.tex Appendix weights.

SOURCE PREPARATION ONLY: not executed, tested, or accepted. This module imports
no repository worker and reads no certificate. A separate bounded checker must
compare it with frozen actual aggregates and documented upper relaxations.

Finite types are (family, rank): A, B, D; E6/E7/E8; F4; H3/H4;
I2 denotes an arbitrary I2(m), m>=5. Its weights are independent of m.
E6 is numbered with center0, arms(1),(2,3),(4,5), long ends3,5.
All other classical numbering is exactly that of the printed manuscript.
"""
from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations

Type = tuple[str, int]
Mark = tuple[int, int]  # component index, node index
ZERO, ONE, HALF = F(0), F(1), F(1, 2)
TEX_SHA256 = "c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0"
E6_LONG_ENDS = frozenset((3, 5))


def _pair(a: int, b: int) -> tuple[int, int]:
    if a == b:
        raise ValueError("two distinct nodes required")
    return (a, b) if a < b else (b, a)


def _check_type(t: Type) -> None:
    k, n = t
    valid = ((k == "A" and n >= 1) or (k == "B" and n >= 2)
             or (k == "D" and n >= 4) or (k == "E" and n in (6, 7, 8))
             or (k == "F" and n == 4) or (k == "H" and n in (3, 4))
             or (k == "I" and n == 2))
    if not valid:
        raise ValueError(f"unsupported finite type {t}")


def _check_mark(t: Type, a: int) -> None:
    _check_type(t)
    if not 0 <= a < t[1]:
        raise ValueError("mark outside component")


def classical_graph(kind: str, order: int, affine: bool = False):
    """Adjacency with native Coxeter labels, not a Gram matrix."""
    g = [dict() for _ in range(order)]
    def edge(a, b, label=3):
        g[a][b] = label
        g[b][a] = label
    if not affine:
        if kind in ("A", "B"):
            for i in range(order - 1):
                edge(i, i + 1, 4 if kind == "B" and i == order - 2 else 3)
        elif kind == "D":
            for i in range(order - 3):
                edge(i, i + 1)
            edge(order - 3, order - 2)
            edge(order - 3, order - 1)
        else:
            raise ValueError("classical family required")
    elif kind == "A":
        if order < 3:
            raise ValueError("cycle graph requires order>=3")
        for i in range(order):
            edge(i, (i + 1) % order)
    elif kind == "B":
        for i in range(order - 3):
            edge(i, i + 1, 4 if i == order - 4 else 3)
        edge(0, order - 2)
        edge(0, order - 1)
    elif kind == "C":
        for i in range(order - 1):
            edge(i, i + 1, 4 if i in (0, order - 2) else 3)
    elif kind == "D":
        for i in range(order - 5):
            edge(i, i + 1)
        for leaf in (order - 4, order - 3):
            edge(0, leaf)
        for leaf in (order - 2, order - 1):
            edge(order - 5, leaf)
    else:
        raise ValueError("classical affine family required")
    return g


def _distances(g, start):
    distances = {start: 0}
    queue = deque((start,))
    while queue:
        i = queue.popleft()
        for j in g[i]:
            if j not in distances:
                distances[j] = distances[i] + 1
                queue.append(j)
    if len(distances) != len(g):
        raise ValueError("connected graph required")
    return distances


def _features(g, finite):
    return {i for i, row in enumerate(g)
            if (len(row) != 2 if finite else len(row) >= 3)
            or any(label != 3 for label in row.values())}


FINITE_REPLACEMENTS = {}
def _finite(kind, rank, value, pairs):
    table = FINITE_REPLACEMENTS.setdefault((kind, rank), {})
    for a, b in pairs:
        p = _pair(a, b)
        if p in table:
            raise ValueError("duplicate finite replacement")
        table[p] = value

_finite("D",12,HALF,[(0,10),(0,11)])
_finite("A",13,ZERO,[(0,12)])
_finite("B",13,ZERO,[(0,12)])
_finite("A",14,ZERO,[(0,12),(1,13),(0,13)])
_finite("B",14,ZERO,[(0,13),(1,13),(2,13),(3,13)])
_finite("D",14,ZERO,[(0,12),(0,13)])
_finite("D",14,HALF,[(0,11)])
_finite("A",15,ZERO,[(0,12),(0,14),(2,14),(0,13),(1,14)])
_finite("B",15,ZERO,[(0,12),(0,13),(0,14),(1,13),(1,14),(2,14),(3,14),(4,14)])
_finite("D",15,ZERO,[(0,11),(1,13),(1,14)])
_finite("D",15,HALF,[(0,12),(0,13),(0,14),(1,12)])
_finite("A",16,ZERO,[(0,12),(0,14),(0,15),(1,14),(1,15),(3,15),(1,12),(3,14)])
_finite("A",16,HALF,[(0,13),(1,13),(2,14),(2,15)])
_finite("A",16,F(1,8),[(0,11),(4,15)])
_finite("B",16,ZERO,[(0,11),(0,12),(0,13),(0,14),(0,15),(1,12),(1,13),(1,14),
    (1,15),(2,13),(2,14),(2,15),(3,14),(3,15),(4,15),(5,15)])
_finite("D",16,ZERO,[(0,10),(0,11),(0,12),(0,13),(1,11),(2,12),(2,14),(2,15)])
_finite("D",16,HALF,[(0,14),(0,15),(1,12),(1,13),(1,14),(1,15),(2,13)])
_finite("A",17,ZERO,[(0,10),(0,11),(0,12),(0,14),(0,15),(1,10),(1,11),(1,13),
    (1,14),(1,16),(2,12),(2,14),(2,15),(2,16),(3,15),(4,14),(4,16),
    (5,15),(5,16),(6,15),(6,16)])
_finite("A",17,HALF,[(0,16),(2,13),(3,14)])
_finite("A",17,F(7,8),[(1,12),(4,15)])
_finite("D",17,ZERO,[(0,10),(0,11),(0,12),(0,13),(1,10),(1,11),(1,12),(1,13),
    (1,14),(1,15),(1,16),(2,12),(2,14),(3,15),(3,16),(4,14),(4,15),(4,16)])
_finite("D",17,HALF,[(0,14),(2,13),(2,15),(2,16),(3,13),(3,14)])
_finite("D",17,F(1,4),[(0,15),(0,16)])


def finite_same_weight(t: Type, a: int, b: int, preliminary=False) -> F:
    _check_mark(t, a)
    _check_mark(t, b)
    p = _pair(a, b)
    kind, rank = t
    if preliminary and t == ("B",13) and p == (0,12):
        return ONE
    if rank <= 11:
        return ONE
    if rank <= 17 and t != ("B",17):
        return FINITE_REPLACEMENTS.get(t, {}).get(p, ONE)
    g = classical_graph(kind, rank)
    da, db = _distances(g,a), _distances(g,b)
    distance = da[b]
    leaf = any(len(g[m]) == 1 for m in (a,b))
    paid = distance <= 8 or (distance == 9 and leaf)
    if t == ("B",17):
        return F(paid or p in ((2,11),(3,12),(4,13)))
    if distance in (6,7,8):
        features = _features(g, finite=True)
        if min(min(da[f] for f in features), min(db[f] for f in features)) >= 6:
            paid = False
    if not paid:
        return ZERO
    if kind == "D" and rank >= 20 and distance == 9 and {a,b} & {rank-2,rank-1}:
        return HALF
    return ONE


def _ends(rank): return {0, rank-1}
def _positions(rank, k): return {k-1, rank-k}


def finite_cross_weight(a: Type, i: int, b: Type, j: int) -> F:
    """Ordinary cross entries: broad rows first, later rows replace."""
    _check_mark(a,i)
    _check_mark(b,j)
    value = ZERO
    for sm,m,la,n in ((a,i,b,j),(b,j,a,i)):
        sk,sr = sm
        lk,lr = la
        if sr <= 5:
            if lk == "D":
                if lr == 10 and n in {lr-2,lr-1}: value = HALF
                if 11 <= lr <= 14 and n in {lr-3,lr-2,lr-1}: value = HALF
                if lr == 15 and ((sr<=2 and n in {12,13,14}) or (3<=sr<=5 and n in {13,14})): value = HALF
                if lr == 16 and sm == ("A",1) and n in {13,14,15}: value = HALF
            if la in (("A",13),("A",15)) and n == (lr-1)//2: value = HALF
            if la == ("A",17) and n in {5,11}: value = HALF
        if la == ("D",11) and n in {9,10}:
            if sm in (("A",1),("A",2)): value = ONE
            if sm == ("A",4) and m in _ends(4): value = ONE
            if sm == ("A",5) and m in _positions(5,2): value = ONE
            if sm == ("E",6) and m in E6_LONG_ENDS: value = ONE
            if sm == ("D",5) and m == 0: value = ONE
        if la == ("D",12):
            if sm == ("A",1) and n in {10,11}: value = ONE
            if sm == ("D",5) and m in {3,4} and n in {9,10,11}: value = ONE
        if la == ("D",13):
            if sm == ("A",2) and n in {11,12}: value = ONE
            if sm == ("A",2) and n == 10: value = F(23,64)
            if (sm == ("A",1) or (sm == ("A",4) and m in _ends(4))) and n == 9: value = HALF
        if la == ("D",14):
            if sm in (("A",1),("A",2)) and n == 10: value = HALF
            if sm == ("A",1) and n == 9: value = F(1,16)
        if la == ("D",15) and sm == ("A",1) and n == 11: value = F(1,16)
        if la == ("D",10) and n == 7 and (sm in (("A",1),("A",2)) or (sm == ("A",3) and m == 1)): value = ONE
        if la == ("B",9) and n == 8 and (sm in (("A",1),("A",2)) or (sm in (("A",3),("A",4)) and m in _ends(sr))): value = HALF
        if la == ("B",10) and n == 9:
            if sm == ("A",1): value = ONE
            if sm == ("A",2): value = HALF
        if la in (("B",11),("B",12)) and n == 9 and sm == ("A",1): value = HALF
        if la == ("D",9) and n in {7,8} and sk == "A" and 1 <= sr <= 8 and m in _ends(sr): value = F(5,32)
        if la == ("A",9) and n == 4 and sk == "A" and 1 <= sr <= 4 and m in _ends(sr): value = HALF
        if la == ("A",11):
            if n == 5:
                if sm == ("A",1): value = ONE
                if sm == ("A",2): value = HALF
            if n in _ends(11) and (sm in (("A",1),("A",2)) or (sm == ("A",4) and m in _ends(4))): value = HALF
            if n in _positions(11,2) and ((sm == ("E",6) and m in E6_LONG_ENDS) or (sm == ("D",5) and m == 0)): value = ONE
        if la == ("A",12):
            if sm == ("A",1):
                if n in {5,6}: value = ONE
                if n in _ends(12): value = F(1,4)
            if sm == ("A",4) and m in _ends(4) and n in (_ends(12)|_positions(12,3)): value = ONE
            if sm == ("D",5):
                if m == 0 and n in _positions(12,2): value = ONE
                if m == 0 and n in _ends(12): value = F(3,8)
                if m in {3,4} and n in _positions(12,3): value = ONE
        if la == ("A",13):
            if n == 6:
                if sm == ("A",1): value = F(41,64)
                if sm == ("A",2): value = F(23,64)
            if sm == ("A",1):
                if n in _ends(13): value = HALF
                if n in _positions(13,2): value = F(1,4)
            if sm == ("A",2) and n in _ends(13): value = ONE
            if sm == ("A",4) and n in _positions(13,2): value = ONE if m in _ends(4) else HALF
        if la == ("A",14):
            if sm == ("A",2) and n in (_positions(14,2)|_positions(14,3)): value = HALF
            if sm == ("A",1) and n in (_positions(14,3)|_positions(14,4)): value = HALF
        if la == ("A",15) and sm == ("A",1) and n in _positions(15,3): value = HALF
        if la == ("A",16) and sm == ("A",1):
            if n in _positions(16,4): value = HALF
            if n in _ends(16): value = F(1,4)
    return value


def finite_corner_weight(forest, first: Mark, second: Mark, *, preliminary=False, final=False):
    """One actual finite marked forest. final=True is dimension58 only."""
    ca,a = first
    cb,b = second
    if not (0<=ca<len(forest) and 0<=cb<len(forest)):
        raise ValueError("component index outside forest")
    ta,tb = forest[ca],forest[cb]
    _check_mark(ta,a)
    _check_mark(tb,b)
    if final and (preliminary or sum(t[1] for t in forest) != 58):
        raise ValueError("final rule requires an actual rank58 forest")
    witnesses = Counter(t for index,t in enumerate(forest) if index not in {ca,cb})
    rank15 = witnesses[("A",15)] + witnesses[("D",15)]
    has1 = bool(witnesses[("A",1)])
    if ca == cb:
        p = _pair(a,b)
        value = finite_same_weight(ta,a,b,preliminary)
        if has1:
            if ta == ("A",13) and p == (0,12): value = ONE
            if ta == ("A",14) and p == (0,13): value = HALF
            if ta == ("A",15) and p in {(0,13),(1,14)}: value = F(1,4)
            if ta == ("D",15) and p in {(1,13),(1,14)}: value = F(1,4)
            if ta == ("A",16) and p in {(1,12),(3,14)}: value = F(5,8)
            if ta == ("D",16) and p in {(0,14),(0,15)}: value = ONE
        if ta == ("D",14) and p in {(0,12),(0,13)} and any(witnesses[("A",r)] for r in (1,2,3)): value = HALF
        if final and ta == ("D",13) and p in {(0,11),(0,12)} and rank15 >= 3: value = HALF
        if final and ta == ("D",28) and p == (0,12) and rank15 >= 2: value += HALF
        return value
    value = finite_cross_weight(ta,a,tb,b)
    addition = ZERO
    for sm,m,la,n in ((ta,a,tb,b),(tb,b,ta,a)):
        if sm == ("A",4) and m in _ends(4) and la == ("A",12) and n in _positions(12,2) and has1: addition += F(1,8)
        if sm == ("A",2) and la == ("A",4) and n in _ends(4) and has1 and witnesses[("A",10)]: addition += F(1,4)
        if sm == ("A",2) and la == ("D",5) and n == 0 and has1 and witnesses[("A",9)]: addition += HALF
        if sm == ("A",1) and la == ("D",5) and n == 0:
            joint = bool(witnesses[("A",9)] and witnesses[("A",2)])
            if witnesses[("A",11)] or joint: addition += F(3,8)
            if joint: addition += F(5,8)
        if final and sm == ("A",12) and m in _ends(12) and la == ("D",16) and n == 0 and rank15 >= 2: addition += HALF
    if ta == tb == ("A",2) and witnesses[("A",10)] and witnesses[("A",1)] >= 2: addition += F(9,32)
    if addition and value:
        raise ValueError("printed contextual additions require ordinary zero")
    return value + addition


def context_sum(forest, *, final=False):
    """Closed symbolic sum of all context increments, in unscaled units."""
    c = Counter(forest)
    a = lambda n: c[("A",n)]
    d = lambda n: c[("D",n)]
    z = 36*a(2)*(a(2)-1)*int(a(1)>=2 and a(10)>0)
    z += 64*a(2)*d(5)*int(a(1)>0 and a(9)>0)
    z += 24*a(1)*d(5)*int(a(11)>0 or (a(9)>0 and a(2)>0))
    z += 40*a(1)*d(5)*int(a(9)>0 and a(2)>0)
    z += 32*(a(15)+d(15))*int(a(1)>0)
    z += 64*d(14)*int(a(1)+a(2)+a(3)>0)
    z += int(a(1)>0)*(80*a(16)+32*a(14)+64*d(16)+64*a(13)
        +32*a(4)*a(12)+64*a(2)*a(4)*int(a(10)>0))
    if final:
        if sum(t[1] for t in forest) != 58:
            raise ValueError("final guards are dimension58 only")
        k = a(15)+d(15)
        z += -64*d(13)*int(k>=3)+(32*d(28)+64*a(12)*d(16))*int(k>=2)
    return F(z,64)


def same_component_sum(t, *, preliminary=False):
    return sum((finite_same_weight(t,a,b,preliminary) for a,b in combinations(range(t[1]),2)),ZERO)


def ordinary_cross_sum(a,b):
    return sum((finite_cross_weight(a,i,b,j) for i in range(a[1]) for j in range(b[1])),ZERO)


def forest_sum(forest, *, preliminary=False, final=False):
    """Literal aggregate, distinct copies included even when types agree."""
    if preliminary and final:
        raise ValueError("preliminary and final stages are distinct")
    z = sum((same_component_sum(t,preliminary=preliminary) for t in forest),ZERO)
    z += sum((ordinary_cross_sum(a,b) for a,b in combinations(forest,2)),ZERO)
    return z + context_sum(forest,final=final)


AFFINE_REPLACEMENTS = {}
def _affine(kind, order, value, pairs):
    table = AFFINE_REPLACEMENTS.setdefault((kind,order),{})
    for a,b in pairs:
        p = _pair(a,b)
        if p in table:
            raise ValueError("duplicate affine replacement")
        table[p] = value
def _product(left,right): return [(a,b) for a in left for b in right]

_affine("B",10,HALF,[(7,8),(7,9)])
_affine("C",10,ZERO,[(0,9)])
_affine("D",10,F(27,32),_product((6,7),(8,9)))
_affine("B",11,ZERO,[(8,9),(8,10)])
_affine("C",11,ZERO,[(0,9),(0,10),(1,10)])
_affine("D",11,HALF,_product((7,8),(9,10)))
_affine("B",12,ZERO,[(0,9),(9,10),(9,11)])
_affine("B",12,HALF,[(8,10),(8,11)])
_affine("C",12,ZERO,[(0,10),(0,11),(1,10),(1,11)])
_affine("D",12,ZERO,_product((8,9),(10,11)))
_affine("B",13,ZERO,[(0,10),(9,11),(9,12),(10,11),(10,12)])
_affine("B",13,HALF,[(1,10),(8,11),(8,12)])
_affine("C",13,ZERO,[(0,10),(0,11),(0,12),(1,11),(1,12),(2,10),(2,12)])
_affine("D",13,ZERO,[(0,8)]+_product((9,10),(11,12)))
_affine("D",13,HALF,[(0,11),(0,12),(8,9),(8,10)])
_affine("B",14,ZERO,[(0,10),(0,11),(1,10),(1,11),(2,11),(9,12),(9,13),(10,12),(10,13),(11,12),(11,13)])
_affine("C",14,ZERO,[(0,10),(0,11),(0,12),(0,13),(1,11),(1,12),(1,13),(2,11),(2,12),(2,13),(3,13)])
_affine("D",14,ZERO,[(0,9),(0,12),(0,13),(9,10),(9,11)]+_product((10,11),(12,13)))
_affine("D",14,HALF,[(1,12),(1,13),(8,10),(8,11)])
_affine("B",15,ZERO,[(0,11),(0,12),(1,11),(1,12),(2,12),(3,12)]+_product((9,10,11,12),(13,14)))
_affine("D",15,ZERO,[(0,10),(0,13),(0,14),(1,13),(1,14),(9,11),(9,12),(10,11),(10,12)]+_product((11,12),(13,14)))
_affine("D",15,F(15,16),[(2,13),(2,14),(8,11),(8,12)])
_affine("D",16,ZERO,[(0,11),(0,14),(0,15),(2,14),(2,15),(9,12),(9,13),(11,12),(11,13)]+_product((12,13),(14,15)))
_affine("D",16,HALF,[(0,9),(0,10),(1,11),(2,11)])
_affine("D",16,F(7,16),[(1,14),(1,15),(10,12),(10,13)])
_affine("D",17,ZERO,[(0,10),(0,12),(0,15),(0,16),(1,11),(2,12),(3,15),(3,16),(9,13),(9,14),(12,13),(12,14)]+_product((13,14),(15,16)))
_affine("D",17,HALF,[(0,11),(1,12),(1,15),(1,16),(11,13),(11,14)])
_affine("D",17,F(1,8),[(2,15),(2,16),(10,13),(10,14)])


def affine_weight(kind, order, a, b):
    if not (0<=a<order and 0<=b<order):
        raise ValueError("affine mark outside factor")
    p = _pair(a,b)
    if order <= 9:
        return ONE  # exceptional factors included; classification is external
    g = classical_graph(kind,order,affine=True)
    da,db = _distances(g,a),_distances(g,b)
    distance = da[b]
    if kind == "A":
        if order == 10 and distance == 5: return HALF
        if order == 11 and distance == 5: return HALF
        return F(distance <= 5) if order >= 12 else ONE
    if order <= 17:
        if kind == "B" and order in (16,17) and distance >= 10: return ZERO
        if kind == "C" and order in (15,16,17) and distance >= 10:
            exceptions = {15:{(2,12)},16:{(1,11),(4,14)},17:set()}
            return ONE if p in exceptions[order] else ZERO
        return AFFINE_REPLACEMENTS.get((kind,order),{}).get(p,ONE)
    if distance <= 5:
        return ONE
    features = _features(g,finite=False)
    near = features and min(min(da[f] for f in features),min(db[f] for f in features)) <= 5
    return F(bool(6 <= distance <= 9 and near))


def affine_sum(kind,order):
    return sum((affine_weight(kind,order,a,b) for a,b in combinations(range(order),2)),ZERO)
