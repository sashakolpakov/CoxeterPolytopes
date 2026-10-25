"""SOURCE ONLY: exact physical-graph to literal-corner adapter.

No repository worker is imported and no file is read at module import. The
caller supplies the independently pinned literal module and one shared fuel
counter. Whole-component geometry and any claimed physical model remain inputs.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')

import ast
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt

MAX_ORDER = 64
MAX_FUEL = 1000000


def require(ok, message):
    if not ok:
        raise ValueError(message)


def charge(fuel):
    fuel[0] += 1
    require(fuel[0] <= MAX_FUEL, 'shared terminal-adapter fuel exhausted')


class Radical:
    """A bounded exact rational linear combination of square-root monomials."""
    def __init__(self, terms):
        self.terms = {d: F(c) for d, c in terms.items() if c}
        require(len(self.terms) <= 16, 'radical term cap')

    @staticmethod
    def scalar(x):
        return Radical({1: F(x)})

    def __add__(self, other):
        out = dict(self.terms)
        for d, c in other.terms.items():
            out[d] = out.get(d, F()) + c
        return Radical(out)

    def __neg__(self):
        return Radical({d: -c for d, c in self.terms.items()})

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        out = {}
        for a, c in self.terms.items():
            for b, e in other.terms.items():
                k = gcd(a, b)
                d = a * b // (k * k)
                require(d <= 10**6, 'radicand cap')
                out[d] = out.get(d, F()) + c * e * k
        return Radical(out)

    def divide(self, other):
        require(set(other.terms) == {1}, 'division is by a nonzero rational only')
        q = other.terms[1]
        return Radical({d: c / q for d, c in self.terms.items()})

    def __eq__(self, other):
        return isinstance(other, Radical) and self.terms == other.terms

    def rational(self):
        require(set(self.terms) <= {1}, 'rational root norm required')
        return self.terms.get(1, F())

    def negative(self):
        for bits in (4, 8, 16, 32, 64):
            upper = F()
            scale = 1 << bits
            for d, c in self.terms.items():
                n = isqrt(d * scale * scale)
                lo = F(n, scale)
                hi = lo if n * n == d * scale * scale else F(n + 1, scale)
                upper += c * (hi if c > 0 else lo)
            if upper < 0:
                return True
        return False


def number(value):
    require(type(value) in (str, int), 'exact scalar schema')
    text = str(value)
    require(len(text) <= 128, 'exact scalar length')
    tree = ast.parse(text, mode='eval')
    require(sum(1 for _ in ast.walk(tree)) <= 40, 'exact scalar syntax cap')

    def visit(node):
        if isinstance(node, ast.Constant):
            require(type(node.value) is int and abs(node.value) <= 10**12, 'bounded integer scalar')
            return Radical.scalar(node.value)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            a = visit(node.operand)
            return -a if isinstance(node.op, ast.USub) else a
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
            if isinstance(node.op, ast.Mult):
                return a * b
            if isinstance(node.op, ast.Div):
                return a.divide(b)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'sqrt':
            require(len(node.args) == 1 and not node.keywords, 'sqrt schema')
            d = visit(node.args[0]).rational()
            require(d.denominator == 1 and 1 <= d <= 30, 'bounded integer square root')
            d = int(d)
            factor = 1
            for k in range(2, isqrt(d) + 1):
                while d % (k * k) == 0:
                    factor *= k
                    d //= k * k
            return Radical({d: F(factor)})
        raise ValueError('unsupported exact scalar syntax')

    return visit(tree.body)


def matrix(encoded):
    require(type(encoded) is list and 1 <= len(encoded) <= MAX_ORDER, 'physical Gram order')
    n = len(encoded)
    require(all(type(row) is list and len(row) == n for row in encoded), 'square physical Gram')
    g = [[number(x) for x in row] for row in encoded]
    for i in range(n):
        require(g[i][i].rational() > 0, 'positive physical root norm')
        for j in range(i):
            require(g[i][j] == g[j][i], 'physical Gram symmetry')
            require(not g[i][j].terms or g[i][j].negative(), 'nonpositive physical wall products')
    return g


def indices(values, n, length=None):
    require(type(values) is list and (length is None or len(values) == length), 'node-list schema')
    require(all(type(i) is int and 0 <= i < n for i in values), 'node-list range/type')
    require(len(values) == len(set(values)), 'duplicate physical node')
    return values


def components(g, ids, fuel):
    remaining, out = set(ids), []
    while remaining:
        charge(fuel)
        seen, queue = {min(remaining)}, [min(remaining)]
        while queue:
            i = queue.pop()
            for j in remaining - seen:
                if g[i][j].terms:
                    seen.add(j)
                    queue.append(j)
        remaining -= seen
        out.append(sorted(seen))
    return out


def labels(g, nodes):
    out = [[0] * len(nodes) for _ in nodes]
    codes = [(Radical.scalar(1), 3), (Radical.scalar(2), 4),
             (Radical({1: F(3, 2), 5: F(1, 2)}), 5),
             (Radical.scalar(3), 6), (Radical.scalar(4), 99)]
    for a, i in enumerate(nodes):
        for b, j in enumerate(nodes[:a]):
            if g[i][j].terms:
                ratio = (g[i][j] * g[i][j] * Radical.scalar(4)).divide(g[i][i] * g[j][j])
                found = [code for value, code in codes if ratio == value]
                require(len(found) == 1, 'unsupported exact native wall label')
                out[a][b] = out[b][a] = found[0]
    return out


def graph(kind, n, weights, affine=False):
    if affine and kind == 'A' and n == 2:
        return [{1: 99}, {0: 99}]
    if affine and kind in ('E', 'F', 'G'):
        out = [dict() for _ in range(n)]
        def edge(i, j, label=3):
            out[i][j] = out[j][i] = label
        if kind == 'E':
            arms = {7: (2, 2, 2), 8: (1, 3, 3), 9: (1, 2, 5)}
            require(n in arms, 'supported exceptional affineE order')
            next_node = 1
            for length in arms[n]:
                previous = 0
                for unused in range(length):
                    edge(previous, next_node)
                    previous, next_node = next_node, next_node + 1
            require(next_node == n, 'exceptional affineE node domain')
        elif kind == 'F':
            require(n == 5, 'affineF4 order')
            for i, label in enumerate((3, 3, 4, 3)):
                edge(i, i + 1, label)
        else:
            require(n == 3, 'affineG2 order')
            edge(0, 1, 3)
            edge(1, 2, 6)
        return out
    if kind in ('A', 'B', 'C', 'D'):
        return weights.classical_graph(kind, n, affine=affine)
    out = [dict() for _ in range(n)]
    def edge(i, j, label=3):
        out[i][j] = out[j][i] = label
    if kind == 'E':
        edge(0, 1)
        edge(0, 2)
        edge(2, 3)
        edge(0, 4)
        for i in range(4, n - 1):
            edge(i, i + 1)
    elif kind in ('F', 'H'):
        for i in range(n - 1):
            edge(i, i + 1, 4 if kind == 'F' and i == 1 else 5 if kind == 'H' and i == 0 else 3)
    else:
        raise ValueError('unsupported target graph')
    return out


def embeddings(left, right, marks, fuel):
    n = len(left)
    require(n == len(right), 'full isomorphism orders differ')
    degree = [sum(bool(x) for x in row) for row in left]
    order = sorted(range(n), key=lambda i: (-degree[i], i))
    used, assignment, result = set(), {}, set()
    def visit(depth):
        charge(fuel)
        if depth == n:
            result.add(tuple(sorted(assignment[i] for i in marks)))
            return
        i = order[depth]
        for j in range(n):
            if j in used or degree[i] != len(right[j]):
                continue
            if any(left[i][a] != right[j].get(b, 0) for a, b in assignment.items()):
                continue
            used.add(j)
            assignment[i] = j
            visit(depth + 1)
            del assignment[i]
            used.remove(j)
    visit(0)
    return sorted(result)


def classify(g, nodes, marks, weights, fuel):
    lab, n = labels(g, nodes), len(nodes)
    if n == 2 and lab[0][1] == 99:
        require(len(marks) == 2, 'parallel affineA1 must contain both black marks')
        return ('A', 2), [(0, 1)], True
    if n == 2 and lab[0][1] in (5, 6):
        return ('I', 2), [(0,), (1,)] if len(marks) == 1 else [(0, 1)], False
    candidates = [('A', n)]
    if n >= 2:
        candidates.append(('B', n))
    if n >= 4:
        candidates.append(('D', n))
    if n in (6, 7, 8):
        candidates.append(('E', n))
    if n == 4:
        candidates.append(('F', 4))
    if n in (3, 4):
        candidates.append(('H', n))
    position = [nodes.index(m) for m in marks]
    found = []
    for typ in candidates:
        im = embeddings(lab, graph(*typ, weights), position, fuel)
        if im:
            found.append((typ, im, False))
    # B2 and its reversed orientation are one type, not two candidates.
    require(len(found) <= 1, 'ambiguous finite graph classification')
    if found:
        return found[0]
    require(len(marks) == 2, 'an affine factor cannot lie wholly in the white support')
    require(n >= 3, 'unsupported small affine factor')
    affine_candidates = [(kind, n) for kind in ('A', 'B', 'C', 'D')
                         if not (kind == 'B' and n < 4 or kind == 'C' and n < 3
                                 or kind == 'D' and n < 5)]
    if n in (7, 8, 9):
        affine_candidates.append(('E', n))
    if n == 5:
        affine_candidates.append(('F', 5))
    if n == 3:
        affine_candidates.append(('G', 3))
    for kind, order in affine_candidates:
        im = embeddings(lab, graph(kind, order, weights, affine=True), position, fuel)
        if im:
            found.append(((kind, order), im, True))
    require(len(found) == 1, 'unsupported or ambiguous full finite/affine component')
    return found[0]


def corner(encoded, white, marks, weights, fuel, extra=(), complete=False,
           sealed_witnesses=()):
    g = matrix(encoded)
    indices(white, len(g))
    indices(marks, len(g), 2)
    require(not set(white) & set(marks), 'black marks occur in white support')
    ids = sorted(set(white) | set(marks))
    parts = components(g, ids, fuel)
    selected = [next(i for i, p in enumerate(parts) if m in p) for m in marks]
    shared = selected[0] == selected[1]
    witness_sets=[tuple(sorted(p))for p in sealed_witnesses]
    require(len(witness_sets)==len(set(witness_sets)),'duplicate sealed-witness component')
    require(all(list(p)in parts and not set(p)&set(marks)for p in witness_sets),
            'a sealed witness must be a distinct unmarked displayed component')
    forest, images, affine, retained, ignored = [], {}, [], {}, []
    for i, part in enumerate(parts):
        local = [m for m in marks if m in part]
        if not local and tuple(part)not in witness_sets:
            ignored.append(part)
            continue
        typ, im, is_affine = classify(g, part, local, weights, fuel)
        retained[i]=len(forest)
        forest.append(typ)
        images[retained[i]] = im
        if is_affine:
            affine.append(retained[i])
    selected=[retained[i]for i in selected]
    if affine:
        require(shared and affine == [selected[0]], 'affine corner marks do not share the unique affine block')
        typ = forest[selected[0]]
        values = {weights.affine_weight(*typ, *pair) for pair in images[selected[0]]}
        require(len(values) == 1 and not extra, 'affine image invariance or unexpected finite spectator augmentation')
        return dict(weight=str(next(iter(values))), kind='ideal', forest=[list(t) for t in forest],
                    marked_images=[list(x) for x in images[selected[0]]], marks=marks,
                    guard_possible=False, complete=False,ignored_unsealed_components=ignored)
    full = forest + [tuple(t) for t in extra]
    rank = sum(t[1] for t in full)
    require(rank <= 58, 'known finite forest exceeds ambient58')
    require(not complete or(rank==58 and not ignored),
            'complete augmentation needs every displayed spectator sealed and exact rank58')
    choices = images[selected[0]] if shared else product(images[selected[0]], images[selected[1]])
    values, possible = set(), False
    for image in choices:
        first = (selected[0], image[0] if shared else image[0][0])
        second = (selected[1], image[1] if shared else image[1][0])
        value = weights.finite_corner_weight(full, first, second, final=complete)
        if not complete and shared and full[selected[0]] == ('D', 13) and image in ((0, 11), (0, 12)):
            # Unsealed displayed spectators may grow. Do not subtract their
            # ranks or use an omitted-component capacity argument here.
            possible=True
            value=min(value,F(1,2))
        values.add(value)
    require(len(values) == 1, 'literal payment not invariant under all native component maps')
    return dict(weight=str(next(iter(values))), kind='finite', forest=[list(t) for t in full],
                marks=marks, marked_component_indices=selected, guard_possible=possible,
                complete=complete,sealed_witness_components=[list(p)for p in witness_sets],
                ignored_unsealed_components=ignored)
