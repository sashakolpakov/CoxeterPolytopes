"""Finite sparse necessary moments, with exact shared pushforward equations."""
from dataclasses import dataclass, field
from fractions import Fraction as F

from . import conic
from .polynomial import Poly, as_poly, var, monomials


@dataclass
class Domain:
    name: str
    variables: set = field(default_factory=set)
    bounded: set = field(default_factory=set)
    upper_bounds: dict = field(default_factory=dict)
    boolean: set = field(default_factory=set)
    fixed: dict = field(default_factory=dict)
    equations: list = field(default_factory=list)
    inequalities: list = field(default_factory=list)
    matrices: list = field(default_factory=list)
    bases: list = field(default_factory=list)
    integer_counts: dict = field(default_factory=dict)
    rank_total: object = None

    def variable(self, name, bounded=True, boolean=False):
        self.variables.add(name)
        if bounded:
            self.bounded.add(name)
            self.upper_bounds[name] = F(1)
        if boolean:
            self.boolean.add(name)
        return var(name)

    def eq(self, name, expression):
        self.equations.append((name, as_poly(expression)))

    def ge(self, name, expression):
        self.inequalities.append((name, as_poly(expression)))

    def matrix(self, name, entries):
        self.matrices.append((name, [[as_poly(p) for p in row] for row in entries]))

    def basis(self, name, names, degree=1):
        self.bases.append((name, monomials(names, degree)))


class Measure:
    def __init__(self, program, name, domain, mass_scale=1):
        self.program, self.name, self.domain = program, name, domain
        self.mass_scale = F(mass_scale)
        self.moments = {}

    def moment(self, monomial):
        seen, reduced = set(), []
        for v in monomial:
            if v in self.domain.boolean:
                if v in seen:
                    continue
                seen.add(v)
            reduced.append(v)
        m = tuple(sorted(reduced))
        if not set(m) <= self.domain.variables:
            raise ValueError("Unregistered domain variable: "+str(set(m)-self.domain.variables))
        if m not in self.moments:
            self.moments[m] = f"{self.name}_m{len(self.moments)}"
        return self.moments[m]

    def reduce(self, polynomial):
        result = Poly()
        for m, c in as_poly(polynomial).items():
            for v, value in self.domain.fixed.items():
                c *= F(value)**m.count(v)
                m = tuple(w for w in m if w != v)
            if not c:
                continue
            # e=[n=58]. On its support t=1/58. This exact quotient avoids
            # introducing disconnected exceptional-dimension moments.
            if "e" in m and "t" in m:
                c /= 58**m.count("t")
                m = tuple(v for v in m if v != "t")
            remaining, seen = [], set()
            for v in m:
                if v in self.domain.boolean:
                    if v in seen: continue
                    seen.add(v)
                remaining.append(v)
            result += Poly({tuple(remaining): c})
        return result

    def linear(self, polynomial):
        result = {}
        for m, c in self.reduce(polynomial).items():
            c /= self.domain.t_scale**m.count("t")
            c *= self.mass_scale
            name = self.moment(m)
            result[name] = result.get(name, F(0))+c
        return {n: c for n, c in result.items() if c}

    def prepare(self):
        d, p = self.domain, self.program
        p.ge(self.name+"_mass", self.linear(1))
        multipliers = (as_poly(1), var("t"), var("e"))
        for v in sorted(d.variables):
            p.ge(self.name+"_nonnegative_"+v, self.linear(var(v)))
            if v in d.bounded:
                p.ge(self.name+"_upper_"+v, self.linear(d.upper_bounds[v]-var(v)))
            if v not in d.fixed and v not in d.boolean:
                x = var(v)
                p.psd(self.name+"_positive_support_"+v,
                    [[self.linear(x), self.linear(x*x)],
                     [self.linear(x*x), self.linear(x**3)]])
                if v in d.bounded:
                    gap = d.upper_bounds[v]-x
                    p.psd(self.name+"_upper_support_"+v,
                        [[self.linear(gap), self.linear(gap*x)],
                         [self.linear(gap*x), self.linear(gap*x*x)]])
        for name, q in d.equations:
            for k, multiplier in enumerate(multipliers):
                p.eq(self.name+"_"+name+f"_mult{k}", self.linear(q*multiplier))
        for name, q in d.inequalities:
            for k, multiplier in enumerate(multipliers):
                p.ge(self.name+"_"+name+f"_mult{k}", self.linear(q*multiplier))
        for name, entries in d.matrices:
            p.psd(self.name+"_"+name, [[self.linear(q) for q in row] for row in entries])
        for name, basis in d.bases:
            # Restricting to a fixed dimension often turns basis rows into
            # duplicates. Remove exact polynomial dependencies, not PSD data.
            selected, pivots = [], {}
            for original in basis:
                q = self.reduce(original)
                for pivot, row in sorted(pivots.items(), key=lambda item: (len(item[0]), item[0])):
                    q -= q.get(pivot, F(0))*row
                if not q: continue
                pivot = min(q, key=lambda a: (len(a), a))
                pivots[pivot] = q/q[pivot]
                selected.append(self.reduce(original))
            if selected:
                p.psd(self.name+"_moment_"+name, [[self.linear(a*b) for b in selected] for a in selected])

    def close_positive_monomials(self):
        # Every local variable is nonnegative. Only bounded factors are
        # removed in upper inequalities; unbounded geometry is never boxed.
        done = set()
        while set(self.moments)-done:
            m = min(set(self.moments)-done, key=lambda a: (len(a), a)); done.add(m)
            if not m:
                continue
            expression = Poly({m: 1})
            self.program.ge(f"{self.name}_monomial_positive_{len(done)}", self.linear(expression))
            for v in sorted(set(m) & self.domain.bounded):
                remaining = list(m); remaining.remove(v)
                self.program.ge(f"{self.name}_monomial_upper_{len(done)}_{v}",
                    conic.plus(conic.times(self.domain.upper_bounds[v], self.linear(Poly({tuple(remaining): 1}))),
                        conic.times(-1, self.linear(expression))))


class Program:
    def __init__(self):
        self.measures, self.blocks, self.equations, self.labels, self.equality_scales = [], [], [], [], []
        self._block_keys, self._equation_keys = set(), set()

    def measure(self, name, domain, mass_scale=1):
        if any(m.name == name for m in self.measures):
            raise ValueError("Repeated measure")
        result = Measure(self, name, domain, mass_scale); self.measures.append(result)
        return result

    def eq(self, name, expression):
        if expression:
            scale = max(abs(c) for c in expression.values())
            row = conic.times(1/scale, expression)
            key = conic.canonical(row)
            if key not in self._equation_keys:
                self._equation_keys.add(key)
                self.equations.append(row); self.labels.append(name); self.equality_scales.append(scale)

    def ge(self, name, expression):
        self.psd(name, [[expression]])

    def psd(self, name, entries):
        if any(entry for row in entries for entry in row):
            scale = max(abs(c) for row in entries for entry in row for c in entry.values())
            entries = [[conic.times(1/scale, entry) for entry in row] for row in entries]
            key = conic.canonical(entries)
            if key not in self._block_keys:
                self._block_keys.add(key)
                self.blocks.append(dict(name=name, entries=entries))

    def compile(self, metadata):
        for m in self.measures:
            m.prepare()
        for m in self.measures:
            m.close_positive_monomials()
        variables = [n for m in self.measures for n in m.moments.values()]
        metadata = dict(metadata, equality_labels=self.labels,
            equality_positive_divisors=self.equality_scales,
            moment_coordinates={m.name: [[list(k), v] for k, v in m.moments.items()] for m in self.measures},
            measure_mass_scales={m.name: m.mass_scale for m in self.measures},
            dimension_moment_variable="t entries encode the scaled variable parity_endpoint/n",
            all_measures_are_necessary_moments_not_realizability_assertions=True)
        return conic.model(variables, self.blocks, self.equations, metadata)


def dimension_domain(name, threshold):
    d = Domain(name)
    d.t_scale = threshold
    t, e = d.variable("t"), d.variable("e", boolean=True)
    d.upper_bounds["t"] = F(1, threshold)
    d.ge("tail_support", 1-threshold*t)
    d.eq("exceptional_dimension58", e*(58*t-1))
    d.ge("integer_gap_at58", (1-e)*(1-57*t)*(1-59*t))
    d.basis("dimension", ("t", "e"), 1)
    return d
