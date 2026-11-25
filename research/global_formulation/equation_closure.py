"""Sparse polynomial-ideal closure of pointwise domain identities.

Aggregate equations between different occupation measures are not pointwise
identities and must never be multiplied by an arbitrary state coordinate.
"""
from . import conic
from .polynomial import Poly, as_poly, var

MAX_IDENTITY_DEGREE = 4
MAX_CLIQUE_VARIABLES = 8
MAX_IDENTITY_TERMS = 24


def variance_cliques(measure):
    domain, program = measure.domain, measure.program
    selected, skipped = [], []
    for label, original in domain.equations:
        q = measure.reduce(original)
        if not q: continue
        degree = max(map(len,q))
        if degree > MAX_IDENTITY_DEGREE or len(q.support) > MAX_CLIQUE_VARIABLES or len(q) > MAX_IDENTITY_TERMS:
            skipped.append(label); continue
        basis = [as_poly(1)]+[var(v) for v in sorted(q.support)]
        basis += [Poly({m: 1}) for m in sorted(q,key=lambda m:(len(m),m)) if m]
        unique = {}
        for b in basis:
            reduced = measure.reduce(b)
            unique.setdefault(conic.canonical(reduced),reduced)
        basis = list(unique.values())
        # An invertible constant basis change exposes q as a zero diagonal
        # after L(q^2)=0. Removing its leading monomial preserves the span.
        pivot = max(q,key=lambda m:(len(m),m))
        basis = [b for b in basis if b != Poly({pivot:1})]+[q]
        program.psd(measure.name+"_ideal_clique_"+label,
            [[measure.linear(a*b) for b in basis] for a in basis])
        program.eq(measure.name+"_ideal_variance_"+label,measure.linear(q*q))
        for i,b in enumerate(basis):
            program.eq(measure.name+"_ideal_clique_kernel_"+label+f"_{i}",measure.linear(q*b))
        selected.append(dict(identity=label,degree=degree,basis_size=len(basis)))
    measure.ideal_closure_report = dict(variance_cliques=selected,
        identities_outside_variance_clique_budget=skipped,
        maximum_identity_degree=MAX_IDENTITY_DEGREE,
        maximum_clique_variables=MAX_CLIQUE_VARIABLES,
        maximum_identity_terms=MAX_IDENTITY_TERMS,
        all_cliques_use_the_same_measure_moment_registry=True)


def existing_products(measure):
    program, domain = measure.program, measure.domain
    available = set(measure.moments)
    added, unknown = 0, 0
    for label, original in domain.equations:
        q = measure.reduce(original)
        if not q: continue
        for v in sorted(domain.variables-domain.fixed.keys()):
            product = measure.reduce(q*var(v))
            if not product: continue
            if all(m in available for m in product):
                before = len(program.equations)
                program.eq(measure.name+"_ideal_times_"+label+"_"+v, measure.linear(product))
                added += len(program.equations)-before
            else: unknown += 1
    measure.ideal_closure_report.update(existing_coordinate_products_added=added,
        products_requiring_unretained_moments=unknown,
        aggregate_incidence_equations_were_not_multiplied=True)


def self_checks():
    from fractions import Fraction as F
    from .moments import Program, dimension_domain
    def construct(enabled):
        d = dimension_domain("closure_fixture",57)
        d.fixed.update(t=F(1,57),e=F(0))
        x,z = d.variable("x",bounded=False),d.variable("z",bounded=False)
        d.eq("product",x*z-1)
        d.basis("joint",("x","z"),1)
        d.equation_closure = enabled
        p = Program(); m = p.measure("fixture",d)
        return p.compile({}),m
    def point(m,values):
        result = {}
        for monomial,name in m.moments.items():
            total = F(0)
            for z in values:
                value = F(1)
                for variable in monomial:
                    value *= F(1) if variable == "x" else z
                total += value
            result[name] = total/len(values)
        return result
    before,m = construct(False)
    conic.check_point(before,point(m,[F(1,2),F(3,2)]))
    after,m = construct(True)
    conic.check_point(after,point(m,[F(1)]))
    mixed = point(m,[F(1,2),F(3,2)])
    if not any(conic.value(q,mixed) for q in after["equations"]):
        raise ValueError("The variance closure failed to reject a cancelling identity mixture")
    return dict(positive_two_point_mixture_passes_mean_identity=True,
        variance_closure_rejects_that_mixture=True,
        genuine_point_passes_all_new_cliques_and_equations=True)
