# Hostile audit: log-Fano cohomological Coxeter carriers

**Verdict (2026-08-10): the narrowed class theorem is valid with its stated
strong carrier hypotheses.  It is not yet a theorem about every maximal
arithmetic reflection chamber.**

## 1. Local discrepancy check

For a smooth surface pair `(S,B=sum b_iD_i)`, the log discrepancy of a
boundary component is `1-b_i`.  If reduced branches of multiplicity one meet
at a smooth point, the exceptional divisor of the point blowup has log
discrepancy `2-sum b_i`.  Hence `epsilon`-klt gives

```text
b_i < 1-epsilon,
r*delta < 2-epsilon.
```

With the canonical mirror coefficient `1/2`, the first point blowup at an
`n`-fold literal vertex has log discrepancy `2-n/2`; the advertised
`klt => n<=3`, `n=4 => lc`, and `n>4 => not lc` trichotomy follows directly.

For an arbitrary divisorial valuation `E`, the exact inequality is

```text
sum b_i ord_E(D_i) < A(E,S,0)-epsilon.
```

Thus the local cuts in the theorem use log discrepancies, not intersection
heuristics.  They apply only when the P-derived facet sheets occur in the
boundary with the declared positive coefficients.

## 2. Global boundedness check

Alexeev's surface theorem bounds the underlying surfaces admitting an
`epsilon`-log-terminal boundary with `-(K+B)` nef, except for the numerically
Calabi--Yau zero-boundary case.  The theorem's nef-and-big hypothesis removes
that exception.  In a bounded smooth projective family one can choose, after
finite stratification, an integral relative very ample `H` and constants

```text
-K.H <= C_epsilon,       rho(S) <= b2(S) <= R_epsilon.
```

Since `-(K+B).H>=0`, one has `B.H<=-K.H`.  Distinct facet components have
positive integral H-degree and coefficient at least `delta`, giving
`delta*f<=C_epsilon`.  No volume/covolume input occurs.

Jiang's explicit optimal bound on `(-(K+B))^2` is compatible with this
argument but does not replace `H`: a nef anti-log-canonical divisor may have
degree zero on boundary curves.

## 3. Trace-rank check

An `O_k` quadratic module of rank `n+1` has rational rank
`[k:Q](n+1)`.  A faithful injection into `NS(S)/torsion` therefore gives

```text
[k:Q](n+1) <= rho(S) <= R_epsilon.
```

This implication is invalid for an external local system or a decorative
Belyi packet.  The theorem expressly requires an isometric injection of the
entire torsion-free restriction-of-scalars trace form, with its distinguished
`(1,n)` real summand and negative-definite conjugate summands, into
`NS(S)/torsion`.  A list of `n+1` divisor roots is not enough.

For a smooth weak del Pezzo surface, rationality and Noether's formula give
`rho=10-K^2<=9`, proving the explicit endpoint.  This endpoint does not apply
to K3, Coble, or external Gram--Belyi carriers.

## 4. Finite cohomological-image check

Bounded projective surfaces, bounded-degree marked effective divisor cycles,
and bounded-degree integral effective correspondence cycles are
parameterized by a finite union of relative Chow spaces.  Signed cycles are
admissible only when their positive and negative parts are separately marked
and bounded.  After finite flat/smooth stratification, cycle intersections
and induced integral cohomological actions are locally constant.  There are
consequently finitely many action polynomials, coefficient algebras,
projected pairings, and returned Grams, provided one global finite return
grammar is fixed across the class.

The following four uniformity clauses are essential:

- cycle degree bounds effective cycles, or separately bounds both parts of a
  signed cycle;
- one global finite return grammar and finite rational alphabet are used,
  rather than one rule per polytope;
- coefficient algebras and projectors are generated functorially by the
  bounded integral correspondence actions; and
- real-place and orientation data are finite selectors, not storage for a
  moving scalar, period, repeated-eigenspace vector, or cross-ratio.

Those possibilities are excluded in the carrier definition.  Carrier
moduli may remain positive-dimensional; only the cohomological Gram image is
claimed finite.  Fields of definition of the moving surfaces are not claimed
finite.

## 5. Fixed-dimensional boundary-complex obstruction

An `epsilon`-klt pair has no coefficient-one boundary, so its standard dlt
dual complex is empty.  Even if one uses the incidence complex of the full
SNC support instead, its dimension on an `s`-fold is at most `s-1`.
Consequently literal recovery of the boundary complex of an `n`-polytope by
boundary strata forces `n<=s`; a surface cannot intrinsically encode the full
facet-incidence sphere for `n>=3` in this way.

The point-blowup inequality gives the complementary local obstruction.  If
`r>=n` facet branches meet at one smooth surface point, then
`r*delta<2-epsilon`; with coefficient `1/2`, klt already forces `n<=3`.
Higher-dimensional surface compilers must therefore separate sheets and put
the missing face/Gram data in a bounded declared decoration.  This is a
restriction on the compiler, not a global dimension theorem for Coxeter
polytopes.

## 6. Exact scope boundary

The following implications are **not** proved:

```text
reflection maximality => epsilon-log-Fano carrier;
reflection maximality => bounded correspondence degree;
primitive Belyi map => bounded degree;
bounded underlying surface => finite cross-ratio fields;
arithmetic Coxeter polytope => intrinsic log-Fano carrier.
```

Prime maps `x -> x^p` disprove the third implication.  Four-pointed projective
lines with moving algebraic cross-ratio disprove the fourth.  The first,
second, and fifth are the genuine completion gates.  Until they are proved
for an exhaustive P-first class, the result is a carrier-relative theorem
and conditional ABSW replacement, not global arithmetic-reflection
finiteness.

Uniformity is substantive: allowing `epsilon` or `delta` to tend to zero
with `n`, or allowing the correspondence degree `E` to grow with the input
Gram, removes the conclusion.  Such input-dependent parameters may still be
useful encodings, but they are not an ABSW replacement.
