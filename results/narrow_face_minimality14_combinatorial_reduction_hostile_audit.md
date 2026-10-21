# Hostile audit: the narrow-face minimality-14 combinatorial reduction

**Verdict (2026-08-10): a self-contained compact conditional lemma survives;
an unconditional finite-volume dimension bound does not.**  The visual-cone
calculation gives the exact normalized threshold `14`.  For a compact Coxeter
chamber, simplicity and the `h`-vector incidence calculation are valid and,
under one explicit component-size hypothesis, force `n<=9`.  A cusped chamber
need not be simple at an ideal vertex, and the component-size hypothesis has
not been independently derived from minimality `14` or reflection maximality.

## 1. The local threshold

Use the square-minus-two convention for outward facet normals:

```text
e.e=f.f=-2.
```

Fix an interior point `O` and choose a facet `E` whose visual cone has minimal
angular aperture.  Let two other closed facets meet `E`, and write
`theta_1,theta_2` for their apertures and `theta_12` for their mutual visual
separation.  The cone-plane calculation is

```text
e.f = 4 X - 2,
X = sin((theta_1+theta_12)/2) sin((theta_2+theta_12)/2)
    / (sin(theta_1/2) sin(theta_2/2)).
```

Minimality gives

```text
0 <= theta_12 <= theta_E <= theta_1,theta_2 < pi.
```

For `0<=b<=a<pi`,

```text
sin((a+b)/2)/sin(a/2) < 2.
```

Indeed, if `a<=pi/3`, the numerator is at most `sin(a)` and the ratio is
`2cos(a/2)<2`; if `a>pi/3`, the numerator is at most one and
`1/sin(a/2)<2`.  Each factor in `X` is therefore below two, so

```text
                         e.f < 14.                       (1)
```

This proof uses the diagonal normalization essentially.  The constant is the
sharp supremum of this cone estimate: setting all three angles equal to
`epsilon` gives `e.f=16cos^2(epsilon/2)-2`, which tends to `14`.  Thus this
argument cannot be tightened to a smaller universal number merely by making
the inequalities non-strict.

The calculation is local and does not use arithmetic classification.  When
facets meet only in the projective closure, the closed-facet incidence and
the three cone angles must remain part of the certificate; an affine or
cuspidal adjacency tag cannot silently replace them.

## 2. Exact compact face average

Let `P` be a **compact** Coxeter chamber in `H^n` and let `Q=E` be its narrow
facet, of dimension `m=n-1`.  Every vertex link of `P` is a spherical Coxeter
chamber and hence a simplex.  Thus `P` and `Q` are simple polytopes.

Let `Q*` be the simplicial dual and let `(h_0,...,h_m)` be its `h`-vector.
Shellability gives `h_i>=0`, Dehn--Sommerville gives `h_i=h_(m-i)`, and
`h_0=h_m=1`.  Duality and the standard `f`--`h` relation give

```text
f_0(Q) = f_(m-1)(Q*) = sum_i h_i,

f_2(Q) = f_(m-3)(Q*)
       = (1/2) sum_i [C(m-i,2)+C(i,2)] h_i.             (2)
```

If `a` is the average number of vertices of a two-face of `Q`, simplicity
also gives

```text
a f_2(Q) = C(m,2) f_0(Q).                              (3)
```

The bracketed weight in (2) is minimized at the middle index or indices.
Because `h_0=h_m=1` occur at strictly larger weights, the resulting bound is
strict:

```text
a < 4 + 4/(m-1) = 4 + 4/(n-2),   m odd  (n even),
a < 4 + 4/(m-2) = 4 + 4/(n-3),   m even (n odd).       (4)
```

This is the whole upper-bound calculation; no polytope census is hidden in
it.

## 3. The precise conditional dimension lemma

Assume additionally that, at every vertex of `Q`, every connected component
of the **full finite vertex Coxeter graph of `P`, including the node `E`,**
has at most two nodes.  The graph is then a matching plus isolated vertices.
After deleting `E`, at most `floor(m/2)` edges remain.  Orthogonal projection
of the other wall normals to `E` creates no extra nonright pair: since `E`
has degree at most one, two different remaining normals cannot both pair
nontrivially with its normal.  Hence at most `floor(m/2)` of the `C(m,2)`
two-face corners at that vertex are nonright.

A compact hyperbolic triangle has at least two nonright angles, and a compact
hyperbolic quadrilateral has at least one: otherwise its angle sum would not
be strictly below `(k-2)pi`.  A `k`-gon therefore has at least `5-k`
nonright corners for `k=3,4`, while for `k>=5` the inequality remains true
with a nonpositive right side.  Summing first over vertices and then over
two-faces yields

```text
(5-a) f_2(Q) <= floor(m/2) f_0(Q).
```

Using (3), this is exactly

```text
a [C(m,2)+floor(m/2)] >= 5 C(m,2).                     (5)
```

For even `n`, the lower bound from (5) is `5(n-1)/n`; at `n=10` it equals
`4.5`, while (4) is strict `<4.5`.  The lower bound increases and the upper
bound decreases thereafter.  For odd `n`, (5) gives `5(n-2)/(n-1)`; at
`n=11` it again equals `4.5`, contradicting the strict upper bound, with the
same monotonicity afterwards.  Therefore

> **Compact conditional lemma.**  If every full finite vertex graph of `P`
> at a vertex of the narrow facet has connected components of size at most
> two, then `n<=9`.

The component-size premise is indispensable.  It does not follow from (1):
a spherical Coxeter component may have three or more nodes while every entry
is far below `14`.  Deriving the premise requires a separate, exhaustive
small-graph/field exclusion.  The alleged five marked four-node families and
their label tuples have not yet been generated in this project from the
endpoint-deletion and inertia predicates alone.  Until that replay exists,
(5) is a checker-ready conditional lemma, not a global dimension theorem.

## 4. Why the cusped export fails

A finite-volume noncompact Coxeter chamber need not be simple at an ideal
vertex.  The regular ideal right-angled octahedron in `H^3` is already a
counterexample: four walls meet at each cusp, and the horospherical link is a
Euclidean square rather than a simplex.  More generally an affine Coxeter
link may be a product of Euclidean simplices.  Consequently the projective
closure of a narrow facet need not be simple, so neither (2)--(3) nor the
matching count at an ideal vertex is available.

Truncating cusps does not repair the proof automatically: it introduces new
facets and changes both the two-face ledger and the nonright-angle counts.
A valid noncompact replacement must certify the full affine link, control its
number of components and incident walls, and derive a new horospherical
incidence inequality.  The separate statement that a noncompact arithmetic
reflection chamber has ground field `Q` is only a field restriction; it does
not restore simplicity or bound `n`.

## 5. Exact scope

The independently retained implications are

```text
minimal visual aperture + normalized cone formula
    => a local Gram window e.f<14;

compactness + simple-polytope h-vector identities
    => the strict two-face average bound (4);

those facts + the explicit component-size-two premise
    => n<=9.
```

What is **not** certified is

```text
minimality 14 or reflection maximality
    => component size at most two;

finite-volume with ideal vertices
    => the compact h-vector/incidence identities;

a finite local field/Gram packet
    => a global dimension bound.
```

The broader dependency record is in
[`nikulin_conceptual_dependency_hostile_audit.md`](nikulin_conceptual_dependency_hostile_audit.md),
and the global local-to-rank obstruction is frozen in
[`integral_packet_global_dimension_bridge_hostile_audit.md`](integral_packet_global_dimension_bridge_hostile_audit.md).
