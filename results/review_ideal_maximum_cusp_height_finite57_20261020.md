# Independent review: maximum cusp height and the ideal five-dimensional obstruction

20 October 2026. Written, adversarial mathematical and primary-source review;
no proof process, census, checker, Lean or build was run.

Reviewed source:
`results/ideal_right_angled_maximum_cusp_height_kissing_candidate_20261020.md`,
SHA256 `95ef1ec3edef6083992aad06dd56e2608b1f120ec098a0d197a2320b55823537`.

**Verdict: PASS as a mathematical implication of the two published theorems
and the standard acute-polytope facts stated below.** I found no sign,
dimension, actual-wall, face-inheritance or equality-propagation gap. The
source has one bibliographic numbering correction: its second external
input is Kolpakov's **Theorem 4**, not Theorem 2. This review does not claim
that a computational certificate or Lean formalization checks the geometry.

## 1. Published inputs checked in their primary sources

Musin's [*The kissing number in four dimensions*](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf),
Annals of Mathematics 168 (2008), 1–32, states its main theorem on printed
page 4. Together with the spherical-code interpretation on printed page 2,
it bounds a set of unit vectors in Euclidean four-space with pairwise
products at most one half by 24. The ambient space here really is
four-dimensional; no three-dimensional kissing bound or uniqueness of the
maximizing configuration is imported.

Kolpakov's [*On the optimality of the ideal right-angled 24-cell*](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf),
Algebraic & Geometric Topology 12 (2012), 1941–1960, **Theorem 4**, printed
page 1945, gives both the 24-facet lower bound for an ideal right-angled
four-polytope and its metric equality characterization. This is the
facet-number theorem, distinct from the minimum-volume theorem immediately
before it. The candidate uses the stronger metric equality conclusion
exactly where it is required.

## 2. The maximizing facet point is actual and correctly normalized

For a fixed future null representative `v` of a cusp, outward unit normals
give `H(m)=-<v,m> >= 0`. Facets through that cusp have height zero. Some
height is positive: in the upper-half-space model with that cusp at
infinity, otherwise all bounding hyperplanes are vertical and the
full-dimensional region has infinite hyperbolic volume near height zero.
Thus the positive maximum `H0` exists because there are finitely many
facets.

For a maximizing normal `n`, put `q=v/H0+n`. Directly,

    q^2=-1,   <q,n>=0,   <q,v>=-H0<0.

The last identity chooses the future component. For every different
actual wall `m`, acute-polytope all-wall acuteness gives `<n,m><=0`, so

    <q,m>=-H(m)/H0+<n,m><=0.

Consequently `q` belongs to the actual facet, not just its supporting
hyperplane. This uses the standard all-wall acuteness theorem for an
acute-angled polytope; assuming it only for adjacent pairs would not
justify this step. Zero-height neighbors and a boundary position of `q`
are permitted at this stage.

## 3. The spherical code has the correct wall set and dimension

Each intrinsic facet of the four-facet `Q` is a ridge with exactly one
actual neighboring facet of the ambient polytope. Its normal `m_i` is
orthogonal to `n`, hence is already a unit intrinsic normal in the
five-dimensional Lorentz space `n^perp`. In a frame centered at `q`,

    m_i=(t_i,sqrt(1+t_i^2)u_i),
    t_i=H(m_i)/H0 in [0,1],   u_i in S^3.

For two different actual walls, all-wall acuteness yields

    u_i.u_j <= t_i*t_j/sqrt((1+t_i^2)(1+t_j^2)) <= 1/2.

The directions cannot coincide: the product of two unit spacelike vectors
with the same spatial direction and nonnegative time coordinates is
strictly positive. Musin therefore applies to the actual intrinsic facet
set, without a visibility, multiplicity or embedded-subdiagram shortcut.

The inherited facet is finite-volume, right-angled and entirely ideal.
Its vertices are vertices of the original face lattice. Right angles
inherit at relative interior points of its ridges; finite volume of the
face is the standard Coxeter-face inheritance fact, or follows here from
the finite ideal vertex hull. Hence Kolpakov applies, and the two bounds
force precisely his metric equality case.

## 4. Metric equality forces the point to be the center

For the regular ideal right-angled 24-cell, the unit outward normal list

    (1,epsilon*e_i+delta*e_j),  i<j, epsilon,delta in {+1,-1}

is correctly normalized: every normal has Lorentz square one. Every
spatial vector occurs with its negative. Transport `q` to `(a,w)`, with
`a>0` and `a^2-|w|^2=1`. The upper bound `t_i<=1` for an opposite pair
gives `a-r.w<=1` and `a+r.w<=1`. Thus `a<=1`, whereas unit future
timelikeness gives `a>=1`. Therefore `a=1`, `w=0`, and every neighboring
height equals `H0`.

This argument needs neither uniqueness of the 24-point spherical code nor
an unproved inference from its equality configuration. It uses only the
published metric classification of the intrinsic four-polytope.

## 5. The global contradiction and scope

The proof repeats at every maximum-height facet. Thus the nonempty set of
such facets is closed under ridge adjacency. The facet adjacency graph of
the finite convex face lattice is connected, including for an ideal
polytope (equivalently, use its Klein-model Euclidean closure). Every
facet would have positive height `H0`, contradicting the actual facets
through the selected cusp.

For dimension greater than five, an actual five-face exists. It inherits
finite volume, right angles and only ideal vertices, so the same
five-dimensional contradiction applies. The conclusion is therefore an
ideal dimension ceiling four, conditional only on the explicitly named
published and standard geometric inputs. It does not concern mixed
polytopes and does not assume a 24-cell decomposition of arbitrary
four-polytopes.

The proof is short but its central mechanism is substantive: a maximum
cusp-height facet simultaneously has an upper facet bound from a spherical
code and a matching lower facet bound with a rigid metric equality case.
I found no remaining local or global mathematical qualification beyond
making the stated all-wall and face-inheritance inputs explicit and
correcting the theorem number above.
