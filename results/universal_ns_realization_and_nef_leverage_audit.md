# Universal Neron--Severi realization versus genuine nef-cone leverage

**Status (2026-08-09): theorem-level audit.**

The lattice-realization, cone, and reservoir claims have been independently
checked in the
[`hostile audit`](polytope_surface_inverse_theorems_hostile_audit.md), with a
qualified-pass verdict after explicit repairs.

This note separates two statements which must not be conflated:

```text
the polytope lattice occurs as NS(X)
```

and

```text
the given Coxeter chamber is Nef(X), or a complete certified nef section.
```

The first statement is automatic in a surprisingly broad range.  The second
is where all useful algebraic geometry enters.

## 1. Automatic exact-Picard K3 compiler through trace rank eleven

Let `P` be an arithmetic Coxeter datum with a nondegenerate integral
Lorentzian trace lattice `(L,B)` of rank `r`.  When an integral-reflection
conclusion is invoked, assume in addition that its specified reflections
preserve the underlying `Z`-module `L`.  This includes a rational
integral-reflection datum with `r=dim(P)+1`, and an order-stable totally-real
admissible datum of field degree `d`, where a specified order `A` in the
field acts on `L`, with

```text
r = d(dim(P)+1).
```

Multiplication of the form by two produces the even lattice `L(2)`.  It does
not alter any hyperplane, angle, oriented halfspace, or reflection:

```text
s_v(x) = x - 2 B(x,v)/B(v,v) v
       = x - 2 (2B)(x,v)/(2B)(v,v) v.
```

Here the underlying abelian group is unchanged.  Consequently common
doubling also leaves the **status** of reflection integrality unchanged: it
preserves integrality when it was present, but cannot repair a nonintegral
reflection.

Morrison's Corollary 2.9, together with Remark 2.11 for existence in rank
eleven, gives:

> **Automatic weak K3 theorem.**  If `r<=11`, there is a smooth projective
> complex K3 surface `X` with
> `NS(X) isometric to L(2)` exactly.  In ranks at most ten the primitive
> embedding in the K3 lattice is unique up to its full orthogonal group; in
> rank eleven existence remains automatic, although uniqueness is not part
> of the statement.

The marked lattice-polarized locus has period dimension `20-r`; outside its
countable Noether--Lefschetz hyperplanes the Picard lattice is exactly
`L(2)`, and its positive direction makes the K3 projective.  For a
totally-real trace construction, every element of the chosen order `A` that
stabilizes `L` acts as an integral self-adjoint endomorphism of `NS(X)`.
This may be promoted to every algebraic integer only when `L` is actually
`O_K`-stable.  Its distinguished real eigenspace and the transported divisor
classes recover `P` as a marked semialgebraic divisor eigenslice.  A rational
divisor-product correspondence realizes that endomorphism on cohomology; no
surface automorphism is asserted or needed.

Consequently:

- every rational integral-reflection Coxeter datum through `H10` has an
  automatic exact-Picard weak K3 carrier after common doubling;
- every admissible order-stable trace-lattice datum of the preceding type
  with `d(n+1)<=11` has the same automatic trace-K3 carrier; in particular,
  the quadratic-field range in this class is automatic through `H4`;
- mere nonemptiness of an NS-realization moduli problem cannot prune a single
  member of this entire class.

For rank `12..20`, exact-Picard K3 existence for the **even lattice
`N=L(2)`** is equivalent to a primitive embedding in

```text
Lambda_K3 = U^3 + E8(-1)^2,
```

or, equivalently, to existence of an even complement of signature
`(2,20-r)` with discriminant form `-q_N`.  Nikulin's local discriminant-form
criterion decides this exactly.  Failure here excludes this K3 carrier
species, not the Coxeter polytope.  Rank above twenty is impossible for a
complex K3, but is not a bound for general surfaces.

Primary references:

- D. Morrison, *On K3 surfaces with large Picard number*, Invent. Math. 75
  (1984), 105--121, [DOI](https://doi.org/10.1007/BF01403093), especially
  Corollary 2.9 and Remark 2.11.
- V. Nikulin, *Integral symmetric bilinear forms and some of their
  applications*, Math. USSR-Izv. 14 (1980), 103--167,
  [DOI](https://doi.org/10.1070/IM1980v014n01ABEH001060), especially the
  primitive-embedding and discriminant-form criteria.
- I. Dolgachev, *Mirror symmetry for lattice polarized K3 surfaces*,
  [arXiv](https://arxiv.org/abs/alg-geom/9502005), for the marked
  lattice-polarized period construction.

## 2. Exact nef realization is a strictly stronger, exact root condition

For a K3 surface with `NS(X)=N`, the closed nef cone is the closure of one
component of the positive cone cut by **all** hyperplanes `delta^perp` with
`delta in N` and `delta^2=-2`.  Riemann--Roch chooses an effective sign, and
the K3 Weyl group acts transitively on the chambers.  Let `C` be a
full-dimensional chamber cone in one component of the positive cone, with
irredundant facets.  Then, after choosing a marking and the corresponding
effective signs, a polytope-first K3 package returns `C` as the whole nef
cone if and only if

1. `N` embeds primitively in `Lambda_K3`;
2. every irredundant facet of `C` is a primitive square-minus-two mirror; and
3. no other square-minus-two mirror meets the interior of `C`.

Thus uniform doubling, while it guarantees exact Picard realization in low
rank, usually destroys the square-minus-two wall condition and deliberately
produces only a weak carrier.  This is the precise boundary between a soft
NS compiler and a compiler with cone-theoretic leverage.

Combined with the direct full-nef compactness theorem already frozen in
[`direct_full_nef_surface_compactness_obstruction.md`](direct_full_nef_surface_compactness_obstruction.md),
this gives a real return theorem: any surface-free class universally promoted
to direct full-nef integral-reflection surface carriers has compact dimension
at most three.  It does not apply to irrational totally-real eigensections.

## 3. The NS lattice and the displayed classical invariants do not determine the nef cone

There is an elementary marked counterexample inside rational surfaces.  Let

```text
X_gen = Bl_{p1,p2,p3}(P2),  with the three points noncollinear,
X_col = Bl_{q1,q2,q3}(P2),  with the three points collinear.
```

Both have the same marked data

```text
NS = Z H + Z E1 + Z E2 + Z E3 = I_(1,3),
K  = -3H + E1 + E2 + E3,
K^2=6,  c2=6,  chi(O)=1,  q=pg=0.
```

On `X_gen`, `-K` is ample: this is the degree-six del Pezzo surface.  On
`X_col`, the strict transform

```text
R = H-E1-E2-E3
```

is an effective curve with `R^2=-2` and `(-K).R=0`; hence `-K` is not ample
and `R` supplies an extra wall.  Therefore even `(NS,K)` together with all
the displayed classical numerical invariants does not determine the ample,
nef, or effective cone.

This proves a formal obstruction to the desired strategy: an inverse
compiler which records only an NS lattice, canonical class, and the displayed
Chern/Hodge numbers cannot have a necessary-and-sufficient exact-nef
conclusion.  It
must also reconstruct period/incidence data and exhaust the effective
negative curves (or use a special category, such as K3 surfaces, where a
root theorem performs that exhaustion).

## 4. Why the obvious larger reservoirs do not repair the gap

### Blowups

Blowing up gives unbounded Picard rank, but each blowup adds an orthogonal
square-minus-one exceptional class, and the effective cone depends on the
incidence of the centers, as the preceding example shows.  This is why the
universal rational-surface lattice embedding is flexible enough to carry
every arithmetic datum weakly and too flexible to return a polytope theorem.

### Products of curves

For `C x D`, the two fiber classes are nef, isotropic, and intersect once.
More precisely, the Picard group has the standard `Z^2` fiber part plus a
correspondence term controlled by `Hom(J_C,J_D)`; the latter is not an
arbitrary negative lattice.  In particular, every full nef cone of a product
contains isotropic boundary rays, so products cannot be compact direct
full-nef carriers.  They can only be considered for proper invariant or
eigen-sections, where cone completeness again needs a separate theorem.

### Finite and branched covers

If `f:X->Y` is finite surjective of degree `e`, then

```text
(f^*D1).(f^*D2) = e(D1.D2)
```

and, by the projection formula,

```text
Nef(X) intersect f^*NS(Y)_R = f^*Nef(Y).
```

The reverse implication follows by lifting a curve on `Y` to a curve on
`X`.  Thus a finite cover transports a cone already controlled downstairs;
it does not promote an arbitrary weak sublattice of `NS(Y)` into a new exact
nef cone.  Noether--Lefschetz results can sometimes prove that a very general
branched cover has no unexpected divisor classes, but then the cover merely
inherits the scaled downstairs NS/cone.  If new classes do occur, their
effective walls are an additional completeness problem.

## 5. Strong conclusion for the inverse program

The strongest currently justified general statement is therefore:

```text
P -> exact NS surface
```

is universal in low trace rank and, for the admissible integral lattice-form
class covered by the rational-blowup construction, universal as a weak
sublattice carrier in all ranks.  It therefore cannot by itself produce a
dimension bound or nonexistence theorem.  The smallest carrier notion with
genuine leverage is

```text
P -> (exact NS + exhaustive effective-root/MBM data + exact nef or
      invariant/eigen-nef equality).
```

K3 surfaces provide an exact implementation through complete square-minus-two
root chambers.  Rational/Looijenga surfaces require period characters and an
exhaustive exceptional-curve test.  Products and branched covers do not give
a missing universal shortcut.  A K3 embedding failure above rank eleven is
only a carrier-relative obstruction unless an independently stated
surface-free polytope class is proved to admit no alternative carrier.
