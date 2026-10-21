# Corrected finite-geometry extensions: a real `P -> surface` theorem and a failed nine-row search

**Status (2026-08-09).**  The false extension implication has been replaced by
a rigorous sufficient theorem.  A moving linear system is now accepted only
after the assigned base cluster is resolved and the residual system gives a
relative closed immersion; constant-rank, base-change, bad-locus and final
Hilbert certificates are part of membership.  Thus every accepted branch is
geometrically nonempty, not merely numerically plausible.

For a raw member which independently passes `D_CobVin`, the existing inverse
theorem then returns the **same marked polytope**, and Kondō's classification
pushes genuine restrictions back to the raw input:

\[
 b\in\{1,2,5,8\},\qquad \dim P\in\{10,11,14,17\}.
\]

In characteristic five this sharpens to `b in {1,2}` and
`dim P in {10,11}`.  These bounds are not inserted in the P-side class.

The simultaneous nine-row search found **no nonzero extension**.  The new
`H10/72` Type-VII construction is integrated as a fifth exact core
calibration, but it remains a zero-step core, not an extension.  Consequently
the honest rating is

```text
L3 + corrected sufficient extension theorem;
zero nonzero-extension witnesses in the nine-row census;
not L4 and no TOOT.
```

> **Hostile quantifier warning.**  `C_FG_HilbExt` is a **P-indexed finite
> geometric-certificate class**, not an unqualified P-only class.  The
> coordinate/Hilbert witness is independently checkable and may not import a
> surface label, but its *existence is part of membership*.  Raw
> Coxeter/arithmetic hypotheses do not yet force such a witness and the class
> does not exhaust all possible carriers.  This is the primary reason the
> theorem is not L4; the absence of a nonzero witness in the nine known rows is
> a second, separate defect.

The executable [certificate](../scripts/certify_coble_corrected_extension_pushback.py),
[frozen JSON](../data/enumerations/coble_corrected_extension_pushback.json),
and [tests](../tests/test_coble_corrected_extension_pushback.py) record the
statement and the exact failed-search denominator.

## 1. The corrected P-indexed certificate class

Let `C_FG` be the five arbitrary-input finite-geometry cores now recognized
from raw oriented wall data.  The fifth, Type VII, also has a separate exact
explicit-section realization certificate:

| core | characteristic | boundary count | rank / dimension |
|---|---:|---:|---:|
| Desargues triangle/cut | 5 | 1 | `11 / H10` |
| Type-VII tetrahedral `K6` factorization | 5 | 1 | `11 / H10` |
| Möbius `S(3,4,10)` | 3 | 2 | `12 / H11` |
| polar/chord | 3 | 5 | `15 / H14` |
| grid/projectivity | 3 | 8 | `18 / H17` |

Define `C_FG_HilbExt` as follows.  A normalized oriented raw wall datum `P`
belongs if its exact all-contraction/dependency compiler finds a `C_FG+`
core and there exists a finite certificate for every subsequent center and
curve operation.  The certificate may contain coordinates and polynomials
as witnesses, but membership may not read a surface name, Kondō type, source
archive label, known existence assertion, or stored output Gram.  Every
degree, multiplicity, incidence and lattice class is checked against the
data derived from `P`.

This is an existential, independently verifiable **P-indexed** class: `P` is
the raw input and a finite algebraic certificate is the witness.  It is not
unqualified P-only because the raw arrays do not force witness existence.  It
is a sufficient class, not a claim that every possible surface carrier admits
this particular certificate.

### Center step

Suppose `X_i -> T_i` is the current smooth projective family with
geometrically integral rational-surface fibers and `T_i` is geometrically
integral.  A new center is accepted only on a certified nonempty open of
`X_i` satisfying all required incidence, distinctness and avoidance
conditions.  The universal center is a relative smooth section.  Blowing up
that section gives the next smooth family.

The ideals of the forbidden loci, their saturations, and a Hilbert-dimension
calculation showing that their union is proper are finite certificate data.

### Moving curve step

After resolving the assigned cluster, let `L` be the strict-transform line
bundle and `V` the proposed relative section space.

1. A finite complex representing `R pi_*L`, together with its minors/Fitting
   ideals, certifies that `V` is locally free of constant rank and commutes
   with base change on the same irreducible parameter component.
2. On the relative Hilbert schemes of length one and length two subschemes,
   the universal evaluation maps

   \[
       V\longrightarrow H^0(L|_Z),\qquad \operatorname{length}(Z)=1,2,
   \]

   are surjective.  Chartwise maximal-minor ideals, saturation, and unit-ideal
   elimination certificates make this a finite check.  Equivalently, the
   induced map

   \[
                    X_i\longrightarrow \mathbf P_{T_i}(V^*)
   \]

   is a relative closed immersion.
3. Smooth Bertini and Jouanolou irreducibility then give a nonempty
   generic-fiber open of smooth geometrically integral hyperplane sections.
   Every additional transversality, exact-incidence and avoidance failure is
   represented by a closed bad locus whose Hilbert certificate proves it
   proper.

The length-two clause matters.  Full first-jet generation at each individual
point does not by itself separate two distinct points.

### Explicit curve step

A rigid or rank-one system can instead supply an explicit coefficient vector
with exact evaluation, Jacobian, saturation, absolute-irreducibility and
strict-transform smoothness certificates.  The Type-VII quartics and boundary
sextic use this branch.  This is not disguised Bertini: the actual unique
sections are checked.

### Final strong-Coble Hilbert gate

The final surface is accepted only after exact matrices certify

```text
anticanonical evaluation/jet matrix:       full column rank, h0(-K)=0;
antibicanonical evaluation/jet matrix:     nullity one;
displayed kernel generator:                reduced boundary B=-2K;
boundary strict transforms:                smooth and pairwise disjoint.
```

Thus `|-K|` is empty and `|-2K|={B}` by an explicit Hilbert calculation.  No
informal appeal to “the expected number of conditions” enters the theorem.

## 2. Corrected extension theorem

> **Theorem.** Every finite dependency branch accepted by
> `C_FG_HilbExt` has a geometric realization over an algebraically closed
> field of characteristic different from two.  The result is a marked strong
> Coble root seed whose boundary and root classes agree coefficientwise with
> the raw P-derived classes.  If the same raw input independently belongs to
> `D_CobVin`, its full nef chamber is literally the original marked `P`.

**Proof.**  A core calibration supplies the initial geometric point.  A
center step is a blowup along a relative smooth section and preserves smooth
geometrically integral fibers.  At a curve step, the Fitting/base-change
certificate gives a constant-rank projective bundle of sections.  Relative
closed immersion gives a dense open of smooth geometrically integral members
by Bertini--Jouanolou.  Each further bad locus is proper by its Hilbert
certificate.  Since a projective bundle over a geometrically integral base is
geometrically integral, finitely many such dense opens have nonempty
intersection.  Induction therefore makes the entire simultaneous parameter
tower nonempty.  An explicit-section step is nonempty by its displayed
verified point.  The final two Hilbert ranks give the strong-Coble property.

For `P in D_CobVin`, the previously proved return theorem adds finite volume,
Vinberg finite index, effective-root exhaustion, primitive-cusp regeneration
and hidden-wall exclusion.  It follows with the original marking and
orientation that `P(Nef S)=P`.  This last step is conditional on `D_CobVin`,
not inferred from finite-sidedness.  `square`

The AG ingredients are standard but the hypotheses are intentionally strong.
Cohomology/base change comes from Grothendieck--Dieudonné
[EGA III, parts I](https://www.numdam.org/articles/10.1007/BF02684273/)
and [II](https://numdam.org/articles/10.1007/BF02684890/); the finite
length-two parameter space is supplied by Grothendieck's primary
[Hilbert-scheme construction](https://numdam.org/item/SB_1960-1961__6__249_0/);
and geometric irreducibility is the Bertini theorem of Jouanolou,
[Théorèmes de Bertini et applications](https://books.google.com/books/about/Th%C3%A9or%C3%A8mes_de_Bertini_et_applications.html?id=Q0vvAAAAMAAJ).
For easy theorem-location cross-checks, see the Stacks treatments of
[smooth Bertini](https://stacks.math.columbia.edu/tag/0FD4) and
[Jouanolou irreducibility](https://stacks.math.columbia.edu/tag/0G4C).

### The finite-field distinction

A nonempty geometric Bertini open need not have an `F_p`-rational point.  The
theorem constructs a point over `overline(F_p)`; when all equations descend
to `F_p`, a closed point gives a member over some finite extension
`F_(p^m)`.  An `F_p`-rational curve is claimed only when its coefficient
vector is explicitly certified.  Type VII happens to supply such vectors
over `F5`; the general moving theorem does not promise them.

## 3. Hostile audit of the old clause

The old condition

```text
base-point-free + generically separable
```

is false.  On `P1 x P1`, `|O(2,0)|` is base-point-free and its morphism to its
image is generically separable, but it is composed with the first projection.
A general divisor is the union of two fibers.  The corrected gate rejects it because the map is not a closed
immersion and cannot separate two points on the same projection fiber.

Three further gaps are also closed explicitly:

- constant rank is imposed on one irreducible parameter component, so later
  choices cannot silently jump components;
- every incidence/avoidance open is proved dense on that same component, so
  separately nonempty opens are not illegitimately combined; and
- geometric nonemptiness is distinguished from rationality over the prime
  field.

## 4. Independent AG pushback onto every accepted raw P

Only after construction and exact return do we invoke Kondō's primary
[classification theorem](https://arxiv.org/abs/2102.00591).  For a finite-Aut
strong Coble surface in characteristic not two, its table gives exactly

| characteristic | possible classified rows | boundary counts `b` |
|---|---|---|
| 3 | I, II, V, VI, MI, MII | `1,2,5,8` |
| 5 | I, II, VI, VII | `1,2` |
| other, not 2 | I, II | `1,2` |

If `B_1,...,B_b` are the disjoint `(-4)` boundaries and
`sum B_i=-2K`, then `K^2=-b`.  A rational surface satisfies
`rho+K^2=10`, hence

\[
                         \rho=10+b,\qquad \dim P=9+b.
\]

This gives the following restrictions on the number `m` of clean additional
blowup centers after each core:

| base core | allowed extension lengths `m` |
|---|---|
| Möbius, characteristic 3, `b=2` | `0,3,6` |
| polar/chord, characteristic 3, `b=5` | `0,3` |
| grid/projectivity, characteristic 3, `b=8` | `0` |
| Desargues, characteristic 5, `b=1` | `0,1` |
| Type VII, characteristic 5, `b=1` | `0,1` |

These are actual `surface -> restriction on raw P` consequences inside the
new P-first certificate class.  For example, any purported two-center
Desargues extension can be rejected before a full Vinberg run, while a
Möbius extension can only jump by three or six centers.

## 5. Exact nine-row extension census

The pushback reduces the direct search to a few possible ranks.  The raw
contraction data then give:

| raw row | simultaneous frame denominator | sequential denominator where relevant | result in the corrected class |
|---|---:|---:|---|
| `H10/22` | max disjoint `4`, no full plane frame | `74` P2 targets / `16` signatures | not a core; nonzero extension cannot keep rank 11 |
| `H10/30` | max disjoint `7`, no full plane frame | `573 / 28` | same |
| `H10/78` | `1` full plane frame | — | zero-step Desargues core |
| `H10/72` | `5` full plane frames | — | fifth zero-step Type-VII core |
| `H11/26` | max disjoint `4`, no 11-frame | `52 / 9` | no clean one-center characteristic-five extension |
| `H11/36` | max disjoint `7`, no 11-frame | `890 / 28` | no clean one-center characteristic-five extension |
| `H11/206` | `154` full 11-frames | — | zero-step Möbius core; direct extension scan zero |
| `H14/60` | `5` full plane frames | — | zero-step polar core; `10<30` roots excludes preserved Möbius core |
| `H17/48` | intrinsic 16-center ruled frame | — | zero-step grid core; `24<30` excludes Möbius and `K4,4` has no polar `K5` boundary subdesign |

Thus the denominator is exactly

```text
9 raw rows;
5 zero-step core calibrations;
0 genuinely nonzero extension witnesses.
```

This is a **failed extension search**, not AG pushback, abstract polytope
nonexistence, or nonexistence of every possible infinitely-near surface
model.

### Smallest concrete missing witness

The closest case is the one-rank jump from the characteristic-five Desargues
core to `H11/206`.  Exhaustion of all 154 full 11-exceptional plane frames
gives

```text
144 frames:  10 lines + 10 conics + 10 cubics;
 10 frames:  12 lines + 18 conics;
```

For each of the 144 promising frames, the prospective twenty Desargues root
curves collectively use **all eleven centers**.  Hence there is no common
eleventh center avoided by the old ten-line/ten-conic core.  The count is

```text
common-avoided-center frames:             0 / 144;
clean one-center Desargues witnesses:     0;
Type-VII 12/4/4 line/conic/quartic frames: 0 / 154.
```

This is the exact obstruction to unfreeze next.  A first genuine extension
must either produce a raw row with the permitted rank jump and such a clean
subcore, or enlarge the theorem to controlled infinitely-near modifications
with a comparably rigorous Hilbert certificate.

## 6. Compact versus cusped scope

This Coble theorem is a cusped finite-volume theorem.  After Kondō
classification, every returned chamber has primitive nef isotropic rays,
equivalently genus-one fibrations, so its hyperbolic polytope has ideal
vertices.  It produces no compact member.

Accordingly:

```text
cusped finite-volume scope:  corrected theorem and exact same-P return;
compact scope:               empty inside this Coble class;
global compact conclusion:   none.
```

A compact polytope-to-surface programme needs a different carrier class—for
example an anisotropic Néron--Severi lattice—not a silent reinterpretation of
these Coble results.

## 7. Replay

```bash
python3 scripts/certify_coble_corrected_extension_pushback.py \
  --check-output data/enumerations/coble_corrected_extension_pushback.json
python3 -m pytest -q tests/test_coble_corrected_extension_pushback.py
```
