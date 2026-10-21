# A maximal reflection group with no Borcherds correction in `H16`

**Status (2026-08-10): proved and exactly certified.**  This is a new
project-level theorem/cross-identification, not a new abstract polytope.

## The result

Let

\[
 S=U\oplus E_8(-1)\oplus E_7(-1),\qquad
 \Gamma=W_2(S).
\]

The complete square-minus-two chamber of `S` is Tumarkin's finite-volume
noncompact Coxeter polytope `Gamma16` in `H16`.  Then:

1. `S` has no proper integral overlattice;
2. every integral reflection of `S` is a square-minus-two reflection, so
   `W_all(S)=W_2(S)`;
3. `Gamma` is maximal among **all discrete reflection groups** in
   `Isom(H16)`, hence in particular among arithmetic reflection groups; and
4. `S` has **no automorphic correction** in Wang's sense.

Consequently

\[
 \boxed{\text{maximal arithmetic reflection group}
        \not\Longrightarrow
        \text{Borcherds automorphic correction}.}
\]

There is an essential warning: `Gamma` is **not** a maximal arithmetic
lattice/group.  Its unique diagram involution gives the proper arithmetic
nonreflection extension

\[
 \Theta_{16}=\Gamma_{16}\rtimes C_2.
\]

Thus four notions that had been drifting together are now separated:

| notion | answer for `U+E8(-1)+E7(-1)` |
|---|---|
| maximal integral lattice in its rational space | yes |
| norm-two group equals the group of all integral reflections | yes |
| maximal discrete/arithmetic **reflection** group | yes |
| maximal arithmetic group | **no**, because of `Theta16` |
| Wang/Borcherds automorphic correction | **no** |

## 1. Exact identification with Tumarkin's `Gamma16`

The archived complete K3 curve diagram has 19 square-minus-two walls, Gram
rank 17 and inertia `(1,16,2)`.  The full Vinberg test gives

```text
finite vertices       125
ideal vertices          3
elliptic edges        1087
facets                  19
bad edge extensions      0
bad cusp neighborhoods   0
```

It is therefore a finite-volume noncompact Coxeter polytope in `H16` with
`n+3=19` facets.  [Tumarkin's theorem](https://arxiv.org/abs/math/0311272)
states that this polytope exists and is unique, which identifies the K3
chamber with `Gamma16` without a diagram-name guess.

[Kellerhals, Example 6](https://doi.org/10.3390/math5030043) records its
complete `f`-vector

```text
(128, 1087, 4768, 14000, 30352, 50960, 67960, 72908, 63204,
 44200, 24752, 10948, 3740, 952, 170, 19, 1),
```

three vertices at infinity, and

\[
 \chi(\Gamma_{16})=
 \frac{2499347}{2360171042879569920000}.
\]

The project computation independently matches the identifying entries:
`128=125+3` vertices, `1087` edges and `19` facets.  The alternating sum of
the published vector is exactly one.

## 2. Lattice maximality and all integral reflections

The first 17 simple roots form a principal Gram minor of determinant two.
Their exact coordinates reconstruct all 19 wall roots integrally.  Since
`|disc S|=2`, they are an `S`-basis, not merely a rational basis.

If `S` admitted an integral overlattice `M` of index `h>1`, then

\[
 |\det M|=\frac{2}{h^2},
\]

which cannot be an integer.  Thus `S` is maximal integral, including against
odd as well as even overclosures.

Now let `r` be primitive and suppose its hyperplane reflection preserves
`S`.  Put `d=div_S(r)`.  Reflection integrality gives

\[
 r^2\mid 2d,
\]

and the exact order of `r/d` in `A_S` is `d`.  Here

\[
 A_S\cong\mathbf Z/2,qquad q_S(1)=\tfrac12\pmod{2\mathbf Z},
\]

so `d` is one or two.  The complete candidate table is

| `div(r)` | `r^2` | `q(r/div(r))` | result |
|---:|---:|---:|---|
| 1 | -2 | trivial class | allowed |
| 2 | -2 | `3/2` | not `1/2` |
| 2 | -4 | `1` | not `1/2` |

Hence every integral reflection has norm `-2`, and
`W_all(S)=W_2(S)`.

## 3. Prasad-free maximality among all reflection groups

This is the cleanest proof and does not use ABSW, Prasad, arithmetic volume,
or surface boundedness.

Suppose `Gamma16<G`, where `G` is any discrete hyperbolic reflection group.
Because `Gamma16` is a lattice, it has finite index in every discrete
overgroup.  The main theorem of
[Felikson--Tumarkin](https://arxiv.org/abs/math/0305093) says that a
finite-index reflection subgroup has at least as many fundamental-chamber
facets as the ambient reflection group.  Therefore a chamber `Q` of `G` has
at most 19 facets.  Any full-dimensional `H16` polytope has at least 17, so
only three cases remain.

- `17` facets would make `Q` a finite-volume Coxeter simplex, impossible in
  `H16`; the complete hyperbolic Coxeter-simplex list stops in `H9`.
- `18=n+2` facets are impossible in `H16` by the complete
  [finite-volume `n+2` classification](https://arxiv.org/abs/math/0301133).
- With `19=n+3` facets, `Q` is Tumarkin's unique noncompact `H16` polytope,
  hence isometric to the original chamber.

The overgroup remains nonuniform because it contains the parabolic elements
of `Gamma16`.  Thus the noncompact classifications apply.  In the last case
the chamber volumes are equal, so `[G:Gamma16]=1`.  This proves maximality
among all discrete reflection groups, a stronger statement than arithmetic
reflection maximality.

## 4. Independent volume and normalizer cross-check

The earlier arithmetic-volume route remains useful as a fully independent
check, but it is no longer a proof dependency.  Kellerhals' Proposition 1
and Table 5, based on
[Belolipetsky's minimal-orientable theorem](https://www.numdam.org/item/ASNSP_2004_5_3_4_749_0.pdf),
identify the minimal orientable cusped arithmetic `H16` orbifold with the
orientation subgroup of `Theta16`.  Its Euler characteristic is exactly
`chi(Gamma16)`.

If an arithmetic reflection overgroup has index `m`, its orientation
subgroup has Euler characteristic

\[
 \frac{2\chi(\Gamma_{16})}{m}.
\]

Minimality therefore gives `m<=2`.  An index-two reflection extension would
make the `Gamma16` chamber the double of a smaller chamber and produce a
hyperplane-reflection symmetry of its wall diagram.

The exact weighted diagram automorphism group is `C2`.  Its unique
involution is

```text
(1 18)(2 15)(3 19)(4 16)(5 17)(6 14)(7 13)(8 12)(9 11)(10).
```

On the free 19-node space it has `(+,−)=(10,9)`.  On the two-dimensional
Gram radical it has `(1,1)`, so on the 17-dimensional Lorentz quotient it
has

\[
                       (+,-)=(9,8).
\]

A hyperplane reflection would have `(16,1)`.  Thus the involution is not a
reflection.  The 17 basis walls also show directly that every wall-permuting
symmetry preserves `S`: its images are integral wall roots and its action on
the basis is unimodular.  This is precisely the nonreflective symmetry that
forms `Theta16`, not a missed reflection extension.

## 5. Wang's obstruction and the completed 18-row experiment

[Wang's Theorem 7.2](https://doi.org/10.4171/JEMS/1358) places
`U+E8(-1)+E7(-1)` among the 18 hyperbolic 2-reflective lattices with no
automorphic correction.  The companion
[exact 18-row ledger](wang18_same_lattice_mixed_root_ledger.md) constructs a
strict same-lattice longer-root reflection overgroup for every one of the
other 17 negative rows:

- `A2` supplies a primitive `(-6,div 3)` reflection;
- `D_n`, `A3`, or two `A1` summands supply a primitive `(-4,div 2)`
  reflection.

Therefore the formerly proposed maximality test is completely decided:

\[
 \boxed{17\text{ nonmaximal negative rows}
        +1\text{ maximal negative row}
        +0\text{ undecided}.}
\]

Equivalently, inside Wang's classified 69-row denominator,

\[
 W_2(S)\text{ reflection-maximal}
 \Longrightarrow
 \bigl(S\text{ has a correction}\bigr)
 \text{ or }
 S\cong U+E_8(-1)+E_7(-1).
\]

This implication is only one-way; no claim is made that every corrected row
is reflection-maximal.

## 6. What maximality means on the algebraic-geometric side

This example gives a concrete negative answer to the proposed
"maximality implies AG finiteness" shortcut.

For a generic exact-Picard carrier `X`,

\[
 NS(X)=U+E_8(-1)+E_7(-1),
\]

so the K3 period dimension is `20-17=3`.  Roulleau's model describes the
coarse family as an open part of `WP3(2,3,5,6)`.  The geometric construction
is especially close to the desired Coxeter-dessin picture.  In the appendix
to Galluzzi--Lombardo,
[Dolgachev](https://sites.lsa.umich.edu/idolga-new/wp-content/uploads/sites/1467/2024/08/appendixGL04.pdf)
starts from a decorated genus-two curve `(W,p)` satisfying

\[
 p\text{ non-Weierstrass},\qquad |3K_W-p|\ne\varnothing.
\]

It produces a two-cuspidal `(3,3)` curve in `P1 x P1`; the double cover
branched over that curve and two rulings resolves to the K3 surface.  For a
generic genus-two curve there are 16 relevant decorated pairs, giving a
generically `16:1` marked cover of genus-two moduli.  The 19 rational curves
then recover the complete `Gamma16` nef chamber.

This is a genuine finite decoration of a three-dimensional moduli problem,
not a zero-dimensional dessin classification.  The marked numerical Gram
image is one point while the carrier varieties and their fields vary in a
three-dimensional family.  Hence:

> reflection maximality is a no-extension property, not a bounded Hilbert
> polynomial, bounded divisor degree, rigidity, or finite-field invariant.

It cannot by itself recover ABSW-style finiteness from surfaces.  A
polytope-to-surface finiteness theorem still needs independent bounded
numerical data plus a zero-dimensional **Gram image** or exact-return locus;
maximality supplies neither.  The direct facet-rank proof above is useful
Coxeter theory, but it is not an AG boundedness theorem.

## 7. Claim boundary and reproduction

The polytope, K3 family, Wang obstruction and source classifications are
published.  New project content is their exact cross-identification, the
three-level maximality separation, the certified counterexample to
`maximal reflection => correction`, and the complete `17+1+0` experiment.
There is no new abstract hyperbolic polytope in this result.

Reproduce the small exact certificate and focused tests with

```bash
python3 -m scripts.certify_wang_h16_maximal_reflection_no_correction --check-archive
python3 -m pytest -q tests/test_wang_h16_maximal_reflection_no_correction.py
```

Frozen artifacts:

- [`wang_h16_maximal_reflection_no_correction.json`](../data/enumerations/wang_h16_maximal_reflection_no_correction.json);
- [`certify_wang_h16_maximal_reflection_no_correction.py`](../scripts/certify_wang_h16_maximal_reflection_no_correction.py); and
- [`test_wang_h16_maximal_reflection_no_correction.py`](../tests/test_wang_h16_maximal_reflection_no_correction.py).
