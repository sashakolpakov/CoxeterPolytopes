# An exact strong-Coble carrier iff on Alexeev's nine genus-zero chambers

**Status (2026-08-09): exact nine-row carrier iff and seven exclusions; `L3`,
not strict `L4` under the strengthened five-part contract, and not an
abstract-polytope nonexistence theorem.**

There is now a genuine non-singleton `POLYTOPE -> SURFACE` theorem.  Its
denominator is fixed independently of the positive answers: Alexeev's nine
complete mixed-reflection chambers on the 2-elementary genus-zero line

\[
                              r+a=22.
\]

With the canonical strong-Coble wall decoration, exactly two of the nine
admit an exact full-nef carrier over an algebraically closed field of
characteristic different from two:

\[
                  \boxed{(15,7,1)\quad\hbox{and}\quad(18,4,0).}
\]

Both occur in characteristic three.  They are respectively Kondō's
five-boundary Type-VI surface and eight-boundary Type-MII surface, and their
full projectivized Néron--Severi nef cones are exactly the marked `H14/60`
and `H17/48` Alexeev chambers.

## 1. Exact statement and decoration

Let `A0` be the following nine raw chambers, specified by their complete
primitive wall Gram matrices:

```text
(11,11,1), (12,10,1), (13,9,1), (14,8,0), (14,8,1),
(15,7,1),  (16,6,1),  (18,4,0), (20,2,1).
```

The class is not defined by asking whether a surface exists.  It is the
independently published `r+a=22` subcollection of Alexeev's even
2-elementary reflection chambers.

We use the decoration forced by a strong-Coble quotient interpretation

\[
                 M=\widetilde{\operatorname{Pic}}(S)(2).
\]

Every primitive square-minus-four, divisibility-two facet is required to be
the image of an actual downstairs irreducible `(-2)`-curve.  Such a curve is
disjoint from the branch boundary and its tame K3 pullback splits into two
exchanged `(-2)`-curves.  Square-minus-two facets are **not** all subjected
to this test: they may encode a boundary half or a downstairs `(-1)`-curve.

For `P in A0`, let `C(P)` mean that there are

1. an algebraically closed field `k` with `char(k) != 2`;
2. a strong Coble surface `S/k`; and
3. a marking respecting the preceding norm/divisibility decoration,

such that the complete projectivized full-`Pic(S)` nef cone, including its
oriented facet set, is exactly `P`.  Then

\[
 C(P_{r,a,\delta})
 \quad\Longleftrightarrow\quad
 (r,a,\delta)\in\{(15,7,1),(18,4,0)\}.                 \tag{1.1}
\]

This is an existence iff over **some** allowed characteristic.  In both
positive cases the characteristic is necessarily three.

## 2. Why the surface classification applies

Suppose `C(P)` holds, and let `b` be the number of boundary components of
`S`.  Since the components are disjoint `(-4)`-curves and sum to `-2K_S`,

\[
                         K_S^2=-b.
\]

Noether's formula for a rational surface gives

\[
                  \rho(S)=10+b.
\]

The chamber `P_(r,a,delta)` has hyperbolic dimension `r-1`, while an exact
full-Picard nef chamber has dimension `rho(S)-1`.  Hence

\[
                             b=r-10.                    \tag{2.1}
\]

The exact finite-sided full-nef hypothesis also forces `Aut(S)` to be
finite.  Indeed, its image on the Picard lattice permutes the finite facet
set.  The facet normals span because the chamber is full-dimensional and
finite-volume, so the image is finite; the kernel of the natural Coble
action on the Coble--Mukai lattice is finite.  Kondō's classification of
strong Coble surfaces with finite automorphism group in arbitrary
characteristic different from two therefore applies.

Its nine rows have boundary counts

```text
1, 2, 1, 2, 1, 5, 1, 2, 8.
```

Consequently an exact carrier in this denominator must have

\[
                           b\in\{1,2,5,8\}.              \tag{2.2}
\]

Combining (2.1) and (2.2) immediately removes the ranks
`13,14,16,20`.

## 3. The polytope-first split-root cut removes ranks 11 and 12

The only survivors of the boundary-count cut are

```text
(11,11,1), (12,10,1), (15,7,1), (18,4,0).
```

For every long wall `v_i`, write its forced K3 splitting as

\[
 v_i=c_i^++c_i^-,\qquad d_i=c_i^+-c_i^-.
\]

If `n_ij=(v_i.v_j)/2`, then, after choosing sheets,

\[
 H_{ii}=2,\qquad
 H_{ij}\in\{-n_{ij},-n_{ij}+2,\ldots,n_{ij}\},
 \qquad H=-\tfrac12\operatorname{Gram}(d_i)\succeq0.    \tag{3.1}
\]

This is a finite condition on the raw wall Gram.  The exact split-root
census proves that neither `(11,11,1)` nor `(12,10,1)` has **any**
positive-semidefinite completion (3.1).  Their unit-intersection graphs are
forests, so every sheet signing is switching-congruent to one fixed
representative; those representatives have inertias `(9,1,0)` and
`(10,1,0)`.  The obstruction works in every characteristic different from
two and precedes periods, effectivity, or a Picard-rank bound.

This step is essential.  Kondō's list does contain other one- and
two-boundary Coble surfaces, but their complete full-nef wall systems are
different chambers (`H10/22`, `H10/30`, and `H11/26` among the certified
examples).  Equality of the invariant lattice does not imply equality of
the nef polytope.

## 4. The two survivors return exactly

The remaining two rows are not merely numerically feasible.

| chamber | `b=r-10` | polytope-first incidence output | exact surface return |
|---|---:|---|---|
| `(15,7,1)` | 5 | conic plus four tangent lines; residual-line equations force `char=3` | Kondō Type VI; all 60 negative curves, all `57*5` cusp lifts, and the full cone give exactly `H14/60` |
| `(18,4,0)` | 8 | `4 x 4` grid on `P1 x P1`; all 24 permutation graphs exist iff `char=3` | Kondō Type MII; all 48 negative curves, all `202*8` cusp lifts, and the marked return give exactly `H17/48` |

Thus existence closes the reverse implication of (1.1).  The split lattice
also independently explains the positive-characteristic jump: Type VI
needs an anti-invariant difference lattice `E6(-2)`, while Type MII needs
`(A2+A2)(-2)` and forces a rank-22 determinant-nine K3 lattice, hence
characteristic three and Artin invariant one.

The complete decision table is therefore

| raw chamber | `b` | decisive carrier test | result |
|---|---:|---|---|
| `(11,11,1)` | 1 | no positive-semidefinite split-root completion | no carrier |
| `(12,10,1)` | 2 | no positive-semidefinite split-root completion | no carrier |
| `(13,9,1)` | 3 | boundary count absent from Kondō's exhaustive list | no carrier |
| `(14,8,0)` | 4 | boundary count absent from Kondō's exhaustive list | no carrier |
| `(14,8,1)` | 4 | boundary count absent from Kondō's exhaustive list | no carrier |
| `(15,7,1)` | 5 | exact Type-VI return | carrier in `char=3` |
| `(16,6,1)` | 6 | boundary count absent from Kondō's exhaustive list | no carrier |
| `(18,4,0)` | 8 | exact Type-MII return | carrier in `char=3` |
| `(20,2,1)` | 10 | boundary count absent from Kondō's exhaustive list | no carrier |

For `(13,9,1)`, both rank-14 rows, and `(16,6,1)`, the split-root census
also gives independent obstructions; the last even has the stronger
rank-22 determinant-twelve K3-lattice contradiction.  These redundancies
are not needed for the short iff proof.

## 5. Why this is `L3`, not strict `L4`

This earns a **finite-class `L3` carrier criterion and classification**:

```text
independently defined nine raw decorated chambers
  -> necessary P-only split test + exhaustive AG classification
  -> seven exclusions
  -> two independently certified exact same-P returns.
```

The class was chosen before, and without reference to, the two accepted
answers.  Both directions of the carrier predicate are closed for every
member, so the seven negative rows and two positive rows are genuine
theorems.

That is nevertheless not strict `L4` under the strengthened five-part
contract.  The carrier compiler/decision procedure is exhaustive, and the
two realized branches return the same input polytopes, but **seven members of
the a priori class have no realized carrier branch**.  Thus class-wide
essential surjectivity fails rather than merely remaining unproved.  The AG
classification pushes back a complete carrier/no-carrier partition, not a
surface construction and same-`P` return for every `P in A0`.

It does **not** say that the other seven arithmetic Coxeter polytopes fail
to exist; Alexeev's chambers are the inputs.  It says that they cannot be
the exact full nef cones of strong Coble surfaces under this canonical
quotient decoration.  It also does not establish a universal
polytope-to-surface theorem beyond the nine-row denominator, nor replace the
still-open uniform no-hidden-wall theorem for the general blowdown/incidence
compiler.

## 6. Frozen evidence

The nine-row polytope-only split computation is frozen in
[`alexeev_g0_split_lattice_census.json`](../data/enumerations/alexeev_g0_split_lattice_census.json),
generated and tested by the files with the same stem in `scripts/` and
`tests/`.  The two exact returns are frozen in
[`coble_type_vi_char3_five_boundary_full_ns_h14.json`](../data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json)
and
[`coble_type_mii_full_ns_h17.json`](../data/enumerations/coble_type_mii_full_ns_h17.json).
The template-free blowdown calculations are in
[`coble_polytope_first_blowdown_incidence_bridge.md`](coble_polytope_first_blowdown_incidence_bridge.md).

Primary references:

- Shigeyuki Kondō,
  [*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
  especially Theorem 1.1 and Table 2.
- Valery Alexeev,
  [*Reflective hyperbolic 2-elementary lattices, K3 surfaces and hyperkähler varieties*](https://arxiv.org/abs/2209.09110),
  with the complete wall matrices in the v4 ancillary file.
