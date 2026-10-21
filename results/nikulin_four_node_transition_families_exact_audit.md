# Exact audit of Nikulin's five marked four-node transition families

**Date:** 2026-08-10  
**Verdict:** **PASS after making the word “type” precise.**  There is no
missing or extra four-node family.  There are, however, **seven** marked
support topologies.  Nikulin's five drawings are five *schematic parameter
families*: `Gamma3` contains the extra support stratum `s=2`, and `Gamma4`
contains the extra support stratum `k=2`.  Since Coxeter label `2` means Gram
entry zero, each of those drawn edges is absent in the special stratum.

This is an independent finite reconstruction from the Gram conditions.  It
does not import a Coxeter classification table, a Lanner list, a field-degree
bound, or a dimension theorem.  It certifies the topology and parameter
grammar used in Section 3.3 of
[arXiv:0708.3991v2](https://arxiv.org/abs/0708.3991) and retained in Section
3.3 of the [final transition paper](https://doi.org/10.1070/IM2011v075n05ABEH002561).
It does **not** certify the later numerical interval extrema.

## 1. Exact support theorem

Number the marked endpoint-exclusive walls `1,2` and the two auxiliary walls
`3,4`.  The marked entry is

```text
a_12 = u > 2,
```

and the five possible finite Coxeter positions are

```text
13, 14, 23, 24, 34.
```

We quotient by the full marking-preserving relabelling group

```text
S({1,2}) x S({3,4}).
```

The endpoint conditions say that deletion of node `1` and deletion of node
`2` are negative definite.  At the support level, either deletion therefore
cannot contain a triangle.  Indeed, every present finite edge has entry at
least `1`; on a three-cycle the vector `(1,1,1)` has quadratic value

```text
-6 + 2(a+b+c) >= 0,
```

so that deletion is not negative definite.  Conversely, every allowed
three-node forest is made negative definite by assigning label `3` to all its
edges.  Thus the support test is exact, not merely necessary.

Exhausting all `2^5=32` masks gives:

```text
32 raw masks
17 connected masks with both endpoint deletions admissible
 7 orbits under S2 x S2.
```

The seven orbits and their placement in the five drawings are:

| Canonical finite support | Printed family | Parameter stratum |
|---|---|---|
| `13,14` | `Gamma4` | `k=2`; drawn edge `23` is absent |
| `13,24` | `Gamma5` | generic |
| `13,34` | `Gamma3` | `s=2`; drawn edge `14` is absent |
| `13,14,23` | `Gamma4` | `k>=3` |
| `13,23,34` | `Gamma2` | generic |
| `13,24,34` | `Gamma3` | `s>=3` |
| `13,14,23,24` | `Gamma1` | generic |

So “five types” is correct only in the schematic sense used by the parameter
lists.  Read as “five literal unweighted support graphs,” it is false: the
literal number is seven.

There is no hidden inertia condition left over.  A negative-definite
three-node endpoint deletion gives at least three negative directions in the
full Gram matrix.  The principal block

```text
[ -2  u ]
[  u -2 ]
```

has one positive direction for `u>2`.  Interlacing therefore forces full
signature `(1,3)` and nonsingularity.

As a non-emptiness check, the certificate sets `u=3` and every present finite
edge to label `3`.  Every one of the seven orbits then has endpoint deletions
negative definite and full inertia `(1,3)`.  Their determinants, in the order
of the table, are

```text
-28, -27, -19, -43, -35, -28, -60.
```

All entries and cyclic products are rational integers, so these are exact
`Q`-valued V-arithmetic witnesses; the nongeometric embedding condition is
vacuous.  Consequently V-arithmeticity does not make any of the seven support
strata empty.

## 2. The parameter lists are complete

For a three-node path with Coxeter labels `p,q`, the Gram matrix is negative
definite exactly when

```text
cos(pi/p)^2 + cos(pi/q)^2 < 1.
```

If both labels are at least `4`, the left side is at least `1`.  Hence one
label is `3`.  With the other label `q`, equality occurs at `q=6`, and
monotonicity gives exactly

```text
(p,q) = (3,3), (3,4), (3,5), up to order.
```

Applying this lemma to both endpoint deletions gives precisely the printed
domains:

| Family | Complete parameter grammar up to the displayed symmetries |
|---|---|
| `Gamma1` | `s=k=3`, `5>=r>=p>=3`; or `s=p=3`, `5>=r>=k>=4` |
| `Gamma2` | `p=3`, `3<=s<=k<=5`; or `s=k=3`, `p=4,5` |
| `Gamma3` | `s=2,k=3,r=3..5`; `s=2,k=4,5,r=3`; `r=3,3<=s<=k<=5`; or `s=k=3,r=4,5` |
| `Gamma4` | arbitrary `k>=2`; `(s,r)=(3,3),(3,4),(3,5),(4,3),(5,3)` |
| `Gamma5` | arbitrary `3<=s<=k` |

The script independently matches every admissible assignment over the exact
logical alphabet

```text
{2,3,4,5,6,101}.
```

Here `2` is absence, `3,4,5` are every possible label in a two-edge
spherical path, and `6,101` are two ordered representatives of the free
`>=6` regime.  Labels `>=6` cannot occur in a two-edge endpoint-deletion
path, so the only unbounded positions are the arbitrary `k` in `Gamma4` and
the two arbitrary leaf labels in `Gamma5`.  The analytic path lemma therefore
upgrades this finite categorical replay to all integer labels.

The replay finds:

```text
7,776 categorical raw assignments
  249 admissible raw assignments
   73 marking-preserving label orbits
    0 uncovered assignments
    0 cross-family ambiguities.
```

The `73` orbits split as

```text
Gamma1  9
Gamma2  8
Gamma3 13
Gamma4 28
Gamma5 15.
```

The latter two counts use the two representatives of the unbounded label
regime and are grammar checks, not finite enumerations of the infinitely many
numeric labels.

There is one harmless redundancy worth recording.  In `Gamma4` with `k=2`,
the auxiliary nodes become interchangeable.  Thus the printed rows

```text
(s,r)=(3,4),(4,3)     and     (3,5),(5,3)
```

represent one marked graph in each pair.  Nikulin later treats both orders in
the same numerical case, so this is duplication, not a missing case.

The exact symbolic determinant of each of the five `4 x 4` Gram matrices is
also regenerated from its support and frozen in the archive.  This catches a
minor source typo in arXiv:0708.3991v2: inside the `Gamma2` determinant
discussion, the sentence introducing the V-arithmetic inequality calls the
graph `Gamma4`.  The displayed matrix, determinant, and subsection are all
`Gamma2`, so the subscript is only a typographical error.

## 3. Ground-field inheritance, without importing Nikulin

Here is the exact general lemma needed by the transition argument.

**Lemma.** Let `A` be a V-arithmetic Gram matrix with cyclic-product ground
field `K`.  Let `A_I` be a principal submatrix which is not negative definite,
and let `K_I` be the field generated by all cyclic products inside `I`.  Then

```text
K_I = K.
```

If `A_I` is hyperbolic, it is itself V-arithmetic over this same field.  The
field equality also holds when `A_I` is parabolic, although the word
V-arithmetic is conventionally reserved for a hyperbolic matrix.

**Proof.** Suppose `K/K_I` were proper.  Since the fields are totally real and
separable, there is a nonidentity real `K_I`-embedding of `K`, and it extends
to the totally real coefficient field.  It fixes every two-cycle
`a_ij^2`, hence changes every nonzero packet entry only by a sign.  It also
fixes every longer cyclic product, so the product of these edge signs around
every cycle is `+1`.  Choose a root in each support component and propagate
vertex signs along a spanning tree.  The cycle condition says that every edge
sign is `d_i d_j`.  Therefore the conjugate packet is

```text
D A_I D,    D=diag(d_i),    d_i in {+1,-1},
```

and has the same inertia as `A_I`.  It is not negative definite.  But the
V-arithmetic condition makes the full Gram matrix at that nongeometric
embedding negative definite, and every principal submatrix of a
negative-definite matrix is negative definite.  Contradiction.

The certificate exhausts this sign-switch statement on all `38` connected
four-vertex supports and all `624` edge-sign assignments.  Exactly `304` are
cycle-consistent, and all `304` are vertex switches.  The proof by a spanning
forest is valid in every rank.

This separates the three marked regimes cleanly:

```text
u>2:  the marked 2-node packet is hyperbolic, so K=Q(u^2);
u=2:  the marked 2-node packet is parabolic with field Q, so K=Q;
u<2:  the marked pair is negative definite, so another non-negative-
      definite packet is required.
```

The last line is precisely where the printed Lanner-only reduction failed.
A minimal indefinite packet can be quasi-Lanner.  Compactness is irrelevant
to inheritance: the full quasi-Lanner packet is hyperbolic and inherits `K`.
Moreover, in the exact `[3,3,6]` counterexample, the affine `[3,6]` principal
deletion is parabolic and has cyclic products `1,3`, so the strengthened lemma
itself gives `K=Q`.

For a general finite-volume noncompact arithmetic chamber, the same `K=Q`
conclusion follows without a quasi-Lanner classification: a cusp makes the
defining `K`-quadratic form isotropic, whereas every nongeometric real
conjugate is definite.  A nonzero isotropic vector remains nonzero after an
embedding, so a nontrivial conjugate is impossible.

The essential carrier warning is equally exact.  Field inheritance requires
the **full cyclic-product algebra** of the principal packet.  An unlabeled
support graph, an affine null vector, or a few Coxeter labels can lose field
data into an undeclared scalar, period, or cross-ratio and does not satisfy the
lemma.

## 4. What is now certified, and what is not

Certified here:

- completeness of the five schematic four-node families;
- the seven literal support strata and both label-`2` degenerations;
- the complete endpoint-deletion parameter grammar;
- exact symbolic Gram determinants and rational V-arithmetic witnesses;
- ground-field inheritance for hyperbolic, quasi-Lanner, and parabolic
  principal packets;
- the `u=2 => K=Q` branch and the field repair of the `[3,3,6]` cusp branch.

Not certified here:

- the later interval extrema or all numerical degree estimates;
- finite-volume realization of every V-arithmetic graph;
- any global dimension theorem;
- any implication from a finite field list to bounded dimension.

For the project-native AG route, the correct four-node “Coxeter dessin” has
five templates but seven support strata, with edge absence part of the data.
It is now a complete local grammar.  It still needs to sit inside one globally
faithful bounded cohomological packet before it can contribute to a dimension
bound.

## 5. Frozen replay

- `scripts/certify_nikulin_four_node_transition_families.py`
- `data/enumerations/nikulin_four_node_transition_families.json`
- `tests/test_nikulin_four_node_transition_families.py`

Pinned SHA-256 values:

```text
7219cf1de954c108a550332b30efb895cab93a49c20ce7c3024ee35ac2a2d3b0  script
f32af7c7b6ecfc00bcc6a2ac4fc0169468f46e715750a93f27370243cf161bb8  archive
59a63f266a8bb95dc85240cb79a972da28a51fac687bb2affc867c0732431f2a  tests
```

The run is deliberately small: its largest finite loop has `7,776` rows and
its only symbolic calculation consists of five `4 x 4` determinants.
