# A four-point constructible-incidence lemma for raw Coble polytopes

**Status (2026-08-09): exact two-member core realization and two label-free
raw certificates; hostile rating `L3`, not `L3+` or flagship `L4`.**

The most naive free-incidence lemma is false on both high-dimensional Coble
positives.  Their residual curves are not produced because independent point
conditions leave a positive-dimensional linear system.  They appear because
square evaluation matrices drop rank on the unique four-point projective
configuration

\[
                         \mathbf P^1(\mathbf F_3).
\]

This note extracts that shared core directly from the raw `H14/60` and
`H17/48` Coxeter data, without consulting surface or Kondō-type labels, and
proves its geometric nonemptiness uniformly.  The two presentations are:

```text
rank 15 / 60 walls:
  four marked tangencies on a conic;
  for every two-subset A, pole(A) lies on chord(complement A);

rank 18 / 48 walls:
  two ordered four-sets on P1;
  every bijection between them is the graph of a projectivity.
```

Both designs exist over an algebraically closed field of characteristic
different from two if and only if the characteristic is three.  In that
characteristic, explicit coordinates construct all boundary curves and all
effective Coble--Mukai root seeds.  The already proved `D_CobVin`
Vinberg/root-exhaustion/Halphen theorem then constructs every integral cusp
lift, excludes all hidden negative curves, and returns literally the input
polytope.

This is a genuine two-member `POLYTOPE -> SURFACE -> same P` calibration.
It is not promoted beyond `L3` because the only two certified members have
no nontrivial extensions, and their
complete four-point cores already encode much of the `H14/H17`
combinatorics and the characteristic.  The class is structural rather than a
union of two archived matrices, but a third input or a genuine extension is
needed before it can plausibly serve as the flagship denominator.

The exact certificate is
[`certify_coble_four_point_constructible_incidence.py`](../scripts/certify_coble_four_point_constructible_incidence.py),
with
[`JSON`](../data/enumerations/coble_four_point_constructible_incidence.json)
and
[`tests`](../tests/test_coble_four_point_constructible_incidence.py).

## 1. The free-incidence proposal fails exactly where the examples live

A natural first proposal was:

> Order the blowup centres and proposed curves acyclically.  At every curve
> step, require the jet-evaluation map to be surjective, the conditions to be
> independent, and the residual system to contain a smooth irreducible
> member by Bertini.

This is a valid sufficient lemma when the kernel still has sections.  It
does not contain either high-dimensional positive row.

### The ruled-grid branch

A `(1,1)` curve on `P1 x P1` has a four-dimensional vector space of
sections.  A proposed graph root must pass through four grid points, one in
each row and column.  Thus its evaluation matrix is `4 x 4`.  For a generic
ordered four-set, exactly four of the 24 permutation conditions have a
nonzero kernel.  At the desired simultaneous solution, all 24 matrices have
rank three, not four:

```text
H0(O(1,1)) dimension       4
point conditions           4
rank at P1(F3)              3
kernel dimension           1
```

Surjectivity would give rank four and zero kernel, hence no curve.

### The conic/tangent branch

A line on `P2` has a three-dimensional vector space of sections.  Each
residual root line is required to pass through three P-derived points: one
intersection of two tangent lines and two complementary points on the
conic.  At the simultaneous solution every `3 x 3` evaluation matrix has
rank two:

```text
H0(O_P2(1)) dimension      3
point conditions           3
rank at P1(F3)              2
kernel dimension           1
```

Again, independent conditions would leave no line.  Therefore the original
free/evaluation-surjective predicate rejects both certified positives.  This
is not a defect in the examples; their special rank drop is precisely the
source of the characteristic-three geometry.

## 2. Label-free dependency extraction from the two raw inputs

No surface label is needed to recover either dependency graph.

### 2.1 The rank-18, 48-wall datum

Start with only the even rank-18 arithmetic lattice and its 48 primitive
mirrors.  The upstream raw compiler exhausts all eight-subsets of the 24
short mirrors:

```text
all eight-subsets                  735,471
pairwise orthogonal frames          21,211
integral Picard de-scalings               1
```

The unique frame reconstructs `I_(1,17)` and numerically partitions the
physical walls as

```text
8 square-minus-four boundaries;
16 square-minus-one curves;
24 square-minus-two curves.
```

The sixteen `(-1)` curves are pairwise disjoint.  Their boundary-incidence
graph, recovered only from intersection number one, is `K4,4`; its
bipartition is intrinsic up to swapping the two sides.  Contracting those
sixteen curves gives two isotropic classes meeting once, hence `P1 x P1`.
The centres are all sixteen transverse intersections of four fibres in each
ruling.

Every `(-2)` wall meets four centres, exactly one in each row and column.
The resulting 24 incidence profiles are all 24 permutations of a four-set,
once each.  Thus the raw wall lattice itself compiles

\[
 C_\sigma\in|\mathcal O(1,1)|,
 \qquad C_\sigma\supset
       \{(a_i,b_{\sigma(i)}):1\le i\le4\},
 \qquad \sigma\in S_4.                                  \tag{2.1}
\]

The centre dependency has depth one and is acyclic.

### 2.2 The rank-15, 60-wall datum

The raw boundary de-scaling reconstructs `I_(1,14)`.  Without reading the
stored curve labels, square and boundary-pairing tests give

```text
 5 boundaries of square -4;
10 actual roots of square -2;
10 crossing curves of square -1 and boundary profile (1,1);
35 cusp-lift walls of square -1 and boundary profile (2).
```

Exact cone enumeration gives five primitive plane rays `H` with

\[
 H^2=1,qquad K.H=-3,qquad
 \{H.B_i\}=\{2,1,1,1,1\}.                              \tag{2.2}
\]

Choose any one; the certificate canonically uses the smallest numerical
conic-boundary index.  The ten crossing curves index every two-subset of the
five boundaries.  Four roots are active on the plane ray, each paired with
the unique crossing curve joining the conic boundary to one of the four
line boundaries.  Together they give four depth-two tangency chains.  The
other six crossing curves are the ordinary intersections of the four lines.

The plane class, the ten crossing classes, and the four total first-tangency
exceptionals form a determinant-one basis with Gram

\[
                         \operatorname{diag}(1,-1^{14}).             \tag{2.3}
\]

For each of the six remaining roots, intersection with this exceptional
basis recovers one ordinary line-line centre indexed by a two-subset `A` of
the four marked tangencies, and the two first-tangency centres indexed by
its complement.  All six two-subsets occur exactly once.  Hence the
dependency grammar, without labels, is

```text
for every A in C({0,1,2,3},2):
  intersection of the tangents indexed by A
  lies on the chord through the two conic points indexed by complement(A).
```

The six ordinary centres have depth one, the four tangency chains have depth
two, and the entire centre dependency graph is acyclic.

## 3. The exact four-point projective-line theorem

Let `k` be algebraically closed of characteristic different from two and
normalize four distinct points on `P1` as

\[
                       (\infty,0,1,\lambda),
                       \qquad \lambda\ne0,1.             \tag{3.1}
\]

### 3.1 Complete graph-projectivity design

Use homogeneous coordinates

\[
 \infty=(1:0),\quad0=(0:1),\quad1=(1:1),\quad
 \lambda=(\lambda:1).
\]

For each permutation `sigma`, evaluate the four basis monomials of
`H0(O(1,1))` at `(a_i,a_sigma(i))`.  The 24 determinant polynomials have the
exact histogram

| determinant | multiplicity |
|---|---:|
| `0` | 4 |
| `-(lambda-1)(lambda+1)` | 4 |
| `2lambda-1` | 4 |
| `lambda(lambda-2)` | 4 |
| `lambda^2-lambda+1` | 4 |
| `-(lambda^2-lambda+1)` | 4 |

If all 24 vanish, the domain `lambda!=1` and two of the equations give

\[
                 \lambda+1=0,\qquad2\lambda-1=0,
\]

and

\[
                 2(\lambda+1)-(2\lambda-1)=3.            \tag{3.2}
\]

Thus `char(k)=3`, and then `lambda=2=-1=1/2`.  Conversely these four points
are `P1(F3)` and

\[
                         PGL_2(\mathbf F_3)\simeq S_4
\]

induces every permutation.  All 24 evaluation matrices have rank three and
give 24 distinct smooth `(1,1)` graph curves.

### 3.2 Complete conic polar/chord design

Use the Veronese conic

\[
                         xz-y^2=0,
\]

with finite point `(t^2:t:1)` and tangent
`x-2ty+t^2z=0`; at infinity the point and tangent are `(1:0:0)` and `z=0`.
For every two-subset `A`, intersect its two tangents and evaluate the chord
through the complementary two points.  The six collinearity polynomials are

```text
(lambda-1)(lambda+1), repeated twice;
lambda(lambda-2),     repeated twice;
1-2lambda,            repeated twice.
```

On the domain (3.1), the first two families force

\[
                   \lambda+1=0,\qquad\lambda-2=0,
\]

and their difference is `3`.  Hence again the design exists only in
characteristic three at `lambda=2`.  Conversely all six polynomials vanish
there.  Each evaluation matrix has rank two and produces one distinct line.

We have therefore proved the common theorem

\[
\boxed{
 \text{either complete four-point design is nonempty}
 \iff \operatorname{char}k=3,
 \quad Q\simeq\mathbf P^1(\mathbf F_3).
}                                                        \tag{3.3}
\]

No Coble surface or classified surface-existence theorem is used in this
calculation.

## 4. Direct construction of the two root-seed surfaces

### 4.1 Grid presentation

Take two copies of `P1(F3)`, four fibres in each ruling of `P1 x P1`, and
blow up their sixteen transverse intersections.  The eight strict transforms
are disjoint `(-4)` curves.  Every element of `PGL2(F3)=S4` gives a distinct
smooth graph, whose strict transform is one of the 24 proposed `(-2)`
curves.  The sixteen exceptional curves are the proposed crossing
`(-1)` curves.  Thus every P-derived Coble--Mukai root seed is effective and
irreducible.

The boundaries satisfy

\[
 \sum_iV_i+\sum_jH_j
   =4F_1+4F_2-2\sum E_{ij}=-2K.                         \tag{4.1}
\]

### 4.2 Conic presentation

Take the Veronese conic and the four tangents indexed by `P1(F3)`.  Blow up
the six distinct line-line intersections and twice along each of the four
conic-line tangencies.  Distinctness is automatic: no three tangent lines to
a smooth conic are concurrent, and a line meets the conic in degree two.
The five strict transforms are disjoint `(-4)` curves and sum to `-2K`.

The ten crossing exceptional curves are smooth `(-1)` curves.  Four strict
transforms of first tangency exceptionals and the six polar/chord lines from
Section 3 give all ten proposed actual `(-2)` curves.  Hence all twenty raw
Coble--Mukai root seeds—ten actual roots and ten roots associated with
crossing exceptionals—are simultaneously effective and irreducible.

### 4.3 The strong-Coble linear systems are automatic

The following argument works in both presentations and for every extension
in Section 6.  Let `B=sum B_i` be the reduced disjoint boundary on the
iterated blowup.  The numerical compiler has already proved

\[
                    B=-2K,\qquad B_i^2=-4.
\]

If `D in |-2K|`, then `D.B_i=-4`, so every `B_i` is a component of `D`.
Their sum exhausts the divisor class; hence `D=B`.  If `A in |-K|`, then
`A.B_i=-2`, so every `B_i` is a component of `A`; but their sum already has
class `-2K`, which makes the residual class `K`, impossible on a rational
surface.  Equivalently, `2A=B` would contradict reducedness.  Therefore

\[
                       |-K|=\varnothing,
                       \qquad|-2K|=\{B\}.                \tag{4.2}
\]

The explicitly constructed surface is strong Coble.  This proves geometric
nonemptiness without importing either known Type-VI or Type-MII surface.

## 5. Why the constructed seed returns the complete same `P`

The construction does not have to draw the 35 rank-15 cusp-lift curves one
by one.  The general `D_CobVin` theorem was designed precisely to turn a
complete effective root seed into the full nef chamber.

Let `Gamma_raw` be the roots compiled from the actual `(-2)` curves and the
crossing `(-1)` curves.  In both inputs:

1. every member of `Gamma_raw` is now effective and irreducible;
2. the P-only Vinberg gate proves
   `[O(CM):W(Gamma_raw)]<infinity`;
3. hence the full effective-root group has finite index and the surface has
   finite automorphism group;
4. Kondō's root-exhaustion theorem makes `Gamma_raw` the entire effective
   irreducible Coble--Mukai root set;
5. primitive isotropic rays are genus-one fibrations, and the Halphen
   theorem realizes `E=f-B_i/2` exactly when the P-derived integral cusp-lift
   test permits it; and
6. every hidden `(-1)` curve would conversely produce one of those same
   primitive cusp lifts.

For the grid input, all `202*8` cusp/boundary tests are nonintegral, so the
48 explicitly drawn curves are already complete.  For the conic input,
exactly 35 of `57*5` tests are integral and Halphen theory produces exactly
those 35 further exceptional curves.  The negative-curve trichotomy leaves
no other wall.  Consequently

\[
                         \boxed{\mathbf P(\operatorname{Nef}S)=P}   \tag{5.1}
\]

with the original integral marking and orientation.  Finite volume is also
independently certified by the literal P-side Vinberg/cone audits; it is not
inferred from a finite wall list.

The crucial separation is:

```text
explicit P1(F3) geometry:       surface and all CM root seeds;
general Coble theorems:         all cusp lifts and no hidden curves;
P-side exact arithmetic:        same marked Coxeter chamber and volume.
```

Known surface existence is nowhere used.

## 6. The finite P-only constructible-incidence predicate

The two cores suggest a precise structural class rather than a union of two
matrices.  Let `D_CobVin^CI4` consist of every `P in D_CobVin` for which at
least one branch of the all-contraction compiler passes the following finite
tests.

1. From only the boundary frame, intersection table, and contraction update,
   derive and canonicalize the centre/curve dependency hypergraph.
2. The non-free dependency core is one complete four-point design of
   Section 3.  This is graph/hypergraph recognition; it does not use rank,
   facet count, a Kondō name, a stored surface, or comparison with a known
   Gram.
3. The complete dependency branch is exhausted by that core: there are zero
   remaining centre or curve operations.
4. The constructed boundary and root-seed intersection matrix agrees
   literally with the P-derived numerical package.

These are finite computations.  The zero-extension clause is essential to
the theorem presently proved.

### Constructible-incidence realization lemma

> Every object of the core-only `D_CobVin^CI4` has a geometric root-seed
> branch over an algebraically closed field of characteristic three.

**Proof.**  Section 3 explicitly realizes either complete core on
`P1(F3)`.  There are no remaining operations.  The exact boundary and root
calculations of Sections 4--5 therefore apply directly.  `square`

The boundary argument (4.2) makes the output strong Coble.  Section 5 then
gives exact same-`P` return for every member, not only for the two archived
calibrations.

The two certified inputs have zero extension steps after their four-point
cores.  They prove a two-member core calibration, not a general extension
theorem.

### Prospective repaired extension verifier

The former extension claim was false as stated.  A base-point-free,
generically separable system can be composed with a pencil; `|O(2,0)|` on
`P1 x P1` is the elementary counterexample, because its general divisor is
reducible.  A future verifier may accept a nonzero curve step only with
either

1. an explicit geometrically integral section certificate, or
2. after resolving the assigned cluster, image dimension at least two,
   separability, and the remaining correct Bertini/avoidance hypotheses.

This repaired predicate is **prospective**.  No nonzero extension is
currently certified, so no nonemptiness theorem for it is claimed here.

## 7. Independent AG pushback

The incidence theorem forces characteristic three but never compares its
output with the surface classification.  Exact same-`P` return makes the
surface finite-Aut, so Kondō's independent classification applies afterward.
It restricts every member of `D_CobVin^CI4` to the characteristic-three
part of the nine-row full-NS list:

| possible surface row in characteristic three | `b` | dimension | facets |
|---|---:|---:|---:|
| I, one boundary | 1 | `H10` | 22 |
| II, one boundary | 1 | `H10` | 30 |
| I, two boundaries | 2 | `H11` | 26 |
| V | 2 | `H11` | 36 |
| MI | 2 | `H11` | 206 |
| VI, five boundaries | 5 | `H14` | 60 |
| MII | 8 | `H17` | 48 |

Thus every object of the structural class is one of at most seven complete
marked full-nef Coxeter types and has dimension in

\[
                         \{10,11,14,17\}.                \tag{7.1}
\]

The two certified raw members occupy the sharp `H14` and `H17` rows.

This pushback is logically independent, but its novelty is qualified.  The
two currently realized four-point cores already determine `b=5` or `b=8`
and characteristic three, so (7.1) gives little new numerical information
on those two members.  It would become materially stronger only after a new
core or a correctly certified nonzero extension produces a dependency
package not already locked to those endpoints.

## 8. Hostile five-gate audit

| gate | verdict | reason |
|---|---|---|
| independent `C` | **pass, narrow** | membership uses only `D_CobVin` and a P-derived complete core with zero extension steps; no surface/type/output Gram |
| finite exhaustive `Car` inside the declared class | **pass for the two certified cores** | each complete four-point core has an explicit point; no extension theorem is invoked |
| exact `Return` | **pass** | constructed root seeds plus Vinberg, root exhaustion, Halphen cusp realization, and hidden-curve exclusion give the identical full nef chamber |
| independent `AG` | **pass** | Kondō is applied only after surface construction and exact return |
| `Pushback` on every a-priori member | **qualified pass** | every member lies in the seven-row characteristic-three list, but current cores already encode most endpoint data |
| broad flagship `L4` | **fail** | only two zero-extension members are certified; the arbitrary-input core recognizer and a valid nonzero extension theorem are absent |

The honest rating is therefore

```text
L3 two-member core calibration;
not L3+ and not the broad class-wide L4 theorem sought by the programme.
```

It is nevertheless directionally stronger than two surface-first lookups:

- both boundary frames and dependency hypergraphs are recovered from `P`;
- one uniform projective-line theorem constructs both surfaces;
- no known surface existence is imported;
- exact same-`P` follows from general Coble theorems; and
- the core recognition does not contain either archived Gram.

## 9. Seven-row result and exact next gate

The requested seven-row extraction is now complete in the
[remaining-seven audit](coble_remaining_seven_contraction_design_audit.md).
All seven rows have rank `11` or `12`, below the minimum ranks `15` and `18`
of the two complete CI4 cores.  Hence no same-core member and no nonzero CI4
extension occurs.

A different design does occur: the raw 206-wall input intrinsically contracts
to the Möbius `3-(10,4,1)` design on `P1(F9)`.  The uniform Frobenius-graph
construction and a `q=3` rigidity formula give a separate structural
`POLYTOPE -> SURFACE -> same P` calibration.  It does not enlarge CI4.

The next finite gate is an arbitrary-input finite-geometry recognizer for the
polar/chord, grid, and Möbius cores.  The naive free predicate must not be
revived, and the repaired extension verifier remains prospective.

## 10. Replay

```bash
python3 scripts/certify_coble_four_point_constructible_incidence.py --check-archive
python3 -m pytest -q tests/test_coble_four_point_constructible_incidence.py
```

The certificate freezes:

- both label-free wall taxonomies and dependency signatures;
- the determinant-one plane exceptional basis;
- all six pair/complement incidences;
- all 24 permutation incidences;
- the complete symbolic determinant alphabets;
- the `24 x rank-3` and `6 x rank-2` characteristic-three evaluations;
- failure of the naive free predicate; and
- the scoped non-overclaim.

The general exact-return input is
[`coble_vinberg_incidence_inverse_transfer.md`](coble_vinberg_incidence_inverse_transfer.md),
and the two upstream P-derived blowdown calculations are in
[`coble_polytope_first_blowdown_incidence_bridge.md`](coble_polytope_first_blowdown_incidence_bridge.md).
