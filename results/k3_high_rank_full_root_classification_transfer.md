# Classification transfer for high-rank complete K3 root chambers

**Frozen scope date:** 2026-08-08.  
**Classified range:** Picard rank `11..20`, equivalently hyperbolic dimension
`10..19`.  
**Positive output:** the finite-volume locus has exactly 19 types, in
dimensions `10..18`.  
**Negative output:** the finite-volume locus has no rank-20 / `H19` type.  
**Novelty boundary:** the lattices and polytopes are published examples.  The
new project result is the surface-free finiteness-criterion formulation, its
structural carrier audit, literal Vinberg certification of all 19 rows, the
oriented arithmetic period/change-of-marking classification, and their
replayable synthesis certificate.

## The surface-free class

Let `K3Root[11,20]` be the class of pairs `(N,C)`, modulo oriented integral
lattice isometry.  Here the orientation is a chosen positive-cone component
and `C` lies in that component.  The defining conditions are the following.

1. `N` is an even lattice of signature `(1,rho-1)`, with
   `11 <= rho <= 20`.
2. `N` admits a primitive embedding in

   \[
   \Lambda_{K3}=U^3\oplus E_8(-1)^2.
   \]

3. If

   \[
   \Delta_2(N)=\{r\in N:r^2=-2\},\qquad
   W_2(N)=\langle s_r:r\in\Delta_2(N)\rangle,
   \]

   then `C` is the closure of one component of the positive cone after all
   hyperplanes `r^perp`, `r in Delta_2(N)`, are removed.

This definition mentions the K3 **lattice**, but no surface, effective divisor,
nef cone, period, automorphism, hyperbolic volume, or finite reflection index.
It is therefore an arithmetic root-chamber class rather than a class defined
by already having an algebraic carrier or by already passing the terminal
volume gate.

Define `K3RootFV[11,20]` to be the finite-volume sublocus.  For a complete root
chamber this is equivalent to

\[
[O^+(N):W_2(N)]<\infty,
\]

but that condition is a theorem output below, not a hypothesis in the ambient
class.

## Classification-transfer theorem

The finite-volume unmarked integral members of `K3Root[11,20]`, equivalently
the members of `K3RootFV[11,20]`, are exactly the following 19 cases.  There
is no finite-volume rank-20 member.

| rho | H-dimension | N | disc length | facets | exact-Picard moduli dimension |
|---:|---:|---|---:|---:|---:|
| 11 | 10 | `U+E8(-1)+A1(-1)` | 1 | 12 | 9 |
| 11 | 10 | `U+D8(-1)+A1(-1)` | 3 | 14 | 9 |
| 11 | 10 | `U+D4(-1)^2+A1(-1)` | 5 | 22 | 9 |
| 11 | 10 | `U+D4(-1)+A1(-1)^5` | 7 | 90 | 9 |
| 12 | 11 | `U+A2(-1)+E8(-1)` | 1 | 13 | 8 |
| 12 | 11 | `U+E8(-1)+A1(-1)^2` | 2 | 14 | 8 |
| 12 | 11 | `U+D8(-1)+A1(-1)^2` | 4 | 19 | 8 |
| 12 | 11 | `U+D4(-1)^2+A1(-1)^2` | 6 | 59 | 8 |
| 13 | 12 | `U+E8(-1)+A3(-1)` | 1 | 14 | 7 |
| 13 | 12 | `U+E8(-1)+A1(-1)^3` | 3 | 17 | 7 |
| 13 | 12 | `U+D8(-1)+A1(-1)^3` | 5 | 39 | 7 |
| 14 | 13 | `U+E8(-1)+D4(-1)` | 2 | 15 | 6 |
| 14 | 13 | `U+D8(-1)+D4(-1)` | 4 | 20 | 6 |
| 14 | 13 | `U+E8(-1)+A1(-1)^4` | 4 | 27 | 6 |
| 15 | 14 | `U+E8(-1)+D4(-1)+A1(-1)` | 3 | 21 | 5 |
| 16 | 15 | `U+E8(-1)+D6(-1)` | 2 | 19 | 4 |
| 17 | 16 | `U+E8(-1)+E7(-1)` | 1 | 19 | 3 |
| 18 | 17 | `U+E8(-1)^2` | 0 | 19 | 2 |
| 19 | 18 | `U+E8(-1)^2+A1(-1)` | 1 | 24 | 1 |

For every row the `N`-polarized coarse moduli space is irreducible, and its
nonempty exact-Picard Noether--Lefschetz complement has the dimension shown.
After choosing an oriented marking on a surface in that complement, its full
nef cone is `C`.  Its hyperboloid slice is a finite-sided, finite-volume,
noncompact hyperbolic Coxeter polytope with the stated number of facets.

Conversely, every **finite-volume** member of the ambient surface-free class
occurs in this table.  Ambient members outside the table still have full-nef
K3 carriers, but their complete root chambers are not finite-volume.  Thus the
Nikulin--Vinberg surface/lattice classification transfers back to a genuine
necessary-and-sufficient finite-volume criterion for this arithmetic class.

## Why the count is really 19

There are three logically different counts here, and all three equal 19.

### 1. Neron--Severi lattice types

Roulleau's pinned source states in its introduction that Nikulin and Vinberg
classified the finite-automorphism cases by a list of 118 Neron--Severi
lattices.  Its final table has exactly

\[
4+4+3+3+1+1+1+1+1=19
\]

rows in ranks 11 through 19.  Artebani--Hausen--Laface, Theorem 2.12,
independently records the same rankwise counts:

\[
\begin{array}{c|cccccccccc}
\rho&11&12&13&14&15&16&17&18&19&20\\ \hline
\#&4&4&3&3&1&1&1&1&1&0.
\end{array}
\]

The rank-20 zero makes hyperbolic dimension 18 sharp for the finite-volume
locus of this full-root K3 class.

### 2. Complete integral chamber types

For fixed `N`, the reflection group `W_2(N)` acts transitively on the
components of the complete root complement.  Hence there is one unmarked
integral chamber type per lattice.  Chamber markings and period points are
not being counted as new chamber types.

### 3. Abstract hyperbolic polytope types

The 19 pairs

\[
(\dim\mathbb H,\#\text{facets})
\]

are pairwise different.  In particular, no two lattices in the same dimension
can collapse to the same unmarked abstract polytope.  Different dimensions
cannot collapse either.  Thus the list also contains exactly 19 distinct
hyperbolic polytope isometry types.  This cheap invariant avoids any unsupported
claim that discriminant alone determines normalized hyperbolic geometry.

The total number of facets over the 19 types is 477.

## Structural `U/R/C/V` calibration

This section proves the unmarked full-root carrier and finite-volume
classification statements.  It is not labelled a completion of the project's
strict M3c milestone.  The literal volume and oriented arithmetic layers are
now classified, but the ambient class is not the bounded raw `IC0` Coxeter
specification demanded there, no generators-and-relations presentations of
the 19 infinite arithmetic period groups are supplied, and geometric family
monodromy has not been independently certified.

### `U` — forward universality: closed

Take any ambient member `(N,C)` of `K3Root[11,20]`.  The primitive embedding
`N -> Lambda_K3` and period surjectivity produce a projective K3 surface with
exact Picard lattice `N`; the complement of the Noether--Lefschetz
hyperplanes is nonempty.  If its root chamber has finite volume, the
Nikulin--Vinberg classification proves that its lattice is one of the 19 rows.

Thus every ambient root chamber has an algebraic carrier, and every
finite-volume member enters the classified surface strata.  This is stronger
than realizing a pinned Gram matrix after the fact.

### `R` — exact marked recovery: closed

The K3 Weyl group acts transitively on root chambers.  Composing an
Neron--Severi marking with a Weyl element sends the closed nef cone to the
prescribed `C`.  Extraction recovers the same integral lattice, every
square-minus-two mirror, the chosen positive-cone component, and the chamber.

After forgetting the integral marking, infinitely many period points in the
same family map to one polytope.  The theorem is an exact recovery theorem for
marked lattice/chamber data, not an injective map from surfaces to abstract
polytopes.

### `C` — cone completeness: closed

The K3 ample-cone theorem identifies the ample cone with a component of the
positive cone cut by **all** square-minus-two root hyperplanes.  The simple
walls are exactly the irreducible smooth rational curves on the corresponding
nef boundary.  Roulleau's configurations are published complete curve lists;
the difficult 90-, 59-, and 39-wall rows were completed by induced chambers
inside `II_(1,17)`, not by an unproved degree cutoff.

Consequently no hidden square-minus-two wall is being suppressed.

### `V` — terminal hyperbolic gate: closed structurally and literally

For each table row, the Nikulin--Vinberg finite-automorphism classification
and the K3 cone theorem imply that `W_2(N)` has finite index in `O^+(N)`.
Equivalently, Artebani--Hausen--Laface identify these rows with polyhedral
effective cone, finite automorphism group, and finitely generated Cox ring.
Since `O^+(N)` is an arithmetic lattice and the complete root chamber is a
fundamental domain for `W_2(N)`, finite index proves finite hyperbolic
covolume.  Finite volume is therefore a classification output rather than an
ambient hypothesis.  The integral intersections of distinct irreducible
`(-2)`-curves give the Coxeter relations:

- intersection 0: angle `pi/2`;
- intersection 1: angle `pi/3`;
- intersection 2: parallel at infinity;
- intersection greater than 2: ultraparallel.

Every row contains `U`.  Hence `N` has a rational isotropic class; Weyl
transport moves a primitive effective representative into the nef boundary,
where it is a genus-one fibre class.  Every chamber therefore has an ideal
cusp and is noncompact.

This structural proof is now independently matched by literal
Vinberg/CoxIter audits for all 19 complete wall systems.  Every elliptic
rank-`(n-1)` edge has exactly two finite or ideal endpoint extensions, and the
full zero-neighbourhood of every connected affine component is parabolic of
rank `n-1`:

| package | cases | facets | elliptic edges | finite vertices | ideal vertices | affine audits |
|---|---:|---:|---:|---:|---:|---:|
| ranks 11--14 | 14 | 375 | 106,630 | 12,854 | 1,119 | 4,177 |
| ranks 15--19 | 5 | 102 | 11,215 | 1,244 | 37 | 78 |
| **all** | **19** | **477** | **117,845** | **14,098** | **1,156** | **4,255** |

Both failure counts are zero in every row: no bad elliptic edge and no bad
full parabolic neighbourhood.  All rows have ideal vertices, so the literal
audit also confirms finite volume **and** noncompactness rather than silently
conflating noncompactness with infinite volume.

One source warning is worth preserving.  For the 22-wall
`U+D4(-1)^2+A1(-1)` row, the literal long-segment-as-chain reading of the
published TeX picture has determinant 128 instead of the named lattice's 32,
and makes the printed `A13` have square `-14`.  That transcription lane is
rejected.  The archived certificate instead reconstructs the complete 22-root
system independently inside the standard lattice by exact norm-two Vinberg
height enumeration; its Gram matrix has rank 11 and inertia `(1,10,11)`, and
its literal audit has 1,343 elliptic edges, 185 finite vertices, 25 ideal
vertices, 72 affine audits, and zero failures.  Nothing is guessed from the
ambiguous drawing.

## Embeddings, markings, and moduli

Let `ell(A_N)` be the minimum number of generators of the discriminant group.
Every row satisfies Nikulin's stable uniqueness inequality

\[
\ell(A_N)+2\le 22-\rho.
\]

The largest lengths by rank are `7,6,5,4,3,2,1,0,1`; at rank 19 the
inequality is the boundary equality `1+2=3`.  Therefore every lattice in the
table has a unique primitive embedding orbit under the full isometry group
`O(Lambda_K3)`.  The oriented refinement is now exact as well: every
orthogonal complement `T` splits a unimodular `U`, whose positive reflection
fixes `N` pointwise and corrects the ambient orientation.  Hence every full
orbit remains one orbit under `O^+(Lambda_K3)`; none of the 19 rows doubles.

Fix a positive-cone component and a complete chamber `C`.  The two positive
cone components form two `O^+(N)`-orbits, while `W_2(N)` acts simply
transitively on the chambers in either fixed component.  Consequently there
is one complete-chamber orbit under `O^+(N)`, even though the raw fixed-marking
chamber set is an infinite `W_2(N)`-torsor.

The arithmetic period/change-of-marking images are also classified.  With
`T=j(N)^perp`, the pointwise `N`-polarized ambient group has `T`-image

\[
\widetilde O^+(T)=\ker\bigl(O^+(T)\longrightarrow O(q_T)\bigr).
\]

If the marking of `N` is forgotten but its nef chamber is retained, write

\[
\Gamma^+_{j,C}=\{g\in O^+(\Lambda_{K3}):g(jN)=jN,\ g(jC)=jC\}.
\]

Then

\[
1\longrightarrow\widetilde O^+(T)
\longrightarrow\Gamma^+_{j,C}
\longrightarrow\operatorname{Stab}_{O^+(N)}(C)
\longrightarrow1,
\]

and the `T`/period image of `Gamma^+_{j,C}` is all of `O^+(T)`.  This last
sentence is an image statement, **not** an isomorphism
`Gamma^+_{j,C}=O^+(T)`.  The discriminant representations on both `N` and `T`
are surjective; the resulting finite quotient orders are computed rowwise and
range from 1 to 40,320.

These are oriented **arithmetic** period/change-of-marking groups.  Identifying
them with geometric family monodromy is conditional on the standard
moduli-stack/deck-group interpretation and is not independently certified by
this package.  Weighted diagram automorphisms, integral chamber stabilizers,
generic surface automorphisms, arithmetic period images, and geometric
monodromy remain distinct objects.

Roulleau's Proposition `The-118-moduli-irred` proves that the coarse
lattice-polarized moduli space is irreducible in all 118 cases, hence in these
19.  The nonempty exact-Picard period complement has complex dimension

\[
20-\rho,
\]

giving four 9-dimensional families, four 8-dimensional families, three
7-dimensional families, three 6-dimensional families, and one family in each
dimension 5, 4, 3, 2, and 1.

The irreducibility assertion applies to the full coarse `N`-polarized moduli
space.  The exact-Picard locus is a nonempty complement of a countable
Noether--Lefschetz arrangement with the same dimension; this report does not
promote that complement to a separate finite-type irreducible algebraic moduli
space.

The certificate records full and oriented embedding uniqueness, chamber
orbits, both arithmetic `T`/period images, and coarse irreducibility.  It does
not pretend that an arithmetic image name is a generators-and-relations
presentation of the corresponding infinite group or an independently proved
geometric-monodromy packet.

## What this proves, and what it does not

The result gives a clean positive and negative theorem inside the ambient
class `K3Root[11,20]`:

\[
\boxed{
\begin{gathered}
\text{exactly 19 complete K3 }(-2)\text{-root chambers in the finite-volume locus}\\
\text{occur in }\mathbb H^{10},\ldots,\mathbb H^{18};\\
\text{all are noncompact, and none occurs in }\mathbb H^{19}.
\end{gathered}}
\]

More explicitly, for `(N,C)` in `K3Root[11,20]`, the following are equivalent:

1. `C` has finite hyperbolic volume;
2. `N` is isometric to one of the 19 table entries;
3. an exact-Picard K3 carrier has finite automorphism group;
4. its effective cone is rational polyhedral; and
5. its Cox ring is finitely generated.

Surface classification has therefore produced a finite-volume decision
criterion and finite list, not merely examples.  The literal `V` and oriented
arithmetic embedding/chamber/period-image layers are now closed for all 19.
Strict M3c remains open because the ambient rank-bounded class does not supply
strict `IC0`'s finite raw Gram/nerve/discriminant bounds, the 19 infinite
arithmetic groups have not been given explicit generators-and-relations
presentations, and geometric family monodromy is not independently certified.

It does **not** classify:

- arbitrary hyperbolic Coxeter polytopes;
- incomplete root subsystems;
- mixed-norm divisor reflections;
- invariant/eigen-nef sections;
- rational-surface or hyperkahler carriers; or
- compact polytopes outside the full-nef K3 class.

No worldwide novelty is claimed for the 19 abstract polytopes.  The theorem is
best viewed as a large, non-singleton unmarked classification/finiteness
calibration of the proposed "classify softer surfaces, infer rigid polytopes"
strategy.  Strict M3c remains open, including its bounded `IC0` input.

## Replayable evidence

- [classification manifest](../data/enumerations/k3_high_rank_full_root_classification.json)
- [synthesis certifier](../scripts/certify_k3_high_rank_full_root_classification.py)
- [regression tests](../tests/test_k3_high_rank_full_root_classification.py)
- [rank-11--14 pinned source inventory](roulleau_rank11_14_nef_stock_scout.md)
- [rank-11--14 literal Vinberg audits](k3_rank11_14_literal_vinberg_volume.md)
- [rank-15--19 exact chambers](k3_high_rank_nef_chambers.md)
- [oriented arithmetic period groups](k3_high_rank_oriented_monodromy.md)
- [strong inverse-carrier theorem](k3_strong_inverse_carrier_transfer_theorem.md)

Pinned source witnesses:

- Roulleau TeX SHA-256
  `e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3`:
  classification statement near lines 203--212, irreducible-moduli proposition
  and proof near lines 1129--1227, and the 19 high-rank rows at lines
  10906--10966;
- Artebani--Hausen--Laface PDF SHA-256
  `44553cad30080e7aa4feb32b13147d5f689261bbee0adbd95d53869f483f6f47`:
  Theorems 2.7, 2.11, and 2.12.

The certifier checks the pinned source hashes, table witnesses, frozen case
data, all 19 literal audit summaries, the oriented arithmetic package, and
synthesis invariants.  It is a consistency/replay certificate; the
Nikulin--Vinberg exhaustion and Roulleau irreducibility theorem remain the
mathematical source inputs rather than being reproved by the script.

Run:

```bash
python3 scripts/certify_k3_high_rank_full_root_classification.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_k3_high_rank_full_root_classification.py
```
