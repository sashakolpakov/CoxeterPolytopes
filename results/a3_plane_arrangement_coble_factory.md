# The A3 plane-arrangement Coble search: an exhaustive one-boundary no-go and a surface-only H19 output

## Frozen verdict

The `A3` reflection arrangement gives **no new complete finite-volume
Coxeter chamber** in this search.

- The complete one-boundary search over all proper intersection points and
  allowed first-near mirror directions is finite and exhaustive: `48,822`
  parent-closed ten-center clusters.  Although `9,925` pass the two raw
  Hilbert inequalities, none produces a reduced, nondegenerate one-boundary
  strong Coble surface.
- The full reduced union of the six mirrors does produce an exact strong
  Coble surface over `Q`, after the four triple concurrences are resolved by
  their forced proper-plus-three-near clusters.  It has ten boundary
  components, Picard rank 20, and hyperbolic dimension 19.
- Its obvious 25-curve wall stock has exact cdd profile

  \[
          2,865\text{ finite}+331\text{ ideal}+790\text{ hyperideal}.
  \]

  Hence that natural stock is not a finite-volume polytope.  It is not the
  full nef wall set.
- More strongly, Kondō's exhaustive finite-automorphism classification has no
  ten-boundary strong Coble surface.  A finite rational-polyhedral full nef
  cone would make the Néron--Severi action finite (and Kondō proves the kernel
  finite), hence the automorphism group finite.  Therefore this `b=10`
  surface cannot have a finite-facet complete nef Coxeter chamber.  There is
  no hidden finite completion waiting beyond the 25 obvious curves.

Thus the positive result is a characteristic-zero **surface in H19**, not a
new finite-volume hyperbolic Coxeter polytope.  No `TOOT` is warranted.

## 1. Exact A3 arrangement

Use the six projective lines

\[
 x=0,\quad y=0,\quad z=0,\quad x-y=0,\quad y-z=0,
 \quad x-z=0.
\]

Their intersection lattice has four triple and three double points, with

\[
                  4\binom32+3\binom22=15=\binom62.
\]

The exact projective arrangement group has order 24.

For a proper point `p` and an incident mirror direction `v`, the program uses
the actual proximity conditions, not merely repeated evaluation:

- `-K`: `F(p)=0`, and at a selected first-near point `D_vF(p)=0`;
- `-2K`: `F,D_vF,D_nF=0` at `p`, followed at the near point by
  `D_v^2F,D_v^3F,D_vD_nF=0`.

These are exactly the multiplicity-two conditions after substituting
`p+u v+u s n` and dividing by `u^2`.

## 2. Exhaustive ten-center Hilbert census

A cluster is parent-closed: a first-near flag may be selected only when its
proper parent is selected.  Each of the four triple points has three allowed
near directions and each of the three double points has two.  The coefficient
of degree ten in

\[
  \bigl(1+x(1+x)^3\bigr)^4
  \bigl(1+x(1+x)^2\bigr)^3
\]

is `48,822`, agreeing with the enumeration.

The exact rational ranks leave `9,925` clusters with
`h0(-K)=0` and `h0(-2K)>0`.  Their sextic-rank split is:

| sextic rank | `h0(-2K)` | marked clusters | A3-group orbits |
|---:|---:|---:|---:|
| 24 | 4 | 25 | 4 |
| 25 | 3 | 948 | 46 |
| 26 | 2 | 3,678 | 162 |
| 27 | 1 | 5,274 | 230 |

Only the rank-27 stratum can have the unique reduced boundary required for a
one-boundary strong Coble surface.  Exact factorization of its 230 orbit
representatives gives:

| degree partition, with multiplicity | orbits |
|---|---:|
| `1+1+1+1+1+1` | 160 |
| `1+1+1+1+2` | 35 |
| `1+1+1+3` | 16 |
| `1+1+2+2` | 2 |
| `1+1+4` | 13 |
| `1+2+3` | 1 |
| `6` | 3 |

So only three orbit representatives are irreducible over `Q`.  Every one of
those three prescribes two or more first-near double points over at least one
proper point.  Exact local jets show what goes wrong: a quadratic tangent cone
cannot support two distinct double infinitely-near points.  The actual proper
multiplicity is at least four, so the first exceptional occurs in `-2K` with
coefficient at least two.  The bi-anticanonical divisor is nonreduced and does
not define a nondegenerate one-boundary strong Coble surface.

This proves the scoped theorem:

> No parent-closed ten-center cluster supported on the `A3` intersection
> points and first-near mirror flags produces a reduced nondegenerate
> one-boundary strong Coble surface.

It is a carrier no-go, not a no-go for arbitrary ten points in the plane.

## 3. The full six-mirror positive surface

At an ordinary triple concurrence, the six-line boundary has multiplicity
three.  Blow up the proper point.  The first exceptional must appear once in
the reduced `-2K` divisor, and it meets the three strict line branches.  Blow
up those three intersections.  Thus each triple point contributes one proper
and three first-near centers.

The full arrangement has four triple and three double points, hence

\[
                 N=4(1+3)+3=19,
                 \qquad b=6+4=10,
                 \qquad N=9+b.
\]

The ten boundary components are the six strict mirrors and the four first
exceptionals.  Every component is a smooth rational `(-4)`-curve, and they are
pairwise disjoint.  Their sum is

\[
                         6H-2\sum_{i=1}^{19}E_i=-2K.
\]

Negative intersection forces this reduced sum to be the unique member of
`|-2K|`; because it is reduced, it cannot equal twice an anticanonical
divisor.  Therefore `|-K|` is empty and the blowup is a strong Coble surface.

This construction is rigid over `Q`; no characteristic specialization is
being hidden.

## 4. Natural wall cone and the completion obstruction

The obvious negative curves are:

- ten boundary `(-4)`-curves;
- three exceptional `(-1)`-curves over the double points;
- twelve terminal exceptional `(-1)`-curves over the first-near points.

Their 25 classes span all of `NS=I_(1,19)`.  Exact cddlib finds `3,986`
extreme rays:

\[
                  2,865+331+790
\]

of positive, null, and negative square respectively.  All 25 facets are
irredundant.  The certificate also expands each sampled cdd generator back to
the physical `H,E_i` basis; the first has square `-1`, pairs nonnegatively with
all 25 walls, and its active wall normals have rank 19.  Thus the 790
negative-square extreme rays are literal Vinberg/cdd hyperideal witnesses, so
this declared cone has infinite volume.

This alone would leave open whether additional effective curves cut every
hyperideal ray.  The algebraic classification closes precisely that loophole
for a **complete finite-facet nef chamber**.  Kondō classifies finite-
automorphism strong Coble surfaces and obtains boundary counts

\[
                         b\in\{1,2,5,8\}.
\]

The present surface has `b=10`.  Since the kernel of the action on the
Coble--Mukai/Néron--Severi lattice is finite, a finite-facet full nef cone
would force finite automorphism group.  Contradiction.  Therefore any full
nef description here necessarily has infinitely many walls or otherwise
fails the finite Coxeter-chamber contract.

That is the desired primal/dual separation:

- dual/example side: a new exact `A3 ->` strong-Coble surface construction in
  dimension 19;
- primal/theorem side: classification prevents it from becoming a complete
  finite-volume Coxeter polytope.

## 5. Lineage and scope

The general reduced rational-sextic Coble factory was already frozen in this
project.  The exact contribution here is the reflection-arrangement version
with forced resolution of ordinary triple points, the complete ten-center
Hilbert/factor/local census, and the rank-20 cone certificate.  Worldwide
priority for this particular resolved complete-quadrilateral presentation was
not established; no new-surface priority claim is made.

The finite-automorphism input is Shigeyuki Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially Theorem 1.1 and Table 2.

## 6. Replay

```bash
python3 scripts/certify_a3_h3_plane_arrangement_coble_factories.py
pytest -q tests/test_a3_h3_plane_arrangement_coble_factories.py
```

Artifacts:

- [`certify_a3_h3_plane_arrangement_coble_factories.py`](../scripts/certify_a3_h3_plane_arrangement_coble_factories.py)
- [`a3_h3_plane_arrangement_coble_factories.json`](../data/enumerations/a3_h3_plane_arrangement_coble_factories.json)
- [`test_a3_h3_plane_arrangement_coble_factories.py`](../tests/test_a3_h3_plane_arrangement_coble_factories.py)
