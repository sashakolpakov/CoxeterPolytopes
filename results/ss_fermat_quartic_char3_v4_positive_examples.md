# Two finite-volume Coxeter polytopes from the characteristic-3 Fermat quartic

## Outcome

The complete Kondo--Shimada chamber of the supersingular Fermat quartic in
characteristic three has two noncyclic chamber-stabilizer fixed sections above
dimension nine.  Both sections are genuine finite-volume, noncompact
hyperbolic Coxeter polytopes:

| fixed subgroup | fixed rank | hyperbolic dimension | facets | finite vertices | ideal vertices | hyperideal rays | decision |
|---|---:|---:|---:|---:|---:|---:|---|
| `V4-A` | 12 | 11 | 328 | 225,408 | 11,554 | 0 | finite-volume, noncompact |
| `V4-B` | 11 | 10 | 123 | 37,451 | 1,732 | 0 | finite-volume, noncompact |

These are fixed sections of the published 5,944-wall chamber, not chambers
obtained by assuming that an arbitrary surface automorphism group is
reflective.  The algebraic surface supplies the exact Neron--Severi lattice,
the complete wall stock and the finite chamber symmetry; the hyperbolic claims
are then checked independently from the restricted inequalities.

## Exact lineage

The primary source is Shigeyuki Kondo and Ichiro Shimada, *The automorphism
group of a supersingular K3 surface with Artin invariant 1 in characteristic
3*, arXiv:1205.6520v2.  The exact data file is Kondo--Shimada's
`compdataAutFQChar3.txt`, archived in this repository under
`data/external/kondo_shimada_fermat_quartic_char3/`.

The 112 published line classes span the rank-22 Neron--Severi lattice of
determinant -9.  Their incidence graph reconstructs the faithful projective
unitary action of order 13,063,680.  Adjoining arithmetic Frobenius gives the
full chamber stabilizer of order 26,127,360, in agreement with the paper's
split-extension theorem.

An exact GAP conjugacy census has 61 element classes.  The only five
nonidentity classes with fixed rank greater than ten are involutions.  Closing
all distinct high-rank pairs gives exactly two conjugacy classes of noncyclic
subgroups:

- `V4-A` has involution fixed-rank profile `{14,16,16}`, common fixed rank 12,
  conjugacy-class size 5,670, normalizer order 4,608 and centralizer order
  2,304.
- `V4-B` has profile `{12,16,16}`, common fixed rank 11, conjugacy-class size
  4,536, normalizer order 5,760 and centralizer order 2,880.

Neither class extends to an elementary abelian group of order eight with
common fixed rank greater than ten.  Since every element of a subgroup fixing
more than ten dimensions must itself lie in the five-class involution list,
the five cyclic casts and these two `V4` casts exhaust all subgroup-fixed
sections above `H9`.

The primitive fixed lattices have the following determinant and Smith data:

- `V4-A` (rank 12): determinant `-768`, Smith diagonal
  `1^6,2^4,4,12`;
- `V4-B` (rank 11): determinant `1024`, Smith diagonal
  `1^2,2^8,4`.

For `V4-B`, all 123 facet reflections act integrally on the fixed lattice.
For `V4-A`, 292 of the 328 facet reflections act integrally on the saturated
fixed divisor lattice and 36 have denominator two there.  This is not a
defect of the Coxeter polytope: the minimal rational reflection-stable
overlattice is obtained after one closure step, contains the fixed divisor
lattice with index four, and is even integral of signature `(1,11)`,
determinant 48 and Smith diagonal `1^8,2^3,6`.  All 328 reflections act
integrally on this stable lattice.  Thus the `V4-A` finite-covolume group is
arithmetic over `Q`; no identification of the determinant-48 lattice with a
named standard lattice is asserted.

## Complete wall and Coxeter reductions

For `V4-A`, the 5,944 walls form 1,792 subgroup orbits and 1,564 distinct
primitive restricted rays.  Exact conic reduction leaves 328 facets and
certifies the other 1,236 inequalities redundant.  All 53,628 facet pairs are
Coxeter-compatible.  The facet-square histogram is
`(-8/3)^180,(-2/3)^96,(-4)^24,(-2)^16,(-8)^12`.

For `V4-B`, the walls form 1,725 subgroup orbits and 1,259 distinct restricted
rays.  Exact conic reduction leaves 123 facets and certifies the other 1,136
inequalities redundant.  All 7,503 facet pairs are Coxeter-compatible.  Its
facet-square histogram is `(-4)^113,(-2)^10`.

Every retained facet has an exact timelike strict relative-interior witness;
every removed inequality has an exact positive conic dependence on retained
facets.  Thus the finite wall lists used by the volume computations are
complete and irredundant.

## Volume certificates

The GMP rational double-description computation for `V4-B` enumerates 39,183
projective extreme rays: 37,451 timelike, 1,732 null and no spacelike ray.  A
separate literal Vinberg audit finds 246,610 elliptic rank-nine edges, all with
exactly two endpoints:

- finite--finite: 140,140;
- finite--ideal: 94,230;
- ideal--ideal: 12,240.

All 1,732 full affine zero-neighborhoods are valid rank-nine parabolic
neighborhoods.  Their active-facet counts are `12:1092`, `13:585`, `14:55`.
This independently proves finite volume and noncompactness.

For `V4-A`, GMP double description enumerates 236,962 projective extreme
rays: 225,408 timelike, 11,554 null and no spacelike ray.  Hence this chamber
too has finite volume and is noncompact.  Its ideal active-facet counts are
`12:576`, `13:9168`, `14:1680`, `15:112`, `16:18`.

## The H11 reflection lattice is larger than the divisor-fixed lattice

For the 328-wall `V4-A` chamber, 292 primitive mirrors reflect the saturated
fixed divisor lattice integrally, but the 36 type `(-24,div6)` mirrors have
reflection denominator two.  This does not affect the Coxeter or finite-volume
proof; it means only that the Coxeter group does not preserve that particular
integral lattice.

The exact one-step closure of the divisor lattice under all 328 reflections
is already stable.  It is the canonical minimal reflection-stable rational
overlattice containing the divisor lattice, has index four, is even integral
of signature `(1,11)`, determinant 48, and Smith diagonal `1^8,2^3,6`.  In
this lattice all 328 reflections are integral, with primitive types

`276*(-6,div3) + 28*(-2,div1) + 24*(-4,div2)`.

Consequently the 328-wall group is arithmetic over `Q`.  An explicit
determinant-minus-one basis map now identifies this stable lattice with
Scharlau--Kirschmer List11/Gram11 entry 20.  Its known full integral
reflection chamber has 22 walls, 543 finite vertices, and 22 cusps, and all
328 mirrors transport to integral roots of that lattice.  Thus the 328-wall
group is a finite-index reflection subgroup of the known 22-wall group.
Priority of the literal subgroup chamber and divisor realization remains
open; see the
[`exact lineage report`](ss_fermat_char3_v4_h11_lineage.md).

## Reproducible artifacts

- Full stabilizer and subgroup representatives:
  `data/enumerations/ss_fermat_quartic_char3_full_stabilizer_exact.json`.
- Complete wall reduction for both casts:
  `data/enumerations/ss_fermat_quartic_char3_v4_fold_screen.json`.
- H10 exact cone, rays and summary:
  `data/enumerations/ss_fermat_quartic_char3_v4_h10_exact_cone.ine`,
  `..._rays.jsonl.gz`, `..._certificate.json`.
- H10 positive Gram, cusp audit and literal Vinberg certificate:
  `data/enumerations/ss_fermat_quartic_char3_v4_h10_coxeter_screen.json`,
  `..._component_audit.json`, `..._vinberg_certificate.json`.
- H11 exact cone, rays and summary:
  `data/enumerations/ss_fermat_quartic_char3_v4_h11_exact_cone.ine`,
  `..._rays.jsonl.gz`, `..._certificate.json`.
- H11 fixed-lattice reflection denominators and the minimal index-four stable
  overlattice:
  `data/enumerations/ss_fermat_quartic_char3_h11_mirror_lattices.json`.
- H11 stable-lattice and full-reflection-group lineage:
  `data/enumerations/ss_fermat_quartic_char3_v4_h11_lineage.json`.

The reconstruction and certification scripts are
`scripts/scratch_fermat_char3_full_stabilizer.g`,
`scripts/screen_ss_fermat_quartic_char3_v4_folds.py`,
`scripts/certify_ss_fermat_quartic_char3_v4_exact_cones.py`, and
`scripts/certify_ss_fermat_quartic_char3_v4_h10_literal_vinberg.py`.  The H11
mirror-lattice reproducer is
`scripts/certify_ss_fermat_char3_h11_mirror_lattices.py`.

## Novelty boundary

Kondo--Shimada publish the ambient Neron--Severi lattice, complete induced
chamber and its stabilizer, but do not state either `V4` fixed Coxeter section
or either finite-volume polytope above.  The pinned official CoxIter corpus
has global maximum 50 facets, so neither the 123-facet H10 graph nor the
328-facet H11 graph occurs literally in that 1,214-file snapshot.  This is a
bounded exact-match statement, not yet a claim of worldwide first
publication, lattice noncommensurability, or a new arithmetic
commensurability class.
