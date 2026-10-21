# Six finite-volume fixed sections from three Alexeev chambers

## Theorem

Consider Alexeev's complete Coxeter chambers for the even hyperbolic
2-elementary lattices `(18,4,0)`, `(15,11,1)`, and `(15,7,1)`.  Their exact
weighted wall-Gram automorphism groups have orders respectively

`1152, 72, 120`,

and have respectively

`221, 26, 19`

subgroup conjugacy classes.  Exactly six proper classes have fixed rank at
least 11.  Every one of their fixed chambers is an integral hyperbolic
Coxeter polytope of finite volume, and every one is noncompact.

| ID | subgroup | space | facets | elliptic edges | finite vertices | cusps |
|---|---:|---:|---:|---:|---:|---:|
| `A18_4_0_C03` | `C2` | `H11` | 28 | 8,851 | 1,252 | 64 |
| `A18_4_0_C04` | `C2` | `H10` | 28 | 6,956 | 1,087 | 74 |
| `A18_4_0_C06` | `C2` | `H13` | 31 | 28,766 | 3,433 | 238 |
| `A18_4_0_C18` | `V4` | `H10` | 21 | 2,088 | 329 | 22 |
| `A15_11_1_C02` | `C2` | `H10` | 21 | 4,159 | 607 | 50 |
| `A15_7_1_C01` | `C2` | `H10` | 27 | 6,221 | 895 | 67 |

The literal priority conclusion is deliberately narrower.  The first five
diagrams have no exact match in all 18 explicit full matrices in Alexeev's
ancillary file or in the pre-existing positive catch ledger.  They are
therefore newly extracted here and unmatched under this bounded comparison,
not proved new in the worldwide literature.  Indeed, Allcock's 2006 theorem
[*Infinitely many hyperbolic Coxeter groups through dimension
19*](https://doi.org/10.2140/gt.2006.10.737) gives infinitely many
finite-volume Coxeter polytopes in every dimension at most 19.  Thus none of
these `H10`, `H11`, or `H13` examples can be first in its dimension, and a
bounded corpus absence cannot establish worldwide literal-diagram priority.
The sixth is not a new abstract polytope: an explicit facet permutation
identifies it with the already certified
`D0004/C007 = D16-D4-C7-H10-27` chamber.

## Why the volume conclusion is complete

For each fixed space, wall orbits were restricted to the saturated fixed
lattice.  Exact nonnegative-dependence certificates remove all conically
redundant inequalities, and every retained inequality has an independent
strict timelike facet witness.  The resulting six facet Grams have no
forbidden Coxeter pair, and all primitive mirror reflections act integrally
on the saturated fixed lattice.

The volume certificate then applies Vinberg's literal criterion, rather than
inferring volume from a finite wall list.  It enumerates

1. every elliptic rank-`n-1` edge diagram;
2. every elliptic rank-`n` finite endpoint;
3. every irreducible affine special subdiagram, including parallel affine
   `A1` pairs; and
4. the full zero-neighborhood of every affine component.

For all six diagrams, every edge has exactly two finite-or-ideal endpoints
and every affine full zero-neighborhood is a rank-`n-1` parabolic cusp.
There are zero bad edges and zero bad affine neighborhoods.  The positive
cusp count in every row proves noncompactness.

The archived `edge_search_prefix_counts` are deliberately labelled as
ordering-dependent search diagnostics.  They are not face numbers.  The
`finite_special_subset_counts_by_rank` are the invariant elliptic-special
subset counts used by the theorem.

## Fixed lattices and completed reflection-stable lineage

| ID | saturated fixed lattice | unique maximal reflection-stable integral overclosure | containing full chamber / index |
|---|---|---|---|
| `A18_4_0_C03` | rank 12, determinant `-256`, Smith `1^6 2^4 4^2` | index `4`, odd determinant `-16`, Smith `1^8 2^4` | Scharlau--Kirschmer `List11` entry 15, 15 walls; finite index, numerical value open because both `H11` characteristics vanish |
| `A18_4_0_C04` | rank 11, determinant `1024`, Smith `1^2 2^8 4` | index `2`, odd determinant `256`, Smith `1^3 2^8` | `List10` entry 11, 13 walls; exact index `130560` |
| `A18_4_0_C06` | rank 14, determinant `-64`, Smith `1^10 2^2 4^2` | index `4`, even determinant `-4`, Smith `1^12 2^2` | `List13` entry 4, 15 walls; finite index, numerical value open because both `H13` characteristics vanish |
| `A18_4_0_C18` | rank 11, determinant `128`, Smith `1^6 2^4 8` | index `2`, even determinant `32`, Smith `1^6 2^5` | `List10` entry 20, 14 walls; exact index `144` |
| `A15_11_1_C02` | rank 11, determinant `8192`, Smith `2^9 4^2` | index `4`, even determinant `512`, Smith `1^2 2^9` | `List10` entry 6, 13 walls; exact index `6720` |
| `A15_7_1_C01` | rank 11, determinant `512`, Smith `1^4 2^5 4^2` | index `4`, even determinant `32`, Smith `1^6 2^5` | `List10` entry 20, 14 walls; exact index `2016` |

Here “maximal” has a precise restricted meaning: maximal among integral
overlattices preserved by every displayed facet reflection.  The exhaustive
discriminant-module computation closes all invariant integral orbit closures
under addition and proves uniqueness in that scope.  Global rescaling and
scaled duality are not silently treated as overlattices: for example,
`A15_11_1_C02` has Gram content two and its globally halved form is explicitly
congruent to `U(2)+<-1>^9`, while the `A18_4_0_C04` identification with
`List10` entry 11 uses the scaled dual.

This closes the arithmetic lineage without producing a new commensurability
class.  Exact local Hasse comparison puts all six in the relevant classical
odd-unimodular rational class, and every maximal stable carrier is identified
with a pre-existing Scharlau--Kirschmer chamber.

For the odd-dimensional `H11` and `H13` cases, finite index needs no volume
ratio: the certified cofinite reflection group is a subgroup of the
arithmetic integral orthogonal group, and both have finite covolume, so the
index is finite.  The index value remains open because the two ordinary
Coxeter Euler characteristics are both zero.

## Exact duplicate and bounded exclusions

For `A15_7_1_C01`, the following one-based permutation sends its 27 facets
to the archived `D0004/C007` facets and makes every signed normalized Gram
entry equal:

```text
[3,1,9,2,4,7,12,17,5,8,13,18,11,16,10,21,14,6,15,23,19,20,22,24,25,26,27]
```

The two 21-wall `H10` diagrams are mutually nonisomorphic.  The 28-wall
`H11` and 31-wall `H13` diagrams have no same-dimension/facet-count entry in
the pinned ordinary CoxIter repository or its index-two doubling archive.
The two 21-wall `H10` diagrams and the 27-wall `H10` diagram are likewise
absent by header.  The doubling archive does contain two `H10/28` records,
but complete edge-coloured relation-graph comparison proves that neither is
`A18_4_0_C04`.  The five bounded-priority diagrams also have no match among
Alexeev's 18 explicit full chamber matrices.  These are useful exclusions,
but no search snapshot is a worldwide priority theorem.

## Algebraic-geometric meaning and its limit

The `(18,4,0)` and `(15,7,1)` ambient lattices lie on Alexeev's `r+a=22`
K3 line and occur as Neron--Severi lattices of complex 2-elementary K3
surfaces.  Hence five of the six constructions genuinely start inside
algebraic divisor lattices, rather than in an arbitrary rational Lorentz
space.

That statement is not automatically the stronger quotient-nef statement.
Alexeev's 48- and 60-wall chambers are mixed-root full-reflection chambers
and include `(-4)` mirrors.  Separate exact discriminant-gluing arguments
now close the geometric lifting problem for all four selected `(18,4,0)`
classes.  `A18_4_0_C06` is induced on `NS(X)` by a genuine non-symplectic
order-four K3 automorphism; see the
[order-four lift theorem](a18_c06_h13_31_k3_order4_lift.md).  Independently,
its equation-level `U+E8(-1)+D4(-1)` carrier now has a complete Weyl framing:
all 26 `(-2)` mirrors are effective, and the primitive square-90666 class
`D506+7H` is ample and strict on all 31 walls.  See the
[nef/effectivity theorem](a18_c06_ue8d4_nef_effectivity.md).  The `C03` and
`C04` actions are realized by symplectic Nikulin involutions at special
periods of Picard ranks 20 and 19, while their fixed Neron--Severi lattices
remain exactly the archived rank-12 and rank-11 lattices; see the
[two-case special-period theorem](a18_c03_c04_special_k3_symplectic_lifts.md).
The `C03 H11/28` chamber now has a second full-Picard carrier.  Its stable
odd index-four closure has doubled dual
`U+D4(-1)+A1(-1)^6`, the even 2-elementary lattice `(12,8,1)`.  An explicit
256-class discriminant graph glue proves a primitive embedding in the K3
lattice, so a very general point in the resulting eight-dimensional period
family has this exact Neron--Severi lattice.  All 28 oriented mirrors
transport and generate it.  This is a full-Picard divisor carrier rather
than another automorphism-fixed or full-nef realization.  It is now attached
to an explicit eight-moduli `I0*+6I2+6I1` Weierstrass family with trivial
very-general Mordell--Weil group; a determinant-minus-one fiber marking
carries all 28 walls.  See the
[scaled-dual theorem](a18_c03_ud4a1x6_scaled_dual_carrier.md) and
[equation-level theorem](a18_c03_ud4a1x6_weierstrass_carrier.md).
The `C04 H10/28` chamber now has a second and arithmetically different
carrier as well.  The doubled dual of its unique reflection-stable odd
index-two closure is explicitly `U+D8(-1)+A1(-1)`, the full Picard lattice
of an actual rank-11 K3 stock.  All 28 oriented halfspaces and integral
reflections transport, with primitive split
`16*(-2,div1)+8*(-2,div2)+4*(-4,div2)`, and the normals generate the full
target lattice.  This is a mixed divisor carrier, not a second automorphism-
fixed or full-nef realization.  It is now attached to an explicit
nine-parameter `I4*+I2+12I1` Weierstrass family with exact very-general
`NS=U+D8(-1)+A1(-1)`.  A deterministic Weyl reduction gives small
fiber-basis coordinates for all 28 mirrors and an explicit square-1408 ample
class lying strictly inside every imported halfspace.  Thus this second
carrier is equation-level and has nonempty intersection with the ample cone,
without being identified with the full nef chamber; see the
[equation-level theorem](a18_c04_ud8a1_weierstrass_carrier.md).
Finally, the `C18` Picard `V4` is the restriction of a genuine `C4 x C2`
action on a special Picard-rank-20 K3; its common fixed `NS` remains the
archived rank-11 lattice.  See the
[C18 finite-group lift theorem](a18_c18_special_k3_c4xc2_lift.md).
Their `H13/31`, `H11/28`, and two `H10` chambers are therefore genuine
automorphism-fixed divisor/reflection-chamber sections.  The source 48-wall
chamber is still not asserted to be the literal nef cone, and no
quotient-surface nef interpretation follows.  Exact seven-bit discriminant
gluing now also realizes `A15_7_1_C01` by a genuine symplectic Nikulin
involution on a special Picard-rank-19 K3.  Its fixed `NS` is exactly the
archived rank-11 lattice, so its `H10/27` section is unchanged; see the
[special-period lift theorem](a15_7_c01_special_k3_symplectic_lift.md).
The source 60-wall chamber is likewise not asserted to be the literal nef
cone.  Thus every selected action from the two K3-admissible ambient stocks
now has a geometric lift.  The `(15,11,1)` source lattice itself remains
outside the complex K3 Picard range, but its rank-11 section no longer does.

The `(15,11,1)` stock has `r+a=26`; it is not a complex K3 Picard lattice in
Nikulin's K3 range.  Nevertheless, the 21-wall section's unique even
reflection-stable index-four overclosure is
`U+E8(-2)+<-2>`.  A very general K3 with this full Néron--Severi lattice
carries all 21 reflections, and the same lattice occurs as the saturated
invariant divisor lattice of a rational involution on Shimada's `X_2`.
Thus the section now has a complex-K3 divisor realization even though its
66-wall source stock does not; see the
[two-tier carrier theorem](a15_11_h10_21_x2_k3_carrier.md).

That Picard carrier is now explicit at equation level.  A normalized
nine-parameter rational elliptic surface with twelve `I1` fibers, followed
by quadratic base change ramified over one `I1` and one smooth fiber, gives
a K3 with `I2+22I1`.  Fiber/section, the pulled-back `E8` Mordell--Weil
lattice, and the new `I2` component give a saturated
`U+E8(-2)+<-2>`; the nine-dimensional moduli calculation makes it the full
very-general Néron--Severi lattice.  Hence all 21 walls attach to an explicit
algebraic family; see the
[base-change theorem](a15_11_ue8m2a1_basechange_family.md).

## Reproduction

The permanent files are:

- [construction and full subgroup denominator](../data/enumerations/alexeev_remaining_high_fixed_section_construction.json), generated by
  [certify_alexeev_remaining_fixed_sections.py](../scripts/certify_alexeev_remaining_fixed_sections.py);
- [literal Vinberg volume certificate](../data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json), generated by
  [certify_alexeev_remaining_fixed_sections_vinberg.py](../scripts/certify_alexeev_remaining_fixed_sections_vinberg.py) and the shared exact C++ enumerator;
- [bounded exact lineage certificate](../data/enumerations/alexeev_remaining_high_fixed_section_lineage.json), generated by
  [certify_alexeev_remaining_fixed_sections_lineage.py](../scripts/certify_alexeev_remaining_fixed_sections_lineage.py); it freezes the rational classes, Solomon characteristics, exact 27-wall duplicate, nonisomorphism of the two 21-wall diagrams, and the pinned Alexeev/CoxIter/prior-ledger exclusions;
- [exhaustive reflection-stable-overlattice certificate](../data/enumerations/alexeev_remaining_high_fixed_section_stable_overlattices.json), generated by
  [certify_alexeev_remaining_stable_overlattices.py](../scripts/certify_alexeev_remaining_stable_overlattices.py); and
- [regression tests](../tests/test_alexeev_remaining_fixed_sections.py) and
  [stable-overlattice tests](../tests/test_alexeev_remaining_stable_overlattices.py).
- [equation-level `A15_11 H10/21` carrier](a15_11_ue8m2a1_basechange_family.md),
  [`A18-C06 H13/31` nef/effectivity framing](a18_c06_ue8d4_nef_effectivity.md),
  [scaled-dual `A18-C03 H11/28` carrier](a18_c03_ud4a1x6_scaled_dual_carrier.md),
  [equation-level `A18-C03 H11/28` carrier](a18_c03_ud4a1x6_weierstrass_carrier.md),
  and [equation-level `A18-C04 H10/28` carrier](a18_c04_ud8a1_weierstrass_carrier.md),
  each with its linked generator, JSON certificate, and focused tests.

Run from the repository root with `PYTHONPATH=src:.`:

```bash
python3 scripts/certify_alexeev_remaining_fixed_sections.py --check-archive
python3 scripts/certify_alexeev_remaining_fixed_sections_vinberg.py --check-archive
python3 scripts/certify_alexeev_remaining_stable_overlattices.py --check-archive
python3 scripts/certify_alexeev_remaining_fixed_sections_lineage.py --check-archive
pytest -q tests/test_alexeev_remaining_fixed_sections.py tests/test_alexeev_remaining_stable_overlattices.py
```

The construction check deliberately re-enumerates all three full weighted
groups and all 266 subgroup conjugacy classes; it is slower than the archive
tests because it is the denominator proof, not a hash-only check.
