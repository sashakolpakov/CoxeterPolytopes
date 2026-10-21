# TOOT: a complete S16 surface transfer and an S17 countercalibration

Status: **full marked `U/R/C/V` theorem for the oriented surface-free S16
class PASS** (2026-08-08).  The parallel S17 class disproves the same theorem
for the minimal `E8(-2)` K3 carrier alphabet.

## The class that actually closes

Let `D16+` be the following raw Coxeter datum.

- There are fourteen labelled facets and the unlabelled nonorthogonality
  support is `S16` from the variable D34 cut denominator.
- Exactly eight support edges have relation `m4`; all other support edges have
  relation `m3`, and all nonedges have relation `m2`.
- Primitive integral wall norm squares lie in `{2,4}`, with exactly ten
  norm-square-four and four norm-square-two facets.
- The Gram target is rank ten with negative-root inertia `(1,9,4)`.
- The positive sheet is the component containing the generated strict exact
  cone interior.

This definition names no Gram, surface, or known polytope.  The
`short4/long10` condition is essential: it is a surface-free arithmetic
selector that distinguishes the geometric lane from its scaled dual.

The variable census starts with all 12,870 eight-`m4` labellings on `S16`,
reduces them to 376 support-automorphism orbits, proves that 136 orbits admit
primitive `{2,4}` normalizations, and checks all 272 scaled-dual exact Gram
lanes.  After the rank/inertia and ten-long selector, `D16+` has exactly one
marked integral Gram orbit.  Its relation pattern is derived to be T2; T2 is
not part of the input.

The finite package list `A(D16+)` has sixteen even reflection-stable fixed
overclosures.  Fifteen fail the distinct-class, `E8/2E8`, or K3-complement
gate.  The unique survivor is the exact `(18,4,0)` package with

```text
fixed lattice       U + E8(-2),
coinvariant          E8(-2),
K3 complement        U(2) + U(2),
source roots         4 fixed + 10 exchanged pairs,
root-span index      1.
```

## The intrinsic surface stratum

Define `M_D16+^alg` to parametrize marked projective K3 surfaces with exact
Picard lattice equal to this rank-18 package, the specified symplectic
Nikulin action, the marked effective root orbits, and

```text
fourteen-wall section = Nef(X) intersect Fix(g).
```

This is an intrinsic algebraic-geometric cone condition.  “Finite hyperbolic
volume,” “passes Vinberg,” and the name T2 are not in its definition.  The
existing K3 period/Torelli and complete invariant-nef theorem proves this
stratum nonempty.

## Why `U/R/C/V` are all closed

- **`U`.**  The surface-free enumeration forces one labelled integral Gram
  and one viable package.  Two spanning labelled root systems with the same
  Gram induce an integral isometry of their generated root lattices; the
  generated strict interior fixes the positive component.  Hence every
  realization of `D16+`, modulo labelled orientation-preserving integral
  isometry, transports to the nonempty marked K3 stratum.  This is the
  missing forward-surjectivity step that the formerly pinned T2 bridge could
  not provide by itself.
- **`R`.**  The exact unimodular bridge transports the rank-18 form, action,
  all ordered 24 roots, all fourteen orbit sums, and the complete 48-wall
  source action.  The generated Gram equals the old T2 input entry-for-entry,
  and the fixed, rank-18, involution, and root hashes agree in exact
  coordinates.
- **`C`.**  The independent strong-promotion theorem gives
  `P_r intersect Fix(g) = Nef(X) intersect Fix(g)` and proves every restricted
  long source wall redundant behind the fourteen active short-root orbit
  walls.
- **`V`.**  The variable compiler, not the old factory transcript, rebuilds
  a strict positive interior and fourteen positive-square facet witnesses.
  Literal Vinberg then gives 683 elliptic rank-eight edges, 90 finite and 20
  ideal endpoints, and zero bad edges or full parabolic neighbourhoods.

Therefore the first fully marked transfer equivalence for this variable raw
denominator is

```text
Pol_fv^mark(D16+) is nonempty
    iff
M_D16+^alg is nonempty.
```

The extraction image has one unmarked Gram/congruence type: the known
finite-volume noncompact T2 polytope in `H9`.  This is a new transfer theorem
and a different surface-first reconstruction of a known polytope, not a new
abstract polytope.

## Why the bare S16 statement must be weaker

If both primitive norm orientations are retained, the complementary
`short10/long4` lane has five stable overclosures but no minimal weak K3
package.  Hyperbolic root rescaling does not change the unmarked polytope, so
the bare S16 class still has the valid unmarked existence corollary

```text
Pol_fv(S16) is nonempty iff M_D16+^alg is nonempty.
```

It does **not** have the advertised marked universality statement for both
scaled-dual integral realizations.  This distinction is required by the
project's rule that scaled duals are separate marked packages.

## S17 is the exact countercalibration

Define `D17+` by replacing `S16` with `S17` and retaining the same
`short4/long10` selector.  Its variable arithmetic is just as rigid:

- exactly one Lorentzian diagram orbit, derived as T1;
- five stable even overclosures;
- four fail before weak K3 realization;
- the unique weak survivor is `(18,2,0)` with complement `U+U(2)`;
- a genuine projective K3 surface and Nikulin action exist;
- exact fold-back `R0` and fresh finite-volume `V` pass, with 112 finite and
  15 ideal endpoints.

But the sole weak survivor has five exchanged-root hidden walls.  The first
already changes sign on two exact candidate-cone rays.  Riemann--Roch forces
one orientation of that wall on every invariant nef section, so the proposed
fourteen-wall cone cannot be the invariant nef section in either effectivity
orientation.  Thus every strong minimal-`E8(-2)` stratum for `D17+` is empty
while its T1 polytope exists:

```text
Pol_fv^mark(D17+) is nonempty,
M_D17+,minimal^alg is empty.
```

This is stronger than saying `U` is merely unknown: it proves that the
minimal K3 carrier alphabet is not universal even in the adjacent two-support
denominator.  It does not exclude a different K3 lift alphabet, another
algebraic surface, or a higher-dimensional algebraic carrier.

The pair `S16/S17` therefore supplies an actual criterion boundary:

```text
cut + Lorentzian + K3 embedding + exact fold + Vinberg
    is not enough;
complete all-root separation is the discriminator between
the strong T2 carrier and the obstructed T1 carrier.
```

## Milestone assessment

This closes **M3b for the sharply scoped oriented class `D16+`**: all four
obligations are independently proved and the algebraic moduli definition
does not contain the hyperbolic conclusion.  It also closes the positive
existence and discrete-Gram-image part of M3c for this class: the image is the
single T2 type.

The broad M3c milestone should remain open.  The present theorem does not
classify every marking/monodromy orbit of the positive-dimensional K3 moduli
stratum, does not cover the scaled-dual marked lane, and does not vary the
support beyond `S16`.  The next honest expansion is to find a surface-free
condition admitting several supports while excluding or repairing the S17
hidden-root mechanism.

## Replay

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_s16_strong_transfer_corollary.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_h9_14_s16_strong_transfer_corollary.py
python3 -m py_compile scripts/certify_h9_14_s16_strong_transfer_corollary.py
```

Artifacts:

- [`h9_14_s16_strong_transfer_corollary.json`](../data/enumerations/h9_14_s16_strong_transfer_corollary.json)
- [`certify_h9_14_s16_strong_transfer_corollary.py`](../scripts/certify_h9_14_s16_strong_transfer_corollary.py)
- [`test_h9_14_s16_strong_transfer_corollary.py`](../tests/test_h9_14_s16_strong_transfer_corollary.py)
- [`variable D34 cut denominator`](h9_14_d34_cut_denominator.md)
- [`T2 exact marked bridge`](a18_c01_t2_marked_bridge_hostile_audit.md)
- [`T2 complete invariant-nef promotion`](a18_c01_c07_k3_nef_promotions.md)
- [`T1 exchanged-root obstruction`](h9_14_t1_strong_k3_nef_obstruction.md)
