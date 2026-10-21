# Four complete next steps in the quartic A9 host

**Superseded status.** The full-host question left open below is now resolved by the [outer edge-walk chamber-stabilizer certificate](fresh_h9_a9_outer_20261013.md): this specified quartic A9 lattice is nonreflective. The finite first steps and geometry recorded here remain valid inputs to that proof.

The full nine-edge star of the previously certified finite corner is now known. Four complete first-mirror calculations, together with the integral A9 reversal, add eight coherent norm-two facets and eight actual finite vertices. The resulting **19-wall H9 prefix still has an exact strict spacelike escape**. This is constructive progress in the full arithmetic chamber, not a new finite-volume polytope or a nonreflectivity theorem.

The host remains

\[
K=\mathbf Q(T),\quad T^4-4T^2+2=0,\quad T=2\cos(\pi/8),
\qquad L=\langle-(T-1)\rangle\oplus A_9.
\]

Its odd inert-prime-five determinant square class distinguishes it from the two previously excluded F4 hosts. The previous [A9 report](fresh_h9_new_host_next_20261013_a9.md) supplies the complete initial A9 stabilizer, the certified first step to the A4+A4+A1 corner, and the structural primitive norm bound. No worldwide novelty assertion is made for this lattice.

## Complete first-mirror calculations

The central corner has incident old facets 0, 1, 2, 3, 5, 6, 7, 8, 10. Releasing facet 10 returns to the initial A9 corner. Reversal pairs the other eight outgoing directions, leaving releases 0–3 to compute.

Every primitive root norm is represented in the safe 18-element pool of positive divisors of 20 modulo unit squares. The finite-angle profile enumeration uses this whole pool, rather than retaining only norm two in advance. For each proposed first contact, all eight angles to the known edge facets must be finite Coxeter angles; the complete profiles reduce the remaining coordinates to a binary relative norm equation.

Put \(\alpha=T-1\). If the old central point is \(p\), the integral vector \(V=zp\) has \(V^2=-20\alpha\). For the released wall, let \(b\) be its residual perpendicular to the eight edge walls. Its norm \(B\) is 5/4 for releases 0 and 3, or 5/6 for releases 1 and 2. Given a finite-angle profile, its fixed component is \(r_0\), and \(R=q-r_0^2\) is positive at all four embeddings. All integral roots in the profile have the form

\[
r=r_0-\frac Z5b+\frac XdV,
\qquad
\begin{cases}
Z^2-\alpha X^2=20R,&d=20,\ B=5/4,\\
Z^2-6\alpha X^2=30R,&d=10,\ B=5/6.
\end{cases}
\]

The saved integral covectors for \(Z\) and \(X\) prove these variables are integral for every lattice root; this is not an assumed coordinate sublattice. Both relative degree-eight extensions have signature (2,3) and certified class number one. The two complete norm/unit interfaces suffice for all four directions.

| Released facet | Complete profiles | Relative targets | Exact interval points | Primitive first mirrors | Endpoint mirror count |
|---|---:|---:|---:|---:|---:|
| 0 | 77 | 37 | 246 | 1 | 21 |
| 1 | 113 | 41 | 204 | 1 | 21 |
| 2 | 113 | 41 | 204 | 1 | 21 |
| 3 | 77 | 37 | 246 | 1 | 21 |

All four sole primitive survivors have norm two. Each endpoint again has complete spherical stabilizer A4+A4+A1. The complete interval contains **900 reconstructed points** in total, including nonintegral, imprimitive, and nonreflective rejections. The saved predecessor and successor inequalities exhaust each relative-unit orbit up to the chosen endpoint. No earlier primitive prefix conflict is discarded.

The finite initial unit-power scout only found an upper endpoint. It does not supply the completeness claim. Completeness uses PARI's certified full integral absolute-norm solution representatives, the full unit bases, the exact relative norm map, and exhaustive ordering of every resulting relative orbit through that upper endpoint. For both extensions the first three unit-norm columns are nonsingular and the fourth has norm one, proving that the latter generates the relative norm-one unit group up to sign. The independent checker replays all of this interface and its exact arithmetic, but does **not** independently recompute PARI's ideal classes or full unit group. The certified PARI interface remains an explicit shared dependency.

Firstness follows from the complete original chamber germ and the finite first-contact argument: any earlier mirror would make a finite spherical vertex on this compact edge segment, hence supply one of the enumerated finite profiles. Every integral solution in each profile has been covered through the endpoint. Exact outward orientations and full reflection closures identify the new endpoint with the same full reflection chamber, rather than merely a compatible reflection subgroup.

## Coherent geometry and remaining queue

The eight new roots and all 171 pairs among the 19 saved walls pass exact integrality and Coxeter-angle gates. Eight rank-nine positive-definite incident Gram matrices and all 19 inequalities certify the new finite vertices. Every new root is perpendicular to the central cap, so these eight steps also advance the geometry of its H8 facet. The shared field here is quartic; this observation is not an identification with a known golden-field H8 polytope.

There are ten certified actual vertices: the initial A9 vertex, the central vertex, and eight new endpoints. Exact lookahead enumerates all 90 outgoing directions. Eighteen are the two orientations of the nine known edges. In each of the other **72 directions**, the light cone is reached before any existing saved wall; reversal leaves **36 unresolved edge orbits**. These are the next arithmetic queue. The lookahead does not say that the full lattice has no intervening mirror.

The previous strict spacelike witness survives every new wall, independently verified with exact Fraction arithmetic. Thus this finite prefix has infinite volume, while the full host remains open. The frontier and the finite transport check below do not imply a global chamber stabilizer.

The exact [oriented frame-transport census](fresh_h9_a9_next_20261013_transports.json) checks all eight A4+A4+A1 diagram identifications for each of the 45 unordered corner pairs including self-pairs: **360 maps**. Of these, 346 have an explicitly saved nonintegral matrix entry. The 14 integral maps are nine identities and five occurrences of the known reversal; both integral directions and the form-preserving identity are checked. This includes each new corner's own diagram automorphisms. It finds no new integral chamber symmetry among these frames. It is not a census of all lattice isometries.

## Files and replay

- [Complete profiles](fresh_h9_a9_next_20261013_profiles.json), [relative orbit intervals and selected roots](fresh_h9_a9_next_20261013_orbits.json), and [independent first-mirror audit](fresh_h9_a9_next_20261013_independent.json).
- Certified PARI output for the [end directions](fresh_h9_a9_next_20261013_end_raw.txt) and [inner directions](fresh_h9_a9_next_20261013_inner_raw.txt), produced by the corresponding [end](../scripts/fresh_h9_a9_next_20261013_end.gp) and [inner](../scripts/fresh_h9_a9_next_20261013_inner.gp) scripts.
- [19-wall geometry and eight actual vertices](fresh_h9_a9_next_20261013_geometry.json), [exact edge queue](fresh_h9_a9_next_20261013_frontier.json), and [independent geometry/queue audit](fresh_h9_a9_next_20261013_geometry_independent.json).
- [Finite frame-transport census](fresh_h9_a9_next_20261013_transports.json) and its [exact producer/replay](../scripts/fresh_h9_a9_next_20261013_transports.py).
- Producers: [profiles](../scripts/fresh_h9_a9_next_20261013_profiles.py), [orbits](../scripts/fresh_h9_a9_next_20261013_orbits.py), [geometry](../scripts/fresh_h9_a9_next_20261013_geometry.py), and [frontier](../scripts/fresh_h9_a9_next_20261013_frontier.py).

```sh
PYTHONPATH=src:. python scripts/fresh_h9_a9_next_20261013_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_a9_next_20261013_geometry_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_a9_next_20261013_frontier.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_a9_next_20261013_transports.py --verify
```

The first checker independently enumerates finite profiles in a different order, solves ten covector equations for their fixed component, and replays every relative-orbit point with separate Fraction field arithmetic. The geometry checker separately verifies the literal roots, all finite-angle gates, eight vertex certificates, the surviving strict escape, and every frontier row.
