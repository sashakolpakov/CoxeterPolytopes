# Certified first mirrors, face cycles, and a remote escape

The edge-walking strategy now gives a stronger arithmetic result than the earlier bounded unit search: complete finite-interval root censuses certify the first two proposed endpoints and find a **previously missed norm-(2+T) mirror before the third proposed endpoint**. Separately, coherent geometric additions produce a 79-facet H9 chamber with sixteen certified finite vertices and five compact two-dimensional faces. All completed chambers still have exact infinite-volume escapes. **No new compact or finite-volume Coxeter dimension has been established.**

This round follows the user's suggestion to distinguish wall enumeration from edge traversal. The previous turn made progress by auditing those algorithms and introducing the explicit frontier; this round implements complete arithmetic first-endpoint checks and exact geometric edge/face closure checks.

## Complete arithmetic edge intervals

The arithmetic host remains over K=Q(T), where T=2cos(pi/8) and T^4−4T²+2=0. For a known corner and a proposed finite endpoint, the new calculation bounds every root whose mirror meets the segment, at all four real embeddings. It enumerates the complete necessary norm pool, integral edge-pairing profiles, and bounded time coordinates, solving the remaining quadratic exactly. Nonsquare residual classes and other binary norm-equation orbits are included.

The first census tests 23 profiles and 2,229 quadratics. It confirms the previously proposed norm-two mirror is first. The second tests eight profiles and 1,206 quadratics and confirms its proposed mirror too. Exact finite F4 block censuses lift the checks from the forced H5 face to the full H9 edge: the edge/start/end stabilizer mirror counts are 30/33/41 and 34/41/49 respectively.

The third census tests 7,105 quadratics and finds a norm-(2+T) mirror at segment parameter T/2<1, before the previous norm-one proposal. Its height is about 120,836, compared with about 3,174,864 for the earlier proposal. That residual class was absent from the earlier square-residual unit grids. Four symmetry images fit all 182 preceding walls and give a corrected 186-wall H9 / 76-wall H5 cone. Independent Fraction replay checks all 17,205 H9 and 2,850 H5 pairs, actual vertices, and surviving strict escapes. This is a concrete correction enabled by exhaustive first-mirror enumeration. See [the arithmetic first-endpoint report](fresh_h9_quartic_first_endpoint_20261013.md) for the exact bounds, corrected chamber, independent replays and limitations.

These are complete local edge steps. They do not exhaust the whole lattice chamber's edge queue or prove reflectivity of the host.

## Coherent geometric face closures

Starting from the previous 74-wall split-seed chamber, five mutually compatible new caps retain every old wall. The exact resulting **79-facet chamber** has sixteen certified finite vertices in its two reached components. Its finite graph has nineteen edges and five cycles: two triangular and three quadrilateral actual compact faces. One new mutual angle is pi/10.

The endpoint census is exact. Every candidate finite vertex has a positive-definite incident Gram matrix and an actual point certificate or a center-projection proof. For every reached spherical eight-wall edge, every remaining wall is tested as a possible finite or ideal endpoint. After dotted-pair exclusions, 216 ninth-wall tests certify the final frontier. It still has **106 open ends**, and twelve strict spacelike witnesses prove infinite volume. See [the second-frontier report](fresh_h9_split_second_frontier_20261013.md).

A separate [two-edge cycle construction](fresh_h9_split_cycle_constraints_20261013.md) gives an exact 73-wall chamber closing a quadrilateral face. Its continuous continuation permits every real lambda>2 while keeping all off-edge wall relations strictly dotted. This exposes a family omitted by an extra-finite-anchor search. However, an all-parameter positive-cone identity makes every member weaker than the lambda=1 cycle cap, and an exact escape survives the whole family.

## Why the larger first-generation bank is exhausted

The existing 3,857 no-removal proposals include 142 caps closing two old directions and 3,715 closing one. Numerical compatible selections close all 32 original open directions, but expose new ones. In a representative 100-wall selection, the reached graph has 36 finite vertices and 248 open directions.

The decisive result is independent of those numerical selections: a rational interval reconstruction of **every one of the 3,857 proposed roots** certifies a single common future strict spacelike point satisfying all 72 old and all 3,857 new inequalities. No subset of that bank can give finite volume with the old walls retained. The norm quadratics, unique branch separation, strict center and escape signs are all certified. The saved 70-digit replay and a separate 90-digit execution pass; all 38,570 new coordinate enclosures nest. See [the bank certificate and numerical scope](fresh_h9_split_multiple_edge_caps_20261013.md).

That exclusion is restricted to the defined bank. It says nothing about arbitrary new caps, other labels or branches, changed old walls, or later edge generations. The opposite second-frontier scouts also miss this remote gap numerically, so repeatedly enlarging those same local cap families is not a justified route to a finite-volume claim.

The continuation is to use the certified first-mirror routine on the arithmetic frontier and to target the remote gap with genuinely different geometric constraints or changed retained walls. The saved escapes and exact open-edge certificates remain the tests that every proposed compact or finite-volume construction must overcome.
