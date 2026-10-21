# Root review of the unified sharp right-angled manuscript

20 October2026. Complete written manuscript review PASS; this record
does not claim typesetting, figure layout inspection or Lean execution.

Source read in full: `manuscripts/right_angled_sharp_dimension_bounds.tex`,
SHA256 `86e77ecdaf85d2d64560ad44a0226b9f5fddfa6b70f06617f93baa64d03170fa`,
375 lines. The general theorem has the same full proof as the reviewed
note9470fa1f...babdbe. The ideal theorem is correctly shortened using
the shared finite-point projection, rather than the old cusp-height
propagation proof. Both earlier manuscripts remain valid predecessors.

## Shared projection and dimension bookkeeping

The paper fixes outward UNIT normals and the future unit Lorentz
model. Every starting point is a finite actual interior point. At
each step the maximum is taken over all actual intrinsic walls;
the all-wall inequality proves strict feasibility at the projected
point against every other wall. Orthogonality of a ridge neighbor
preserves the intrinsic unit normalization. These checks justify
iteration inside the selected actual facet.

The height update A->A/sqrt(1+A²) gives exactly
A_j<=A0/sqrt(1+j A0²)<1/sqrt j. A finite positive A0 makes the
inequality strict, including tied maximizing facets. Three steps
start in dimension9 and end in6, whose intrinsic Lorentz space has
dimension7 and whose radial sphere is S5 in R6. No step descends
to dimension2. The one-step ideal argument starts in5 and ends in4,
with radial sphere S3 in R4.

## The universal eight-dimensional ceiling

The tensor-square proof is complete in the manuscript. I independently
checked the36-entry rank-two and216-entry rank-three contractions,
S(u).S(v)=24(u.v), the cubic factorization and exact positive square
coefficients. They prove the real sum inequality directly, without
an imported spherical linear-programming theorem.

Strict products below1/4 give at most27 directions first. At27,
every off-diagonal cubic value must vanish and every product must
be-1/2. Four such unit vectors have squared sum-2, giving the strict
cap26. The numerical boundary is proved, not rounded away.

The numerical geometric input is explicitly the universal finite-volume
RA6 facet minimum27, stated in the primary mathematical text of
Alexandrov2023 Section4 and attributed there to Dufour Lemma1 and
Proposition4. It is not confused with the mere existence of a27-facet
polytope. Potyagailo--Vinberg2005 is cited separately for the dimension8
examples establishing sharpness. The proof does not use the old
mixed10 theorem, the ideal4 theorem, arithmeticity or an enumeration.

## The entire-ideal four-dimensional ceiling

One projection makes every intrinsic height of the actual ideal4
facet strictly less than1. The radial code bound is consequently
strictly below1/2, so Musin's non-strict24 cardinality bound applies.
Kolpakov's24-facet minimum gives equality and his METRIC theorem
identifies the ideal24-cell. This is not a kissing-code uniqueness
claim. The explicit opposite normals force any point with all heights
at most1 to be its center with all heights1, contradicting strictness.

Actual five-face or nine-face inheritance supplies the higher tails.
The published regular ideal24-cell supplies sharpness in dimension4.
The three numerical external statements and standard acute-face
foundations are isolated at the beginning of the paper.

## Presentation and remaining checks

The two results share one focused paper. The ladder figure correctly
represents successive actual facet choices and projection formulas,
not ordinary edges or an unproved geometric realization. Its final
code has one vector per actual wall. Final raster layout is to be
reviewed after a capped build. The old ideal-cusp graph proof and
12-page mixed10 proof need no expansion and are not additional
current proposed publications.

This is a mathematical proof review with named published inputs;
external refereeing and worldwide priority are not asserted. Exact
SOS/model checks and the conditional Lean logic require their own
actual execution records, with all real and geometric interfaces
kept explicit.
