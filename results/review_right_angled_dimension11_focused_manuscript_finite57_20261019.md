# Independent full source review of the focused dimension-eleven paper

19 October 2026. **PASS for mathematical translation and completeness.**
This review binds `manuscripts/right_angled_dimension11.tex`, SHA256
`6c455e9e36420581a1987d4fab8e24e2bfa745c8c4648b4c81ff73f89576ddcf`
(759 lines). I read the entire source, including all three vector figures,
against the new positive-part proof and the inherited opposition,
weighted-row and polygon/three-face arguments. This was a written/source
review only: no build, script replay, Lean job, PDF inspection or scientific
calculation was run in this lane. The proof author and manuscript author
are different agents; this is a translation review, not an independent
discovery of the new analytic argument or external refereeing.

## Geometric and analytic coverage

* Sections 1--2 identify the published inputs separately from the new
  proof. The projection argument proves the normal-product gap with an
  actual relative-interior point. The bad-triangle wall, all-wall
  feasibility of its three-wall opposition, convex-hull closure and the
  two-finite-vertex multiplicity are given, rather than assuming a formal
  Gram diagram is a chamber face.
* The invariant `h`, its spacelike square, and one fixed future null `Z`
  per actual finite opposition component are established before ratios
  are compared. The energy signs and the equations `sum p=3`,
  `sum p^2=1` are consistent. Individual coordinates need not be positive.
  Isolated components cause no exception.
* Section 3 distinguishes all 220 formal triples from the actual subset.
  Their formal vectors are future unit timelike, the ratio interval is
  `[1/2,2]`, and the inclusion probabilities and centered moment identities
  have the correct normalization and signs at every finite vertex.
* Section 4 contains the complete polynomial obligation. Its positive
  and negative halves, reciprocal remainder identity, degree-six
  expansion, Bernstein definition, three coefficient rows, central
  convexity argument and endpoint values cover the entire interval.
  The negative cubic coefficient uses an upper bound on `-E d^3`;
  positive quadratic and quartic coefficients use upper bounds in the
  ordinary direction. No sampled positivity or unrecorded root location
  is substituted for the interval proof.
* Section 5 inserts absent triples only after proving `P>=0`. It sums
  over ordinary finite vertices. Each actual undirected edge has the
  same componentwise `Z` at both ends and contributes
  `r^-2+r^2-2>=0`. Thus the passage to
  `2B <= (1420561/12500) f` uses neither a maximum-height estimate nor a
  weighted-to-unweighted averaging assumption.

## Counting coverage

Section 6 supplies the truncation convention, both relevant
Dehn--Sommerville rows and the strict weighted row. The cubical-link
numbers 2048, 11264 and 28160 occur in their correct roles. The
eight-face bound of 82 and the elementary ten-face minimum are explicit
external/counting inputs; the coface multiplicities five and three
produce the stated strict signs.

Section 7 proves the required small-face statements rather than invoking
a census: `F>=6`, the unique `P3` equality at `F=6` or `nu=9`, the
two-bad-triangle attachment, the zero-excess `Q8` structure and its
one-ideal-ending-edge property. The attachment normalizations and the
cross products `-t,-1/t` are consistent; the final null sum is made from
actual walls and convexity closes the six-wall body. The exceptional
incidence bound uses actual ideal-ending directions, and its displayed
13-value integer row covers every possible value of `k`.

Section 8 uses precisely the same `B` and finite-vertex count `f` as the
transport bound. It distinguishes nonnegative excess `E` from signed
defect `R`, accounts for removal of each `P3` contribution, derives the
Euler/coface identity, and retains strictness in both lower bounds. The
two positive multipliers cancel the ambient cusp count, giving
`2B > (184460/1619) f`; the final positive cross difference is the same
5861741 as the exact coefficient packet. There is no global facet-type,
finite-edge or arithmetic restriction hidden in this elimination.

## Figures, references and verification boundary

Figure 1 is explicitly an incidence drawing and marks opposition as a
non-edge. Figure 2 shows certificate intervals rather than numerical
samples. Figure 3 separates the formal row from the actual graph and
shows both reciprocal directions with the correct exponents. No figure
is used as a missing geometric proof; this review makes no claim about
the rendered PDF.

The cited Dufour ceiling twelve supplies the higher-dimensional tail;
the compact and entirely ideal published ceilings four and six are used
only for the final unrestricted corollary. The dimension-eight example
remark makes no sharpness claim. The published eight-dimensional facet
minimum and standard acute-wall/link/Dehn--Sommerville inputs remain
explicit external dependencies.

The verification paragraph accurately separates the already executed
exact coefficient replay from written universal geometry. I also read
`docs/right_angled_dimension11_proof_map_20261019.md` and
`formal/CoxeterProofs/RightAngledMixed11.lean` (SHA256
`f0b1cffede6210dbddf7cc7d95f29d14baa5984bcd2e19641e029afb1e7f864d`).
The conditional Lean interface exposes the count identities,
exceptional-incidence inequalities and aggregate analytic upper bound;
it does not disguise real geometry as integer ratios or claim to
formalize the polynomial interval proof. Its stated axiom/audit boundary
agrees with the recorded root verification. No source correction is
requested.
