# Independent source audit of the frozen compact17 foundation translation

24 October 2026. Reviewer: `/root/audit_sharp_ra`. **Source review only.**

The frozen translation should not be accepted verbatim. It has four local
statement/wording defects, described below. None of the actual applications
examined requires the invalid form of those statements. The new A12, A15,
N17–19 and exact18 arguments have a coherent original-vector dependency
route, and their displayed payments fit their individual anchors. This is
scoped mathematical review of the written foundation batch, conditional on
the stated classification/geometry inputs. It is not acceptance of a complete
compact17 manuscript, a runtime receipt, or a formalization.

## 1. Exact inputs and amount read

The immutable input directory is
`results/source_archive/compact17_foundation_transcription_20261024/`.

| Input | Bytes | SHA-256 |
|---|---:|---|
| `compact_farthest_vertex_bounds_consolidated.tex` | 193119 | `6ee917c01903a4843ada6e43b9f37ecb1930fafec84f71d7eb84c69e38e761b2` |
| `source_ledger.md` | 38486 | `3ecfc31c0c5de4fddf5cbe7bb6fe745d556c7ea078677362c25412420c82b082` |
| `dependency_snapshot.json` | 1122137 | `dcc7b00debacc5f585a7eb67fbce42dc6a6f011fe2ccfe62ea66e8a86cba3b38` |
| `results/compact17_foundation_transcription_freeze_20261024.json` | manifest | `55803c059eb4ca9e54533bc51bef15153f22a827da90d7a9d86f3900ef5b0fb7` |

The hashes above were checked by ordinary file hashing. They identify bytes;
they do not certify the mathematics.

I read all 5093 lines of the frozen TeX and all 209 lines of its ledger.
For the large JSON snapshot I read the opening source/status metadata,
the final470 input boundaries, the selected section DAG, and the closing
foundation/freeze records; I did not reread its approximately 20000-line
archival literal-pin inventory or re-audit every historical source to which
it points. The ledger is the provenance map, not a substitute for the new
proof bodies. The contradictory-looking older queue entries in that
snapshot are explicitly superseded by its final freeze record.

I also consulted the required repository README, manuscript index and
current-bounds guide. I did not write or edit the frozen compact TeX or its
ledger. Earlier reports authored by this reviewer concern the right-angled
and finite57 lanes; their conclusions do not establish this compact batch.

No worker was imported; no enumeration, test, scientific process, Lean job,
TeX build or Sphinx build was launched. No canonical TeX/PDF pair, source
archive, shared documentation or bound was changed. This report is the only
new deliverable of this review.

## 2. Required local repairs

### F1. The point-pressure lemma needs the spatial hypothesis

Frozen lines182–199 write `xi = A p - U` but do not require `U` to lie in
`p`'s Euclidean orthogonal complement. That hypothesis is essential both to
the proof and to `U^2=A^2-k` immediately afterwards. Literally the lemma is
false: take the feasible actual vertex `xi=p`, write `A=2,U=p`, and use
`<x,p>=0`. The asserted conclusion would be `0>=1`.

Correct statement: require `U in p^perp`, `xi=A p-U` future timelike with
square `-k`, and actual feasibility. Then `A=-<p,xi>>0`, and the all-point
comparison gives exactly the displayed pressure. At a different actual
vertex use that vertex's spatial complement and its separately stated
comparison factor.

Every actual application read in this batch has the required spatial form:
old normals, original component duals and their Euclidean linear
combinations lie in the relevant complement. The local-A18 argument uses
the complement of its actual vertex `q` and a different vector `y`; it
does not add `y`'s energy to the original `x` budget.

### F2. Reference payment needs a nonnegative new weight, and generally gives a weak inequality

Frozen lines330–350 require the retained old coefficients
`b_i-mu*t_i>=0`, but also need `mu>=0`. Otherwise the new lower pressure
cannot be multiplied by `mu` in the indicated direction. For multiple new
pressures all new weights must be nonnegative.

The generally valid first inequality in (payment) is `>=`. It is strict
if some strictly constrained original anchor retains a positive coefficient,
or a positively weighted new pressure is itself strict. If all old
coefficients cancel and the new pressure is only weak, the claimed `>`
does not follow. For example `w=P`, `mu=1`, `Q=0-P` and the new pressure
`<x,0>>=0` produce `W=0` and equality on both sides.

All applications in the reviewed new bodies use positive weights and have
strictness available:

- A12/A15/A16 retain positive margins at the individually debited strong
  anchors; their conservative deficit estimates are also strict.
- N18/N19's changed endpoint references use strict actual endpoint
  pressures and maintain positive coordinate junctions.
- Exact18's `1/7` core consumers use K8–10, with terminal price strictly
  greater than `1/7`.
- The final D18 payments retain the positive unused terminal16 coefficient;
  the local-A18 graft leaves a strictly positive price at its actual strong
  anchor.

Thus this is a repair to the general payment statement, not a discovered
failure of one of those applications.

### F3. Two actual-deficit monotonicity assertions are too strong

At lines3604–3612 (N17) and4441–4444 (A12), the derivative computation
establishes that `A_far` minus the linear reference pairing increases.
It does not establish monotonicity after also subtracting the varying
core radical. In the dotted case `sqrt(1-gamma^-2)` has infinite right
derivative at1, so the actual deficit initially decreases.

The written numerical bounds remain available by the method the displayed
data actually support: subtract a fixed global radical cap from the
increasing time-minus-linear envelope. Use `1/4` for finite owns at least7,
and `.338` for dotted owns in `(1,17/16)`. N17 then gives
`2.391-.25-5H6/6=983/3000>.327` and a dotted bound exceeding `.3343`.
A12 similarly gives `2.266-.25-1.60125=.41475`.

Replace “the deficits increase” with this conservative-envelope argument.
The A15 body already uses the narrower, valid wording.

### F4. The local `1/7` payment is overbroad at K13

Lines4786–4791 correctly distinguish terminal capacity `>1/7` at K8–12
from the weaker `>7/50` at K13. Lines4805–4813 then offer weight `1/7`
without retaining that restriction. An ordinary reference with the
preceding knot retained has terminal price
`sqrt(13)-sqrt(12)<1/7`. Consequently the advertised arbitrary-K13
payment is not funded from the stated information.

A safe correction is to restrict that displayed `1/7` formula to K8–12,
and state the general formula with `0<=mu<=beta_K` otherwise:

`e > q_K + 2*mu*(G-f_K) - mu^2*(R-K)`.

Alternatively K13 requires an additional price proof, for example a
specified preceding omission; the existing `>7/50` statement alone is
insufficient.

The actual displayed consumers are not affected. D16's two core stages
use only K8/9 and K9/10. B17's K13 extension pays a new slope less than
`.39/3=.13<7/50`; it does not use the `1/7` local-core formula at K13.
The K12 extension has new slope less than `.53/4<1/7`. These are distinct
individual-price checks, not an inference from the total reference energy.

## 3. Geometry and the dependency order

The foundational route is acyclic at the following exact interfaces:

1. Compactness, acute wall products, spherical/affine classification and
   Lanner's order-five bound supply the external geometry/classification
   boundary. Strict interiority and one actual farthest vertex supply the
   all-point comparison.
2. Positive inverse matrices give actual positive projections. Reciprocal
   shadows give full-tangent comparisons and prohibit an opposite wall
   from closing orthogonal original directions.
3. Ordinary small components and actual I/J/N7/O7 defects give the
   spectator-uniform debits. The rank-free I argument and typed propagation
   use only that actual geometry, not a numerical dimension cover.
4. The A-chain branch uses the common-pair inventory and positive L/H/F/C
   feet. The B/D branch uses the actual large-component floor and the
   exact-rank prefix inventory. These are different finite-type analyses.
5. The new N13/N16/N17–19 windows rederive their actual-rank inventories,
   nonadjacency and reference data. Their paid helpers use the common
   small actual core matrices, not a larger-rank energy theorem.
6. Exact18's fresh B3/D4/D3 boundary layers extend local coverage only after
   their additional actual roles have been handled. Forest composition
   comes afterwards.

The basic positive projection adds nonnegative multiples of actual wall
normals. Its selected products vanish, its unselected products stay
negative, and its norm is negative. For a positive actual-normal core,
the selected products must instead be checked directly; all remaining
wall signs then follow from acuteness. This includes unknown walls,
foreign walls, both genuine D forks, and a wall released at an adjacent
vertex. A spacelike core is never itself normalized: the finite-point
construction uses `p+(A/delta)*eta`.

The reciprocal-shadow derivation uses the inverse of the retained old
Gram matrix, so its nonnegative Schur remainder is legitimate. The
orthogonal-opposite obstruction uses strict interiority of the same
original `z` and both original half-height inequalities. Its foreign-cost
consequence is strict `q<A^2`, not an assumed lattice bound. The native
gap `(-1/2,0)` concerns products of actual unit wall normals.

The local-A18 main proof and Appendix localAinputs point to one another in
exposition. The actual dependency is not circular: its `3/2` ratio
calibration is the explicit positive polynomial before the endpoint
assumption; the appendix uses that primitive, and the `.78` endpoint
argument uses the resulting inventory. Move that calibration before the
inventory statement in a revised exposition to make this order obvious.
The appendix does not use the original B18 consumer or the new `.78`
conclusion to establish its own inventory.

## 4. Ordinary components, defects and uniform B/D floor

I read the complete ordinary appendix1510–2262, including the strengthened
B3/H3 costs and both propagation arguments. The common-wall forcing
inequality keeps the full foreign cost in `h=q+Q-1`; its small-height
isolation follows from each positive foreign cost being at least1/4 and
strictly below `h`. The actual finite endpoint, not the full cap diagram,
is the spherical object being classified.

The short A1 split keeps the isolated dotted I and the two-A1 J distinct.
J pays both of its factors; its light factor is not charged again. The
two typed rank-two exceptions are isolated in their whole original blocks.
Their products with another defect or short wall lie in the forbidden
native interval. This supplies uniqueness/compatibility independently of
the later dimension-cover enumeration.

The rank-free I payments use the actual other opposites and the original
orthogonal projection. The scalar `.767` barrier is reproduced as a real
quartic argument; the mixed-title historical source contributes no field
condition. The typed-seven propagation checks the full native families
and all foreign costs, including the old-I2(7) common-wall alternative.

In2263–2470, the B/D floor first establishes distinct4/5 and6/7, uses
those constraints to eliminate common1/2 and2/3, and only then combines
the three distinct-pair sums. It does not use its final energy contradiction
to justify an earlier inventory. The B9 short-end successor retains the
removed own6 check; D9 retains its right-fork exception. Positive actual
points in the `4@2+3@3` and isolated7 cases have explicit positive old-root
coefficients, selected products, negative norm and positive time.

The many small finite inverse/contact tables here are mathematical finite
algebra, not newly replayed script outputs. This review checked the stated
matrix/Schur mechanism and the branch coverage; it did not run an algebra
verifier or create a new execution receipt for every scalar entry.

## 5. Shared A-chain bodies and endpoint supports

The actual common1/2 and3/4 inventory is stated before the energy-dependent
H/F/central-double discards. Consequently the later4/7 and5/7 proofs can
reuse the classification and redo the correct individual cutoff, rather
than inherit an exclusion proved only at energy1/2.

The H/F/L constructions retain every foreign coefficient. L clears the
entire spectator forest positively, gives a pure A-spatial vector and
uses the explicit `8/81<z0<1/3` range. At rank7 the N8 data are recomputed
rather than inferred from N9; in particular the two central-common and
far-end alternatives are separately handled.

The paired-end Gram is entrywise decreasing in N, permitting the stated
rank-uniform seals. Raw support length matters: L has length3, H length4,
F length5, and C length5. In the4/7 proof the only overlapping special
case is F/F at N9; it is handled by the exact actual-wall product
`1-phi^2/2+phi^2/N` in the forbidden native interval, not by a false
disjoint-support identity. For the5/7 proof N>=14 makes all length-at-most5
end supports disjoint. Reversal changes the raw sum's sign and gives the
positive centered cross `2BB'/N`, which is retained.

The central4/3 case, no longer eliminated by the new individual cutoff6,
is explicitly restored as foot C. Its positive actual correction and
selected products are written out. No inherited half-bound discard is
silently reused at the larger energy.

## 6. N17–19 unification and individual prices

The central window has retained side order at least5. The finite-white
classification is complete at the spherical-classification boundary:
one permitted attachment per retained path, a joined ordinary-path bridge,
and all foreign attachments represented by their Schur cost. Two contacts
on one retained path create a cycle. Joined ranks exceed8; exceptional
H/F/E joined completions cannot be inserted at these ranks.

Zero own is treated separately. Isolated owns cover `u<1`, `u=1`, foreign
`u>1`, the native5 H4 core, native7 A2 core, finite>=8 singleton and
affine/parallel cases. The direct ratio comparisons retain foreign Q.
Every residual weak pure helper has Q0 and `gamma<17/16`, so an unbounded
dotted input has been cut down by an actual geometric inequality, not by
testing finitely many labels.

For adjacent pure helpers, the mutual native alphabet is determined before
the energy target. Same-end/outward/reversed cases are separated. The
reversed determinant polynomial is increasing in each own and concave in
the position. Its N17 endpoint bounds increase with N at the specified
comparison pairs. Native5/7 uses a genuine degree-three obstruction over
Q(sqrt5); two native5 owns use the nonzero rational coefficients
`N(N+6)` and `-N^2` of the raw determinant. Those coefficients are signs
of literal cosine algebra, not a restriction on the polytope's field.

The surviving two-cap Schur matrices have positive inverses. For opposite
helpers the ratio of positive weights lies in `(5/6,6/5)` simply from
diagonals in `(1/2,3/5]` and off-diagonal `-1/2`. This gives
`D/L^2<=3/55`, `L>40` and the printed radical-loss bound. The spatial
rationalization uses adjacent dual products and the actual lower bounds
on both lengths; all original spectator blocks remain in the projection.

Near-fours are removed only after the appropriate target is assumed.
At the other end either a second near-four supplies a mate or helper
nonadjacency supplies a strong anchor. The correlations are bounded at
the actual N17/N18/N19 indices. The proof does not assume a near-four
is itself a pure helper in order to invoke pure-helper nonadjacency.

The common endpoint cores have the advertised all-wall signs:

- Native5: selected old roots are cleared; the spacelike residual is used
  only inside the actual future finite point.
- Finite>=7: the core products on the cap and own root vanish, and the
  helper product is `-1/2+2(1-gamma^2)<0`.
- Dotted: the cap product vanishes and the own product is
  `-gamma+gamma^-1<0`.

Their times are positive, all corrections have nonnegative coefficients,
and their spatial vectors are pure A-component vectors. Removing the
positive old-normal products only weakens their valid pressures.

I checked the reference omission counts directly. For N17 the internal
counts are3/3/2 after zero/one/two weak ends; for N18 they are4/3/3;
for N19 they are4/4/3. Nonadjacency makes the loss intervals disjoint.
The printed maximal second differences bound every eligible position.
For example the N17 two-end energy is
`H6^2/6+.083^2+.041^2-2*.042^2/4=1195276/1500000`.

In the N17 paid proof the separate prices are `839/6000` and
`349/2500`; a near foot costs at most31/320, and two feet cost at most
101/960 at either actual end. Those are smaller than the actual individual
capacities. The six/six residual formulas are fixed-block variances;
near orientation adds `gamma/3`, and six/five far variance decreases1/60.
The full two-foot cross is bounded by
`||YL+YR||^2<=2||YL||^2+2||YR||^2<32/5`.
No separately optimized energy is added.

At N18/N19 the raw endpoint pressures are fresh. Finite7 and8 use their
own radical/time brackets; finite>=9 uses a monotone envelope with the
fixed1/4 cap; the entire dotted interval uses its exact rationalized
time correction. The resulting `.442`, `.432` and near `.58` pressures
modify one reference by preserving its end-block sums. The written
coordinate junction bounds keep all required new coefficients positive.
For two far helpers at N19 the two changed end blocks are disjoint.

## 7. Fresh A15 and A12 bodies

### A15 at N16

The fresh N16 isolated finite>=8 singleton estimate uses its exact finite
projection and checks loss monotonicity at the actual minimum D=55/8.
It does not import the slightly inadequate larger-rank coarse estimate.
The polynomial derivative comparison is positive throughout the stated
real box. All foreign, native5/native7, affine and parallel branches remain.

The reversed-helper case has its fresh endpoint value
`P(.809,.923)=.0907947698>0`; the native5/native7 cubic obstruction and
two-native5 determinant handle the remaining cases. Mixed-point bounds
are then re-established with times>2 and the actual N16 dual minima.
Near-fours are removed using their actual mates and the opposite10/11
anchor. No N17 energy conclusion has been shifted down a rank.

The zero/one/two-end omission counts are3/2/2; for two weak ends the
largest eligible second difference is `.046`. The one-end separate
prices are43/300 and179/1250; two-end prices are43/300. The stated
31/320,1/100 and101/960 debits fit. The core variances are recomputed at
N16 and the full two-foot norm is retained. The final two-end total
`1242171/1500000=.828114` exceeds413/500. The rank-three spectator
consequence uses a disjoint proved ordinary debit>17/70, not a historical
consumer imported as an extra premise.

### A12 at N13

The loss is freshly `.65`; near-four is directly strong at its new ratio
threshold and is not silently omitted. Actual isolated pressures and
positive foreign-cost inequalities close all other roles. The weak helper
bound `gamma<.99` excludes dotted owns geometrically.

For reversed adjacent helpers, zero mutual would make P vanish. Two
non-five owns and two-five owns are separately excluded. With exactly
one native5, P is a nonzero quadratic over Q(phi), giving degree at most4
over Q for the other literal cosine. The complete totient list is
`m=7,8,9,10,12,15`. The reasons for excluding it are explicit:

- m7,9 retain degree3 after a quadratic extension.
- The m8 real cyclotomic quartic has unique quadratic subfield Q(sqrt2).
- The m12 biquadratic field has only quadratic subfields Q(sqrt2),
  Q(sqrt3), Q(sqrt6).
- For m10 the irreducible quadratic is even while P has nonzero linear term.
- For m15 the specified irreducible quadratic cannot be proportional to P,
  since the resulting rational-linear expression in phi cannot vanish.

This finite-label reasoning uses only a standard real-cyclotomic degree
formula. It does not assume arithmeticity, conjugate positivity, ground-field
degree or integrality for the reflection group.

The reference omission counts0/1/2 outer sites give exactly the stated
energies; with both outer sites weak,6 and7 are both strong. The fresh
six-block variances and separate prices are explicit. One selected helper
suffices even in the worst two-end reference. The finite-five far pairing
is bounded by1.46, preserving the correction to the older false1.438.
After F3's envelope wording repair, the uniform deficit>.4 and weight.1
give the stated total `3959881/6000000>13/20` without a second helper or
an overlapping energy debit.

## 8. Generic prefix and fresh exact18 continuation

The generic theorem is used only for `6<=M<=16`, B rank>=M+5 or D
rank>=M+6. Its deleted-node right ranks and all old long-root neighbors
are explicit. The isolated-own, complete left/right/bridge and zero-own
cases retain all foreign Schur costs; exceptional D-fork attachments are
kept with their allowed Q. The foreign far-helper estimates cover the
entire permitted index interval. They are not a bounded own-label scan.

The two-cap nonadjacency point clears actual old roots, has positive
inverse weights and appends untouched spectator blocks. All four weak5
pressures use actual selected roots and the correct finite-point
normalization. Their finite and dotted bounds cover the complete own
intervals. The strong5 and weak5 interpolation branches are both proved;
the latter is not supplemented by an overlapping separately optimized
first-coordinate energy.

The exact18 local continuation has three distinct classification changes:

- B16 at node13 retains A12+B3. Its extra short3 is the sealed F4
  white extension; the zero-own affine path is retained and positive owns
  are directly strong. This supplies the additional actual layer before
  interpolation is extended.
- D16 at node12 keeps all three D4 leaves, including both signed forks.
  The remaining native5/fork3 point uses its positive left H4 correction,
  pure scaled displacement flat through9 and legal K8/9 payment.
  At node13 the full D3=A3 list is middle3, end3, end4, end5. Both full
  bridges and the explicit end4 null are present. The remaining native5
  and finite>=7 end3 cases have actual H4 or two A2 cores, respectively;
  the two negative-core spatial pressures are averaged, not their energies.
  Their legal K9/10 payments retain both forks and every spectator wall.
- B17's B3 layer is separately supplied before the B2 layer. The last
  exceptional short3/own3 cap produces an actual A17 endpoint using the
  true old short root. Positive inversion checks the released old16 as
  well as all unknown walls. The new S16 floor is appended using the
  distinct K12/K13 terminal capacities described in F4. An actual I is
  not charged the nonshort-A1 debit.

Thus the actual written exact18 consumers do not invoke the generic
prefix theorem outside its rank hypotheses. Neither a sign choice at a
D fork nor an extra finite white role has been removed by arithmeticity.

## 9. Acceptance and execution boundaries

This review supports the new foundation architecture and the listed actual
applications **after F1–F4 are repaired in a reviewed successor**. The
frozen bytes remain untouched and are not approved verbatim. Root acceptance
and review of the exact repair delta are still required. The local repairs
must not be credited retroactively as if present in this source hash.

The external assumptions retained here are compact acute simplicity/wall
acuteness, spherical and affine Coxeter classification, the Lanner bound and
the elementary real-cyclotomic degree formula. I did not newly authenticate
the original Vinberg classification pages in this review. The published29
ceiling belongs to the separate high-dimensional tail, not to these local
component proofs. All displayed internal finite tables remain written
mathematical claims; no new exact checker run is asserted.

The frozen manuscript is still explicitly incomplete. Its own pending
queue includes the other exact18 row consumers and wholeA18; the full
higher-dimensional general covers/consumers; wholeB18 terminal proofs;
and the earlier wholeD18 head/.58/reference/U5/strong16/finite-head/fork/pair
bodies. Those are not converted into external axioms by being absent here.
The separate wholeA18 working fragment is outside this freeze.

I read the main final470 translation and checked its use of the repaired
spatial/payment primitives and its distinction between original x and
local y. This is not a new transitive audit of all of its historical raw
reference or terminal-inventory inputs. The global compact17 conclusion
continues to require those separately identified inputs and the entire
missing transcription queue. No full-manuscript mathematical PASS, Lean
completion, numerical replay, PDF build or publication claim follows from
this source assessment.
