# Released-R norm239 checker: full root source review

22 October 2026. This is a complete SOURCE and mathematical-interface
review, not execution. I read the entire 670-line worker and 158-line
specification, reread the accepted finite-domain argument and inspected
the complete diff against the preceding successful worker. No import,
plan, test, syntax compilation or scientific run of this new worker was
used in this review. A separate author's complete review is pending.

## Bound sources

| Source | SHA-256 |
| --- | --- |
| `scripts/check_compact_h9_prime239_e7_a2_released_r_batches_20261022.py` | `c2dc88dc5f355e631317c8c1f829c2ec2413e1f742a9ff1aebaab0b8f1e7d380` |
| `results/compact_h9_prime239_e7_a2_released_r_checker_spec_20261022.md` | `040a4ebc56bda02c7dd43399ee03dc8231abf256bfa42f9c8acb6e0bef6c9da0` |
| `results/compact_h9_prime239_e7_a2_released_r_small_folded_domain_candidate_20261022.md` | `dd5b74807f6604cbcbaeb71e40d326691da7c7be82989418607e25ac464bdc38` |
| Accepted third-edge geometry | `a545ec06561d49e52d5df8483beabce09b76f87caffcb18eeb9255f5eaa85a68` |
| Previous worker, preserved unchanged | `6569c6083ce781400fa2052f5a863af526605e4d27e367b4f172f45a6d28c9e2` |

The accepted mathematical contract and its previous root/outside reviews
provide the original module, actual starting E7+A2 germ, primitive norm
alphabet, finite spherical-extension list and local-finiteness argument.
These are explicit inputs, not propositions newly established by this
implementation. A header recording that a source was initially prepared
or pending is historical; subsequent acceptance records supply its status.

## Exact field arithmetic and geometric reconstruction

The field, serialization, square-root, floor/ceiling and linear-algebra
primitives are copied literally from the pinned old source. I checked
their entire use in the new code. The sign algorithm compares rational
squares in Q(sqrt5), including opposite signs and zero. The field-square
routine includes rational roots, pure sqrt5 roots, both signs and zero;
its general case tests both rational norm signs and rechecks the square.
There is no floating-point ordering or imported scientific worker.

`geometry()` reconstructs the old boost from original norm2 reflections,
then constructs the accepted S and P*. It checks R²=S²=2, R.S=-1,
P*.R=P*.S=0 and the correct negative future ray. Reflection K in R+S
has norm2 and sends R to -S, S to -R, fixing P*. It is used as an
integral isometry; no chamber symmetry is inferred.

The e7 check is essential. The original e7 vector has H projection H/2,
so the correction `(V1^j H-H)/2`, together with X and Y images, checks
the actual original basis rather than just the active H basis. Both
signed powers through three are checked quarter-integral and powers
plus/minus three integral. The conjugate V1=KVK fixes S and translates
in the outgoing direction because K sends S+R/2 to -(R+S/2).

The future point V1³P* is decomposed exactly as lambda P*+eta W1,
with lambda>0, eta<0. The upper is -eta/lambda, positive and strictly
before the ideal endpoint. V1³R is an original-integral norm2 reflector
incident there, with its sign and collision checked. This proves a
finite mirror upper, not that this particular mirror is the next facet.
All eight old native walls and r0,m,R,S occur in the known-wall checks.

## Completeness of the new scalar domain

The five profiles are exactly those in the accepted contract. E7 is
fixed by every phase, and Z0=J-3H/2 is in that fixed summand. Accordingly
the active three-coordinate representation omits only epsilon Z0; its
positive norm 3 epsilon²/2 is restored explicitly in the norm equation.

For g=seed.R and seed.S=-kappa, inversion of the A2 Gram gives
the coefficients (2g-kappa)/3 and (g-2kappa)/3. Since P* is negative,
the implemented T² is the positive-spatial-square-minus-2d expression
divided by -P*². The worker retains both exact field-square signs.
It imposes quarter membership on the ENTIRE reconstructed original
seed. It does not assume quarter membership of T or of either A2
coefficient, which would incorrectly discard valid seeds.

For norm2 the original contact lies in quarter O, |g|<46 and
|g'|<=2. The code generates all c=-87,...,87 and all integer a in
the exact conjugate interval, followed by the exact physical and
conjugate filters. For norm2b it uses g/b in quarter O, the proved
covariance divisibility, c=-14,...,14, |g/b|<7 and conjugate absolute
value<43/50. The rational interval generation is an overinterval;
strict inside tests follow. The row caps 2975 and 203 therefore apply
before root solving, yielding at most 48412 scalar rows over all units.
There are at most two field-square roots per row. This is termination
and coverage, not a measured runtime guarantee.

The centered-window denominator is g-kappa/2=W1.seed. Zero gives
no isolated collision and is skipped. Both signed nonzero centered
parameters, including boundary equality, remain. Reconstruction with
V1^j then checks original integrality, exact norm and, for norm2b,
every reflection-matrix coefficient. For norm2, integral coordinates
in the original even integral module imply integral reflection.

Finally the worker checks the outgoing collision signs, finite upper
and acuteness to EVERY known actual wall, including the released R.
It imposes no acuteness test on the folded seed. Thus all actual first
simple facets are covered by some unit, and every admitted candidate
is a genuine mirror of the original full arrangement. Combined with
the starting full germ, these two implications identify the first
collision after all units finish. An empty single unit has no such
conclusion.

## Exact minima, tie cone and execution boundary

Each unit keeps its exact least collision and all distinct original
root ties. The tie bound raises an error rather than truncating. The
combiner requires all twenty distinct NEW phase/profile kinds with
identical current bindings; it rejects old worker outputs. It decodes
each root and recovers its E7 covectors, original decomposition,
signed inverse phase, quarter seed, scalar profile, centered window,
original reflection and collision. This revalidation supplements the
execution records; it does not authenticate them by itself.

For a proposed endpoint, the code checks negativity and every known
wall. For each tied root it tests positive definiteness of the actual
nine-wall Gram E7+S+root and exact incidence. Every other tie must be
an exact nonnegative linear combination of that basis, with all ten
original coordinates reconstructed. A full germ is certified only
when exactly one candidate has this complete cone property. No row
order, predetermined root or presumed diagram supplies a winner.

The eight fixed inputs and self are streamed-hashed before and after;
the combiner also freezes every part's bytes through completion. Output
destinations must be new direct results JSON files. Storage is bounded
by small scalar state and the explicitly capped tie dictionary; there
are no saved scalar arrays or loaded search catalogues.

Root must separately match terminal actual bounded-runner ledgers and
their stdout/output hashes for all twenty units and the combiner.
`PASS_COMPLETE_DOMAIN` in an untrusted or incomplete part collection
does not replace those runtime records. No execution, new wall, new
germ, chamber period, finite volume or polytope example is credited by
this source review.

## Verdict

The frozen worker correctly implements the accepted complete finite
domain at the declared mathematical interfaces. I found no source-level
blocker. It is eligible for guarded admission ONLY after the separate
author's full review also passes. Every launch must retain the unchanged
shared lock, fresh 40GiB disk/6GiB available-RAM admission, 8GiB/3GiB
abort floors, 59-second/1GiB/one-thread watchdog and new-only outputs.
This review itself grants no resource admission or actual runtime PASS.

## Separate review received and admission clearance

I subsequently read the entire separate-author review
`results/review_prime239_released_r_worker_and_contract_full_scope_foundations_20261022.md`,
SHA-256 `430f4870af75869bf12d52e3d6aec4923a803fc447d380d92c2bc720486d284c`.
It binds these exact unchanged worker/spec bytes, independently checks
the complete new scalar/geometry/implementation boundary and discloses
its inherited first-edge authorship. It requests no correction and
credits zero execution. Its literal source comparison counts twenty-nine
unchanged generic blocks, supplementing its full review of their new uses.

With that full independent review received, the frozen implementation is
cleared for root-only guarded execution, one phase/profile at a time.
Fresh resource admission and actual terminal runner evidence remain
mandatory for EACH unit and the separate all-twenty combiner. No wall,
germ, example or runtime conclusion is credited by this clearance.
