# Bounded short-profile search at the actual thirty-one-wall edge

19 October 2026. **Source prepared; not executed.** No scientific import,
root search, finite census, Lean run or build was performed in preparing
this worker. The coordinator owns the serial execution queue and its
resource gates.

## Fixed source and mathematical input

The worker is
`scripts/search_h19_thirty_one_edge_short_profiles_20261019.py`, SHA256
`d92086f97ab207412c8910c795b3fe34f180bcdfd6b06b0ccb2462110da9453b`.
Its mathematical input is
`results/h19_thirty_one_escape_finite_profile_reduction_20261019.md`,
SHA256 `2be3d769f288ab25bfdb5831115c3002c1f8022890d30b4940d527ed4c01336c`.
The coordinator has independently reviewed that reduction. The worker
does not replace its finite-diagram completeness proof with a guessed
bond cutoff.

The current chamber is the pinned thirty-one-wall certificate
`results/h19_thirty_one_wall_seed_20261019.json`, SHA256
`5a3b25aa139793f304d535f56d6033ad103d88cffe07c20278bc7b42c2712d57`.
The worker reconstructs the actual native roots, reproduces their full
465-pair ledger, and matches the recorded escape from vertex3=p_s,
freed wall old3, first hit q. Its eight whites are precisely

    old0, old1, old4, old5, old6, old7, old8, s.

They form A4+4A1. All certificates and pure helper files are pinned
before and after the calculation. Only the existing exact arithmetic,
field-square, matrix and native-angle helpers are imported; none of
their old input, producer, checker or search entrypoints is called.

## Declared finite scope

The worker reconstructs all63 finite short contact profiles, checks
their nine-root Gram matrices positive definite and their exact Schur
constants, and removes only the zero profile excluded at every height
by the proved modulo5 obstruction. It searches the remaining62 profiles
for native integral norm-two roots in the box

    0 < x <= 600,  |x'| <= 2.

For x=a+b*phi, this box implies b>=0 and b<270. The implementation
deliberately uses the larger b range0..300. With
m=floor((sqrt(5)-1)*b/2), the six choices a=m-2..m+3 contain every
possible integer a satisfying |x'|<=2. The exact box test then rejects
the extra values. The implementation asserts1806 coefficient pairs
and1073 accepted x values, hence66526 profile/x pairs. These are
source assertions, not claimed execution results.

For each profile the worker first tests the necessary native congruences

    n=(x+u*epsilon)/t in O,  y=(beta-phi*x)/2 in O,

then tests whether 5D*n^2+K is a square in Q(sqrt5). Both signs of its
square root A are retained. It reconstructs U,T,y,z using the reviewed
Pell identity and tests all native coordinates, the norm and reflection
integrality. An integral norm-two vector in this even O-lattice is
primitive: a nonunit common divisor would have its square divide1.
There is no long-root branch and no assertion outside the declared box.

## Geometric acceptance and retained failures

Every integral root receives all31 signed native-pair tests, rather
than just its eight finite-corner tests. It must cut the actual recorded
spacelike escape. The worker then checks the exact first endpoint on
that edge: positive parameter before the old first hit, future timelike
square, every old inequality, the exact active set, and a positive
definite nine-wall Gram matrix. Old-center positivity is not a rejection;
a common strict center is constructed from the edge vertex and the old
strict center, and projects to a strict witness on the new facet.

The seven named old finite vertices are checked separately. For each
old compact cap r,w, an ultraparallel pair and a strict point on its
hyperplane certify preservation of the entire cap hyperplane. A right
cap pair is not falsely declared lost: it is retained among compatible
cuts with cap preservation explicitly uncertified. Only roots carrying
all seven strict vertex tests and both full-cap certificates enter the
`survivors` list. A compatible cut at a nonfinite endpoint, or a finite
endpoint without the requested retention certificates, remains visible
under `other_compatible_cuts`.

Every integral outcome contributes to a deterministic digest. Rejected
previews are capped at128 records with an explicit omitted count. Rich
retained outcomes are capped below850000 bytes, and the entire result
must be below1MiB; hitting either cap fails rather than certifies an
incomplete search. Candidate roots are individually compatible with
the current chamber; no mutual compatibility between different
survivors is inferred.

## Execution handoff

Use the coordinator's bounded parent-lock research wrapper, one thread,
45-second worker alarm,59-second outer cap and1GiB memory ceiling,
after its ordinary disk/RAM start gates pass. The worker owns no inner
lock, launches no process, and writes only its exclusive-create result

    results/h19_thirty_one_edge_short_profiles_20261019.json

under `--write`. A separate `--verify` invocation recomputes and compares
the complete result. Neither invocation has occurred in this lane.
Static AST parsing and streamed source-hash inspection are preparation
checks, not execution or mathematical verification.

Any returned root still needs the coordinator's separate thirty-two-wall
certificate and global edge-walk/volume investigation. No new mirror,
finite-volume H9 example, all-height reflection census, nonreflectivity
statement or dimension bound is claimed by this preparation note.
