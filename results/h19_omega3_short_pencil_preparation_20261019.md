# Bounded omega3 short-pencil source and independent obstruction review

19 October 2026. **SOURCE PREPARED, NOT EXECUTED.** This lane performed
only source reading, input metadata/hash inspection and an AST syntax parse.
No root search, matrix calculation program, checker body, solver, Lean or
build was executed. The coordinator owns review and any guarded execution.

Source: `scripts/search_h19_omega3_short_pencil_20261019.py`,255 lines,
11,540 bytes, SHA256
`02750aed919d5d507f6e32816718484e81cf8f5f1570ad783295aa591ae1edbd`.

## Exact finite search and completeness of its domain

The sole pencil is

    r=(2u,-phi*u,t*omega3),  omega3^2=20,
    10t^2-Du^2=1,  D=1+5phi,
    u,t in Z[phi], 0<u<=300, |u'|<=2, t>0.

For u=a+b*phi, b=(u-u')/sqrt5 lies strictly between-1 and151. Hence
0<=b<=150. Since phi' is between-1 and0, a=u'-b*phi' lies in
[-2,b+2]. The source traverses all12,080 pairs in this conservative
integer box and applies the requested inequalities exactly. It never
replaces this box by an approximate Pell recurrence or a numerical cutoff.

For each u it tests whether q=(1+Du^2)/10 is a square in Q(sqrt5).
Writing q=A+B*sqrt5 and t=C+E*sqrt5, the rational norm of q must be a
rational square and C^2=(A plus or minus sqrt(N(q)))/2. Both signs are
tested by integer numerator/denominator square roots; E=B/(2C), with the
C=0 case handled separately. These cases exhaust every field square.
The selected positive square root is then tested for membership in Z[phi].
There is no independent bound on t and no loss from a negative choice of
its conjugate. The elementary necessary congruence u=plus or minus2
modulo sqrt5 is asserted only after a genuine integral solution is found.

Every integral solution is tested against ALL26 roots. Products and the
exact squared normalized products distinguish right angles, requested
native labels3,4,5,6, strict ultraparallel pairs, and rejection cases.
Parallel pairs are explicitly rejected. Every candidate is also required
to cut the saved old4 spacelike first-hit ray. No first positive wall
product short-circuits the other25 tests. Native finite-angle exclusions
beyond this explicit label list are supplied by the written field argument
below, rather than by an approximate angle comparison.

The old strict centre c is **not a rejection filter**. Since

    (p,r)=-2Du<0,

the source halves a positive rational epsilon from1 until
(p,r)+epsilon*(c,r)<0. It then checks the exact square, future orientation
and all27 strict wall products of p+epsilon*c. This procedure terminates
even when (c,r)>=0. The point p is feasible and future timelike, c is
strict and in the same time cone, so every positive mixture is timelike
and strict on all26 old walls. The certificate records the old-centre
product/status, epsilon and the new strict-centre norm/product; surviving
roots also carry the complete centre coordinates and26-product table.

## Input reconstruction and provenance

The10-dimensional native form and its inverse are reconstructed exactly.
The19 retained roots come from the22-root native prefix by excluding
indices10,20,21. The seven newer roots are recovered from the saved
26-wall certificate's `new_covectors` using G inverse; no duplicate
scientific producer is imported. The source checks their saved norms,
coordinate integrality and reflection divisibility. The omega3 column
and the saved old4 first-hit ray are checked against that reconstructed
form. The only imported helper is the pinned standard-library quadratic
field and inverse implementation.

The four direct pins are:

|Input|Bytes|SHA256|
|---|---:|---|
|`results/h19_twenty_six_wall_seed_20261019.json`|10,261|`e310723817ce78bce5c8ed1c11f2a210925d6edc3c6fcb8faf443d5458fc011b`|
|`scripts/check_h19_twenty_six_wall_seed_20261019.py`|12,402|`f0ac1d28225a37be17e9b61f2e436ea2fc3b1511efb32b02212b33539fecb363`|
|`data/candidates/h19_vinberg_prefix.json`|12,259|`91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1`|
|`scripts/check_arithmetic_dn_opposed_fork_models_20261019.py`|10,507|`185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9`|

All seven additional `sources` entries of the26-wall certificate are
also rehashed, including its12KB native prefix, pure helper, constructive
reports and7.5KB complete H8 archive. This gives nine distinct pinned
dependency paths overall. The archived H8 facets/41-vertex counts are
inherited only and are not enumerated here. No giant result file is read.

## Independent written review of the excluded pencils

I independently checked the congruence and angle arguments in
`results/h19_initial_edge_pencil_obstructions_20261019.md`, SHA256
`e86f18aa7a8a64b81a869f70ccef2f1e4d9d0e0ebec167933fcb475b600f58e3`.
They pass within the stated canonical norm/pencil scope:

* O/3O is F9. Its square set is exactly
  {0,1,2,1+phi,2+2phi}. For D=1+2phi modulo3, D inverse is2phi.
  Thus the omega2 short equation forces u^2=phi, and its long equation
  forces u^2=1+2phi. Both are nonsquares. For a canonical long root,
  t=b*tau follows from the spatial reflection covector t*e2/b, so the
  long reduction does not omit an integral-root possibility.
* Modulo the ramified prime over5, phi=3. The omega3 long equation forces
  u^2=3, impossible in F5; the short equation instead requires
  u^2=4 and remains unexcluded. Hence there is no reason to run an
  omega2 or omega3-long pure-pencil search.
* A canonical long/short finite nonright pair would require
  g^2=b*q_m. The only possible labels from cos^2(pi/m) in Q(sqrt5)
  are3,4,5,6,10; the rational norms19,76,19,171,95 are nonsquares.
  A short/short pair permits only products-1 and-phi: labels4 and6
  require sqrt2 or sqrt3 in this field, and label10 has q_m=phi+2
  with rational norm5. Thus the script's finite label test is sufficient
  for this canonical norm-two candidate against the current norm2/2b
  alphabet. The proof uses the inherited completeness of that norm
  alphabet; this small search does not prove it anew.

These facts exclude the whole stated omega2 pencil, not every short
mirror meeting the omega2 edge. A new short wall may contact fixed edge
walls nonorthogonally, and those profiles remain outside this search.
No whole-lattice nonreflectivity or completion impossibility follows.

## Generation and separate replay, coordinator only

The script owns **no inner flock**. It must run under the reviewed wrapper
with `--parent-lock`, retaining the single shared lock,6GiB RAM/40GiB disk
start gates,3GiB/8GiB abort gates,59-second outer deadline and1GiB live
process-tree RSS limit. The script itself has a45-second alarm. Only the
fixed12,080 coefficient pairs, exact rational square tests and26-wall
checks occur; there is no optimizer or scientific dependency. Timing is
not yet measured. The small fixed inputs and lack of a growing search
queue give no expected memory pressure near the1GiB limit.

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/search_h19_omega3_short_pencil_20261019.py --write
```

Separate replay uses the same argv with `--verify`. The certificate path
is `results/h19_omega3_short_pencil_20261019.json`; exclusive creation
refuses replacement, output is capped below1MiB, and replay recomputes
and compares canonical JSON without timestamps or tuple/list differences.

Any PASS means completeness only for this bounded integral pencil with
the stated cut/native-pair requirements. It is neither a general root
census nor a finite-volume, compactness, vertex-preservation, reflection
group completion, nonreflectivity or dimension-bound theorem.
