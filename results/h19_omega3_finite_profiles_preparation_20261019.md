# Prepared bounded omega3 finite-profile search: not executed

19 October 2026. **SOURCE ONLY.** I wrote and read the new source,
inspected metadata/hashes, and parsed its AST without importing it.
No scientific body, root enumeration, optimizer, Lean process or
build ran. The coordinator owns review and any later guarded launch.
No new search certificate or scientific result exists from this lane.

Source: `scripts/search_h19_omega3_finite_profiles_20261019.py`,
290 lines,15,350 bytes, SHA256
`da6be4a933eba17b74701b1aac0699123410427cab4440045e4822c285b636e0`.
Proposed output: `results/h19_omega3_finite_profiles_20261019.json`.

## Mathematical domain and filters

The pinned input proof is
`results/compact_h9_discriminant19_omega3_profile_congruences_20261019.md`,
SHA256 `0c9a29fc83a7c1abac1330b72b771612423c0f0ea9009a63daf6f71d47a038fd`.
The fixed initial old4-edge has white support A1(binary),
A4(E8 nodes0,1,2,7), and A3(E8 nodes4,5,6).

The source considers every choice of at most one0/1/phi contact in
each white component:3*9*7=189 labeled choices. The exact Schur
test gives29 positive-definite profiles including zero. It then
applies the proof's prime-D square test, and the divided-by5 test
when R belongs to5O. The asserted counts are29->22->18, retaining
the actual node labels. The latter test uses A=20T+R in5O to
reduce -2x²=K/5 modulo sqrt(5); it is not just a bare K-square test.
No finite classification computation is substituted for that proof.

For a selected profile, the exact equations are

    beta=phi*x+2y,
    z=C^{-1}(T*e3+eta), eta_3=0,
    R=(C^{-1})_3*eta, Q0=eta^T*C^{-1}*eta,
    A=20T+R,
    A²-10D*x²=K=40-10beta²-20Q0+R².

For every field-square solution both A signs are tried. The source
requires T=(A-R)/20 and y=(beta-phi*x)/2 in O, reconstructs the
entire native vector, and rechecks its norm2, coordinate integrality,
and every prescribed signed white contact. It never treats the
scalar equation alone as a valid root.

## Complete stated bounded domain

The only height restriction is

    x=a+bphi, 0<x<=600, |x'|<=2.

Exactly as in the independently read pinned omega2 worker, it uses
b=0..300 and a=m-2,...,m+3, with
m=(isqrt(5*b*b)-b)//2=floor(b/phi). These1,806 conservative pairs
contain the complete requested domain, followed by exact field-sign
filters. The conjugate bound is conservative for canonical norm2
in this positive-conjugate host. There is no unproved Pell-period,
T bound, approximate square test or floating coefficient floor.

The previous executed worker used this identical x domain and
retained1,073 x, suggesting19,314 exact field-square tests across
18 profiles. This is an estimate inherited from that domain,
not a timing, count or outcome measured for the new source.

## All27 retained walls and actual endpoint checks

The source reconstructs the native form, its inverse, and every
current root. It uses the19 literal old prefix roots, the seven
new26-wall covectors, and the smaller q from the fixed27-wall
certificate. It verifies all coordinates, norms and full-lattice
reflection ratios. It does not run the old or new checker body.

Every integral candidate, even one already failing a wall or cut,
receives ALL27 exact pair tests. Positive, parallel and forbidden
finite-angle products are rejected. The pinned field-angle theorem
makes this finite menu complete for the native norm classes. The
old center is only diagnostic; a smaller p+epsilon*c supplies a
strict future timelike center if the old c does not pass the new root.

The saved27-wall old4 first-hit point is spacelike. A candidate must
cut it strictly. For each survivor the source reconstructs

    tau=D*x/A,
    endpoint=p+tau*omega3,
    endpoint²=-2D*K/A²<0.

It verifies0<tau<tau_old, the future time cone, EVERY27 old-wall
inequality, the new root equality, and exactly the eight old white
walls active. The Schur-positive profile therefore corresponds to
an actual finite endpoint only after all these checks pass.

The r,w pair kinds are retained as diagnostics; no new H8-facet,
vertex-preservation, full-cone or finite-volume claim is embedded
in the search's status. Further consequences need separate review.

## Pins, common arithmetic type, output and resource limits

The imported pinned omega2 source supplies ONLY Q and pure arithmetic
helpers. Its inputs(), profiles() and run() are never invoked. Using
its same Q class avoids mixing separately generated quadratic-field
classes. Its underlying exact field/inverse helper is also a direct
runtime pin.

The27-wall certificate and checker,26-wall certificate, native prefix,
pure omega2 source/helper, and the18-profile proof are direct pins.
The27/26/omega2 saved certificate source maps are all rehashed before
use. The four small JSON inputs total78,325 bytes. No giant inherited
JSON, optimizer, external package or reflection-orbit queue is used.
The27-wall generation certificate is pinned independently of its
still-pending separate replay at source-preparation time; no replay
success is assumed or fabricated by this worker.

The source has no inner flock. It must run under the coordinator's
single-lock bounded wrapper, with the required40GiB disk/6GiB RAM
start gates and all ongoing resource/owner-death guards. The child
adds a45-second alarm, and the parent must keep its59-second/1GiB
caps. The parent must review and schedule this, not infer permission
to launch from the command below:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/search_h19_omega3_finite_profiles_20261019.py --write
```

Generation uses exclusive output creation; replay uses --verify and
recomputes the entire result before comparing canonical JSON. Retained
solution records stop at512KiB compact JSON, and the final indented
certificate is capped below1MiB before any file is written. Reaching
a cap fails explicitly rather than declaring a partial search complete.

Any eventual PASS or empty result covers only these18 profiles in
the specified x domain against the fixed27-wall seed. It does not
prove all-height nonexistence, other-star completeness, nonreflectivity,
finite volume, compact H9, or a new dimension record.
