# Prime811/k5: complete bounded transverse-edge checker specification

21 October 2026. UNREVIEWED AND UNEXECUTED source preparation. No checker
import, generation, replay, scientific job, Lean or TeX was performed. This
specification and its worker require independent full source review before
root admission. They change no current host status, example, bound, paper or
accepted receipt. A completed future execution must be recorded separately.

Worker: `scripts/check_compact_h9_prime811_first_edge_20261021.py`,
SHA256 `b72f3274cd5cb953f424d24bf358e22b355cc7ff005c02a0edcff4c6597174d0`
(24121 bytes,599 lines). The worker does not import another repository
module. Only Python standard-library arithmetic, hashing and bounded I/O
are used. Its `--verify` means a full fresh computation into a NEW JSONL
file, not comparison with an existing certificate.

## Mathematical statement and explicit inputs

Let O=Z[phi], phi^2=phi+1, and use the literal lattice

    L=E8(O) orthogonal-sum [[-2phi,a],[a,2]],
    u=3+5phi, a=5+8phi, b=59+89phi, D=phi*b=89+148phi.

The native E8 edges are02,13,23,34,45,56,67; outward old normals are
`-e_j`. The actual starting point is p1=(2,-a;0), of square -2D, with
whole germ E8+A1(c), c=(0,-1;0). Release -e7, retaining -e0,...,-e6,c.
The outgoing edge is p1+tau*omega7. The explicitly integral upper is

    Vb=5556+8989phi, Xb=62229+100689phi,
    b*Vb^2-phi*Xb^2=1,
    mb=Xb*p1+b*Vb*omega7, tau_b=phi*Xb/Vb<sqrt(D).

The worker derives the exact least collision at or before tau_b among
ALL possible first-simple profiles below. It then derives the complete
oriented new simple germ using finite cone witnesses. A literal-host
nonreflectivity conclusion is permitted ONLY if that unique simple wall
is exactly mb and the original-lattice integral period passes.

The following are explicit external mathematical inputs, supplied by the
pinned sources, not encoded as computationally proved geometry:

* The full Bugaenko H8 base calibration, the actual first binary edge,
  and completeness of their actual ambient germs.
* The golden ring/PID/unit facts and primitive reflection norm result;
  since N(b)=811 is prime, norm classes are exactly2 and2b.
* Complete finite Coxeter classification at E7+A1, finite native angle
  classification, prime-norm angular isolation, and local finiteness of
  the full integral arrangement. Finite chamber corners are simple.
* An integral future isometry matching full oriented actual germs
  preserves their full arrangement chamber. The reflection group acts
  freely on its chambers; the ambient arithmetic group is discrete and
  a lattice when the infinite-volume conclusion is drawn.

These inputs are distinct from a formal positive definite Gram, an
arbitrary root list, a selected invariant lattice or a similarity class.

## Complete finite domain and termination

Write z=A+Bphi with A,B integers. Positive orientation and the conjugate
strip force B>=0. The five first-simple profiles are:

| Profile | Normal and exact equation | Conservative scalar strip | Maximum B |
|---|---|---|---:|
| U0 | X p1+Y omega7, Y^2-DX^2=1 | 0<X<16000, abs(X')<2/3 |8000|
| U1 | (U/2)p1-c/2+Y omega7,4Y^2-DU^2=3 |0<U<28000,abs(U')<9/8,U=1 mod2|14000|
| Uphi | (U/2)p1-phi*c/2+Y omega7,4Y^2-DU^2=3-phi |0<U<19000,abs(U')<5/4,U=phi mod2|9500|
| E8 | X p1+omega6+Y omega7,(2Y+3)^2-4DX^2=1,Y>=0 |0<X<8000,abs(X')<1/3|4000|
| B | X p1+b V omega7,bV^2-phi X^2=1 |0<V<20102,abs(V')<251/500|10051|

The U congruences express integrality of the ORIGINAL binary second
coordinate; a=1 modulo2 for k5. They are not the k4 parity recipe.
Collision comparison with tau_b gives the sharper square bounds in the
candidate source; the worker independently checks that every displayed
conservative constant contains those bounds, as well as the coefficient
maxima. It also verifies the fixed upper is actually admitted by the
domain and square test; it is not enough that some unrelated wall survives.

For a strip |z'|<c, its integer A lies in B(phi-1)+(-c,c). The exact
floor of B(phi-1) is `(isqrt(5*B*B)-B)//2`, including B=0. Trying this
floor plus -1,0,1,2 covers every strip because c<=5/4. Exact field signs
then retain precisely the strip members. At most2,3,3,1,2 integers A
per B survive, respectively for U0,U1,Uphi,E8,B. Therefore the number of
strip pairs counted across the FIVE profiles is at most

    16002+42003+28503+4001+20104=110613.

The raw four-integer trials before strip filtering are a slightly larger
constant-size loop; they are not mislabeled as the110613 admitted domain.
Every loop is finite with explicit endpoints, every integer operation is
exact, and no elapsed time, Pell-prefix length or root height is a domain
cutoff. Termination is thus mathematical, while successful completion
inside the external45/59-second limits remains unmeasured.

### E8/U0 sharing, including the actual ties

The E8 equation is the U0 equation with X0=2X,Y0=2Y+3. The worker searches
U0 once. Each even X0 is also counted as an E8 strip member; E8 admission
additionally requires Y0-3 even and nonnegative. Thus it covers the full
E8 scalar domain without a second square computation. Both candidates
have collision D X0/Y0, but are DISTINCT mirrors. They are retained as
ties, not merged by row order. The identity

    m_U0=2 m_E8+sum_{j=0}^6 (2,3,4,6,5,4,3)_j*(-e_j)

is checked whenever this pair is reconstructed. It makes U0 redundant
in the E8 simple cone. Their positive product forbids treating both as
simple walls of the acute chamber.

## Exact arithmetic, mirror validity and firstness

`F(a,b)` means a+b phi with rational coefficients. Multiplication uses
phi^2=phi+1; inverse is conjugate divided by rational norm. Physical
signs reduce to the sign of 2a+b+b sqrt5 by rational square comparisons.
The conjugate is (a+b)-b phi. No floating arithmetic is used.

For w=A+Bphi in O, a square root v=a+bphi must have norm s with
s^2=N(w), and

    (2a+b)^2=Tr(w)+2s,  b^2=(Tr(w)-2s)/5.

Testing BOTH rational norm signs, both integer square-root signs and
parity reconstructs every possible integral square root. Exact squaring
validates the witness. This includes the rational, pure-radical and
trace-zero cases. If w is nonintegral, it has no required O-valued root.
The Ueta division by4 is checked exactly. Positive physical orientation
selects the one positive square root.

Every admitted mirror is checked in the original 10-dimensional basis:
all coordinates integral, squared norm2d, the reflection covector
`Gm/d` integral, expected retained-wall products, negative product at p1,
nonpositive product against released -e7, exact collision and
0<tau<=tau_b. Primitivity follows from the inherited even-lattice
norm result and squarefree ideals d=1,b: a nonunit common coordinate
divisor would have its square divide d. The worker is not claiming to
perform an unrelated golden Euclidean gcd algorithm.

Every admitted mirror is streamed. Only the least tau and all its ties
are retained. For fixed profile, its positive norm equation and collision
determine the scalar uniquely, so at most one mirror per profile can tie
at a fixed tau, hence at most five minimum ties. An earlier unlisted
arrangement mirror would yield an earlier actual simple opposite, which
the finite-profile theorem puts in this exhausted domain. Thus the
minimum is the ACTUAL endpoint under the stated geometric inputs; it
is not inferred from the upper alone.

## Full simple germ and conditional period

At the minimum, the point's square is checked negative. For each tied
candidate q, let the nine-element basis be the eight retained normals
plus q. The worker checks positive exact Gram pivots and expresses every
tied normal in this basis by exact Gaussian elimination. It reconstructs
each normal and records all coefficients and their signs. A proposal is
accepted only if EVERY tied normal has nonnegative coefficients.

The actual new simple is among the ties. If q's simple cone contains all
these roots, it in particular contains that actual simple root. Since q
is itself an arrangement mirror with the correct orientation, the two
chamber halfspace descriptions coincide. This identifies the full germ.
Exactly one accepted proposal is required for an actual-simple-germ
conclusion. Zero/multiple accepted proposals are explicit unresolved
results, not heuristically repaired.

Only if q is exactly mb at tau_b does the period branch run. It constructs
the candidate's T on the ORIGINAL native basis and checks:

* every matrix entry belongs to O, determinant1 and exact preservation
  of the entire original Gram;
* T(r)=mb,T(g)=c,T(Pold)=Pnew with future timelike Pnew;
* all nine old simple normals map to the nine new simple normals;
* exact trace253415+410024phi>10.

Determinant1 over O gives an integral adjugate inverse. The trace rules
out finite order. The pinned full-arrangement/germ theorem supplies
chamber stabilization; its reflection group has trivial intersection
with the chamber stabilizer. This proves infinite index and, under the
explicit arithmetic lattice input, infinite chamber volume. This branch
concerns only the literal k5 host. The worker cannot transfer that result
to its full form class or claim a new compact H9 example.

If a different actual first wall wins, its exact full germ is saved and
no period/exclusion conclusion is asserted. Further actual edges would
then remain necessary. A new compact example would require an independent
complete facet/cone-ray compactness certificate, not merely this edge.

## Pins, bounded execution and certificate schema

The worker hard-pins candidate f814a811, contract c91abbb5, all six
candidate input sources, the unchanged bounded runner4d8a2a74, its
monitor27f14501, the signal-safe shim c0ddc98c and the base resource
implementation `run_lean_batches_20261018.py` at97822ecf. The last
file implements the actual40/6 start and8/3 abort gates and is part of
the reviewed transitive resource closure. The caller additionally supplies
the exact reviewed worker SHA and this specification's SHA; all are
stream-hashed BEFORE and AFTER the mathematical computation. There is
no large JSON input and no repository code import.

An invocation requires `COXETER_CHECK_OWNER` to equal the live direct
parent and verifies the shared exclusive lock is already held. The
unchanged wrapper must use `--parent-lock`, retain its RSS/process
watchdog and apply the fresh40GiB-disk/6GiB-RAM admission gate. The
child additionally has45-second wall and CPU limits and periodically
checks its owner. No weakened gate or uncapped fallback is authorized.
No runtime is predicted by this source-only specification.

After full source review and ROOT admission only, the single-module
command shape is:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py \
  --parent-lock --ledger results/NEW_prime811_edge_execution.json -- \
  scripts/check_compact_h9_prime811_first_edge_20261021.py --verify \
  --expected-source-sha b72f3274cd5cb953f424d24bf358e22b355cc7ff005c02a0edcff4c6597174d0 \
  --spec-sha REVIEWED_SPEC_SHA \
  --out results/NEW_prime811_first_edge.jsonl
```

The destination is an exclusively created new `.jsonl` directly under
results. It is never overwritten or used as an input. Records are:

1. Header: explicit incomplete status, all pins, encoding and domain bound.
2. Mirror: deterministic ID, profile, scalar, positive square witness,
   original-basis normal, half-norm and exact collision.
3. Final result: counters, exact minimum and ties, all Gram/cone witnesses,
   unique-simple flag, conditional period data or its failure reason,
   and confirmation of after-computation pin checks.

Working storage is independent of domain size: counters, one scalar,
at most five ties and fixed9/10-dimensional matrices. Certificate bytes
are streamed with a256MiB hard output cap. No candidate array is loaded
or accumulated. A failure/missing final result is not a completed
certificate, even if earlier records happen to contain the upper.
The complete-result statuses distinguish unresolved germ, a new actual
germ without exclusion, and the fully checked conditional literal-k5
obstruction. An exit0 alone is not a declaration that the obstruction
branch succeeded. Output contains no owner/timestamp-dependent values;
a fresh replay must use another new destination and can then be compared
byte-for-byte externally. No such execution or comparison has occurred.

The pre-review14f633 worker and454af24f specification are byte-preserved
in `source_archive/prime811_first_edge_before_resource_closure_20261021.py.gz`
and `source_archive/prime811_first_edge_spec_before_resource_closure_20261021.md.gz`.
The successor adds ONLY the previously omitted base resource-source pin
and corresponding documentation; no mathematical domain or algorithm
is changed, and neither version has been executed during preparation.
