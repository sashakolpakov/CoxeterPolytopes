# Independent hostile audit: golden Leech norm-128 manifest

## Verdict

**PASS, with the claim restricted to the archived first necessary prefixes.**

An independent backend that does not import the owner manifest reconstructs
all 28 prime-ideal atoms, all 47 squarefree coefficient-ideal directions, all
primitive local mirror-norm classes and spatial gates, every first necessary
binary-Leech packet, and every strict fixed-direction witness field by field.
The tight row is again `p41-`, with exactly the archived coordinate, margin,
and hyperideal excess.

This is not a compact-polytope example and not a global nonreflectivity
theorem.  It proves only that each finite wall prefix through its first
necessary packet still contains an open hyperideal patch.  Later mirrors may
cut that patch.

## Independent implementation

The clean-room replay is in
`src/coxeter_polytopes/search/golden_leech_manifest_cleanroom.py`.  It imports
only exact arithmetic in `Q(phi)`, not the owner manifest.  In particular it
independently:

1. detects splitting by factoring `X^2-X-1` modulo each rational prime;
2. finds and canonicalizes the two narrow-principal generators over every
   split prime;
3. forms squarefree ideal products by their generator lineage;
4. derives mirror classes directly from that lineage rather than refactoring
   the owner coefficient;
5. enumerates controls from separately certified bounds at both real
   embeddings;
6. solves the binary Leech Gram equations and exhausts every integral binary
   combination that could have square below four; and
7. recomputes the fixed-direction `3x3` Gram bounds and strict witnesses.

The four Leech direction norms are not accepted from a shell-count assertion:
the audit stores explicit standard-MOG vectors of norms `4,6,8,10`.

## Arithmetic census

For `O=Z[phi]`, `phi^2=phi+1`, the ring is norm-Euclidean, hence has class
number one.  Its units are `+/-phi^n`.  Because `N(phi)=-1`, the narrow and
ordinary class groups agree.  The totally positive units are precisely
`phi^(2n)`, and every one is a unit square.  Thus each squarefree ideal product
has one canonical sign-`(+,-)` coefficient in the selected unit-square
fundamental domain.

The independent factorization gives:

- inert primes `(2),(3),(7)`, with ideal norms `4,9,49`;
- the ramified prime over `5`, with ideal norm `5` and generator `2+phi`;
- both primes over each of
  `11,19,29,31,41,59,61,71,79,89,101,109`.

That is `3+1+2*12=28` atoms.  Direct squarefree-product enumeration through
ideal norm 128 gives 47 lineages.  The two primes over 11 remain independent:
their coefficients are `-1+3phi` and `-2+3phi`, while their product is the
rational-content direction `-11+11phi` of ideal norm 121.

These 47 rows are squarefree representatives of coefficient directions in
`K^*/K^{*2}`.  Removing nonunit square content does not generally preserve an
integral lattice, so 47 is not a count of integral diagonal hosts.

## Primitive local mirror classes

Let `r=(x,y)` be primitive and `q=r^2`.  Reflection integrality says

\[
(q)\mid 2\,(a x,(y,\Lambda)).
\]

At a prime away from `2a`, a positive valuation of `q` would force both `x`
and `y` into that prime, contradicting primitivity.  At an odd prime `P`
dividing squarefree `(a)`, a positive valuation forces `y in P Lambda` and
`x` to be a `P`-unit.  Then

\[
q=-a x^2+y^2\equiv-a x^2\pmod {P^2},
\]

so the valuation is exactly one.  The choices at distinct odd prime ideals
are independent, including conjugate primes over the same rational prime.

The hostile dyadic check confirms the owner rule.  The prime `(2)` is inert,
with `v_(2)(2)=1`, while Leech evenness gives `y^2 in 2O`.

- If `(2)` does not divide `(a)` and `x` is a unit, then `v(q)=0`.
- If `(2)` does not divide `(a)` and `v(q)>0`, the norm equation forces
  `x in (2)`; primitivity makes `y` primitive, and reflection integrality
  gives `v(q)<=1`.  Hence `v(q)=1`.
- If `(2)` divides squarefree `(a)`, evenness gives `v(q)>=1`.  When `y` is
  primitive, reflection integrality gives `v(q)<=1`; when `y` is not
  primitive, `x` is a unit and the norm equation gives `v(q)=1` directly.

Therefore the dyadic exponent is `{0,1}` off `(a)` and exactly `{1}` on
`(a)`.  It imposes no spatial divisibility gate.  This reproduces `q=1,2`
for `a=phi-1`, only `q=2` for `a=2(phi-1)`, four classes for the norm-five
row, and eight independent classes for the product of the two primes over 11.

## Height zero and first necessary packets

For a primitive height-zero vector, unimodularity makes the pairing ideal of
`y` equal to `O`.  Evenness and reflection integrality then force
`(y^2)=(2)`.  Since every totally positive unit is a square, unit rescaling
would give `y^2=2` and trace four.

Writing `y=u+phi*v` gives

\[
\operatorname{Tr}(y^2)=u^2+(u+v)^2+2v^2.
\]

The Leech minimum is four.  If `v=0`, the first two terms coincide and sum to
at least eight; if `v` is nonzero, the last term alone is at least eight.
Thus height-zero mirrors do not exist.

For positive height, the clean-room search independently enforces

\[
2ax/q\in O,
\qquad T=q+ax^2\succeq0,
\qquad H=ax^2/q,
\]

and, for every odd spatial divisor `d`, solves `z^2=T/d^2` with `y=dz`.
All 47 first nonempty necessary packets occur below the exact height cap
5001.  Each has one local component.  The independent rows agree with the
owner rows in coefficient, lineage, mirror classes, height, control, gated
target, and every surviving Gram triple.

The packets are necessary super-packets: their triples need not be realized
by Leech vectors, need not be primitive roots, and need not survive Vinberg
acceptance.  Enlarging the true wall set in this way is safe for the negative
prefix obstruction below.

## Witness normalization and the tight row

For a component with spatial gate `y=dz`, the affine wall inequality is

\[
c\,(D_0,z)\le ax/d.
\]

This confirms that the wall bound is `B=ax/d`, not `ax` and not an extra
multiple of `a`.  Canonicalizing `q` by a totally positive unit square merely
unit-rescales the root and does not change the ideal gate `d Lambda`.

For each explicit Leech direction `D_0` of square `D in {4,6,8,10}`, the
clean-room code exhausts all integral pairings allowed by positive
semidefiniteness of the `3x3` Gram matrix and the Leech minimum.  Its maximum
`M` is therefore an upper bound even when a necessary binary Gram type is not
Leech-realized.  Taking

\[
c_{\max}=\min_i (a x_i/d_i)/M_i
\]

gives a uniform wall-safe coordinate.  A strict shrink has both positive wall
margin and `D c^2-a>0` in all 47 rows.

The minimum limiting excess occurs uniquely at

\[
a=-5+6\phi,\qquad \text{lineage }p41-,\qquad H_1=122+198\phi.
\]

Its sole component has 12 necessary Gram types and

\[
D=4,\quad M=17+26\phi,\quad B=16+30\phi,
\]

\[
c_{\max}=B/M=-\frac{92}{55}+\frac{94}{55}\phi,
\]

\[
4c_{\max}^2-a=
\frac{3373}{121}-\frac{10398}{605}\phi>0.
\]

The strict shrink `c=(141/142)c_max` gives

\[
B-cM=\frac8{71}+\frac{15}{71}\phi>0
\]

and

\[
4c^2-a=
\frac{16807457}{609961}
-\frac{51937482}{3049805}\phi>0.
\]

These are exactly the owner artifact values.

## Artifact freshness and reproducibility

Before the scope clarification, regenerating the owner JSON from the current
backend produced a byte-identical file.  The owner generator and JSON were
then updated together to state explicitly that the rows are K-square
directions rather than integral-host classes.  The clean-room JSON separately
records that scope and `owner_backend_imported: false`.

Reproduce the independent certificate and both test suites with

```console
PYTHONPATH=src python3 scripts/audit_golden_leech_low_norm_manifest_cleanroom.py
PYTHONPATH=src pytest -q tests/test_golden_leech_manifest.py \
  tests/test_golden_leech_manifest_cleanroom.py
```

The companion SHA-256 manifest is
`data/enumerations/golden_leech_low_norm_manifest_independent_audit.sha256`.
