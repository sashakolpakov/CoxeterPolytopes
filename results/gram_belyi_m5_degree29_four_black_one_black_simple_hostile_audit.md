# Hostile audit: final degree-29 four-black `BW/WB` no-go

**Verdict (2026-08-10): PASS.**  No correctness blocker was found.  The
canonical normalization, the two endpoint orientations, the complete
resultant/tower census, the characteristic-zero reduction argument, and the
four-black/global count boundary all survive independent review.  The global
degree-`29` problem remains open.

## Frozen artifacts

The audit pins the following SHA-256 hashes:

| artifact | SHA-256 |
|---|---|
| `scripts/certify_gram_belyi_m5_degree29_four_black_one_black_simple.py` | `db7075b67d3bd696d95d4bea60d65144d4e36a92fe0d2e1a4d5df86efce5bda8` |
| `data/enumerations/gram_belyi_m5_degree29_four_black_one_black_simple.json` | `21b06f6714561937b494713d6fef36638fea4dd4be36ceb2bd75b294c0732bef` |
| `results/gram_belyi_m5_degree29_four_black_one_black_simple.md` | `8324f4b08d92d6892427a201db60bdb5736fa61bfafb62e9f545f19eb54a1f70` |
| `tests/test_gram_belyi_m5_degree29_four_black_one_black_simple.py` | `ed9d17ff750f393b40839fc312904c3482caaa5990e9a84b1148e9274e65099f` |

The canonical JSON encoding is `9,308,646` bytes.  The independently encoded
`.rows` array, including its final newline, has SHA-256
`5470c827a362663eb690641b57fff6cdd5a738c1d24e8eff4e5bb2f46caa8ea7`.

## 1. Normalization and orientation coverage

The ledger convention `B(lambda)=0` makes the marked point black.  For a
degree record `(a,b,c,d)`, where `a` is the marked-root multiplicity, `b` is
the black-endpoint multiplicity, and `c<=d` are the two remaining unlabelled
multiplicities, the two orientations are exactly:

- `BW`: use `y=x/lambda`; the black endpoint is `0`, the marked root is `1`,
  and the white endpoint is `z_BW=1/lambda`;
- `WB`: first reflect `x -> 1-x`, then use
  `y=(1-x)/(1-lambda)`; the black endpoint is again `0`, the marked root is
  `1`, and the white endpoint is `z_WB=1/(1-lambda)`.

Thus in both orientations

\[
 G(y)=y^b(y-1)^a(y-u)^c(y-v)^d,
\]

so the script's canonical weight order `(b,a,c,d)` is correct.  Independent
enumeration of all positive tuples summing to `29` with `c<=d` gives exactly
`1,729` distinct records.  The two ordered orientations therefore give
`3,458` systems, with no missing endpoint ordering.

## 2. Critical equations and finite-field coverage

After removing the forced black-root powers from `G'`, the residual derivative
is a cubic `D`.  If

\[
 G \bmod D=R_2y^2+R_1y+R_0,
\]

then `R1=R2=0` is precisely the scheme-theoretic condition that every
non-root critical point has the same value.  It remains valid when `D` has
repeated roots.  The endpoint condition is only `G(z)-R0=0`; no simplicity
or derivative condition is used, so the calculation legitimately contains
both simple and critical white endpoints.

The homogeneous degree bounds are `deg R1<=28` and `deg R2<=27`, hence the
projective Bezout bound for
`A(u)=Res_v(R1,R2)` is `28*27=756`.  Every selected reduction retains exact
total degrees `28,27`, bidegrees `(27,27)`, and resultant degree `756`.

For every irreducible factor `f(u)` of `A`, the producer works in
`K_f=F_p[u]/(f)`, computes the complete specialized gcd of `R1,R2` in `v`,
and factors that gcd.  For every irreducible factor `g(v)`, it evaluates both
endpoint differences in the field `K_f[v]/(g)`.  A nonzero field element
stays nonzero under every algebraic-closure embedding, so this checks entire
orbits rather than only rational points.  Resultant and gcd multiplicities do
not affect point coverage, but their weighted degrees are retained.

The archive-wide audit found:

- `1,729` rows in exactly the independently enumerated order;
- weighted `u`-factor degree `756` in every row;
- at least one specialized critical equation of `v`-degree `27` on every
  `u`-factor;
- `28,159` irreducible `(u,v)` tower orbits, each with endpoint flags
  `[true,true]`;
- prime histogram `101:1529, 109:174, 131:23, 139:2, 149:1`;
- maximum `u`-factor degree `600`, maximum relative `v`-factor degree `2`,
  and maximum absolute tower degree `600`.

As a calculation independent of the producer's tower helper, one genuine
component for weights `(1,1,1,26)` was replayed directly in
`F_101[T]/(T^2+61)` with `u=-1`.  Direct polynomial division gave
`R1=R2=0` and `R0=54T`.  At `z_BW=18` and `z_WB=7`, the two differences were
respectively `48+89T` and `75+3T`, both nonzero.

## 3. No escape under characteristic-zero reduction

All selected primes are split for `Q(sqrt(5))` and avoid `2` and `29`.
The critical equations therefore have coefficients integral at the selected
place.  Exact resultant degree `756` modulo the prime makes the
characteristic-zero leading coefficient a local unit.  Dividing by it makes
`A` monic, so any putative solution coordinate `u` is integral.

For the residue factor containing the reduction of `u`, the stored
degree-`27` condition says that a leading coefficient of `R1` or `R2` in `v`
is nonzero modulo the residue field.  It is consequently a local unit, and
the corresponding equation makes a common root `v` integral.

Finally,

\[
 N(\lambda)=\frac1{29},\qquad N(1-\lambda)=\frac{16}{29}.
\]

Thus `lambda` and `1-lambda`, and hence both fixed endpoint coordinates, are
local units at every selected prime.  A characteristic-zero solution of the
critical and endpoint equations would therefore reduce to one of the tower
orbits checked in the archive, where the relevant endpoint difference is
nonzero.  This contradiction excludes both orientations.  There is no
uncovered branch at `u`, `v`, or an endpoint infinity.

## 4. Full deterministic replay

An instrumented invocation of the producer's unchanged `main(--check)`
delegated every prime attempt to the original `certify_one` and only emitted
progress/RSS records.  It rebuilt all `1,729` successful rows from `1,959`
prime attempts, then compared the complete canonical byte string with the
frozen archive.

The producer printed:

```text
status = EXACT_FOUR_BLACK_COMPLETE_NO_GO
rows replayed = 1,729
critical orbits checked = 28,159
archive sha256 = 21b06f6714561937b494713d6fef36638fea4dd4be36ceb2bd75b294c0732bef
wrapper status = SCRIPT_MAIN_SUCCESS
```

Elapsed wall time was `1,094.107` seconds (`1,094.22` seconds from the outer
timer).  Peak in-process RSS, including construction and final archive-byte
comparison, was `133,627,904` bytes = `127.438 MiB`, well below `0.5 GiB`.

The outer `/usr/bin/time -l` command returned status `1` only after the
successful child run because the sandbox denied its auxiliary
`sysctl kern.clockrate` query.  The producer itself completed, printed its
success JSON, matched the expected hash, and the instrumentation printed
`SCRIPT_MAIN_SUCCESS`; the timer-wrapper status is not a certificate failure.

## 5. Tests, bytes, and count boundary

The focused test file completed with `6 passed in 2.56s`.  Independent scans
confirmed valid UTF-8, final LF, and absence of NUL, backspace, vertical-tab,
form-feed, and carriage-return bytes in the script, archive, theorem report,
and test.

Cross-ledger arithmetic also closes exactly:

\[
 6734+2560+3458=12752
\]

four-black coefficient systems, and

\[
 23660+8960+10374=42994
\]

four-black marked passports.  The final block contributes three simple-white
partitions per oriented system, so `3*3,458=10,374`.  Globally,

\[
 37969+10374=48343,
 \qquad 42610678-48343=42562335.
\]

It is therefore correct to set `complete_four_black_no_go_claimed=true` and
to conclude that any degree-`29` witness must have at least five distinct
black vertices.  It is equally important that the archive leaves
`complete_degree_29_no_go_claimed=false`: all strata with at least five black
vertices remain open.

## 6. Auditability note

The JSON records monic resultants, component degrees, multiplicities, and
endpoint booleans, but not every irreducible factor polynomial or a Bezout
witness for each nonvanishing assertion.  Archive-only verification is
therefore weaker than deterministic producer replay.  This is not a
correctness blocker here because the complete exact replay succeeded within
the stated resource bound.

```text
hostile_audit_verdict = PASS
mathematical_blockers = 0
implementation_blockers = 0
four_black_status = IMPOSSIBLE
global_degree_29_status = OPEN
```
