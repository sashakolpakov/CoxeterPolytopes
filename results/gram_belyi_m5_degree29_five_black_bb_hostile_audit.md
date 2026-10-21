# Hostile audit: degree-29 five-black quotient and `BB` no-go

**Verdict (2026-08-10): PASS.**  No mathematical, implementation, coverage,
or claim-boundary blocker was found.  The exact quotient contains `11,011`
five-black `BB` master systems, and the projective finite-field certificates
exclude all of them.  The other five-black blocks and the global degree-`29`
problem remain open.

## Frozen artifacts

| artifact | SHA-256 |
|---|---|
| `scripts/certify_gram_belyi_m5_degree29_five_black_quotient.py` | `deaaed9f90e550ecd8a2256dc65a2bdc4ca57aea7e18354b7dabfff0d344936f` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json` | `89e2d1fafc10ac5ad91edd6e20ed5e96f6c238b404b8bb303a054f330cbc0e87` |
| `results/gram_belyi_m5_degree29_five_black_quotient.md` | `aea1d5c480ca78ad50d214b60c8b68d1ce196bc4ae25be3bfebc7157c445ac9d` |
| `scripts/certify_gram_belyi_m5_degree29_five_black_bb.py` | `5426a3ea63a0b7bc2f05ece8fec5f8920727da9419f5559e362734fe3c7cc389` |
| `data/enumerations/gram_belyi_m5_degree29_five_black_bb.json` | `817ae3eff77012e7808ecefedc5fafff368e2631ef3a783229a44221901206cd` |
| `results/gram_belyi_m5_degree29_five_black_bb.md` | `1c0ac50efbb292076d669e5ff9da90a0872766be4bb11b385d8305cf48d3e96c` |
| `tests/test_gram_belyi_m5_degree29_five_black_bb.py` | `7477f3c301db0c5082c6698694f0eff832d6757801c560e431fc1d6244cbe882` |

The auxiliary preflight script, which is not used as proof, has SHA-256
`71022c98d345515e864071e07b870fe2c0f9328c930cbe0bde4c605d65b2d4fe`.
The canonical `.rows` array of the `BB` archive, including its final newline,
has independently reproduced SHA-256
`2e8fa0a261b7b1227ff2b5018465d550031b70280887f46ee9c46225c2d6a29a`.

## 1. Independent quotient and canonical coverage

Five black vertices force `25` white vertices and white excess ramification
four.  The five white partitions are

```text
(5,1^24), (4,2,1^23), (3^2,1^23),
(3,2^2,1^22), (2^4,1^21).
```

They give one-white status counts `S/C=5/7` and ordered two-white counts
`SS/SC/CS/CC=5/7/7/7`.  Independent marking of length-five black partitions
and direct ordered/unordered composition enumeration both give

```text
M0 = 1,292   (marked lambda root; four other degrees unordered)
M1 = 4,225   (plus one black endpoint; three degrees unordered)
M2 = 11,011  (plus two ordered black endpoints; two degrees unordered)
```

Consequently the exact raw endpoint-colour census is

| block | raw marked records | binary-status masters |
|---|---:|---:|
| `BB` | 55,055 | 11,011 |
| `BW` | 50,700 | 8,450 |
| `WB` | 50,700 | 8,450 |
| `WW` | 33,592 | 5,168 |
| **total** | **190,047** | **33,079** |

For `BB`, each canonical record is

\[
 (a,b,c,d,e),\qquad a+b+c+d+e=29,qquad d\le e,
\]

where `a` is the marked-root degree and `b,c` belong to the ordered endpoints
`0,1`.  A fresh direct enumeration produced `11,011` unique records in
exactly the archive order.  The raw multiplier is the five possible white
partitions, giving `5*11,011=55,055`.

## 2. Critical equations and fixed projective degrees

For

\[
 G(x)=(x-\lambda)^a x^b(x-1)^c(x-u)^d(x-v)^e,
\]

remove the forced black-root powers from `G'`.  The residual derivative `D`
is a quartic with leading coefficient `29`.  If

\[
 G\bmod D=R_3x^3+R_2x^2+R_1x+R_0,
\]

then a Shabat polynomial necessarily satisfies `R1=R2=R3=0`.  This is the
scheme-theoretic common-critical-value condition and remains necessary when
`D` has repeated roots.  Omitting distinctness inequations only enlarges the
closed system, which is safe for a no-go proof.

Before fixing roots, scaling all root coordinates and `x` shows that `R_j`
is homogeneous of root degree `29-j`.  After fixing `lambda,0,1`, this gives
the initial bounds `28,27,26` for `R1,R2,R3`.  The degree-`28` part of `R1`
in `u,v` is obtained by setting all three fixed roots to zero.  In that
specialization, `G` and `D` are both divisible by `x^2`; Euclidean division
therefore leaves a remainder divisible by `x^2`, whose `x` coefficient is
zero.  The valid fixed bounds are exactly

\[
                 (\deg R_1,\deg R_2,\deg R_3)\le(27,27,26).
\]

These bounds hold in characteristic zero and are used as fixed
homogenization degrees.  Reduction is coefficientwise because every selected
prime avoids `29`; a modular leading term may vanish without invalidating the
argument.  Across the archive, the only actual total-degree patterns are
`(27,27,26)` for `10,903` rows and `(26,27,26)` for `108` rows.  Every row
has bidegrees `((26,26),(26,26),(26,26))`, consistent with the fixed bounds.

Independent scalar specializations at three `(u,v)` pairs reproduced direct
univariate division of `G` by `D`, checking the multivariate remainder
construction separately from the stored certificate.

## 3. Affine Singular unit criterion

On the affine chart `w=1`, each certificate sends the three exact reductions
to Singular over `F_p[u,v]` and computes a reduced standard basis.  The
condition

```text
size(J) = 1 and J[1] = 1
```

is precisely the assertion that `(R1,R2,R3)` is the unit ideal.  By the weak
Nullstellensatz this excludes common affine zeros even over the algebraic
closure, not merely over `F_p`.

For an independent representative at `p=101`, an alternate Singular program
returned basis size `1`, generator `1`, and normal form
`reduce(1,J)=0`.  Thus the producer's exact string test cannot be a false
positive for that component.  All frozen rows store both
`affine_standard_basis=[1]` and `singular_stdout="UNIT"`, and the full replay
recomputed them.

The audited environment is Singular `4.4.1` (Homebrew build `4.4.1p5_3`),
FLINT `3.6.0`, and python-flint `0.9.0`.

## 4. Complete infinity-chart coverage

Let `L_i(u,v)` be the degree-`27,27,26` part used at `w=0`.  Projective
infinity is a copy of `P^1`.

- On `u!=0`, set `u=1`.  The producer takes the gcd of every nonzero
  `L_i(1,v)`.  Gcd one is equivalent to absence of a common root over the
  algebraic closure.  Omitting an identically zero leading form is correct,
  since it imposes no equation on this chart.
- The sole point missed by that chart is `[u:v]=[0:1]`.  The producer
  evaluates all three fixed-degree forms there and requires at least one
  nonzero value.

Every one of the `11,011` records has infinity gcd `[1]`.  In fact the
nonzero pattern at `[0:1]` is uniformly `(false,false,true)`, so `L_3`
excludes the missing point in every row.  A separate SymPy calculation on a
representative reproduced the chart gcd as `1` and the point values as
`[0,0,23]`.  These two charts exhaust projective infinity.

## 5. Characteristic-zero implication by properness

Fix the selected prime of `Q(lambda)` and its local valuation ring.  The
elements `lambda` and the critical-equation coefficients are integral there:
the only division in the residual polynomial reduction is by `29`, which is
a unit at every selected prime.  Homogenizing to the fixed structural degrees
defines a closed subscheme `X` of projective plane over that valuation ring.

The affine unit calculation and the two infinity checks say exactly that the
special fibre of `X` has no geometric point.  If the characteristic-zero
generic fibre were nonempty, it would have a closed point over a finite field
extension.  Properness of projective space (equivalently, the valuative
criterion after extending the valuation) would specialize that point to a
nonzero projective point of the special fibre.  This contradicts its
emptiness.  Equivalently, the image of the proper morphism `X -> Spec(O)` is
closed: containing the generic point would force it to contain the closed
point.

No affine integrality or leading-coefficient monicity assumption is needed;
a coordinate that is nonintegral simply specializes onto `w=0`, which has
already been checked completely.  Thus no characteristic-zero branch can
escape at infinity.

## 6. Archive-wide invariants and full replay

The canonical archive has `11,011` rows and the first-successful-prime
histogram

```text
101: 10,589
109:    401
131:     21
```

Every row has the correct marked degree record, split-prime embedding,
structural degree bounds, affine unit basis, infinity gcd, nonzero missing
point evaluation, and a 64-hex-character exact-equation digest.

The producer was then run unchanged through its `main(--check)` path with
progress every `500` systems.  It regenerated all rows, serialized the full
archive, and matched SHA-256
`817ae3eff77012e7808ecefedc5fafff368e2631ef3a783229a44221901206cd`.
The successful producer summary was

```text
status = EXACT_FIVE_BLACK_BB_NO_GO
eliminated master systems = 11,011
eliminated marked passports = 55,055
remaining master systems = 22,068
remaining marked passports = 134,992
wrapper status = SCRIPT_MAIN_SUCCESS
```

Elapsed in-process wall time was `915.291` seconds (`915.36` seconds from the
outer timer).  Peak RSS, including final serialization and archive comparison,
was `371,490,816` bytes = `354.281 MiB`, below the requested `0.6 GiB` cap.

The outer `/usr/bin/time -l` wrapper returned status `1` only after the
successful child because its auxiliary `sysctl kern.clockrate` query was
denied by the sandbox.  The unchanged producer completed, printed the exact
success JSON, matched the archive hash, and the wrapper printed
`SCRIPT_MAIN_SUCCESS`; the timer status is not a certificate failure.

## 7. Tests, bytes, and claim boundary

The focused test completed with `6 passed in 0.77s`.  Independent scans of
the quotient and `BB` scripts, archives, reports, and focused test confirmed
valid UTF-8, final LF, valid Python/JSON syntax where applicable, and no NUL,
backspace, vertical-tab, form-feed, or carriage-return bytes.  The quotient
producer's independent `--check` also passed.

The count boundary is exact:

\[
 33079-11011=22068,\qquad 190047-55055=134992.
\]

The preceding exact reductions eliminated `48,343` global marked records,
so the `BB` theorem gives

\[
 48343+55055=103398,
 \qquad 42610678-103398=42507280.
\]

It is correct to mark the five-black `BB` block impossible.  It is also
correct that both `complete_five_black_no_go_claimed` and
`complete_degree_29_no_go_claimed` remain false.

## 8. Auditability note

The archive stores equation digests rather than the full equations or
explicit Singular Bezout combinations.  Archive-only review therefore does
not independently prove each unit-ideal assertion.  This is not a correctness
blocker here because the unchanged producer reconstructed every equation and
every Singular/infinity certificate in the successful bounded replay.

```text
hostile_audit_verdict = PASS
mathematical_blockers = 0
implementation_blockers = 0
five_black_BB_status = IMPOSSIBLE
five_black_status = PARTIAL
global_degree_29_status = OPEN
```
