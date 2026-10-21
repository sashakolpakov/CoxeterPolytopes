# Hostile audit: exact Fekete probe on `[-2,1]`

**Verdict (2026-08-10): PASS.** The frozen certificate proves the stated
degree-six minimality and uniqueness theorem, the exact degree-twenty
supremum, optimality in the explicitly bounded 1,280-row factor grammar, and
the one-exceptional-generator degree bound `23`. It does not promote these
conditional statements to a global Coxeter theorem.

## Pinned artifacts

| artifact | SHA-256 |
|---|---|
| `scripts/certify_fekete_interval_minus2_1.py` | `531e6aa1750e57e2d033285b53aa4f84fec8e8e039ab508272308ab0eee9e144` |
| `data/enumerations/fekete_interval_minus2_1_certificate.json` | `509d16580e98075bd762467a23e9bf61514c6f027582fcbac68d7822ab51579e` |
| `results/fekete_interval_minus2_1_probe_theorem.md` | `bcab709ca69362902fc5bcfc4e9582e9f660663b95fc8f6b1073bbb41c437977` |
| `tests/test_fekete_interval_minus2_1_certificate.py` | `013bdeb48c231c9e1d82162a6de3154cdd8be0dfb1d09ba0d2f35f89d34f2b6e` |

## Independent algebra

Writing `z=(x+1/2)^2` sends `[-2,1]` to `[0,9/4]` and gives

```text
P6=(z-1/4)(z-5/4)(z-9/4).
```

Its only non-root critical points are `5/4 +/- sqrt(3)/3`; their squared
magnitudes are `4/27`. The endpoint magnitudes are `45/64` and `0`, and
`4/27 < (45/64)^2`. Hence the exact supremum is `45/64`.

Any integer polynomial of norm below one vanishes at `-2,-1,0,1`, hence is
divisible in `Z[x]` by `A=x(x-1)(x+1)(x+2)`. At the two roots of
`q=x^2+x-1`, one has `A=-1` and root difference `sqrt(5)`. For a linear
quotient the two strict bounds force its slope and constant to vanish. For a
quadratic quotient, reduction modulo `q` forces it to be an integer multiple
of `q`; since `2(45/64)>1`, only `+P6` and `-P6` survive. This proves both
degree minimality and degree-six uniqueness without a field classification.

For

```text
P20=[x(x+1)]^4[(x-1)(x+2)]^3[x^2+x-1]^3,
```

the exact derivative in `z` factors as

```text
(4z-9)^2(4z-5)^2(4z-3)(4z-1)^3(20z-37)/131072.
```

The nonzero endpoint/critical magnitudes are

```text
91125/1048576,  27/1024,  884736/9765625.
```

Thus the archived supremum `884736/9765625` is exact. This also identifies
the repair behind the interval probe: cubing, rather than squaring, the
factor `x^2+x-1` produces the advertised `0.0905969664`.

## Bounded-family and norm checks

The grammar has exactly `4^4*5=1,280` exponent tuples. For each of the
other 1,279 rows, the producer evaluates at a rational point and checks
exactly

```text
|P(r)|^20 > (884736/9765625)^deg(P).
```

This is a sound one-sided certificate: any such value is a lower bound for
the true supremum. Root isolation is used only to find a useful rational
point; it is not trusted as numerical evidence. I also replayed the family
with an independent critical-point search and exact rational comparisons;
all 1,279 competitors were eliminated and `(4,4,3,3,3)` was the unique
winner. Optimality is asserted only in the frozen factor grammar.

For one distinguished conjugate in `[2,14]`, monotonicity in
`w=x(x+1)` gives

```text
M=210^4*209^3*208^3=159773975758768250880000.
```

Exact cross-multiplication verifies `M*q^22>1` and `M*q^23<1`, where
`q=884736/9765625`. Since a degree-`d` norm contains `d-1` ordinary
conjugates, `d>=24` is impossible and `d<=23` follows.

The report correctly preserves the logical boundary: this bounds
`Q(alpha)`, not an arbitrary overfield; the sharper ambient-field bound
requires a primitive generator, while `t` generating coordinates give only
the compositum bound `23^t`. Integrality, total reality, the conjugate
envelope, bounded `t`, and compact/cusp-aware Gram coverage remain explicit
hypotheses.

## Replay

```text
python3 scripts/certify_fekete_interval_minus2_1.py --check
PASS: degree-6 norm 45/64; unique bounded-family degree-20 winner;
      one-exceptional-generator degree <= 23
real 8.08s; user 7.87s; sys 0.10s

python3 -m pytest -q tests/test_fekete_interval_minus2_1_certificate.py
5 passed in 7.91s
```

Both Python files compile. A control-byte scan of all four frozen artifacts
was empty. No frozen producer, archive, theorem, or test was modified by
this audit.
