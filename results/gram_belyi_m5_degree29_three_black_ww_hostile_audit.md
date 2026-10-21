# Hostile audit: degree-29 three-black/`WW` no-go

**Status (2026-08-10): PASS.**

This audit is scoped to the theorem that all `1,372` marked degree-`29`
three-black/`WW` records are impossible.  It does not claim that degree `29`
is globally impossible: all strata with at least four black vertices remain
open.

## Pinned artifacts

```text
bf3eddb5867e43c5710b0fd69dd97bd80a478a213aae6d7a9c18c7fe9fb0b47e  results/gram_belyi_m5_degree29_three_black_ww.md
757e882dfedd822609f3346b6b4362847c261813bdcc3e1afad1cbb2282f4229  scripts/certify_gram_belyi_m5_degree29_three_black_ww.py
7f0a3110c8e05ecabd18314df96f70bdc796465fba791476fbdb6ead7db5234e  data/enumerations/gram_belyi_m5_degree29_three_black_ww.json
2132bd51d534e6f2eb37de93b5aff81a177201b96356ea2e0920dd0953451b6c  tests/test_gram_belyi_m5_degree29_three_black_ww.py
```

## Independent checks

The hostile audit independently verified:

1. The residual logarithmic-derivative quadratic
   `D(x)=29x^2-Ax+B`, both endpoint eliminations, the exceptional `C=K=0`
   branches, all denominator-clearing factors, and the asserted exact
   univariate degrees.
2. Every one of the `588` endpoint-critical modular gcd certificates.
3. The canonical identity
   `q=29 lambda-13`, `q^2+12q+16=0`, the canonical residual quadratic
   `D_t`, and the criterion `R_1(t)=0`, including a double non-root critical
   point.
4. All `196` simple/simple systems with a separate pure-Python verifier.
   It reconstructed all `1,159` irreducible factors in the archive and found
   every extension-field gcd constant.
5. The exact census: `196` marked black-degree records, `784` distinct
   coefficient systems, and `1,372` marked passport records.
6. The local-integrality argument.  Full degree makes the Shabat equation
   locally monic, hence `t` is integral; monicity of `H_t` then makes `y`
   integral.  Reduction therefore cannot discard a characteristic-zero
   solution by sending it to infinity.

No timeout, rational-point-only search, numerical coefficient solve, or
plane-tree census enters the proof.

## Replay

The independent audit obtained:

```text
7 passed
archive --check succeeded
```

The larger focused Gram--Belyi suite subsequently obtained:

```text
19 passed
```

## Claim boundary

The audited conclusion is exactly:

```text
three-black/WW status = IMPOSSIBLE
remaining three-black/WW records = 0
global degree-29 status = OPEN
remaining marked records with at least four black vertices = 42605329
```
