# Record 17: all eight residual nerves are compact-Coxeter impossible

## Theorem

Under the frozen compact-H9/14 reductions, second-layer record 17 is not the
missing-edge graph shadow of a compact hyperbolic Coxeter polytope.

The exact whole-shadow certificate leaves eight and only eight residual
incidence-complex orbits.  Complete order-four/order-five Lanner placement
joins and every actual rank-nine vertex restriction leave zero label states
on each orbit.  This upgrades the earlier result from a no-go for one selected
ordinary-polytopal nerve to a **whole-shadow compact-Coxeter no-go**.

## Exact input

The input is the frozen certificate
`data/enumerations/h9_14_record17_certificate.json`:

```text
file SHA-256        db613bd1b3f6606f7b140f4183f5fe2897db6ef8292d72670e38523f7c397d94
certificate SHA-256 b98e58b390cbf8d9b40e08da16da29f78a2f04805bb0b13d7f2efb623b9afb51
residual orbits      8
```

For each orbit the verifier extracts the complete minimal-nonface clutter and
reconstructs all nine-element faces directly.  The resulting actual-vertex
counts are

```text
108, 108, 106, 106, 110, 108, 108, 106.
```

The fifth entry is the previously selected `a=0` rational
ordinary-polytopal nerve.  No PL or convex realization is assumed for the
other seven.

## Complete local label search

Every order-four or order-five minimal nonface of a compact Coxeter nerve is
a Lanner diagram.  The search uses the complete classical tables:

```text
order 4: 9 isomorphism types, 108 labelled placements
order 5: 5 isomorphism types, 420 labelled placements
```

At each stage a further forced Lanner support is joined to the current label
states.  Every distinct restriction obtained by intersecting the current
cluster with every actual rank-nine face must be a finite Coxeter diagram.
Within Burcroff's globally exhaustive H9/14 alphabet `m=2,3,4,5`, the finite
test is the complete `A/B/D/E/F/H/I_2(m)` classification, hence an exact
positive-definiteness test rather than floating-point numerics.

The table records `compatible Lanner extensions -> elliptic survivors`.
Support indices refer to the concatenation of the orbit's order-four and
order-five minimal nonfaces in the frozen certificate.

| orbit | `a` | `X` | rank-nine faces | decisive supports | exact stages |
|---:|---:|---|---:|---|---|
| 0 | 1 | `H???FB_` | 108 | `1,2,0` | `420->420; 11634->24; 432->0` |
| 1 | 1 | `H???FB_` | 108 | `2,3,4,0` | `420->420; 1332->96; 384->18; 1944->0` |
| 2 | 2 | `H???F?^` | 106 | `2,3,0` | `420->420; 1332->96; 10368->0` |
| 3 | 2 | `H???Fbo` | 106 | `2,3,4` | `420->420; 1332->96; 384->0` |
| 4 | 0 | `H??EEBo` | 110 | `0,1,2` | `420->420; 1332->96; 384->0` |
| 5 | 1 | `H??EEBw` | 108 | `3,4,5,0` | `420->420; 11634->24; 88->8; 32->0` |
| 6 | 1 | `H??FFB_` | 108 | `4,5,6,0` | `420->420; 1332->96; 384->18; 828->0` |
| 7 | 2 | `H??FF?^` | 106 | `1,2,0` | `108->108; 1542->84; 9072->0` |

Thus five orbits die after three supports.  The remaining three have exactly
18, 8, and 18 states after three supports and die after the fourth.  The
archive freezes the full cluster, restriction-size histogram, raw extension
count, compatible count, survivor count, and SHA-256 of the complete survivor
state set at every stage.

## Why sparse incremental pruning is unrestricted

For an actual vertex restriction, normalize the Coxeter Gram matrix as

\[
G=I-A,
\]

where `A` is symmetric and entrywise nonnegative, with
`A_ij=cos(pi/m_ij)` on a finite non-right edge and zero on a right edge.
Because `A` is symmetric and nonnegative,

\[
G>0 \quad\Longleftrightarrow\quad
\lambda_{\max}(A)=\rho(A)<1.
\]

Changing any unforced right angle to any finite non-right angle gives
`A' >= A` entrywise.  Perron--Frobenius monotonicity gives
`rho(A') >= rho(A)`.  Therefore a sparse restriction with `rho(A) >= 1`
cannot be made elliptic by labels added at a later support or anywhere else
in the nerve.  This proves that every incremental rejection is lossless.
Dotted pairs never occur together in an actual vertex restriction.

Since all eight searches end at zero before continuous dotted Gram entries
are relevant, no global rank-ten Lorentzian Gram search is required.

## Independent hostile replay

`scripts/independent_h9_14_record17_all_coxeter_audit.py` imports neither the
production gate nor the project's Coxeter/Lanner module.  It independently:

1. closes hard-coded classical representatives under all permutations and
   recovers the `108/420` labelled tables and `9/5` isomorphism counts;
2. checks that every row is nonfinite and every vertex-deleted row finite;
3. implements a separate finite Coxeter classifier;
4. reconstructs the eight nerves and all rank-nine faces from the frozen
   clutters; and
5. reproduces every support, restriction histogram, extension count,
   survivor count, and survivor-state hash.

Its frozen summary is

```text
e19c45bb7134690b3649b9476c69829072076a777335c124216aeee5b0e3b609
```

## Scope

This proves a compact theorem for the **entire record-17 graph shadow**, not
merely for the displayed ordinary-polytopal completion.  It does not prove
that the seven other abstract residuals are PL spheres or convex polytopes;
those properties are unnecessary for excluding them as Coxeter nerves.

It does **not** exclude a noncompact finite-volume chamber with different
ideal-face data.  Such a claim would require a completed Lorentzian diagram
and Vinberg's full parabolic zero-neighborhood criterion, not the compact
elliptic-vertex condition used here.

## Artifacts and reproduction

```text
src/coxeter_polytopes/combinatorics/h9_14_record17_all_coxeter_labels.py
scripts/certify_h9_14_record17_all_coxeter_labels.py
scripts/independent_h9_14_record17_all_coxeter_audit.py
data/enumerations/h9_14_record17_all_coxeter_labels.json
tests/test_h9_14_record17_all_coxeter_labels.py
```

```bash
PYTHONPATH=src:. python3 scripts/certify_h9_14_record17_all_coxeter_labels.py --verify
python3 scripts/independent_h9_14_record17_all_coxeter_audit.py
PYTHONPATH=src:. pytest -q tests/test_h9_14_record17_all_coxeter_labels.py
```

Frozen digests:

```text
certificate 59a0ebdd63a28914412cda34d60289e06a539dcf3ae4c8d65debc9e7b2f7d820
archive     592fcdaf863704f97e718bc927e8e0cf051cef5d9c373591010f822499e807cf
hostile     e19c45bb7134690b3649b9476c69829072076a777335c124216aeee5b0e3b609
```
