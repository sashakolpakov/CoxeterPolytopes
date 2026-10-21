# Hostile audit: raw adjunction/log Stage 0

**Verdict (2026-08-10): terminal PASS; no repair to the theorem or compiler
was required.**  The raw solver really uses only an integral unimodular
Lorentzian Gram matrix and primitive oriented full-span negative wall
vectors.  Its finite adjunction enumeration is exhaustive and independent
of the computational wall basis, and its inferred adapted-module test is an
exact lattice test.  The six advertised calibration outputs all replay.

Audited theorem:
[`raw_adjunction_log_stage0.md`](raw_adjunction_log_stage0.md).

## Pinned audited artifacts

| artifact | SHA-256 |
|---|---|
| `results/raw_adjunction_log_stage0.md` | `6a37b85fc2f762feb4b8b4b89b1d4fc47a85a5f3f0d06682b809b3e692dfe93e` |
| `scripts/certify_raw_adjunction_log_stage0.py` | `a60d4780c9897daf2c60d3bf432c3dc41983449c531ffada1871ba5c118171b5` |
| `data/fixtures/raw_adjunction_log_stage0_input.schema.json` | `f007d407886b729259796d85891bbce49662e4631231e9444cebb48b19db97d9` |
| `data/fixtures/raw_adjunction_log_stage0_inputs.json` | `c197624eb2986e44c2a991f3433cfe81d5b16ca00593282b481844721dcc4113` |
| `data/enumerations/raw_adjunction_log_stage0.json` | `55700e272eaeb527b0bcc913ab401a45898b9d076cf9de50f5159986f6403730` |
| `tests/test_raw_adjunction_log_stage0.py` | `d80c35143395b7474ca07ea421a6ad1103dc4fa010da0940c29fe532efc38348` |
| `tests/test_raw_adjunction_log_stage0_hostile_audit.py` | `db02c3efc87ec14a1f2abd7fc68dbe5ab6ac6b6b1ee2d4a5a131ab0583a0a62c` |

## 1. The finite search is complete

For a physical smooth rational class (C=dr), with primitive supplied wall
(r) and (r^2=-q), adjunction gives

\[
        K.r=qd-2/d.
\]

Because (K.r) and (qd) are integers, (d\mid2).  Characteristic parity
(K.r\equiv r^2\pmod2) then leaves exactly

\[
 (d,K.r)=(1,q-2),
 \qquad (2,2q-1)\text{ only when }q\text{ is odd}.
\]

There is no omitted local scaling.  Any (ho) independent wall vectors
give an invertible system for (K) over (mathbf Q).  Enumerating the
local rows on one such basis and then testing the resulting (K) against
every wall is therefore exhaustive; it does not assume that the selected
basis is an integral basis of the wall span.  The final requirements that
(K\in L), (K) be globally characteristic, and (K^2=10-ho) are all
checked independently after solving.

The weighted greedy choice is sound.  An odd (q) contributes one binary
choice and an even (q) contributes none, so its weight is exactly the
logarithm base two of the number of local rows.  The vector-matroid greedy
algorithm minimizes the sum of these zero-one weights and hence minimizes
the number of assignments.  This optimization affects running time, not the
candidate set.

The hostile suite recomputes every candidate from eight independently
tie-broken minimum-weight bases for each of the five larger fixtures and
from both available bases for the (F_4) fixture.  Every basis returns the
same complete ((K,(d_i))) set.  It also applies three unrelated wall
permutations to every fixture and compares the entire candidate data after
reindexing, not merely the number of branches.  All comparisons pass.

## 2. The inferred module calculation is exact

For inferred boundary components (C_i^2=-2m_i), the code constructs

\[
       M_\Delta=L+\sum_i\mathbf Z(C_i/m_i).
\]

After clearing a common coordinate denominator (D), the columns used by
the compiler generate exactly the integral lattice (DM_\Delta).  SymPy's
column Hermite normal form is therefore a basis of that lattice, and
integrality of the inverse-HNF coordinates is equivalent to membership.
Testing the images of the displayed generators proves
(s_r(M_\Delta)\subseteq M_\Delta).  It proves equality because every
reflection satisfies (s_r^2=1); no separate index comparison is missing.

The hostile test does not reuse HNF.  It enumerates the finite residue group

\[
       M_\Delta/L
       =\left\langle C_i/m_i+L\right\rangle
       \subset (mathbf Q/mathbf Z)^\rho
\]

directly, reflects every (L)-basis vector and every fractional boundary
generator, and tests their residues.  It separately recomputes reflection
integrality on (L) and integrality of the (ell)-scaled form.  This agrees
byte-for-boolean with the compiler for all five accepted branches and for
the rejected (\operatorname{Bl}_p(\mathbb F_4)) log-shaped candidate.

The distinction between the two square-minus-four cases survives this
independent calculation.  For `H10/30`, (C=2r) with (r^2=-1), so
(C/2=r\in L), the inferred module has coordinate denominator one, and all
reflections are already integral on (L).  For `H11/26`, `H14/60`, and
`H17/48`, the physical boundary is a primitive vector of square (-4), the
module has denominator two, and some reflections fail on (L) but preserve
the inferred half-lattice.  The compiler does not collapse these two
geometrically different inputs.

## 3. The six calibrations survive

The exact replay gives:

| case | rank | walls | assignments | numerical (K)'s | accepted |
|---|---:|---:|---:|---:|---:|
| `H10/30` | 11 | 30 | 2 | 2 | 1 Coble branch |
| `H11/26` | 12 | 26 | 1 | 1 | 1 Coble branch |
| `H14/60` | 15 | 60 | 1 | 1 | 1 Coble branch |
| `H17/48` | 18 | 48 | 1 | 1 | 1 Coble branch |
| toric `H6/9` | 7 | 9 | 2 | 2 | 1 weak branch |
| `Bl_p(F_4)` | 3 | 3 | 4 | 1 | 0 |

The accepted physical-square histograms are respectively
`1*(-4)+12*(-2)+17*(-1)`, `2*(-4)+10*(-2)+14*(-1)`,
`5*(-4)+10*(-2)+45*(-1)`, `8*(-4)+24*(-2)+16*(-1)`, and
`6*(-2)+3*(-1)`.  The hostile test independently exhausts the four local
pairing systems for the (F_4) basis, proves that all possible coordinate
solutions lie in its brute-force box, and recovers only
(K=(-2,-6,1)) after the global gates.  The candidate passes boundary
shape, disjointness, and the half-module reflection test but has residual
((-3/2,-6,1)), so rejection at (K+\Delta\ne0) is exact.

No canonical class or boundary marking reaches `solve_raw_input`: its exact
key check permits only `schema`, `picard_gram`, and `wall_vectors`.  The
fixture file stores source canonical classes only under a separate
post-solve `expected` object.  Inspection of the solver confirms that it
does not read case identifiers, source files, or expected outputs.

## 4. Scope that did not pass—and is not claimed

This PASS is numerical.  Pairwise nonnegative wall intersections are a
necessary effective-curve orientation condition, not proof of effectivity.
The solver does not construct a rational surface, prove irreducibility or
smoothness of the inferred curve classes, prove that they generate the full
Mori cone, prove finite hyperbolic volume, or return the same Coxeter
polytope.  Consequently the conditional `9`, `19`, and `21` chamber bounds
do not yet constitute the desired raw-polytope dimension theorem.  No
Nikulin--Vinberg or ABSW/Prasad input is hidden in this audit.

## Replay

```bash
/Users/sasha/.venvs/global/bin/python \
  scripts/certify_raw_adjunction_log_stage0.py \
  --check-fixtures --check-archive

/Users/sasha/.venvs/global/bin/python -m pytest -q \
  tests/test_raw_adjunction_log_stage0.py \
  tests/test_raw_adjunction_log_stage0_hostile_audit.py
```

The first command exactly rebuilds both frozen fixtures and the 30,245-line
archive.  The combined focused suite has 18 passing tests: 11 owner tests
and 7 independent hostile tests.
