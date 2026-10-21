# Hostile audit of the wrong-log orbifold rank gate

**Status (2026-08-10): PASS after two substantive repairs and one proof
completion.**  The audited theorem is correct on its declared input: a
smooth projective rational surface over \(\mathbf C\) whose complete
rational-polyhedral Mori cone has the physical \((-1),(-2),(-4)\) alphabet,
with every \((-4)\)-curve included in the standard half-boundary host and
satisfying (W1)--(W3).  It proves

\[
                         \rho(S)\leq 26,
 \qquad \dim \mathbf P(\operatorname {Nef}S)\leq25.
\]

The theorem does **not** prove raw Coxeter-to-surface existence, equality of
an arbitrary raw polytope with a nef chamber, or a bound in the hostile
incidence locus (H1)--(H3).  Its finiteness corollary is only for abstract
reflection-maximal chambers already carrying the stated literal surface
geometry.

## 1. Audited artifacts

The terminal audit pins the following owner artifacts:

| artifact | SHA-256 |
|---|---|
| `results/rational_wrong_log_orbifold_rank_gate.md` | `36ef3a14968300119f4f0178e9dad11003220d4ee2879cc49c9f0494dc961afa` |
| `scripts/certify_rational_wrong_log_orbifold_rank_gate.py` | `9d5d4657ae5699f6ff044738e094b00b57caab462491fb7835d5c529a9b0c22c` |
| `data/enumerations/rational_wrong_log_orbifold_rank_gate.json` | `62caed0ac6c73021b0c61d4155511028c9d712b003b7b34ae3b3fa8ff3213913` |
| `tests/test_rational_wrong_log_orbifold_rank_gate.py` | `c4abb82716f7b0fd48afa4a2bbea3d566c1f3f8409fe8675d58980aaa107bc77` |

The independent replay is
`tests/test_rational_wrong_log_orbifold_rank_gate_hostile_audit.py`.  It does
not import the production certifier.  It reconstructs the relevant
intersection matrices, continued fraction, discrepancy, orbifold Euler
number, rank elimination, coefficient inequalities, ruled-surface controls,
and \(\mathbb F_4\) calibration from fresh formulas.

## 2. Repairs forced by the hostile audit

### 2.1 H3 had to be aggregate, not componentwise

The first draft treated (H3) as the existence of one component with
\(B_i.E\geq3\).  That is not the negation of (W3), which is

\[
                          B.E=\sum_i B_i.E\leq2.
\]

The omitted pattern was \((B_1.E,B_2.E,B_3.E)=(1,1,1)\).  It is geometric:
three disjoint fibers of one ruling, each made into a \((-4)\)-curve by four
private blowups, meet a transverse fiber once each; one private blowup makes
the transverse fiber a \((-1)\)-curve.

The repaired proof puts

\[
 t_i=B_i.E,\qquad T=\sum_i t_i,\qquad Q=\sum_i t_i^2,
 \qquad c=a_E.
\]

Nefness of \(-(K+\Theta)\) on the \(B_i\) and on \(E\) would imply

\[
 a_i\geq\frac12+\frac{t_ic}{4},\qquad
 c\geq\sum_i t_i a_i-1,
\]

hence

\[
                    \left(1-\frac Q4\right)c\geq\frac T2-1. \tag{2.1}
\]

If \(T\geq3\) and \(Q\geq4\), the left side is nonpositive and the right
side positive.  If \(Q<4\), integrality leaves only \((T,Q)=(3,3)\), and
(2.1) forces \(c\geq2\).  Thus the corrected aggregate H3 is exact and also
catches the distributed \((1,1,1)\) case.

### 2.2 The rank ceiling had to be symbolic

The first certificate inferred the largest Picard rank from a bounded search
with an arbitrary cutoff.  The repaired certificate instead eliminates
\(b\) exactly:

\[
 D^2\geq0\Longrightarrow b\geq\rho-10,
 \qquad
 \widehat c_2(Y)\geq0\Longrightarrow7b\leq4\rho+8,
\]

so

\[
             7(\rho-10)\leq4\rho+8,qquad3\rho\leq78.
\]

The finite replays through ranks \(27,100,500\), and the hostile replay
through rank \(1000\), are now poison checks only.  The unique extremal
integer pair is \((\rho,b)=(26,16)\).

### 2.3 The contraction proof needed its base-point-free hypotheses

The final proof no longer invokes an unnamed contraction assertion.  Given
an ample Cartier divisor \(A\), it uses

\[
 H=A+\sum_i\frac{A.B_i}{4}B_i,
 \qquad H.B_i=0,
 \qquad H^2=A^2+\sum_i\frac{(A.B_i)^2}{4}>0.
\]

Every distinct irreducible curve has positive \(H\)-degree, so
\(\operatorname {Null}(H)=\coprod_iB_i\).  For a Cartier multiple \(L=mH\)
and \(\Delta=B/2\), the pair \((S,\Delta)\) is klt and

\[
                  aL-(K_S+\Delta)=aL+D
\]

is nef and big.  The surface base-point-free theorem makes \(L\) semiample;
the connected-fiber part of its Stein factorization is birational and
contracts exactly the \(B_i\).  This supplies a projective normal target and
also a nef-and-big class on it, as required by the later IJZ input.

## 3. Local singularity and orbifold audit

For a contracted component, every positive exceptional cycle satisfies

\[
 p_a(mB_i)=1-2m^2+m\leq0.
\]

Artin's criterion therefore gives a rational singularity.  Its complete
minimal-resolution string is the one-vertex string \([4]=4/1\), so the
one-string Hirzebruch--Jung theorem identifies it as
\(\frac14(1,1)\).  Solving the discrepancy equation independently gives

\[
                 2=K_S.B_i=(-4)a_i,qquad a_i=-\frac12>-1.
\]

Thus the target is klt and its local uniformizing group has order four.
The relative exceptional matrix is \(-4I_b\), so the exceptional classes
are independent and \(\rho(Y)=\rho(S)-b\).

The orbifold Euler computation was checked in both equivalent forms.  Since
\(e(S)=2+\rho\), deleting the \(b\) exceptional projective lines gives

\[
 e(Y_{\rm reg})=2+\rho-2b.
\]

Adding \(1/4\) for each local quotient chart yields

\[
 \widehat c_2(Y)
 =e(Y_{\rm reg})+\frac b4
 =e(Y)-\frac{3b}{4}
 =2+\rho-\frac{7b}{4}.                                    \tag{3.1}
\]

No ordinary Euler contribution was accidentally counted twice.

## 4. External theorem hypothesis audit

The pinned primary source is Iwai--Jinnouchi--Zhang,
arXiv:2607.18699v1, Theorem 1.3, source label
`thm-nef-anti-canonical`.  The PDF, source tarball, and theorem-containing
TeX member were independently hash-checked against the archive.  In the
source, the theorem assumes:

1. a compact normal complex analytic variety;
2. at most klt singularities;
3. existence of a nef and big class; and
4. nef anticanonical divisor.

The constructed \(Y\) is normal and projective, hence compact and carrying
an ample class; the local discrepancy is \(-1/2\); and

\[
       f^*(-K_Y)=-\left(K_S+\frac12B\right)=D
\]

is nef.  Thus every hypothesis is supplied.  In dimension two the product
of the \(n-2\) auxiliary nef-and-big classes is empty, so the conclusion is
exactly \(\widehat c_2(Y)\geq0\).  No K-semistability assumption from a
different theorem is being imported, and no positive-characteristic claim
is made.

The source pins used in the audit are:

| primary object | SHA-256 |
|---|---|
| PDF | `2739f4489b0a4cf0722141ea1954f871fae65269c8aec63c02f66bfad398c501` |
| source tar | `a18775cdb3f99507da702375c03e983760c894ab9ddd36eaebc9ee96910e8442` |
| `semipositivity_chern_2026_0717.tex` | `6cb03d24e199e4524ed0e2cec4e7bce6371c2aee09f910bde9deb4a297877517` |

## 5. Supported-boundary and geometric controls

The repaired Proposition 5.1 was checked without the production helpers.
For (H1), the standard host forces \(B_i.B_j\) even; at the least value two,
the two curve inequalities already give \(a_i+a_j\geq2\).  For (H2), the
Coxeter-angle exclusion removes intersection one; at intersection two the
two inequalities force \(a_B\geq1\).  Larger intersections only strengthen
the contradictions.  Aggregate (H3) is handled by (2.1).

Fresh ruled-surface calculations reproduce all four hostile controls:

| control | strict-transform squares | retained degree | canonical degrees |
|---|---:|---:|---:|
| two \((1,1)\)'s, six private blowups each | \((-4,-4)\) | \(2\) | \((2,2)\) |
| two \((1,1)\)'s, six and four blowups | \((-4,-2)\) | \(2\) | \((2,0)\) |
| \((1,1)\), \((1,2)\), six and five blowups | \((-4,-1)\) | \(3\) | \((2,-1)\) |
| three boundary fibers and a transverse fiber | \((-4,-1)\) | \((1,1,1)\) | \((2,-1)\) |

These examples certify that the local obstruction types are nonempty.  They
do not assert that the resulting full nef cones are Coxeter, maximal, or
reflection-saturated.

The independent \(\mathbb F_4\) matrix replay also gives

\[
 (K+B/2).(B,E,G)=(0,-1,-1/2),\quad D^2=8,
 \quad\widehat c_2(Y)=13/4.
\]

Thus that maximal triangle lies inside the weak incidence lane and is a
positive calibration, while still refuting the stronger false claim that
maximality forces \(B=-2K\).

## 6. Finiteness scope

The corollary begins only after the literal full-nef carrier and the complete
Mori alphabet are present.  The rank theorem leaves finitely many abstract
rational-surface Neron--Severi lattice types.  For each, a standard host is
an intermediate exponent-two group between the fixed lattice and its
half-lattice, hence there are finitely many possible abstract hosts.
Relative reflection maximality then identifies the chamber with a chamber
of the complete reflection arrangement of that host, whose chambers form one
reflection-group orbit.

This proves finiteness of the **already-carried abstract maximal chambers**.
It does not prove finiteness of surface isomorphism classes, manufacture a
carrier from raw Coxeter data, establish a same-polytope return theorem, or
bound carriers containing (H1), (H2), or aggregate (H3).  Section 7 of the
owner report states exactly this remaining dichotomy: either maximality must
exclude the hostile locus above bounded rank, or a different bounded AG
carrier must handle its explicitly non-klt contractions.

## 7. Terminal replay

With the stable project interpreter:

```text
$ /Users/sasha/.venvs/global/bin/python -m pytest -q \
    tests/test_rational_wrong_log_orbifold_rank_gate.py \
    tests/test_rational_wrong_log_orbifold_rank_gate_hostile_audit.py
..........                                                               [100%]
10 passed
```

The hostile suite contains no production-helper import and includes a
rank-1000 poison replay, the formerly omitted distributed H3 pattern, both
orbifold Euler formulas, the HJ/discrepancy reconstruction, conditional
primary-source hash verification, the four ruled-surface controls, the
\(\mathbb F_4\) matrix, and explicit negative scope gates.

**Verdict:** PASS, carrier-relative.  The wrong-log weak-incidence lane now
has a valid project-native absolute dimension bound.  The result is an
ABSW-type theorem for that declared AG carrier class, not a replacement for
the still-missing raw polygon-to-surface existence theorem and not a theorem
about the aggregate hostile incidence locus.
