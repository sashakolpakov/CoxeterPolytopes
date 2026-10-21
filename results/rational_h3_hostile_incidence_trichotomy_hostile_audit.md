# Hostile audit of the rational H3 incidence trichotomy

**Status (2026-08-10): PASS after one substantive repair.**  The audited
memo gives a correct exact reduction of the H3 branch under its declared
complete rational Mori alphabet.  It does not prove an absolute Picard or
Coxeter-dimension bound.  Its affine conclusions remain conditional on
semiampleness, and every branch leaves the length of a negative-definite
exceptional or vertical tail uncontrolled.

The initial version contained one genuine proof-import error in the
hyperbolic branch: it referred to the weak-incidence contraction proof even
though the weak-lane nef divisor is not nef in H3.  Section 5 of the owner
memo now supplies the correct finite-Mori-generator base-point-free argument.
No other numerical or lattice assertion required repair.

## 1. Frozen artifacts

| artifact | SHA-256 |
|---|---|
| `results/rational_h3_hostile_incidence_trichotomy.md` | `c0bae8e9447f0793d21a059fd43c13b23bd4b26ff75d5947925965fcab2e5d44` |
| `scripts/certify_rational_h3_hostile_incidence_trichotomy_hostile_audit.py` | `941c3dddcbd2ad7306531443721394397636344b65b6d371446e18e7bcff5e6d` |
| `data/enumerations/rational_h3_hostile_incidence_trichotomy_hostile_audit.json` | `d5769ceb509f7c962cf8082fc9434b8e8ba4dc127ae965ed6965546c5ad890bf` |
| `tests/test_rational_h3_hostile_incidence_trichotomy_hostile_audit.py` | `d5196b1e82c3c521378c6acb6f2a3db5e57b42cedf628c28d2399550f91bf1fc` |

The certifier and test replay use exact integer/Fraction arithmetic.  The test
does not import the certifier.  The JSON archive explicitly records
`semiampleness: conditional`, `absolute_picard_bound: false`, and
`tail_bound: false`.

## 2. The repaired contraction step

The weak-incidence report contracts the disjoint boundary by applying the
base-point-free theorem with

\[
 D=-\left(K_S+B/2\right).
\]

That proof cannot be copied into H3: for the offending square-minus-one
curve,

\[
                         D.E=1-\frac T2<0.
\]

The contraction itself nevertheless exists projectively under the H3 memo's
complete rational-polyhedral Mori hypothesis.  For an ample Cartier divisor
\(A\), put

\[
 H=A+\sum_i\frac{A.B_i}{4}B_i.
\]

Then \(H\) is nef and big, its null locus is exactly the disjoint union of
the \(B_i\), and it is positive on every other irreducible curve.  Clear
denominators, \(L=mH\).  On the finite list of non-boundary Mori generators
the ratios

\[
             -\,\frac{D.C}{L.C}
\]

have a finite maximum, while \(D.B_i=L.B_i=0\).  Hence an integer
\(a\gg0\) makes \(aL+D\) nef; increasing \(a\) also makes its square
positive, so it is big.  The pair \((S,B/2)\) is klt, and the surface
base-point-free theorem makes \(L\) semiample.  Its connected-fiber
morphism contracts precisely the \(B_i\).

This repair is important: the contraction is valid, but not because H3
inherits W3 nefness.

## 3. Canonical star and the exact profile split

Write \(t_i=B_i.E\), \(T=\sum t_i\), \(Q=\sum t_i^2\), and

\[
                         F=4E+\sum_i t_iB_i.
\]

Using \(E^2=-1\), \(B_i^2=-4\), \(K.E=-1\), \(K.B_i=2\), and
the disjointness hypotheses gives independently

\[
 F.B_i=0,\quad F.E=Q-4,\quad
 F^2=4(Q-4),\quad K.F=-4+2T.
\]

Every other complete Mori generator has nonnegative intersection with
\(F\), so \(F\) is nef exactly in the asserted range \(Q\geq4\), and it is
big for \(Q>4\).

The profile classification needs no bounded search.  For nonnegative
integral \(t_i\),

\[
                   Q-T=\sum_i t_i(t_i-1)\geq0.
\]

Thus \(T\geq3\) implies \(Q\geq3\).  Equality \(Q=3\) forces three entries
equal to one.  If \(Q=4\), the only square partitions are \(2^2\) and four
ones; the former has \(T=2\), so H3 forces four ones.  All remaining H3
profiles have \(Q>4\).

In the standard host, with \(e=E\) and \(h_i=B_i/2\), the Gram matrix has
central and leaf squares \(-2\), leaf pairings zero, and \(e.h_i=t_i\).
The leaf-block Schur complement is

\[
                         \frac{Q-4}{2},
\]

so the spherical/affine/Lorentzian split is exact.  If a relevant \(B_i\)
is twice an integral class, then \(t_i=E.B_i\) is even; hence neither of the
all-one low profiles can use such a leaf.

## 4. Spherical packet

For \((t_i)=(1,1,1)\), order the physical curves as
\((E,B_1,B_2,B_3)\).  Their matrix and its leading principal determinants
are

\[
 Q_{\rm phys}=
 \begin{pmatrix}
 -1&1&1&1\\
 1&-4&0&0\\
 1&0&-4&0\\
 1&0&0&-4
 \end{pmatrix},
 \qquad (-1,3,-8,16).
\]

It is negative definite.  The full-support anti-nef cycle

\[
                         Z=(3,1,1,1)
\]

has intersection vector \((0,-1,-1,-1)\).  No coefficientwise smaller
full-support cycle is anti-nef: central coefficient at most two gives
\(-a+\sum b_i>0\).  Thus \(Z\) is the fundamental cycle.  Directly,

\[
                         Z^2=-3,\qquad K.Z=3,\qquad p_a(Z)=1.
\]

Consequently any normal contraction of this packet is non-rational.  If the
target is additionally \(\mathbf Q\)-Gorenstein, the discrepancy equation
\(Q_{\rm phys}a=(-1,2,2,2)^t\) has the unique solution

\[
                         a=(-2,-1,-1,-1).
\]

The central discrepancy proves non-log-canonicity.  The
\(\mathbf Q\)-Gorenstein qualification is indispensable and is present in
the owner memo.

Blowing down \(E\) changes each boundary square from \(-4\) to \(-3\) and
each pairwise intersection from zero to one.  Because the three original
intersection points on \(E\) are distinct, the image is an ordinary triple
point.  The pullback identity

\[
 \pi^*(K_{S'}+B_1'+B_2'+B_3')
       =K_S+2E+B_1+B_2+B_3
\]

is correct.

### 4.1 Highest root and leaf swap

The host Gram matrix is the negative \(D_4\) Cartan matrix.  Starting at
\(h_1\), reflections in

\[
                         e,\ h_2,\ h_3,\ e
\]

produce \(2e+h_1+h_2+h_3\).  Its square is \(-2\), so the apparent
highest-root reflection already belongs to the facet group.

For \(v=h_i-h_j=(B_i-B_j)/2\), one has \(v^2=-4\) and

\[
                  s_v(x)=x+\frac{(x,v)_M}{2}v.
\]

The half-boundary generators pass the integrality test.  For
\(x\in\operatorname {NS}(S)\),

\[
                         (x,v)_M=x.(B_i-B_j).
\]

Since the Neron--Severi lattice of a smooth rational surface is unimodular,
all these pairings are even exactly when
\(B_i-B_j\in2\operatorname {NS}(S)\).  Conversely, if one pairing is odd,
the reflection correction contains \((B_i-B_j)/4\), which is not in the
standard host \(M\subset\frac12\operatorname {NS}(S)\).  This proves the
memo's exact equivalence.

The condition is genuinely additional, not a formal consequence of the
star.  On \(\mathbf P^1\times\mathbf P^1\), take three distinct fibers of
one ruling, blow up four private points on each, and take a transverse fiber
with one private blowup.  Numerically,

\[
 B_i=f-\sum_{j=1}^4e_{ij},\qquad E=g-e_0.
\]

Then \(B_i^2=-4\), \(E^2=-1\), \(B_i.B_j=0\), and \(E.B_i=1\), but
\(B_1-B_2\) has odd exceptional coefficients.  Pairing it with \(e_{11}\)
gives one.  Thus the local standard host need not contain the leaf-swap
reflection.  This is only a lattice/gluing control, not a claim that this
particular blown-up surface is a full reflection-maximal carrier.

## 5. Affine packet

For four unit incidences, the physical null cycle and host null root are

\[
             F=(4,1,1,1,1),\qquad r=F/2=(2,1,1,1,1),
\]

and direct multiplication by the two Gram matrices gives zero.  Also

\[
                         F^2=0,\qquad K.F=4,\qquad p_a(F)=3.
\]

No step in the hypotheses proves semiampleness.  Conditional on
semiampleness, write \(F\equiv mG\) for the primitive general-fiber class.
Adjunction gives

\[
                         m(2g-2)=4,
\]

whose only positive integral solutions are

\[
                         (m,g)=(1,3),\ (2,2).
\]

After vertical relative minimalization, standard Arakelov nefness for a
relatively minimal genus-\(\geq2\) fibration gives

\[
 0\leq K_f^2
   =K_{S_{\min}}^2+8(g-1)
   =10-\rho(S_{\min})+8(g-1).
\]

Therefore the core bounds are exactly \(26\) for genus three and \(18\) for
genus two.  They are not bounds on \(\rho(S)\): each vertical blowdown
changes \(\rho\) by one, and

\[
                         \rho(S)=\rho(S_{\min})+k
\]

contains the uncontrolled tail length \(k\).  The original \(E\) is
vertical because \(F.E=0\), so even the displayed star contributes at least
one such contraction.

## 6. Hyperbolic packet

After the repaired projective contraction of the disjoint boundary, each
point has type \(1/4(1,1)\) and

\[
 K_S=f^*K_Y-\frac12B,\qquad
 f^*\bar E=E+\sum_i\frac{t_i}{4}B_i.
\]

Equivalently \(f^*(4\bar E)=F\).  Projection formula and the already checked
star identities give

\[
                 \bar E^2=\frac{Q-4}{4}>0,\qquad
                 K_Y.\bar E=\frac T2-1>0.
\]

Since its pullback is nef and big, \(\bar E\) is nef and big.  Its positive
canonical degree proves that \(-K_Y\) is not nef.  Adding an effective
boundary cannot reverse the sign because every effective divisor has
nonnegative intersection with the nef curve class \(\bar E\).

The smooth-carrier coefficient elimination is also exact.  If
\(-\left(K+\Theta\right)\) were nef and \(a_i,a_E\geq0\) were the
coefficients of \(B_i,E\), testing the two kinds of curve gives

\[
 a_i\geq\frac12+\frac{t_i a_E}{4},\qquad
 a_E\geq\sum_i t_i a_i-1.
\]

Hence

\[
             \left(1-\frac Q4\right)a_E\geq\frac T2-1.
\]

For \(Q>4\), the left side is nonpositive and the right side is positive.
Thus no effective boundary, klt or otherwise, makes the same carrier
log-anticanonical nef.

Finally,

\[
 D.F=4-2T,\qquad F^2=4(Q-4)
\]

and Hodge index yields

\[
 10-\rho+b=D^2\leq\frac{(T-2)^2}{Q-4}.
\]

Rearrangement gives a lower bound on \(\rho\), never the upper bound needed
for a dimension theorem.

## 7. What survives the audit

The three outputs are therefore:

1. \(Q=3\): a non-log-canonical spherical packet whose highest-root
   reflection is already present and whose leaf-swap extension needs
   independent mod-two gluing;
2. \(Q=4\): a nef isotropic packet with conditional genus-two/three core
   bounds, but no semiampleness theorem and no vertical-tail bound; and
3. \(Q>4\): a nef-big marked quotient curve of positive canonical degree,
   for which both effective-boundary repair and Hodge index point in the
   wrong direction.

The local tail warning is real.  A chain of length \(\ell\) with squares

\[
                         (-2,\ldots,-2,-1)
\]

and adjacent intersection one has leading determinants
\((-1)^j(j+1)\) for \(j<\ell\) and full determinant
\((-1)^\ell\).  Hence arbitrarily long negative-definite unimodular finite
tails exist at the level of local intersection theory and can lie
orthogonal to a pulled-back nef class.  This does not construct an
unbounded family satisfying finite volume and reflection maximality.

Accordingly the owner memo's missing theorem is stated at the right level:
finite-volume full-nef geometry plus reflection maximality must bound, or
equivariantly eliminate, the exceptional/vertical tail.  Another local
discrepancy or Chern computation will not supply an absolute rank bound.

## 8. Terminal replay

With the stable project interpreter:

```text
$ /Users/sasha/.venvs/global/bin/python \
    scripts/certify_rational_h3_hostile_incidence_trichotomy_hostile_audit.py \
    | diff -u \
      data/enumerations/rational_h3_hostile_incidence_trichotomy_hostile_audit.json -
[no output]

$ /Users/sasha/.venvs/global/bin/python -m pytest -q \
    tests/test_rational_h3_hostile_incidence_trichotomy_hostile_audit.py
......                                                                   [100%]
6 passed
```
