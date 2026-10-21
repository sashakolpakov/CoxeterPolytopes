# A genuine polytope-to-K3 theorem: complete 2-root chambers with a Weyl vector

**Status (2026-08-08).** This is an exact, independently recognizable
POLYTOPE → SURFACE → THEOREM ABOUT POLYTOPES route. It is not a new
classification: the mathematical input is Gritsenko--Nikulin, especially
Definition 3.2, Theorem 3.1, and Remarks 3.1 and 3.3 of
[arXiv:1602.08359v2](https://arxiv.org/abs/1602.08359). What is new here is
the strict inverse-transfer formulation, the count/rank audit, and the exact
separation test on Bugaenko's compact H8 chamber.

The route is powerful but narrow. It classifies marked rational integral
chambers for the **full** norm-two reflection arrangement and having a lattice
Weyl vector. It does not classify arbitrary arithmetic Coxeter polytopes and
does not extend from \(\mathbb Z\)-lattices to lattices over totally real
number fields.

## 1. The independently defined polytope class

Use the sign convention of Gritsenko--Nikulin. Let \(S\) be an even integral
lattice of signature \((n,1)\), \(n\geq 2\), and put

\[
R_2(S)=\{r\in S:r^2=2\},\qquad
W^{(2)}(S)=\langle s_r:r\in R_2(S)\rangle .
\]

Choose one component of the negative cone. Let
\(\mathcal M^{(2)}(S)\subset\mathbb H^n\) be a chamber of the complete
\(R_2(S)\)-mirror arrangement, and let
\(P(\mathcal M^{(2)})=\{r_1,\ldots,r_N\}\) be its primitive outward wall
roots. Define \(\mathcal C_{\mathrm{K3W}}\) to be the marked pairs
\((S,\mathcal M^{(2)})\) satisfying:

1. \(\mathcal M^{(2)}\) has finite hyperbolic volume, equivalently
   \(W^{(2)}(S)\) has finite index in \(O(S)\);
2. every facet is a norm-two mirror and the chamber is for the **full**
   group \(W^{(2)}(S)\), not a proper reflection subgroup; and
3. there is a lattice Weyl vector \(\rho\in S\otimes\mathbb Q\) such that

   \[
   (\rho,r_i)=-1=-r_i^2/2\quad(1\leq i\leq N).
   \]

This class is defined without choosing a surface. The ambient integral
lattice is part of the marking. An unmarked real Coxeter polytope is not
enough: the same real chamber can support different finite-index integral
lattices, and the full-root condition depends on which lattice is chosen.

Finite volume implies that the wall roots span \(S\otimes\mathbb Q\), so a
Weyl vector, if it exists, is unique. It is timelike and its ray is the center
of a sphere tangent to every facet. Thus the Weyl-vector condition is already
a strong equal-inradius constraint.

## 2. Exact P-first recognition test

For a marked arithmetic Coxeter chamber \((S,P)\), membership in
\(\mathcal C_{\mathrm{K3W}}\) is decidable in the following order.

1. **Integral carrier.** Verify that \(S\) is even integral of signature
   \((n,1)\), every marked wall root is primitive of square \(2\), and the
   wall roots span \(S_{\mathbb Q}\).
2. **Actual finite volume.** Apply Vinberg's finite-volume criterion, or an
   equivalent exact extreme-ray/elliptic--parabolic face audit. Finite
   sidedness alone is insufficient.
3. **Full 2-root chamber.** Certify, normally by complete Vinberg root
   enumeration, that no other norm-two mirror cuts the open chamber and that
   the supplied walls bound a connected component of the full
   \(R_2(S)\)-arrangement. A chamber for a proper subgroup does not pass.
4. **Weyl equations.** If \(G=((r_i,r_j))\) is the wall Gram matrix, solve

   \[
   Gc=-\mathbf 1,\qquad \rho=\sum_i c_i r_i.
   \]

   Since the walls span, existence is equivalent to
   \(\operatorname{rank}G=\operatorname{rank}[G\mid-\mathbf 1]\). This test
   works directly from the exact marked Gram data.
5. **Classification lookup.** Test the even lattice \(S\), up to integral
   isometry, against Theorem 3.1's finite list. Determinant alone is only a
   preliminary filter.

Passing steps 1--4 already places the datum under the classification theorem;
step 5 identifies its row. Conversely, every row of the theorem comes with
the complete chamber and Weyl vector. This is therefore a total recognition
procedure for the stated class, conditional only on carrying out the standard
exact full-root and finite-volume computations in steps 2--3.

An immediate cheap filter is useful. With all roots normalized to square
\(2\) in a \(\mathbb Z\)-lattice, finite off-diagonal Gram entries are
integral; an acute Coxeter chamber can therefore have finite labels only
\(2\) and \(3\). A \(5\)-edge has product
\(-2\cos(\pi/5)=-\phi\), so it cannot occur in this rational all-2 class.

## 3. The exact classification and its polytope bounds

Gritsenko--Nikulin Theorem 3.1 classifies **59 even lattice isomorphism
classes**. This is a count of integral carriers, not automatically 59
distinct unmarked real Coxeter polytopes: finite-index lattice variants can
share a real root chamber.

The theorem list has the following audited distribution. Compactness is
equivalent to anisotropy of \(S_{\mathbb Q}\), because the full reflection
group has finite index in the arithmetic orthogonal group.

| rank \(S\) | hyperbolic dimension | lattice classes | compact | noncompact finite-volume |
|---:|---:|---:|---:|---:|
| 3 | 2 | 15 | 4 | 11 |
| 4 | 3 | 11 | 1 | 10 |
| 5 | 4 | 7 | 0 | 7 |
| 6 | 5 | 8 | 0 | 8 |
| 7 | 6 | 3 | 0 | 3 |
| 8 | 7 | 6 | 0 | 6 |
| 9 | 8 | 4 | 0 | 4 |
| 10 | 9 | 4 | 0 | 4 |
| 18 | 17 | 1 | 0 | 1 |
| **total** |  | **59** | **5** | **54** |

The compact split is explicit in Theorem 3.1: four anisotropic rank-three
lattices and the anisotropic rank-four lattice \(S_{4,28}\). All listed
lattices of rank at least five are isotropic; this is also forced abstractly
by Meyer's theorem.

> **Exact class bound.** A chamber in \(\mathcal C_{\mathrm{K3W}}\) has
> \(n\leq17\). If it is compact, then \(n\leq3\).

There are no rows in dimensions \(10,\ldots,16\); the last row jumps from
\(\mathbb H^9\) to \(\mathbb H^{17}\). These are necessary and sufficient
bounds for this class, not universal bounds for Coxeter polytopes.

## 4. The actual POLYTOPE → SURFACE transfer

Given a passing pair \((S,\mathcal M^{(2)})\), reverse the sign and set

\[
N=S(-1),\qquad \operatorname{sign}(N)=(1,n).
\]

Remarks 3.1 and 3.3 of the primary source, using the global Torelli theorem,
give a complex projective K3 surface \(X\) with
\(\operatorname{NS}(X)\cong N\). Choose the marking and positive component
compatibly. Then:

* the positive-cone cross-section of the closed nef cone is the
  sign-reversed closed chamber \(\mathcal M^{(2)}\), while its interior is
  the projectivized ample cone;
* the full set \(R_2(S)\) is generally infinite and its reflections tile the
  positive cone;
* the **simple boundary roots** \(P(\mathcal M^{(2)})\), not all roots in
  \(R_2(S)\), become precisely the classes of all irreducible smooth rational
  curves on this marked K3 surface;
* other \((-2)\)-roots lie in Weyl orbits of simple roots and, in the chosen
  nef chamber, need not themselves be irreducible curves; they can represent
  signed or reducible root divisors;
* the finitely many simple curve classes generate \(\operatorname{NS}(X)\)
  up to finite index, and \(\operatorname{Aut}(X)\) is finite; and
* after clearing denominators, \(h=t\rho\) is ample and every smooth rational
  curve \(E_i\) has the same degree

  \[
  h\cdot E_i=t.
  \]

The assertion about simple roots is exactly the content of Remark 3.1:
\(P(\mathcal M^{(2)})\), after sign reversal, gives the classes of all
nonsingular rational curves. It must not be misread as saying that every
norm-two root of the infinite reflection arrangement is a distinct
irreducible curve on the chosen surface.

Conversely, Remark 3.3 says that every complex K3 surface of Picard rank at
least three with finite automorphism group and an ample rational ray giving
the same degree to every smooth rational curve produces one of these 59
lattices and its complete chamber. Thus, within this independently defined
class, one has a genuine equivalence:

\[
\begin{array}{c}
\text{complete finite-volume norm-two chamber}\\
\text{with a rational lattice Weyl vector}
\end{array}
\quad\longleftrightarrow\quad
\begin{array}{c}
\text{finite-Aut K3 Picard lattice}\\
\text{with all smooth rational curves of one ample degree.}
\end{array}
\]

The K3 side is what turns the finite list into a polytope theorem. It is not
merely a surface realization of an already selected example: any marked
polytope satisfying the P-first hypotheses is forced into the classified K3
Picard list, and the class bounds follow back on the polytope.

For a row of rank \(r\), the generic exact-Picard lattice-polarized K3 locus
has dimension \(20-r\). The correspondence is therefore with a moduli
family, not with a canonical individual equation for a surface.

## 5. The high-dimensional positive endpoint

The unique top-rank row is

\[
S_{18,1}=U\oplus 2E_8.
\]

It gives a known finite-volume, noncompact Coxeter polytope in
\(\mathbb H^{17}\). Section 3.4 describes its diagram as

\[
\operatorname{St}(\widetilde E_8,\widetilde E_8).
\]

There is one central wall, two eight-node \(E_8\) bases, and two affine
walls. Hence the chamber has exactly \(1+2(8+1)=19\) facets. It is not
compact because the lattice contains \(U\).

Writing \(\varpi_{E_8}\) for the ordinary \(E_8\) Weyl vector and using the
paper's \(U\)-basis, its lattice Weyl vector is

\[
\rho=(31,-\varpi_{E_8}-\varpi_{E_8},30),\qquad \rho^2=-620.
\]

The transferred K3 has
\(\operatorname{NS}(X)\cong U\oplus2E_8(-1)\), Picard rank 18, a
two-dimensional generic moduli locus, exactly 19 smooth rational curves, and
an integral ample class of square \(620\) meeting every one of them in degree
one.

This is a substantial positive high-dimensional example in the requested
direction, but it is old: neither the abstract \(\mathbb H^{17}\) polytope nor
the Gritsenko--Nikulin K3 interpretation is claimed new here.

## 6. Resolving the paper's rank-19/rank-18 discrepancy

The introduction of arXiv v2 says that the 59 Weyl-vector lattices have
maximal rank 19. That sentence is a typo.

Three independent internal checks settle the point.

1. Theorem 3.1 lists ranks
   \(3,4,5,6,7,8,9,10,18\) with counts
   \(15,11,7,8,3,6,4,4,1\), whose sum is exactly 59. Its last row is
   \(U\oplus2E_8\) of rank 18.
2. The rank-19 lattice \(U\oplus2E_8\oplus A_1\) appears one subsection
   earlier only in the **larger** list of all elliptically 2-reflective
   lattices. It is removed by the Weyl-vector filter. On the star subdiagram
   \(U\oplus K\), the equations force the same \(U\)-coordinate to equal the
   Coxeter number of every irreducible component of \(K\). Here the two
   \(E_8\) components have Coxeter number 30 while \(A_1\) has Coxeter number
   2, an immediate contradiction.
3. The authors' contemporaneous summary
   [Lorentzian Kac--Moody algebras with Weyl groups of 2-reflections](https://ems.press/content/serial-article-files/46626)
   states the corrected maximum, rank 18.

Therefore the sharp hyperbolic-dimension bound supplied by this theorem is
\(18-1=17\), not 18.

## 7. Automorphic corrections: what is proved and what is not

The automorphic layer is useful extra structure, but it is not needed for the
classification or the dimension bound.

For \(T=U(m)\oplus S\) of signature \((n+1,2)\), the paper conjectures for
all 59 rows the existence, for some \(m\), of a holomorphic 2-reflective
automorphic form whose divisor is the complete norm-two Heegner arrangement
with multiplicity one. Such a form gives the Borcherds denominator and the
automorphic correction of the associated Lorentzian Kac--Moody algebra.

The paper's headline proven coverage is **36 of 59**, with more than 70
reflective modular forms considered. Its detailed bookkeeping must be kept
separate from that headline:

* Theorem 4.3 constructs 24 strongly reflective cusp forms
  \(\Phi_{k,K}\). Proposition 4.1 supplies a second correction for the
  already occurring \(U\oplus D_4\) root system. The nine Leech pullbacks in
  Theorem 4.4 are parabolic rather than rows of the elliptic 59-list.
* Section 6 constructs corrections for fourteen elliptic root systems and
  four parabolic systems. Two of those fourteen are the rank-three
  \(\Delta_{11}\) and \(\Delta_5\) cases explicitly identified there as
  previously known.
* Thus simply adding every displayed modular form is not a count of distinct
  new lattice rows. The source does not give a clean 59-row yes/no table.
  Its statement about a planned “rest 10” concerns rank-four and rank-five
  cases and must not be read as the numerical complement \(59-36\).

The certified coverage statement used here is therefore the authors' stated
\(36/59\), together with the explicit fact that existence for all 59 remains
conjectural in that paper. No claim is made here that the other 23 rows all
fail, or that later literature has not resolved individual cases.

Most importantly, the automorphic discriminant belongs to an orthogonal
type-IV moduli space for the extension \(T\), and the K3 surfaces discussed in
Remark 6.7 are arithmetic mirrors with transcendental lattice \(T(-1)\).
This is not automatically a principal divisor on, or an automorphic function
of, the original K3 surface with Picard lattice \(S(-1)\).

Thus automorphic correction is additional coverage on a subset, not a
necessary condition used in the P-first classification and not a converse
theorem for arbitrary arithmetic Coxeter polytopes.

## 8. Exact Bugaenko H8 separation

Bugaenko's compact eleven-wall chamber is over
\(\mathcal O_F=\mathbb Z[\phi]\), \(F=\mathbb Q(\sqrt5)\), not over
\(\mathbb Z\), and it has two \(5\)-edges. It is therefore outside
\(\mathcal C_{\mathrm{K3W}}\) before any classification lookup.

For completeness, the generalized Weyl-vector equations were nevertheless
tested over the full field \(F\). All eleven archived normalized walls have
square two. The first nine wall vectors form an \(F\)-basis, with coordinate
determinant

\[
-3+2\phi\neq0.
\]

Their nine equations force a unique candidate \(\rho_9\). Its pairings with
the eleven walls are

\[
(-1,-1,-1,-1,-1,-1,-1,-1,-1,3-8\phi,3-8\phi).
\]

The last two residuals are therefore

\[
(\rho_9,r_i)+1=4-8\phi,\qquad
N_{F/\mathbb Q}(4-8\phi)=-80\neq0.
\]

Equivalently, the coefficient matrix has rank 9 and the augmented eleven-wall
system has rank 10.

> **Bugaenko obstruction.** The compact Bugaenko H8/11 chamber has no
> vector \(\rho\in L\otimes_{\mathcal O_F}F\) satisfying
> \((\rho,r_i)=-r_i^2/2\) for all eleven walls. In particular it has no
> integral \(\mathcal O_F\)-Weyl vector.

This is a direct field-linear calculation, not an application or extension
of the Gritsenko--Nikulin \(\mathbb Z\)-lattice theorem. It proves that the
most important compact calibration lies genuinely outside this leverage
class. A successful compact high-dimensional theory must relax
equal-norm/equal-inradius, the rational ground field, or both.

Replay:

    PYTHONPATH=src:. python3 scripts/certify_bugaenko_h8_weyl_vector_obstruction.py --check-archive
    PYTHONPATH=src:. pytest -q tests/test_bugaenko_h8_weyl_vector_obstruction.py

Frozen artifacts:

* script: scripts/certify_bugaenko_h8_weyl_vector_obstruction.py;
* certificate: data/enumerations/bugaenko_h8_weyl_vector_obstruction.json;
* tests: tests/test_bugaenko_h8_weyl_vector_obstruction.py.

## 9. What this route earns, and its exact boundary

This route earns all four things demanded of a useful inverse theory:

1. a surface-free, P-first recognition test;
2. an exhaustive K3 carrier classification for the passing class;
3. a sharp transferred finite-volume bound \(n\leq17\) and compact bound
   \(n\leq3\); and
4. positive finite-volume examples through \(\mathbb H^{17}\), together with
   explicit K3 Picard/nef/curve interpretations.

It does **not** imply a bound for any of the following larger classes:

* arithmetic Coxeter polytopes over totally real fields;
* integral chambers with roots of several norms;
* chambers of proper reflection subgroups of \(W^{(2)}(S)\);
* finite-volume chambers without a common lattice Weyl vector; or
* arbitrary surface eigenslices whose full nef cone is larger than the given
  chamber.

Those exclusions are not technicalities. Bugaenko already violates the
ground-field, \(5\)-edge, and Weyl-vector gates. The theorem is therefore a
successful model of P → surface → bound/examples, and at the same time a
precise guide to which hypothesis must next be generalized.
