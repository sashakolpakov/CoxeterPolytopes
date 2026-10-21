# Facet sheets inject into automorphic divisor degree after exact level separation

**Status (2026-08-09):** qualitative Missing Lemma C is proved for the
compact anisotropic determinant-line carrier, after passing to the already
required exact-support level and then to a neat finite cover. The inequality
is exact but is not yet a numerical polytope bound: the automorphic weight and
top automorphic degree still need the uniform upper bounds in Missing Lemma A+.

## 1. Setup

Let \(P\) be a compact marked arithmetic Coxeter polytope with facets
\(F_1,\ldots,F_s\). Apply Proposition 5.5 of
[the main transfer theorem](polytope_to_automorphic_divisor_transfer_theorem.md)
at a direct determinant level where the target special divisor \(Z_i\) has
closed real-slice support exactly \(F_i\). Pass, if necessary, to a finite
neat subgroup contained in that discriminant kernel. Exact support persists
under this further pullback.

Assume the orthogonal quotient \(X\) is compact (in particular, the
anisotropic totally-real case) and let \(\lambda=\mathcal M_1\) be the
weight-one automorphic line bundle on the neat quotient. Bruinier constructs
\(\mathcal M_1\) and gives it a Petersson metric with first Chern form the
invariant Kähler form \(\Omega\). Thus \(\lambda\) is positive and hence
ample by Kodaira. It is a genuine line bundle on the neat cover, not merely
an orbifold class.

Let a holomorphic automorphic section have effective, nonzero divisor

\[
 D=\sum_Y a_Y Y,\qquad a_Y\in\mathbf Z_{>0}.
\tag{1.1}
\]

Assume every target \(Z_i\) occurs in \(D\) positively. Let \(q_i\) be the
minimum multiplicity in \(D\) of a component of \(Z_i\) containing the
marked sheet \(F_i\), and put \(q_{\min}=\min_i q_i\). In the positive
compensated Borcherds carrier these are positive integers after clearing
denominators and killing the finite multiplier by a power.

## 2. The sheet-to-degree theorem

### Theorem 2.1

The marked facets inject into distinct irreducible components of \(D\), and

\[
 \boxed{
 s q_{\min}\leq \deg_\lambda D
 :=D\cdot c_1(\lambda)^{n-1}.
 }
\tag{2.1}
\]

If the section has integral weight \(w>0\) and finite character, then

\[
 \deg_\lambda D
 =w\,c_1(\lambda)^n
 =w\,\operatorname{vol}_\Omega(X)
\tag{2.2}
\]

in Bruinier's normalization. Consequently

\[
 \boxed{
 \#\operatorname{Facets}(P)
 \leq \frac{w\,c_1(\lambda)^n}{q_{\min}}.
 }
\tag{2.3}
\]

For rational weight or a finite multiplier, take a common positive power.
Both weight and divisor multiplicities scale by that power, so (2.3) is
unchanged after the corresponding normalization.

#### Proof

For each \(i\), choose an irreducible component \(Y_i\subset Z_i\) whose
upstairs mirror contains \(F_i\). Suppose \(Y_i=Y_j\) for \(i\ne j\). The
inverse image of this quotient component is one orbit of mirror components
under the neat group. Hence the support of \(Z_i\), pulled back to the
canonical real slice, contains \(F_j\) as well as \(F_i\). This contradicts
exact closed-slice support, which says \(Z_i\cap P=F_i\). Thus the \(Y_i\)
are pairwise distinct. This is a statement about quotient components after
level separation, not raw vector representatives or divisor classes.

Because \(\lambda\) is an ample integral line bundle and \(Y_i\) is a
nonzero integral \((n-1)\)-cycle,

\[
 \deg_\lambda(Y_i)
 =Y_i\cdot c_1(\lambda)^{n-1}\in\mathbf Z_{>0}.
\tag{2.4}
\]

Additivity gives

\[
 \deg_\lambda D
 =\sum_Ya_Y\deg_\lambda(Y)
 \geq\sum_{i=1}^s q_i\deg_\lambda(Y_i)
 \geq s q_{\min}.
\tag{2.5}
\]

An automorphic form of weight \(w\) and finite character is a section of
\(\lambda^{\otimes w}\otimes T\), with \(T\) torsion. Hence
\([D]=w c_1(\lambda)\) rationally. Intersecting with
\(c_1(\lambda)^{n-1}\) gives the first equality in (2.2); Bruinier's degree
formula gives the second. This proves (2.1)--(2.3). \(\square\)

## 3. What this closes, and what it does not

After exact level separation, every marked real facet consumes at least one
unit of global automorphic degree, or \(q_{\min}\) units with multiplicities.
No Kudla intersection formula, Hirzebruch proportionality, or real-facet
volume estimate is needed for this injection.

This is not yet a numerical bound intrinsic to the original polytope class.
The right side changes on the level/neat cover: \(c_1(\lambda)^n\) scales by
cover degree, and the construction does not yet bound the required level,
weight \(w\), or cover degree. Those are the remaining quantitative tasks in
Missing Lemma A+. Once

\[
 w\,c_1(\lambda)^n\leq U(P)
\tag{3.1}
\]

is bounded intrinsically, (2.3) immediately returns
\(s\leq U(P)/q_{\min}\).

The level hypothesis is essential. On a coarser quotient, several chamber
walls can be sheets of one arithmetic Heegner component, so facets need not
inject into components. Redoubling makes the failure structural: larger
Coxeter chambers inside one arithmetic mirror arrangement can have
unboundedly many walls while a fixed coarse quotient retains its fixed mirror
orbits. The separating cover records this growth through its index and hence
through \(c_1(\lambda)^n\). A cover-independent inequality using only the
unseparated coarse divisor degree is therefore false in this generality.

For \(F=\mathbf Q\), the same component count works on a projective
compactification only after controlling boundary divisors and an ample
extension of a power of \(\lambda\). No such boundary version is asserted
here. The theorem is the compact anisotropic version needed for the principal
high-dimensional totally-real route.

## 4. Primary sources

1. J. H. Bruinier,
   [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076),
   J. Reine Angew. Math. **672** (2012), 177--222. Section 2 constructs
   \(\mathcal M_1\), identifies its Petersson curvature with \(\Omega\), and
   records projectivity in the anisotropic case; the degree remark in
   Section 4 gives
   \(\deg c_1(\mathcal M_w)=w\operatorname{vol}(X)\); Theorem 6.8 gives the
   automorphic section and its special-divisor cycle.
2. K. Kodaira,
   [*On Kähler varieties of restricted type (an intrinsic characterization
   of algebraic varieties)*](https://doi.org/10.2307/1969701),
   Ann. of Math. **60** (1954), 28--48.
3. D. Allcock,
   [*Infinitely many hyperbolic Coxeter groups through dimension 19*](https://arxiv.org/abs/0903.0138),
   Geom. Topol. **10** (2006), 737--758. Its doubling/redoubling
   construction supplies the coarse-quotient warning: finite unions of
   adjacent Coxeter chambers acquire growing boundary sets inside one
   underlying mirror arrangement. This warning is not used in the proof of
   Theorem 2.1.

## 5. Novelty and lineage

Positivity of \(\lambda\), integrality of ample intersection degrees, and
Bruinier's weight--degree identity are standard. The project observation is
that the already proved exact closed-slice level separation converts them
into an injection from marked facets to irreducible divisor components.
Thus qualitative Missing Lemma C is not an independent arithmetic-volume
theorem; its quantitative content reduces to bounding the carrier weight and
the automorphic volume of the required cover.
