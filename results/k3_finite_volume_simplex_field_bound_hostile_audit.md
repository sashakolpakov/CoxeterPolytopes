# Hostile audit: finite-volume K3 Coxeter simplices

**Verdict (2026-08-11): PASS for ordinary finite-volume simplices.**  The
owner theorem proves a finite field set in every dimension and field degree
compatible with a primitive full-trace complex-K3 carrier.  It uses neither
a simplex classification nor a hidden imposed height.  The word
``ordinary'' is essential: truncated hyperideal simplices may have
ultraparallel facet pairs and are outside the proof.

## 1. Rank-two face audit

An ordinary hyperbolic simplex has \(m=n+1\) facets and is combinatorially
a simplex.  Every pair of facets therefore has a nonempty common face in
the finite or ideal closure.  Equivalently, their two reflecting
hyperplanes meet in hyperbolic space or at infinity.  The normalized
two-normal Gram block is positive semidefinite, so for

\[
 p_{ij}={4h(e_i,e_j)^2\over h(e_i,e_i)h(e_j,e_j)}
\]

one has

\[
                              0\le(p_{ij})_0\le4.
\tag{1.1}
\]

At every nondistinguished embedding the form is positive definite, hence
strict Cauchy--Schwarz gives \(0\le\sigma(p_{ij})<4\).  This proves the
claimed all-place pair bound.  No list of simplex diagrams or dihedral
angles is being used.  A hyperideal vertex without its truncation would
give infinite volume; adding truncation facets leaves the ordinary-simplex
scope and can introduce ultraparallel pairs, exactly as the owner report
warns.

## 2. Integral simple-cycle audit

Put \(P_i=1-s_i\).  Since every facet reflection preserves the finite
projective \(\mathcal O\)-module, \(P_i\) is an \(\mathcal O\)-linear
rank-one endomorphism.  With

\[
 a_{ij}={2h(e_i,e_j)\over h(e_i,e_i)},
\]

direct rank-one composition gives

\[
 \operatorname {Tr}_{\mathcal O}(P_{i_1}\cdots P_{i_\ell})
 =a_{i_1i_2}\cdots a_{i_\ell i_1}=c(C)\in\mathcal O.
\tag{2.1}
\]

Thus no division or freeness assumption is hidden in cycle integrality;
the canonical trace of an endomorphism of a finite projective module lies
in \(\mathcal O\).

Symmetry of \(h\) makes each directed cycle equal to its reversal.  Hence

\[
                         c(C)^2=\prod_{(i,j)\in C}p_{ij}.
\tag{2.2}
\]

Equations (1.1)--(2.2) imply

\[
                           |\sigma(c(C))|\le2^{|C|}\le2^m
\tag{2.3}
\]

at every embedding.  Cutting a closed walk at a repeated vertex factors
its scalar Cartan product into products of simple-cycle coordinates.
Consequently, the simple cycles generate every cyclic product and, by the
literal minimal-Vinberg-field hypothesis, generate \(k\).

## 3. Degree and primitive-grid audit

The primitive full trace has rank

\[
                              md=\rho(X)\le20.
\tag{3.1}
\]

For \(d>1\), this gives \(3\le m\le10\) and \(d\le6\).  A minimal
field-generating simple-cycle packet has at most two elements: three strict
adjunctions would force degree at least eight.  Among
\(\theta_1+q\theta_2\), \(0\le q\le15\), one element is primitive, because
the at most \(\binom62=15\) unordered embedding pairs each forbid at most
one coefficient.  Its conjugates have magnitude at most \(2^{m+4}\), so

\[
 |\operatorname {Disc}(k)|
 \le2^{(m+5)d(d-1)}.
\tag{3.2}
\]

Exact maximization subject to \(m\ge3\), \(d\ge2\), and \(md\le20\) gives
the exponent table

\[
 240, 180, 120, 66, 24, 26, 28, 30
 \quad\text{for }m=3,\ldots,10.
\]

The maximum is \(240\), attained at \((m,d)=(3,6)\).  Therefore

\[
                  |\operatorname {Disc}(k)|\le2^{240}=256^{30},
\]

and Hermite--Minkowski proves the finite field set.

## 4. Scope verdict

The theorem assumes the returned field is the minimal cyclic Cartan/Vinberg
ground field; a gratuitous correspondence-field extension is not covered.
Completeness, saturation, and reflection maximality may be retained but are
not needed for the arithmetic bound.  The theorem does not construct a K3
carrier, prove essential surjectivity, treat truncated hyperideal
simplices, or close general higher-facet polytopes.  Within those stated
limits, the proof is complete.
