# Independent hostile audit: nondotted field-generating cores

**Verdict (2026-08-11): PASS under the exact generating-core hypothesis.**
The theorem proves a finite set of fields with the safe bound
\(\lvert\operatorname{Disc} k\rvert\le 2^{300}\). It permits dotted edges elsewhere in
the chamber, but not on a pair or chordless perimeter used to generate the
field.

## 1. Integrality on a finite projective order-module

For a returned reflection \(s_i\), put
\(P_i=1-s_i\in\operatorname{End}_{\mathcal O}(M)\). Over \(k\), this is a
rank-one operator \(P_i(x)=\alpha_i(x)e_i\). Hence

\[
 \operatorname{Tr}_{\mathcal O}(P_{i_1}\cdots P_{i_\ell})
 =\alpha_{i_1}(e_{i_2})\cdots\alpha_{i_\ell}(e_{i_1})
 =c(C).
\]

The canonical trace of an endomorphism of a finite projective module lies
in \(\mathcal O\), so every cycle coordinate is an algebraic integer. This
argument does not require a free module, primitive normals, or integrality
of the individual directed Cartan entries. The executable audit also
checks that arbitrary rescaling of the normals changes those directed
entries but leaves every closed product invariant.

## 2. The chordless minor and its exact length

At a definite embedding, the Gram matrix on a chordless \(\ell\)-cycle is
positive semidefinite and has only diagonal and cyclic-perimeter entries.
Take rows \(2,\ldots,\ell-1\) and columns \(3,\ldots,\ell\). This
\((\ell-2)\times(\ell-2)\) minor is triangular with diagonal

\[
                  G_{2,3},G_{3,4},\ldots,G_{\ell-1,\ell}.
\]

Every factor is nonzero because it is a perimeter edge of the nonzero
Cartan graph; injectivity of a field embedding preserves nonvanishing.
Thus

\[
             \ell-2\le\operatorname{rank}G_C\le m,
             \qquad \ell\le m+2.
\]

The independent replay evaluates this exact minor for lengths \(3\) through
\(12\), including the triangular \(3\)-cycle endpoint.

## 3. Uniform height at every embedding

Symmetry of \(h\) identifies a directed cycle with its reverse and gives

\[
                         c(C)^2=\prod_{(i,j)\in C}p_{ij}.
\]

At the distinguished place, every pair on the selected generating
perimeter satisfies \(0\le p_{ij}\le 4\) by hypothesis. At every definite
place, Cauchy--Schwarz gives \(0\le\sigma(p_{ij})<4\). Therefore

\[
                         \lvert\sigma(c(C))\rvert\le 2^\ell\le 2^{m+2},
\]

at every embedding. A selected pair generator satisfies the same safe
final bound directly.

## 4. Primitive grid and exponent \(300\)

The full-trace K3 equality gives \(md=\rho(X)\le 20\), so for \(m\ge 3\) a
nonrational field has \(d\le 6\). A minimal generating subpacket has at most
two coordinates: three strict adjunctions would force degree at least
eight. If two coordinates \(x,y\) are needed, each unordered pair of the
at most six embeddings forbids at most one \(q\) in \(x+qy\). The at most
\(\binom{6}{2}=15\) bad values leave a primitive element for some
\(q\in\{0,\ldots,15\}\).

Its conjugates have magnitude at most \(16\cdot 2^{m+2}=2^{m+6}\), so
conjugate differences are at most \(2^{m+7}\). Vandermonde yields

\[
 \lvert\operatorname{Disc} k\rvert\le 2^{(m+7)d(d-1)}.
\]

Exact maximization under \(m\ge 3\), \(d\ge 2\), and \(md\le 20\) gives
\(300\), uniquely at \((m,d)=(3,6)\). Hermite--Minkowski therefore gives a literal finite
field set.

The scope is narrow and correctly stated.  The theorem assumes that the
minimal Vinberg field admits a generating packet of pair/chordless cycles
whose perimeter edges are all finite or parabolic.  Dotted edges may occur
outside that packet.  A dotted pair generator or a chordless generator with
an ultraparallel perimeter edge remains open, as do carrier construction,
canonical mixed-height bounds, chamber/host finiteness, and essential
surjectivity; neither finite volume nor reflection maximality is used.

The independent executable audit is
[audit_k3_nondotted_field_generating_core_bound_independent_hostile.py](../scripts/audit_k3_nondotted_field_generating_core_bound_independent_hostile.py),
with separate tests in
[test_k3_nondotted_field_generating_core_bound_independent_hostile.py](../tests/test_k3_nondotted_field_generating_core_bound_independent_hostile.py).
