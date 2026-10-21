# Hostile audit of the canonically colored affine cusp-factor no-go

## Verdict

**PASS after two exact repairs.**  The marked-packet theorem is correct: the
canonical colors impose one weighted difference equation, while primitive
physical \(\widetilde D_n\) packets satisfying it exist for every even
\(n\ge4\).  The coordinate mirror really bisects the split affine chamber,
and its integrality on the full adapted host is exactly an ambient mod-two
condition.

The repairs were:

1. \(\widetilde D_4\) has one four-branch exceptional cluster, not two
   three-branch end clusters; its proximity equation is
   \(1+1+1+1=4\).
2. the adapted-host warning was sharpened to the equivalence

   \[
   q_M(M,v)\subset2\mathbf Z
   \quad\Longleftrightarrow\quad
   B_n-B_{n-1}\in2\operatorname {NS}(S).
   \]

No global surface realization, reflection-maximal counterfamily, or
dimension theorem follows.  The result is a local obstruction to deriving
a rank bound from marks and canonical colors alone.

## 1. The weighted canonical-color identity

Put

\[
 M=\langle\operatorname {NS}(S),B_\alpha/2\rangle,
 \qquad q_M=2q_{\operatorname {NS}(S)}.
\]

The three physical wall roots have the exact colors

\[
 q_M(K,B/2)=2,\qquad q_M(K,R)=0,\qquad q_M(K,E)=-2.
\]

For an irreducible affine factor at the cusp \([F]\), let

\[
 u=\sum_i a_ir_i=\lambda F/2
\]

be its primitive positive integral null relation.  The quotient
\(M/\operatorname {NS}(S)\) has exponent two, hence

\[
 D=2u=\sum_Ba_iB_i+2\sum_Ra_iR_i+2\sum_Ea_iE_i
     =\lambda F
\]

lies in \(\operatorname {NS}(S)\).  Since \(F\) is primitive, a
\(\mathbf Z\)-linear functional on \(\operatorname {NS}(S)\) takes \(F\)
to one; applying it to \(D=\lambda F\) proves
\(\lambda\in\mathbf Z_{>0}\).  Pairing the null relation with \(K\) then
gives, without a hidden factor of two,

\[
 2\sum_Ba_i-2\sum_Ea_i=4\lambda,
 \qquad \sum_Ba_i-\sum_Ea_i=2\lambda.
\]

This controls a difference of marked masses, not either mass or the number
of vertices.

## 2. Independent reconstruction of the even-rank packets

Take the vectors

\[
 \alpha_0=-e_1-e_2,\quad
 \alpha_i=e_i-e_{i+1}\ (1\le i<n),\quad
 \alpha_n=e_{n-1}+e_n.
\]

The negative Euclidean Gram matrix has diagonal \(-2\), the displayed
affine \(D_n\) edges, rank \(n\), and radical vector

\[
 (1,1,2,\ldots,2,1,1).
\]

For even \(n\), color the four leaves and odd internal vertices by \(B/2\),
and the even internal vertices by \(E\).  Every edge is bipartite.  Replacing
the roots by physical curves therefore gives diagonal squares \(-4\) and
\(-1\), edge intersection one, and no other intersection.  The positive
and negative mark masses are exactly

\[
                         \sum_Ba_i=n,
 \qquad                  \sum_Ea_i=n-2.
\]

The physical null divisor has coefficients one on leaf \(B\)'s, two on
internal \(B\)'s, and four on every \(E\).  At each component its
intersection is zero, its coefficient gcd is one, and adjunction gives

\[
 F_n^2=0,\qquad K.F_n=2n-2(n-2)=4,
 \qquad p_a(F_n)=3.
\]

The proximity calculation is also exact after the repair.  For \(n\ge6\),
the two end clusters have branch multiplicities \((1,1,2)\), while each
internal cluster has \((2,2)\).  For \(n=4\), the single cluster has
\((1,1,1,1)\).  Every sum is four, the exceptional multiplicity.  This is
ordinary point-blowup compatibility only; it does not manufacture the
downstairs projective branches or a second pencil member.

## 3. The coordinate refinement and the exact glue gate

In

\[
 D_n=\{x\in\mathbf Z^n:\sum x_i\equiv0\pmod2\},
\]

the vector \(v=2e_n=\alpha_n-\alpha_{n-1}\) is primitive of norm four.
Its reflection changes the sign of \(x_n\), preserves \(D_n\), and is not
in the even-sign-change Weyl group.  It interchanges the two terminal
simple roots.  Both have canonical color \(+2\), so \(q_M(K,v)=0\).

The standard affine alcove

\[
 x_1>x_2>\cdots>x_{n-1}>|x_n|,qquad x_1+x_2<1
\]

contains interior points with \(x_n>0\), \(x_n<0\), and \(x_n=0\).
Thus \(v^\perp\) genuinely crosses and bisects its interior; this is not
merely a symmetry of its boundary diagram.

On the physical adapted host,

\[
                         v=(B_n-B_{n-1})/2,
 \qquad s_v(x)=x+q_M(x,v)v/2.
\]

The pairwise-disjoint square-minus-four boundaries make \(M\) an even
integral lattice.  Since \(q_M(v,v)=-4\), \(v\) remains primitive in
\(M\): a nontrivial divisor would have nonintegral or odd square.  Hence
\(s_v\) preserves \(M\) exactly when every \(q_M(x,v)\) is even.

For half-boundary generators those pairings are \(0\) or \(\pm2\).  For
\(x\in\operatorname {NS}(S)\),

\[
 q_M(x,v)=x.(B_n-B_{n-1}).
\]

The Neron--Severi lattice of a smooth rational surface is unimodular, so
this is even for every \(x\) precisely when
\(B_n-B_{n-1}\in2\operatorname {NS}(S)\).  The affine factor does not
record that congruence.  If it holds and the packet occurs in the cusp of
a literal complete reflection chamber, the crossing integral mirror gives
a proper discrete reflection refinement.  If it fails, this particular
coordinate reflection is not integral on \(M\).  Neither outcome alone
settles all possible overgroups.

## 4. Frozen claim boundary

The independent hostile tests reconstruct the Gram matrices, marks,
physical divisors, all proximity clusters through rank 60, explicit points
on both sides of the coordinate mirror, and both outcomes of the exact
adapted-host parity test.  They import no producer function.

The audited theorem proves only:

* arbitrary even affine rank at the marked root-packet level;
* the weighted canonical-color identity and its inability to bound rank;
* ordinary proximity compatibility; and
* the exact conditional coordinate-refinement parity gate.

It does **not** prove arbitrary-rank projective rational genus-three
realization, a full-nef or reflection-maximal family, failure of every
possible maximality argument, or any dimension theorem.

The frozen producer artifacts are:

    183f9bca555e38f5d7ef41255078835636e99b8ce6428756ba6339d8b616710e  results/rational_affine_h3_colored_cusp_factor_no_go.md
    81b4a94588aa2326973031b30df32a00475c5e09e6af3472578afe1d0cfd1b16  scripts/certify_rational_affine_h3_colored_cusp_factor_no_go.py
    2ae937a404f69bffafd2f030b58ea95b70bb7236e433ba70e97e298be79038b0  tests/test_rational_affine_h3_colored_cusp_factor_no_go.py
