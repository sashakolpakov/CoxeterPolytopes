# Hostile audit: the dominant-correspondence field-height gate

**Verdict (2026-08-11): PASS after scope repairs; two were required.**  A uniform
bound on one mixed intersection per member of a bounded field-generating
pair-spectral packet gives a finite field set.  The theorem does not derive
that height bound from reflection maximality.

The two repairs were genuine but did not weaken the theorem's formal
hypotheses:

1. preservation of an order-module does not force the individual Cartan
   coefficients of nonprimitive displayed normals into the order.  The
   cyclic product is nevertheless integral by the invariant trace identity
   \(\operatorname {tr}_k(s_is_j)=n-3+\xi_{ij}\);
2. a bounded pair packet bounds the ambient field only when it generates
   that field.  The opening and final conclusion now say
   “field-generating” literally.  A duplicated `Assume:` marker was also
   removed.

## 1. Signed-place inequality

For the orthogonal trace decomposition, write

\[
 H=\sum_\sigma H_\sigma,\qquad
 w_\sigma=q(H_\sigma,H_\sigma).
\]

All nondistinguished \(w_\sigma\le0\).  If
\(0<\sigma(\theta)\le\sigma_0(\theta)\), then

\[
 H.\tau_\theta(H)-\sigma_0(\theta)H^2
 =\sum_{\sigma\ne\sigma_0}
   (\sigma(\theta)-\sigma_0(\theta))w_\sigma\ge0.
\]

Thus \(H^2\ge s_0\) and \(H.\tau_\theta(H)\le C\) imply that every
conjugate of a dominant totally positive \(\theta\) lies in
\((0,C/s_0]\).  This is signed Hodge linear algebra, not an operator-norm
bound.

The identity applies to the full trace space.  If the trace lattice is a
proper Néron--Severi summand, the positive class must itself lie in that
summand, or the complement must receive a separately certified signed-place
decomposition.  The K3 specialization avoids this issue by assuming
\(\operatorname {NS}(X)=N\).

## 2. Primitive positive combination

For \(d\) embeddings there are \(d(d-1)/2\) unordered embedding pairs.
Each pair cuts out one proper coefficient hyperplane because the joint
generators separate embeddings.  A proper hyperplane contains at most
\(M^{r-1}\) points of \(\{1,\ldots,M\}^r\).  With

\[
 M_d={d(d-1)\over2}+1,
\]

the union contains fewer than \(M_d^r\) grid points.  A positive integral
combination outside it has \(d\) distinct conjugates and is therefore a
primitive generator.  If each input coordinate has all conjugates in
\((0,R_0]\), the combination has all conjugates at most
\(R=rM_dR_0\), and the general discriminant bound is
\(R^{d(d-1)}\).

### Dominant generators are not arbitrary generators

The theorem therefore requires the field to be pair-generated; it never
replaces arbitrary generators by dominant elements.

## 3. Cartan integrality and zero conjugates

At every nondistinguished definite place,

\[
 \sigma(\xi_{ij})
 ={4h(e_i,e_j)^2\over h(e_i,e_i)h(e_j,e_j)}\in[0,4]
\]

by Cauchy--Schwarz.  At the distinguished Lorentzian place it is at most
four for intersecting/asymptotic walls and greater than four for
ultraparallel walls.  A nonzero field element cannot vanish under an
embedding, so nonzero pair coordinates are totally positive.  The
ultraparallel coordinates are dominant; the remaining coordinates already
have the uniform bound four.

For integrality, write a reflection as \(s_i=I-e_i\alpha_i\), where
\(\alpha_i(e_i)=2\).  Then on the \((n+1)\)-dimensional field space

\[
 \operatorname {tr}_k(s_is_j)
 =(n+1)-4+\alpha_i(e_j)\alpha_j(e_i)
 =n-3+\xi_{ij}.
\]

Since the product preserves the returned order-module, its trace is in the
order and so is \(\xi_{ij}\).  This proof is scale invariant and remains
valid for nonprimitive displayed normals; it imports no integrality from an
abstract Coxeter diagram.

## 4. K3 rank, generator count, and discriminant

Exact full-Picard trace gives

\[
 d(n+1)=\rho(X)\le20.
\]

For \(n\ge2\), \(d\le6\).  Starting from \(\mathbf Q\), each selected pair
coordinate outside the current subfield enlarges its degree by at least two.
Thus \(2^r\le d\) and a pair-generating packet contains a generating
subpacket with \(r\le\lfloor\log_2d\rfloor\le2\).  The case \(d=1\) is
trivially \(\mathbf Q\).

Uniformly using \(D=6\), \(M_D=16\), and
\(R_0=\max\{4,C/s_0\}\), the primitive element has conjugates at most
\(R=32R_0\).  Its Vandermonde contains \(d(d-1)/2\) differences; after
squaring and enlarging \(d(d-1)\le30\),

\[
 |\operatorname {Disc}k|
 \le|\operatorname {disc}(\theta)|
 \le(32R_0)^{30}.
\]

The elementary symmetric bounds also put every coefficient of the monic
minimal polynomial in a finite integer box.  Therefore this is a finite
field set, not merely a degree bound or discriminant slogan.

## 5. Geometric height and remaining implication

For a returned correspondence,

\[
 H.\tau_\theta(H)
 =\deg\bigl(Z_\theta.p_1^*H.p_2^*H\bigr).
\]

The theorem assumes a uniform upper bound for this mixed bidegree on a
field-generating spectral packet.  An effective bounded-degree return, or a
signed grammar bounding both positive and negative degrees, supplies it.
Fixed \(H^2\) alone does not: in the exact rank-four K3 family,
\(H^2=14\) while \(H\cdot\tau_t(H)=32a-10\) is unbounded.

The still-open implication is precisely

\[
 \text{reflection maximality plus exact eigen-nef return}
 \Longrightarrow
 \text{bounded canonical mixed bidegree}.
\]

No such implication is claimed: the theorem does not derive the scalar
height bound from reflection maximality.  Unbounded mixed degree, non-pair
cyclic generation, or failure to return a positive trace-block class remain
genuine escape routes.

## 6. Frozen owner bytes

| artifact | SHA-256 |
|---|---|
| theorem report | `3c1f779bbb6a765c368ff928c41287ddb221b1c3c56b60a4f8c920dec39c01b8` |
| certifier | `552b598393e4f5f17616d5a04e75bb841d2cf838c5175e672d96b59fdd6b1ee6` |
| archive | `938496ca23ad8da741d3bfa7bff9735f3a3d37e97492d3da1c045003d2b9e423` |
| owner tests | `237389fbc61c40ebaa378125e38de34bd35f65522c69b5ed23667fa18628ae15` |

The independent hostile replay is
[here](../tests/test_k3_dominant_correspondence_height_field_bound_hostile_scope.py).
