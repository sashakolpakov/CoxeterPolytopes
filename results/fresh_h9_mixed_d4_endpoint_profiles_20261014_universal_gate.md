# Every finite endpoint cap on the 52 old mixed-D4 edges misses one strict gap

The complete **1,825 finite profiles and 72 continuous dihedral tails** on the saved 52 outgoing edges of the 12-wall mixed-D4 chamber cannot cut the rational spacelike point recorded in the [certificate](fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.json). This is an exact all-parameter, all-label statement for those finite-endpoint families. It does not exclude caps without such an endpoint, caps based on later edges, affine endpoint profiles, or changed cores.

The [producer](../scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py) uses integer outward-rounded dyadic intervals at 128 bits. Its complete closed-box cover has 9,941 boxes, 5,919 terminal boxes, and maximum depth 9. All 72 tails are treated continuously, including their limiting parameter. No floating-point acceptance, sampled-label conclusion, or numerical optimizer is used.

The point has rational coordinates in the original ten-dimensional Lorentz basis. Its square is

\[
E^2=\frac{330006865140157191295442899174944207
-147195324722660821894201306369903199\sqrt5}{10^{34}}>0.
\]

It has strictly negative pairing with all 12 old walls and negative pairing with the old future center. These signs are checked by separate nested-quadratic rational arithmetic. The chosen decimal strings are interpreted as exact rationals; the numerical proposal files are not dependencies.

For each saved positive-definite rank-eight support, normalize its roots to square 2 and write its inverse Gram matrix as (J). If the nonright contact magnitudes are (c_i=2\cos(\pi/m_i)), then

\[
p=-\sum_i(Jc)_i s_i,\qquad \rho=2-c^TJc>0.
\]

The frozen [complete local-profile certificate](fresh_h9_mixed_d4_endpoint_profiles_20261014.md) supplies every finite possibility and every isolated (I_2(m)) tail, (m\ge6). Its source manifest is verified. The new checker reconstructs physical normalized Gram entries by exact signs and squares, verifies (BJ=I), positive definiteness, the profile costs, and the oriented residual frames.

Let (U^2=-1,V^2=1,U\cdot V=0), with (U) the old finite start and (V) pointing along the outgoing edge. Every genuine new finite endpoint retaining the old walls uses

\[
r=p+hU+\sqrt{h^2+\rho}\,V,\qquad h>0.
\]

The opposite residual branch meets the backward half-line. Acuteness and a common strict interior exclude a new facet through the old finite start: the inverse Gram of the full spherical corner is entrywise nonnegative. Thus the stated forward branch is necessary; the interval proof additionally includes (h=0).

For a tail, put (k=2\cos(\pi/m)). The checker uses

\[
\sqrt3\le k\le2,\qquad p=kp_0,\qquad \rho=2-k^2/2.
\]

Unlike a parameter divided by (\sqrt\rho), this representation retains the nonzero null-residual limits at (k=2). That endpoint is an enlarged closure case, not an assertion that (m=\infty) gives a finite endpoint.

Every scalar pairing is

\[
f(h)=b+ah+c\sqrt{h^2+\rho}.
\]

For each old nonsupport wall, the checker establishes (a\le0) and (a+c<0), so (f) decreases with (h\ge0). For (E), it also establishes (a<0<c). Hence

\[
r\cdot E\le b+c\sqrt\rho+(a+c)h.
\]

This yields a certified finite integer bound on every possible escape-cutting (h), uniformly over each tail interval. All larger (h) have strictly negative escape pairing.

On that compact domain, each closed dyadic box is rejected either by a strictly negative upper bound for (r\cdot E), or by an old-wall pairing interval disjoint from

\[
\{-2\cos(\pi/m):2\le m\le15\}
\ \cup\ (-\infty,-2\cos(\pi/16)].
\]

This is a relaxation containing **every** acute Coxeter pairing, including all finite labels at least 16, parallel walls, and dotted walls. The integer 16 is a partition boundary, not a label cutoff. Strict rejection comparisons preserve every closed endpoint. Exact Machin arctangent bounds enclose \(\pi\); the degree-128 cosine Taylor polynomial has remainder bounded by \(2^{130}/130!<2^{-140}\). All rounding is outward.

Consequently, any finite coherent collection of caps from these old finite-endpoint families leaves (E) strictly inside every halfspace. A common strict future interior point then gives an open ideal-boundary patch and infinite volume. This conclusion does not cover a construction that first adds other kinds of caps and then uses newly available endpoint supports. In particular, arbitrary all-dotted caps lie outside the finite-endpoint mechanism considered here.

Validation: full generation and saved replay both pass; the saved certificate records deterministic exhaustive-tree digests. The separate [proof and arithmetic review](fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate_review.md) states its independent and shared-engine scopes.
