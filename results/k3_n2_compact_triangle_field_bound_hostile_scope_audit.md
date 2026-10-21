# Hostile audit: the finite-area K3 Coxeter-triangle field bound

**Verdict (2026-08-11): PASS in finite-area cyclic-ground-field scope,
including ideal vertices.** Exact full trace gives
\(3d=\rho(X)\le 20\), hence \(d\le 6\). Every two facets of a finite-area
triangle meet at a finite vertex or are asymptotic at an ideal vertex, so
there is no ultraparallel wall coordinate.

Put \(T_i=1-s_i\) on the finite projective order-module. The canonical
module traces

\[
 \operatorname{Tr}(T_iT_j)=p_{ij},\qquad
 \operatorname{Tr}(T_0T_1T_2)=c_3
\]

put every pair product and the directed triangle product in
\(\mathcal O\), without requiring primitive normals or integral individual
directed Cartan entries. At every real embedding,

\[
 0\le \sigma(p_{ij})\le 4,\qquad
 \lvert\sigma(c_3)\rvert\le 8.
\]

These are deliberately weak endpoint bounds: an ideal distinguished
vertex may give \(p_{ij}=4\). Right-angle zeros simply delete the
corresponding coordinate. Every closed directed walk splits recursively at
a repeated vertex into simple cycles; on three vertices those are constants,
pair cycles, or the triangle cycle and its equal reverse. Thus the bounded
simple packet generates the full cyclic Cartan ground field stipulated by
the carrier.

A minimal field-generating subpacket has at most two elements because
\(d\le 6\). Among the coefficients \(0,\ldots,15\), some linear
combination is primitive. Its conjugates have magnitude at most
\(16\cdot 8=128\), so

\[
 \lvert\operatorname{Disc}(k)\rvert\le 256^{30}.
\]

Hermite--Minkowski therefore gives a literal finite field set, not only a
degree bound.

Complete saturation and declared-host reflection maximality are retained
but unused in the height argument. The theorem covers compact triangles
and ideal-vertex finite-area triangles. It does not cover truncated
ultraideal polyhedra, gratuitous extensions beyond the cyclic ground
field, carrier construction, or arbitrary \(n=2\) polygon incidence. No
reflective-lattice classification enters.
