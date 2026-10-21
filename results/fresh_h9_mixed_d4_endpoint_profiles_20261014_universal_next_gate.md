# The 14-wall chamber has the same gap for every current finite-edge endpoint family

Every finite-endpoint cap on **all 71 open directions from the finite vertices of the 14-wall mixed-D4 chamber** strictly misses the same rational spacelike point used in the [52-edge gate](fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.md). The exact certificate covers **2,221 finite profiles and 92 continuous dihedral tails**, with no angle-label cutoff.

In particular, it covers all 21 open directions at the three vertices introduced by the [mutually orthogonal cap pair](fresh_h9_mixed_d4_endpoint_pair_20261014.md): these contribute 459 finite profiles and 21 tails. Thus those new corners do not provide a finite-endpoint cutter of this point while retaining the 14 walls.

The [producer](../scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_next_gate.py) and [saved certificate](fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_next_gate.json) use the frozen dyadic engine of the preceding gate. They independently reconstruct every physical edge frame from the complete 14-wall cone. For old finite ray (X), outgoing spacelike cone ray (q), and (A=-X^2>0), the second frame vector is

\[
Y=Aq+(X\cdot q)X,\qquad
U=X/\sqrt A,\qquad V=Y/\sqrt{Y^2}.
\]

Exact guards establish the signs, support orthogonality, forward orientation, positive-definite normalized support Gram, its inverse, and every profile cost. The rational point has strictly negative pairing with **all 14** retained walls; its norm is the same strictly positive exact value recorded by the preceding gate.

The parameterization is

\[
r=p+hU+\sqrt{h^2+\rho}\,V,\quad h\ge0.
\]

Each tail is treated on the entire closed interval (\sqrt3\le k\le2), with (p=kp_0\) and (\rho=2-k^2/2\). The exact norm-two Coxeter pairing set is enlarged to the isolated labels 2 through 15 together with the whole ray ( (-\infty,-2\cos(\pi/16)]\). This includes every larger finite label and all parallel and dotted contacts. The proof and interval-error details are those of the frozen [base gate](fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.md).

The exhaustive closed-box cover uses **11,093 boxes**, including **6,703 terminal boxes**, with maximum depth 9. The new-vertex profiles account for 2,930 boxes. Every terminal box either has strictly negative escape pairing or violates one retained wall's enlarged Coxeter pairing set. There are no unresolved boxes or numerical acceptance tests.

A finite collection of such caps admitting a common strict future interior point therefore leaves a strict spacelike gap and has infinite volume. The scope is the complete saved **finite-corner edge queue**. It does not exclude arbitrary caps, affine endpoint profiles, endpoints on subsequently created supports, or changed-core constructions. It also does not classify lines with no finite endpoint.

Validation: full generation and saved replay both pass. The separate [extension review](fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_next_gate_review.md) states which arithmetic is shared with the already reviewed interval engine.
