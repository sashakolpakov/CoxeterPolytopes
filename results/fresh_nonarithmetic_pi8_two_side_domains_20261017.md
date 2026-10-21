# Continuous endpoints for both new r8 edges

**Generation and fresh saved replay PASS.** This package gives complete continuous endpoint domains for the two open edges supported by \(S_6+[7,r_8]\) and \(S_6+[30,r_8]\) in the58-wall chamber. It does not claim that every parameter satisfies the remaining Coxeter contacts or that two compatible closing sides have already been selected.

The [source](../scripts/fresh_nonarithmetic_pi8_two_side_domains_20261017.py) and [data](fresh_nonarithmetic_pi8_two_side_domains_20261017.json) bind the frozen [compact tetrahedral-face construction](fresh_nonarithmetic_pi8_tetrahedron_completion_20261017.md) and [all-label triangle-facet obstruction](fresh_nonarithmetic_pi8_triangle_facet_gate_20261017.md). Arithmetic and finite-profile enumeration are shared with those certified dependencies; this is a new domain construction, not their independent review.

The common current finite endpoint is right to \([0,1,2,3,5,6,7,30,55]\). Its future timelike norm and all58 wall inequalities are checked explicitly. At either new endpoint, a further wall must be right to H4. The remaining supports are \(A_3+A_1\) and \(A_2+A_1+A_1\). Their exact inverse-cost tests give16 and36 finite contact profiles, respectively, and one and two isolated-A1 high-label tails. These domains include nonright A2 and r8 contacts. Every label at least6 is retained through its full parameter \(c=2\cos(\pi/m)\in[\sqrt3,2)\), with no upper cutoff.

For either support, write the contact magnitudes as \(c_i\), the positive-space duals as \(d_i\), and
\[
p=-\sum_i c_i d_i,\qquad \rho=2-p^2>0.
\]
A certified future unit vector \(k_0\) and positive unit vector \(k_1\) span its Lorentz2 complement. For either sign \(s\), put
\[
A=(k_0+s k_1)/\sqrt2,\quad B=(-k_0+s k_1)/\sqrt2.
\]
Every compatible norm2 direction with that prescribed contact row has the form
\[
r(t)=p+\sqrt{\rho/2}(tA+B/t),\qquad t>0.
\]
Both orientations are stored. The corresponding future finite intersection point is
\[
x(t)=tA-B/t,\qquad x(t)^2=-2,\quad x(t)\cdot r(t)=0.
\]
For every one of the58 old walls, the data retain the full projection offset and both null coefficients. Thus its actual-endpoint inequality is exactly \(a_i t^2-b_i\le0\), and its new-wall product is
\[
p\cdot n_i+\sqrt{\rho/2}(a_i t+b_i/t).
\]
The latter must still be zero, a finite Coxeter product, or at most−2. In particular, the equations for walls49,50,51 and the two tetrahedron walls56,57 have not been omitted or replaced by a sampling rule.

The same persistent finite point \(V\) is included as a common-interior certificate: if a finite collection of proposed acute additions is strict at \(V\), one sufficiently small common perturbation toward the old strict center remains strict for the whole collection. Mutual acuteness and mutual Coxeter contacts are separate requirements.

Generation29628 and fresh saved replay19752 exited0. Syntax, direct source hashes, and report links passed before freezing. Numerical finite-event searches using these domains are separate discovery work and are not completeness certificates.
