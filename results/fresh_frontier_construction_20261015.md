# A discriminant-31 arithmetic host over Q(sqrt2)

This is a new bounded construction experiment, not a new finite-volume polytope or a nonreflectivity theorem. The complete Vinberg prefix through height 20 has **13 actual walls, three finite vertices, one finite edge, and a strict spacelike escape**. An exact order-two chamber symmetry supplies one additional actual wall of height 
\((193+137\sqrt2)/4\); the same escape remains strict for all 14 walls.

The host is
\[
 O=\mathbb Z[\sqrt2],\qquad
 L=E_8\perp\begin{pmatrix}-2\sqrt2&-1\\-1&2\end{pmatrix}.
\]
Write \(s=\sqrt2\), \(\delta=1+4s\), \(\varepsilon=1+s\), and
\(\beta=\varepsilon\delta=9+5s\). Then
\(\det L=-\delta\), \(N(\delta)=-31\), and \(N(\beta)=31\).
The physical signature is \((9,1)\); the other embedding is positive definite.
The discriminant ideal has prime norm 31. This is a different field from the
previous discriminant-11 and discriminant-55 hosts and is not the rank-23
SilverLambda22 lattice. A literal repository search found no prior instance.
A focused primary-literature search did not identify a theorem excluding this
literal lattice; this is not a complete literature classification. The existing
[rank-three classification over this ring](https://arxiv.org/abs/1512.01133)
does not classify this rank-ten lattice.

## Complete primitive norm alphabet

The ring is norm-Euclidean: rounding both coefficients leaves
\(a,b\in[-1/2,1/2]\), and \(|a^2-2b^2|\le1/2<1\). In particular it is a PID.
Its units are \(\pm\varepsilon^k\). One elementary verification of the last
statement is to reduce a positive unit by powers of \(\varepsilon\) into
\([1,\varepsilon)\). If its norm is 1, its integer coefficient of \(s\) is
\((u-u^{-1})/(2s)\in[0,1)\); if its norm is -1 that coefficient is
\((u+u^{-1})/(2s)\in(0,1)\) for \(1<u<\varepsilon\). Neither gives a nontrivial
unit in that interval. Totally positive units are therefore exactly
\(\varepsilon^{2k}\), all squares.

Let a primitive integral reflection normal have norm \(q=2b\), using the
evenness of the displayed lattice. Reflection integrality says \(Gr\in bO^{10}\).
Multiplication by the adjugate gives \(\det(G)r\in bO^{10}\). Since the coordinates
of a primitive vector generate the unit ideal, Bezout gives \(b\mid\det(G)\).
The prime norm 31 makes \((\delta)=(\beta)\) prime. Moreover \(b\) is totally
positive, by the physical spacelike condition and positive definite conjugate.
After unit rescaling, every reflection direction consequently has a representative
of norm **2 or \(2\beta\)**. Conversely all enumerated vectors are explicitly
checked for reflection integrality. Their primitivity follows from evenness:
a nonunit common factor \(a\) would give \(a^2\mid q/2\), impossible for
\(q/2=1\) or the prime \(\beta\).

This argument concerns all integral lattice reflectors, rather than just a chosen
norm-two subgroup. For nonzero mixed-norm pairing \(g=\beta u\), the Cartan
cycle is \(g^2/\beta=\beta u^2\), with norm at least 31. Conjugate positivity
bounds its conjugate by 4. Thus such a pair cannot meet or be parallel, since
that would bound its physical value by 4 and its norm by 16. Mixed-norm simple
walls are right or strictly dotted.

## Zero sector and bounded enumeration

Coordinates are \((x,y_1,\ldots,y_8,t)\), with basepoint
\(H=(2,0^8,1)\), \(H^2=-2\delta\), and \((r,H)=-\delta x\).
At height zero, write \((y,t)=u+s v\) in the rational even positive lattice
\(E_8\perp A_1\). A norm-two vector satisfies
\(u^2+2v^2=2\), so \(v=0\). The exact rational shell enumeration finds all
240 E8 roots, with the two A1 roots giving precisely 242 zero roots.
A norm-\(2\beta\) zero root would have its E8 covector and \(-2t\) divisible by
\(\beta\); since the E8 inverse is integral and \(\beta\) is coprime to 2, the
whole vector would be divisible by \(\beta\), contradicting primitivity.
The complete initial sector is therefore \(-E_8\) simple roots and the positive
last coordinate root, with diagram \(E_8+A_1\).

Put \(\lambda=E_8y\ge0\), \(l=x-2t\ge0\), and \(J=E_8^{-1}\). The search uses
\[
 2\lambda^T J\lambda+l^2=2q+\delta x^2,
 \qquad h(r)=\frac{\delta x^2}{2q},
 \qquad l\equiv x\pmod{2O}.
\]
For the long norm, every \(\lambda_i\) and \(l\) is also divisible by \(\beta\).
The remaining binary reflection condition follows from these conditions and
\(\beta\sim\delta\); it is also checked directly on every accepted vector.

The finite control boxes are rigorous at both embeddings. Since
\(q<33\), \(q' <4\), \(\delta>6\), and \(-\delta'>4\), height \(T\) gives
\(0<x<\sqrt{11T}\), \(|x'|<2\). Every positive search coordinate has
\(z^2<66(1+T)\), \(|z'|<3\). The physical and conjugate bounds give explicit
integer coefficient boxes, followed by exact filtering. Rational prefix Schur
minima of the literal positive form provide complete pruning at both embeddings.
Prior-wall pruning is valid because all remaining coefficients
\((2y_8,\ldots,2y_1,l)\) are physically nonnegative. No floating-point decision is
used. The ordering is the usual distance ordering of
[Vinberg's construction](https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=3032&wshow=paper).

The run checks 14 exact control levels and 3,723 recursion nodes, accepting four
positive-height walls. The 78 pairs of the 13 walls are 57 right, ten pi/3,
one pi/4, and ten dotted. A strict timelike center and 13 explicit relative
interior facet points certify actuality.

## Complete small cone and finite symmetry

The cone checker tests all \(\binom{13}{9}=715\) possible extremal-ray supports:
70 are dependent, 589 give no feasible ray, and 56 distinct feasible rays remain.
There are three finite vertices, no ideal vertices, and 53 hyperideal rays.
Exactly one pair of finite vertices shares eight facets; the finite-vertex
frontier has 25 remaining open germs. A hyperideal ray plus a small positive
multiple of the strict center gives the saved exact strict spacelike escape.
This proves infinite volume of this finite prefix only.

The third finite vertex has incident walls
\([0,1,2,3,4,5,6,10,11]\), again with the literal norm-colored diagram \(E_8+A_1\).
An exact integral Lorentz matrix maps initial walls \([0,\ldots,8]\) to that list
in order and maps \(H\) to \((7+5s)V\). Its determinant is 1, trace is 6, and
square is the identity. Both the matrix and its inverse preserve the lattice.

For completeness, an actual finite corner of an acute complete Vinberg prefix
persists in the full reflection chamber. Write
\(H=aV+\sum b_i n_i\), where \(a>0\) and
\(b=G_S^{-1}((H,n_i))\le0\), using the entrywise nonnegative inverse of the
spherical vertex Gram. Every later wall \(r\) satisfies
\((H,r)<0\) and \((n_i,r)\le0\); hence
\((H,r)\ge a(V,r)\), forcing \((V,r)<0\). The checker stores and verifies this
decomposition. The integral isometry maps the complete outward local sectors
of these full corners, so it stabilizes the full reflection chamber.

The involution fixes walls 0 through 6 and 12, swaps 7 with 10 and 8 with 11,
and transports wall 9 to a new wall 13. The 14-wall set is invariant under it.
An invariant strict timelike center and all 14 facet points are checked exactly.
The old escape remains strict against wall 13. The involution has finite order;
it supplies **no nonreflectivity conclusion**. No complete 14-wall vertex census
is claimed.

## Artifacts and next bounded decision

- [Complete height-20 producer](../scripts/fresh_frontier_construction_20261015.py)
  and [certificate](fresh_frontier_construction_20261015.json).
- [Complete cone checker](../scripts/fresh_frontier_construction_20261015_geometry.py)
  and [certificate](fresh_frontier_construction_20261015_geometry.json).
- [Finite symmetry and transport](../scripts/fresh_frontier_construction_20261015_symmetry.py)
  and [certificate](fresh_frontier_construction_20261015_symmetry.json).

Generation and saved replay pass for all three packages; source hashes, syntax,
and local links also pass. Independent review is pending. A useful next experiment is one further complete bounded
height interval followed by exact maps between the new norm-colored corners,
using the known involution to avoid duplicated symmetry work. Merely adding
translates that preserve the displayed escape is not finite-volume progress.
