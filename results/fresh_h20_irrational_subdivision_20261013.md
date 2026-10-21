# All-real subdivision obstruction for the 59 twin-leaf section planes

Every full-dimensional section in the **59 frozen twin-leaf planes** fails to admit a finite **face-to-face subdivision into full-dimensional hyperbolic Coxeter polytopes**. This includes irrational section parameters and arbitrary Coxeter labels. It does not exclude other section planes, changed source polytopes, or a construction that is not a face-to-face subdivision of the whole section. No new Coxeter polytope is claimed.

The [exact certificate](fresh_h20_irrational_subdivision_20261013.json) checks 118 orientations: 114 use an actual spherical codimension-three link, and four use two actual ridges and an algebraic-integer contradiction. They give 57 and two unoriented planes, respectively. The [producer and verifier](../scripts/fresh_h20_irrational_subdivision_20261013.py) uses only exact rational arithmetic and the elementary matrix helpers of the [frozen section certificate](fresh_h20_cusp_sections_20261013.md). It imports no numerical scout. Earlier finite-angle numerical experiments are unnecessary to this all-parameter result.

A [separate Fraction reconstruction](../scripts/fresh_h20_irrational_subdivision_20261013_review.py)
also [passes](fresh_h20_irrational_subdivision_20261013_review.json). It
imports no producer or repository arithmetic helper, binds all 118 rows
to the earlier section packet, and independently checks 187 distinct
face points, all 122 interpolations, every projected Gram polynomial,
and both exceptional ridge formulas. The spherical tiling argument was
reviewed by two other agents; the cyclotomic argument was separately
reviewed by the root agent. These are scoped arithmetic and mathematical
checks, not formal proof-assistant verification.

## Parameter range and actual face certificates

Use the frozen norm-two orthogonal twin walls \(a,b\), their future null cusp \(\delta\), and
\[
 e_0=(a-b)/2,\qquad e(t)=e_0+t\delta.
\]
The checker verifies \(e_0^2=1\), \((e_0,\delta)=\delta^2=0\), and both orientations of each plane. Reversing the twins changes \(t\) to \(-t\). The normal \(e(1/2)\) is a strictly positive combination of the component walls other than \(b\). The reverse orientation gives the analogous identity for \(-e_0+\delta/2\), so \(\delta\) is also a positive wall combination. Consequently \(|t|\geq1/2\) cannot give a section meeting the strict interior of the source polytope. A pure null normal \(\delta\) does not cut hyperbolic space. It suffices to treat \(0\leq t<1/2\).

Every claimed face is certified against all **210** source inequalities, rather than inferred from a small positive-definite Gram block alone. For the specified left and right ambient faces, projection of the strict timelike center gives future timelike points \(z_L,z_R\), with exactly the saved zero sets and all other inequalities strictly negative. Put
\[
 A=-(\delta,z_L)>0,\quad B=-(\delta,z_R)>0,
 \quad x(t)=(1/2-t)Bz_L+tAz_R.
\]
The identities \((e_0,z_L)=0\) and \((e(1/2),z_R)=0\) imply \((e(t),x(t))=0\). For \(0<t<1/2\), this is a positive combination of future timelike points. Its zero walls are exactly the intersection of the two ambient face zero sets; every other inequality is strict. At \(t=0\), the extra zero wall \(b\) has the same projected normal as \(a\), since both project to \((a+b)/2\). Thus the stated ridge or codimension-three face is still actual, with no extra distinct incident facet.

## The 57 spherical-link planes

In these planes the twins have a norm-two common neighbor \(c\), with \((a,c)=(b,c)=-1\). The selected auxiliary wall \(j\) satisfies \(j^2=2\), \((c,j)=0\), \((\delta,j)=-1\), and \((e_0,j)=\nu\in\{-1/2,1/2\}\). The ambient left face is \(\{a,b,c,j\}\); the right face contains the component except \(b\), together with \(j\). Their intersection gives exactly the actual section face \(\{a,c,j\}\).

In the order \(a,c,j\), the projected Gram matrix is
\[
 \begin{pmatrix}
 1&-1&t-1/2\\
 -1&2&0\\
 t-1/2&0&7/4+2\nu t-t^2
 \end{pmatrix}.
\]
It is positive definite on the whole parameter interval. The corresponding spherical normal link is a triangle with angles
\[
 \pi/4,\quad\pi/2,\quad\theta,
 \qquad
 \cos^2\theta=\frac{(1/2-t)^2}{7/4+2\nu t-t^2},
 \qquad \pi/3<\theta<\pi/2.
\]

Here is a local obstruction that applies to **any** such triangle. All three side lengths are strictly less than \(\pi/2\), by the spherical cosine law. The whole triangle has diameter less than \(\pi/2\): its vertices have positive mutual dot products, and every point is a normalized positive combination of them.

Suppose a finite face-to-face Coxeter subdivision of the section existed. Choose a generic point of the original codimension-three face, avoiding all subdivision walls that do not contain that face. Its induced spherical tiling is face-to-face. A tile of local rank less than three would contain a line in the transverse cone and could not fit in the pointed triangular cone. Thus the tiles are rank-three spherical Coxeter chambers. Reducible chambers have a side of length \(\pi/2\), so cannot fit in this triangle. The remaining possibilities are precisely
\[
 A_3:(2,3,3),\quad B_3:(2,3,4),\quad H_3:(2,3,5).
\]
Indeed, a spherical Coxeter triangle has positive area, so its reciprocal label sum exceeds one; after ordering the labels, the irreducible cases are exactly those three.

At the \(\theta\)-corner one tile angle is impossible, and three tile angles are too large because each is at least \(\pi/5\). The only possible sums of two angles strictly between \(\pi/3\) and \(\pi/2\) are \(2\pi/5\) and \(9\pi/20\). In either case an \(H_3\) tile occurs. But the exact squared cosines of the side lengths are

| Tile | Squared side cosines |
|---|---|
| \(A_3\) | \(1/9,1/3,1/3\) |
| \(B_3\) | \(1/3,1/2,2/3\) |
| \(H_3\) | \(1/3+2\sqrt5/15,1/2+\sqrt5/10,1/2+\sqrt5/6\) |

An \(H_3\) tile cannot share a full side with either other type. Face-to-face tiling and connectedness of the tile adjacency graph therefore force every tile to be \(H_3\). Such tiles cannot give the \(9\pi/20\) corner. At \(\theta=2\pi/5\), the whole triangle has area \(3\pi/20\), whereas an \(H_3\) tile has area \(\pi/30\); the ratio is \(9/2\), not an integer. Both cases are impossible. This argument also covers \(t=0\).

## The two remaining planes

The remaining orientations are source rows **4, 5, 18, 19**. Their common neighbor is long, so the preceding finite link is not used. The certificate instead constructs two actual ridges by the same strict-face interpolation. Their squared cosine invariants are
\[
 A(t)=4\cos^2\theta_A=
 \frac{4(1/2-t)^2}{7/4+t-t^2},\qquad
 D(t)=4\cos^2\theta_D=
 \frac{4(1-4t)^2}{4-(1-4t)^2}.
\]
The second ridge may be obtuse; the argument only uses its squared cosine. A finite face-to-face Coxeter subdivision would make each original actual ridge angle a finite sum of angles \(\pi/m\), hence a rational multiple of \(\pi\). Therefore both \(A,D\) would be algebraic integers, and every conjugate would lie in \([0,4]\).

The equation
\[
 (1-4t)^2=\frac{4D}{D+4}
\]
first shows that \(t\) is algebraic, and then shows that every embedding of \(\mathbb Q(t)\) takes \(t\) to a real number in
\[
 [(1-\sqrt2)/4,(1+\sqrt2)/4]\subset(-1/8,5/8).
\]
On that rational outer interval, the denominator of \(A\) is positive and
\[
 (7/4+t-t^2)-4(1/2-t)^2=3/4+5t-5t^2>0.
\]
This concave quadratic has endpoint values \(3/64\) and \(123/64\), both positive. Thus every conjugate of \(A\) lies in \([0,1)\). Physically \(A>0\), since \(t<1/2\), so \(A\ne0\). Every embedding of \(\mathbb Q(A)\) extends to \(\mathbb Q(t)\); its nonzero algebraic-integer norm would consequently have absolute value strictly between zero and one. This is impossible. No angle-label bound or parameter sampling enters this argument.

## Replay and scope

```sh
python scripts/fresh_h20_irrational_subdivision_20261013.py --verify
python scripts/fresh_h20_irrational_subdivision_20261013_review.py
```

The exact replay checks the source hashes, all 118 oriented face constructions, the 59 reversal pairs, positive endpoint combinations, all off-face inequalities, the spherical side arithmetic and the exceptional interval inequalities. The proof above supplies the finite-tiling and algebraic-number implications; a replay is not a formal proof-assistant verification of those arguments.

This closes the proposed irrational-parameter subdivision route in these 59 specific planes. It is stronger than the earlier statement that the sections themselves are not Coxeter polytopes, but remains a local obstruction within this fixed section family.
