# Power-coupled K3 right pentagons have only five possible fields

**Status (2026-08-11): proved scoped \(n=2\) obstruction.**  In the first
genuinely nonorthogonal compact chamber ansatz, exact order-integrality and
rank-three pentagon closure already bound a field-generating cyclic Cartan
coordinate.  More precisely, an exact rank-six full-trace K3 return whose
prescribed eigen-nef chamber is a compact right-angled pentagon cannot run
through infinitely many real quadratic fields if two chosen nonadjacent
pair coordinates satisfy the additional power coupling

\[
                              C=A^q,
                 \qquad q\in\mathbf Z,\quad q\ge2,
\tag{0.1}
\]

and \(A\) is an algebraic-integer generator of the returned field.  The
field must lie in

\[
 \left\{
 \mathbf Q(\sqrt2),\ \mathbf Q(\sqrt3),\ \mathbf Q(\sqrt5),
 \ \mathbf Q(\sqrt{13}),\ \mathbf Q(\sqrt{17})
 \right\}.
\tag{0.2}
\]

The field-generating pair cycle \(\xi=4A\) satisfies the uniform height
bound

\[
                         4<\xi_0<2(1+\sqrt{69}).
\tag{0.3}
\]

This is not the general \(n=2\) theorem.  For a right pentagon with
independent coordinates \(A,C\), the closure denominator is
\(AC-A-C\); there is no fixed divisor identity analogous to the one below.
The arbitrary nonorthogonal coupling problem remains open.

## 1. The exact carrier package retained

The theorem is conditional on the following complete geometric return.
These hypotheses are deliberately listed separately from the extra ansatz
in Section 2.

Let \(k\) be a real quadratic field, \(\mathcal O\subset k\) a returned
order, and \(M\) a finite projective rank-three \(\mathcal O\)-module with
an admissible symmetric \(k\)-form \(h\).  At the distinguished embedding
\(\sigma_0\), the form has signature \((2,1)\); at the other embedding
\(\sigma_1\), it is positive definite.  Assume all of the following.

1. There is a minimal complex projective K3 surface \(X\) with
   \(\operatorname {NS}(X)\) of rank six equal to the primitive full trace
   lattice returned from \((M,h)\).  The \(\mathcal O\)-action is
   self-adjoint and is realized by divisor-product correspondences.
2. The prescribed projectivized eigen-nef chamber is exactly a compact
   right-angled pentagon \(P\subset\mathbf H^2\), not merely a divisor-wall
   subchamber.
3. Its five facet reflections preserve the returned order-module.  The
   restricted K3 mirror packet is complete and saturated, and its group is
   reflection-maximal among the integral \(k\)-linear reflections of the
   declared order host.

Thus finite volume, exact same-chamber return, order integrality, complete
reflection saturation, and relative reflection maximality are all retained.
Only order integrality is needed for the numerical obstruction below; no
step derives the conclusion by silently weakening the other hypotheses.
The projective nef dimension is \(\operatorname {rank}_k(M)-1=2\).

## 2. The additional power-coupled ansatz

Choose cyclically ordered outward normals
\(e_0,\ldots,e_4\in M\otimes_{\mathcal O}k\).  Adjacent normals are
orthogonal and every nonadjacent pair is ultraparallel.  Put

\[
\begin{array}{lll}
 A={h(e_0,e_2)^2\over h(e_0,e_0)h(e_2,e_2)},
 & C={h(e_1,e_3)^2\over h(e_1,e_1)h(e_3,e_3)},
 & B={h(e_0,e_3)^2\over h(e_0,e_0)h(e_3,e_3)},\\[3mm]
 D={h(e_1,e_4)^2\over h(e_1,e_1)h(e_4,e_4)},
 & E={h(e_2,e_4)^2\over h(e_2,e_2)h(e_4,e_4)}.
\end{array}
\tag{2.1}
\]

Thus all five quantities belong to \(k\).  Formally rescaling the normals
to unit length gives the same ratios, without requiring the square roots of
their norms to belong to \(k\).  At the distinguished place the five ratios
are strictly greater than one.

In addition to the carrier package of Section 1, assume

\[
        A\in\mathcal O,\qquad k=\mathbf Q(A),\qquad C=A^q
        \quad\hbox{for an integer }q\ge2.
\tag{2.2}
\]

The algebraic integrality of \(A\) is an actual extra hypothesis.  Integral
facet reflections give \(4A\in\mathcal O\), but divisibility by four does
not by itself give \(A\in\mathcal O\).  Equation (2.2) is therefore not
being inferred from the exact carrier package.

This ansatz is nonorthogonal: the five nonadjacent wall pairs of a compact
right pentagon are ultraparallel.  It is a coupling condition inside one
rank-three Lorentzian form, not an orthogonal product stabilization.

## 3. Exact right-pentagon closure

Every four of the five normals lie in a three-dimensional \(k\)-space.
The corresponding principal Gram determinants vanish.  Applied to the
normalization (2.1), three of these identities are

\[
\begin{aligned}
 B&=(A-1)(C-1),\\
 E&=(A-1)(D-1),\\
 D&=(C-1)(E-1).
\end{aligned}
\tag{3.1}
\]

Solving the last two and writing

\[
                              \Delta=AC-A-C
\tag{3.2}
\]

gives

\[
             D={A(C-1)\over\Delta},\qquad
             E={C(A-1)\over\Delta}.
\tag{3.3}
\]

Because \(D,E>1\), the denominator is positive.  On substituting
\(C=A^q\), define

\[
                  u=A^q-A^{q-1}-1
                    =A^{q-1}(A-1)-1.
\tag{3.4}
\]

Then \(\Delta=Au\), and (3.3) becomes

\[
          E=1+{1\over u},\qquad
          D=1+{A^{q-1}\over u}.
\tag{3.5}
\]

These are necessary closure identities for any carrier in the declared
ansatz; no existence of such a K3 pentagon is asserted.

## 4. Reflection integrality produces a divisor of four

For two facets, the cyclic Cartan product is

\[
 a_{ij}a_{ji}
 ={4h(e_i,e_j)^2\over h(e_i,e_i)h(e_j,e_j)}.
\tag{4.1}
\]

It belongs to \(\mathcal O\).  Indeed, if \(s_i,s_j\) are the two
reflections, then \(\operatorname {tr}_k(s_is_j)\) is the module trace of
an endomorphism of the finite projective \(\mathcal O\)-module, and it
differs from (4.1) by an integer.  Applied to the \(E\)-pair, this yields
\(4E\in\mathcal O\).  Hence

\[
             v:=4E-4={4\over u}\in\mathcal O.
\tag{4.2}
\]

By (2.2) and (3.4), also \(u\in\mathcal O\), and therefore

\[
                                  uv=4.
\tag{4.3}
\]

This fixed divisor identity is the decisive feature of the power-coupled
locus.

## 5. The signed-place bound

Write \(A_0=\sigma_0(A)\) and \(A_1=\sigma_1(A)\).  Ultraparallelity gives
\(A_0>1\).  At the definite place, strict Cauchy--Schwarz gives

\[
                                  0<A_1<1.
\tag{5.1}
\]

The inequalities are strict because the two wall normals are independent
and their pairing is nonzero over \(k\).

At the distinguished place, \(E_0>1\), so (3.5) gives \(u_0>0\) and
\(v_0>0\).  At the definite place,

\[
 u_1=A_1^q-A_1^{q-1}-1
     =-\bigl(1+A_1^{q-1}(1-A_1)\bigr)\in(-2,-1),
\tag{5.2}
\]

and consequently

\[
                                  v_1={4\over u_1}\in(-4,-2).
\tag{5.3}
\]

Since \(v\) is an algebraic integer,
\(N_{k/\mathbf Q}(v)=-n\) for a positive integer \(n\).  Taking norms in
(4.3) gives

\[
             N(u)N(v)=16,\qquad n\mid16.
\tag{5.4}
\]

Moreover \(n=v_0|v_1|<4v_0\), whence

\[
                 u_0={4\over v_0}<{16\over n}\le16.
\tag{5.5}
\]

Equations (3.4) and (5.5) imply

\[
 A_0^{q-1}(A_0-1)<17.
\tag{5.6}
\]

As \(q\ge2\) and \(A_0>1\), this gives

\[
 A_0(A_0-1)<17,
 \qquad
 A_0<{1+\sqrt{69}\over2}.
\tag{5.7}
\]

This is an automatic bound on the field-generating cyclic product
\(\xi_0=4A_0\).  Reflection maximality is not used in deriving it.

## 6. Exact quadratic enumeration

Because \(A\) is a quadratic algebraic integer, write

\[
 T=A_0+A_1\in\mathbf Z,
 \qquad N=A_0A_1\in\mathbf Z.
\tag{6.1}
\]

Equations (5.1) and (5.7) give

\[
                         2\le T\le5,
                 \qquad 1\le N\le4.
\tag{6.2}
\]

Enumerating these sixteen integer pairs, requiring positive nonsquare
discriminant \(T^2-4N\), and imposing

\[
 {T+\sqrt{T^2-4N}\over2}<{1+\sqrt{69}\over2},
 \qquad
 0<{T-\sqrt{T^2-4N}\over2}<1,
\tag{6.3}
\]

leaves exactly

\[
\begin{array}{c|c|c|c}
T&N&T^2-4N&k\\ \hline
3&1&5&\mathbf Q(\sqrt5)\\
4&1&12&\mathbf Q(\sqrt3)\\
4&2&8&\mathbf Q(\sqrt2)\\
5&2&17&\mathbf Q(\sqrt{17})\\
5&3&13&\mathbf Q(\sqrt{13}).
\end{array}
\tag{6.4}
\]

This proves (0.2) and (0.3).

### Theorem 6.1 (power-coupled right-pentagon obstruction)

Every return satisfying the exact rank-six K3/maximal carrier package of
Section 1 **and** the additional power-coupled algebraic-integrality ansatz
of Section 2 has ground field in the five-element set (0.2).  In particular,
there is no infinite-field family inside this ansatz.

## 7. Exact scope and surviving \(n=2\) problem

The theorem does not:

- construct a K3 surface or a right-pentagon carrier satisfying the ansatz;
- prove a bound for arbitrary compact Coxeter polygons;
- treat a right pentagon whose two initial coordinates \(A,C\) vary
  independently, without the relation \(C=A^q\);
- derive a uniform mixed correspondence degree
  \(H\mathbin\cdot\tau_\theta(H)\) from maximality; or
- prove that reflection maximality alone supplies any height bound.

For independent \(A,C\), equations (3.2)--(3.3) only give

\[
 E={C(A-1)\over AC-A-C}.
\tag{7.1}
\]

Integrality of \(4E\) no longer pairs an algebraic integer with its inverse
to produce a fixed rational integer.  The exact surviving problem is to
control this moving denominator, or the equivalent length-at-most-five
field-generating cycle heights, using the same-chamber K3 geometry and
reflection maximality.  No reflective-lattice classification is used here.
