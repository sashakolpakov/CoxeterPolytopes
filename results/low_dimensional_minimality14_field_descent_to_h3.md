# Minimality-14 field descent in dimensions four and five

**Status (2026-08-12).**  This note independently reconstructs the exact
codimension-one/codimension-two descent behind the low-dimensional part of
the narrow-face method.  It uses no arithmetic-reflection finiteness theorem.
Combined with the finite transition integral locus and the exact
two-dimensional orbifold-curve theorem, it shows that the only remaining
ground-field obstruction in dimensions \(3,4,5\) is compact dimension three.

## 1. The theorem

Let \(P\subset\mathbb H^n\), \(n\geq4\), be a finite-volume arithmetic
Coxeter chamber with minimal cycle field \(k\).  Let \(E\) be a
minimal-aperture facet, so every normalized pairing between closed facets
meeting \(E\) is less than \(14\).

If \(P\) is noncompact, its cusp makes the defining \(k\)-form isotropic;
definiteness at every nonidentity real place then forces \(k=\mathbf Q\).
Thus the only case needing descent is compact, and all faces used below are
compact.

A compact hyperbolic Coxeter chamber is simple: the link of each finite
vertex is a spherical Coxeter simplex.  Consequently two meeting facets
share a ridge, every face is simple, a wall cutting a face cuts it in a
facet, and an edge leaving that facet has a unique wall at each endpoint
which does not contain the whole edge.  These elementary link facts justify
the face dimensions and the endpoint-exclusive packet used below.

> **Theorem 1.1 (codimension descent or transition).**  At least one of the
> following holds:
>
> 1. \(k=\mathbf Q\), or \(k\) belongs to the finite Lanner/triangle field
>    packets;
> 2. \(k\) is the full field of an exact finite-volume arithmetic Coxeter
>    chamber of minimality \(14\) in dimension \(n-1\);
> 3. \(k\) is the full field of such a chamber in dimension \(n-2\); or
> 4. \(k\) belongs to the finite integral locus of one of the first four
>    four-node transition arrangement charts.

Consequently

\[
 \begin{aligned}
 \mathscr K_4&\subset
   \mathscr K_3\cup\mathscr K_2\cup\mathscr K_{\rm tr}\cup\mathscr K_{\rm small},\\
 \mathscr K_5&\subset
   \mathscr K_3\cup\mathscr K_2\cup\mathscr K_{\rm tr}\cup\mathscr K_{\rm small},
 \end{aligned}                                             \tag{1.1}
\]

where \(\mathscr K_2\), \(\mathscr K_{\rm tr}\), and
\(\mathscr K_{\rm small}\) are finite by the independently proved
orbifold-curve, contracting-cycle, and Lanner results.  Therefore a finite
field theorem in compact \(\mathbb H^3\) would close the field problem in
dimensions four and five; no other low-dimensional field input is missing.

This does not say that every pentagonal two-face is an arithmetic Coxeter
pentagon.  The all-right pentagon is precisely why the codimension descent is
needed instead of a naive polygon-face argument.

## 2. First descent: the narrow facet itself

The minus-one eigenspace \(\ell_e\subset V\) of the \(k\)-rational
reflection in \(E\) is a \(k\)-line. Choose a nonzero \(k\)-vector \(e\)
on \(\ell_e\); use a square-\((-2)\) representative only after extending
to the distinguished real place when normalized Gram entries are discussed.
Let
\(\mathcal S(E)\) be the closed star: the walls whose closed facets meet
\(E\).

Suppose every \(\delta\in\mathcal S(E)\setminus\{e\}\) is orthogonal to
\(e\).  Put \(U=e^\perp\).  The intrinsic normal of the corresponding wall
on \(E\) is ordinarily its orthogonal projection to \(U\); here it is
literally \(\delta\).  Hence all norms, pairings, angles, parallelisms, and
dotted distances are unchanged.

The restricted reflections preserve \(U\) and \(L\cap U\) for an ambient
arithmetic lattice \(L\).  Their finite-volume chamber is exactly \(E\).
The wall packet spans the Lorentzian \(U\), so its Gram has a positive
eigenvalue and the full-field inheritance lemma gives cycle field \(k\).
At every nonselected place the restricted form is definite.  Thus \(E\) is
an exact arithmetic Coxeter chamber over \(k\), of dimension \(n-1\) and
minimality \(14\).  This is alternative 2.

## 3. Second descent: a codimension-two face

Otherwise choose a neighboring wall with \(e.f>0\), and choose a nonzero
\(k\)-vector \(f\) on its \(k\)-rational minus-one eigenspace \(\ell_f\).
Since the two facets meet, their normal lines span a negative-definite plane.
Put

\[
                    F=P\cap H_e\cap H_f,
 \qquad U=(\ell_e+\ell_f)^\perp
          =(\operatorname {span}_k\{e,f\})^\perp.         \tag{3.1}
\]

If every other wall \(\delta\) meeting \(F\) is orthogonal to both \(e\)
and \(f\), its projection to \(U\) is again literal.  The argument of
Section 2 makes \(F\) an exact arithmetic Coxeter chamber over the full
field \(k\), now of dimension \(n-2\) and minimality \(14\).  This is
alternative 3.

## 4. Failure of descent forces four transition walls

Suppose neither descent occurs.  Then there is a wall \(g\) meeting \(F\)
with \(g.e>0\) or \(g.f>0\).  The three walls \(e,f,g\) meet on a compact
face, so their principal Coxeter diagram is negative definite; because
\(e.f>0\) and \(g\) meets at least one of them, it is connected.

The face \(F\) has dimension at least two.  Choose a one-dimensional edge
\(r\subset F\) with one endpoint on \(H_g\) and the other endpoint off
\(H_g\).  Let \(h\) be the wall containing the latter endpoint but not the
whole edge.  The walls \(g,h\) are the endpoint-exclusive pair of the finite
edge chamber of \(r\).

If \(g.h=2\), their independent parabolic packet has the full field, forcing
\(k=\mathbf Q\).  If \(g.h<2\), choose an inclusion-minimal principal packet
with a positive direction inside the complete wall packet of the finite edge
\(r\).  It has rank greater than two.  If all its proper principal packets
are elliptic, it is a triangle/Lanner packet; if a proper principal packet is
parabolic, it is the quasi-Lanner noncompact-simplex case and its rational
cusp forces \(k=\mathbf Q\).  This is alternative 1.  Outside those small
sets one therefore has

\[
                             g.h>2.                         \tag{4.1}
\]

Now \(\{g,h;e,f\}\) is a connected four-node finite-edge packet: \(g,h\)
are its dotted endpoint pair and \(e,f,g\) are a connected spherical triple.
Every wall meets \(E\), so its minimality is \(14\).  The exact seven-support
enumeration places it in transition family 1, 2, 3, or 4.  Family 5 is
impossible because, after deleting \(h\), the component containing \(g\)
contains the full three-node packet \(e,f,g\), whereas both endpoint-deletion
components in family 5 have only two vertices.

The dotted pair gives the four-node Gram a positive eigenvalue.  Full-field
inheritance therefore identifies its cyclic field with \(k\), and the
independently proved contracting tail/fork theorem makes the possible fields
finite.  This is alternative 4.

## 5. Exact remaining obstruction

Iterate Theorem 1.1 once for \(n=5\).  Every dimension-five field descends to
dimension four or three, or enters the finite transition/small packets.  A
dimension-four descendant again lands in dimension three or two, or in those
finite packets.  Dimension two is closed by the exact real-DM orbifold-curve
theorem.  Thus only dimension three remains.

The obstruction in \(\mathbb H^3\) is not cosmetic.  A finite edge packet
may split as a two-node Lorentz block (the dotted endpoint pair) orthogonal
to spherical blocks.  Its primitive integral coordinate has

\[
                 \alpha_0>4,\qquad0<\alpha^\sigma<4,        \tag{5.1}
\]

and \([0,4]\) has integer Chebyshev capacity one.  The compact absolutely
maximal tetrahedron \(H_6\) also shows that no automatic rank-three
Lorentzian subdiagram, Coxeter facet, or genus-zero mirror quotient is
available.  Closing compact dimension three requires a new determinant
coupling, a richer bounded AG carrier, or an imported three-dimensional
arithmetic-reflection theorem.  None is supplied here.

## 6. Relation to the all-dimensional goal

This descent improves the location of the arithmetic gap, but it does not
prove chamber finiteness.  Even after finite fields are known, one still
needs a dimension ceiling, a global facet/dotted-height bound, and native
lattice/glue finiteness.  The local transition charts and descended faces do
not recover the transverse Gram extension.
