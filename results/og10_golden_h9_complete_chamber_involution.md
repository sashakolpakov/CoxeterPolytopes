# A complete-chamber involution of the golden fixed \(\mathbb H^9\) slice

> **Subsequent decision.**  This involution becomes decisive when paired with
> the extended reversal `A` of the complete Bugaenko facet:
> `tr(A*g)=8+8*phi>10`.  Their product has infinite order, so the complete H9
> chamber and its H10 parent have infinite volume, and the enhanced birational
> centralizer is infinite.  See
> [`og10_golden_h9_infinite_stabilizer.md`](og10_golden_h9_infinite_stabilizer.md).

## Result

Let

\[
L_9=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle,
\qquad (x,z,u)\in L_9,
\]

and let \(P_9\) be the complete norm-two chamber obtained as the fixed slice
of the golden OG10 chamber.  On the \((x,u)\)-plane consider the Pell
translation

\[
T=\begin{pmatrix}\phi&1\\ \phi&\phi\end{pmatrix},
\qquad
T^{-1}=\begin{pmatrix}\phi&-1\\-\phi&\phi\end{pmatrix}.
\tag{1}
\]

Exact chamber folding sends both maps to

\[
g=\begin{pmatrix}\phi&-1\\ \phi&-\phi\end{pmatrix}\oplus I_8.
\tag{2}
\]

More precisely, \(T\) crosses only the wall

\[
r_+=(1,0^8,\phi),
\]

whereas \(T^{-1}\) crosses only

\[
r_-=(0,0^8,-1),
\]

and

\[
s_{r_+}T=s_{r_-}T^{-1}=g.
\tag{3}
\]

The endpoint is strict for every one of the 61 archived fixed walls.  The
global Lorentzian estimate below proves that no later simple wall is
violated.  Consequently

\[
                         gP_9=P_9
\]

for the **complete** chamber, not merely for the known prefix.  Moreover

\[
g^2=1,
\qquad
g=s_\eta,
\qquad
\eta=(1,0^8,\phi^2),
\qquad
\eta^2=4+4\phi.
\tag{4}
\]

Thus the infinite Pell motion does not yield an infinite-order chamber
symmetry: its two directions fold to the same non-root reflection.

## Complete-membership bound

Use the invariant interior point inherited from the ambient chamber and put
\(q'=gq\).  In fixed coordinates it is

\[
q'=\left(
\frac{197}{400}+2\phi,
\frac{23}{50},\frac{91}{100},\frac{27}{20},\frac{11}{10},
\frac{21}{25},\frac{57}{100},\frac{29}{100},\frac{17}{25},
\frac{199}{400}+\frac{399}{200}\phi
\right).
\tag{5}
\]

Write \(q'=ae_0+v\), with \(v\perp e_0\).  A norm-two root
\(r=xe_0+w\) has

\[
w^2=2+2\phi x^2.
\]

If \((q',r)\geq0\), positive-definite Cauchy--Schwarz on \(e_0^\perp\)
gives

\[
x^2\leq
\frac{v^2}{\phi(2\phi a^2-v^2)}
=
\frac{202727309617}{200773372141}
-\frac{1021288932}{200773372141}\phi
\approx1.00150148034.
\tag{6}
\]

The first coordinate unprocessed at the time of this proof was
\(3+6\phi\), whose square is

\[
(3+6\phi)^2=45+72\phi\approx161.498.
\tag{7}
\]

The ambient enumeration is exhaustive below that coordinate, and every
fixed simple root is an ambient simple root with \(y=0\).  Since the largest
pairing of (5) with the 61 fixed walls is already \(-1/100\), (6)--(7)
exclude every unseen violator.  This is the complete-chamber step.

## Why this does not lift to the golden \(\mathbb H^{10}\) chamber

The ambient field lattice has coordinates \((x,y,z,u)\), fixed hyperplane
\(y=0\), and orthogonal normal

\[
\nu=(\phi-1,2,0^9).
\]

Any \(K=\mathbb Q(\sqrt5)\)-linear isometry preserving \(y=0\) and
restricting to \(g\) must send \(\nu\) to \(\pm\nu\).  These are the only
two possible lifts.  In either one,

\[
u'=\phi x-\frac12y-\phi u.
\tag{8}
\]

Therefore neither lift preserves the \(\mathbb Z[\phi]\)-lattice.  In
particular, \(g\) is an intrinsic symmetry of \(P_9\), but it is **not** a
symmetry of the golden \(\mathbb H^{10}\) lattice or its complete Coxeter
chamber.

## Quarantined full-lattice observation

The obstruction above is specifically an obstruction to a golden-field
\(\mathbb H^{10}\) lift.  It is not an obstruction in the abstract rank-24
OG10 Beauville lattice.  Exact discriminant gluing gives a positive,
infinite-order integral isometry of that lattice whose restriction to the
rank-20 \(h\)-fixed lattice is \(g\).  In a basis of the transverse lattice

\[
\langle2\rangle^2\oplus\langle-2\rangle\oplus\langle-6\rangle
\]

its transverse block is

\[
B=\begin{pmatrix}
1&1&-1&0\\
-1&2&-1&3\\
-1&1&0&3\\
0&-1&1&-1
\end{pmatrix},
\qquad
\chi_B(t)=(t-1)^3(t+1).
\tag{9}
\]

This extension does **not** preserve the original period plane,
Néron--Severi lattice, golden \(\mathbb H^{10}\) carrier, or ambient Coxeter
chamber.  It is consequently not an automorphism or birational map of the
archived OG10 period model and does not itself advance the covolume question.
It is recorded only to prevent the stronger, false claim that no abstract
Beauville-lattice extension exists.  The later covolume conclusion comes from
the intrinsic H9 stabilizer and the established primal--dual equivalence, not
from this abstract extension.

## Scope

This theorem supplies one complete-chamber symmetry of \(P_9\).  By itself it
does not establish finite volume or compactness, and it creates no new
finite-volume polytope.  The subsequent infinite-stabilizer theorem combines
it with the Bugaenko-facet reversal and proves \(P_9\) and the parent H10
chamber infinite-volume.  The explicit non-lift statement remains essential:
the enhanced birational conclusion is indirect.

The machine-readable certificate is
[`og10_golden_h9_complete_chamber_involution.json`](../data/candidates/og10_golden_h9_complete_chamber_involution.json),
recomputed by
[`certify_og10_golden_h9_complete_chamber_involution.py`](../scripts/certify_og10_golden_h9_complete_chamber_involution.py).
