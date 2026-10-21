# Marked minimal-base single-chain simplices have a doubled-wall escape

**Status (2026-08-10): exact scoped no-go for the direct inverse-diagonal
depth argument.**  No reflective-lattice classification, chamber census, or
Nikulin--Vinberg bound is used.

The plane result has two parts.

1.  In the standard plane root frame the visible simplex has finite volume
    exactly through chain length nine.
2.  That threshold does **not** bound the chain.  At the first hyperideal
    ray, physical adjunction says that only a smooth rational \((-4)\)-wall
    can cut.  On the direct full Neron--Severi host its integrality forces it
    to be doubled and supported at the escape stage.  The same one wall then
    cuts every later member of the original hyperideal-ray family.

The second statement is an obstruction to the proposed proof method, not an
unbounded geometric counterfamily.  After the doubled wall is inserted, new
mixed extreme rays appear; the first is already hyperideal at chain length
twelve.  Bounding those mixed rays remains the genuine global gate.

The surviving direct minimal ruled base has the same mechanism one stage
earlier.  On \(\mathbb F_2\), the universal negative-section/affine-fibre
simplex has visible threshold \(m=8\), first escape
\(-K_{S_9}\), one-wall repair only for \(m=9,10\), and first mixed escape at
\(m=11\).  Section 8 proves that companion calculation independently.

## 1. Exact scope

Let

\[
 S_m\longrightarrow \mathbf P^2
\]

be a succession of \(m\) point blowups, and use orthogonal total-transform
coordinates

\[
 H^2=1,\qquad e_i^2=-1,\qquad H.e_i=e_i.e_k=0\quad(i\ne k).
\]

This note studies the **marked standard root simplex**

\[
 \alpha _0=H-e_1-e_2-e_3,\qquad
 \alpha_i=e_i-e_{i+1}\ (1\le i<m),\qquad
 \alpha_m=e_m.                                      \tag{1.1}
\]

When used geometrically, the intended lane is the exact free-chain frame:
each \(e_i-e_{i+1}\) is the strict exceptional curve, \(e_m\) is the last
exceptional curve, and the line class \(\alpha _0\) follows exactly the
first three marked centres and is not blown up again.  These conditions are
**not automatic for an arbitrary curvilinear cluster**.  Nothing below
claims that full-NS reflection maximality by itself forces \(\alpha _0\) to
be a facet.

Any additional irreducible curve in this free-chain frame is the strict
transform of a plane curve and can be written

\[
 C=dH-\sum_{i=1}^m\mu_i e_i,\qquad \mu_i\ge0.        \tag{1.2}
\]

The nonnegativity in (1.2) is essential in the hidden-wall argument.  An
arbitrary satellite chain or an arbitrary exceptional class lies outside
that exact statement.

## 2. Closed-form inverse Gram matrix

The roots (1.1) are an integral basis.  Their dual basis is

\[
 v_0=H,
 \qquad
 v_i=\min(i,3)H-e_1-\cdots-e_i\quad(1\le i\le m).    \tag{2.1}
\]

Indeed, direct intersection gives

\[
                         v_i.\alpha_j=\delta_{ij}.    \tag{2.2}
\]

Consequently

\[
 (G_m^{-1})_{00}=1,
 \qquad
 (G_m^{-1})_{ii}=\min(i,3)^2-i,                       \tag{2.3}
\]

and hence

\[
 \operatorname {diag}(G_m^{-1})
   =(1,0,2,6,5,4,3,2,1,0,-1,-2,\ldots,9-m).          \tag{2.4}
\]

A hyperbolic simplex has finite volume exactly when each deleted-wall dual
ray is finite or ideal, which here is exactly nonnegativity of every entry
in (2.4).  Thus

\[
 \boxed{\text{the marked standard simplex has finite volume}
        \quad\Longleftrightarrow\quad m\le9.}         \tag{2.5}
\]

Starting at \(j=10\), the escaping rays have the uniform form

\[
 v_j=3H-e_1-\cdots-e_j,\qquad v_j^2=9-j.              \tag{2.6}
\]

For \(K_{S_m}=-3H+e_1+\cdots+e_m\),

\[
 K_{S_m}.v_j=j-9=-v_j^2,
 \qquad p_a(v_j)=1.                                  \tag{2.7}
\]

The first ray is

\[
                         v_{10}=-K_{S_{10}},
 \qquad v_{10}^2=-1.                                 \tag{2.8}
\]

Thus the first failure is intrinsically genus one.  Treating \(v_{10}\)
itself as a smooth rational \((-1)\)-wall would violate adjunction.

## 3. Which physical wall can cut an escape?

Let \(C\) be a smooth rational plane strict transform (1.2).  Since

\[
 v_j=-K_{S_m}+\sum_{i>j}e_i,
\]

adjunction gives the exact identity

\[
 \boxed{
 C.v_j=C^2+2+\sum_{i>j}\mu_i.}                        \tag{3.1}
\]

Here \((-K).C=C^2+2\) because \(p_a(C)=0\).  For the complete smooth
rational reflection alphabet,

\[
\begin{array}{c|c}
C^2&C.v_j\\ \hline
-1&1+\sum_{i>j}\mu_i\\
-2&  \sum_{i>j}\mu_i\\
-4&-2+\sum_{i>j}\mu_i.
\end{array}                                           \tag{3.2}
\]

Therefore a \((-1)\)- or \((-2)\)-curve cannot remove any \(v_j\).
Every finite-volume completion of the marked simplex at \(m\ge10\) must
contain a hidden smooth rational \((-4)\)-wall with

\[
                         \sum_{i>j}\mu_i\le1          \tag{3.3}
\]

for each original escape which it removes.

Equation (3.1) is about an actual curve \(C\), whereas the reflection mirror
depends only on its primitive half-root when one exists.  This distinction is
decisive in the direct full-NS lane.

## 4. Full-NS doubling kills the diagonal depth argument

Let \(L=\operatorname {NS}(S_m)\cong I_{1,m}\).  It is unimodular.  If a
square-minus-four class \(C\) acts by an integral reflection on \(L\), then

\[
 s_C(x)=x+{x.C\over2}C
\]

is integral for every \(x\in L\).  Unimodularity forces the divisibility of
\(C\) to be even, hence

\[
                    C=2r,\qquad r^2=-1.               \tag{4.1}
\]

In the plane coordinates (1.2), all \(\mu_i\) are consequently even.
Combining this with (3.3) sharpens the cutting condition to

\[
 \boxed{C.v_j<0\quad\Longrightarrow\quad
        \sum_{i>j}\mu_i=0,qquad C.v_j=-2.}           \tag{4.2}
\]

Thus a direct full-NS physical wall which cuts \(v_j\) is pulled back from
stage \(j\): every later centre lies off that curve.  But then for every
\(j'\ge j\),

\[
                             C.v_{j'}=-2.              \tag{4.3}
\]

One doubled wall cuts the **entire later family** of old simplex escapes.
Checking \(v_{10},v_{11},\ldots\) one at a time can therefore never yield a
chain-length bound.

On an adapted half-boundary host, (3.1) still holds but the full-NS evenness
need not.  The exact conclusion there is only
\(\sum_{i>j}\mu_i\le1\).  A wall already supported at stage \(j\) still
cuts all later \(v_{j'}\), so the same inverse-diagonal no-go survives;
the stronger parity conclusion (4.2) is direct-full-NS only.

## 5. Exact numerical countermodel and the first mixed debt

The simplest doubled numerical wall is

\[
 B_{10}=2v_{10}
       =6H-2e_1-\cdots-2e_{10}.                       \tag{5.1}
\]

After any number of later blowups away from it,

\[
 B_{10}^2=-4,qquad K.B_{10}=2,qquad p_a(B_{10})=0,
 \qquad B_{10}.v_j=-2\quad(j\ge10).                  \tag{5.2}
\]

Hence (5.1) satisfies the exact numerical smooth-rational adjunction and
full-NS reflection-integrality tests in every rank.  This note does **not**
assert that an irreducible smooth representative exists for an arbitrary
chosen cluster.  The reflecting root is the genus-one numerical class
\(v_{10}\); a reduced physical curve, if it exists, would instead have the
doubled class \(B_{10}=2v_{10}\) and arithmetic genus zero.  Thus (5.1) is
the promised numerical countermodel to the diagonal argument, not a
geometric carrier.

There is also a precise warning against overreading (5.2).  Orient every
physical wall by \(x.C\ge0\), and add the one mirror
\(v_{10}^{\perp}=B_{10}^{\perp}\) to the standard simplex.  Write a vector
in the old dual basis as \(x=\sum t_i v_i\), \(t_i\ge0\).  The new
inequality has coefficients

\[
 p_i=v_i.v_{10}=(3,2,4,6,5,4,3,2,1,0,-1,-1,\ldots). \tag{5.3}
\]

Besides the surviving \(v_i\), the new extreme rays are

\[
                         v_i+p_i v_j
 \quad(p_i>0,\ j\ge10).                               \tag{5.4}
\]

This list is exhaustive: an extreme ray either retains one old coordinate
with \(p_i\ge0\), or lies on the new equality with exactly one positive and
one negative old coordinate.  For every mixed ray in (5.4),

\[
 (v_i+p_i v_j)^2=v_i^2+p_i^2(11-j).                  \tag{5.5}
\]

All of these squares are nonnegative for \(j\le11\).  At \(j=12\), the
choice \(i=1\), for which \(v_1^2=0\) and \(p_1=2\), is already negative.
Thus checking **all** surviving and pair rays proves that the one-wall
augmented Coxeter cone has finite volume exactly for

\[
                              m=10,11.                 \tag{5.6}
\]

At \(m=12\), for example,

\[
                         w=v_1+2v_{12},
 \qquad w^2=4(11-12)=-4.                              \tag{5.7}
\]

So the doubled wall removes every old \(v_j\) while creating a new mixed
hyperideal direction.  This is exactly why the remaining theorem must
control mixed-root degree or a global contraction orbit, rather than repeat
the inverse-diagonal calculation.

## 6. What is and is not proved

Proved:

- the exact inverse diagonal (2.4) and the sharp visible threshold \(m=9\);
- the genus-one nature of every old tail ray;
- the physical-wall identity (3.1);
- the direct full-NS doubling/support theorem (4.2)--(4.3);
- the exact first mixed escape (5.7).

Not proved:

- that the standard line root \(\alpha_0\) occurs for every curvilinear
  cluster;
- effectivity or irreducibility of the numerical class \(B_{10}\);
- an unbounded family of complete nef chambers;
- a uniform bound, or its failure, under the full geometric hypotheses;
- closure of the adapted-host parity lane.

The usable conclusion is deliberately narrow but global: **the original
simplex escape family is exhausted as a possible depth invariant.**  Any
successful single-chain theorem must bound the new mixed rays, force a
bounded-degree physical truncator, or prove a finite equivariant contraction
orbit.

## 7. Exact replay

- [certifier](../scripts/certify_rational_single_chain_standard_simplex_gate.py)
- [JSON archive](../data/enumerations/rational_single_chain_standard_simplex_gate.json)
- [regression tests](../tests/test_rational_single_chain_standard_simplex_gate.py)

The certifier constructs the dual basis directly, checks every Kronecker
pairing, verifies adjunction and doubling through rank 32, derives all
one-wall extreme rays without importing another project module, and freezes
the first mixed hyperideal ray.

## 8. The surviving minimal ruled base \(\mathbb F_2\)

The direct minimal ruled endpoint has an equally rigid universal frame.  Let
\(s\) be the negative section and \(f\) a fibre on \(\mathbb F_2\):

\[
 s^2=-2,\qquad f^2=0,\qquad s.f=1.                    \tag{8.1}
\]

Blow up first a point of a fibre and then the intersection of its strict
transform with the first exceptional curve.  Continue thereafter at free
points of the newest exceptional curve, always off \(s\), and use the walls

\[
 s,\qquad
 \delta_0=f-e_1-e_2,\qquad
 \delta_i=e_i-e_{i+1}\ (1\le i<m),\qquad
 \delta_m=e_m.                                       \tag{8.2}
\]

Geometrically, (8.2) again means the exact marked frame.  The second centre
is the fibre-following free infinitely-near point where the strict fibre
meets \(e_1\); it is proximate to the single exceptional component \(e_1\),
so it is not satellite in the standard proximity terminology.  Its special
position on the strict fibre is what makes both \(f-e_1-e_2\) and
\(e_1-e_2\) physical \((-2)\)-curves.  Every later exceptional difference
shown is the strict exceptional curve, and no centre lies on \(s\).  It is
not a claim that an arbitrary chain on a ruled surface has all these
physical walls.

### 8.1 Duals and visible finite volume

In the wall order in (8.2), the dual rays are

\[
 u_s=f,\qquad u_0=s+2f,
\]

and

\[
 u_i=\min(i,2)(s+2f)-e_1-\cdots-e_i
 \quad(1\le i\le m).                                 \tag{8.3}
\]

Direct intersection gives the full Kronecker pairing.  Since
\((s+2f)^2=2\), the inverse diagonal is

\[
 \operatorname {diag}(G_{\mathbb F_2,m}^{-1})
   =(0,2,1,6,5,4,3,2,1,0,-1,-2,\ldots,8-m).          \tag{8.4}
\]

Therefore

\[
 \boxed{\text{the visible \(\mathbb F_2\) simplex has finite volume}
        \quad\Longleftrightarrow\quad m\le8.}         \tag{8.5}
\]

For \(j\ge9\), the old escaping rays are

\[
 u_j=2s+4f-e_1-\cdots-e_j,
 \qquad u_j^2=8-j.                                   \tag{8.6}
\]

Since

\[
 K_{S_m}=-2s-4f+e_1+\cdots+e_m,
\]

one has

\[
 K.u_j=j-8=-u_j^2,
 \qquad p_a(u_j)=1,
 \qquad u_9=-K_{S_9}.                                \tag{8.7}
\]

Thus the first ruled escape is again a genus-one numerical root, not a
smooth rational \((-1)\)-curve.

### 8.2 The same physical-wall and doubling theorem

Apart from the already displayed negative section and exceptional curves,
write a hidden irreducible strict transform as

\[
 C=as+bf-\sum_i\mu_i e_i,\qquad \mu_i\ge0.           \tag{8.8}
\]

The equality

\[
 u_j=-K_{S_m}+\sum_{i>j}e_i
\]

and smooth-rational adjunction give exactly

\[
                         C.u_j=C^2+2+\sum_{i>j}\mu_i. \tag{8.9}
\]

Hence only a \((-4)\)-curve can cut an old \(u_j\).  The lattice
\(\operatorname {NS}(\mathbb F_2)\perp\langle-1\rangle^m\) is unimodular,
so direct full-NS reflection integrality again forces

\[
                         C=2r,\qquad r^2=-1.          \tag{8.10}
\]

All \(\mu_i\) are even.  Negative pairing in (8.9) therefore forces
\(\sum_{i>j}\mu_i=0\), and the same wall has

\[
                              C.u_{j'}=-2
 \qquad(j'\ge j).                                    \tag{8.11}
\]

On an adapted half-boundary host, only the weaker tail bound
\(\sum_{i>j}\mu_i\le1\) survives, exactly as in the plane frame.

The canonical numerical doubled wall is

\[
 B_9=2u_9=4s+8f-2e_1-\cdots-2e_9.                   \tag{8.12}
\]

After later blowups away from it,

\[
 B_9^2=-4,\qquad K.B_9=2,\qquad p_a(B_9)=0,
 \qquad B_9.u_j=-2\quad(j\ge9).                      \tag{8.13}
\]

As before, this is an exact lattice-and-adjunction countermodel.  No
effective, reduced, irreducible representative is asserted.  The mirror is
that of the half-root \(u_9\); the sign of a reflecting root is immaterial,
while the chosen positive class \(B_9\) is the one with rational adjunction
color \(K.B_9=2\).

### 8.3 One wall and the first ruled mixed escape

Orient the added wall by \(x.B_9\ge0\).  In the old dual basis, its
coefficients are

\[
 (u_s.u_9,u_0.u_9,u_1.u_9,\ldots)
  =(2,4,3,6,5,4,3,2,1,0,-1,-1,\ldots).               \tag{8.14}
\]

The exhaustive extreme rays are the old rays with nonnegative coefficient
and

\[
                         u_i+p_i u_j
 \quad(p_i>0,\ j\ge9).                               \tag{8.15}
\]

For every pair ray,

\[
                   (u_i+p_i u_j)^2=u_i^2+p_i^2(10-j). \tag{8.16}
\]

All extreme rays are finite or ideal through \(j=10\).  At \(j=11\), use
\(u_s=f\), with \(u_s^2=0\) and \(p_s=2\), to get

\[
                    (f+2u_{11})^2=4(10-11)=-4.        \tag{8.17}
\]

Consequently

\[
 \boxed{\text{the \(B_9\)-augmented \(\mathbb F_2\) cone has finite volume}
        \quad\Longleftrightarrow\quad m=9,10.}        \tag{8.18}
\]

This independently confirms the same no-go on the only surviving direct
minimal ruled base: the old inverse diagonal cannot bound chain depth, and
the unresolved information is carried by new mixed directions.
