# Hostile audit: bounded Coxeter degree and the ideal-polygon tower

Date: 2026-08-09

## Verdict

The bounded-degree finiteness theorem is **valid after four repairs**.  Its
global form is existential; its form relative to a specified integral maximal
carrier is effective.  The proposed rational $H^2$ family is also valid, but
it is an ideal-$(N+2)$-gon tower, not a triangle family except at $N=1$.

## 1. Correct bounded-degree statement

Let $n\geq 2$, let $P\subset H^n$ be a finite-volume arithmetic Coxeter
polytope, and let $W=W(P)$.  Let $W_0$ range over maximal arithmetic
hyperbolic reflection groups containing $W$, and put

\[
 \deg_{\rm Cox}(P)=\min_{W_0\supseteq W}[W_0:W].
\]

ABSW prove both that such a maximal reflection overgroup exists and that
there are only finitely many conjugacy classes of maximal arithmetic
hyperbolic reflection groups.  Consequently, for every $N$, there are only
finitely many $P$, **up to ambient hyperbolic isometry**, with
\(\deg_{\rm Cox}(P)\leq N).

Indeed, choose a representative of each maximal conjugacy class.  Each
$W_0$ is finitely generated and therefore has only finitely many subgroups
of index at most $N$.  Restricting to reflection subgroups only decreases
this finite set.  For a fixed embedded reflection group $W$, the reflecting
hyperplanes cut $H^n$ into one $W$-orbit of chambers, so $W$ determines
its unmarked Coxeter polytope up to isometry.

The minimum is invariant under ambient conjugacy, but a minimizing carrier
$W_0$ need not be unique.  Thus the integer is canonical; a chosen carrier
or a choice among ties is extra decoration.

Reference: I. Agol, M. Belolipetsky, P. Storm, K. Whyte, *Finiteness of
arithmetic hyperbolic reflection groups*, arXiv:math/0612132.

## 2. Correct carrier-relative height estimate

Fix one maximal carrier $W_0=\langle s_1,\ldots,s_f\rangle$, a base chamber
$C_0$, a $k$-model $V$, and fixed normals $r_i$ to its walls.  Choose
a $W_0$-stable $\mathcal O_k$-lattice (a projective lattice and a
pseudobasis are sufficient).  If coordinates are used, record once and for
all the fixed fractional ideals, or a fixed denominator clearing them.

Suppose $W\leq W_0$ is a reflection subgroup of index $m\leq N$, and let
$P$ be the chamber of $W$ containing $C_0$.  Then $P$ is the connected
union of $m$ $W_0$-chambers.  Every boundary facet contains a panel of a
tile $wC_0$.  In the connected $m$-vertex tile graph one may choose that
tile with

\[
 |w|_{\{s_i\}}\leq m-1\leq N-1.
\]

The boundary reflection is $w s_iw^{-1}$, and a fixed normal to its mirror
is $\pm wr_i$.  For a submultiplicative norm at every embedding
$\sigma:k\hookrightarrow\mathbb C$, set

\[
 M=\max\!\left(1,\max_{\sigma,i}\|\sigma(s_i)\|\right),\qquad
 R=\max_{\sigma,i}\|\sigma(r_i)\|,
\]

and let $B$ bound the embedded bilinear forms.  Then every boundary root has
all-embedding height at most

\[
 R M^{N-1},
\]

and every boundary Gram entry has all-embedding height at most

\[
 C M^{2N-2},\qquad C=BR^2.
\]

This exponent is for roots and Gram entries.  It must not be asserted without
adjustment for the entries of the conjugated reflection matrix
$w s_iw^{-1}$.

After the one fixed denominator is cleared, the height-to-level norm argument
gives an explicit ideal level which injects the finite candidate set.  This
recovers $P$ from a **carrier-marked dessin that records the labelled tile
graph or its boundary roots**.  It is false, without an additional congruence
subgroup property, that the bare finite reduction image of $W$ determines an
arbitrary finite-index subgroup: noncongruence subgroups are the obstruction.

The constants above are explicit once $W_0$, its lattice, and its base roots
are supplied.  ABSW gives a finite global carrier list existentially, not an
enumerated effective list; therefore a uniform constant obtained by taking a
maximum over all carriers is presently existential unless that list is made
effective.

## 3. Exact audit of the rational $H^2$ tower

Use the Lorentz form

\[
 \langle x,y\rangle=-x_0y_0+x_1y_1+x_2y_2
\]

and, for $N\geq1$, the spacelike roots

\[
 \ell=(0,-1,0),\qquad u_N=(N,1,N),
\]

\[
 a_j=(1+j(j+1),\,2j+1,\,j(j+1)-1),
 \qquad 0\leq j<N.
\]

Direct expansion gives

\[
 \ell^2=u_N^2=a_j^2=1,
\]

\[
 \langle\ell,u_N\rangle=-1, 
 \langle\ell,a_j\rangle=-(2j+1), 
 \langle u_N,a_j\rangle=2j+1-2N,
\]

and

\[
 \boxed{\langle a_i,a_j\rangle=1-2(i-j)^2}.
\]

Thus the cyclically adjacent pairs

\[
 (\ell,a_0),\qquad (a_j,a_{j+1})\ (0\leq j<N-1),\qquad
 (a_{N-1},u_N),\qquad(u_N,\ell)
\]

all have pairing $-1$, while every nonadjacent pairing is strictly below
$-1$.  The ideal vertices are explicitly

\[
 v_\infty=(1,0,1),\qquad
 v_j=(j^2+1,\,2j,\,j^2-1),\quad0\leq j\leq N.
\]

They satisfy $v_j^2=0$, and consecutive roots annihilate the corresponding
$v_j$.  Hence these roots bound a genuine ideal $(N+2)$-gon $P_N$.
Its area is

\[
 \operatorname{area}(P_N)=((N+2)-2)\pi=N\pi.
\]

All reflections

\[
 s_r(x)=x-2\langle x,r\rangle r
\]

are integral matrices in $O(-1,1,1;\mathbb Z)$.  Since their polygon has
finite area,
$W_N=\langle s_\ell,s_{a_0},\ldots,s_{a_{N-1}},s_{u_N}\rangle$
is a finite-index subgroup of this arithmetic lattice and hence arithmetic.

There is also an exact index statement.  Put $u=u_1$,
$W_1=\langle s_\ell,s_{a_0},s_u\rangle$, and $T=s_us_\ell$.  Induction
gives

\[
 T^q a_0=a_{2q},\qquad T^q s_u(a_0)=a_{2q+1},
\]

\[
 T^q u=u_{2q+1},\qquad T^q\ell=-u_{2q}.
\]

Signs do not change mirrors.  Therefore every boundary reflection of $P_N$
belongs to $W_1$, so $W_N\leq W_1$.  The covolume ratio now proves

\[
 [W_1:W_N]
 =\frac{\operatorname{area}(P_N)}{\operatorname{area}(P_1)}=N.
\]

This family demonstrates why the degree cutoff is substantive: a single
arithmetic carrier has finite-index reflection subgroups of unbounded index
and unbounded wall/Gram height.

## 4. Required wording repairs

1. Say “up to ambient isometry,” not finitely many embedded copies.
2. Say “carrier-marked tile/boundary-root congruence dessin,” not a bare
   congruence image.
3. Attach $CM^{2N-2}$ to Gram height, and include fixed denominator or
   projective-lattice data.
4. Label global finiteness existential and fixed-carrier reconstruction
   effective.
