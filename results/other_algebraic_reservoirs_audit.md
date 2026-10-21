# Other algebraic reservoirs: seven exact rational nef chambers and a sharp design ceiling

## Outcome

The broader algebro-geometric search now gives seven strict divisor-cone
examples above `H9`.  The first two are literal cases of a February 2026
cone theorem.  A controlled off-fibre continuation, followed by explicit
quartic and sextic branches, gives exact surface chambers through `H16`.

Let $C\subset \mathbf P^2$ be a smooth cubic, let $p\in C$ be a flex, and
start with eleven successive infinitely-near points along the strict
transform of $C$ above $p$.  For $n=10,11$ truncate that chain.  Higher
members leave $C$ after the eleventh point and follow explicitly chosen
branches on the exceptional tail.  Then

\[
n=10,11,12,13,14,15,16,
\]

the complete nef cone of the resulting rational surface $X_n$, after
projectivizing the positive cone, is Vinberg's classical
odd-unimodular chamber

\[
P_n\subset\mathbb H^n.
\]

All seven are finite-volume and noncompact.  Exact cone enumeration gives

| surface | lattice | walls | finite vertices | cusps | volume |
|---|---|---:|---:|---:|---|
| $X_{10}$ | $I_{1,10}$ | 12 | 17 | 2 | finite, noncompact |
| $X_{11}$ | $I_{1,11}$ | 13 | 26 | 2 | finite, noncompact |
| $X_{12}$ | $I_{1,12}$ | 14 | 35 | 2 | finite, noncompact |
| $X_{13}$ | $I_{1,13}$ | 15 | 43 | 3 | finite, noncompact |
| $X_{14}$ | $I_{1,14}$ | 17 | 89 | 5 | finite, noncompact |
| $X_{15}$ | $I_{1,15}$ | 18 | 139 | 6 | finite, noncompact |
| $X_{16}$ | $I_{1,16}$ | 20 | 313 | 12 | finite, noncompact |

The chambers $P_{10}$ and $P_{11}$ are classical.  Moreover, Galindo,
Monserrat, Moreno-Ávila, and Pérez-Callejo's
[*Nef divisors of surfaces given by pencils at
infinity*](https://arxiv.org/abs/2602.15703), Theorems 3 and 7, literally
covers the blow-up construction: the first nine points resolve the pencil
$\langle C,3L\rangle$, and the tenth and eleventh are a free chain of length
one or two on the selected fibre $C$.  The contribution of this audit is the
exact Coxeter identification, independent finite-volume certificate, and a
newly extracted **sharp $H^{11}$ ceiling for this entire published factory**.
No worldwide priority claim is made for any of those deductions.

The surfaces $X_{12}$ through $X_{16}$ lie outside that published factory.
Their effective-cone equalities are proved below by exact rational
polyhedral inclusions: proximity and Bézout put every irreducible curve in a
small necessary cone, and every extreme ray of that cone decomposes
nonnegatively into the displayed wall curves.  The abstract polytopes
$P_{12},\ldots,P_{16}$ are classical; the bounded search found no prior
source for these particular surface realizations, but worldwide priority
remains open.  A divisor-restriction argument proves that the same standard
marked construction cannot realize the required $P_{17}$ wall set.

The three novelty questions must be kept separate:

| question | status |
|---|---|
| new abstract hyperbolic Coxeter polytope? | **no**: every $P_{10},\ldots,P_{16}$ is Vinberg's 1972 odd-unimodular chamber |
| first example in a previously empty dimension? | **no**: these dimensions were already covered abstractly |
| new algebro-geometric realization? | **proved here** for the displayed $X_{12},\ldots,X_{16}$; worldwide priority remains open |

The audit also repairs a major omission in the abstract lineage ledger.
Borcherds' Reed--Muller construction gives a known 960-wall finite-volume
noncompact chamber in $\mathbb H^{17}$.  It does not give a new example now,
and a discriminant-capacity theorem below excludes its full wall set from
being a $K3^{[3]}$-type movable/divisor chamber.

The exact finite calculations are frozen in two independent layers.  The
$P_{10}$--$P_{13}$ and Reed--Muller calculations are replayed by
[`certify_other_algebraic_reservoirs.py`](../scripts/certify_other_algebraic_reservoirs.py)
and its
[`certificate`](../data/enumerations/other_algebraic_reservoirs_certificate.json),
whose payload SHA-256 is

```text
a0f43f751eaf9c75d687e919e395adb6a400553e2d87635c9d3809c8e98cec6e
```

The $P_{12}$--$P_{16}$ cone decompositions, volume audits, and the
marking-independent $P_{17}$--$P_{19}$ rational-carrier obstructions are replayed by
[`certify_rational_i1n_inverse_design.py`](../scripts/certify_rational_i1n_inverse_design.py),
with full theorem notes in the
[`inverse-design report`](rational_i1n_inverse_design.md) and exact
[`certificate`](../data/enumerations/rational_i1n_inverse_design.json).
Its payload SHA-256 is

```text
17ab489b76eba0042c10d3dddce32466020595470ce3d64223e169c9266a03b7
```

## 1. A general rational-polyhedral surface factory

### Theorem 1.1

Let $X$ be a smooth projective surface and put
$N=\operatorname{NS}(X)/\text{torsion}$, of rank $\rho$.  Suppose:

1. $\operatorname{Nef}(X)$ is rational polyhedral; and
2. every irreducible negative curve whose orthogonal hyperplane is a genuine
   facet of the projectivized nef cone has square $-1$ or $-2$.

Then

\[
\mathbf P\bigl(\operatorname{Nef}(X)\cap\mathcal C_X^+\bigr)
\subset\mathbb H^{\rho-1}
\]

is a finite-volume integral hyperbolic Coxeter polytope.  It is compact if
and only if no extremal nef ray is isotropic.

### Proof

The Hodge index theorem gives $N_{\mathbb R}$ signature
$(1,\rho-1)$.  Distinct irreducible curves $A,B$ have
$k=A\cdot B\ge0$.  Writing $a=-A^2$ and $b=-B^2$, the normalized
wall product is $k/\sqrt{ab}$.  Because $a,b\in\{1,2\}$, every pair is
one of the following:

| curve squares | intersection | wall relation |
|---|---:|---|
| any | 0 | $m=2$ |
| $(-2,-2)$ | 1 | $m=3$ |
| $(-1,-2)$ | 1 | $m=4$ |
| $(-1,-1)$ | 1 | parallel |
| $(-2,-2)$ | 2 | parallel |
| any larger ratio |  | ultraparallel |

Thus every dihedral angle is a Coxeter angle.  Moreover

\[
s_A(x)=x-\frac{2(x\cdot A)}{A^2}A
\]

is integral for $A^2\in\{-1,-2\}$.  A nef class has nonnegative square, so
every extremal ray of the rational-polyhedral nef cone lies in the closed
future cone.  Its Klein-model section is therefore a finite Euclidean
polytope all of whose vertices are finite or ideal.  Such a hyperbolic
polytope has finite volume.  It is compact exactly when every extremal ray
has positive square.  Strict convexity of the Lorentz light cone shows that
any isotropic point of a finitely generated nef cone already lies on an
isotropic extremal ray.  This proves the last assertion.

Null or positive-square curve classes do not define hyperbolic facets:
their orthogonals are respectively tangent to, or disjoint from, the
positive cone.  Thus the theorem does not silently omit another kind of
interior wall.

## 2. The flex-cubic surfaces $X_{10}$ through $X_{16}$

### The published pencil-at-infinity cone theorem

Let $L$ be the flex line.  Then the smooth cubic $C$ has one place at
infinity relative to $\mathbf P^2\setminus L$: it meets $L$ only at $p$,
with local intersection multiplicity three.  Consider the pencil

\[
\mathcal P_C=\langle C,3L\rangle.
\]

Theorem 3 of Galindo--Monserrat--Moreno-Ávila--Pérez-Callejo (quoting the
earlier theorem of Campillo--Piltant--Reguera) says that after resolving all
basepoints of this pencil, the cone of curves is generated by the strict
infinity line and the strict exceptional curves.  Their Theorem 7 permits
one or more free chains along selected fibres of the pencil and says exactly
that

\[
\operatorname{NE}(\widetilde X)
=\operatorname{Cone}
\bigl(\widetilde L,\widetilde E_i,\widetilde C_j\bigr),
\]

while also giving an explicit finite generating set for the dual nef cone.
This is an equality of complete cones, not merely a list of visible curves.

In the smooth degree-$d$ case the base cluster can be read off without a
black box.  Every centre lies on the smooth strict transform of $C$, hence
has multiplicity one and is free.  A resolved pencil fibre has square zero,
so

\[
0=(dH-E_1-\cdots-E_n)^2=d^2-n;
\qquad n=d^2.
\tag{2.1}
\]

The infinity line and $C$ have intersection multiplicity $d$ at $p$, so the
line follows precisely the first $d$ centres.  For a cubic, therefore, the
base cluster consists of **nine successive free points**, and
$\widetilde L=H-E_1-E_2-E_3$.  Choosing $C$ itself as a selected fibre and
adding a free chain of length one or two gives exactly $X_{10}$ and
$X_{11}$.  Thus the 2026 theorem supplies the requested algebro-geometric
completeness proof; the calculation below identifies its cone with the
classical Coxeter chamber.

Let $H$ be the pullback of a line and let $E_i$ denote the total-transform
basis for the successive exceptional divisors.  Then

\[
\operatorname{NS}(X_n)=
\mathbb ZH\oplus\bigoplus_{i=1}^n\mathbb ZE_i,
\qquad
H^2=1,\quad E_i^2=-1,
\]

so the lattice is $I_{1,n}$.  The following irreducible curves exist:

\[
\begin{aligned}
R_i&=E_i-E_{i+1} &&(1\le i<n),\\
R_n&=E_n,\\
L&=H-E_1-E_2-E_3,\\
\widetilde C&=3H-E_1-\cdots-E_n=-K_{X_n}.
\end{aligned}
\tag{2.2}
\]

The $R_i$ are the strict exceptional chain.  The flex line has contact
order three with $C$, hence passes through precisely the first three
centres and gives $L$.  The last class is the smooth strict transform of
the cubic.  Their squares are

\[
R_i^2=L^2=-2\ (i<n),\qquad R_n^2=-1,
\qquad \widetilde C^2=9-n.
\]

Consequently the norm histograms are

\[
n=10:\quad 10(-2)+2(-1),
\qquad
n=11:\quad 12(-2)+1(-1).
\]

The nonzero intersections are a chain

\[
R_1-R_2-\cdots-R_n-\widetilde C,
\]

together with $L\cdot R_3=1$; every displayed adjacency has intersection
one.  For $n=10,11$, the exact Gram matrices are respectively the standard
simple-root Grams of the full odd-unimodular reflection chambers
$P_{10}$ and $P_{11}$.  Vinberg's construction of the odd-unimodular
groups is the primary lineage source; the maintained census identifies the
same full chambers as Scharlau--Kirschmer `List10` and `List11`, entry 1.
See Vinberg's
[*Groups of units of certain quadratic forms*](https://doi.org/10.1070/SM1972v016n01ABEH001346).
That 1972 paper is also the primary source for $P_{12}$ through $P_{17}$;
the later Vinberg--Kaplinskaja continuation is needed only in dimensions 18
and 19.

### An independent completeness inequality

Let $D\ne C$ be an irreducible plane curve of degree $d$, and let
$m_1,\ldots,m_n$ be its multiplicities at the successive cluster.  Its
strict transform has class $dH-\sum m_iE_i$.  Since $D$ and $C$ are
distinct effective curves,

\[
\widetilde C\cdot\widetilde D=3d-\sum_i m_i\ge0.
\tag{2.3}
\]

The arithmetic genus of the integral strict transform is nonnegative, so

\[
(d-1)(d-2)-\sum_i m_i(m_i-1)\ge0.
\tag{2.4}
\]

Adding (2.3) and (2.4) gives the decisive bound

\[
\widetilde D^2
=d^2-\sum_i m_i^2
\ge 3d-2-\sum_i m_i
=\widetilde C\cdot\widetilde D-2
\ge-2.
\tag{2.5}
\]

The exceptional curves are already listed, and
$\widetilde C^2=-1,-2$.  Therefore **every** irreducible negative curve on
$X_{10}$ or $X_{11}$ has square $-1$ or $-2$.

### Coxeter identification and nef-cone equality

The curves (2.2) orient the full all-root chamber $P_n$ toward any ample
class.  Any further negative curve is an integral norm-$-1$ or norm-$-2$
root.  Its mirror cannot cross the interior of the full root chamber, and
the ample class selects its nonnegative side.  An effective curve of
nonnegative square lies in the closure of the same future cone, so it also
has nonnegative intersection with every point of $P_n$.  Conversely, the
curves (2.2) are effective and hence every nef class satisfies their
inequalities.  Thus

\[
\boxed{\operatorname{Nef}(X_n)\cap\mathcal C^+_{X_n}=P_n,
\qquad n=10,11.}
\tag{2.6}
\]

The certificate independently enumerates the extreme rays of both rational
cones.  It finds `17+2+0` and `26+2+0` finite, ideal, and hyperideal rays.
This is a direct finite-volume/noncompact audit, not an inference from the
facet count.

For the single cubic-fibre chain, the construction stops sharply at
$n=11$: at $n=12$ the irreducible cubic wall has square $-3$, and its
reflection is not integral in the unimodular divisor lattice.  The next
subsection shows that this is actually the sharp ceiling for the full
one-place-at-infinity factory of Theorem 7, not just this one chain.

### Theorem 2.1: sharp ceiling for the published one-place factory

Within the Galindo--Monserrat--Moreno-Ávila--Pérez-Callejo construction,
require every negative Mori generator to define an integral divisor-lattice
reflection.  Then the projectivized nef cone has hyperbolic dimension at
most eleven.  Equality forces the smooth cubic base pencil, exactly one
selected fibre, and a free chain of length two.  It is the $P_{11}$ case
above.

Indeed, first suppose that the one-place curve $C$ is singular at infinity.
The embedded resolution of a singular plane branch contains a satellite
point.  One quick proof uses the proximity equality: if every point were
free, the only point proximate to $p_i$ would be $p_{i+1}$, so the branch
multiplicities would satisfy $m_i=m_{i+1}$ all the way to the eventual value
one; this would force $m_1=1$, contrary to singularity.  The pencil base
cluster contains the embedded-resolution cluster.  If $p_j$ is proximate
both to its predecessor and to an older point $p_i$, then at least two later
centres are proximate to $p_i$.  Hence the
strict exceptional curve

\[
\widetilde E_i=E_i^*-\sum_{q\to p_i}E_q^*
\]

has square at most $-3$.  It is an irreducible negative extremal curve in
the exact cone of Theorem 7.  Its class is primitive in the unimodular
plane-blow-up lattice, so it has divisibility one and its reflection cannot
be integral.  Thus the integral-root lane forces $C$ to be smooth at
infinity.

Now let $C$ be smooth of degree $d$, let $k$ selected fibre chains be
attached to the terminal base exceptional, and write their positive lengths
as $s_1,\ldots,s_k$.  The complete list supplied by Theorem 7 contains
curves with squares

\[
\widetilde L^2=1-d,\qquad
\widetilde E_n^2=-1-k,\qquad
\widetilde C_j^2=-s_j.
\tag{2.7}
\]

All three types are primitive.  A primitive vector in a unimodular lattice
defines an integral reflection only when its negative square divides two.
Consequently

\[
d\le3,\qquad k\le1,\qquad s_1\le2.
\tag{2.8}
\]

(The degree-one case is harmless and low-dimensional; its infinity line is
isotropic rather than a hyperbolic wall.)  Since the number of blowups, and
hence the hyperbolic dimension, is

\[
d^2+\sum_j s_j,
\]

the maximum is $3^2+2=11$.  Conversely the cubic with one length-two chain
has only square-$-1$ and square-$-2$ negative generators and is precisely
the certified $P_{11}$ chamber.  This proves sharpness.

### Theorem 2.2: leave the fibre and obtain $P_{12}$ and $P_{13}$

The preceding ceiling is a ceiling for the configuration in the 2026 cone
theorem, not for rational surfaces.  Starting with $X_{11}$, choose $p_{12}$
on the terminal exceptional $E_{11}$ away from both $\widetilde C$ and the
older exceptional curve.  If desired, choose $p_{13}$ freely on the new
$E_{12}$, again away from the older exceptional.  Let $X_n$ be the resulting
surface for $n=12$ or $13$.  Put

\[
R_i=E_i-E_{i+1}\quad(1\le i<n),\qquad R_n=E_n,
\]

and retain

\[
L=H-E_1-E_2-E_3,
\qquad C=3H-E_1-\cdots-E_{11}.
\tag{2.9}
\]

Then

\[
\boxed{
\operatorname{NE}(X_n)
=\operatorname{Cone}(R_1,\ldots,R_n,L,C),
\qquad n=12,13.}
\tag{2.10}
\]

In particular, the dual nef cones projectivize to the complete classical
odd-unimodular chambers

\[
P_{12}\subset\mathbb H^{12},
\qquad P_{13}\subset\mathbb H^{13}.
\tag{2.11}
\]

Here is the completeness proof.  Let $D$ be an irreducible plane curve
different from $L$ and $C$, with strict-transform class
$dH-\sum_{i=1}^n m_iE_i$, and set $A_j=\sum_{i=1}^j m_i$.  The centres form
a free chain, so the proximity inequalities give

\[
m_1\ge m_2\ge\cdots\ge m_n\ge0.
\tag{2.12}
\]

Bézout with $L$ and $C$ gives

\[
A_3\le d,\qquad A_{11}\le3d.
\tag{2.13}
\]

Let $\ell=n-11\le2$ and define

\[
y=\max\left(0,\frac{A_n-3d}{2}\right),\qquad
x=d-3y,
\]

\[
\mu_j=x\min(j,3)+y\min(j,11)-A_j.
\tag{2.14}
\]

Since $A_n-A_{11}\le\ell m_{11}$, one has
$2y\le2m_{11}$ and hence $0\le y\le m_{11}\le d/3$; thus $x\ge0$.
The remaining coefficients are nonnegative as follows.

For $j=1,2,3$, the three inequalities

\[
m_1+2m_{11}\le A_3\le d,
\quad A_2+4m_{11}\le2A_3\le2d,
\quad A_3+6m_{11}\le3A_3\le3d
\]

give $\mu_j\ge0$.  For $4\le j\le8$,

\[
A_j+(9-j)y
\le d+(j-3)m_4+(9-j)m_{11}
\le d+6m_4\le3d,
\]

so again $\mu_j\ge0$.  For $9\le j\le11$, use
$\mu_j=3d+(j-9)y-A_j\ge3d-A_{11}$; on the exceptional tail use
$\mu_j=3d+2y-A_j\ge3d+2y-A_n$.  Therefore

\[
[\widetilde D]=x[L]+y[C]+\sum_{j=1}^n\mu_j[R_j]
\tag{2.15}
\]

is a nonnegative decomposition.  The exceptional prime curves, $L$, and
$C$ are already among the displayed generators, proving (2.10) for every
irreducible curve and hence for the whole cone.

Every $R_i$ with $i<n$, as well as $L$ and $C$, has square $-2$; $R_n$
has square $-1$.  Their exact Gram matrices are the standard simple-root
Grams of $P_{12}$ and $P_{13}$.  Exact dual-cone enumeration gives

| surface | walls | finite vertices | cusps | hyperideal vertices |
|---|---:|---:|---:|---:|
| $X_{12}$ | 14 | 35 | 2 | 0 |
| $X_{13}$ | 15 | 43 | 3 | 0 |

Thus both volumes are finite and both chambers are noncompact.  No inference
from wall count is being used.

Vinberg's 1972 odd-unimodular chambers are the abstract lineage, so
these are not new hyperbolic polytopes or new arithmetic groups.  The 2026
pencil theorem does not cover the off-fibre blowups.  No matching published
surface realization appeared in the bounded source and phrase search, but
that is not a worldwide priority proof; the safe status is **newly proved in
this project, external priority open**.

### Theorem 2.3: quartic and sextic branches give $P_{14}$ through $P_{16}$

The off-fibre construction continues, but the new Vinberg walls must now be
real curves rather than merely lattice vectors.  Use coordinates
$[X:Y:Z]$, put $p=[0:1:0]$, and define

\[
L=(Z=0),\qquad
C=(F=0),\qquad F=Y^2Z-X^3-Z^3,
\]

\[
Q=(Z^4+F(X+Z)=0),
\qquad
S=(F^2+Z^2Q=0).
\tag{2.16}
\]

The cubic $C$ is smooth and $L|_C=3p$.  The quartic $Q$ has an ordinary
node at $p$, with one branch transverse to $C$ and the other having contact
eleven with $C$.  Equivalently, $Q|_C=Z^4$, so the total intersection is
$12p=1p+11p$.  Follow the latter branch through $p_{11}$ and then through
$p_{12},\ldots,p_{15}$.  The sextic satisfies

\[
S|_C=Z^6,\qquad S|_Q=F^2.
\]

Its selected branch follows $Q$ through $p_{15}$ and separates there; take
$p_{16}$ on that branch.  The multiplicity sequences along the resulting
free chain are

\[
L:(1^3),\qquad C:(1^{11}),\qquad
Q:(2,1^{14}),\qquad S:(2^7,1^9).
\tag{2.17}
\]

An exact absolute-factorization replay over characteristic zero returns one
nonconstant absolute factor, with exponent one, for each of $F,Q,S$.
Thus the required curves are integral, not formal divisor classes.

For $X_{14}$ and $X_{15}$ truncate the cluster after the corresponding
point and use $R_1,\ldots,R_n,L,C,Q$.  For $X_{16}$ use
$R_1,\ldots,R_{16},L,C,Q,S$.  Then

\[
\boxed{
\begin{aligned}
\operatorname{NE}(X_n)
 &=\operatorname{Cone}(R_1,\ldots,R_n,L,C,Q),&&n=14,15,\\
\operatorname{NE}(X_{16})
 &=\operatorname{Cone}(R_1,\ldots,R_{16},L,C,Q,S).
\end{aligned}}
\tag{2.18}
\]

The completeness certificate is finite and exact.  For a plane curve
$D=dH-\sum m_iE_i$, put $A_j=\sum_{i\le j}m_i$.  Proximity and Bézout force

\[
d\ge0,\quad m_1\ge\cdots\ge m_n\ge0,\quad
A_3\le d,\quad A_{11}\le3d,\quad
A_{\min(n,15)}+m_1\le4d,
\tag{2.19}
\]

For $n=16$ there is additionally

\[
A_{16}+A_7\le6d.
\tag{2.20}
\]

The rational cones cut out by these inequalities have respectively

\[
94,\qquad145,\qquad325
\]

primitive extreme rays.  Every ray passes an exact nonnegative-decomposition
test in the visible curve cone.  Concretely, for $n=14,15$ the certificate
solves

\[
D=xL+yC+zQ+\sum_i\mu_iR_i,qquad
x=d-3y-4z,
\]

\[
\mu_j=x\min(j,3)+y\min(j,11)+z(j+1)-A_j,
\tag{2.21}
\]

with every coefficient nonnegative.  For $n=16$ it adds $wS$, replaces
$x$ by $d-3y-4z-6w$, and adds

\[
w\,\begin{cases}2j,&j\le7,\\j+7,&j\ge8\end{cases}
\]

to $\mu_j$.  There are zero failed rays in all three dimensions.  Since
every actual irreducible plane curve satisfies (2.19)--(2.20), and the
exceptional curves are already listed, this proves (2.18); it is not a
bounded search for negative curves.

The wall Grams are exactly Vinberg's full odd-unimodular
chambers.  The complete volume data are

| surface | chamber | walls | finite vertices | cusps | hyperideal vertices |
|---|---|---:|---:|---:|---:|
| $X_{14}$ | $P_{14}$ | 17 | 89 | 5 | 0 |
| $X_{15}$ | $P_{15}$ | 18 | 139 | 6 | 0 |
| $X_{16}$ | $P_{16}$ | 20 | 313 | 12 | 0 |

Thus all three are finite-volume and noncompact integral Coxeter nef
chambers.

### Theorem 2.3a: non-isotrivial carriers and geometric cusp cycles

The positive construction is not rigid.  Replacing the flex cubic by

\[
C_a:Y^2Z=X^3+aXZ^2+Z^3
\]

and putting `Q_a=Z^4+(X+Z)C_a`, `S_a=C_a^2+Z^2Q_a` preserves the contact-11
and contact-15 branch data, and hence the same complete cone proofs, on a
nonempty open parameter set.  At `a=1` the unique genus-one Mori facet has
`j=6912/31`, while it has `j=0` on the archived `a=0` surface.  Thus the two
explicit `P16` carriers are nonisomorphic, and nonconstancy of
`j(a)=6912a^3/(4a^3+27)` produces infinitely many pairwise nonisomorphic
complete carriers of every `P12`--`P16`.  The exact cluster and
irreducibility replay is [the `j`-moduli theorem](rational_i1n_j_moduli_carriers.md).

An independent exact dual-cone calculation gives a classical-invariant
interpretation of finite volume.  All 28 ideal rays are primitive nef
isotropic divisor classes.  If a cusp has `c` connected affine components,
its complete orthogonal negative-curve configuration has `n+c-1` nodes and
rank `n-1`; the primitive positive null cycle of every component equals the
cusp class itself.  Sixteen cusps have at least two such disjoint cycles,
which give a base-point-free pencil.  Five line rulings and the `P13`
quartic pencil give the first 22.  Exact complete-ideal jet matrices and
coprime second plane sections resolve five further single-component systems:
genus five on `P14`--`P16`, and genus ten and eleven on `P16`.  Thus 27 of
the 28 cusps move in their primitive classes.  The last `P16` genus-three
class has exact first moving multiple `|4F|`: the kernel dimensions for
`F,2F,3F,4F` are `1,1,1,2`, and the coprime degree-16 member proves
base-point-freeness.  Hence all 28 cusp rays are semiample.  See the exact
[cusp-fibration theorem](rational_i1n_cusp_fibrations.md).

### Theorem 2.4: no smooth rational surface carries the complete $P_{17}$--$P_{19}$ facet cones

The next full odd-unimodular chambers require, in their standard
$B_n$ markings, the additional wall classes

\[
A_n=4H-E_1-\cdots-E_{\min(n,18)},\qquad n=17,18,19.
\tag{2.22}
\]

Thus the literal class follows 17 centres for $P_{17}$ and 18 centres for
$P_{18}$ and $P_{19}$.

It cannot be represented by the required irreducible curve on the cluster
above.  Indeed, the flex condition gives
$\mathcal O_C(1)=\mathcal O_C(3p)$.  The quartic branch leaves $C$ after
$p_{11}$, so a putative $A_n$ following the full free chain has

\[
A_n|_C=11p+r
\]

for a residual point $r\ne p$.  On the other hand,

\[
\mathcal O_C(A_n)=\mathcal O_C(4)=\mathcal O_C(12p),
\]

and therefore $\mathcal O_C(r)=\mathcal O_C(p)$.  The degree-one Abel map
on the smooth locus of the integral cubic is injective, forcing $r=p$, a
contradiction.  The literal class changes at $P_{18}$, but the same
restriction obstruction persists in all three packets.

The Abel calculation first proves that $P_{16}$ is a sharp endpoint for the
standard flex-cubic/free-chain marking.  The marking hypothesis can now be
removed for all three subsequent classical chambers.  If a smooth
projective rational surface had complete Mori cone generated by the $P_n$
facet classes, every $K$-negative extremal contraction would contract one of
the current facet images.  Enumerating the exact update
$f^*f_*D=D+(D\cdot e)e$ gives respectively 8, 24, and 120
$\mathbf P^2$ endpoints for $n=17,18,19$, and the same numbers of minimal
$\mathbb F_2$ endpoints.  In each dimension their canonical
reconstructions coincide and form exactly the full weighted-diagram orbit
of the standard plane-blowup canonical class; the group orders are 8, 24,
and 120.

After conjugating any putative marking by the corresponding integral
diagram isometry, adjunction makes its $B_n$ chain smooth rational.
Castelnuovo contraction recovers the standard free chain, so the Abel
contradiction applies.  Hence **no smooth projective rational surface has
the complete $P_{17}$, $P_{18}$, or $P_{19}$ Mori/nef carrier**.  This is
not an abstract-polytope nonexistence theorem: the Vinberg and
Vinberg--Kaplinskaja finite-volume chambers are classical.  The exact
replays are
[`certify_rational_p17_mmp.py`](../scripts/certify_rational_p17_mmp.py) and
[`certify_rational_p17_canonical_cusps.py`](../scripts/certify_rational_p17_canonical_cusps.py),
together with the
[`P18`/`P19` MMP theorem](rational_p18_p19_mmp.md).

## 3. The missing known $\mathbb H^{17}/960$ chamber

### Exact lineage

Borcherds' [*Reflection groups of Lorentzian
lattices*](https://arxiv.org/abs/math/9909123), Duke Math. J. **104**
(2000), Section 12, case `N=2`, pp. 342--344 (arXiv PDF pp. 24--26), treats

\[
L=II_{1,17}(2^{+6})=(18,6,0),
\qquad \operatorname{sign}(L)=(1,17),
\qquad |\det L|=2^6.
\]

He writes $L=K\oplus II_{1,1}$, where $K$ is obtained by Construction A
from the length-16 first-order Reed--Muller code.  The complete fundamental
domain has

\[
64\text{ norm-}(-2)\text{ walls}
+896\text{ norm-}(-4)\text{ walls}=960\text{ walls}.
\tag{3.1}
\]

The 896 count is the displayed MacWilliams calculation in that section:
the dual code has 448 words of weight 6 and 448 of weight 10.  Borcherds
states that this case has an arithmetic reflection group and exhibits null
vectors $z,z'$ in the domain.

There is also a later literal cusp audit.  Alexeev's
[*Reflective hyperbolic 2-elementary lattices, K3 surfaces and hyperkahler
manifolds*](https://arxiv.org/abs/2209.09110), Theorem
`thm:table` and subsection `(18,6,0)`, attributes the diagram to Borcherds
and computes exactly

\[
90,897,634
\]

maximal parabolic subdiagrams, in 28 diagram-symmetry orbits.  Therefore the
volume status is not guessed:

\[
\boxed{\mathbb H^{17}/960\text{ is finite-volume and noncompact}.}
\]

It was a new example in Borcherds' 2000 paper.  It is a **known** abstract
example now and was simply missing from this repository's ledger.

### Algebraic-carrier obstructions

Borcherds already notes that $L$ cannot be the Picard lattice of a complex
K3 surface: its discriminant 2-rank is 6, larger than its codimension 4 in
the K3 lattice.

It also cannot be the full movable/divisor chamber of a $K3^{[2]}$-type
variety.  Markman's [Theorem 1.2](https://arxiv.org/abs/0912.4981) says that
a primitive negative class whose multiple is an irreducible divisor on a
$K3^{[n]}$-type variety has square $-2$ or $2-2n$.  For $n=2$ both
values are $-2$, excluding the 896 square-$-4$ walls.

The apparently better $K3^{[3]}$ match also fails for the **full** wall
set, for a subtler reason.  Here square $-4$, ambient divisibility 2 is a
permitted prime-exceptional numerical type.  Let

\[
\Lambda_3=3U\oplus2E_8(-1)\oplus\langle-4\rangle,
\qquad A_{\Lambda_3}\simeq\mathbb Z/4.
\]

Suppose $L\hookrightarrow\Lambda_3$ were primitive and let $T=L^\perp$,
of rank 5.  If the gluing projection $H_L\subset A_L=(\mathbb F_2)^6$
has order $2^d$, the determinant formula gives

\[
|A_T|=\frac{|A_{\Lambda_3}|\,|H_L|^2}{|A_L|}=2^{2d-4}.
\tag{3.2}
\]

Injection into $A_T$ forces $d\ge4$.  The rank-five complement excludes
$d=6$.  If $d=4$, then $H_T=A_T\simeq(\mathbb Z/2)^4$, and the residual
discriminant group is the 2-elementary order-four orthogonal complement of
$H_L$, contradicting $A_{\Lambda_3}\simeq\mathbb Z/4$.  Hence

\[
\dim_{\mathbb F_2}H_L=5.
\tag{3.3}
\]

A primitive norm-$-4$ root $r\in L$ has divisibility exactly 2 in $L$,
and its class $r/2\in A_L$ is anisotropic.  It retains ambient divisibility
2 precisely when

\[
[r/2]\in H_L^\perp.
\tag{3.4}
\]

But (3.3) makes $H_L^\perp$ one-dimensional.  The Reed--Muller certificate
shows more: modulo the 32-word code, the 896 weight-6/10 dual words fill all
28 anisotropic classes of $A_L$, exactly 32 roots per class.  Thus at most
32 long roots can satisfy (3.4), and at least

\[
896-32=864
\]

long roots acquire ambient divisibility 1.  Their reflections are not
integral on $\Lambda_3$, and Markman's theorem excludes them as
prime-exceptional divisor classes.  Therefore the complete 960-wall chamber
cannot be a $K3^{[3]}$-type movable/divisor chamber.  This statement is
scoped to the full wall set and this deformation type; it does not exclude
other hyperkahler types or a non-divisorial lattice occurrence.

## 4. Three sharp surface-reservoir filters

### 4.1 Coble surfaces stop at $\mathbb H^9$

For a complex Coble surface $S$, the raw Picard group can grow with the
number of bi-anticanonical boundary components.  The reflective object is
nevertheless the Coble--Mukai lattice, the orthogonal complement of those
boundary components after adjoining their half-classes.  Kondō's
Proposition `CME10` in [*Coble surfaces with finite automorphism
group*](https://arxiv.org/abs/2102.00591) proves

\[
\operatorname{CM}(S)\simeq E_{10},\qquad
\operatorname{sign}\operatorname{CM}(S)=(1,9)
\]

over $\mathbb C$.  His Proposition `NefCone` identifies the nef intersection
with the corresponding root chamber.  Thus the natural Coble--Mukai factory
has a sharp `H9` ceiling.  The surface automorphism group preserves the
negative bi-anticanonical boundary and acts on its orthogonal Coble--Mukai
lattice, so increasing the raw Picard rank by adding boundary components
does not enlarge this **standard reflective quotient**.

This is only a no-go for the classical Coble--Mukai mechanism.  Kondō's
theorem does **not** exclude a different Coxeter chamber built directly in
the larger raw Picard lattice, and no claim about the divisibility of every
raw boundary component is made here.

### 4.2 The canonical Cremona simplex also stops at $\mathbb H^9$

For the blow-up of $\mathbf P^2$ at $n$ marked points, McMullen's
[*Dynamics on blowups of the projective
plane*](https://www.numdam.org/article/PMIHES_2007__105__49_0.pdf),
Sections 2--3, uses the roots

\[
\alpha_0=H-E_1-E_2-E_3,
\qquad \alpha_i=E_i-E_{i+1},
\]

whose diagram is $T_{2,3,n-3}$ in $K_X^\perp$.  At $n=10$ this is the
finite-volume `H9` simplex $T_{2,3,7}$.  For every $n\ge11$, deleting the
long-arm terminal leaves that Lorentzian $T_{2,3,7}$ as a proper vertex
subdiagram.  The canonical simplex therefore has a hyperideal vertex and
infinite volume.

This does not conflict with the flex-cubic examples: those use the **full**
odd Néron--Severi lattice and add the effective anticanonical wall; they are
not the canonical $K_X^\perp$ Cremona simplex.

### 4.3 Smooth toric surfaces fail above $\mathbb H^9$

Let $X$ be a smooth complete complex toric surface with $r$ boundary
divisors $D_i$.  Then $\rho=r-2$, $K_X^2=12-r$, and

\[
K_X^2=\left(\sum_iD_i\right)^2=\sum_iD_i^2+2r,
\qquad
\sum_iD_i^2=12-3r.
\tag{4.1}
\]

If $\rho\ge11$, then $r\ge13$.  Were every negative invariant curve of
square at least $-2$, (4.1) would give simultaneously

\[
\sum_iD_i^2\ge-2r
\quad\hbox{and}\quad
\sum_iD_i^2=12-3r<-2r,
\]

a contradiction.  Hence some invariant negative curve has square at most
$-3$.  Such a curve is a genuine Mori/nef facet.  Its primitive class lies
in the unimodular divisor lattice of a rational surface, so its reflection
is not integral.  Smooth toric rational surfaces therefore cannot satisfy
Theorem 1.1 in Picard rank at least 11.

## 5. Ranked reservoir inventory

| rank | reservoir | hard status | exact reason |
|---:|---|---|---|
| 1 | flex-cubic/quartic/sextic branch factory | **non-isotrivial positive through `H16`; every smooth-projective-rational `P17`--`P19` carrier excluded** | Theorems 2.2--2.4 prove complete nef cones $P_{12},\ldots,P_{16}$ by exact curve-cone decompositions; the `j`-family supplies infinitely many pairwise nonisomorphic carriers and all 28 cusp fibrations (27 primitive, one exact index four); exact facet-cone MMPs conjugate every possible rational `P17`, `P18`, and `P19` marking to the cubic Abel obstruction |
| 2 | pencils at infinity on rational surfaces | **positive and sharply closed in that published lane** | Theorems 3 and 7 of arXiv:2602.15703 give the complete cones; exact root-facet analysis gives $P_{10}\subset\mathbb H^{10}$, $P_{11}\subset\mathbb H^{11}$, and no integral-root case above $\mathbb H^{11}$ inside that construction |
| 3 | Borcherds Reed--Muller lattice | **known abstract positive; divisor-carrier no-go** | `H17/960`, cofinite noncompact; K3, $K3^{[2]}$, and full $K3^{[3]}$ divisor lanes excluded as above |
| 4 | other rational Mori-dream surfaces | open but sharply filtered | arbitrary Picard rank exists, but published redundant-blowup families visibly create square-$\le-3$ curves |
| 5 | Coble--Mukai surfaces | standard mechanism closed above `H9`; raw-Picard alternative open | the natural reflective orthogonal lattice is always rank 10 over $\mathbb C$ |
| 6 | standard Cremona $K^\perp$ chamber | closed above `H9` for finite volume | a Lorentzian proper vertex subdiagram appears at $n=11$ |
| 7 | smooth toric surfaces | closed above `H9` for the integral facet factory | toric Noether identity forces a square-$\le-3$ facet |
| 8 | $K3^{[2]}$ Alexeev chambers | lattice-only for mixed roots | their square-$-4$ walls are not prime-exceptional/movable walls |
| 9 | cone-conjecture/MBM orbit finiteness | insufficient by itself | finitely many wall orbits modulo birational automorphisms is not a finite wall list or a cofinite reflection chamber |

Hwang--Park's [*Cox rings of rational surfaces and redundant
blow-ups*](https://arxiv.org/abs/1303.2274) constructs Mori-dream rational
surfaces of arbitrarily large Picard rank, but its explicit high-rank lane
passes through configurations containing $(-3)$-curves.  This makes those
families valuable negative controls rather than immediate reflection
factories.  It does not prove that every high-rank rational Mori-dream
surface fails.

## 6. Next executable case

The next rational case should be designed, not sampled blindly.

1. Freeze the non-isotrivial Theorems 2.2--2.3a and search their
   surface-realization lineage.  The abstract $P_{12}$--$P_{16}$ lineage is
   closed; only priority for these exact rational nef-carrier families is
   open.
2. Extract a reusable finite facet-cone MMP criterion from the completed
   `P17`--`P19` certificates and apply it to other rational candidates, with
   terminal `P2/F_n` reconstruction and diagram symmetry built into the proof.
3. Search a genuinely nonrational or singular algebraic reservoir, or a
   different finite-volume chamber in `H17`--`H19`; the complete classical
   `P17`--`P19` facet cones are now excluded on every smooth projective
   rational surface, and merely Cremona-relabeling them cannot evade the
   theorem.
4. The isotropic nef rays of $X_{12}$ through $X_{16}$ are now complete over
   the whole nonempty carrier-open family: all 28 are semiample fibration
   classes, with 27 primitive pencils and one exact index-four pencil.  The
   coefficientwise `ZZ[a]` calculation has no additional exceptional
   parameter because its paired chart gates have gcd one.  Combine this
   completed symbolic theorem with the already published $X_{10},X_{11}$
   pencil description to decide the full Mori-dream status.

For modular forms, the corresponding executable rule is also now clear:
use a reflective Borcherds form to generate a candidate lattice and mirror
stock, but require an independently complete simple-root description and a
literal Vinberg volume audit.  Borcherds himself gives reflective modular
forms whose Lorentzian reflection groups are not cofinite; the modular form
is a search oracle, not a volume certificate.
