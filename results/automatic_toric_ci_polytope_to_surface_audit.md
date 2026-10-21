# Automatic toric complete-intersection carriers: theorem, data loss, and BMY no-go

**Status (2026-08-08).**  This note audits a genuinely automatic

\[
  \text{compact Coxeter polytope}
  \longrightarrow
  \text{projective toric pair}
  \longrightarrow
  \text{algebraic surface}
\]

construction.  The incidence theorem below is rigorous and uniform.  The
conclusion is negative for the intended leverage: the construction is a useful
**incidence carrier**, but neither a Coxeter-metric carrier nor a source of a
new dimension/facet bound.  The failure is structural rather than a missing
computer search.

## 1. Executive verdict

Let $P\subset\mathbb H^n$ be a compact Coxeter polytope.  Its face lattice
automatically produces a projective simplicial toric variety, its canonical
smooth toric Deligne--Mumford stack, and general complete-intersection
surfaces.  The ambient toric boundary recovers the whole facet nerve, while
the induced boundary curves on the surface recover exactly the facet-adjacency
graph.

This does **not** give the desired inverse carrier theorem.

1. Rationalizing the Klein-model polytope discards the Lorentzian Gram matrix,
   the arithmetic lattice/field, and the hyperbolic volume.
2. The normal fan does not contain the affine support constants.  Even before
   rational perturbation, those constants enter the Lorentzian wall normals
   and hence the hyperbolic angles.
3. A general surface section avoids all toric strata of codimension at least
   three.  It therefore retains pairwise facet incidence but loses the higher
   face lattice.
4. Toric DM-stack stabilizers and all boundary-root-stack stabilizers are
   diagonalizable, hence abelian.  At a transverse crossing, order-two roots
   give $\mu _2\times\mu _2$, which geometrizes only $m_{ij}=2$.
5. Faithful nonabelian $I_2(m)$ inertia is locally algebraic, but for
   $m>2$ its coarse branch discriminant is not two transverse boundary
   curves.  Thus it cannot enhance the same toric SNC surface pair.
6. For the reduced toric boundary, logarithmic BMY is an explicit positive
   identity depending only on the complete-intersection classes.  It contains
   no fan, face, or Coxeter-label term.
7. Even an artificially engineered label-sensitive orbifold correction is two
   powers lower in the freely scalable complete-intersection polarization and
   is eventually swamped by the positive unlabeled term.

Accordingly this route is classified as

\[
  U_{\mathrm{incidence}},
  \qquad
  \text{not }U_{\mathrm{Gram}},
  \quad
  \text{not }U_{\mathrm{nef}},
  \quad
  \text{not }U_{\mathrm{bound}}.
\]

## 2. The exact automatic theorem

### Theorem 2.1 (automatic toric incidence surface)

Let $P\subset\mathbb H^n$, $n\ge 3$, be a compact hyperbolic Coxeter
polytope with facets $F_1,\ldots,F_f$.  Then there are:

- a rational simple Euclidean $n$-polytope $Q$ with the same labelled face
  lattice as $P$;
- its projective simplicial toric variety $X_Q$;
- the canonical smooth toric DM stack
  $\pi:\mathcal X_Q\rightarrow X_Q$;
- labelled toric boundary divisors
  $\mathcal D_1,\ldots,\mathcal D_f$; and
- for any sufficiently positive ample line bundles
  $A_1,\ldots,A_{n-2}$, a general smooth complete-intersection DM surface

  \[
    \mathcal S=A_1\cap\cdots\cap A_{n-2}\subset\mathcal X_Q
  \]

  with projective coarse surface,

such that the following statements hold.

1. For every set $I\subseteq\{1,\ldots,f\}$,

   \[
     \bigcap_{i\in I}\mathcal D_i\ne\varnothing
     \quad\Longleftrightarrow\quad
     \bigcap_{i\in I}F_i\ne\varnothing .
     \tag{2.1}
   \]

   Thus the labelled ambient toric boundary has exactly the facet nerve of
   $P$.

2. Put $\mathcal C_i=\mathcal D_i\cap\mathcal S$.  Then

   \[
     \mathcal C_i\cap\mathcal C_j\ne\varnothing
     \quad\Longleftrightarrow\quad
     F_i\cap F_j\ne\varnothing .
     \tag{2.2}
   \]

   The intersections on the left are finite and positive when nonempty.

3. A general $\mathcal S$ misses every boundary stratum of codimension at
   least three.  Consequently no three distinct $\mathcal C_i$ meet, and the
   surface boundary records only the one-skeleton of the facet nerve.

4. The ambient toric log canonical class is trivial:

   \[
      K_{\mathcal X_Q}+\mathcal D=0,
      \qquad
      \mathcal D=\sum_i\mathcal D_i.
      \tag{2.3}
   \]

If an honest smooth projective variety is preferred, take a projective regular
toric subdivision and then a general complete-intersection surface.  The
subdivision introduces extra boundary components; the proper map to $X_Q$
and the labels on the original rays are needed to recover (2.1).

### Proof

A compact Coxeter polytope is acute angled.  Its finite vertex links are
spherical Coxeter chambers and hence simplices, so $P$ is simple.  In the
Klein model it is a bounded real simple polytope.  Translate an interior point
to the origin and dualize.  The polar is simplicial.  A sufficiently small
rational perturbation of its vertices preserves its face lattice; polarizing
back gives the required rational simple $Q$.  This is the same rational
perturbation used in the toric proof of the $g$-theorem.

The normal fan $\Sigma_Q$ is rational, complete, projective, and simplicial.
Its rays correspond to the facets of $Q$, and a set of rays lies in a common
cone exactly when the corresponding facets have nonempty intersection.  This
proves (2.1) for the toric orbit stratification.  A simplicial toric variety has
a canonical smooth toric DM stack; equivalently one may use the stacky fan
whose ray maps are the primitive ray generators.

The orbit closure $\mathcal D_i\cap\mathcal D_j$, when nonempty, has dimension
$n-2$.  Intersecting it with $n-2$ general ample divisors gives a positive
zero-dimensional cycle.  If the two boundary divisors are disjoint, their
restrictions remain disjoint.  This proves (2.2).  A boundary stratum of
codimension $k\ge3$ has expected intersection dimension $2-k<0$ with a
general surface section, proving (3).  Finally, the sum of the invariant prime
divisors is the toric anticanonical divisor, proving (2.3).  $\square$

### What is canonical and what is not

The existence is automatic, but the object is not canonical.  The rational
realization $Q$, its lattice, its support numbers, and the ample classes
$A_a$ are choices.  Only the labelled Stanley--Reisner incidence ideal is
forced by the face lattice.  This distinction is fatal for numerical
inequalities intended to transfer back to $P$.

## 3. Exact data-loss ledger

| datum of $P$ | ambient toric pair | complete-intersection surface | reason |
|---|---:|---:|---|
| labelled facet nerve | yes | no, except its graph | toric strata give all cones; a general surface avoids codimension $\ge3$ strata |
| facet adjacency | yes | yes | equation (2.2) |
| finite labels $m_{ij}$ | no | no | normal-fan incidence has no Lorentzian angle |
| right angle versus $m>2$ | no | no | every adjacent pair has a positive curve intersection |
| parallel versus ultraparallel | no | no | all nonadjacent boundary divisors are disjoint |
| Lorentzian Gram and signature | no | no | fan intersections are not the normalized wall pairing |
| arithmetic field, lattice, discriminant | no | no | rational perturbation ignores arithmeticity |
| hyperbolic covolume | no | no | toric intersection volume is ordinary lattice volume |
| ordinary polarized polytope volume | yes, after choosing $Q,H$ | partially | $H^n=n!\operatorname{Vol}_{\rm lattice}(Q)$, but $H$ is noncanonical |
| Hodge-index signature | n/a | yes, abstractly | it applies to every projective surface, not to the Coxeter wall Gram |
| full nef/Mori cone equal to $P$ | no | no | the construction supplies boundary curves, not a cone equality theorem |

There is also a rank mismatch.  Rationally, the toric divisor-class group
generated by the $f$ boundary divisors has rank $f-n$.  The Coxeter wall
Gram has rank $n+1$.  For a simplex the former rank is one while the latter
is $n+1$.  Hodge index therefore supplies a free Lorentzian signature on
$\operatorname{NS}(\mathcal S)$, but not the required wall lattice or its
metric.

In affine Klein coordinates a facet has an equation
$a_i\cdot x\le b_i$.  The normal fan remembers the ray of $a_i$, not the
support constant $b_i$.  The associated Lorentzian hyperplane normal uses
both $a_i$ and $b_i$.  Thus even the normal fan of the unperturbed real
polytope would not determine the Coxeter Gram matrix.

## 4. The log-BMY calculation is a tautology

The strongest clean conclusion of this audit is an exact identity.

### Theorem 4.1 (logarithmic Chern collapse)

Let $X$ be a smooth projective toric $n$-fold, let $D=X\setminus T$ be
its reduced toric boundary, put $r=n-2$, and let

\[
  S=H_1\cap\cdots\cap H_r
\]

be a smooth complete-intersection surface transverse to $D$.  Write
$\Delta=D|_S$ and $h_a=H_a|_S$.  Then

\[
 \begin{aligned}
  c_1\bigl(\Omega_S^1(\log\Delta)\bigr)
    &=\sum_{a=1}^r h_a,\\
  c_2\bigl(\Omega_S^1(\log\Delta)\bigr)
    &=\sum_{a=1}^r h_a^2+\sum_{a<b}h_ah_b.
 \end{aligned}
 \tag{4.1}
\]

In particular

\[
 3\bar c_2-\bar c_1^2
   =2\sum_a h_a^2+\sum_{a<b}h_ah_b>0.
 \tag{4.2}
\]

The same formula holds in rational Chow groups on a smooth toric DM stack.

### Proof

For a smooth toric pair, the logarithmic cotangent bundle is trivial:

\[
  \Omega_X^1(\log D)\simeq M\otimes_{\mathbb Z}\mathcal O_X.
\]

Logarithmic transversality gives the conormal sequence

\[
  0\longrightarrow\bigoplus_{a=1}^r\mathcal O_S(-h_a)
   \longrightarrow\Omega_X^1(\log D)|_S
   \longrightarrow\Omega_S^1(\log\Delta)
   \longrightarrow0.
\]

Therefore

\[
 c\bigl(\Omega_S^1(\log\Delta)\bigr)
   =\prod_{a=1}^r(1-h_a)^{-1}.
\]

Expanding through degree two gives (4.1).  Equation (4.2) follows by direct
subtraction.  Every $h_a$ is ample, so all displayed surface intersection
numbers are positive.  $\square$

Since $K_S+\Delta=\sum h_a$ is ample, logarithmic BMY applies, but (4.2)
shows that its inequality is already a formal positivity identity.  The fan,
the number of boundary components, the face vector, and all labels
$m_{ij}$ have canceled.  For $H_a=H$,

\[
 \frac{\bar c_1^2}{\bar c_2}
   =\frac{2r}{r+1}<2,
\]

so these pairs remain uniformly far from the BMY boundary (3), even as the
dimension grows.

## 5. Why a toric/root stack cannot remember $m_{ij}>2$

Smooth toric DM stacks are built from stacky fans; their stabilizers are
finite diagonalizable groups.  The bottom-up classification by canonical
stacks, roots of boundary divisors, and gerbes makes the same point
geometrically.  Consequently all toric inertia is abelian.

If one takes a square root along every $\mathcal C_i$, the generic inertia of
a wall curve is $\mu _2$, while a transverse crossing has

\[
  I_p=\mu _2\times\mu _2.
\]

The product of the two distinguished involutions has order two.  This is the
Coxeter group $I_2(2)$, not $I_2(m)$ for $m>2$.

The obstruction persists outside the toric category if the same coarse SNC
boundary is retained.

### Proposition 5.1 (local SNC/dihedral incompatibility)

Let $\mathcal Y$ be a smooth effective tame complex DM surface with smooth
coarse surface $Y$.  Suppose that near $p\in Y$ its branch divisor consists
of two smooth transverse branches, each with generic reflection inertia of
order two.  If the two distinguished reflections generate the full inertia
$I_2(m)$ of order $2m$, then $m=2$.

### Proof

The local linear reflection representation can be written

\[
 r(x,y)=(\zeta x,\zeta^{-1}y),
 \qquad
 s(x,y)=(y,x),
 \qquad \zeta^m=1.
\]

Its invariant ring is

\[
 \mathbb C[x,y]^{I_2(m)}
   =\mathbb C[u,v],
 \qquad
 u=xy,
 \quad
 v=x^m+y^m.
\]

The image of the reflecting lines is the discriminant

\[
  v^2-4u^m=0.
  \tag{5.1}
\]

For odd $m>1$, (5.1) is irreducible and singular.  For $m=2k>2$, it has
two smooth branches $v=\pm2u^k$ with intersection multiplicity $k>1$, so
they are tangent.  Only $m=2$ gives two transverse branches.  $\square$

The Chevalley--Shephard--Todd theorem explains why the coarse quotient is
nevertheless smooth: $I_2(m)$ is generated by complex reflections.  The
price is the singular discriminant (5.1).

### Global consequence

A faithful chart $[\mathbb A^2/I_2(m)]$ exists algebraically for every $m$,
but it cannot be glued in while leaving the toric SNC crossing unchanged.
One must replace that crossing by (5.1).  For odd $m$, the two wall branches
are even identified in the coarse discriminant; for even $m>2$, they become
tangent.  Globalizing all such replacements is a separate branched-cover and
orbifold-monodromy problem, not a consequence of the normal fan.

Putting a noneffective residual gerbe $BI_2(m)$ at a point would merely store
the label as decoration.  It would not act by the two wall reflections on the
normal directions and therefore is not geometric recovery of the Coxeter
angle.

Finally, a surface section avoids rank-three and higher toric strata.  Even a
successful pairwise dihedral decoration would not retain the compatibility of
the finite spherical Coxeter subgroups occurring at vertices of $P$.

## 6. The hypothetical orbifold inequality and the scaling obstruction

It is useful to record exactly what would happen if the local problem were
ignored.  Suppose hypothetically that a compact orbifold surface had generic
order-two inertia along curves $C_i$ and full inertia $I_2(m_p)$ at their
crossings.  Weighted topological stratification would give

\[
 e_{\mathrm{orb}}
  =e(S\setminus C)
   +\frac12\sum_i e(C_i^\circ)
   +\sum_p\frac1{2m_p}.
 \tag{6.1}
\]

The transverse double-root stack has $1/4$ in the last summand.  Thus the
formal label-sensitive correction is

\[
  \sum_p\left(\frac1{2m_p}-\frac14\right),
  \tag{6.2}
\]

which is negative for $m_p>2$.  This is the only plausible BMY entry point
for the finite labels in this construction.

It still cannot yield an intrinsic obstruction.  Take
$H_1=\cdots=H_r=kH$, where $r=n-2$.  The positive unlabeled leading term in

\[
 3c_2^{\mathrm{orb}}-(K_{\mathrm{orb}})^2
\]

is

\[
  \frac{r(r+3)}2\,k^{r+2}H^n.
  \tag{6.3}
\]

Boundary terms have order at most $k^{r+1}$, while the number of induced
pairwise crossing points, and hence (6.2), has order only

\[
  D_iD_j(kH)^r=O(k^r).
\]

For every fixed labelled polytope, the positive term (6.3) therefore dominates
as $k\to\infty$.  Since the normal-fan construction supplies no canonical
primitive polarization and permits this scaling freely, orbifold BMY cannot
reject a polytope or bound its dimension/facet number.

## 7. Relation to the $g$-theorem

This toric association applies to every simple polytope, whether Coxeter,
arithmetic, or neither.  Toric hard Lefschetz for the dual simplicial polytope
is precisely the classical algebraic-geometric source of the necessity half
of the $g$-theorem.  Surface Hodge index in the present construction gives
ordinary intersection-form or Alexandrov--Fenchel/Khovanskii--Teissier
positivity inherited from the same toric intersection algebra.

For the particular proposed inequality there is an even stronger answer:
log-BMY reduces to the identity (4.2), so it does not produce an additional
face-number condition at all.  Any Coxeter-specific inequality would have to
reintroduce data that this normal-fan functor discarded.

## 8. Go/no-go decision and the precise missing lemma

The automatic theorem is worth retaining as a cheap algebraic incidence
compiler and as a source of explicit projective test surfaces.  It should not
be counted as progress toward a universal divisor-lattice/nef-cone transfer or
a new compact-Coxeter bound.

There is no useful next enumeration on this lane.  The exact go/no-go lemma
needed to revive it is the following.

> **Canonical label-sensitive surface lemma.**  Starting from the arithmetic
> Lorentzian Gram datum of $P$, not merely its face lattice, construct
> functorially a projective effective DM surface whose distinguished divisors
> carry faithful reflection inertia $I_2(m_{ij})$, whose higher spherical
> parabolics remain visible, and whose polarization/log-canonical class has a
> canonical normalization forbidding $H\mapsto kH$.  Prove an orbifold-Chern
> identity in which the label term occurs in the same leading degree as the
> positive term.

Proposition 5.1 proves that such a lemma cannot use a general toric
complete-intersection surface with its SNC boundary.  It must instead use a
different carrier, such as an arithmetic reflection/Heegner stack where the
dihedral inertia and the polarization are intrinsic.  Until that lemma is
available, this toric lane is closed for theorem-boundary leverage.

## 9. Primary sources

- E. B. Vinberg, [*Discrete groups generated by reflections in Lobachevskii
  spaces*](https://www.mathnet.ru/eng/sm4099), Math. USSR-Sb. **1** (1967),
  429--444.
- V. I. Danilov, [*The geometry of toric
  varieties*](https://doi.org/10.1070/RM1978v033n02ABEH002305), Russian Math.
  Surveys **33** (1978), 97--154.
- R. P. Stanley, [*The number of faces of a simplicial convex
  polytope*](https://klein.mit.edu/~rstan/pubs/pubfiles/46.pdf), Adv. Math.
  **35** (1980), 236--238.
- L. A. Borisov, L. Chen, and G. G. Smith, [*The orbifold Chow ring of toric
  Deligne--Mumford stacks*](https://arxiv.org/abs/math/0309229), J. Amer. Math.
  Soc. **18** (2005), 193--215.
- B. Fantechi, E. Mann, and F. Nironi, [*Smooth toric DM
  stacks*](https://arxiv.org/abs/0708.1254), J. Reine Angew. Math. **648**
  (2010), 201--244.
- G. C. Shephard and J. A. Todd, [*Finite unitary reflection
  groups*](https://doi.org/10.4153/CJM-1954-028-3), Canad. J. Math. **6**
  (1954), 274--304; C. Chevalley, [*Invariants of finite groups generated by
  reflections*](https://doi.org/10.2307/2372597), Amer. J. Math. **77**
  (1955), 778--782.
- A. Langer, [*Logarithmic orbifold Euler numbers of surfaces with
  applications*](https://arxiv.org/abs/math/0012180), Proc. London Math. Soc.
  **86** (2003), 358--396.
