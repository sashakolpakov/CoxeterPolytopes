# The K3-cover bridge from a Coxeter chamber to a Coble surface

**Status (2026-08-09): quotient theorem and complete lift taxonomy proved;
template-free finite-root promotion theorem proved; the unrestricted Coble
L4 claim remains open.**  The invariant 2-elementary lattice alone produces
a strong Coble quotient, but it does **not** return the input polytope.  The
missing datum is an equivariant unfolding of every divisor wall into actual
K3 roots, followed by a period/ramification and projected-hidden-root test.
In characteristic three this must be a crystalline, not a complex-period,
test.

This distinction is decisive for the two positive endpoints.  The
characteristic-three Type-VI and Type-MII covers are supersingular K3
surfaces.  Their invariant lattices have Nikulin types `(15,7,1)` and
`(18,4,0)`, but the same invariant types also occur over the complex
numbers.  The generic complex quotients are strong Coble surfaces with five
and eight boundary components; their nef cones are **not** the finite
`H14/60` and `H17/48` chambers.  Thus

```text
correct invariant lattice + correct fixed-locus topology
    does not imply
effectivity of the proposed walls + exact same-P nef return.
```

The results below give the non-tautological compiler that is needed to close
that implication and isolate the gates which are not yet certified.

## 1. A fixed rational locus automatically gives a strong Coble quotient

### Theorem 1.1 (fixed-locus-to-Coble theorem)

Let `k` be algebraically closed of characteristic different from two.  Let
`X` be a K3 surface and let `tau` be a tame anti-symplectic involution whose
fixed locus is

\[
       X^\tau=R_1\sqcup\cdots\sqcup R_b,
       \qquad R_i\simeq\mathbf P^1,\qquad b>0.
\]

Then the quotient `pi:X -> S=X/<tau>` is a smooth rational surface.  If
`B_i=pi(R_i)`, then

\[
 B_i^2=-4,\qquad B_iB_j=0\ (i\ne j),\qquad
 |-K_S|=\varnothing,\qquad
 |-2K_S|=\{B_1+\cdots+B_b\}.                         \tag{1.1}
\]

In particular `S` is a strong Coble surface with exactly `b` boundary
components.

### Proof

The local action along a fixed curve is `(u,v) -> (u,-v)`.  Hence `S` is
smooth and `pi` is the double cover branched along
`B=sum B_i`, with

\[
       \pi^*B_i=2R_i.
\]

Adjunction on the K3 surface gives `R_i^2=-2`; the projection formula then
gives

\[
 2B_i^2=(\pi^*B_i)^2=(2R_i)^2=-8,
\]

so `B_i^2=-4`.  Let `L` be the double-cover line bundle, so that
`2L~B` and

\[
        \pi_*O_X=O_S\oplus L^{-1},\qquad
        K_X\sim\pi^*(K_S+L).
\]

Since `K_X~0`, the line bundle `pi^*(K_S+L)` is trivial.  Taking its norm
gives `2(K_S+L)~0`; this conclusion does not require pullback on `Pic(S)`
to be injective.  Using `2L~B` therefore gives, before rationality or any
torsion argument,

\[
                         B\sim-2K_S.                    \tag{1.2}
\]

For a tame quotient,
`H^i(O_S)=H^i(O_X)^tau`.  Thus `q(S)=0`; and `p_g(S)=0` because tau is
anti-symplectic, so its action on `H^0(omega_X)` (equivalently, by Serre
duality, on `H^2(O_X)`) is `-1`.  Since an ample divisor has negative
intersection with `K_S=-B/2`, no positive multiple of `K_S` is effective;
`S` has Kodaira dimension `-infinity`.  The classification of ruled
surfaces, together with `q(S)=0`, makes `S` rational.

There are two equivalent proofs of the linear-system assertions.

* Downstairs, if `D` is effective and `D~-K_S`, then
  `D.B_i=-2`.  Hence every `B_i` is a component of `D`, and `D-B` is an
  effective canonical divisor.  This contradicts `p_g(S)=0`.  If
  `D~-2K_S`, then `D.B_i=-4`; again every `B_i` is a component, and
  `D-B` is effective and linearly trivial.  Therefore `D=B`.

* Upstairs, put `R=sum R_i`.  For every `m>=1`, the negative-definite
  support shows `|mR|={mR}`: intersect a putative member with each `R_i`,
  remove one copy of every component, and induct on `m`.  Pulling back an
  anticanonical or bi-anticanonical member gives the same conclusion.

This proves (1.1).  Notice that neither finite automorphisms nor Kondō's
nine-class list is used.  This is the robust direction

```text
anti-symplectic K3 pair with rational fixed locus -> strong Coble surface.
```

## 2. What the raw boundary frame forces in the invariant lattice

Let a raw marked wall system reconstruct an odd unimodular lattice

\[
 N\simeq I_{1,9+b},
\]

pairwise orthogonal vectors `B_i` of square `-4`, and

\[
             K=-\frac12\sum_{i=1}^bB_i\in N.            \tag{2.1}
\]

Put `h_i=B_i/2` and

\[
 \widetilde N=\langle N,h_1,\ldots,h_b\rangle,
 \qquad M=\widetilde N(2).                              \tag{2.2}
\]

Here `(2)` means that the form is multiplied by two.  The module `M` is
integral and even.  Indeed `2(x,h_i)=x.B_i` is integral before scaling,
and all scaled squares are even.  Also `2M^* subset M`: if
`x in M^*`, pairing with the unimodular sublattice `N` shows `2x in N`.
Thus `M` is 2-elementary.

There is one further saturation condition which must be tested, not
assumed.  Require that the only relation among the boundary halves modulo
`N` is

\[
                 h_1+\cdots+h_b=-K\in N.                \tag{2.3}
\]

Equivalently `[widetilde N:N]=2^(b-1)`.  Then, with
`r=rank(M)=10+b`, unimodularity of `N` gives

\[
 |\det M|
 =\frac{2^{10+b}}{2^{2b-2}}
 =2^{12-b}.                                             \tag{2.4}
\]

Consequently the 2-elementary invariants have

\[
                         (r,a)=(10+b,12-b),              \tag{2.5}
\]

while `delta` is computed from the discriminant quadratic form.  Condition
(2.3), the discriminant form, and primitive K3-lattice embeddability are all
finite Smith/discriminant-form calculations on the polytope marking.

Over the complex numbers, suppose `M` is an admissible invariant lattice of
a non-symplectic K3 involution.  Nikulin's fixed-locus formula is

\[
 g=\frac{22-r-a}{2},\qquad k=\frac{r-a}{2},              \tag{2.6}
\]

where the fixed locus is one genus-`g` curve and `k` additional rational
curves, outside the two rank-ten exceptional cases.  Substitution of (2.5)
gives

\[
                         g=0,\qquad k=b-1.               \tag{2.7}
\]

Thus the invariant lattice forces exactly `b` rational fixed components.
Theorem 1.1 then gives a strong Coble quotient.  Before using any finite-nef
classification, this also gives the coarse K3 bound

\[
                  b\le10,\qquad \dim P=9+b\le19.        \tag{2.8}
\]

The parity check has geometric content: `delta=0` if and only if the class
of the total fixed locus is divisible by two in the K3 lattice.  It does not,
by itself, identify which effective roots are the fixed components.

## 3. Why `M` alone cannot return the polytope

For a generic complex 2-elementary K3 pair of type `M`, one may take
`NS(X)=M`; the involution is `+1` on all of `NS(X)`.  This realizes the
fixed-locus topology (2.7), but it loses precisely the split curves needed
by the Coble wall system.

If `C subset S` is a `(-2)`-curve, then `C.B=0` and

\[
 \pi^*C=C^++C^-,\qquad
 (C^\pm)^2=-2,\qquad C^+C^-=0,\qquad \tau C^+=C^- .     \tag{3.1}
\]

Thus the invariant wall class has square `-4` in `M`, but its two K3 roots
require anti-invariant algebraic classes.  They cannot occur when
`NS(X)=M`.

There is a similar issue for a `(-1)`-curve meeting one boundary component
twice.  If the two intersections are distinct, its pullback is one
tau-stable `(-2)`-curve.  If the intersection divisor is twice one point,
the pullback splits into two exchanged `(-2)`-curves meeting once.  The
intersection numbers in `N` do not distinguish these alternatives.

This explains the endpoint phenomenon without hand-waving:

* Type VI in characteristic three has ten actual `(-2)` walls, hence at
  least those ten exchanged-root pairs upstairs.
* Type MII has twenty-four actual `(-2)` walls, hence twenty-four such
  pairs.
* Kondō proves that the K3 covers of Types V, VI, VII, MI, and MII are
  supersingular of Picard rank `22`.  In particular the H14 and H17
  endpoints are not generic complex period points.

So Nikulin existence produces many strong Coble surfaces, but not the
claimed finite-nef chamber.  This is an actual obstruction to the shortcut

```text
raw P -> M -> choose a generic M-polarized K3 -> exact P.
```

## 4. The exhaustive K3 lift alphabet

The correct polytope-first object is an **orbit-root unfolding**.  Every
physical negative wall has one of the following exhaustive lift types.

| downstairs wall | boundary profile | roots upstairs | internal product |
|---|---|---|---:|
| `B_i^2=-4` | boundary | one pointwise fixed root `r_i`, with `pi^*B_i=2r_i` | -- |
| `E^2=-1` | `(1,1)` | one tau-stable root `e=pi^*E` | -- |
| `E^2=-1` | `(2)` at two distinct points | one tau-stable root `e=pi^*E` | -- |
| `E^2=-1` | tangency divisor `2p` | two exchanged roots `e,e'`, `pi^*E=e+e'` | `1` |
| `C^2=-2` | zero | two exchanged roots `c,c'`, `pi^*C=c+c'` | `0` |

The first, second, and fifth rows are forced.  For the fourth/fifth
distinction inside a numerical `(2)` profile, the compiler must enumerate
both possibilities.  Riemann--Hurwitz, or the square/non-square branch
section on `P1`, proves that the table is exhaustive in characteristic not
two.

All cross-pairings have a finite arithmetic enumeration.  Let `D,E` be
downstairs non-boundary curves with `D.E=d`.

1. Two tau-stable roots have product `2d`.
2. A stable root and an exchanged pair have products `(d,d)`.
3. For exchanged pairs `{u,u'}` and `{v,v'}`, tau symmetry gives

   \[
     u.v=u'.v'=x,\qquad u.v'=u'.v=y,\qquad x+y=d,        \tag{4.1}
   \]

   so `x` ranges over the finite set `0,...,d`.
4. A boundary root `r_i` has product `B_i.D` with a stable lift.  With
   an exchanged lift it has equal products `B_i.D/2` with the two roots;
   the split rows have profiles for which this is integral.

All these products are nonnegative because they are intersections of
distinct irreducible curves.  Therefore a raw finite wall Gram gives only
finitely many orbit-root Gram matrices:

\[
  2^{\#\{(2)\text{-profile walls}\}}
  \prod_{\{D,E\}\text{ both split}}(D.E+1),             \tag{4.2}
\]

before rank, signature, and symmetry pruning.  This is the promised finite
compiler denominator; it does not import one of Kondō's answers.

For every candidate Gram, quotient by its saturated integral radical, take
all finite-index even tau-stable overlattices, and impose:

```text
L is nondegenerate of signature (1,rho-1);
the orbit-sum embedding identifies L^tau with the reconstructed M;
the physical wall classes generate the prescribed descent lattice N;
L has a primitive equivariant embedding in the K3 lattice.
```

There are finitely many overlattices because they are isotropic subgroups of
a finite discriminant group.  For each fixed nondegenerate `L`, primitive
embedding in the K3 lattice is an effective integral-lattice problem: one
must enumerate complements with discriminant form `-q_L`, the finite gluing
anti-isometries, and the compatibility of the two involutions on those
discriminant groups.  Finiteness follows from finite class number at fixed
rank/signature/discriminant and from finiteness of the discriminant groups.
This is a proof that `U3` has a finite mathematical search space, **not** a
claim that the equivariant embedding algorithm has already been implemented
or certified for every compiler output.  Merely writing “by discriminant
forms” would leave an unjustified gap.  Conversely, the actual lifts on
every genuine Coble carrier occur in this enumeration, so the orbit-root
Gram stage itself is essentially exhaustive.

What it does **not** yet prove is that one of these finite lattices has the
right period, ramification roots, and invariant nef slice.

## 5. A genuine template-free complex promotion theorem

The following restricted theorem has the desired logical direction.  It is
stated separately from the unfinished claim that all Coble candidates pass
its hypotheses.

### Theorem 5.1 (finite-root cover promotion)

Let a raw marked Coxeter datum `P` pass the boundary-frame construction of
Section 2 and the exhaustive lift construction of Section 4.  Suppose one
candidate `(L,tau)` has all of the following properties.

1. `L` has rank at most `20` and admits a primitive tau-equivariant
   embedding in the complex K3 lattice.  The fixed lattice of the extended
   involution is exactly `M`, and the orthogonal period domain is nonempty.
2. There is a tau-invariant chamber `A` of the **complete** `(-2)`-root
   arrangement of `L`, with a finite certified set `Delta` of simple roots,
   and `A` contains a tau-fixed timelike interior vector lifting an interior
   point of the raw chamber.
3. Among subsets of `Delta^tau`, the proposed `r_1,...,r_b` are, up to an
   automorphism of the marked datum, the unique ramification frame.  Here a
   ramification frame means `b` mutually orthogonal roots such that, for
   `F=sum r_i`:

   * every other tau-fixed simple root `c` has `c.F=2`;
   * for every two-element simple-root orbit relevant to a negative
     quotient wall, `c.F=c.tau(c)` (hence `0` or `1`);
   * `F/2` is integral in the K3 lattice exactly when `delta=0`.

4. Intersecting `A` with `M_R`, quotienting the simple roots by tau, and
   minimizing redundant inequalities gives exactly the oriented raw chamber
   `P`.  Its orbit walls have precisely the squares and integral descent
   lattice `N` prescribed in Section 4.

Then there is a complex K3 surface `X` with anti-symplectic involution tau
whose quotient is a strong Coble surface `S`, and

\[
                         \operatorname{Nef}(S)=P         \tag{5.1}
\]

under the recovered marking.

### Proof

Choose a generic period in the anti-invariant orthogonal complement of `L`.
Surjectivity of the K3 period map gives a marked K3 surface with
`NS(X)=L`.  The chosen chamber is made the Kähler chamber by changing the
marking by the Weyl group.  Since tau fixes an interior vector of the
chamber, it preserves the chamber.  Global Torelli therefore realizes the
lattice involution by an anti-symplectic automorphism.

Nikulin's formula and (2.5) show that the fixed locus consists of `b`
rational curves.  Their classes are tau-fixed simple roots, are mutually
orthogonal, satisfy the two intersection conditions in item 3, and have the
stated divisibility parity.  Uniqueness of the ramification frame therefore
identifies them with `r_1,...,r_b`.

Theorem 1.1 makes the quotient a strong Coble surface.  A tau-fixed
non-ramification simple root descends to a `(-1)`-curve.  A two-root orbit
with internal product `m=0` or `1` descends to a curve of square

\[
 \frac{(c+\tau c)^2}{2}=-2+m.                           \tag{5.2}
\]

Finite pullback identifies nefness downstairs with nefness upstairs on the
invariant subspace.  Item 4 is exactly the finite polyhedral equality of
that invariant slice with `P`, proving (5.1).  Integral generation gives the
claimed Picard marking rather than merely a rational isometry.

This theorem is not a lookup in the nine Coble templates: its hypotheses are
K3-lattice, root-arrangement, involution, and polyhedral calculations.  But
the present project has **not** yet certified a non-singleton raw-P family
passing all four items.  In particular, calling it the completed Coble L4
theorem now would be premature.

The existing complex Type-I/Type-II `H10/22`, `H10/30` and Type-I
`H11/26` surface constructions prove that plausible positive inputs exist,
and the hostile `H10/30` replay reconstructs its boundary/CM/cusp package
from the polytope.  They do not yet provide raw-P certificates for the full
K3 lift lattice, every projected K3 root, and uniqueness of the
ramification frame.  Consequently even the char-zero **non-singleton**
version of Theorem 5.1 is a concrete next computation, not a theorem already
earned by those forward constructions.

### Why the ramification test is necessary

Nikulin's invariants count the fixed components but do not name their root
classes.  A tau-stable rational curve need not be pointwise fixed: if it is
not, the involution on `P1` has two fixed points, and the curve meets the
fixed divisor transversely in two points.  This proves the second numerical
condition above but also shows why merely selecting `b` orthogonal invariant
roots is insufficient.  The actual fixed frame must be forced by a complete
root-chamber calculation or by an independent geometric model.

## 6. The projected hidden-wall gate

If the complete K3 simple-root set is not finite, the right condition is
still purely arithmetic but is no longer automatically a finite
calculation.  Fix an invariant timelike vector `h` in the interior of `P`.
For every K3 root `alpha` choose the sign with `h.alpha>0` and put

\[
             p_+(\alpha)=\frac{\alpha+\tau\alpha}{2}.   \tag{6.1}
\]

There must be no anti-invariant root (`p_+(alpha)=0`), since an invariant
ample class would be orthogonal to it.  Exact same-P return is equivalent to

\[
 \{x\in M^+:x.w\ge0\text{ for every raw wall }w\}
 =
 \{x\in M^+:x.p_+(\alpha)\ge0\text{ for every positive root }\alpha\}.
                                                               \tag{6.2}
\]

Equivalently every projected positive-root inequality is in the dual cone
of the raw walls and every raw facet is attained by a root orbit.  This is
the precise hidden-wall statement.  Writing “full root completeness” without
proving (6.2) merely renames the problem.

For a finite full K3 Weyl chamber, Theorem 5.1 checks (6.2) by exact
polyhedral comparison.  For supersingular covers or an infinite K3 chamber,
one needs a separate finite reduction, usually by a Conway/Vinberg chamber,
a group fundamental domain, or a direct classification of curve orbits.
No such template-free reduction is presently frozen for all nine Coble
rows.

## 7. The characteristic-three gate is crystalline

The H14 and H17 endpoint covers have Picard rank `22`, so they cannot be
constructed by the complex period argument of Theorem 5.1.  A supersingular
compiler must instead output

```text
the supersingular K3 lattice Lambda_(p,sigma);
a strictly characteristic Ogus subspace in its discriminant space;
an involution preserving that characteristic subspace;
the required anti-symplectic action on the period line;
a preserved ample chamber;
the complete projected-root equality (6.2).
```

Crystalline Torelli then produces the K3 automorphism.  An abstract isometry
of the rank-22 lattice is not enough: failure to preserve the characteristic
subspace is a genuine period obstruction.  Even after Torelli, the selected
ramification frame and the tame fixed locus must be checked.

For Artin invariant one the period space is rigid, which makes the Type-MII
lane especially promising, but this does not remove the ample-cone,
anti-symplectic, ramification, or hidden-wall checks.  Type VI and MII are
therefore positive `SURFACE -> POLYTOPE` endpoints and exact template-based
inverse calibrations, but not yet products of a general crystalline
polytope-to-surface compiler.

## 8. Exactly when Kondō's classification pushes back to polytopes

The following implication is unconditional once any cover construction has
passed the ramification, descent, and exact-nef gates.

### Theorem 8.1 (classification pushback)

Let `P` be a finite-sided Coxeter polytope and suppose a template-free
polytope-first construction produces a strong Coble surface `S` in
characteristic different from two with

\[
                       P=\mathbf P(\operatorname{Nef}S).\tag{8.1}
\]

Then `Aut(S)` is finite.  Consequently Kondō's classification applies and

\[
 b\in\{1,2,5,8\},\qquad
 \dim P=9+b\in\{10,11,14,17\}.                          \tag{8.2}
\]

More precisely the surface belongs to exactly one of these nine rows:

| type | characteristic | `b` | `dim P` |
|---|---:|---:|---:|
| I | any `!=2` | 1 | 10 |
| I | any `!=2` | 2 | 11 |
| II | any `!=2` | 1 | 10 |
| V | 3 | 2 | 11 |
| VI | 5 | 1 | 10 |
| VI | 3 | 5 | 14 |
| VII | 5 | 1 | 10 |
| MI | 3 | 2 | 11 |
| MII | 3 | 8 | 17 |

### Proof

Every automorphism preserves the nef cone and hence permutes its finite set
of facets.  Use the primitive integral outward normal on each oriented
facet.  An integral isometry taking one oriented facet to another takes its
primitive normal to the corresponding primitive normal: no nontrivial real
scaling remains.  The normals span `NS(S)_R` because (8.1) is a
full-dimensional nondegenerate hyperbolic polytope.  Consequently the image
of `Aut(S)` in `O(NS(S))` acts faithfully on this finite normal set and is
finite.  The kernel of `Aut(S)->O(NS(S))` acts trivially on the
Coble--Mukai lattice; it is contained in the finite kernel of Kondō's
natural map `Aut(S)->O(CM(S))`.  Thus `Aut(S)` is finite, and his main
classification theorem gives the table.  Finally

\[
 K_S^2=-b,\qquad \rho(S)=10+b,
\]

so the projectivized positive/nef cone has dimension `9+b`.

This is real algebraic-geometric leverage.  On any independently defined
polytope class for which the cover compiler is essentially surjective, it
immediately proves:

```text
no Coble-nef members in H12, H13, H15, H16, or above H17;
the only possible dimensions are H10, H11, H14, H17;
the H14 and H17 members can exist only in characteristic 3.
```

It also supplies positive boundary strata: the certified Type-VI
`H14/60` and Type-MII `H17/48` full-nef chambers.  What is missing is not
Kondō's pushback; it is the template-free proof that a non-singleton raw
polytope denominator automatically passes the cover compiler.

All these polytopes are noncompact.  A smooth rational Coble surface can be
contracted to a minimal rational surface; pulling back a ruling (using the
penultimate `F_1` model if the minimal model is `P2`) gives a nonzero nef
isotropic class.  Hence its projectivized nef cone meets the ideal boundary.

## 9. Completion ledger

The genuine `POLYTOPE -> SURFACE` route now has the following exact gates.

| gate | content | status |
|---|---|---|
| `U0` | enumerate every boundary frame and odd-unimodular de-scaling | proved in principle; frozen for H10 and Type-VI H14 |
| `U1` | compute `M=2<N,B_i/2>`, `(r,a,delta)`, and K3 admissibility | theorem above; finite |
| `U2` | enumerate every stable/split K3 wall lift and every cross-pairing | theorem above; finite |
| `U3` | enumerate even tau-stable closures, K3 complements/glue, and equivariant embeddings | finite search space proved above; full algorithm/certificates not yet completed |
| `U4C` | complex period plus complete chamber and unique ramification frame | Theorem 5.1; no non-singleton raw family yet certified |
| `U4p` | supersingular characteristic subspace plus crystalline Torelli | open for the raw H14/H17 inputs |
| `U5` | prove projected all-root equality (6.2), excluding hidden cuts | open uniformly; surface/template proofs exist row by row |
| `U6` | quotient, strong-Coble lemma, Picard descent, exact same-P return | proved once `U4/U5` pass |
| `AG` | finite Aut, Kondō pushback, `dim<=17`, nine rows | proved once `U6` passes |

Therefore the defensible rating is

```text
strong-Coble quotient from a rational K3 fixed locus:       proved;
template-free finite lift compiler:                         proved;
restricted finite-root complex promotion theorem:           proved;
Kondo dimension/type pushback after exact return:            proved;
non-singleton char-0 raw denominator passing all gates:       not yet proved;
crystalline H14/H17 essential surjectivity:                  not yet proved;
unrestricted nine-class L4 polytope theorem:                 not yet proved.
```

## 10. Sources

1. V. V. Nikulin, classification of K3 involutions and 2-elementary
   lattices; the fixed-locus formulas are recalled in the modern literature
   on 2-elementary K3 surfaces.
2. Valery Alexeev,
   [*Reflective hyperbolic 2-elementary lattices, K3 surfaces and
   hyperkähler manifolds*](https://arxiv.org/abs/2209.09110), especially
   the distinction between `W_2` (K3 nef walls) and the larger mixed
   `W_r` generated by `(-2)` and divisibility-two `(-4)` roots.
3. Shigeyuki Kondō,
   [*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
   especially the main nine-class theorem, the finite-kernel statement for
   `Aut(S)->O(CM(S))`, the K3-cover/R-invariant discussion, and the Type-VI
   and Type-MII constructions.
4. A. Ogus, crystalline Torelli for supersingular K3 surfaces; this is the
   required positive-characteristic replacement for the complex period
   argument.
