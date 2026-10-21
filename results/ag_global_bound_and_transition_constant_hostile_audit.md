# Hostile audit: can algebraic geometry supply a global Coxeter dimension bound?

**Date:** 2026-08-10  
**Verdict:** **no unconditional global bound follows from any presently located
surface theorem.**  There is a precise conditional route through Alexeev's
boundedness theorem, but its uniform discrepancy hypothesis is not implied by
reflection maximality or arithmeticity.  Assuming it is equivalent to assuming
the missing uniform Coxeter-label/complexity bound.  Nikulin's surface papers do
not repair this: their only unconditional boundedness theorem fixes the Picard
rank, and their rank bounds import hyperbolic reflection-lattice results.

The successful part of the surface analogy is instead a **method transfer**:
exceptional-curve diagrams led to narrow-face, four-vertex V-arithmetic Gram
graphs.  Those small graphs can be treated as exact arithmetic dessins.  This
does give effective field-degree funnels, but it is direct hyperbolic
geometry/number theory, not an object-level `polytope -> surface -> theorem on
the same polytope` functor.  It also does not give a dimension cutoff.

This note treats Nikulin's and Vinberg's numerical theorems as regression
targets, not project dependencies.  In particular, the final journal version of
Nikulin's transition-constant paper repairs the polynomial defect in arXiv v4;
the obsolete preprint defect is not a current theorem gap.  A complete
independent exact replay of all remaining numerical branches is still required
before this project imports `N(14)<=25`.

## 1. The answer in one implication diagram

The tempting surface proof would be

```text
arithmetic maximal Coxeter chamber P in H^n
       |
       | universal exact carrier, not presently known
       v
surface pair (X,B),  n+1 <= rk Pic(X),  -(K_X+B) nef
       |
       | Alexeev, but only if b_j <= 1-epsilon_0 for one absolute epsilon_0>0
       v
rk Pic(X) <= A_1(epsilon_0)
       |
       v
n <= A_1(epsilon_0)-1.
```

The last two arrows are valid.  The first arrow and the absolute
`epsilon_0` are not known.  If the natural wall coefficient is
`b_j=1-1/m_j`, then

```text
epsilon_0 <= 1/max_j(m_j).
```

Thus a uniform discrepancy is already a uniform bound on the Coxeter labels.
No theorem located here derives that bound from maximality.  Alexeev gives an
explicit countermechanism when `epsilon -> 0`: repeated blowups produce smooth
log Fano pairs with strictly increasing Picard rank.

For an ABSW-style finiteness result one needs still more: a bounded marked
family and a uniform bound on the number/degree of the divisors encoding the
facets.  Boundedness of unmarked surfaces over a finite-type base does **not**
bound degrees of their closed points over `Q`; a finite-type variety has closed
points of unbounded number-field degree.  A bounded Belyi/dessin degree would
make the marked combinatorics finite, but no surface theorem supplies that
degree from reflection maximality.

## 2. What the primary algebraic-geometry theorems actually say

### 2.1 Nikulin's finite-polyhedral Mori-cone theorem fixes the rank

In V. V. Nikulin, *On algebraic varieties with finite polyhedral Mori cone*,
arXiv:math/9806047, Theorem 1.1, the input is a smooth projective surface with
finite-polyhedral Mori cone and **fixed**

```text
rho = rk NS(X),
delta_E = max{-E^2 : E is an exceptional curve},
p_E = max{p_a(E) : E is an exceptional curve}.
```

The conclusion is the existence of an ample effective divisor `h` with
`h^2 <= N(rho,delta_E)` and, over `C`, a very ample divisor `h'` with
`(h')^2 <= N'(rho,delta_E,p_E)`.  This is bounded polarization after the
dimension variable `rho-1` has already been fixed.  It cannot give a global
rank bound.

The same paper makes this failure literal.  Example 1.2 constructs surfaces
with finite-polyhedral Mori cone, fixed `delta_E>=3` and fixed `p_E`, but

```text
rho(X_k)=3+k.
```

The rank bounds in Examples 1.3 and 1.4 are not independent terminal AG
inputs:

- `delta_E=1` gives `rho<=9` by the classification/finiteness of hyperbolic
  lattices reflective in square-minus-one vectors;
- `delta_E=2` gives `rho<=22`, and the pure square-minus-two case gives
  `rho<=19`, by Vinberg/Esselmann/Nikulin reflective-lattice results.

Using those conclusions to prove a general Coxeter dimension bound would be
circular.

There is a genuine geometric statement in the proof: the projectivized nef
cone of a finite-polyhedral surface has finite hyperbolic volume.  It recognizes
finite volume once an exact complete-nef carrier is already present; it does
not show that an arbitrary Coxeter chamber has such a carrier or bound its
rank.

### 2.2 Alexeev gives exactly the right-shaped conditional theorem

In V. Alexeev, *Boundedness and `K^2` for log surfaces*,
arXiv:alg-geom/9402007, Theorem A1 (source label `thm A1`), let `X` be a
smooth surface and

```text
B = sum b_j B_j,       0 <= b_j <= 1-epsilon < 1,
-(K_X+B) nef.
```

Then there is a constant depending only on `epsilon` such that

```text
rk Pic(X) <= A_1(epsilon).
```

The next theorem says that, after excluding the case `K_X` numerically trivial
and `B=0`, these smooth pairs form a bounded class.  The singular version fixes
`epsilon>0`, assumes `K_X+B` is MR `epsilon`-log terminal and
`-(K_X+B)` nef, and excludes the numerically trivial `B=0`, Du Val case; it
again concludes boundedness.

The dependence on `epsilon` is essential, not cosmetic.  Alexeev's example
immediately before the theorem starts with two lines in `P^2`, repeatedly
blows up the intersection of the first line with the newest exceptional
divisor, chooses its coefficient arbitrarily close to one, and obtains smooth
pairs with ample `-(K+B)` and strictly increasing Picard number.

There is a second noncircularity warning.  The proof of A1 explicitly reduces
to the conditions of “Nikulin's main theorem” on exceptional-curve graphs.
So Alexeev is an algebraic-geometric theorem, but its rank mechanism is the
same diagram technology transferred from hyperbolic geometry.  That is fine
for a new proof if its hypotheses can be derived geometrically; it is not an
independent black-box replacement for the missing reflection complexity bound.

### 2.3 A literal surface-boundary dessin is too small

For an snc divisor on a surface, the dual complex has dimension at most one.
It can encode components and intersection points, but not the full face poset
of an arbitrary `n`-polytope when `n>2`.  This rules out a literal uniform
“facet/face strata = boundary strata” construction on a surface.

The Neron--Severi lattice avoids this dual-complex obstruction, because its
rank is unbounded on rational blowups.  It then meets the opposite obstruction:
the effective/nef cone and its hidden walls are not determined by the abstract
lattice, and no uniform rank boundedness remains.

### 2.4 Complete-nef compact carriers have a sharp elementary obstruction

If a compact arithmetic Coxeter chamber is exactly the projectivization of
the **entire** nef cone of a smooth projective surface and all facet reflections
are integral on `NS(X)`, then its dimension is at most three.  Hodge index gives
signature `(1,rho-1)`, finite volume makes the reflection subgroup finite index
in the arithmetic orthogonal group, compactness makes the rational quadratic
space anisotropic, and Meyer makes every indefinite rational form in at least
five variables isotropic.  Hence `rho<=4` and `n<=3`.

This is a genuine `P -> surface -> theorem on the same P`, but only for the
direct full-nef carrier class.  It is also a no-go for a universal direct
surface compiler: known compact arithmetic examples in `H^8` cannot be exact
full-nef chambers of a surface.  An eigenspace/folded slice can evade Meyer,
but then exact return and the absence of hidden restricted walls become new
problems.

### 2.5 Replacing the surface does not presently close the gap

For `n>=3`, real hyperbolic space for `O(n,1)` is not a Hermitian symmetric
domain.  Its arithmetic quotient has no canonical Baily--Borel algebraic
variety analogous to an orthogonal `O(2,n)` Shimura variety.  Stabilizing the
lattice to signature `(2,n+1)` changes the problem and does not automatically
produce a reflective Borcherds form.  Wang's classification supplies explicit
hyperbolic 2-reflective lattices of rank at least six with no automorphic
correction, so an object-level Borcherds return cannot be universal.

One can canonically attach an algebraic group, quadratic space, arrangement,
or congruence tower to `P`.  Current boundedness theorems for these objects
require a bounded rank, volume, level, height, Hilbert polynomial, or divisor
degree.  Reflection maximality has not been translated into any one of those
bounds without importing the target arithmetic-volume machinery.

## 3. The successful “AG analogy”: small Gram graphs as arithmetic dessins

The real leverage inherited from surfaces is the diagram method.  It has four
separate pieces.

1. A closed finite-sided hyperbolic polyhedron has a facet of minimality `14`.
2. Edges in this facet give V-arithmetic edge chambers.
3. A field of sufficiently high degree forces the hyperbolic connected
   component of such an edge graph to have fewer than four vertices.
4. Four-vertex components have five finite graph shapes; their Gram
   determinants and all nongeometric embeddings constrain one algebraic
   parameter to explicit intervals.

The five marked graph shapes, their Coxeter labels, cyclic products, embedding
intervals, and determinant signs are a viable “Coxeter dessin.”  Fekete or a
purpose-built integer polynomial converts the interval dessin into a degree
bound.  This is effective and finite after the graph has been reduced to four
vertices.

It is important not to misname the output:

```text
narrow-face graph theorem  -> field-degree/field-list funnel
separate face-count theorem -> forces such a graph in high dimension
separate dimension theorem  -> rules out all sufficiently high dimensions.
```

The first two do **not** imply the third.

## 4. Exact hostile audit of the transition constant

### 4.1 Frozen primary artifacts

The obsolete arXiv artifact is:

```text
/private/tmp/nikulin_transition.pdf
SHA-256 e5391ab018e6e3f92ba9e2af8c35aff79332f973f3021f9c085148f52336f54b

/private/tmp/nikulin_transition.tar
SHA-256 2222824ada4f3d6ece3cab458da85f91b4b516cd76d754ae69492829c23bed85

/private/tmp/nik0910_src/transconst0.tex
SHA-256 02902c2e3a3bb745c7312cc7cdf119bd1354a2d0de4a7cb4bfaee5f9480b4b82
```

This is arXiv:0910.5217v4, 23 January 2010.

The final peer-reviewed artifact is:

```text
/private/tmp/nikulin_transition_journal.pdf
SHA-256 716e6e7333efbfaa22060b60dc97e8303363a4914409ccdd56ef7a90c923c7af
```

V. V. Nikulin, *The transition constant for arithmetic hyperbolic reflection
groups*, Izvestiya: Mathematics **75**:5 (2011), 971--1005,
<https://doi.org/10.1070/IM2011v075n05ABEH002561>.

### 4.2 The arXiv-v4 defect is real

On arXiv-v4 p.15, equation (24) prints

```text
P_bad(x)=(x+2)^3 (x+1)^4 x^4 (x-1)^3 (x^2+x-1)^2
```

and assigns it the maximum `delta=0.0905969664` on `[-2,1]`.  Exact
differentiation gives

```text
P_bad'(x)=
x^3(x-1)^2(x+1)^3(x+2)^2(2x+1)(x^2+x-1)
 * (9x^4+18x^3-10x^2-19x+8).
```

The actual maximum is approximately

```text
max[-2,1] |P_bad| = 0.1598126956594688...
```

and the displayed logarithmic estimate at `t=14` is

```text
1 + log|P_bad(14)|/(-log delta_bad)
  = 27.2225709073...,
```

so that argument gives degree at most `27`, not `23`.

The literal arXiv-v4 dependency is

```text
eq. (24), p.15
  -> Theorem 11, p.19: Gamma_1^(3), s=k=3, degree <=23
  -> Proposition 12, p.20: Gamma_1^(4), degree <=23
  -> Proposition 13, p.20: Gamma_2^(4), degree <=23
  -> Theorem 4, pp.7--8: N(14)<=25.
```

Thus the printed arXiv-v4 proof of that branch is invalid.  This observation
alone would not prove the theorem false; at worst it would replace the
preprint's advertised global bound by `27`.

### 4.3 The final journal version repairs it

Journal equation (4.19), p.986, changes the last exponent from `2` to `3`:

```text
P_good(x)=(x+2)^3 (x+1)^4 x^4 (x-1)^3 (x^2+x-1)^3.
```

This is exactly the translate `x -> x+2` of journal equation (4.8).  Its
derivative factors as

```text
x^3(x-1)^2(x+1)^3(x+2)^2(2x+1)(x^2+x-1)^2
 * (2x^2+2x-1)(5x^2+5x-8).
```

At both roots of `5x^2+5x-8`,

```text
P_good(x) = -884736/9765625,
max[-2,1] |P_good| = 884736/9765625 = 0.0905969664.
```

The bound becomes

```text
1 + log|P_good(14)|/(-log(884736/9765625))
  = 23.2493142903...,
```

and hence the integer degree is at most `23`.

Exact journal locators are:

- Theorem 3.2, p.977: `Gamma_1^(4),Gamma_2^(4)` degree at most `23`,
  `Gamma_3^(4),Gamma_4^(4),Gamma_5^(4)` degree at most `25`, hence
  `N(14)<=25`;
- equation (4.19), p.986: the corrected exponent-three polynomial;
- Theorem 4.6, p.989: the `Gamma_1^(3)` bounds, including `s=k=3`;
- Proposition 4.1, p.990: propagation to `Gamma_1^(4)`;
- Proposition 4.2, p.991: propagation to `Gamma_2^(4)`;
- Theorems 4.7, 4.8, 4.9, pp.994, 996, 998: the remaining three
  four-vertex families;
- closure of the proof of Theorem 3.2, p.998;
- Remark 5.1, p.1000: the older claim `N(14)<=56` had actually proved only
  `N(14)<=120`, and the present paper corrects it.

The phrase in Theorem 3.2 “one can take `N(14)=25` (or `N(14)<=25`)” should
be reported as an upper bound.  The definition makes `N(14)` an actual maximum;
equality would require a degree-25 witness, which is not supplied by this
sentence.

**Frozen verdict:** arXiv v4 contains a real defect; the journal proof repairs
that defect.  It is wrong to present the exponent error as a present gap in the
published theorem.

### 4.4 Why the project still must replay the journal theorem

Repairing (4.19) certifies only one dependency branch.  Theorem 4.6 contains
several other explicit-polynomial maxima over algebraic intervals, and Theorems
4.7--4.9 use further decimal extrema and finite parameter lists.  The article
also says that the five four-node graph shapes and their parameter lists are
“easy to draw.”  Those are finite claims, but they are nontrivial claims.

Therefore the project's acceptable independent certificate is:

1. enumerate all five four-node marked edge-graph shapes from the two spherical
   vertex-deletion conditions;
2. enumerate every allowed integer-label tuple;
3. derive every determinant and ground-field identity symbolically;
4. isolate every real embedding interval exactly;
5. certify every polynomial maximum by Sturm sequences/algebraic interval
   arithmetic, never by an unguarded decimal;
6. reproduce the integer degree bound for every branch;
7. take the exact maximum.

Until that replay is complete, `N(14)<=25` is a strong published benchmark,
not an imported axiom of this project.

## 5. Structural dependency audit: what is solid, false, repaired, and separate

### 5.1 Narrow face of minimality 14: short independent proof

Nikulin's 1981 paper, *On the classification of arithmetic groups generated by
reflections in Lobachevsky spaces*, proof of Fundamental Theorem 4.1.1,
journal p.113, chooses an interior point `O` and a facet whose visual angular
opening is minimal.  For two facet normals `e,f` of square `-2`, formula (2.1)
is

```text
e.f = g( sin((theta_1+theta_12)/2) sin((theta_2+theta_12)/2)
           / (sin(theta_1/2) sin(theta_2/2)) ),
g(x)=4x-2.
```

If both closed facets meet the minimal facet, then
`theta_12 <= theta_min <= theta_i`.  Each sine ratio is strictly below two,
so their product is below four and `e.f<g(4)=14`.  This is a clean,
surface-free lemma and can be adopted after writing the elementary visual-cone
argument carefully.

### 5.2 Hyperbolic subgraphs inherit the ground field

The 1981 Theorem 2.1.3(2) says that a hyperbolic principal subgraph of a
V-arithmetic Gram graph is V-arithmetic with the same ground field.  The key
reason can be replayed without authority: if the subgraph generated a proper
subfield, a nongeometric real embedding of the full ground field agreeing on
that subfield would leave its hyperbolic cyclic products geometrically
equivalent, while V-arithmeticity requires the conjugate principal form to be
negative definite.  Signature gives a contradiction.  The embedding-extension
and diagonal-rescaling steps should be included explicitly in a formal proof.

### 5.3 Four-node reduction: one false sentence with a repairable conclusion

Journal Proposition 3.1 reduces every edge chamber with a hyperbolic component
of at least four vertices to:

```text
FL^4, FT, or one of five families F Gamma_i^(4)(t).
```

For `u>2`, a connected four-node subgraph containing the distinguished pair
exists and the five-shape enumeration applies.  For `u=2`, any nongeometric
embedding would have to make a principal `2x2` parabolic matrix negative
definite, which is impossible; hence the ground field is `Q`.

The published proof of the `u<2` branch says that a minimal hyperbolic
subgraph is Lanner.  This implication is false.  The marked arithmetic
orthoscheme `[3,3,6]` has square-minus-two Gram matrix

```text
[-2, 1, 0, 0; 1,-2,1,0; 0,1,-2,sqrt(3); 0,0,sqrt(3),-2],
```

marked middle pair `u=1`, full signature `(1,3)`, and three elliptic vertex
deletions, while its fourth deletion is the affine `[3,6]` chain with null
vector `(1,2,sqrt(3))`.  It satisfies the finite-edge hypotheses and is
minimal indefinite but quasi-Lanner, not Lanner.  The exact replay is pinned
in the
[`[3,3,6]` dependency audit](nikulin_conceptual_dependency_hostile_audit.md).

The proposition's field-set conclusion is repairable.  A V-arithmetic
quasi-Lanner simplex is finite-volume noncompact arithmetic and therefore has
ground field `Q`; the compact branch is Lanner and lands in `FL^4` or `FT`.
Later dependencies use only this repaired field-set conclusion, so no new
numerical branch appears.  Any diagram- or AG-level return must nevertheless
retain the affine null root and cusp cluster rather than replacing it by a
compact negative-definite configuration.

### 5.4 High-dimensional face count: exact and independently checkable

For a compact arithmetic chamber in `H^n`, the chamber and its narrow facet are
simple.  The exact h-vector inequality used in arXiv:0708.3991 gives the average
number `a` of vertices of a two-face of the `(n-1)`-dimensional narrow facet:

```text
a < 4+4/(n-2),  n even;
a < 4+4/(n-3),  n odd.
```

If the field avoids the small graph sets, every connected component of each
finite vertex Coxeter graph has at most two nodes.  Counting nonright two-angles
both by vertices and by polygonal two-faces gives

```text
a * ((n-1)(n-2)/2 + floor((n-1)/2))
  >= 5(n-1)(n-2)/2.
```

For `n=10` in the even case and `n=11` in the odd case, the lower bound is
`4.5` while the upper bound is strict `<4.5`; the contradiction only gets
stronger afterwards.  Hence, for `n>=10`, the field must occur in the finite
small-diagram sets.  This calculation is short enough to certify independently.

For `n>=6`, the same average bound is `<5`, so the narrow facet has a triangle
or quadrilateral two-face.  The triangle and quadrilateral lemmas reduce the
field to the Lanner/triangle/five-core lists or to an arithmetic quadrilateral
field.  The printed lemma saying `P(M)={delta_1,...,delta_4}` is plainly meant
to say the polygon's wall set; the ambient `n`-polytope generally has more than
four facets.

### 5.5 This chain bounds fields, not dimensions

The output is:

```text
n>=10: field in FL^4, FT, or five four-node families;
n>=6:  additionally arithmetic quadrilateral fields;
n=4,5: reduce to n-1/n-2 or a small core;
n=3,2: separate low-dimensional inputs.
```

The published transition paper then yields degree at most `25` in dimensions
at least six and, with low-dimensional work, a global field-degree bound.

There is no implication

```text
bounded [k:Q]  => bounded n
```

without another theorem.  The classical `n<30` result is a separate
Vinberg theorem.  Likewise, global finiteness of maximal arithmetic reflection
groups uses three independent inputs:

1. a dimension cutoff;
2. a finite field/degree funnel;
3. fixed-`(n,[k:Q])` finiteness of maximal arithmetic reflection groups.

The narrow-face transition constant supplies only item 2.  It cannot be
advertised as an AG-derived or dessin-derived global dimension theorem.

## 6. The alleged “deeper Belolipetsky error”

No public primary statement by M. Belolipetsky alleging a structural failure of
the Nikulin narrow-face/transition chain was located in this audit.  The primary
record found points in the opposite direction:

- Belolipetsky, *Arithmetic hyperbolic reflection groups*,
  arXiv:1506.03111, Section 3, states Nikulin's edge-chamber theorem, the
  minimality-14 reduction, the correction history, and `N(14)<=25`, and calls
  the resulting all-dimensional degree-25 bound the best known general bound;
- Belolipetsky--Linowitz, *On fields of definition of arithmetic Kleinian
  reflection groups II*, arXiv:1211.4443, Corollary 5.1, explicitly combines
  their three-dimensional result with Maclachlan and Nikulin to obtain the
  all-dimensional degree bound `<=25`;
- Agol--Belolipetsky--Storm--Whyte's finiteness proof is independent, but its
  published introduction describes Nikulin's proof as another proof, not as a
  refuted one.

The frozen primary artifacts checked here are:

```text
/private/tmp/belolipetsky_survey.pdf
SHA-256 28f6a68799028258e023964c88f5d70cf01df79d45a96ab526c8e9e01c530200

/private/tmp/absw.pdf
SHA-256 974ac7c5141ca1c45b515f145215ad115aba867e42551a9b103233ea52930ee4

/private/tmp/nikulin_finiteness.pdf
SHA-256 9a820fef0b88158b283057f7d97acc89130610a9a4aacc6c53a8215e23ce874a
```

Exact locators make the attribution question decisive.  Belolipetsky's survey,
arXiv:1506.03111, printed p.11, Section 3, states Theorem 3.1, says that a
minimality-14 face exists without a cocompactness assumption, identifies the
published correction trail, records `N(14)<=25`, and then states Theorem 3.2,
the all-dimensional degree bound `<=25`.  ABSW, arXiv:math/0612132, printed
p.2, Introduction, says that Nikulin's contemporaneous proof uses his earlier
work and the dimension-two and dimension-three finiteness inputs; it calls the
two methods independent.  Nikulin's final finiteness paper, Izv. Math. 71:1
(2007), pp.53--56, puts the narrow-face induction in Theorem 2.1 and explicitly
uses Vinberg's separate `n<30` theorem in Theorem 2.2.

Accordingly, there is no source-backed Belolipetsky allegation to identify.
Our independent audit nevertheless found the false `u<2 => Lanner` sentence
and repaired its field conclusion in Section 5.3.  A separate conceptual
danger is a possible misinterpretation by us: the narrow-face/small-graph
chain ends in a finite field/degree funnel, not in a dimension theorem.
Appending

```text
finite field set -> global dimension bound
```

is invalid.  Nikulin's 2007 global finiteness argument does not make that
implication; it appends Vinberg's independent dimension cutoff.  This is exactly
why the present project may use the graph calculation as a benchmark but may not
advertise it as the desired AG-derived dimension theorem.

This does not disprove the existence of a private email, referee report, or
unindexed comment.  Any such allegation needs its exact text or source before
it can be attributed.  The remaining places that deserve hostile
certification are the five-shape completeness claim, the unreplayed algebraic
interval extrema, and a standalone formalization of the field-inheritance
lemma.  The quasi-Lanner defect now has an exact matrix certificate.

Do not confuse this with arXiv:math/0407313, M. N. Prokhorov,
*Boundedness of Degree of the basic field of an arithmetic reflection group in
Lobachevsky space*.  Its own arXiv record says that its main theorem has a
serious proof error and is not proved.  That is a different author and a
different proposed degree theorem.

## 7. A minimal independent completion plan

The useful theorem expedition now has two lanes.

### Lane A: certify the arithmetic-dessin field funnel

1. Reprove the visual-cone minimality-14 lemma exactly.
2. Implement an exhaustive generator for finite edge graphs with marked
   endpoints and spherical endpoint deletions; recover exactly five four-node
   shapes.
3. Prove field inheritance under hyperbolic principal subgraphs, including the
   `u=2` and quasi-Lanner cases.
4. Recompute every Gram determinant and embedding interval in the 2011 paper.
5. Replace every decimal maximum by an exact algebraic certificate.
6. Recompute `N(14)` as a certified upper bound.
7. Reprove the `n>=10` and `n>=6` face-count reductions from the exact h-vector
   inequality.

This lane can genuinely produce an effective finite list of high-dimensional
field dessins.  It still will not prove a dimension cutoff.

### Lane B: search for a genuinely new terminal geometric invariant

A surface/compiler theorem capable of replacing a dimension or global
finiteness theorem must prove, rather than assume, at least one of:

```text
absolute epsilon-lt gap;
absolute marked-facet/divisor degree;
absolute Picard rank;
absolute Borcherds compensation weight/height;
absolute level/monodromy degree;
absolute log volume bounded away from zero and above.
```

Reflection maximality also needs a rigorous AG translation, for example:
“there is no nontrivial reflection-compatible finite cover/contraction of the
marked pair.”  No boundedness theorem located here turns that minimality
condition alone into any of the absolute bounds above.

The most concrete falsifiable target is therefore:

> **Maximal-carrier discrepancy conjecture.** Construct, from every maximal
> arithmetic Coxeter chamber, an exact marked log pair with `-(K+B)` nef and
> prove an absolute lower bound on its log discrepancies directly from the
> absence of reflection-compatible AG extensions.

If true, Alexeev immediately gives a dimension bound.  If the construction
uses `1-1/m` coefficients, the conjecture must genuinely prove an absolute
label bound and cannot hide it in the definition.  A counterfamily with
maximal chambers and discrepancies tending to zero would kill this route
decisively.

## 8. Final status

- **New unconditional AG global dimension theorem:** none found.
- **Unconditional theorem for a meaningful carrier class:** an exact compact
  full-nef integral surface carrier has dimension at most three.
- **Best conditional AG route:** universal exact log-pair carrier plus an
  absolute discrepancy gap gives `n<=A_1(epsilon_0)-1`.
- **Missing noncircular input:** derive that gap, or another absolute marked
  complexity bound, from arithmetic reflection maximality.
- **Successful present leverage:** narrow-face four-node Gram dessins give an
  effective field funnel.
- **Current numerical status:** the journal correction validates the formerly
  defective `(4.19)` branch; the whole `N(14)<=25` theorem remains a published
  benchmark until all finite branches receive project-native exact
  certificates.
- **Dimension boundary:** still separate.  Neither surfaces nor the
  transition constant currently replace an independent dimension theorem.
