# Bounded universal mixed AG carriers: the exact no-go and the Mukai next species

**Status (2026-08-12).**  The requested unconditional global transfer theorem
is **not proved**.  There is, however, a non-tautological mixed K3/derived
carrier which passes both compulsory controls `X3` and `H21`.  What is missing
is not another choice of geometric vocabulary: it is a uniform
essential-surjectivity theorem, including a full-trace embedding and bounded
decoration complexity.  Such a theorem would already imply a stronger
rank/field bound before any moduli boundedness is invoked.

This report proves three precise statements.

1. A quasi-compact finite-type **decorated** carrier stack with exact,
   constructible extraction has finite Coxeter image.
2. Boundedness of the underlying varieties alone cannot be substituted for
   boundedness of the curve, graph, correspondence, or Fourier--Mukai-kernel
   decorations.
3. Vinberg's `X3` chamber has a mixed curve/Fourier--Mukai-action carrier, and
   Borcherds' absolutely maximal `H21` chamber has one as well.  Thus these
   tests do not kill the derived species; they isolate its unsolved uniform
   embedding and boundedness gates.

No Nikulin, Vinberg, Esselmann, or ABSW dimension/degree/finiteness theorem is
used to claim an AG finiteness result here.  The named `X3` and `H21` chamber
descriptions are calibration inputs only.

## 1. What a bounded carrier theorem would actually prove

Let `Cox` be the discrete groupoid whose objects are finite marked Coxeter
data (equivalently, exact wall Gram data with their chamber orientation), with
isomorphisms the marking changes.  Give every isomorphism class the discrete
topology.

### Theorem 1 (bounded decorated extraction)

Let \(\mathcal C\) be a quasi-compact algebraic stack and suppose an exact
carrier construction gives a map

\[
                         e:|\mathcal C|\longrightarrow |\mathrm {Cox}|.
\]

Assume either

1. \(e\) is locally constant, or
2. \(\mathcal C\) has a finite constructible stratification on every stratum
   of which \(e\) is constant.

Then only finitely many marked Coxeter data occur in \(\mathcal C\).

#### Proof

In the first case the inverse images of points are open.  They cover the
quasi-compact space \(|\mathcal C|\), hence a finite subcollection covers it.
In the second case there are only finitely many strata and each contributes
one value.  This proof uses no arithmetic reflection finiteness.  \(\square\)

The theorem is deliberately about the **decorated** stack.  For a mixed
carrier the datum to be bounded includes, at a minimum,

- a polarization or projective Hilbert polynomial for the variety;
- the number and Hilbert polynomials of marked wall curves/divisors;
- the number and graph Hilbert polynomials of regular or birational maps;
- for derived actions, a finite set of kernel numerical classes and a uniform
  Tor-amplitude/regularity bound; and
- a finite-type, constructible certificate that the selected hyperplanes are
  exact: no geometric wall or reflective action cuts the proposed chamber.

With these bounds, the usual Hilbert schemes and fixed numerical strata of
the moduli stack of perfect complexes provide the expected finite-type
parameter spaces.  The last, exact-cone condition is not automatic and is
therefore retained as a hypothesis.  Hacon--McKernan--Xu boundedness of stable
pairs with fixed dimension, coefficients and volume is an important source of
bounded underlying log pairs, but it does not supply bounded graphs/kernels or
the exact-cone predicate.

### Proposition 2 (the decoration barrier)

Boundedness of the unmarked underlying varieties does not imply the
conclusion of Theorem 1 for any carrier category which allows unrestricted
rational Chow correspondences, unrestricted birational graphs, or arbitrary
perfect kernels.

#### Proof

Already on a fixed rational surface whose N\'eron--Severi lattice has rank at
least two, divisor-product correspondences contain the full rational matrix
algebra on that lattice.  Matrices of unbounded entries therefore give
infinitely many numerical correspondence classes on the same surface.
Birational self-maps of the plane have arbitrarily large degree, so their
graph Hilbert polynomials are likewise unbounded.  Even for kernels of honest
autoequivalences, tensoring by \(L^m\) gives diagonal kernels with unbounded
numerical class as \(m\) varies.  Finally, Lieblich's moduli of complexes is
locally of finite presentation after numerical/boundedness data are fixed;
the union over all such data is only a locally finite-type union, not a single
quasi-compact bounded stratum.

Thus the forgetful image in the moduli of varieties can be one point while
the decoration lies in infinitely many numerical components.  A boundedness
theorem for the forgotten object contains no statement about those
components.  This is exactly the mechanism in the universal weak Chow-action
carrier: matrices and wall cuts encode the answer, while their complexity is
unbounded.  \(\square\)

This proposition is a logical no-go for **underlying-object boundedness**, not
a claim that no ingenious bounded decorated stack can exist.

## 2. `X3` has a genuine mixed curve/derived-action carrier

Let

\[
 L_{X_3}=U\oplus2E_8(-1)\oplus A_2(-1)=\operatorname {NS}(X_3).
\]

The frozen `X3` packet has 36 facets:

\[
 24\times(-2,\operatorname {div}1),\qquad
 12\times(-6,\operatorname {div}3).
\]

The 24 short-norm classes are smooth rational curve walls.  The 12 remaining
classes are internal divisor cuts, so `X3` fails a literal full-nef K3 return.
Nevertheless their reflections are not merely formal Chow operators.

### Proposition 3 (Fourier--Mukai lift of every long `X3` reflection)

For every primitive \(v\in\operatorname {NS}(X_3)\) with
\(v^2=-6\) and \(\operatorname {div}(v)=3\), the integral reflection

\[
                  r_v(x)=x+\frac{(x,v)}3v
\]

on \(\operatorname {NS}(X_3)\) is induced by a Fourier--Mukai
autoequivalence of \(D^b(X_3)\).

#### Proof

The discriminant group of \(L_{X_3}\) has order three.  The nonzero class
\(v/3\) therefore generates \(A_{L_{X_3}}\), and

\[
                         r_v(v/3)=-v/3.
\]

Hence \(r_v\) acts as \(-1\) on the discriminant group.  The transcendental
lattice \(T(X_3)\) has rank two and discriminant three.  The actions \(r_v\)
on \(L_{X_3}\) and \(-\mathrm{id}\) on \(T(X_3)\) agree through the
anti-isometry of discriminant forms defining the unimodular K3 lattice, so
they glue to an integral isometry of \(H^2(X_3,\mathbf Z)\).

Extend it by the identity on the Mukai hyperbolic plane
\(H^0\oplus H^4\).  The extension is a Mukai Hodge isometry: it acts by
\(-1\) on \(H^{2,0}\), and it preserves the Hodge decomposition.  It also
preserves the natural orientation of the positive four-plane.  Indeed,
\(-\mathrm{id}\) has determinant \(+1\) on the positive real period plane,
while \(r_v\) fixes a positive two-plane inside the algebraic orthogonal
\(v^\perp\).

The Mukai--Orlov derived Torelli theorem, with the orientation refinement of
Huybrechts--Macr\`i--Stellari, identifies the
cohomological image of Fourier--Mukai autoequivalences with the
orientation-preserving Mukai Hodge isometries.  The extension is therefore
induced by a Fourier--Mukai autoequivalence and restricts to \(r_v\) on
\(\operatorname {NS}(X_3)\).  \(\square\)

For a \((-2)\)-curve \(C\), the spherical twist in
\(\mathcal O_C(-1)\) has Mukai vector \((0,[C],0)\) and induces its root
reflection.  Consequently all 36 `X3` facet reflections have actual derived
actions; 24 also have intrinsic curve walls.  The complete pinned chamber is
recovered from their \((-1)\)-eigenlines and the chosen positive component.

This improves the weak “12 arbitrary divisors” interpretation.  It does not
prove absolute reflection maximality of `X3`; the frozen source only proves
completeness for its pinned integral host.

## 3. `H21` passes the same mixed species

For

\[
 L_{21}=U\oplus D_{20}(-1)
\]

the 210-wall Borcherds chamber is finite-volume and absolutely
reflection-maximal.  Its characteristic-two Artin-invariant-one
supersingular K3 realization has

\[
 42\times(-2,\operatorname {div}1)
       +168\times(-4,\operatorname {div}2).
\]

The 42 short classes are smooth rational curves.  Dolgachev--Kond\=o's 168
involutive Cremona transformations are regular automorphisms of the K3
surface and act by the 168 long reflections.  Pullback/pushforward by a
regular automorphism is a Fourier--Mukai autoequivalence (its kernel is the
structure sheaf of the graph).  The 42 curve reflections are again spherical
twists.  Hence `H21`, like `X3`, is a mixed K3 curve/Fourier--Mukai-action
carrier.

This also explains why the known OG10 route is not the universal repair.
Although all 210 input vectors embed as OG10 wall-divisor types, an additional
\((-4,\operatorname {div}1)\) wall cuts the strict interior.  The OG10
K\"ahler chamber is therefore not the Borcherds chamber.

### Corollary 4 (answer to the two-control question)

There is a non-tautological mixed species covering both controls without any
imported arithmetic dimension bound: K3 surfaces in arbitrary characteristic,
with curve walls and actual Fourier--Mukai actions.

For these two individual objects one may choose polarizations and the finite
sets of curves and kernels, then take the finite union of their Hilbert/kernel
strata; that union is bounded.  This observation is **not** a bounded
universal theorem.  The stack of all polarized K3 surfaces over all degrees,
decorated by arbitrary Fourier--Mukai kernels, is not the required single
quasi-compact stratum.

## 4. Why Mukai rank 24 does not yet give the global theorem

For a projective K3 surface,

\[
 \operatorname {rk}K_{\mathrm {num}}(X)
       =\rho(X)+2\le
       \begin{cases}
          22,&\text{characteristic zero or finite height},\\
          24,&\text{supersingular characteristic }p.
       \end{cases}
\]

Let \(P\subset\mathbb H^n\) have totally real ground field \(k\) of degree
\(d\), and let \(T_P\) be its integral full-trace lattice.  Any faithful
full-trace Mukai realization forces

\[
                  \boxed{d(n+1)=\operatorname {rk}_{\mathbf Z}T_P\le24.}
\]

Thus essential surjectivity into this species would itself prove a strong
dimension/field inequality.  It cannot be justified by first importing a
classical rank bound without making the proposed AG proof circular.

Rank is only the first gate.  One must also prove a primitive embedding of
the trace/order module with the correct discriminant gluing, extend the
\(\mathcal O_k\)-action to the relevant Hodge or crystalline realization,
lift every side reflection to an actual autoequivalence, and exclude all
hidden reflective actions/walls.  None follows from the inequality alone.
The supersingular K3 lattices are \(p\)-elementary; in particular, the frozen
golden-field audit proves that no supersingular K3 N\'eron--Severi lattice
admits the required integral self-adjoint golden full-NS action.  A Mukai
extension may evade that particular lemma, but no general construction is
known.  Positive characteristic is therefore a genuine extra branch, not a
universal solution to non-rational ground fields.

The same issue remains for motives and higher-dimensional hyperk\"ahler
varieties.  Rational Chow motives permit full matrix algebras and so fall
under Proposition 2 unless correspondence complexity is bounded.  For the
known hyperk\"ahler deformation types the cohomological ranks are small, but
there is neither a classification of all deformation types nor a total
trace-lattice/MBM-cone compiler.  Coherent full-\(H^2\) real multiplication
even satisfies the stronger period-slot inequality \(d(n+3)\le b_2\).

## 5. The exact missing theorem

The strongest viable next species is the following **bounded mixed Mukai
cone/action stack**.  An object consists of

\[
 (X,H,\iota; C_1,\ldots,C_a;Phi_1,\ldots,\Phi_b),
\]

where \(X\) is a projective K3 surface (allowing all characteristics), \(H\)
is a polarization, \(\iota:T_P\hookrightarrow K_{\mathrm {num}}(X)\) is a
primitive trace/order marking, the \(C_i\) are curve walls, and the
\(\Phi_j\) are Fourier--Mukai autoequivalences whose cohomological
\((-1)\)-eigenlines supply the action facets.

The required theorem, presently unavailable, is:

> **Uniform Mukai essential-surjectivity conjecture.**  There is one tuple of
> bounds \(B=(D,F,\mathcal H,\mathcal V,[a,b])\), independent of \(P\), such
> that every absolutely reflection-maximal finite-volume arithmetic Coxeter
> chamber has such an object with \(H^2\le D\), at most \(F\) decorations,
> divisor/graph Hilbert polynomials in the finite set \(\mathcal H\), kernel
> numerical classes in the finite set \(\mathcal V\), Tor amplitude in
> \([a,b]\), and exact chamber return after exhausting every curve, movable,
> and cohomological reflection wall meeting the marked positive cone.

If this conjecture were proved and exact extraction were constructible, its
stack would be quasi-compact and Theorem 1 would give unconditional global
finiteness without classical reflection finiteness.  Conversely, merely
calling the unbounded union a “stack of K3 categories” supplies none of the
conjecture.

There are five currently independent missing assertions:

1. `TRACE`: \(T_P\) embeds in the numerical Mukai lattice, equivariantly for
   the full ground-field order;
2. `LIFT`: every simple side reflection is induced by a curve twist,
   automorphism, birational action, or Fourier--Mukai kernel;
3. `EXACT`: their fundamental chamber is exactly \(P\), with no hidden
   curve/MBM/derived reflection;
4. `BOUND`: polarizations, number of facets, graphs and kernel classes have
   uniform bounds independent of \(P\); and
5. `CONSTR`: exact extraction is constant on a finite constructible
   stratification.

`X3` and `H21` prove that `LIFT` is plausible and genuinely geometric in two
hard endpoint cases.  They do not prove `TRACE`, `EXACT`, or `BOUND` for an
arbitrary maximal chamber.  In particular, no presently known bounded
K3/Mukai rank-24 category is total on the requested input class.

## 6. Candidate audits and claim boundary

| candidate | `X3` | absolute-maximal `H21` | class-wide bounded/total? |
|---|---|---|---|
| full K3 nef cone | fails: 12 internal cuts | fails: 168 action facets | no |
| mixed K3 curve + FM action | passes (Proposition 3) | passes | no theorem |
| rational/Enriques surface | no universal rank/return | fails known complex routes | no |
| OG10 full K\"ahler cone | not needed | fails by an explicit hidden wall | no |
| unrestricted Chow/motive/kernel | formally passes | formally passes | no leverage |
| bounded mixed Mukai stack | would pass | would pass | exact missing conjecture above |

Accordingly the literal user acceptance criterion is not met.  The strongest
rigorous conclusion is a positive two-control derived carrier theorem plus a
sharp bounded-decoration reduction.  A global finite-P theorem still awaits
the uniform Mukai essential-surjectivity conjecture or another species with
equally strong `TRACE/LIFT/EXACT/BOUND/CONSTR` gates.

## 7. Primary sources and frozen local dependencies

- Huybrechts--Macr\`i--Stellari,
  [*Derived equivalences of K3 surfaces and orientation*](https://arxiv.org/abs/0710.1645),
  Duke Math. J. 149 (2009), 461--507: orientation and the cohomological image
  of Fourier--Mukai equivalences.
- Orlov,
  [*Equivalences of derived categories and K3 surfaces*](https://arxiv.org/abs/alg-geom/9606006):
  representability and the K3 derived Torelli theorem.
- Bridgeland,
  [*Stability conditions on K3 surfaces*](https://arxiv.org/abs/math/0307164):
  the complex K3 stability/derived setting.
- Lieblich--Olsson,
  [*Fourier--Mukai partners of K3 surfaces in positive characteristic*](https://arxiv.org/abs/1112.5114),
  Ann. Sci. \`Ec. Norm. Sup\'er. 48 (2015), 1001--1033: positive-characteristic
  derived Torelli context.  No characteristic-two stability-manifold theorem
  is imported here.
- Lieblich,
  [*Moduli of complexes on a proper morphism*](https://arxiv.org/abs/math/0502198):
  locally finite-presentation moduli of fixed derived objects.
- Hacon--McKernan--Xu,
  [*Boundedness of moduli of varieties of general type*](https://arxiv.org/abs/1412.1186):
  bounded stable pairs with fixed data.
- Dolgachev--Kond\=o,
  [*A supersingular K3 surface in characteristic 2 and the Leech lattice*](https://arxiv.org/abs/math/0112283):
  the 42 curves and characteristic-two action realization.

The exact `X3`, `H21`, OG10, weak-Chow, supersingular-field and coherent-Hodge
claims are imported from the following frozen local reports and pinned by the
replay certificate:

```text
results/universal_weak_chow_action_carrier_and_leverage_no_go.md
results/projected_k3_universality_rank_saturation_no_go.md
results/borcherds_h21_outer_routing_counterexample.md
results/borcherds_d4_h21_og10_carrier_obstruction.md
results/supersingular_k3_golden_h10_obstruction.md
results/coherent_full_h2_real_multiplication_bound.md
```

The replay checks dependency hashes and the theorem/gate matrix.  It does not
machine-prove derived Torelli or download external PDFs.
