# Exact Coxeter local systems exist, but rigidity supplies no global AG bound

**Status (2026-08-12): sharp no-go for the standard orbicurve/nonabelian-Hodge
route, plus an exact conditional finiteness theorem.  No Vinberg, Nikulin,
ABSW, or arithmetic-reflection finiteness theorem is used.**

## 1. Result

There is a simple exact local-system carrier for every labelled arithmetic
hyperbolic Coxeter chamber.  If the chamber has `f` facets in `H^n`, put
`N=n+1`.  On a genus-zero orbicurve with `2f` order-two points, send two
consecutive orbifold meridians to each simple reflection.  The product
relation is then automatic.  The monodromy traces recover the normalized
Gram exactly:

\[
  g_{ij}^{2}=\frac{\operatorname {tr}(s_i s_j)-N+4}{4},
  \qquad g_{ij}\leq 0\quad(i\ne j).
\tag{1.1}
\]

After retaining the coefficient order, its invariant lattice, and the
distinguished real embedding, this is an integral orthogonal local system.
Tame nonabelian Hodge theory turns its semisimplification into a parabolic
Higgs bundle.  Thus this route really does give

\[
       P\longmapsto (\mathcal C_f,\mathbb V_P,Q_P,\sigma,
                       \{p_{i,+},p_{i,-}\})
       \longmapsto P .
\tag{1.2}
\]

The obstruction is not the existence or exactness of the carrier.  It is
the complete absence of the uniform boundedness and rigidity needed to push
back a finiteness theorem.

In fact the standard exact carrier is **never cohomologically rigid**.  For
an irreducible rank-`N` reflection representation its parabolic deformation
space has dimension

\[
 \dim H^1\!\left(\mathbf P^1,j_{!*}\operatorname {End}\mathbb V_P\right)
       =2(N-1)(2f-N-1).
\tag{1.3}
\]

Since every full-dimensional hyperbolic polytope has `f>=N`, this is at
least `2(N-1)^2`.  The compulsory controls give

\[
\begin{array}{c|c|c|c}
  \text{carrier}&N&f&\dim H^1\\ \hline
  \text{Borcherds }H^{21}/210&22&210&16674\\
  \text{Vinberg }X3\ H^{19}/36&20&36&1938\\
  \text{Lambert quadrilateral}&3&4&16.
\end{array}
\tag{1.4}
\]

Esnault--Groechenig therefore does not apply to this exact carrier.  Their
theorem proves **integrality from cohomological rigidity**; it does not prove
rigidity, bound rank, bound the curve, or bound a coefficient field.  In the
quasiprojective counting step they explicitly fix the rank and the orders of
the determinant and boundary eigenvalues before obtaining a finite count.

Nor would an unspecified `rigidification' finish the problem.  Polarized
geometric hypergeometric local systems on the fixed curve
`P^1-{0,1,infinity}` give:

* cohomologically rigid integral systems of unbounded rank (the Dwork
  hypergeometric systems); and
* cohomologically rigid rank-two systems whose coefficient fields have
  unbounded degree (cyclic-cover Gauss systems).

Consequently none of Simpson correspondence, Katz rigidity,
Esnault--Groechenig integrality, or an Arakelov inequality supplies the
missing uniform rank/field/presentation bound.  A successful version of this
route would require a new theorem saying that arithmetic reflection
maximality produces a uniformly bounded **rigid** exact local-system package.
That is not a consequence of the named AG results, and the Lambert family
shows that maximality alone cannot imply it.

## 2. The exact orbicurve carrier

Let `P` be a labelled arithmetic Coxeter chamber.  Write

\[
       s_i(v)=v-2\frac{q(v,e_i)}{q(e_i,e_i)}e_i,
       \qquad 1\le i\le f,
\tag{2.1}
\]

for its simple reflections on the arithmetic quadratic space `(V,q)` over
its totally real field `k`.  Arithmeticity gives a `W(P)`-stable lattice
`Lambda` over an order in `k`; replacing the order by its normalization and
the lattice by a commensurable one gives an `O_k`-integral representation.
For a strict functor one either includes this invariant lattice in the source
datum or retains the groupoid of commensurable invariant lattices.  This
minor integral-normalization choice is not used in the obstruction below.

Choose distinct points

\[
 p_{1,+},p_{1,-},\ldots,p_{f,+},p_{f,-}\in\mathbf P^1
\tag{2.2}
\]

and let `C_f` be the order-two root stack at these points.  Its orbifold
fundamental group has the presentation

\[
 \left\langle\gamma_{i,\pm}\ \middle|\
       \gamma_{i,\pm}^{2}=1,\ 
       \prod_{i=1}^f\gamma_{i,+}\gamma_{i,-}=1\right\rangle .
\tag{2.3}
\]

Define

\[
             \rho_P(\gamma_{i,+})=\rho_P(\gamma_{i,-})=s_i.
\tag{2.4}
\]

Every factor `s_i s_i` is the identity, so (2.4) is a representation.
All boundary monodromies have order two, hence are quasi-unipotent with one
`-1` eigenvalue.  Restriction of scalars gives a `Z`-local system of rank
`[k:Q]N`; retaining the `O_k` action and the distinguished embedding
`sigma:k->R` prevents restriction of scalars from mixing the geometric and
nongeometric conjugates.

### Proposition 2.1 (exact trace return)

The labelled local system, its coefficient action, and `sigma` recover the
normalized signed Gram of `P`.

### Proof

Normalize `q(e_i,e_i)=1` at the distinguished real place and put
`g_ij=q(e_i,e_j)`.  If `P_i(v)=q(v,e_i)e_i`, then `s_i=I-2P_i` and

\[
\begin{aligned}
 \operatorname {tr}(s_i s_j)
 &=\operatorname {tr}(I-2P_i-2P_j+4P_iP_j)\\
 &=N-4+4g_{ij}^{2}.
\end{aligned}
\tag{2.5}
\]

This proves (1.1).  The outward acute convention gives `g_ij<=0`, so the
distinguished real square root recovers `g_ij`, including finite angles,
parallel pairs, and dotted distances.  The labelled Lorentzian Gram and the
chosen time sheet recover the chamber.  \(\square\)

The construction is functorial on the labelled Coxeter category: a
label-preserving isometry conjugates the monodromy and its quadratic form.
For the finiteness question this loses nothing, since an unlabelled chamber
has only finitely many facet labellings.  A strict unlabelled version should
be formulated as a finite quotient groupoid; no claim that every permutation
of a fixed point configuration is induced by an automorphism of `P^1` is
being made.

This carrier is genuinely algebraic on the Betti side: it is a point of a
relative character stack with prescribed order-two conjugacy classes.  The
Corlette--Simpson--Mochizuki correspondence associates to a reductive such
point a tame parabolic Higgs bundle with vanishing parabolic characteristic
numbers.  This last correspondence is real-analytic, not an algebraic map
between the Betti and Dolbeault moduli spaces.  Also, a flat indefinite
quadratic form is not by itself a Hodge polarization: the resulting harmonic
bundle need not be a `C^*`-fixed system of Hodge bundles.

## 3. The exact carrier is never rigid

Let `U=P^1-{p_{i,+},p_{i,-}}` and `j:U->P^1`.  For an irreducible rank-`N`
local system on a punctured sphere, Katz's rigidity index is

\[
 \operatorname {rig}(\mathbb V)
 =\chi\!\left(\mathbf P^1,j_{!*}\operatorname {End}\mathbb V\right)
 =(2-r)N^2+\sum_{a=1}^r\dim Z(M_a),
\tag{3.1}
\]

where `r` is the number of punctures and `Z(M_a)` is the matrix centralizer
of the local monodromy.  Irreducibility gives `H^0=H^2=C`, hence

\[
        \dim H^1=2-\operatorname {rig}(\mathbb V).
\tag{3.2}
\]

A reflection has eigenspace multiplicities `(N-1,1)`, so

\[
                 \dim Z(s_i)=(N-1)^2+1=N^2-2N+2.
\tag{3.3}
\]

For the carrier (2.4), `r=2f`, and substitution gives

\[
\begin{aligned}
 \operatorname {rig}(\mathbb V_P)
  &=(2-2f)N^2+2f(N^2-2N+2)\\
  &=2N^2-4f(N-1),\\
 \dim H^1
  &=2-2N^2+4f(N-1)\\
  &=2(N-1)(2f-N-1).
\end{aligned}
\tag{3.4}
\]

This proves (1.3).  Passing from `End` to `End^0` does not change `H^1`,
because the scalar summand extends as the constant sheaf on `P^1`.

The size in (1.4) is not an artifact of a loose estimate.  It measures the
relations lost when the Coxeter representation is placed on a punctured
curve.  The curve group remembers that each `s_i` is an involution and that
the ordered product is one.  It does **not** impose all Coxeter relations
between the different `s_i`.  Deformations in (3.4) are allowed to leave the
Coxeter locus.

One can instead cut out the Coxeter relations inside a representation scheme
of the presented Coxeter group.  For a fixed lattice in dimension at least
three, hyperbolic local rigidity may then isolate the distinguished point.
That does not produce global boundedness: the rank, the number of generators,
the Coxeter matrix, and even the representation scheme all vary with `P`.
In dimension two the maximal Lambert family gives a positive-dimensional
Coxeter locus itself.

## 4. The compulsory controls

### 4.1 Borcherds `H21`

Here `N=22` and `f=210`.  Formula (3.4) gives

\[
             \dim H^1=2\cdot21\cdot(420-23)=16674.
\tag{4.1}
\]

The carrier is exact and all `420` local monodromies have order two.  Thus
the absolute-maximal high-rank control passes the construction but fails the
rigidity hypothesis by a huge, explicit margin.  No argument that starts
`apply Esnault--Groechenig to the natural Coxeter local system' can pass this
example.

### 4.2 Vinberg `X3`

For the rank-twenty, thirty-six-wall chamber,

\[
             \dim H^1=2\cdot19\cdot(72-21)=1938.
\tag{4.2}
\]

Mixed square-minus-two and square-minus-six wall lengths cause no problem
for the monodromy carrier: reflection matrices and (2.5) see both.  They also
do not improve rigidity.

### 4.3 The maximal Lambert family

For the compact four-sided family in `H^2`, `N=3`, `f=4`, so (3.4) gives
`dim H^1=16`.  More sharply, the exact trace coordinates include

\[
 \operatorname {tr}(s_1s_3)=-1+4x^2,
 \qquad
 \operatorname {tr}(s_2s_4)=-1+4y(x)^2,
\tag{4.3}
\]

and therefore trace a nonconstant algebraic curve in the character variety.
Outside finitely many parameters the reflection groups are absolutely
maximal.  Hence absolute reflection maximality does not imply local-system
rigidity, even when all local conjugacy classes and the base orbicurve are
fixed.

The generic Lambert member is not arithmetic.  The arithmetic members are
cut out by the additional integral equation

\[
                    (4x^2-4)(4y^2-4)=4.
\tag{4.4}
\]

Thus this control does not refute a future theorem using the full arithmetic
condition.  It does prove that the required discreteness cannot come from
maximality, the character stack, or nonabelian Hodge theory alone.

## 5. What the named theorems actually provide

### 5.1 Simpson and tame nonabelian Hodge theory

For a **fixed** smooth projective variety and a **fixed** rank, Simpson
constructs the Betti/de Rham/Dolbeault moduli spaces and identifies their
semisimple or polystable points through nonabelian Hodge theory.  Simpson's
noncompact-curve theorem and Mochizuki's tame correspondence add prescribed
parabolic data.  These are correspondences and fixed-numerical-data moduli
constructions, not theorems bounding the rank, number of punctures, or
parabolic denominators.

Moreover, only `C^*`-fixed Higgs points are systems of Hodge bundles.  The
orthogonal harmonic bundle attached to (2.4) is not automatically a
polarizable variation of Hodge structure, so an Arakelov inequality cannot
be applied without a new Hodge-lifting theorem.

### 5.2 Esnault--Groechenig

Their Theorem 1.1 says that an irreducible cohomologically rigid complex
local system with finite determinant and quasi-unipotent boundary monodromy
is integral.  Its logical direction is

\[
 \text{rigid + finite determinant + quasi-unipotent boundary}
      \Longrightarrow \text{integral}.
\tag{5.1}
\]

Arithmetic Coxeter monodromy already starts on the integral side; (5.1)
does not give the converse.  In particular it cannot turn the nonzero group
in (3.4) into zero.

The introduction of the same paper makes the bounded inputs especially
clear.  For the quasiprojective counting argument it first bounds the rank
`r` and fixes integers `d` and `h` controlling the determinant order and the
orders of the boundary eigenvalues.  Only then is the number
`N(r,d,h)` finite.  Those are precisely the kinds of uniform inputs the
Coxeter programme is trying to derive.

### 5.3 Katz rigidity

Katz's algorithm classifies irreducible cohomologically rigid local systems
on a punctured projective line from their local conjugacy data, using tensor
operations and middle convolution.  It is extraordinarily effective after
the rank and local data are given.  It does not bound those data as the
local system varies.  Section 6 gives explicit rigid geometric systems of
arbitrarily large complexity.

### 5.4 Arakelov inequalities

Arakelov--Jost--Zuo inequalities bound degrees or slope differences of Hodge
bundles by expressions of the form

\[
       \text{constant(weight)}\cdot
       \operatorname {rank}(E)\cdot
       \deg\Omega^1_C(\log D).
\tag{5.2}
\]

They require a polarized VHS and are homogeneous in the rank.  They do not
give an upper bound for the rank.  Even their base complexity grows in the
exact carrier, since `|D|=2f`.

Deligne-type finiteness over a finite field likewise fixes the rank and
bounds ramification.  Spreading a Coxeter local system out does not create
those bounds, nor does it control the degree or height needed to return the
distinguished real Gram.

## 6. Rigid hypergeometric systems are themselves unbounded

These two standard families prevent `rigid' from being used as a synonym
for `uniformly bounded'.

### 6.1 Dwork: unbounded rank on one three-punctured curve

For every integer `m>=3`, put `r=m-1` and consider the hypergeometric
differential equation

\[
 \left[\theta^r-t\prod_{a=1}^{r}
              \left(\theta+\frac a m\right)\right]u=0,
 \qquad \theta=t\frac d{dt},
\tag{6.1}
\]

on `P^1-{0,1,infinity}` (rescaling `t` only moves the finite singular
point).  It is the rank-`r` invariant primitive piece of the degree-`m`
Dwork pencil.  Hence it is of geometric origin and carries an integral
polarized variation of Hodge structure.

Its local monodromies have the following centralizer dimensions:

\[
       r,\qquad r,\qquad (r-1)^2+1.
\tag{6.2}
\]

The first two correspond to a regular block and a semisimple matrix with
distinct eigenvalues; the third is a pseudoreflection.  Therefore

\[
  (2-3)r^2+2r+(r-1)^2+1=2,
\tag{6.3}
\]

so the system is cohomologically rigid for every `r`.

Integrality is also visible in Levelt's companion-matrix description: the
two characteristic polynomials can be chosen as

\[
       \frac{X^m-1}{X-1}=X^{m-1}+\cdots+X+1,
       \qquad (X-1)^{m-1}.
\tag{6.4}
\]

Both have integral coefficients.  Thus a fixed algebraic curve supports
rigid polarized integral local systems of unbounded rank.

### 6.2 Cyclic covers: unbounded coefficient fields in rank two

For `m>=3`, consider

\[
        {}_2F_1\!\left(\frac1m,1-\frac1m;1;t\right).
\tag{6.5}
\]

It is a rank-two local system on the same three-punctured projective line.
All three non-scalar local centralizers have dimension two, and hence

\[
                    (2-3)\,2^2+2+2+2=2.
\tag{6.6}
\]

So it is cohomologically rigid.  It is a weight-one geometric system: an
eigenspace of the cohomology of the cyclic-cover family

\[
             C_{m,t}:\quad y^m=\frac{x(1-tx)}{1-x}.
\tag{6.7}
\]

Indeed the Euler period of `dx/((1-x)y)` is (6.5).  The cup-product
polarization on `H^1(C_{m,t})` induces its complex polarization.

The monodromy at infinity has eigenvalues `zeta_m,zeta_m^{-1}`, so its trace
generates

\[
                K_m=\mathbf Q(\zeta_m+\zeta_m^{-1}).
\tag{6.8}
\]

The companion polynomial

\[
       X^2-(\zeta_m+\zeta_m^{-1})X+1
\tag{6.9}
\]

is integral over `O_{K_m}`.  Along prime `m`,
`[K_m:Q]=(m-1)/2`, which is unbounded.  Restriction of scalars gives
polarized `Z`-local systems of rank `2[K_m:Q]`.

This family satisfies the qualitative hypotheses used in the integrality
theorem: irreducibility, cohomological rigidity, finite determinant, and
quasi-unipotent local monodromy.  The conclusion is integrality, not a bound
on `[K_m:Q]`.  A field bound requires a uniform bound on the local
eigenvalue order `m` (the `h` in Section 5.2) or some genuinely new input.

## 7. The exact conditional finiteness statement

There is a valid AG-native theorem after making the missing uniform package
explicit.

### Theorem 7.1 (fixed rigid orbicurve package)

Fix finitely many punctured curves `U_a`, a rank `r`, an integer `d`, and an
integer `h`.  For each `U_a`, consider irreducible complex local systems of
rank `r` whose determinant has order dividing `d`, whose boundary
eigenvalues have order dividing `h`, and which are cohomologically rigid.
There are only finitely many such isomorphism classes.  Consequently any
exact-return Coxeter compiler whose entire image lies in this package
returns only finitely many normalized Grams, fields, and labelled chambers.

### Proof

For fixed `U_a,r,d,h`, there are only finitely many possible boundary Jordan
types.  The corresponding relative character stacks are of finite type.
Cohomological rigidity makes the relevant points smooth and isolated in the
fixed-local-data stacks, so only finitely many occur.  Equivalently, this is
the finite count `N(r,d,h)` used by Esnault--Groechenig.  A finite union over
the `U_a` remains finite.  Exact return is injective on normalized labelled
Grams.  \(\square\)

This theorem does not solve the global Coxeter problem.  For the natural
carrier, the curve has `2f` orbifold points and (3.4) disproves rigidity.
Bounding the restriction-of-scalars rank would already bound
`[k:Q](n+1)`.  Thus treating Theorem 7.1 as a global proof would simply hide
the desired dimension/degree/presentation bounds in its hypotheses.

## 8. Exact missing theorem and verdict

To make the rigid-local-system strategy global, one needs a construction
different from (2.4) and a theorem of the following form:

> Every maximal arithmetic hyperbolic Coxeter chamber admits an exact
> orbicurve local-system carrier such that (i) the base orbicurves belong to
> a fixed finite list, (ii) the ranks and all parabolic denominators are
> uniformly bounded, and (iii) the carrier points are cohomologically rigid.

No cited theorem proves any of (i)--(iii) from reflection maximality.
Moreover:

* (3.4) disproves (iii) for the universal direct monodromy carrier;
* `H21` shows the failure at the mandatory absolute-maximal high-rank input;
* the Lambert curve disproves `maximality implies rigidity';
* Dwork systems disprove `rigid polarized integral implies bounded rank';
* cyclic-cover systems disprove `rigid polarized geometric implies bounded
  coefficient-field degree'; and
* Arakelov inequalities cannot repair either failure because rank is an
  input and the exact Coxeter harmonic bundle is not automatically a VHS.

Therefore the nonabelian-Hodge route presently supplies a clean exact
encoding and the scoped Theorem 7.1, but not an AG-native proof of global
finiteness of maximal arithmetic Coxeter polytopes or reflection groups.

## Primary references

1. C. Simpson, *Higgs bundles and local systems*, Publ. Math. IHES 75
   (1992), 5--95,
   <https://pmihes.centre-mersenne.org/articles/10.1007/BF02699491/>.
2. C. Simpson, *Moduli of representations of the fundamental group of a
   smooth projective variety I*, Publ. Math. IHES 79 (1994), 47--129,
   <https://numdam.org/articles/10.1007/BF02698887/>.
3. C. Simpson, *Harmonic bundles on noncompact curves*, J. Amer. Math. Soc.
   3 (1990), 713--770, DOI `10.1090/S0894-0347-1990-1040197-8`.
4. T. Mochizuki, *Kobayashi--Hitchin correspondence for tame harmonic
   bundles II*, Geom. Topol. 13 (2009), 359--455,
   <https://arxiv.org/abs/math/0602266>.
5. H. Esnault and M. Groechenig, *Cohomologically rigid local systems and
   integrality*, Selecta Math. 24 (2018), 4279--4292,
   <https://arxiv.org/abs/1711.06436>.
6. N. Katz, *Rigid Local Systems*, Annals of Mathematics Studies 139,
   Princeton University Press (1996),
   <https://www.degruyter.com/document/doi/10.1515/9781400882595/html>.
7. F. Beukers and G. Heckman, *Monodromy for the hypergeometric function
   nF(n-1)*, Invent. Math. 95 (1989), 325--354,
   <https://doi.org/10.1007/BF01393900>.
8. N. Katz, *Another look at the Dwork family*, in *Algebra, Arithmetic,
   and Geometry*, Progress in Mathematics 270 (2009), 89--126,
   <https://web.math.princeton.edu/~nmk/dworkfamilyfinal.pdf>.
9. P. Deligne and G. Mostow, *Monodromy of hypergeometric functions and
   non-lattice integral monodromy*, Publ. Math. IHES 63 (1986), 5--89,
   <https://www.numdam.org/item/PMIHES_1986__63__5_0/>.
10. J. Jost and K. Zuo, *Arakelov-type inequalities for Hodge bundles over
    algebraic varieties. I. Hodge bundles over algebraic curves*, J.
    Algebraic Geom. 11 (2002), 535--546,
    <https://arxiv.org/abs/math/0007102>.
