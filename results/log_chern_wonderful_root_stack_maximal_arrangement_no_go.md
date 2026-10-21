# Log Chern, wonderful-model, root-stack, and K-moduli tests for the exact Coxeter arrangement carrier

**Status (2026-08-12): sharp scoped no-go for the standard log-Chern/BMY,
BAB, K-stability, wonderful-model, and root-stack routes; plus a new exact
uniformly K-stable log-Fano countercontrol.  This is not a counterexample to
finiteness of maximal arithmetic reflection groups and not a proof of that
finiteness.**

## 1. Outcome

For the exact carrier

\[
 (\mathbf P(V),Q_q,H_1,\ldots,H_f,\sigma,C_P),
\]

put

\[
 B_{a,b}=a\sum_{i=1}^fH_i+bQ_q.
\]

Its canonical class is

\[
 K_{\mathbf P^n}+B_{a,b}
   =\bigl(af+2b-(n+1)\bigr)H.                         \tag{1.1}
\]

No standard log-geometric inequality turns arithmetic reflection maximality
into the missing uniform bounds.  There are four independent reasons.

1. The mandatory controls occupy incompatible positivity regimes.  For the
   natural reflection coefficient \(a=1/2\), Borcherds'
   \(\mathbf H^{21}/210\) carrier is deep in log general type, while the exact
   Vinberg \(X3\) carrier is log Calabi--Yau after the quadric is included.
2. Log Chern and arrangement Miyaoka--Yau expressions are intersection-poset
   invariants.  They do not see the quadric-relative Gram coordinates, the
   Vinberg field, or arithmetic height.
3. Crepant wonderful and root-stack modifications preserve (1.1), rather
   than changing a general-type carrier into a log-Fano one.  Noncrepant
   wonderful models store the unbounded arrangement-flat count in their
   exceptional boundary.
4. A simple root stack has abelian inertia at a crossing and therefore does
   not faithfully encode a Coxeter label \(m>2\).  A faithful dihedral stack
   stores \(m\) in stabilizer order, so bounded stack complexity would already
   imply the label bound one is trying to prove.

There is a particularly sharp new obstruction.  The maximal Lambert family
admits an **exact** klt log-Fano carrier

\[
 \left(\mathbf P^2,
       \frac12\sum_{i=1}^4H_{i,x}+\frac14Q_x\right)       \tag{1.2}
\]

such that every member is uniformly K-stable, all BAB data are fixed, the
reflection group is absolutely maximal off a finite parameter set, and the
exact Gram nevertheless moves in a curve.  Rational maximal members have
unbounded Gram height.  Thus even exact K-stability and BAB boundedness do not
make the return image finite.

The missing result remains an arithmetic theorem:

> **MAX--ARR.** Absolute arithmetic reflection maximality forces the exact
> arrangement carriers into finitely many bounded strata and forces the
> integral admissible Gram map on those strata to have finite image (or gives
> uniform degree and height bounds).

None of the log-geometric tools tested here proves MAX--ARR.

## 2. The canonical-class trichotomy is birationally persistent

Since \([Q_q]=2H\), \([H_i]=H\), and
\(K_{\mathbf P^n}=-(n+1)H\), formula (1.1) is immediate.  In particular,
\(B_{1/2,0}\) is of log-Fano numerical type only if

\[
                         f<2(n+1).                         \tag{2.1}
\]

More generally, if every facet coefficient is at least
\(\delta>0\) and the quadric has coefficient \(b\geq0\), log-Fano
positivity gives only the conditional estimate

\[
                         \delta f+2b<n+1.                  \tag{2.2}
\]

This is a useful scoped facet-ratio lemma, but it does not produce its own
hypotheses.  Choosing coefficients of order \((n+1)/f\) merely rewrites the
missing bound as the missing positive coefficient floor.

The four controls give

| carrier | \(n\) | \(f\) | \(K+B_{1/2,0}\) | \(K+B_{1/2,1}\) |
|---|---:|---:|---:|---:|
| Borcherds \(\mathbf H^{21}\) | 21 | 210 | \(83H\) | \(85H\) |
| Vinberg \(X3\) | 19 | 36 | \(-2H\) | \(0\) |
| Lambert | 2 | 4 | \(-H\) | \(H\) |
| Farey \(P_N\) | 2 | \(N+2\) | \((N/2-2)H\) | \((N/2)H\) |

The Borcherds side group is arithmetic and absolutely reflection-maximal.
Hence maximality does not imply (2.1), even before checking klt
singularities.  The exact \(X3\) pair is log Calabi--Yau.  The Farey sequence
is arithmetic and finite-area but nonmaximal; it shows that arithmeticity and
finite volume alone allow the boundary count and stable volume to escape.

Now let \(\pi:Y\to\mathbf P^n\) be any wonderful-model blowup and choose the
crepant boundary \(B_Y\) whenever it exists.  By definition,

\[
                         K_Y+B_Y=\pi^*(K_{\mathbf P^n}+B). \tag{2.3}
\]

For example, if a codimension-\(c\) flat contained in boundary components of
coefficients \(a_i\) is blown up, the crepant exceptional coefficient is

\[
                       \sum_i a_i-(c-1).                   \tag{2.4}
\]

Thus a crepant wonderful model preserves the sign, Iitaka type, and top
log volume of (1.1).  It cannot put the Borcherds carrier into BAB.  If one
instead gives all exceptional divisors a convenient noncrepant coefficient,
the number and incidence of those divisors are controlled by the building
set of arrangement flats.  No uniform flat-count bound follows from the
construction or from reflection maximality.

## 3. An exact uniformly K-stable Lambert carrier

Use the compact Lambert family \(P_x\subset\mathbf H^2\), \(x>1\), with
angle sequence

\[
                    (\pi/2,\pi/2,\pi/2,\pi/3).
\]

Its exact carrier consists of a smooth conic \(Q_x\subset\mathbf P^2\) and
four polar lines \(H_{i,x}\).  The conic is transverse to every line, no
three lines are concurrent, and no conic--line--line triple point occurs.
The normalized Gram varies nontrivially with (x).  Only finitely many
isometry classes in the family have a proper discrete reflection overgroup,
so absolute reflection maximality is cofinite on the parameter curve.

Set

\[
 \Delta_x^{\mathrm F}
   =\frac12\sum_{i=1}^4H_{i,x}+\frac14Q_x.                \tag{3.1}
\]

> **Theorem 3.1 (exact uniformly K-stable moving carrier).**  Every pair
> \((\mathbf P^2,\Delta_x^{\mathrm F})\) is \(1/2\)-lc, klt log Fano, and
> uniformly K-stable.  It has fixed coefficient set, boundary degrees,
> incidence data, polarization, and anti-log-canonical volume
>
> \[
>       -\bigl(K+\Delta_x^{\mathrm F}\bigr)=\frac12H,
>       \qquad
>       \bigl[-(K+\Delta_x^{\mathrm F})\bigr]^2=\frac14.  \tag{3.2}
> \]
>
> With the usual real-place and side marking it returns \(P_x\) exactly.
> Its image in log-Fano K-moduli is positive-dimensional modulo a finite
> relabelling action, and its absolutely reflection-maximal locus is
> cofinite.

### Proof

The SNC statement gives klt.  The smallest log discrepancy along a boundary
component is \(1-1/2=1/2\), and blowups of the transverse double points have
larger log discrepancy, so the pair is \(1/2\)-lc.  Equation (3.2) follows
from

\[
 -3H+\frac12(4H)+\frac14(2H)=-\frac12H.
\]

First omit the conic and put

\[
       \Delta_0=\frac12\sum_iH_{i,x},
       \qquad L=-\bigl(K+\Delta_0\bigr)=H.
\]

Fujita's hyperplane-arrangement criterion shows that this four-line pair is
uniformly K-stable.  Alternatively, its two strict inequalities are
(2>3/2) and (4>3).

Apply Fujita's perturbation proposition to

\[
                   B=\frac14Q_x\sim_{\mathbf Q}\frac12L.
\]

Here \(r=1/2\).  Since \(Q_x\) is smooth and transverse to the line boundary,

\[
       \operatorname {lct}(\mathbf P^2,\Delta_0;B)=4
            >2=r^{-1}.                                  \tag{3.3}
\]

The proposition makes \((\mathbf P^2,\Delta_0+B)\) uniformly K-stable.

The coefficient \(1/4\) uniquely identifies the conic and the coefficient
\(1/2\) components identify the four lines.  Quadric polarity therefore
recovers their root lines, and the real marking recovers the chosen chamber.
An isomorphism of the pairs is an isomorphism of exact carriers.  The
Lambert Gram parameter has only a finite relabelling ambiguity, so the
moving Gram gives a positive-dimensional K-moduli image.  Cofinite maximality
is the Lambert overgroup theorem.  \(\square\)

For \(x=N\in\mathbf Z_{\geq2}\), the carrier has a rational presentation and
the invariant two-cycle \(4x^2=4N^2\) has unbounded height.  All but finitely
many such members are absolutely reflection-maximal.  These groups are not
arithmetic: the other dotted cycle is

\[
                      4+\frac1{N^2-1},
\]

which violates Vinberg integrality.  This qualification is essential.  The
family proves that maximality plus exact BAB/K-moduli data gives neither
discreteness nor a height bound; arithmeticity must add a genuinely
Diophantine input.

## 4. The log-Chern inequalities are compatible with the escapes

For an SNC weighted surface pair \(B=\sum a_iD_i\) with no triple points,
the stratified orbifold Euler number is

\[
 c_2^{\rm orb}(X,B)
   =e(X)-\sum_i a_i e(D_i)
       +\sum_{i<j}a_ia_j(D_i.D_j).                         \tag{4.1}
\]

This is the second parabolic/log Chern number in the SNC BMY formula.
Take a smooth conic of coefficient \(b\) and \(f\) general-position lines of
coefficient \(1/2\), all transverse and with no triple point.  Then

\[
\begin{aligned}
 c_2^{\rm orb}
   &=3-2b-f+bf+\frac{f(f-1)}8,\\
 (K+B)^2
   &=\left(2b+\frac f2-3\right)^2.                       \tag{4.2}
\end{aligned}
\]

For the stable exact Lambert pair \(b=1,f=4\),

\[
              (K+B)^2=1,
              \qquad c_2^{\rm orb}=\frac52,              \tag{4.3}
\]

so logarithmic BMY is strict, \(1<3(5/2)\), at every parameter.  More
generally, for \(b=1\),

\[
 3c_2^{\rm orb}-(K+B)^2
       =\frac{f^2+5f+16}{8}>0.                            \tag{4.4}
\]

Thus this inequality is compatible with arbitrarily many transverse lines;
it supplies no upper bound on (f).  For the exact log-Fano carrier (3.1),

\[
        (K+\Delta_x^{\mathrm F})^2=\frac14,
        \qquad c_2^{\rm orb}=1,                           \tag{4.5}
\]

again independently of (x).  BMY is not invoked in the anti-ample
direction; (4.5) records that the Chern packet is constant while the exact
Gram moves.

The higher-dimensional arrangement inequality of de Borbon--Panov is even
more explicit about what it sees.  For \(f\) generic hyperplanes in
\(\mathbf P^n\) with equal weight \(a\), its quadratic form is

\[
 Q_{\mathcal H}(a,\ldots,a)
   =2(n+1)fa^2-2(fa)^2
   =2fa^2(n+1-f).                                         \tag{4.6}
\]

When (f>n+1), the equal-weight arrangement is stable and (4.6) is
negative.  In particular, taking (f=n+2) gives

\[
                       Q_{\mathcal H}=-2(n+2)a^2<0        \tag{4.7}
\]

in every dimension.  On the four Lambert lines, (n=2,f=4,a=1/2), it is
the constant value (-2).  The theorem depends only on the intersection
poset; it cannot see the moving conic-relative Gram.  Equations
(4.6)--(4.7) also show directly that its Miyaoka--Yau inequality gives no
absolute dimension or facet upper bound.

Base field degree and arithmetic height are invisible for the same reason.
Log Chern classes are unchanged by extension of the ground field.  Applying
the same inequality to every Galois conjugate repeats the same numerical
inequality; summing over embeddings multiplies both sides by the field
degree and cancels it.  A separate norm contraction or height inequality,
such as the one available on the Lambert **arithmetic** locus, is required.

## 5. BAB and K-stability stop at boundedness

BAB fixes dimension and a positive discrepancy gap before it concludes
boundedness.  The universal Coxeter carrier supplies no dimension bound.
Even within fixed dimension, (2.1) excludes the required
\(\mathbf H^{21}/210\) control from the natural log-Fano lane.

Theorem 3.1 shows the independent second failure.  Its family has

* fixed dimension (2);
* a fixed (1/2)-lc gap;
* fixed rational coefficient set ({1/2,1/4});
* fixed anti-log-canonical volume (1/4);
* uniform K-stability member by member; and
* exact Coxeter return.

BAB correctly places the family in a bounded parameter space.  K-moduli
correctly records a positive-dimensional family of nonisomorphic stable
pairs.  Neither theorem says that a bounded moduli space has finitely many
algebraic points, bounded field degree, bounded height, or a zero-dimensional
Gram image.

Consequently there are two logically separate hypotheses in any valid
log-Fano pushback theorem:

1. show that all maximal arithmetic Coxeter carriers enter finitely many
   fixed BAB/K-moduli species; and
2. show that arithmetic maximality cuts their exact Gram maps to finite
   image.

Current BAB and K-stability results prove neither implication.

## 6. Root stacks do not supply the Coxeter labels

Let \(D_i,D_j\) meet transversely and take order-two roots along both.  The
inertia at the crossing is

\[
                         \mu_2\times\mu_2.                \tag{6.1}
\]

The two local involutions commute.  This faithfully models only
\(I_2(2)\).  For a Coxeter angle \(\pi/m\), the local reflection group is the
dihedral group \(I_2(m)\) of order \(2m\), and the product of its generating
involutions has order \(m\).  For \(m>2\), (6.1) forgets the label.

A faithful dihedral quotient has, in suitable coarse coordinates, branch
discriminant

\[
                         v^2-4u^m=0.                     \tag{6.2}
\]

For odd \(m\) it is irreducible; for even \(m>2\) it has tangent rather than
transverse branches.  It is therefore not the simple SNC root construction
on the two projective carrier hyperplanes.  Ultraparallel facets have no
finite dihedral local group at all; their dotted distance must remain a
separate quadric-relative scalar.

There is also a precise boundedness circularity.  In a quasi-compact
finite-type family of Deligne--Mumford stacks with finite inertia, geometric
stabilizer lengths are uniformly bounded.  If a chart faithfully contains
\(I_2(m)\), this immediately bounds \(2m\).  But no construction proves that
the faithful Coxeter stacks lie in such a bounded family.  Establishing
that assertion is already the desired uniform label/stack-complexity theorem.

For an SNC root stack \(\rho:\mathcal X\to Y\), the canonical bundle formula
is

\[
 K_{\mathcal X}
  =\rho^*\!\left(K_Y+\sum_i(1-1/r_i)D_i\right).           \tag{6.3}
\]

Thus order-two roots recover the half-boundary in (1.1) and preserve its
positivity type.  They do not turn the (83H) Borcherds half-boundary into a
Fano object.  Adding the full mirror orbit would encode the group, but an
arithmetic hyperbolic reflection group has infinitely many mirrors; that is
not a finite-type boundary.

## 7. What survives and what is still missing

The valid positive consequences are conditional.

> **Proposition 7.1 (scoped log-geometric bound).**  Suppose a class of exact
> carriers has dimension \(n\leq N\), facet coefficients at least
> \(\delta>0\), quadric coefficient \(b\geq0\), and is uniformly
> \(\epsilon\)-lc log Fano.  Then (2.2) bounds \(f/(n+1)\), and BAB puts the
> underlying log pairs in finitely many bounded families.  If, in addition,
> the scheme-theoretic exact-Gram image on those families is
> zero-dimensional, or the Gram is returned by bounded integral
> cohomological data, then only finitely many normalized Grams occur.

The first conclusion is elementary positivity plus BAB.  The last conclusion
uses that a zero-dimensional finite-type image has finitely many geometric
points, or that bounded integral cohomological data are locally constant
after finite stratification.  Theorem 3.1 proves that the last hypothesis
does not follow from K-stability, boundedness, exactness, or absolute
reflection maximality alone.

For the global maximal arithmetic class, the missing inputs are exactly:

1. a uniform dimension/facet/flat/stabilizer-complexity bound;
2. a finite collection of coefficient and volume/Hilbert species covering
   both the Borcherds and (X3) regimes; and
3. a Diophantine theorem converting Vinberg integrality, conjugate
   definiteness, and reflection maximality into a zero-dimensional integral
   Gram locus or uniform field-degree and height bounds.

Log Chern, BMY, BAB, K-stability, wonderful compactification, and root stacks
do not establish any of these three implications.  They remain useful after
MAX--ARR is proved, but invoking them before then assumes the needed
boundedness in their hypotheses.

## 8. Primary inputs and internal controls

Primary sources:

* C. Birkar, [*Singularities of linear systems and boundedness of Fano
  varieties*](https://arxiv.org/abs/1609.05543).
* M. de Borbon and D. Panov, [*A Miyaoka--Yau inequality for hyperplane
  arrangements in \(\mathbf {CP}^n\)*](https://arxiv.org/abs/2411.09573),
  v3 (2026), especially Definition 6.1, Theorem 6.12, and Example 6.20.
* K. Fujita, [*K-stability of log Fano hyperplane
  arrangements*](https://arxiv.org/abs/1709.08213), especially Theorem 1.5
  and Proposition 1.9.
* A. Langer, [*Logarithmic orbifold Euler numbers of surfaces with
  applications*](https://arxiv.org/abs/math/0012180).
* C. Cadman, [*Using stacks to impose tangency conditions on
  curves*](https://arxiv.org/abs/math/0312349).
* C. De Concini and C. Procesi, [*Wonderful models of subspace
  arrangements*](https://doi.org/10.1007/BF01589496).

Internal exact controls:

* [universal exact quadric-arrangement theorem](universal_quadric_arrangement_carrier_and_log_boundedness_no_go.md);
* [maximal Lambert stable-pair theorem](maximal_lambert_exact_stable_pair_and_vhs_obstruction.md);
* [Vinberg integral-point and Farey audit](vinberg_integral_points_on_quadric_arrangement_moduli.md);
* [Borcherds \(\mathbf H^{21}/210\) absolute-maximality theorem](borcherds_h21_outer_routing_counterexample.md); and
* [Vinberg \(X3\) exact carrier](vinberg_x3_h19_explicit_carriers.md).
