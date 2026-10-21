# Independent hostile audit: the anticanonical genus-one facet bound

**Verdict (2026-08-11): PASS under the literal anticanonical-facet
hypothesis.**  The intrinsic bound is characteristic-free and uses no
Nikulin--Vinberg or chamber classification.  Both endpoints are realized:
(P_{10}) has (F^2=-1), while (P_{11}) has (F^2=-2) and attains the
sharp projective-nef dimension (11).

## 1. Why the anticanonical class is primitive

Let (S) be smooth projective rational over an algebraically closed field,
with (\rho(S)\ge3).  Minimal rational surfaces have Picard rank one or two,
so (S) is nonminimal and has a ((-1))-curve (E).  Adjunction gives

\[
                              K_S.E=-1.
\]

If (K_S=mD) in (\operatorname{NS}(S)), then (mmid K_S.E=-1).
Thus (K_S), and hence the literal class (F=-K_S), is primitive.  This
is exactly where (\rho\ge3) is used.

## 2. Unimodularity and the reflection norm

The Neron--Severi lattice of a rational surface over an algebraically closed
field is unimodular.  A primitive vector in a unimodular lattice has
divisibility one, so some (x\in\operatorname{NS}(S)) satisfies (x.F=1).
If

\[
                  s_F(x)=x-{2(x.F)\over F^2}F
\]

is integral, primitivity of (F) forces (2/F^2\in\mathbf Z).  A genuine
hyperbolic facet has negative normal square: a positive normal has
negative-definite orthogonal complement and misses hyperbolic space, while
an isotropic normal is tangent only at infinity.  Therefore

\[
                              F^2\in\{-1,-2\}.
\]

For a rational surface, (K_S^2+\rho(S)=10).  Since (F^2=K_S^2),

\[
 \rho(S)\in\{11,12\},\qquad
 \dim\mathbf P(\operatorname{Nef}S)=\rho(S)-1\in\{10,11\}.
\]

Adjunction for (F=-K_S) gives (p_a(F)=1), whether or not the irreducible
curve is geometrically smooth.

## 3. Exact (P_{10}) and (P_{11}) endpoints

I independently regenerated the walls in (I_{1,n}):

\[
 E_i-E_{i+1},\quad E_n,\quad H-E_1-E_2-E_3,
 \quad F=3H-E_1-\cdots-E_n.
\]

The walls other than (F) form a unimodular basis.  Exact rational
active-set enumeration checks all (66) sets for (n=10) and all (78)
sets for (n=11).  It recovers

\[
\begin{array}{c|ccc|c}
 &\text{finite}&\text{ideal}&\text{hyperideal}&F^2\\ \hline
P_{10}&17&2&0&-1\\
P_{11}&26&2&0&-2.
\end{array}
\]

Every displayed inequality is a genuine facet.  Direct basis tests show
that reflection in (F) is integral in both lattices, and the pinned
pencil-at-infinity source supplies the complete-nef equality for the two
flex-cubic carriers.  Their strict cubic is smooth genus one.  Hence both
allowed norms occur and (P_{11}), not (P_{10}), is the sharp endpoint.

The theorem requires the literal primitive curve class (F=-K_S).  It does
not apply merely because a nonprimitive curve spans the same hyperplane,
does not identify the twice-primitive rational ((-4)) decoration with this
genus-one branch, and does not classify all chambers in dimensions ten or
eleven.

The independent executable audit is
[audit_rational_pfirst_anticanonical_genus_one_facet_dimension_bound_independent_hostile.py](../scripts/audit_rational_pfirst_anticanonical_genus_one_facet_dimension_bound_independent_hostile.py),
with separate tests in
[test_rational_pfirst_anticanonical_genus_one_facet_dimension_bound_independent_hostile.py](../tests/test_rational_pfirst_anticanonical_genus_one_facet_dimension_bound_independent_hostile.py).
