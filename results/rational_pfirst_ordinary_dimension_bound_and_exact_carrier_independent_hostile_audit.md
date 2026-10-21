# Independent hostile audit: ordinary P-first cutoff and exact carrier

**Verdict (2026-08-11): PASS after making canonical escape automatic in
the owner theorem.**  This audit was written independently of the owner
certifier.  It checks the live refrozen bytes, reconstructs the standard
simplex and its dual rays over the integers, and separately checks the
timelike, cusp, and flex-carrier compositions.

## 1. The wall-cone argument has exactly the needed scope

Let

\[
 \mathcal W=\{x:x.r_j\geq0\}
\]

for the finite list of facet normals, and take a strict future point
\(A\).  If \(w\in\mathcal W\) is spacelike, then \(A+tw\) remains strict
for every wall while

\[
 (A+tw)^2=A^2+2t(A.w)+t^2w^2
\]

has one positive zero.  The zero is nonzero and future-null.  Because the
wall list is finite and every slack is still strict, a positive-measure
ideal patch and a radial collar lie in the chamber.  Their ball-model
volume diverges.  Thus a finite-volume finite-sided chamber cannot have a
spacelike vector in its homogeneous wall cone.

Finite-sidedness is essential at precisely the ideal-patch step.  It is
present in Lemma 1.1 and Theorem 2.1; throughout the combined theorem,
“Coxeter chamber” must retain its stated finite-facet convention.

For an ordinary wall, \(F=-K\) has slack \(1\) on a \((-1)\) facet and
\(0\) on a \((-2)\) facet.  Hence \(F\in\mathcal W\), while

\[
                             F^2=9-n.
\]

The independent replay constructs a strict point with every wall slack
equal to one for each \(10\leq n\leq15\).  Along its ray toward \(F\),
the light-cone polynomial has negative root product and the positive
crossing retains all strict slacks.  This is an exact sample of the
general symbolic argument, not a census input.  The cutoff \(n\leq9\) uses
neither a surface nor reflection maximality.

## 2. Canonical escape is a conclusion, not a hypothesis

There was one mismatch in the first frozen presentation: the opening
claimed the result from ordinary finite volume and relative maximality,
while Theorem 3.1 additionally assumed \(K\notin P\).  The extra assumption
is unnecessary.

After the spacelike cutoff, \(F\) is a nonzero causal vector in
\(\mathcal W\).  The second part of the wall-cone lemma excludes the past
sheet: a past causal vector not proportional to \(-A\) produces a
spacelike positive combination, and proportionality contradicts the wall
slacks.  Therefore \(F\) is future-directed and \(K=-F\) is past-directed,
so automatically

\[
                              K\notin P.
\]

The owner theorem now states this as equation (3.0) and no longer assumes
canonical escape.  Its opening and theorem hypotheses therefore agree.

## 3. The timelike orientation order is non-circular

The physical signs alone do not choose a positive sheet.  The valid order
is:

1. maximality inside the declared full host \(O^+(I_{1,n})\) puts every
   integral reflection into the side group, because adjoining one remains
   a discrete reflection subgroup of the arithmetic host;
2. signed permutations and the Cremona reflection reduce every primitive
   reflecting vector of square \(-1\) or \(-2\) to a standard side;
3. an oriented reflection isometry carries the chamber to the standard
   simplex;
4. only then is the image \(G\) of \(F\) placed in the homogeneous
   standard inequality cone; and
5. the reduced-characteristic identity forces \(G=F\), so the isometry
   fixes \(K\).

For \(3\leq n\leq8\), writing

\[
 k_1=k+A,\qquad k_2=k+B,\qquad k_3=k
\]

gives the exact gap

\[
 (9-n)k(k+1)+4k(A+B)+2AB+2A+2B.
\]

Every term is nonnegative and the first coefficient is positive.  The
gap vanishes only at \(k=A=B=0\), after which ordering and the norm equation
force every remaining coordinate and \(\ell\) to vanish.  Thus \(F\) is
the unique reduced characteristic vector of square \(9-n\).  No initial
nefness or sheet choice occurs in this step.

The host qualification matters: this reasoning uses maximality among all
reflection overgroups inside the declared full \(O^+\), not merely
maximality in a smaller hand-picked family.

## 4. The rank-ten cusp premise is proved internally

At \(n=9\), \(F\) is primitive, characteristic, null, and future.  If its
ray were not extreme in \(\mathcal W\), it would split into two
nonproportional vectors of \(\mathcal W\).  Spacelike exclusion and time
orientation make both summands future causal, while reverse Lorentz
Cauchy--Schwarz says a null sum is possible only for proportional null
summands.  Therefore \(F\) is an actual ideal vertex; the characteristic-
cusp theorem receives its cusp premise rather than assuming it.

Independent dualization of the standard ten-facet simplex gives dual
squares

\[
                  0,0,2,6,5,4,3,2,1,1.
\]

The two primitive null rays are

\[
 F=(3,-1,\ldots,-1),\qquad Q=H-E_1.
\]

Only \(F\) is characteristic.  Hence a chamber conjugacy cannot exchange
the two cusps and must fix \(F\), exactly as required by the rank-ten
composition.

The nine square-\(-2\) standard sides satisfy

\[
 F=2r_1+4r_2+6r_3+5r_4+4r_5+3r_6+2r_7+r_8+3r_9.
\]

Their incidence graph is a tree on nine vertices with eight edges, so its
reduced fibre has Euler number \(18-8=10\).

## 5. The flex carrier is an exact geometric certificate

At the common flex the pencil ideal is

\[
                 (z^3,\ z-x^3-xz^2).
\]

Put \(w=z-x^3-xz^2=z(1-xz)-x^3\).  Since
\(\partial w/\partial z(0)=1\), \(w\) is a regular parameter.  On
\(w=0\), \(z\) has \(x\)-adic order three, so \(z^3\) has order nine.
Thus the ideal is analytically \((w,x^9)\), proving colength nine.  Each
successive point blow-up changes \((w,x^e)\) to
\((w_1,x^{e-1})\).  After \(n\) blow-ups the residual ideal is exactly

\[
                         (w_n,x^{9-n}).
\]

The base branch has contact order three with the flex line.  Consequently
the first three centers lie on its successive strict transforms and give
the effective class \(H-E_1-E_2-E_3\); the exceptional chain gives all
classes \(E_i-E_{i+1}\), with \(E_n\) the last exceptional.  For
\(n\leq8\) the residual base scheme is zero-dimensional, hence the pencil
has no fixed curve and its class \(-K\) is nef by intersection with a
member not containing a chosen irreducible curve.

For \(n=9\), the discriminant \(4+27\lambda^2\) has two simple roots, so
the only reducible fibre is the displayed affine nine-component fibre.
The determinant of the lattice generated by \(F\), \(E_9\), and the eight
components not meeting \(E_9\) is \(-1\).  Hence the section quotient is
trivial and \(E_9\) is the unique section.  Adjunction then exhausts the
negative curves: a negative horizontal curve must be a section, and a
negative vertical curve must be one of the displayed fibre components.

For \(3\leq n\leq8\), nefness of \(F\) similarly restricts any negative
irreducible curve to square \(-1\) or \(-2\).  Its integral reflection is
in the full side group, so its mirror cannot cross the chamber interior;
an ample class selects the nonnegative side.  Effective curves of
nonnegative square are handled by the closed future-cone Lorentz
inequality.  This proves both inclusions in the exact-nef return.

## 6. Scope and replay

The audited theorem is exactly for finite-sided integral Coxeter chambers
in the standard marked \(I_{1,n}\), \(n\geq3\), with only ordinary
\((-1,-1)\) and \((-2,0)\) canonical colors.  Relative reflection
maximality is needed for the carrier classification, but not for the raw
dimension bound.  Nothing here covers twice-primitive \((-4)\) facets,
nonstandard characteristic orbits, or arbitrary adapted hosts.

The independent executable audit is
[audit_rational_pfirst_ordinary_dimension_bound_and_exact_carrier_independent_hostile.py](../scripts/audit_rational_pfirst_ordinary_dimension_bound_and_exact_carrier_independent_hostile.py),
with separate tests in
[test_rational_pfirst_ordinary_dimension_bound_and_exact_carrier_independent_hostile.py](../tests/test_rational_pfirst_ordinary_dimension_bound_and_exact_carrier_independent_hostile.py).
