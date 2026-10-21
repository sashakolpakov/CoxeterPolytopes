# Rational parameters cannot be rescued by finite Coxeter subdivision

For rational values of the affine parameter, the sections in the [59 excluded fork-leaf normal planes](fresh_h20_cusp_sections_20261013.md) cannot admit a finite face-to-face subdivision into full-dimensional Coxeter polytopes. The same statement holds for the interior rational parameters of the [pure A1 plane](fresh_h20_cusp_sections_20261013_absorption.md). This is a consequence of their certified actual ridge angles; it does not address irrational parameters or different section planes.

The elementary angle test is the following. If `θ/π` is rational, then `exp(2iθ)` is a root of unity. Thus `2cos(2θ)` is an algebraic integer. If `cos²θ` is also rational, that algebraic integer is rational and therefore an integer between −2 and 2. Consequently

\[
\cos^2\theta\in\{0,1/4,1/2,3/4,1\}.
\]

At every interior rational parameter of each fork-leaf plane, the displayed actual ridge has rational squared cosine strictly between 0 and 1/4. At the central parameter it is 1/7. The ridge angle is therefore an irrational multiple of π.

For the pure A1 plane, let `u=t²`, with `0≤u<1/4`. Its actual obtuse ridge has

\[
\cos^2\theta=\frac{(1/4-u)^2}{49/16-(9/2)u+u^2}.
\]

This is positive and at most 1/49: the denominator minus 49 times the numerator is `4u(5−12u)≥0`. For rational t it is again rational and outside the displayed five-value list.

Finally, take a generic point in the relative interior of the certified ridge. In a finite face-to-face subdivision, the finitely many tile faces of codimension at least three cannot cover that ridge. Every tile contributing an angular sector at the chosen point therefore contributes one of its dihedral angles, `π/m`. Those finitely many sectors partition the original wedge, making its angle a finite sum of rational multiples of π. The irrational ridge angle contradicts this necessary condition.

This rules out finite Coxeter subdivisions of the stated rational sections. It is not a general obstruction to gluing other chambers, changing the section normal, or using an irrational parameter.
