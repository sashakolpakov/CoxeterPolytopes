# Physical exclusion of the low-S/low-R continuous high-N family

On the two fixed parent chains, no wall system in the declared low-S/low-R
domain can have a finite high N/X bond and the required finite disk corners.
This upgrades the [numerical event search](fresh_oblique_low_s_continuous_n_events_20261018.md)
to directed physical interval exclusion. Every S/R norm sheet is covered;
no numerical native-wall rejection is used as a proof. It does not exclude
high R/M or high R/S, other parent chains or other incidences, and supplies
no new finite-volume example.

Files: [producer](../scripts/fresh_oblique_low_s_continuous_n_interval_20261018.py),
[certificate](fresh_oblique_low_s_continuous_n_interval_20261018.json),
[standalone checker](../scripts/check_oblique_low_s_continuous_n_interval_saved_20261018.py).

## Domain and physical inputs

The two complete-instruction parents retain all four original selected
origins. Fixed new walls are J67,L68,X69,M71. Remove ancestral R70 before
any new R construction. The original deletion set remains30/44/63/64/65;
old63/64 are recipe witnesses, never retained inequalities.

N has profile0 and frame10/J/X labels(2,2,m) or(3,2,m), for every finite
m>=6. Write n=2cos(pi/m); the proof uses the larger closed interval
[sqrt3,2]. A finite X/N/R corner forces RX=RN=2: a norm2 high bond has
squared magnitude at least3, and any other nonright native bond already
prevents positive definiteness. The seven possible R white profiles are
retained. S has any of the seven profiles and all three frame52/L/M labels
in2–5. Its actual duties are52/L/S and L/M/S, without a52/M corner.
R has frame X/M/S labels(2,r,s), r,s in2–5, with actual duties X/M/R
and M/R/S, without an X/S corner.

The exact inherited finite-duty classifier builds labeled diagrams and
matches their connected components against finite spherical diagrams;
it is not a floating eigenvalue filter. This gives42 S recipes on each
parent and1,238 R recipes over the76 surviving S sheets. The standalone
checker separately reproduces these lists with the complete connected
rank-at-most-five catalogue: A/B, D4/D5, F4, H3/H4 and I2. The four fixed
H4 whites form a separate positive component, leaving the five-node
A2-plus-three-wall test. The seven-profile premise itself remains inherited.

The [original63-wall certificate](fresh_nonarithmetic_release49_62_replacement_20261017.json)
provides physical facet witnesses and the exact Coxeter scaffold. Its
normalization makes every original wall norm2. The
[prior physical projection](fresh_oblique_high_s_polar_interval_20261018.md)
encloses the Gram of the six fixed whites and old8/10/52/63/64. Using
facet witnesses directly introduces interval dependency widths that grow
under repeated extensions. The new producer sharpens only identities
already forced by that original Coxeter scaffold:

* norm2 diagonal entries;
* an original product whose enclosure isolates zero;
* an original finite product whose enclosure isolates exactly one native
  label, with strict rational comparisons against the adjacent labels.

There are60 such upper-triangular identities, involving norm2 and labels
2,3,5,27. In particular old8 is orthogonal to10 and52, old52 is orthogonal
to64, and old63/64 has label27. Each new enclosure is intersected with the
original enclosure. Dotted products are unchanged. Floats may select a
candidate integer label, but the exact isolation inequalities certify it;
no new S or R product is rounded to an allowed angle.

The six-white inverse and projected Gram are recomputed. All five parent
extensions are then reconstructed using their frozen frames, profiles,
labels, witness roots and witness signs. This preserves the exact parent
definitions, while reducing the interval dependency widths. A projected
parent contact needed as a right finite-duty edge is also enclosed strictly
between the first nonzero native magnitudes; any admissible Coxeter parent
must therefore have that contact zero. This is a necessary condition,
not an assumption that every fixed parent gate has independently been proved.

## Complete physical sheets and native rejections

For each three-frame, its determinant is strictly negative and its first
diagonal is positive. It therefore has Lorentz signature(2,1). The residual
of the selected witness has strictly positive square, and the requested
norm remainder is strictly positive. The two roots are exactly the two
signs of this positive one-dimensional residual line.

Both signs relative to the same witness are checked for every recipe.
Numerical products name which saved branch corresponds to each sign; the
proof includes both signs regardless of those names. This avoids depending
on the old absolute comparison tolerance for rejected numerical roots.
All native products restore the full six-white profile offset.

A rejection uses just one retained wall and one strict interval condition:
positive product, a product strictly between−1 and0, or a product strictly
between two consecutive negative native magnitudes. These are necessary
Coxeter gates. The92 S rejections consist of84 positive products and eight
native gaps. The2,416 R rejections consist of1,472 positive products and
944 native gaps; their witnesses are old8, old10 or J67.

The remaining76 S sheets and60 R sheets are an overdomain. The proof does
not assume they pass every retained native or center gate: retaining extra
systems makes the subsequent exclusion stronger. No arbitrary supporting
wall or extra inequality is silently added to the frames.

## Continuous N event exclusion

The frozen [continuous-N proof](fresh_oblique_continuous_high_n_high_s_interval_20261018.md)
excludes the other N sheet by a positive old8 product on all66 covering
intervals. Its remaining sheet has

\[
N(n)=u+n v+\sqrt{\rho_0+\rho_1n+\rho_2n^2}\,b.
\]

For each remaining R, the right contact RN=0 implies

\[
(A+Bn)^2-C^2(\rho_0+\rho_1n+\rho_2n^2)=0.
\]

Every leading coefficient and discriminant has a strictly positive lower
bound. Directed quadratic formulas enclose both roots of all120 equations.
Thus constant polynomials, tangencies and missing norm sheets do not occur.
Squaring gives an overdomain: exclusion of every squared root suffices,
so no unsquared existence claim is needed.

Of240 roots,236 lie strictly outside [sqrt3,2]. The remaining four, two on
each reflected parent, lie strictly in native-label gaps. Rounded outward
display intervals are:

| N/10 label | Enclosed positive root | Native gap |
|---|---|---|
|2|[1.9127046425,1.9135468897]|2cos(pi/10) < n < 2cos(pi/11)|
|3|[1.8764703175,1.8775602451]|2cos(pi/8) < n < 2cos(pi/9)|

The minimum proved separation from either neighboring native magnitude
is greater than0.00182499. General cos(pi/m) bounds use directed rational
Machin arctangent bounds and a Taylor remainder, including m=9,11,27 and
the neighboring labels used for isolation. Monotonicity of2cos(pi/m)
for integer m>=2 excludes every native finite label in these gaps.

| Complete stage | Count |
|---|---:|
| S recipes / physical sheets |84 /168|
| Strict S exclusions / remaining overdomain |92 /76|
| R recipes / physical sheets |1,238 /2,476|
| Strict R exclusions / remaining overdomain |2,416 /60|
| Quadratics / enclosed roots |120 /240|
| Outside-domain roots / native-gap roots |236 /4|
| Unresolved cases / accepted systems / new cones |0 /0 /0|

## Scope and validation

The exact old debt remains19 lost finite supports,12 touched compact
two-faces and one compact three-face, including three lost supports on8.
The numerical predecessor independently reconstructs this from all199
original finite supports. The disk remains V/E/F=10/21/12 and boundary
8–10–N–R–S–52–8. Since every system is excluded before the cone stage,
no lost compact face is claimed restored and no full ambient H9 region
is certified. This is a complete necessary obstruction for the stated
low-S/low-R domain on the two exact parents, not an overall search result.

Final source-bound generation96091/ffff44 passed. The independent checker
uses standard-library rational interval operations, reconstructs the
original physical projection, checks every native-identity intersection,
replays the parent and S/R extensions, checks common-witness opposite sheets,
independently enumerates the exact finite-duty lists, and verifies every
strict native rejection and event-root gap. It does not reprove the original
scaffold or seven-profile classification and does not reconstruct cones.
Separate fresh producer replay51919/ad00e8 and final checker34481/ac7cbe
both passed. Root completed the full producer/report/helper mathematical
read, including the final source-binding and exact-catalogue/common-witness
deltas (1f4d54,0996d3,3bf3aa,782a2e,cd82f9,c087f3,dcaade). The read checked
original-scaffold native isolation, physical parent reconstruction, complete
finite-duty lists, both sheets with a common witness, every strict rejection,
quadratic-root overdomains, native gaps and the inherited opposite-N-sheet
exclusion. It did not execute the code or reprove the original scaffold,
seven-profile theorem or geometric cones. Earlier drafting attempts stopped
at conservative comparison and native-gap guards; no failed attempt wrote
this final certificate. The packet is final and frozen.
