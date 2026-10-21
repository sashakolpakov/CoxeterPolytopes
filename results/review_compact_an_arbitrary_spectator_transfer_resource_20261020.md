# Outside review: A21--A29 components with arbitrary spectators

20 October 2026. Full paper-only source audit by resource_cleanup.
No scientific execution, imports, optimization, Lean or TeX; no accepted
source or status edits.

**Verdict: PASS.** No mathematical gap or requested correction was
found in
`results/compact_an_large_components_arbitrary_spectator_transfer_20261020.md`,
SHA256 `8d3b999761be7dd46e7990973b3ec4d4a746d4d682b044d53695f4aae405e3e5`.
The conclusion excludes a WHOLE A_r COMPONENT, 21<=r<=29, at a
farthest vertex of a compact Coxeter polytope, with arbitrary finite
orthogonal spectator components. It is not a new ambient dimension
bound or an unbounded-rank local result.

I had previously reviewed the whole-A plateau transfer and contributed
to checking some underlying actual-projection interfaces. This review
does not claim independent discovery of those imported arguments. The
new all-spectator clearing, zero-own-contact repair and complete weak-
profile transfer were independently checked here from the displayed
formulas and actual-wall hypotheses.

## Spectator clearing and isolated caps

At the actual endpoint all spectator normals plus c form a positive
definite subdiagram. Hence Q<1, and its positive-inverse correction
w_S is a nonnegative combination of actual spectator normals. The
single-contact-per-component rule and the cost floor1/4 follow from
finite trees and the unit inverse-Gram diagonal bound. They require
no rank or type restriction on the spectator forest.

For isolated component contact gamma, I rederived both clearing
formulas. With d0=1-Q and u=d0/gamma^2, the positive residual
d0-gamma^2 gives the finite-face pressure on
`gamma*(omega_i-n_i)` after every spectator term cancels. For u>1,
putting B=sqrt(D_i-1), t=sqrt(u-1) gives the stated bound
`h_i>sqrt(B^2-t^2)/(1+t/B)>B-t`. The last squared difference is
exactly `(B-t)*t^2/(B+t)>0`.

For 0<u<1 the separate vector
`xi=c+w_S+(d0/gamma)n_i` has the three asserted products, square
`-d0*(1-u)` and positive time coefficient A. Its coefficients are
positive combinations of ACTUAL normals, with the selected c, n_i
and spectator inequalities checked separately. Thus it respects
every other wall of the full polytope. Farthest comparison gives
`h_i-u*(x_A.n_i)>=sqrt(D_i-u)-sqrt(u*(1-u))`, and strict old
interiority supplies the displayed strict bound. At u=1 it is a
nonzero feasible future null ray, incompatible with compactness.
There is no promotion of an orthogonalized wall to an actual facet.

For Q>0, gamma>=phi/2 gives u<=6-3phi. The two deficit estimates
follow respectively from q+B>21/4 and from maximizing
`(4/21)u+sqrt(u*(1-u))` on [0,1]. Their strict upper bounds
17/28 and127/210 are below61/100. The affine principal exclusions
for isolated3,4,6 and parallel contact survive spectators. With Q=0
the imported isolated pressure proofs apply in the original p+A
span. Therefore no isolated central cap is weak.

## Exhaustive zero-own and foreign-contact repairs

The two retained paths have rank at least four. The finite one-node
extension list is exactly the one imported: a single3 at an allowed
A/D/E position, an endpoint4, or the two endpoint3 contacts. A
two-path connection has rank at least21, excluding exceptional
finite diagrams; its branch-arm lengths exclude D, and an internal4
excludes B. Spectators cannot enlarge this restriction of an actual
positive-definite Gram matrix.

With own contact zero, two white3 contacts contain the all3 affine
cycle. A single internal3 or4 contains the stated affine principal
star. The remaining far-end and next-to-end cases are all covered:

* End3 has full component cost `(N-1)/(2N)` and retained cost
  `(m-1)/(2m)`. Its positive height requires a spectator, but is
  less than1/(2m)<=1/10. Every contacted spectator has squared
  derivative at least1/4, so actual shadow recognition forces a
  second-component closure.
* Node2/simple3 has Q<2/m<=2/5, so there is at most ONE contacted
  spectator component. Its full cost Q exceeds `A^2=Q-2/N`, and
  its squared derivative is at least Q. This establishes the
  required actual second-component endpoint, not merely a small
  aggregate height. Without that component the height is negative.
* End4 has Q<1/m<=1/5, forcing Q=0 and negative height.

These closures contradict the accepted farthest orthogonal-component
theorem. The argument is symmetric at the other global end and does
not use the whole-vertex signature shortcut that would fail in the
presence of spectators.

For positive own contact, adding Q only decreases the actual white
Schur residual in the edge formula. The accepted nonexceptional
edge-only strong bounds therefore persist. The two exceptional
foreign-contact estimates were checked separately: the far-end3
ratio `(2fm+1)/sqrt(m*(m+2))` increases for m>=5, with derivative
sign `(2f-1)m-1>0`, and its endpoint squared comparison has margin
342365. In the near4/simple3 role, OWN4 and WHITE3 give
`sigma=1/sqrt2+(m-1)/(2m)>11/10` and residual at most7/20 after
Q>=1/4. Its squared ratio is greater than121/35. Both exceed the
uniform F(q)/q upper761/500. Hence EVERY weak central cap has
Q=0, with exactly the same literal profile and height as before.

## Complete certificate transfer

The required immutable inputs were hash-checked:

* whole A21 certificate: `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb`;
* A21--29 plateau transfer: `8a58b7747e1a3d8c9affe7df17038781a62983ddd70ccb670502c6f40b1f9cee`;
* finite-face pressure: `abfcb96b6eb74aefdc6498cda362d81b1d4a1f9c32ad17fe3b464cff26fba565`;
* actual mate/packing interface: `8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b`.

The new repair supplies precisely what the old proof lacked for
arbitrary spectators. Every weak near4 or far-helper cap has zero
spectator projection. Its mate recognition, height, mutual products
and native-gap calculations are therefore unchanged. All selected
finite-face and timelike pressure constructions use these actual
caps and old A normals with the same positive coefficients. They
pair zero with old spectators and nonpositively with every other
actual wall. For selected full-rank finite replacements, adjoining
the spectator block restores the full vertex rank by an orthogonal
positive-definite direct sum. Actual point comparison consequently
uses only x_A, even though the original interior point has arbitrary
spectator coordinates.

The old A concavity, norm identity, nonadjacency, positive-dual
weights, isometric plateau lift and six-hole losses therefore apply
to x_A itself. The inherited strict lower bound4013/4000 contradicts
`||x_A||^2<=||x||^2<1`. No certificate step was replaced by scalar
edge floors alone. No intrinsic Coxeter property or full reflection
arrangement of a metric face is assumed.

The source's rank21--29 theorem is established relative to its
explicit accepted classical and actual-geometry inputs. Its broader
statement about all relevant large A ranks also explicitly imports
the existing universal compact ceiling29; no new ceiling follows
from this review.
