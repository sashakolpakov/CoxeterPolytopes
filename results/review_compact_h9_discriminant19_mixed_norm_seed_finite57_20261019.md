# Independent written review of the mixed-norm discriminant19 seed

**PASS, within the stated archival dependencies.** No mathematical correction
is requested to
[the construction report](compact_h9_discriminant19_mixed_norm_seed_20261019.md),
SHA256 `bc47f8c2253ac61308c38c497cecaa7ef161f9f825ad4c3e2d5ffa1249c63b13`.
This verdict binds that original20-wall version only. A subsequently
announced21st-wall/Bugaenko-facet addendum is outside this review and
remains unreviewed here, even if later saved at the same report path.
The report supplies a changed integral20-wall seed, an obstruction to its
entire finite binary-only completion family, and a universal lower height
for long mirrors with E8 contact. It claims neither a finite-volume example
nor a complete mixed-norm census. Those scope restrictions are necessary.

This review used source reading and paper algebra. No matrix-calculation
program, scientific checker, search, Lean or build was run. Small JSON
coordinate inspection and streamed source hashing are metadata work, not
a replay of the archived Vinberg search or ray certificate.

## Primitive reflection norm alphabet

The form is even and integral over O=Z[phi], with determinant-phi*(4+phi).
The element b=4+phi is totally positive, with norm19. If a primitive
lattice-preserving reflection vector has square2d, reflection integrality
and a Bezout combination of its coordinates give Gr in dO^10. Multiplying
by the adjugate and using Bezout again gives d dividing detG. The determinant
ideal is prime, so d is a unit or a unit times b. Positivity at both real
embeddings makes that unit totally positive. The existing norm-Euclidean
and unit argument for O applies unchanged: totally positive units are
phi^(2k), and a unit rescaling of the vector leaves canonical norm2 or2b.

This is complete for reflection directions preserving this particular
free lattice. It is not an enumeration of either class, and does not
restrict an arbitrary geometric completion that need not preserve it.
Both classes occur: the archived simple roots have norm2 and the new
r=(1,2,0^8) is primitive of norm2b.

## Literal wall replacement and actual facets

For r one obtains Gr=(0,b,0^8), so its reflection is exactly
S(x,y,z)=(x-y,-y,z), an integral involution. Its height is
D/(4b)=phi/4. The old literal coordinates verify

    old10=old0+r,
    old20=old12+phi*r,
    old21=old13+phi*r.

For example old12 has binary coordinates(1+phi,-phi), while old20 has
(1+2phi,phi), with the same E8 coordinates. Thus deleting these three
walls and adjoining r gives exactly oldcone intersect{y<=0}; the dropped
inequalities follow from positive combinations of retained inequalities.
The deleted walls have positive product with r, so an append-only claim
would have been invalid. The report correctly makes the replacement.

The retained y values are0,-1,-phi,-phi^2. Their products with r are by;
the nonzero normalized magnitudes are at least sqrt(b)/2>1. Hence each
new nonright pair is strictly ultraparallel. Other products retain their
archived Coxeter status. The archived strict timelike center remains
strict for r because its y=-phi+1/10 is negative.

The projection argument proves that every one of the20 inequalities is
an actual facet. For a wall m, project the strict center c to
x=c-(c·m)m/m^2. Then x^2=c^2-(c·m)^2/m^2<0. For every other wall n,
x·n=c·n-(c·m)(m·n)/m^2<0 because both c·m and m·n are nonpositive.
Thus the projected point is timelike, on m, and strict on every other
inequality. There is no hidden redundancy inferred from the wall count.

The initial p is still incident exactly to E8+A1 and strict on r and all
positive-height walls. The proposed v_new=(1,0,0^8) has square-2phi and
is incident exactly to E8+r. At v_new the old A1 gives-phi<0 and every
retained positive-height wall gives phi*(y_i-2x_i)<0. These are actual
finite vertices. The report correctly declines to preserve the whole
old vertex ledger across the cut.

## Whole binary-only extension family

The archived spacelike ray v has binary part p/2. Its y is negative, so
it survives r, and all its nine archived independent active walls remain
retained. The archived ray therefore remains extreme. Extremality is not
needed for the stronger completion obstruction.

Every wall w orthogonal to E8 satisfies v·w=(p·w)/2. Retaining p necessarily
retains v. Consider any finite full-dimensional extension with this property.
It has a timelike strict interior point u. For sufficiently small epsilon>0,
(1-epsilon)v+epsilon*u is strict for all its finitely many walls and remains
spacelike. The segment from u to that strict spacelike point crosses the
future light cone inside every wall halfspace. An open ideal-boundary patch
is therefore exposed, forcing infinite hyperbolic volume.

This proves the stated family obstruction without a norm, height, label
or integrality restriction on the added binary walls. Its finite-extension
and retained-p hypotheses are explicit. It is not a proof that no completion
using E8-contacting mirrors exists, and it does not prove nonreflectivity.

## Long-mirror congruences and universal E8 trace bound

For a canonical norm2b reflector, the E8 block's unimodularity and
Gr in bO^10 give z_actual=bz. The binary congruence is y=2x modulo b.
The second binary coordinate of Gr gives integral
L=(phi*x+2y)/b; reducing modulo2, where b=phi, gives L=x modulo2O.
Substitution into the exact diagonalized norm identity gives

    b*(L^2+2q(z))=4+phi*x^2,
    h=phi*x^2/4.

For x>0 the conjugate equation gives
0<=q(z)'<2/b'<2(phi')^2. The strictness is valid because x is nonzero
and therefore its conjugate is nonzero. Also b'=5-phi>phi^2.

Writing z=a+phi*d with ordinary E8 vectors gives
T=Trq=2A+2J+3E=2||a+d/2||^2+(5/2)E. This is positive even when z!=0.
Under T<=12, evenness and positivity restrict E to0,2,4. The report's
subcases are exhaustive:

* E=0 gives q'>=2.
* E=4 gives J<=-A and q'>=4(phi')^2.
* E=2 gives A+J<=3. For A=0, q'=2(phi')^2. For A=2, J<=1 gives
  q'>=2+2phi'+2(phi')^2=4(phi')^2. For A>=4, J<0 gives q'>=A>=4.

Every case contradicts the strict conjugate upper bound, so T>=14.
Taking the trace of the norm identity yields
T<18/19+2h/b: the conjugate phi' term is strictly negative and the L^2
trace is nonnegative. Thus h>124b/19. This is an all-vector argument,
not a bounded shell census. At x=0, reflection congruences force the
whole vector into bLattice; its divided norm2/b is not an even integral
lattice norm, so no long height-zero case escapes the argument.

For theta=(2,4,6,5,4,3,2,3), the native E8 graph gives theta^2=2: its
coefficient-square sum is119 and its edge-product sum is118. With
old9=(1,0,phi*theta), reflection in old0 sends r to(1,-2-phi,0).
That vector has product-D with old9. Applying the latter reflection
therefore gives exactly
w_*=(2+5phi,-2-phi,(5+6phi)*theta). Its norm2b and height(45+74phi)/4
follow, but its product with old9 is+D. The report correctly records
this as an inadmissible reflected mirror, not as a successful cap.

## Archival boundary

All three report source pins match the current files:

* `data/candidates/h19_vinberg_prefix.json`:
  `91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1`.
* `src/coxeter_polytopes/search/h19_vinberg.py`:
  `991eb7a6078e5880c8d33778a2f425e5d3d8d29023f61bec772ebcc36d46be38`.
* `results/coupled_discriminant19_candidate.md`:
  `5242d0b8333a357164ae671525882762a58834bebb5430d1a6447f67a66d3b16`.

I read the literal archived wall coordinates, the norm-two search's declared
scope and algebra, the E8 Cartan definition, the existing independent
prime-norm alphabet proof, and the ray-diagnostic source. The historical
prefix enumeration, all old pair tests, the strict-center norm and the
ray's full coordinate/norm/active-rank certificate are inherited rather
than freshly replayed. Hash agreement identifies those inputs; it does
not itself prove their mathematics. No complete mixed-norm first-height
search, final vertex/cusp census, finite-volume chamber or K3 nef marking
is supplied by this review.
