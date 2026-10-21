# A universal short-path alternative for deep finite distance-eight pairs

Let an actual quadrilateral face of a finite-sided finite-volume real
hyperbolic Coxeter polytope have a finite corner whose two marked roots
lie in one whole finite A/B/D component at graph distance8. Suppose both
marks have distance at least17 from every endpoint, branch and endpoint
of a nonsimple bond of that whole component. Then the **full actual
supporting star contains a path of length at most7 between the two
opposite black pairs, with white internal roots**.

This statement covers all low/high, disjoint/overlapping minimal-core
cases, all native labels and arbitrary additional supporting roots.
It does not assume the old no-short8 condition. It is a conditional
short-path theorem, not a claim that every low or disjoint configuration
is impossible. Only in the explicit full-star no-short7 branch does it
exclude a deep finite-distance8 pair.

The [source](../scripts/fresh_universal_deep_finite8_corridor_20261018.py)
and [certificate](fresh_universal_deep_finite8_corridor_20261018.json)
check the exact rational inequalities below. No contact-profile or
low/high joining census is needed. A complete changed payment rule,
its finite-forest budget and any dimension consequence require their
separate composition. The unrestricted ceiling stays65 in this report.

## Actual-face inputs and the forced path

Write the star as W+a+b+c+d, with opposite black pairs {a,c} and {b,d}.
Every root of W is white. The actual white diagram and each W plus one
black are positive definite. Each W plus an opposite pair is indefinite;
these are the inherited strict-edge and closed-opposite-side facts of
the [complete face reduction](fresh_selective_distance9_dimension108_20261017.md).
The corner under consideration is W+a+b. Roots are normalized to
squared norm2, so every nonzero native acute contact has magnitude at
least1, including arbitrarily large finite labels, parallel or dotted
contacts. All distinct-root products are nonpositive.

Assume there is no full-star path of the stated length. In particular
c.d=c.b=d.a=0. The products c.a and d.b are unrestricted nonpositive
native products; they may vanish or have arbitrarily large magnitude.

The whole finite A/B/D component is a tree. Feature depth17 gives a
simple eight-edge path from a to b and two disjoint simple outward
arms of17 whites. A feature exactly at outward distance17 is allowed:
every bond inside the corresponding17-white arm is still simple.
Keep15 outward whites on each side. The retained39 roots form an
induced positive definite A39 path with positions0,...,38, and

    a at15,     b at23.

Its other37 roots are supporting whites. The two unused whites on
each outer arm remain available for actual-edge inverse bounds. This
is the same forced geometry as the frozen
[A39 high-domain theorem](fresh_full_a39_corridor_high_domain_20261017.md),
but no high-domain enumeration is imported here.

## Four possible contact positions on each side

For a norm-two A_n path the positive inverse entries are

    J_n[i,j]=(min(i,j)+1)(n-max(i,j))/(n+1),  0<=i,j<n.

All its leading principal determinants are k+1, so the path is positive
definite. The checker verifies the determinant recurrence and inverse
formula exactly at orders17 and39.

The omitted root c and any17-white outer arm belong to an actual edge
diagram. If c meets that arm at outward distance j=3,...,15, the
projection square is at least

    (A17 inverse)[j-1,j-1]=j(18-j)/18 >=5/2>2.

The native contact is at least1; all further contacts have nonnegative
magnitudes and the inverse is positive, so they cannot lower this bound.
It contradicts strict positive definiteness of that actual edge.
Thus c has no retained contact at those distances on either outer arm.

If c meets an inner white at distance j>=2 from a, following the inner
path to b gives a c-to-b path of length9-j<=7. Meeting either of the
first two whites beyond b gives a path of length2 or3. These contacts
are excluded by the full-star hypothesis. The argument for d is the
reflection. Therefore, inside the retained A39 path, the only possible
nonzero contact positions are

    c: 13,14,15,16;       d: 22,23,24,25.

No restriction on the magnitudes at these positions is needed. In
particular this proof does not impose a finite-label alphabet, rule
out simultaneous contacts, or bound the two own-mark products.
Contacts to roots outside A39 are allowed throughout.

## Why both retained contact vectors are nonzero

This is the step that covers the low-core and disjoint branches without
their old no-short8 hypothesis. Choose a minimal indefinite principal
diagram K inside W+a+c. Both a and c belong to K: deleting either
would leave a principal diagram of a positive-definite actual edge.
Minimal indefiniteness makes K connected.

Use the published all-label fact that a minimal indefinite Coxeter
diagram has order at most10. The precise classification input and
primary source are recorded in the frozen
[minimal-indefinite catalogue report](fresh_high_rank_quasilanner_census_20261017.md)
and its [local primary PDF](../data/references/fresh_finite_edge_cusp_tradeoff_20261017_chein.pdf).
The current proof uses only this order ceiling. It neither infers the
ceiling from a bounded search nor uses the26 high-order diagrams or
the finite low-order census.

Here is the precise bridge to that published statement. If K has
order2, a dotted pair already satisfies the ceiling. At order at
least3, no dotted bond is possible: its two-root indefinite principal
would contradict minimality. The remaining finite and parallel
products are exactly the Coxeter cosine-matrix alphabet, including
the infinite label, used by Chein. Because K with c deleted is strictly
positive definite, interlacing and indefiniteness give full signature
(|K|-1,1), with no zero eigenvalue. Every other deletion is positive
semidefinite and has finite or affine connected components. On printed
page4 (PDF page3), Chein explicitly allows such disconnected deletions
in his graph criterion, and states the order-ten bound for irreducible
hyperbolic matrices, attributing it to Vinberg. Thus the bound applies
to these actual minimal cores. No finite bound on order3 black-black
labels is assumed. The [primary paper on Numdam](https://www.numdam.org/item/M2AN_1969__3_3_3_0.pdf)
and the bound local copy contain the same printed-page statement.

Let B be the connected component containing a in K with c deleted.
It has at most9 roots. It is contained in W+a, hence also in the
actual finite corner W+a+b. Every root of B lies in the same whole
finite component as a, and is at distance at most8 from a. Feature
depth17 forces B to be a simple subpath of that component. It cannot
pass through b, since b does not belong to K. Consequently B fits
inside the retained A39: at most8 outward steps from a and at most7
steps toward b, while15 outward whites were retained.

Since K is connected, c has a nonzero contact with B. Otherwise B
would already be a separate connected component of K. Thus c has a
nonzero contact inside the retained A39. The same argument with a
minimal indefinite core in W+b+d proves the analogous assertion for d.
No hypothesis about whether the two cores intersect has entered.

Extra supporting roots, old spectator blocks, or additional cross-core
contacts do not affect this reasoning. B is an actual connected
principal subdiagram of the whole finite tree; its edges and paths
cannot be removed by restoring other roots. Only a necessary subset
of actual roots is used in the final Gram obstruction.

## Uniform correlation and a strictly positive determinant

Let x,y be the nonnegative native contact-magnitude vectors of c,d
against A39, and put J=J_39. Their supports lie respectively in
I={13,14,15,16} and R={22,23,24,25}, and both vectors are nonzero.
Every nonzero coordinate is at least1. Define

    q=x^T J x,    r=y^T J y,    s=x^T J y.

The positive entries of J give

    q,r >= min_{i in I union R} J[i,i] =91/10.

Put U=sum_{i in I} x_i(i+1), V=sum_{j in R} y_j(39-j).
The cross block of a path inverse has rank one, and here

    s=U*V/40.

For i,k in I, assuming i<=k without loss of generality,

    J[i,k]/((i+1)(k+1))=(39-k)/(40(k+1)) >=23/680.

The right block has the reflected inequality. Nonnegative contact
coordinates therefore imply

    q >=(23/680)U^2,    r >=(23/680)V^2,
    s^2 <=(289/529)q*r.

The source checks every one of the32 within-block inequalities, all16
cross-block factorizations and the eight inverse diagonals. These
finite identities prove the bound for arbitrary real nonnegative
contact magnitudes, not just the values appearing in a catalogue.

Because c.d=0, eliminating the positive definite A39 block from the
actual41-root Gram matrix leaves the two-by-two Schur matrix

    [ 2-q   -s  ]
    [ -s    2-r ].

Both diagonals are at most-71/10. Its determinant has lower bound

    (240/529)q*r -2(q+r)+4.

Write q=91/10+u and r=91/10+v, where u,v>=0. This lower bound becomes

    13674/2645 +(1126/529)(u+v) +(240/529)u*v >0.

Thus the Schur matrix is negative definite. By congruence, the full
Gram matrix has two negative directions. Actual wall roots lie in a
Lorentz space with only one negative direction, a contradiction.
Additional roots or contacts outside this principal matrix cannot
repair it. The argument does not require a nonsingular old core,
a rational gauge for that core, or any finite cutoff on supporting
components or contact labels.

This proves the universal full-star short7 alternative. It supplies
the missing conditional low-core and high-disjoint coverage for a
future deep-distance8 payment change. It does not prove a general
finite7 quadrilateral theorem, and it does not itself assign a full
unit to an ideal short6/7 pair: the existing ideal-cycle/tree transfers
must still be composed with it.

## Verification and limits

Final generation28326f and a separate fresh saved replay58e035 both
terminated exit0/PASS. The source imports no project arithmetic module
and uses exact Fraction arithmetic. It checks path positivity/inversion,
actual-edge projection gates, support/correlation identities and the
strict determinant polynomial. The proof's published order-ten bound
and actual-face facts remain explicit inherited inputs.

Two other agents completed full source/report mathematical reads.
Both passed the minimal-core deletion/component argument, its placement
inside the retained whole finite path, actual-edge support exclusions,
all-label nonzero contact bounds, arbitrary spectators, the uniform
rank-one correlation inequality and the Lorentz index-two contradiction.
They separately checked Chein's printed-page4 statement and the explicit
bridge from the actual minimal cores to its hypothesis, including
disconnected positive-semidefinite deletions and order2 dotted cores.
Neither duplicated the exact arithmetic run or any original census,
and neither review is claimed as a review of the later payment/global
composition. No mathematical gap was identified in the stated scope.

All nine direct bindings, Python syntax, six local links and the stated
primary-paper link were checked. Source, certificate and report are
final/frozen.

No original all-label classification or high-overlap census is
re-enumerated here. No complete changed payment system, finite budget,
new polytope, new ambient dimension theorem, external refereeing or
formal verification is asserted.

    python3 scripts/fresh_universal_deep_finite8_corridor_20261018.py
    python3 scripts/fresh_universal_deep_finite8_corridor_20261018.py --verify
