# PB-Atlas Stage B2b: a generated strong-Coble component with all facets

**Status (2026-08-10): proved on one structurally selected raw plane packet.**
This is the first PB-Atlas stage that starts from the numerical contraction
packet, generates the incidence equations and collision obstruction, and ends
with a strong Coble surface carrying every proposed physical facet.  It is a
real geometric compiler positive, not yet the universal raw-chamber theorem.

The production geometry core is
`scripts/pb_atlas_geometry_stage_b2b.py`; the exact calibration is
`scripts/certify_pb_atlas_h10_geometry_stage_b2b.py`, with frozen archive
`data/enumerations/pb_atlas_h10_geometry_stage_b2b.json`.

## 1. Input and selection

The construction consumes one certified Stage-B1 `P2` packet.  The
calibration selects it as the unique packet whose proximity-forest component
sizes are

```text
[2,2,2,1,1,1,1].
```

The construction never reads a surface name, a Coble/Kondō type, stored point
coordinates, an accepted factor, or a stored curve equation.  Its ten centers,
divisor classes, multiplicities, and physical boundary weight come from the
Stage-B1 arrays; the chart and every jet entry come from Stage B2a.

There are three free centers, hence eight `A/B` blowup-chart choices on the
deterministic four-point `PGL3` gauge open.  For every chart, the program
selects precisely the nine noncontracted degree-one square-minus-two divisor
rows and computes their `3x3` evaluation determinants.  Seven charts contain
a nonzero constant determinant and therefore cannot realize all nine root
lines.  The sole surviving chart has modes

```text
center 3: A,   center 6: B,   center 9: B.
```

This is a finite chart calculation, not a remembered orientation.

## 2. The generated incidence ideal and its bad component

On the surviving chart, exact Gröbner reduction gives

```text
q5_0-q8_1, q5_1-1, q7_0-1, q7_1-q8_1, q8_0-q8_1,
(q8_1-1)(q8_1+1), t3, t6, t9.
```

Thus the line-root equations have two reduced rational components.  The
compiler then evaluates the homogeneous coordinates of every proper center
and generates the pairwise collision locus.  At `q8_1=+1`, centers
`4,5,7,8` all coincide, producing all six pairwise collisions.  At
`q8_1=-1`, the seven proper centers are distinct.  Hence collision avoidance,
not a hard-coded sign, selects the unique good component.

The compiler also substitutes both components into the actual selected
affine-patch coordinates and all four (3\times3) determinants of the PGL3
gauge quadruple.  Every declared chart-open factor is nonzero; in particular,
the accepted point is not an artifact of clearing a chart denominator.

This closes the exact false-positive defect of the old prime-binomial pilot:
both signs are saturated binomial components, but only one belongs to the
ordered-cluster open.

## 3. Exact sections, not virtual dimensions

At the collision-free component, the generated matrices have

```text
|-K|:   10 x 10, rank 10, kernel 0;
|-2K|:  30 x 28, rank 27, kernel 1.
```

The B2b core extracts the primitive kernel polynomial of every relevant
matrix and independently transports that polynomial through the chronological
blowups.  At each center it computes the actual order of vanishing and checks
that it equals, rather than merely exceeds, the multiplicity in the divisor
table.  This excludes hidden exceptional components.

The physical facets split as follows.

* Ten contracted walls are the smooth irreducible strict exceptional curves
  produced by reversing the ten blowups.
* Nineteen noncontracted nonboundary rows have one-dimensional section spaces.
  Their degrees are one, two, or three.  Every actual cluster multiplicity is
  exact, and exact absolute factorization proves every plane image
  geometrically integral.
* The remaining wall row is the square-minus-one root (-K), marked by
  physical boundary weight two.  The compiler therefore constructs the
  degree-six physical curve (B=2(-K)=-2K), rather than pretending the
  ineffective wall root itself is a curve.  Its kernel is one-dimensional,
  every cluster multiplicity is exactly two, and exact absolute factorization
  proves it geometrically integral.

The three wall-index sets above are checked to be disjoint and to partition
the complete thirty-row input list; the ten contracted indices are also
checked coefficient-for-coefficient against the source exceptional index of
every reverse blowup.

Each nonboundary class has arithmetic genus zero in the compiled intersection
table.  Geometric integrality plus exact strict-transform class therefore
makes it a smooth rational curve.  For the boundary, `B=-2K`, `K^2=-1`, so
`B^2=-4`, `K.B=2`, and again `p_a(B)=0`; its strict transform is smooth
rational.  Consequently all thirty proposed physical facets are effective
and irreducible.

The rank statements and boundary geometry prove

```text
|-K_S| = empty,
|-2K_S| = {B},
```

with `B` reduced and smooth.  The output is therefore a strong Coble surface,
derived from the raw numerical packet rather than supplied as a witness.

## 4. The dimension theorem carried by this compiler gate

The reusable implication is now unconditional on every packet that passes
these same finite closed-component gates.  If its Picard rank is `n+1`, the
canonical double cover of the resulting strong Coble surface is a smooth K3
surface and contains the adapted rank-`n+1` lattice.  Therefore

```text
for the present Q/C compiler: n <= 19.
```

There is a separate characteristic-not-two analogue: an independently
constructed strong Coble carrier has a K3 cover and hence (n\leq21).  The
current exact chart, Gröbner, and Singular certificate lives over
\(\mathbf Q\) and does not itself construct a positive-characteristic
carrier, so that analogue is not part of this calibration's output.

This rank pushback does not assume a K3 lattice or a rank ceiling in the raw
input.  It also does not require equality with the nef cone: carrier existence
and the full Picard marking already bound `n`.  If the raw input is additionally
reflection-maximal, the independently proved maximal-Coble return theorem
turns the thirty constructed facets into the complete nef chamber and returns
the literal same polytope.

For this calibration, `rho=11` and `n=10`; the value is a nontrivial end-to-end
test of the compiler, not a new numerical exclusion.  The theorem-shaped
payoff is the rank bound for the variable raw class accepted by the same
finite algebraic gate.

## 5. Exact boundary of the result

This theorem does **not** yet prove:

* a second raw packet passes the complete B2b grammar;
* the ruled `F_e` chart compiler or all other plane chart strata;
* reflection maximality of this particular calibration;
* that every maximal arithmetic reflection chamber enters the rational/Coble
  compiler.

Thus it is a strict generated Stage-B2b positive and a real AG dimension
bridge, but not global L4 or a replacement for the remaining essential-
surjectivity theorem.
