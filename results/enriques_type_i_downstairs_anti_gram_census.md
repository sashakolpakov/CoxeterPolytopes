# The Type-I downstairs anti-Gram compiler collapses to one generic orbit

**Status:** exact exhaustive criterion; hostile-audit strengthened and frozen.

## Theorem

Start only with Shimada's Type-I (`12A`) downstairs curve-intersection
matrix.  For two quotient curves with intersection `m`, write `a,c>=0` for
the same-sheet and cross-sheet intersections of chosen K3 lifts.  Then

\[
a+c=m,
\qquad
\delta_i\cdot\delta_j=2(a-c)=4a-2m.
\]

Thus an `m=1` edge allows `{-2,2}` and an `m=2` edge allows `{-4,0,4}`.
The Type-I graph has ten edges of the first kind and three of the second, so
there are exactly

\[
2^{10}3^3=27,648
\]

raw anti-Gram assignments.

Exact enumeration and quotient by independent sheet swaps
`delta_i -> +/-delta_i` proves:

1. the 27,648 assignments form exactly **18** sheet-switch orbits;
2. precisely one orbit is negative semidefinite of rank nine;
3. that unique orbit contains exactly 1,024 raw assignments and is exactly
   the geometric `E8(-2)+<-4>` package reconstructed from OSCAR Number 172;
4. precisely one further orbit is negative semidefinite of rank ten, with
   512 raw assignments, but the half-sum glue forces a square-minus-two
   vector in its saturated anti-invariant span, so `L2` eliminates it; and
5. every other orbit has at least one positive direction and therefore
   cannot be the algebraic anti-invariant span of an Enriques K3 cover.

Consequently, among **all** Type-I two-sheet assignments, the downstairs
intersection alphabet, negative-semidefiniteness, half-sum glue, and
saturated rootlessness leave exactly the rank-nine geometric anti-Gram
**up to swapping the two sheets over individual curves**.  In particular,
the conclusion no longer needs a prior Picard-rank-19 assumption.

This is the first finite projection result for the variable Enriques lift
theorem.  It decides the entire Type-I anti-Gram layer.  It does not yet
reconstruct the surviving rank-nine embedding/glue orbit from downstairs
data alone, nor independently replay its period and `L4` gates without the
known OSCAR carrier.

The exact archive is
[`enriques_type_i_downstairs_anti_gram_census.json`](../data/enumerations/enriques_type_i_downstairs_anti_gram_census.json),
rebuilt by
[`certify_enriques_type_i_downstairs_anti_gram_census.py`](../scripts/certify_enriques_type_i_downstairs_anti_gram_census.py).

## Exact switching quotient

For a fixed support graph, sheet swaps act by diagonal sign conjugacy.  The
certifier chooses a canonical representative without a graph library:

1. split the nonzero support into connected components;
2. in each component take the sorted breadth-first spanning tree; and
3. switch vertex signs so that every tree edge is negative.

The signs on non-tree cycle edges are invariant.  This produces a complete
canonical representative.  As an independent check, the script explicitly
generates every diagonal sign conjugate of every representative.  If a
support has `c` connected components, its orbit contains exactly
`2^(12-c)` matrices; this agrees with the multiplicity of every one of the
18 canonical classes.

The orbit sizes are

```text
size  512:  2 orbits
size 1024:  6 orbits
size 2048: 10 orbits
```

and sum with multiplicity to 27,648.

## Exact inertia census

Exact rational inertia is computed on all 18 orbit representatives.  The
raw-assignment and orbit counts are:

| inertia `(positive,negative,zero)` | raw assignments | switch orbits |
|---|---:|---:|
| `(0,9,3)` | 1,024 | **1** |
| `(0,10,2)` | 512 | **1** |
| `(1,9,2)` | 2,048 | 1 |
| `(1,10,1)` | 9,216 | 7 |
| `(1,11,0)` | 512 | 1 |
| `(2,10,0)` | 14,336 | 7 |

The unique `(0,9,3)` representative is compared entry-for-entry, after its
canonical sheet switch, with the anti-difference Gram independently produced
by the explicit Type-I lift package.  They agree exactly.  Its support has
two components, explaining its `2^(12-2)=1,024` raw sheet choices.

The rank-ten representative has three support components and hence
`2^(12-3)=512` choices.  It is eliminated before any period calculation.
In the fixed wall labelling its eleventh and twelfth generators satisfy

\[
\delta_{11}^{2}=\delta_{12}^{2}=-4,
\qquad \delta_{11}\cdot\delta_{12}=0,
\]

while Shimada's corresponding roots obey

\[
r_{11}+r_{12}=(8,4,8,14,12,10,8,6,4,2)\in2L_{10}.
\]

For arbitrary sheet choices `eps_11,eps_12 in {+1,-1}`, the half-sum glue
therefore makes

\[
e=\frac{\epsilon_{11}\delta_{11}+\epsilon_{12}\delta_{12}}2
 =b_{11}^{\epsilon_{11}}+b_{12}^{\epsilon_{12}}
  -\frac{u_{11}+u_{12}}2
\]

an integral anti-invariant vector in the saturation of the generated span.
Its square is

\[
e^2=\frac{-4-4}{4}=-2.
\]

Thus the rank-ten orbit violates `L2` saturated rootlessness.  The formula
works for all four sign choices, so no sheet switch escapes the obstruction.

## What this changes

Before this census, the Enriques lift theorem retained arbitrary
anti-invariant vectors as auxiliary input.  For Type I those vectors no
longer constitute a fishing problem at the Gram level:

```text
downstairs Gram + nonnegative curve intersections
    + negative semidefiniteness + half-sum glue + saturated rootlessness
    => one anti-Gram orbit
    => the known E8(-2)+<-4> package.
```

The remaining compiler steps are now sharply separated:

1. reconstruct the saturated anti-lattice and its glue from that forced Gram;
2. enumerate primitive embeddings into the fixed Enriques anti-lattice;
3. decide the exact period and ample orientation; and
4. run the finite `L4` and complete-cone/Vinberg gates.

For the rank-nine orbit the explicit package has already passed all four.
There is no surviving rank-ten branch.  This is a scoped criterion for one
known chamber, not a new abstract hyperbolic polytope.
