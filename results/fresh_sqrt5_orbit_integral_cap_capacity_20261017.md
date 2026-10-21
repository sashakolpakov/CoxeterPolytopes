# Core orbits and arithmetic two-cycles give 815 cap capacities

The [exact certificate](fresh_sqrt5_orbit_integral_cap_capacity_20261017.json)
and [producer](../scripts/fresh_sqrt5_orbit_integral_cap_capacity_20261017.py)
extend the preceding241 cross-coface cap capacities to815, retaining every
old column. Of these,167 have positive load at the new full-model dual,
and19 are violated. This invalidates that witness in the enlarged cap
model; it does not compute a new optimum or prove a dimension bound.

## A core symmetry gives finitely many orientations

Retain the previous actual maximal-simplex/prism incidence interpretation.
The full simplex star spans the ambient Lorentz space and is the orthogonal
sum of a nonsingular Lorentzian core and finite spectators. Restrict to
prism stars whose deletion of the cap b leaves exactly this core and one
isolated A1 root q internally. All appended spectators are orthogonal to b.
Normalize every outward wall normal to squared norm2.

Fix the complete contacts of b with the core, its base-vertex subset B,
and tau=(b,q)^2>4. The core contacts determine a unique projection v and

    b = v - sqrt(tau) q/2,
    h = (v,v) = 2 - tau/2 < 0.

Choose an identification of an actual core with its saved canonical model.
When the core has symmetries, a selected incidence may have any orientation
in the simultaneous orbit of its contact vector and B under the complete
coloured core automorphism group. Let k be the size of that orbit. There
are k distinct pairs of contacts and base subsets, even when the group
itself is larger. In the present census k is1 or2.

Within one orientation, caps at distinct isolated A1 copies q,q' have
mutual normalized product h. If h is forbidden, that orientation occurs
at most once across all those copies. At the same q the spanning star
determines b uniquely; fixed B gives one base face, and the inherited
vertex/coface capacity prevents counting it twice. Therefore the sum of
all selected incidences in this orbit is at most k per actual simplex S.
This is a separate valid signed capacity for each profile. No assertion
is made that the capacities of different profiles are disjoint resources.

The old trivial-symmetry capacities are exactly the k=1 subfamily.
Symmetric cores can also have k=1: the contacts and B may be fixed by
every symmetry. This is responsible for most new violations.

## Two independent reasons why a repeated orientation is forbidden

Here h belongs to Q(sqrt(5)). We use either of the following necessary
conditions on a pair of actual walls.

First, if -2<h<0, the pair must have h=-2cos(pi/m) for an integer m>=3.
The degree of 2cos(pi/m) is phi(2m)/2. Degree at most two restricts m
to3,4,5,6. Indeed phi(N)<=4 forces prime divisors2,3,5 with exponents
at most3,1,1; the surviving N are1,2,3,4,5,6,8,10,12. The checker
reconstructs this list. The values sqrt(2) and sqrt(3) do not belong
to Q(sqrt(5)), so only h=-1 and h=-(1+sqrt(5))/2 remain. Every other
h in this open interval is forbidden. This is a geometric condition
and does not require arithmeticity of the whole chamber.

Second, an arithmetic reflection group satisfies Vinberg's cyclic
integrality criterion. For two norm-2 walls its length-two cycle is the
square of their product. Thus h^2 must be an algebraic integer. This
condition applies also when h<-2, where the walls could geometrically
be disjoint. It is a necessary arithmetic condition, not a sufficient
test for a pair or a chamber to exist.

For the norm-2 convention and the integrality condition, see Theorem13
and the length-two cycles in section5 of Antolín-Camarena, Maloney and
Roeder, [Computing arithmetic invariants for hyperbolic reflection groups](https://arxiv.org/pdf/0708.2109).
That paper applies Vinberg's criterion in dimension3; the necessary
cyclic-integrality condition used here is the same general criterion
already used in the inherited arithmetic census.

For z=a+b sqrt(5), integrality is checked in two equivalent exact forms:
the ring coordinates a-b and2b are integers, and the characteristic
polynomial has integral trace2a and norm a^2-5b^2. A nonintegral trace
or norm proves that z is not an algebraic integer. Passing to a larger
number field cannot change that fact. The criterion is applied to h^2,
not to h or to the original dotted parameter by mistake. The same
normalization is used by the inherited exact arithmetic star census.

We retain a selected occurrence when either the finite-angle condition
fails or h^2 is nonintegral. Parallel products and permitted finite or
dotted products are left out unless the integrality test excludes them.
No omitted profile receives a positive coefficient in these columns.

## Reconstruction and exact diagnostic

Among1,827 isolated-A1 base occurrences,1,100 have the required sole
internal A1 extra. Of these,815 pass at least one obstruction test:
727 fail the finite-angle condition and635 fail cyclic integrality;
these sets overlap. The old open interval(-1,0) accounts for690.
The resulting815 columns have572 capacities1 and243 capacities2.

For every selected occurrence, the producer checks the raw core map,
orthogonality of q to the core, the complete contact data, the base
subset, and the action of each saved automorphism. It computes the
simultaneous orbit exactly and groups only identical orbit profiles.
Completeness of the saved automorphism groups is inherited from the
frozen exhaustive coloured-core census. Every old241 column is compared
coefficient by coefficient to its retained successor.

The diagnostic uses the exact dual from the
[new full model](fresh_sqrt5_cross_coface_cap_full_model_20261017.md),
whose bound is12.69277669282928. Seventeen of the19 violated capacities
involve nontrivial core symmetry. The other two use genuinely new
arithmetic information: tau=(13+5sqrt(5))/2 gives

    h = -5(1+sqrt(5))/4,
    h^2 = (75+25sqrt(5))/8.

The latter is not integral; its trace is75/4. These are base types87981
and87997, with containing-simplex rows278 and42. Their caps could have
been mutually dotted by the real-signature test alone, so the finite-angle
gap argument would not have excluded their repetition.

All arithmetic is rational or exact Q(sqrt(5)) arithmetic. The producer
invokes no optimizer. Generation and fresh saved replay both pass.
Repairing the face masses
requires recomputing affected vertex and local/joint charges. No old
fixed-mass cover is imported as a certificate at changed masses.
