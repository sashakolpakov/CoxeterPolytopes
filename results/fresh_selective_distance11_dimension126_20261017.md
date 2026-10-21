# Selective distance 11 and the unrestricted dimension 126 bound

Every finite-sided, finite-volume real
hyperbolic Coxeter polytope has dimension at most **126**. No arithmeticity,
compactness, simple-cusp, cusp-factor-rank, or facet-count restriction is
imposed. The compact-facet consequence gives ambient dimension at
most **95**.

The change from the [completed distance-12 proof](fresh_selective_distance12_dimension136_20261017.md)
is to pay pairs at distance at most 11, while still paying every pair in
a connected corner component of order at most 17. Two new overlap types
appear. One has an actual radical relation that isolates the component.
The other admits a supporting-root extension, but the extension lies in
the original core span and can add at most one node.

The [primary source](../scripts/fresh_selective_distance11_dimension126_20261017.py)
and [certificate](fresh_selective_distance11_dimension126_20261017.json)
record the finite classification and scalar inequalities. The
[separate review](fresh_selective_distance11_dimension126_20261017_review.md)
checks the geometry and reconstructs the arithmetic independently.

## Actual face geometry and exact obstruction certificates

Use the actual-face and cusp-factor incidence statements of the preceding
proof. For a quadrilateral, let W be its supporting white roots, and let
U and V be its two opposite pairs of black side roots. The diagram W,
and W together with any one black root, is positive definite. The diagram
W with an adjacent black pair is positive semidefinite. Each of W+U and
W+V is indefinite. Choose minimal indefinite cores L and M in these two
diagrams. Each contains its two marked black roots, and deleting either
black root leaves a positive-definite diagram.

The inherited complete classification has 26 diagrams of order at least
6, giving 384 marked contexts. Their maximum distance from a marked pair
is 7. Cores of order at most 5 have marked radius at most 3, with no
restriction on their finite labels. Assume there is no path from U to V
of length at most 11 with white internal nodes.

For a positive-definite black-deletion Gram A, an additional supporting
white contacting its node i has squared projection at least (A inverse)ii.
The inverse is entrywise nonnegative, all contacts are nonpositive, and
every nonzero norm-two Coxeter contact has magnitude at least 1. If either
black deletion has inverse diagonal at least 2, its extension contradicts
positive definiteness of W with the remaining black. Equality 2 is
excluded as well. The complete test of 200 ports at marked distance at
least 4 leaves 55 allowed ports.

Some rejections use an explicitly nonnegative vector x with
x^T A x <= 0 in a diagram required to be positive definite. These are
exact certificates with a monotonicity property: additional nonpositive
contacts or larger bond magnitudes can only decrease x^T A x. Thus the
certificate remains valid when unspecified contacts are restored.

The certificates are obtained by symmetric elimination of a Z-matrix
(a symmetric matrix with nonpositive off-diagonal entries). At every
positive pivot q, the remaining basis vectors change by subtracting
nonpositive multiples of the pivot vector. Starting from coordinate
vectors, all residual basis vectors stay nonnegative. The first
nonpositive pivot therefore supplies the required vector. The source
checks its entries and its quadratic value exactly. Labels 4 are
rationalized by a positive diagonal rescaling with entries 1 and sqrt(2);
this preserves inertia and nonnegativity of the witness. No approximate
eigenvalue is used for an exclusion.

## Disjoint cores

Two disjoint indefinite cores cannot be orthogonal in an ambient space
with negative index one. A cross bond must connect white ports whose
marked distances sum to at least 11. Both distances are at least 4,
because each core radius is at most 7. This also excludes cores of order
at most 5 from the present disjoint case.

The 55 allowed ports give 425 ordered marked joins. In 119 joins, the
minimum Gram with a simple cross bond already has a nonnegative
nonpositive-square vector in the wholly white block. Those joins are
impossible for every allowed magnitude of that bond. Among the remaining
306 joins, each ordered marked-context pair admits exactly one cross
edge. Every other potential cross edge is already excluded by a short
path, an inverse test, or the preceding positive-vector certificate.

The remaining cross bond must be simple, including when the native
diagrams contain labels 4. Each endpoint has a three-node wholly white
path starting at that endpoint. Give the consecutive nodes coefficients
1, 2/3, 1/3. Its square is at most 4/3, because each native bond magnitude
is at least 1. Combining the two paths gives square at most

    8/3 - 2 beta,

where beta is the cross-bond magnitude. White positive definiteness
implies beta < 4/3. The norm-two Coxeter alphabet has no nonright bond
between the simple magnitude 1 and sqrt(2), so beta=1. Parallel and dotted
bonds are also excluded. This argument has no finite-angle cutoff.

The exact full Grams for all 306 remaining joins have two negative
directions. Rational congruence certificates verify this inertia. They
exhaust 13 underlying unmarked joins, all with one endpoint L25/node 9:

| Other diagram and port | Ordered marked occurrences |
|---|---:|
| L10/node 5 | 6 |
| L13/node 6 | 6 |
| L14/node 6 | 6 |
| L16/node 7 | 6 |
| L18/node 7 | 42 |
| L20/node 8 | 6 |
| L21/node 8 | 6 |
| L22/node 8 | 54 |
| L22/node 4 | 42 |
| L23/node 9 | 6 |
| L24/node 9 | 6 |
| L25/node 9 | 114 |
| L25/node 3 | 6 |

Since no additional cross bond remains possible, these are the actual
core Grams. Negative index two contradicts the ambient Lorentzian space.

## Exhaustive overlap classification

A shared white root has its two marked distances sum to at least 12.
Each distance is therefore at least 5. In the complete marked catalogue
each context has at most three such nodes. Enumerating injective partial
identifications of these nodes exhausts all possible overlaps.

There are 2,137 candidate identifications before injectivity and common
Gram-entry consistency, and 1,010 after those checks. In 856 cases a
shared node has an external white neighbor from the other core and a
black-deletion inverse diagonal at least 2. The inverse test excludes
them. Of the remaining 154 cases, 40 have a nonnegative
nonpositive-square vector in W with one black root, which must be
positive definite. Restoring other nonpositive contacts cannot repair
that obstruction.

The remaining 114 ordered marked contexts are precisely:

| Overlap | Contexts |
|---|---:|
| Two L25 cores sharing their far node 9 | 60 |
| Two L25 cores sharing the terminal edge in reverse | 27 |
| Two L25 cores sharing their terminal three-node path in reverse | 9 |
| L22/node 8 identified with L25/node 9, in either orientation | 18 |

Every unspecified bond between the exclusive parts is checked: 8,280
create a marked path of length at most 11, and the other 45 have exact
nonnegative nonpositive-square witnesses in W with one black root. All
45 occur in shared-vertex contexts. This determines
the actual core Gram in all 114 surviving contexts. The following
extension arguments therefore apply to the actual normals, rather than
to a freely chosen zero-completion.

## Shared vertex: the larger marked domain

The 19-root L25/L25 shared-vertex Gram is the same as in the completed
distance-12 proof. Number the first copy 0 through 9 and map the second
copy's 0 through 8 to 10 through 18, with node 9 shared. Its actual
radical relation is

    z=(-6,-3,-4,-2,-5,-4,-3,-2,-1,0,
        6,3,4,2,5,4,3,2,1).

To see that this is an actual root relation, its Gram-null combination
is orthogonal to the core span, which contains a timelike vector, and
has square zero. The orthogonal complement of that timelike vector is
positive definite, so the combination vanishes.

In each copy choose the preferred marked root as follows: choose 4 if
it belongs to the marked pair, then 0 if it belongs; otherwise choose 2
from {1,2} or {2,3}, and 1 from {1,3}. The complete preferred-corner
calculation gives a marked component of order at most 17 in all 60
contexts.

For each finite corner component, its inverse diagonal restricts a new
supporting white's possible contacts to entries at most 2. For an affine
component, its positive kernel forces all contacts to that component to
vanish. Include the two omitted black roots separately. Every allowed
contact has nonzero z-coordinate; all negative-z contacts are within
distance 3 of U, and all positive-z contacts within distance 3 of V.

Pairing an additional supporting white with z forces contacts of both
signs if there is any contact. That creates a white-interior marked path
of length at most 8. Hence every additional supporting white is
orthogonal to the core, and the preferred marked component is isolated.

## Shared edge: the larger end sets

The 18-root shared-edge core maps the second L25 copy's 0 through 7 to
10 through 17, its 8 to 9, and its 9 to 8. It is nonsingular Lorentzian,
with determinant -1. The two-center Schur complement is

    [[-5/66,-1/11],[-1/11,-5/66]].

Use the same preferred-root choices. All 27 preferred marked components
have order at most 16. The finite inverse and affine kernel restrictions
confine contacts of an additional supporting white to

    E_left={0,1,2,3}, E_right={10,11,12,13}.

Every allowed end node is within distance 2 of its own marked pair.
Contacts at both ends give a marked path of length at most 6. If the
contacts are supported at just one end, the corresponding principal
block of the inverse of the full core Gram is

    [[30,15,20,10],
     [15, 8,10, 5],
     [20,10,14, 7],
     [10, 5, 7, 4]].

It has strictly positive entries and diagonal minimum 4. Thus a nonzero contact vector b has
b^T G^{-1}b >= 4.

The core span contains the unique negative direction, so its ambient
orthogonal complement is positive definite. A norm-two additional root
would need squared orthogonal residual 2-b^T G^{-1}b >= 0, a contradiction.
All supporting whites outside the core are orthogonal to it.

## Three-node overlap: the extension is controlled

The new 17-root core maps the second L25 copy by

    (0,1,2,3,4,5,6,7,8,9) -> (10,11,12,13,14,15,16,9,8,7).

Its determinant is -2 and its signature is (16,1). Both black pairs are
one of {1,2}, {1,3}, {2,3}. Choose the preferred node 2, 1, 2 respectively.
The [separate exact extension lemma](fresh_distance11_extension_review_20261017.md)
proves that each preferred marked component initially has order 15, and
that a new supporting white's contacts can lie only in

    {1,2,3} or {11,12,13}.

A contact at both ends creates a marked path of length at most 4. The
full inverse Gram on either one-end set is

    [[7/2,4,2],[4,6,3],[2,3,2]].

Its strictly positive entries and the norm-two residual inequality imply
that the only nonzero possibilities are a sole contact -1 at node 3,
or a sole contact -1 at node 13. Each attains squared projection 2, so
its orthogonal residual vanishes. The actual roots are uniquely
r=-G^{-1}e_3 and s=-G^{-1}e_13, in the core span. Their mutual product is
-2. They cannot both belong to the positive-definite supporting diagram
W. Duplicate normalized roots do not define distinct facets.

Consequently at most one extra root occurs. Every other supporting white
orthogonal to the core is also orthogonal to that extension, because it
lies in the core span. The preferred marked component remains isolated
and has order at most 16. The argument explicitly permits this extension;
it does not assume that every surviving core must be isolated unchanged.

## Mixed L22/L25 overlap

Number L22 by 0 through 8, and map L25 by

    (0,1,2,3,4,5,6,7,8,9) -> (9,10,11,12,13,14,15,16,17,8).

The [mixed extension lemma](fresh_distance11_extension_review_20261017_mixed.md)
gives inertia (16,1,1) and the actual radical relation

    z=(-4,-2,-3,-2,-1,-3,-2,-1,0,6,3,4,2,5,4,3,2,1).

Each black pair is one of {1,2}, {1,3}, {2,3} in its original copy. Use
the preferred choices 2, 1, 2. The nine contexts, and their nine reversed
orientations, have marked corner component order at most 16. The full
corner restrictions exclude the zero-radical node from the contact
domain. Allowed negative-z contacts lie within distance 2 of the left
black pair; positive-z contacts lie within distance 1 of the right pair.

The actual radical relation forces both signs for any nonzero additional
supporting-white contact, giving a marked path of length at most 5.
Every additional supporting white is therefore orthogonal to the core.

## From local isolation to actual corner payment

At finite corners, the resulting isolated marked component is a component
of the actual ambient vertex diagram. At ideal corners, the inherited
factor-fiber incidence theorem places the entire varying affine factor
inside W plus the two side roots. An isolated finite component cannot be
that varying factor; an isolated affine component is already the whole
factor. Omitted spectator cusp walls do not enlarge the marked component
relevant to the face payment. Thus the proof includes nonsimple cusps.

Pay weight 1 on a marked pair if its distance is at most 11 or its
connected corner component has order at most 17; otherwise pay 0. The
complete classification above pays every quadrilateral at least 1.
The inherited minimal-core argument gives every triangle two paths of
length at most 8, hence payment at least 2. Nonnegative weights suffice
for polygons with at least five sides. Every k-gonal two-face therefore
receives total weight at least 5-k.

## Exact dimension inequalities

Components of order at most 17 have at most 8 times their order in pairs.
For larger finite and affine components, the inherited all-rank
path/cycle/fork counts give at most the component order at each distance.
The finite vertex cost is at most 11n. At an ideal vertex with affine
factor ranks r_i and link vertex count V=product(r_i+1), the factor-fiber
multiplicity gives

    C(v)/V <= sum_i min(11,r_i/2).

The inherited shape inequality has normalized ideal curvature

    kappa/V = 3d/4 - (5/2) E|sum_i X_i|,   d=n-1,

where X_i is centered uniform on r_i+1 consecutive integers. The proved
all-rank convex-order transfer applies with cap 22. Concentrating two
ranks below 22, or two above 22, preserves the cost and cannot increase
curvature. Increasing the sum of squared ranks proves termination. The
reduced profiles consist of copies of 22, at most one part from 1 through
21, and at most one part at least 23.

For d=126 through 135, all 1,095 reduced profiles have positive
curvature-minus-cost. Integer convolution and a separate moment formula
give the same exact values. The minimum is

    40603/559682,

at ranks (22,22,22,22,22,16). The finite slack beta_n-11n is positive for
n=127 through 136, and the global shape inequality has its strict
positive term 5(n-1). These dimensions are therefore impossible. The
completed 136 theorem excludes all higher dimensions, proving 126.

At ambient n=126 the reduced profile (22,22,22,22,22,15) instead has
margin -238535/12872686. These same scalar weights do not prove 125.
This is a local link profile, not a finite-volume construction.

For a compact facet of dimension d=n-1, the finite payment gives
beta_d < 11(d+1). For d=95+2t and d=96+2t, eight times the opposite
difference has coefficient lists (12,192,4) and (104,196,4), respectively.
Both are positive for t>=0. Thus d<=94 and the ambient dimension is
at most 95.

The finite/quasi-Lanner classification, actual incidence statements,
distance counts, convex-order theorem, and shape inequality remain
explicit inherited inputs. Exact finite arithmetic does not by itself
formally verify these geometric arguments. No finite-volume realization,
external refereeing, proof-assistant formalization, or publication-priority
claim follows from this package.

Primary generation 85083 and fresh saved replay 10135 both terminated
with exit code zero and PASS. The separately implemented Fraction
reconstruction agrees on all 200 tested ports, 425 joins, 1,010 consistent
overlaps, 114 preferred-corner cases, 8,325 extra-contact exclusions, and
1,095 profiles. Its written mathematical review finds no gap. These are
internal agent-workflow checks; the report does not claim external
mathematician refereeing or formal verification.

The primary source and certificate are frozen at SHA256
`8fcbea5f00ae2b63cd3e0b7fe3fa8616d77b7ed559c321af0648eccfc31bfb66`
and `f6340fcf99d3884b4a22113fc446cef5ac8bd6880f0b1354f525f2bcaab47646`,
respectively. All 25 direct bindings and source syntax passed. The
separate review records its own saved replay and final bindings.
