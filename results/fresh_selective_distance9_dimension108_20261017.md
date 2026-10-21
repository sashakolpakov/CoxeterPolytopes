# Selective distance nine and controlled marked components

**Complete internally checked computer-assisted proof.** Primary generation
and a fresh saved replay pass. A separate mathematical/source review and
independent exact reconstruction agree with the complete high-order census,
all preferred corners, all 4,371 saved obstruction witnesses and the scalar
bounds. The linked review records its own final binding and replay status.
No gap has been identified. External refereeing and Lean formalization
have not been performed.

The argument below gives dimension at most **108** for every finite-sided,
finite-volume real hyperbolic Coxeter polytope. It has no arithmetic,
compactness, simple-cusp, affine-factor-rank or facet-count hypothesis.
Its compact-facet consequence is ambient dimension at most **79**.

The [primary source](../scripts/fresh_selective_distance9_dimension108_20261017.py)
and [certificate](fresh_selective_distance9_dimension108_20261017.json)
record the high-order finite split and scalar inequalities. The
[complete small-core lemma](fresh_distance9_low_core_20261017.md)
supplies the additional low-order cases. The new supporting-component
argument has a [separate local review](fresh_distance9_supporting_extensions_review_20261017.md).
The [independent full review](fresh_selective_distance9_dimension108_20261017_review.md)
checks the new step. The completed
[dimension-116 proof](fresh_selective_distance10_dimension116_20261017.md)
supplies the inherited actual-face and cusp incidence, all-rank graph
counts, shape inequality and convex-order transfer, and the excluded
higher-dimensional tail.

## Actual-face setup and complete core domain

For an actual quadrilateral, let W be its supporting white roots and
U,V its opposite black pairs. Normalize roots to square two. W and W
plus any one black root are positive definite; W plus an adjacent black
pair is positive semidefinite. Each opposite pair gives an indefinite
diagram. Choose minimal indefinite cores L in W+U and M in W+V. Each
contains both its black roots; deleting either black root is positive
definite. Their proper principal subdiagrams are positive semidefinite.

Assume no U-to-V path of length at most nine has white internal nodes.
The complete high-order catalogue, of orders six through ten, has 384
marked contexts and maximum marked radius seven. Orders at most three
cannot participate in a disjoint unpaid join or overlap. Two cores of
order at most five cannot do so either. The remaining small/high cases
require a separate classification; they are not discarded by the old
distance-ten radius argument.

The small-core lemma proves that all order-four labels are at most six,
and all order-five labels at most five. It exhausts all weighted graphs,
obtaining 32 and 14 canonical minimal-indefinite types, respectively.
Black-deletion and port tests leave 104 low-first joining-label cases.
Eleven fail white positive definiteness; the other 93 have two negative
directions. Six multiple-edge choices have a white affine cycle. Of 30
low-first overlaps, 24 fail the inverse test and six have an isolated
marked corner of order at most 12. Reversing the cores covers the other
orientation. Both fork/radical and label-five/nonsingular overlap cores
are handled over all further supporting whites.

## High-order disjoint cores

An external supporting white contacting a node of a positive-definite
black deletion has squared projection at least that normalized inverse
diagonal. Thus a simple contact requires both deletion inverse diagonals
strictly below two. Among 1,262 high-order ports at marked distance at
least two, exactly 378 pass. Each distant cross edge has endpoint
distances summing to at least nine, giving 4,442 ordered possibilities.

For a simple edge between white ports x,y, eliminating the two white
blocks gives the necessary condition

    (G_WL inverse)xx (G_WM inverse)yy < 1.

It excludes 1,103 joins. For every remaining join, either one endpoint
black-deletion inverse diagonal is at least one, or twice the displayed
product is at least one. Therefore every magnitude at least sqrt(2)
violates an actual-edge or white positive-definiteness condition. All
labels at least four, parallel contacts and dotted contacts are excluded.
Every remaining cross edge is simple, without a label cutoff assumption.

The 3,339 remaining joins have 85 underlying single-edge Gram types,
each with exactly two negative directions. Some context pairs admit
several candidate edges: 2,007 admit one, 348 admit two, 180 admit three
and 24 admit four. All 1,332 subsets with at least two edges have an
explicit nonnegative nonpositive-square witness in W plus one black.
Extra negative bonds only strengthen these witnesses. Every nonempty
cross-edge set is therefore impossible.

The producer rationalizes each required positive-definite principal
diagram separately. A full joined indefinite graph may have an odd
label-four cycle even though such a principal diagram is a forest. If
a principal graph itself has no rational sqrt(2) gauge, replacing its
nonright bonds by simple bonds gives an entrywise larger Z-matrix; a
nonnegative obstruction for that relaxation remains valid for the native
Gram. The certificate records the gauge or relaxation used for each
witness. No irrational Gram is silently treated as rational.

## High-order overlaps and fixed core Grams

Every shared white has marked-distance sum at least ten. Both distances
are at least three. The source exhausts all injective partial
identifications of eligible nodes, checking every contact on shared
nodes. Of 312,316 candidates, 20,718 are consistent; inverse tests
exclude 18,721 and actual-edge witnesses exclude 1,048 more.

The remaining 949 marked contexts have 54 oriented literal core types.
All four adjacent corner decompositions are checked as finite/affine.
An extra bond between exclusive core parts either creates a path of
length at most nine (60,831 tests), or contradicts actual-edge positive
definiteness (888 tests). Thus the full Gram of each surviving core is
fixed before the extension argument is used.

At a selected corner, a further supporting white is orthogonal to every
affine component by its positive kernel. On each finite component, a
contact is possible only at normalized inverse diagonal at most two.
Omitted black roots are also retained in the possible contact domain.
This accounts for contacts anywhere in the core.

The ordinary argument handles 880 contexts. In a singular indefinite
core, its radical is an actual root relation in a Lorentzian realization.
The chosen contact domain avoids zero radical coordinates; balancing its
two signs forces a forbidden short path. In nonsingular cores, projection
onto the core leaves a positive-definite residual space. Positive inverse
end blocks and short-path tests force a single simple contact. Pairwise
residual Cauchy tests allow at most one core-contacting white; the entire
indirect tail is bounded by the lemma below.

## A rank bound for all indirectly attached roots

Let E be a connected positive-definite Coxeter Gram of order m, normalized
to diagonal two. Every diagonal entry of E inverse is at least m/(m+1).
To prove this, write E=2I-A with A entrywise nonnegative. Positive
definiteness implies spectral radius(A)<2. Removing bonds to a spanning
tree and replacing its bond magnitudes by one decreases A and therefore
decreases the inverse entrywise, by the convergent Neumann series.

At a chosen root of that unit tree, let its branches have orders m_i.
Induction and Schur complementation give

    (E inverse)root,root >= 1 / (2 - sum m_i/(m_i+1)).

Since sum m_i=m-1 and sum m_i/(m_i+1)>=(m-1)/m, the required lower
bound follows. The base case is the one-node Gram [2].

If a root of residual square rho contacts several orthogonal supporting
components, with total order m, its projection cost is at least the sum
of their order/(order+1) bounds, hence at least m/(m+1). Nonpositive
multiple contacts cannot decrease this cost because inverse entries are
nonnegative. Therefore, for rho<1,

    m <= floor(rho/(1-rho)).

This bounds the entire contacted components, including every indirect
tail. It is not a fixed-rank enumeration. Equality can occur for a
simple A_m path contacted at an endpoint; rho=1 must be treated separately.

## The new relevant-component argument

The remaining 69 contexts have a preferred marked component of order at
most 16 after all extensions. Some roots outside that component can have
residual square one, so bounding every external root by the preceding
lemma would be invalid. Only attachments that can enlarge the marked
component need to be controlled.

The exact inverse tests split each contact domain into two positive
inverse end blocks. Any contact at both ends gives a short marked path;
two contacts at one end already have projection square greater than two.
Thus each external core-contacting white has one contact. All candidate
projection squares are at least one. A projection-one node can have a
simple or label-four contact; every larger label exceeds the norm budget.
Every node that can directly enlarge the selected marked component has
projection square strictly greater than one, hence only a simple contact.

If there is no such relevant contact, the marked component is isolated,
regardless of external attachments elsewhere. Otherwise call its new
white G, with residual vector g of square rho<1. Cauchy inequalities and
the Coxeter gap between -1 and zero exclude every other relevant root.
They also exclude all label-four external contacts in the presence of G.
At most one possible outside simple-contact node remains. Call its root C
and residual vector c. It lies outside the marked component and must be
orthogonal to G: the allowed inner-product interval has lower endpoint
greater than -1, while acuteness requires a nonpositive Coxeter product.

Write P for the positive product of the two core projections. Then
g.c=-P. If c has square less than one, two roots of its type cannot
coexist by the same positive-product test. If c has square one, two such
roots would have to be orthogonal with opposite unit residuals. Their
common product -P with g makes this impossible. Thus at most one C
coexists with G.

The exact surviving residual cases are:

| G projection square | C projection square | P | det S | (S inverse)GG |
|---:|---:|---:|---:|---:|
| 4/3 | 4/3 | 2/3 | 0 | — |
| 3/2 | 1 | 1/2 | 1/4 | 4 |
| 7/4 | 1 | 1/2 | 0 | — |
| 5/4 | 1 | 1/2 | 1/2 | 2 |

Here S is the Gram of g,c, with negative off-diagonal entry -P. If C is
absent, every remaining white is core-orthogonal and the rank lemma
bounds G's entire tail. If det S=0, g and c have a positive linear
relation; nonpositive products force every core-orthogonal white to be
orthogonal to both, so no tail can join them.

If det S>0, its inverse is entrywise positive. A core-orthogonal norm-two
white contacting G has projection cost at least (S inverse)GG. The value
four excludes it. At the value two, equality forces contacts (-1,0) to
G,C and puts that white entirely in their residual span, with negative
coefficients. Its product with every other core-orthogonal supporting
white is then nonnegative; acuteness forces zero. It is also orthogonal
to C. Thus this one possible tail root cannot lead to a longer tail or
to C's other attachments. This explicitly excludes an indirect bridge
through an otherwise unbounded projection-one attachment.

All preferred marked components in the 949 high-order cases have order
at most 17, and the small-core lemma gives at most 12. At ideal corners,
the inherited factor-fiber theorem places the entire varying affine
factor inside W plus the two marked sides. Omitted spectator cusp walls
cannot enlarge this component. No simple-cusp assumption is made.

## Payments and the dimension conclusion

Pay one on a corner pair when its distance is at most nine, or when both
marks lie in a connected component of order at most 17. The complete
classification pays every quadrilateral at least one. The inherited
triangle argument supplies two paths of length at most eight, hence pays
each triangle at least two. Nonnegative weights suffice for larger
polygons, giving total payment at least 5-k on a k-gonal two-face.

The inherited all-rank finite/affine path, cycle and fork counts give
finite cost at most 9n. Small components of order at most 17 have total
pair cost at most eight times their order, so the extra payment does not
break that bound. At a cusp with affine ranks r_i and link vertex count
V=product(r_i+1), cost divided by V is at most sum min(9,r_i/2).

The shape inequality and convex-order transfer apply with cap 18. For
link dimensions d=108 through 115, the reduced profiles consist of copies
of 18, at most one smaller positive part and at most one larger part.
All 756 profiles have positive curvature-minus-cost, with minimum
45873/2476099 at six copies of 18. The finite slack is also positive in
ambient dimensions 109 through 116. The completed 116 theorem excludes
the higher tail, giving dimension at most 108.

At ambient 108, the profile (18,18,18,18,18,17) has margin
-4659/260642. The same scalar weights do not prove 107; that profile is
not a finite-volume realization.

For a compact facet of dimension d=n-1, finite payments imply
beta_d<9(d+1). At d=79+2t and 80+2t, eight times the opposite difference
has coefficient lists (12,160,4) and (88,164,4), respectively. Both are
positive for all t>=0. Therefore d<=78 and ambient dimension is at most 79.

The finite classifications and numerical inequalities are exact. The
written geometric and all-rank arguments remain explicit mathematical
inputs. This package makes no construction, external-refereeing,
proof-assistant or publication-priority claim.

Primary generation (session 41240) and the fresh `--verify` replay
(session 2981) both terminated with exit code zero and
`PASS_SELECTIVE_DISTANCE9_CONTROLLED_EXTENSIONS_DIMENSION108`.
The replay recomputes the complete certificate and compares it with the
saved JSON. The low-order package has its own completed generation and
fresh replay; its separate mathematical/source review does not claim a
second independent implementation of that entire enumeration.
