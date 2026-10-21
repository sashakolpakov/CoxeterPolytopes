# Realizing a singular acute star as a compact quadrilateral face

Let G be a symmetric real Gram matrix of n spacelike roots, with positive
diagonal entries and nonpositive off-diagonal entries. Distinguish four
black roots in cyclic order0,1,2,3 and let W be the remaining n−4 white
roots. Suppose:

1. G has inertia(n−2,1,1).
2. Each principal block W together with a cyclically adjacent black pair
   is positive definite.
3. The kernel of G contains no nonzero vector with all coefficients
   nonnegative.

Then the roots in the nondegenerate Lorentz quotient define a chamber with
a strict timelike interior point. Its common-white face is a **compact
quadrilateral with exactly the four prescribed finite vertices**. Every
listed wall is an actual facet near at least one of those vertices. This
is a local realization statement, not a finite-volume assertion about the
whole chamber.

## A common strict center exists

Decompose the nonzero-off-diagonal graph into connected components. Exactly
one Gram block has a negative eigenvalue. In that block, a unit vector
minimizing the Rayleigh quotient may be replaced by its coordinatewise
absolute value, since all off-diagonal entries are nonpositive. A zero
coordinate of a nonnegative minimizing eigenvector would force every
neighbor coordinate to vanish; connectedness rules this out. Thus its
least-eigenvalue eigenvector v is strictly positive. The root combination
x with coefficients v has negative square and strictly negative products
with every root of this block.

Every other block is positive definite. Otherwise a singular positive
semidefinite connected Z-matrix has a nonnegative kernel vector by the same
least-eigenvector argument, contrary to hypothesis3. In each positive
block use the root combination with coefficients −G_block^{-1}1; its
products are all−1. Add sufficiently small positive multiples of these
orthogonal combinations to x. The result h stays timelike and has strictly
negative products with every root. Choose its time orientation as future.

## Projection gives all four actual finite corners

A positive-definite Z-matrix has entrywise nonnegative inverse. One proof
chooses a scalar t at least its largest eigenvalue and writes
G^{-1}=t^{-1}Σ_{k≥0}(I−G/t)^k; the matrix in parentheses is entrywise
nonnegative and has spectral radius below1.

Project h orthogonally to any principal spherical block S. Its projection
onto span(S) has nonpositive coefficients, since its contact row is
negative and G_S^{-1} is nonnegative. Subtracting that projection makes
the products with every other wall strictly more negative or unchanged.
The squared norm decreases by a nonnegative quantity; its product with h
is negative. Hence the residual vector is future timelike, lies on S,
and is strict at every unselected wall.

Apply this first with S=W and then with S=W plus each adjacent black pair.
The latter blocks have rank n−2 in a Lorentz space of dimension n−1, so
the resulting vectors are the four actual finite corner rays. The active
normals at each corner are independent, and no other root is incident.

## The whole face cone has only those four rays

Eliminate the positive white block. The four projected black normals have
a symmetric Z-Gram of rank3, inertia(2,1,1), positive diagonal, and
positive-definite adjacent two-by-two principal blocks. Its kernel lifts
to the full kernel by the nonnegative matrix −G_W^{-1}G_WB. Thus a
nonnegative residual relation is impossible.

In a mixed relation, move the positive and negative terms to opposite
sides. A sign class consisting of one root has positive squared norm,
whereas its product with the other side is nonpositive; that is
impossible. With four roots, the two sign classes must therefore each
have size two and no coefficient can vanish. An adjacent pair has
positive-definite Gram, so the same argument excludes an adjacent
same-sign pair. The signs consequently alternate: +,−,+,−, up to reversal.

In the three-dimensional common-white subspace, taking the four black
products identifies the space with the hyperplane k·t=0, where k is this
alternating kernel vector. The face cone is exactly its intersection with
the four inequalities t_i≤0. Its extreme rays have one nonzero coordinate
from each sign class of k. There are exactly four such rays, and their
zero coordinates are precisely the cyclically adjacent black pairs.
They are the four future timelike rays already constructed. Their
projective convex hull is a compact hyperbolic quadrilateral. The cone
description also excludes hidden ideal or exterior rays.

## Arithmetic realization from integral cycles

Suppose additionally that G is defined over a totally real field K, its
diagonal entries are totally positive, and every other real embedding
gives a positive-semidefinite matrix of corank one. The K-radical is
common under embeddings, so the quotient is an admissible K-form. For
the present applications K=Q(sqrt5); this form is anisotropic over K
because its other real embedding is positive definite.

Set a_ij=2G_ij/G_jj. Reflection in root j sends root i to
r_i−a_ij r_j. Assume all directed cyclic products of these coefficients
are algebraic integers. Form the O_K-module generated by

    a_i0,i1 a_i1,i2 ... a_i(k−1),ik r_ik

over every simple directed path, including all length-zero roots. There
are finitely many such generators. Every walk with repeated vertices
reduces to a simple path times closed-walk products, and each closed walk
splits into simple cycles. The cyclic-integrality assumption therefore
puts every walk-scaled root into this module. Applying a root reflection
appends a path edge, so preserves the module. It is finitely generated,
torsion-free, spans the quotient, and is a full O_K-lattice. No principality
or free-basis assumption is needed.

Thus integral cycles suffice for a stable integral reflection module;
one must not require or silently assume that every individual a_ij is
integral. The existing exact cycle checker multiplies
2G_ij/G_ii around each undirected simple cycle, including two-cycles.
This gives the same cyclic product as the destination-denominator
convention above, because both denominators are the product of the cycle's
vertex norms. Reversing a cycle also leaves the product unchanged.

Orthogonal positive-definite finite Coxeter blocks may be added as white
spectators. They preserve the old four-ray face cone, admissibility, and
the stable module; small strict-center components in these blocks give
a common timelike interior. This permits exact local quadrilateral
families at every larger ambient dimension allowed by the spectator rank.

This lemma and its arithmetic module argument were developed and reviewed
within the agent workflow. Matrix hypotheses and cycle coverage must still
be checked for each application; no arithmeticity or global compact
realization is inferred from a floating-point Gram test.
