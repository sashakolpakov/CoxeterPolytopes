# Internal root review of the exact cusp-link curvature inequality

The [shape-curvature proof](fresh_cusp_link_shape_curvature_20261017.md)
passes this internal mathematical/source review. It applies to every
finite-sided finite-volume real hyperbolic Coxeter polytope of dimension
n>=3. It is a counting inequality, not by itself a dimension theorem.

The reviewed final files have SHA256 values:

* source: `98e501cc2d61fbe34f7f63aa2644d1421ce2ec7a4b5197eb336e07d144739ee7`;
* data: `cd1fea243dd0cd56d6116aa65ca53c40478ae9f1b509c2b01c12d4768060ebe9`;
* report: `0946759e18de93d13ae4729b95026b6cb92c59e3c48c0f2ef33c61032b8ebd3e`.

Put d=n−1. The projective closure is a bounded polytope, simple at edges;
its vertex figures are simple. For a generic linear function, the exact
vertex-figure reduction gives I=sum E_w and2f2=sum E_w^0, where E_w^0
counts uncut edges of the vertex figure. These are the identities in
[Khovanskii, Theorem9](https://www.math.utoronto.ca/askold/1986-Faa-1-english.pdf).
They count original two-face corners, including ideal ones.

The lower bound E_w^0>=sum h_j min(j,d−j) is justified vertex by vertex.
Below the cutting hyperplane, all downward edges at a vertex are uncut
and are assigned to their upper endpoint. Above it, all upward edges
are uncut and assigned to their lower endpoint. Every uncut edge is
counted once. The index formula for a simple polytope gives exactly h_j
vertices of index j, independent of the chosen generic linear function.
No assertion that an ideal vertex itself is simple is needed.

At the distinct global minimum and maximum of the original polytope,
the induced cuts miss their entire links. Thus E_w^0=E_w at both.
The h-vector symmetry gives
E_w−sum h_j min(j,d−j)=sum_(j>d/2)h_j(2j−d)>=d,
because h_d=1 and all coefficients are nonnegative. Each endpoint
therefore contributes at least(5/2)d beyond its local lower bound.
The claimed additive5(n−1) is valid even if one or both endpoints are
ideal. It is additional to, not already included in, the local kappa terms.

For a finite vertex the link is a simplex, giving beta_n exactly. For
an ideal vertex the link is a product of simplices; its H-polynomial is
the product of the individual H-polynomials, its vertex count is H(1),
and its edge count is dH(1)/2. Substitution proves the full displayed
inequality. The probabilistic expression is the normalized coefficient
convolution of independent discrete uniforms; the variance relaxation
uses Cauchy–Schwarz with the stated inequality direction.

The exact negative margin at ranks(30,30,30,30,30,13) is scoped correctly:
it rules out a particular per-cusp domination by the coarse distance15
cost. It neither constructs a finite-volume H164 polytope nor contradicts
the existing163 theorem. The positive additive global slack alone cannot
turn a negative per-cusp coefficient into uniform pointwise domination.

The all-rank reasoning, the elementary scalar identities, the primary
reduction input and the scope of the finite corroborating checks have
been reviewed. No mathematical gap was identified. This is an internal
agent-workflow review, not external refereeing or formal verification.
