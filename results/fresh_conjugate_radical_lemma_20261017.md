# Conjugate positivity excludes affine principal diagrams

Let K be a number field with two real embeddings, the chosen physical one
and sigma. Let G be a symmetric matrix over K, and suppose sigma(G) is
positive semidefinite. If a principal block G[I] is singular and v is in
its kernel over K, extend v by zeros outside I to obtain vbar. Then

    sigma(vbar)^T sigma(G) sigma(vbar) = 0.

For a real positive semidefinite matrix, a vector of quadratic norm zero
lies in its kernel. Therefore sigma(G) sigma(vbar)=0, and injectivity of
sigma gives G vbar=0. This argument needs no positivity of G[I] at the
physical embedding and no enumeration of principal subsets.

For a physical affine Coxeter principal block, its kernel has a strictly
positive vector on that connected block. This vector can be chosen over K:
row reduction gives a K-vector spanning its one-dimensional kernel, and
its physical sign can be chosen positive. Consequently G has a nonzero
nonnegative kernel vector, extended by zero outside the affine block.

Such a relation is incompatible with a common strict chamber interior X:
if G represents wall normals r_i and sum(v_i r_i)=0 with v_i>=0 and some
v_i>0, pairing with X gives a strictly negative sum equal to zero. In a
nondegenerate ambient Lorentz space with the wall span containing a
Lorentzian subspace, a Gram-kernel relation is an actual vector relation.
The enumerated cores have negative index one and full physical span modulo
their kernel, which is the setting used here.

Thus the existing exact tests already exclude every affine principal block
once they require sigma(G) positive semidefinite and either reject corank
greater than one or explicitly reject a nonnegative generator of the full
kernel. Both conditions hold in the new dotted attachment queue. The old
subset search for affine diagrams is redundant for its accepted states.
Removing that preliminary rejection changes some rejection explanations,
but must leave the retained attachment set unchanged.

This is a linear-algebra lemma used for computation, not a new dimension
bound. It does not apply without the conjugate-positivity hypothesis and
must not be used in the unrestricted nonarithmetic construction searches.
It also does not make the remaining rank and face-payment checks optional.

The [new queue](../scripts/fresh_sqrt5_global_b1_dotted_fast_20261017.py)
starts from the same377 exact parameter seeds and uses the same full
finite-corner domains, physical/conjugate arithmetic, marked isomorphisms
and hereditary weight pruning as its
[predecessor](../scripts/fresh_sqrt5_global_b1_dotted_closure_20261016.py).
The old queue is retained as a stopped partial comparison input. It is not
called complete, and no previously frozen source is modified.
