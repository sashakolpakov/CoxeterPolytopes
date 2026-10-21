# Vertex-weighted capacity for simplex bases of maximal prisms

Let P be a compact simple n-polytope, S a maximal d-simplex face, and
let I(S) be the actual incidences of maximal prisms Q along one of their
simplex bases B contained in S. Use the conventions and scope of the
[shared-base theorem](fresh_prism_simplex_shared_capacity_20261017.md):
a k-prism has a (k-1)-simplex base with k vertices; a square uses all four
edges as bases. Distinct incidences count distinct actual bases, with
their unique maximal prism, rather than overlapping local diagram
descriptions.

For every nonnegative function w on the vertices of S,

\[
 \boxed{\quad
 \sum_{(Q,B)\in I(S)}\ \sum_{v\in V(B)}w(v)
 \ \le\ (n-d)\sum_{v\in V(S)}w(v).
 \quad}                                                   \tag{1}
\]

In particular, every individual vertex of S belongs to at most n-d
such bases, with all prism dimensions counted together. This refines the
unweighted bound by retaining where the bases lie inside S.

## Proof

The shared-base theorem associates each incidence to an immediate
(d+1)-coface U of S. Within a fixed U, bases of distinct maximal prisms
have disjoint vertex sets in S. Otherwise the shared external direction,
triangular base faces and quadrilateral side faces force their union to
belong to one larger prism, contrary to maximality. Its existing proof
constructs the relevant facets and establishes this for actual faces.

There are n-d immediate cofaces U because P is simple along S. At a
fixed vertex v of S, at most one base assigned to a given U can contain
v. Therefore its total base multiplicity is at most n-d. Multiply this
pointwise inequality by w(v)>=0 and sum over all v. This proves(1).
There is no independence assumption between the different vertices.

The identical argument applies to compact S and compact counted prisms
inside a finite-volume Coxeter polytope: the
[finite-volume transfer](fresh_finite_volume_prism_capacity_20261017.md)
supplies the same simple finite vertices and immediate cofaces. It does
not assert a capacity at nonsimple ideal vertices.

## Aggregation by an invariant vertex label

Choose a label ell of a vertex within its maximal simplex star.
Examples are its full finite Coxeter forest, or an automorphism orbit
of vertices in the complete simplex star. The label must be preserved
under the star equivalences used to aggregate incidences. Let

* c(S,ell) count vertices of S with label ell;
* c(B,ell) count the vertices with that label in an incident simplex base B.

Putting w(v)=1 when v has label ell gives

\[
 \sum_{(Q,B)\in I(S)}c(B,\ell)\le(n-d)c(S,\ell).       \tag{2}
\]

For a global finite-face model with simplex type S of frequency f_S
and full prism-provenance frequencies f_rho, let m(rho,S,ell) be the
sum of these vertex counts over all bases of that actual provenance
whose maximal containing simplex has type S. Then

\[
 \sum_\rho m(\rho,S,\ell)f_\rho
 \le(n-d)c(S,\ell)f_S.                              \tag{3}
\]

Here all representatives of type S have the same c(S,ell). If
c(S,ell)=0, every actual coefficient m(rho,S,ell) must be zero;
this is an exact incidence consistency condition. Different labels
must use the same full-provenance frequency f_rho. Taking independent
coordinatewise minima does not supply such a simultaneous split.

Summing(3) over a complete partition of labels recovers the previous
shared-base inequality because sum_ell c(B,ell)=k and
sum_ell c(S,ell)=d+1. The separate labelled inequalities can be stronger:
the summed inequality alone can hide an excess concentrated at one
kind of vertex. No strict improvement is asserted for every domain.

## Scope

This is a direct weighted refinement of the previously reviewed
pointwise geometric proof. It uses no arithmetic classification,
optimization, finite-label cutoff or new computational certificate.
An application still must reconstruct the actual base vertices and
their labels and verify consistent provenance masses. It is not a
dimension theorem and does not itself exclude the current weighting
witness. Predecessor files were not changed.

**FINAL/FROZEN.** The exact_hub_counts agent independently reviewed the
pointwise argument, invariant-label aggregation, square-base convention,
structural zero condition and finite-volume scope. That internal proof
review passed; it is not an arithmetic replay or external refereeing.

