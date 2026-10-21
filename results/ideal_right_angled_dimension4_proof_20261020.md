# Maximum cusp height and the ideal right-angled dimension-four bound

20 October2026. Complete written proof, submitted for independent review.
This note supersedes the citation numbering in the preserved candidate
`ideal_right_angled_maximum_cusp_height_kissing_candidate_20261020.md`
(SHA256 `95ef1ec3edef6083992aad06dd56e2608b1f120ec098a0d197a2320b55823537`).
No scientific computation, formal proof compilation, or TeX build was run.
No status table or previously frozen source was modified.

## Statement and external theorems

**Theorem.** An entirely ideal finite-volume right-angled hyperbolic
polytope has dimension at most four. The bound is attained by the
regular ideal right-angled24-cell.

Two published theorems are isolated as inputs:

* Musin's Main theorem, printed page4: a finite set in `S^3` with
  pairwise products at most `1/2` has at most24 members.
  [O. R. Musin, The kissing number in four dimensions, Annals of
  Mathematics168 (2008),1–32](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf).
  Only cardinality is used, not uniqueness of a kissing configuration.
* Kolpakov's **Theorem4**, printed page1945: an entirely ideal
  right-angled four-polytope has at least24 facets, with equality only
  for the regular ideal right-angled24-cell.
  [A. Kolpakov, On the optimality of the ideal right-angled24-cell,
  Algebraic & Geometric Topology12 (2012),1941–1960](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf).
  Theorem2 is the preprint numbering; the linked published version
  uses Theorem4. Its final equality argument is on printed page1958.

We also use standard Coxeter-polytope foundations. Distinct actual
outward unit normals of an acute-angled hyperbolic polytope have
nonpositive products, including nonadjacent facets. This is the
all-wall condition in Vinberg, **printed page42, §2.1** of
[Hyperbolic reflection groups, Russian Mathematical Surveys40:1
(1985),31–75](https://doi.org/10.1070/RM1985v040n01ABEH003527).
Vinberg explicitly notes Andreev's result that nonobtuse actual ridge
angles imply that nonadjacent supporting hyperplanes do not intersect.
Thus the condition applies to a right-angled polytope from its actual
dihedral angles; it is not an extra assumption about invisible angles.
Faces of dimension at least three inherit finite volume and right
angles. The facet-adjacency graph is connected. For a finite-volume
ideal polytope this is the ordinary facet graph of its compact
projective convex closure, whose vertices are the finitely many ideal
vertices. Positive-dimensional ridges have actual hyperbolic interior.

## Proof in dimension five

Assume `P` is an entirely ideal right-angled finite-volume five-polytope.
Use Lorentz signature `(-,+,+,+,+,+)`, outward unit normals `m`, and
halfspaces `<x,m> <= 0`. Fix a future null representative `v` of an
actual ideal vertex and set `H(m)=-<v,m> >= 0`. Incident facets have
height0, and some other facet has positive height. Choose a facet
normal `n` of maximum height `H0>0`, and put

    q = v/H0 + n.

Then `q^2=-1`, `<q,n>=0`, and `<v,q>=-H0<0`, so `q` is future.
For every actual normal `m != n`, all-wall acuteness gives

    <q,m> = -H(m)/H0 + <n,m> <= 0.

Hence `q` is an actual point of the four-facet `Q=P intersect n-perp`.
It may initially be on its boundary; no strictness is assumed here.

Let `m_i` run over the actual ridge neighbors of `n`. These give
exactly the facets of `Q`, and are intrinsic unit normals since
`<m_i,n>=0`. In the Lorentz five-space `n-perp`, use `q` as the time
axis and write

    m_i=(t_i,sqrt(1+t_i^2)*u_i),
    t_i=-<q,m_i>=H(m_i)/H0 in[0,1],    u_i in S^3.

For distinct neighbors, `<m_i,m_j> <= 0` implies

    u_i.u_j <= t_i*t_j/sqrt((1+t_i^2)*(1+t_j^2)) <= 1/2.

Equal directions would instead give positive product for two unit
normals, so the directions are distinct. Musin's theorem bounds the
number of facets of `Q` by24. Face inheritance and Kolpakov's theorem
give the reverse bound and make `Q` the regular ideal24-cell.

We now use the metric equality case. Its standard unit facet normals
in signature `(-,+,+,+,+)` are

    (1,epsilon*e_i+delta*e_j),
    1<=i<j<=4,    epsilon,delta in{+1,-1}.

The center is `(1,0)`, and these normals occur in opposite spatial
pairs `(1,r),(1,-r)`. Write the transported point `q` as `(a,w)`,
where `a>0` and `a^2-|w|^2=1`. Its heights on an opposite pair are
`a-r.w` and `a+r.w`. Both are at most1, hence `a<=1`. Unit future
timelikeness gives `a>=1`, forcing `a=1,w=0`. Thus `q` is the center,
every `t_i=1`, and every facet adjacent to `n` has height `H0`.

Repeat the same argument at every maximum-height facet. The nonempty
set of such facets is closed under adjacency, so connectivity makes
it the entire facet set. This contradicts the height0 facets at the
chosen cusp. There is no such `P`.

For any entirely ideal right-angled polytope of higher dimension, a
five-face would inherit entire ideality, finite volume, and right
angles. Thus none exist above dimension four either. The regular
ideal24-cell gives equality.

## Scope and reproducibility

The proof is independent of arithmeticity, field, facet spectrum,
face-number estimates, search completeness and numerical optimization.
It uses Musin's cardinality theorem and Kolpakov's geometric equality
case exactly as stated above. It does not assert uniqueness of a
24-point kissing arrangement or derive finite volume from local data.
Every projection is checked against all actual walls. The canonical
unit-normal scaling is retained throughout, so the maximum and all
height ratios refer to the same normalization.

The standalone human-readable source is
`manuscripts/ideal_right_angled_dimension4.tex`. At this handoff it is
unbuilt, and no Lean verification is asserted. Subsequent reviews and
executions must be recorded separately and must not be inferred from
this source-only proof note.
