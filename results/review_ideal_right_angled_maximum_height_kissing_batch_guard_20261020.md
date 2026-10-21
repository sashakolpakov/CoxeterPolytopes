# Independent adversarial review: maximum cusp height and the ideal ceiling

20 October2026. **Written mathematical review PASS.** No computation,
scientific check, Lean compilation or build was run. This is an independent
read and derivation of the whole argument, not an execution certificate.

Reviewed source:
`ideal_right_angled_maximum_cusp_height_kissing_candidate_20261020.md`,
SHA256 `95ef1ec3edef6083992aad06dd56e2608b1f120ec098a0d197a2320b55823537`.
Its bytes were read completely and remain unchanged. The conclusion below
uses the published external theorems identified in Section1. No claim of
literature priority or exhaustive novelty checking is made.

## 1. External theorem statements independently checked

I opened both primary journal PDFs during this review.

* Musin, *The kissing number in four dimensions*, Annals168(2008),
  printed p4, Main theorem: the four-dimensional kissing number is24.
  Equivalently a finite set of unit vectors in R4 with distinct-pair
  product at most1/2 has cardinality at most24. Its formulation as a
  spherical pi/3-code is explicit in the proof.
  [Primary PDF](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf).
* Kolpakov, *On the optimality of the ideal right-angled24-cell*,
  AGT12(2012), printed p1945, **Theorem4**, states both the minimum
  of24 facets among ideal right-angled H4 polytopes and metric uniqueness
  at equality. [Primary PDF](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf).

The candidate calls the second result Theorem2. That is a citation-number
correction, not a gap in the stated input: use journal Theorem4 in a
successor proof. No uniqueness theorem for spherical kissing codes is
needed. Only Kolpakov's metric equality theorem is used.

## 2. Height signs and all-wall feasibility

For an actual ideal point v of P, choose its future null representative.
Every outward unit facet normal m satisfies v.m<=0 by continuity of its
supporting inequality. Hence H(m)=-v.m>=0; precisely the cusp-incident
facets have height0. Not every facet contains v: a full-dimensional
finite-volume polytope cannot be the intersection solely of vertical
halfspaces at one ideal point. In upper halfspace such an intersection
would have infinite volume toward the finite boundary. Thus some height
is strictly positive. Finiteness of the actual facet set supplies Hmax.

The ALL-WALL sign n.m<=0 for distinct facets is essential and is valid
for an acute/Coxeter chamber; it is not being inferred from the products
only at ridge neighbors. One can also check the right-angled case from
the standard Coxeter chamber theorem: two nonadjacent side reflections
generate an infinite dihedral standard parabolic subgroup. Their support
hyperplanes therefore cannot intersect (an irrational elliptic product
would violate discreteness; a rational one would give finite order).
For disjoint hyperplanes their outward orientations toward an actual
nonredundant chamber give n.m<=-1. Adjacent right-angled facets give0.
This is the standard all-wall acuteness input, with its orientation fixed.
It is also the Gram sign convention explained in Vinberg, *Hyperbolic
reflection groups*, Russian Math. Surveys40:1(1985), Section2.1
([primary PDF](https://www.mathnet.ru/links/df899e7a8aeed4a3951f34fd2ab8427c/rm2140_eng.pdf)).
No exact printed page is asserted here.

With n maximal and H0=H(n)>0, independently expanding

    q=v/H0+n

gives q.q=-1 and q.n=0. Furthermore q.v=-H0<0, so q is future.
For every other wall

    q.m=-H(m)/H0+n.m<=0.

Consequently q is in the ACTUAL facet Q=P intersect n-perp. A positive
maximum, rather than a minimum, is the correct choice for the later
ratios. A height-zero neighbor merely puts q on the boundary of Q;
this causes no problem below and is not silently excluded.

## 3. Intrinsic facets, radial directions and cardinality

Every facet of Q is a ridge of P lying on n, hence is supported by an
actual ridge neighbor m_i of n. These are in bijection with the facets
of Q. The finite-volume right-angled face correspondence ensures there
are no extra intrinsic facets from projecting unrelated ambient walls.
Each m_i.n=0, so its ambient unit normal is already the intrinsic unit
normal in the Lorentz five-space n-perp.

Idealness is inherited: a vertex of a facet is a vertex of the whole
polytope, and an ambient ideal point on n-perp is intrinsically ideal.
The facet is the convex hull of its finitely many ideal vertices and
has finite intrinsic volume, for example by a triangulation into
finite-volume ideal4-simplices. Intrinsic right angles follow from the
orthogonal ambient facet normals. Thus Kolpakov's hypotheses apply to Q.

With q as time origin, write

    m_i=(t_i,sqrt(1+t_i^2)*u_i),
    t_i=-q.m_i=H(m_i)/H0, 0<=t_i<=1, u_i in S3.

For every DISTINCT pair, all-wall acuteness yields

    u_i.u_j<=t_i*t_j/sqrt((1+t_i^2)(1+t_j^2))<=1/2.

This is a code in FOUR-dimensional Euclidean space, not a code in
three-dimensional space. Repeated directions are impossible: equal
directions would make the unit Lorentz product positive. Therefore
Musin gives f(Q)<=24, while Kolpakov gives f(Q)>=24 and identifies Q
isometrically with the regular ideal right-angled24-cell. There is no
average-facet substitution or abstract corner completion in this step.

## 4. The equality-center deduction

The standard regular-cell normals (1,r), r=+/-e_i+/-e_j, have square1
and occur in opposite spatial pairs. An arbitrary future unit point
transported from q is (a,w), with a>=1 and a^2-|w|^2=1. The two height
inequalities at (1,r),(1,-r) read

    a-r.w<=1, a+r.w<=1.

Adding gives a<=1, so a=1,w=0. All24 intrinsic t_i are consequently1.
This metric argument works even though q was not assumed interior.
It uses neither a stronger unpublished minimum theorem nor the
unproved uniqueness of a maximal kissing arrangement.

Thus EVERY actual neighbor of ANY globally maximal-height facet has
the same maximum positive height. Every use of the local argument
still meets the same hypotheses; no induction on an unverified
description of neighboring facets is required.

## 5. Global contradiction and higher-dimensional scope

The adjacency graph of actual facets of a convex finite-volume ideal
polytope is connected. This can be seen from the connected1-skeleton
of the dual of its compact projective closure; adjacency is across
actual ridges, exactly as used above. Propagating maximality along a
finite path therefore reaches every facet, including the height-zero
facets at the originally chosen cusp. This contradicts H0>0.

The proof excludes an entirely ideal finite-volume right-angled
H5 polytope without an arithmetic assumption, a facet-spectrum
restriction, or an enumeration hypothesis. Higher dimensions are
excluded by actual face inheritance: a polytope of dimension>=5
has an entirely ideal finite-volume right-angled five-face. The
known ideal24-cell in H4 establishes sharpness of the resulting
ceiling4. This review finds no missing local-to-global implication.

## 6. Verification boundary

PASS is for the complete written mathematics with the named external
theorems and standard acute-polytope/face facts. It does not mean
Musin's proof, Kolpakov's classification, the Lorentz reduction or
face inheritance has been formalized in Lean. A conditional kernel
may check the cardinality squeeze and finite graph propagation while
keeping those inputs explicit. No compilation or axiom audit is
claimed in this review, and no sealed theorem record was changed.
