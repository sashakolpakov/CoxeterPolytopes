# A degree-two positive facet forces excess at least twenty

19 October 2026. Written corollary, no execution. This note preserves the
reviewed W=16 exclusion source and the primary manuscript unchanged.
Use W,d(F),G,H as in the reviewed compact RA4 proofs.

**Theorem.** If G has a vertex of unweighted degree two, then W>=20.

The inherited exact degree-two interface makes its facet F a Loebell
L_k, k>=6, with excess2(k-5), two k-gonal bases, and2k lateral
pentagons. The distinct-facet excess sum gives W>=5k-15. If k>=7
this is already at least twenty. It remains to treat k=6.

Let B_1,B_2 be the two base neighbors, and let A_(j,1),...,A_(j,6)
be the six lateral neighbors meeting base B_j in F. These twelve
lateral neighbors, the two base neighbors and F are fifteen distinct
facets: a convex pair of facets cannot share two different ridges.
Each has excess at least two, giving baseline thirty in sum d=2W.

The full metric lemma and its parent-facet incidence justification are
in `compact_ra4_sixteen_excess_equality_exclusion_research_20261019.md`,
SHA256
`8caf24e085b204375e3a36b108513e55d8fb0cc66559b44f3e3caf1bd17977e8`,
Sections2--3. Every lateral pentagon in L6 has cyclic lengths
(b,a,c,c,a) with cosh(b)=2, cosh(a)=sqrt(2), cosh(c)=sqrt(3).
Its unique longest edge lies on its hexagonal base. Consequently,
if d(A_(j,i))=2, that facet is also L6 and its unique longest edge
on F intersect A_(j,i) must lie on the same parent B_j. Compact
simplicity identifies the third parent facet along that edge uniquely.
Thus A_(j,i) intersect B_j is a hexagon, giving a unit G-edge.

Put m_j=number of the six A_(j,i) with excess greater than two.
Their total excess ABOVE baseline is at least m_j. The remaining
6-m_j facets give that many distinct unit G-edges at B_j, in
addition to F--B_j. Hence

    d(B_j)>=max(2,7-m_j),
    d(B_j)-2>=max(0,5-m_j).

The total extra excess from this row and its base is therefore at least

    m_j+max(0,5-m_j)>=5.

The two rows and bases are disjoint selected facet sets, so these
extras add, regardless of other positive facets elsewhere. Including
F's baseline gives

    2W=sum d >=30+5+5=40.

Thus W>=20. The proof does not suppose that all lateral neighbors
have excess two, that every positive facet was selected, or that their
ridge contributions are disjoint. It sums distinct FACET excesses in
the exact complete sum, which already counts each ridge twice.

The only metric inputs are those of the reviewed L6 lemma: compact
3D Andreev uniqueness realizing the explicit reflection and the
elementary pentagon identities. The result is not a minimum-volume
or four-dimensional classification theorem. No TeX/source integration,
new count threshold, scientific run, formalization or novelty claim is
made here.
