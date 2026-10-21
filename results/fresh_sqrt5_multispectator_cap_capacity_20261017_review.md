# Independent review of the two-spectator cap capacities

**PASS.** The [independent checker](../scripts/fresh_sqrt5_multispectator_cap_capacity_20261017_review.py) and [saved certificate](fresh_sqrt5_multispectator_cap_capacity_20261017_review.json) pass generation and a fresh saved replay. They bind the frozen [producer source](../scripts/fresh_sqrt5_multispectator_cap_capacity_20261017.py) and [data](fresh_sqrt5_multispectator_cap_capacity_20261017.json), with SHA256 respectively
`c0438f6942324cc54997d80394059d433c33c34f0ef819c1b1487dbedd9cdcca` and
`e7482aa690f1e9e01663de6e967f487e128a13c740306be2e37774aea558a092`.

The arithmetic review imports no producer or project arithmetic. SymPy independently computes the exact sign and the monic minimal polynomial over the rationals of each squared pair product, for all eleven distinct values of \(\tau\) and all four pair classes: **44 tests**. All signs, integrality outcomes, and permitted-pattern labels agree. This includes cases where the pair square has degree four, and cases where it drops to degree two or one. Integrality is tested on the independently obtained minimal polynomial, rather than by reproducing the producer's relative-norm algorithm.

The normalization and counting received a separate source/proof read. With all roots normalized to square two, a fixed full colored core orientation has

\[
b_{q,p}=v-\frac{\sqrt\tau}{2}q-\frac12p,
\qquad q\ne p,\qquad v^2=h=\frac{3-\tau}{2}.
\]

Here \(q,p\) are distinct isolated \(A_1\) roots of the same actual containing simplex star. The full star spans the ambient vector space, so fixed contacts determine the cap normal. Fixing the actual base-vertex subset \(B\) also removes multiple counting of the same base incidence. This uses the inherited selected dotted maximal-prism scope; it does not extend the statement to arbitrary product diagrams.

For two ordered spectator choices, disjoint supports give product \(h\), a common \(p\) gives \(h+1/2\), one crossed spectator gives \(h+\sqrt\tau/2\), and swapped choices give \(h+\sqrt\tau\). A common \(q\) with different \(p\)'s gives \(h+\tau/2=3/2>0\), hence is always impossible. These are necessary acuteness and arithmetic-integrality gates only: a permitted pair is not asserted to be geometrically realizable.

Represent \((q,p)\) by an arrow \(q\to p\) on the \(k\) isolated copies. The common-\(q\) exclusion gives outdegree at most one. The five observed patterns have the following all-\(k\) bounds.

| Permitted pair types | Maximum number of arrows |
|---|---:|
| None | \(1\) |
| Disjoint only | \(\lfloor k/2\rfloor\) |
| Disjoint and swapped | \(2\lfloor k/2\rfloor\) |
| Crossed chain only | \(3\) for \(k\ge3\), and \(1\) for \(k=2\) |
| Common head only | \(k-1\) |

For disjoint-only pairs the underlying edges are a matching. Allowing swapped pairs permits at most both orientations on each matched edge. In the chain-only case every two underlying edges meet, so they form a star or a triangle: a star admits at most two arrows because repeated tails and repeated heads are forbidden, and a triangle admits at most its directed three-cycle. In the last case all arrows have the same head and distinct other endpoints. These arguments are uniform in \(k\). Independent exhaustive maximum-clique checks for \(2\le k\le7\) give **30 corroborating finite checks**; they are not the justification for unbounded \(k\).

The producer groups by the simultaneous orbit of the entire core contact vector and the actual subset \(B\), with \(\tau\) fixed. The cap for one literal orientation is multiplied by the number of distinct such orientations, not by the automorphism-group order. All identical isolated spectator copies are already included in \(k\). The tree-gauge scale ratios correctly reconstruct normalized multi-contact projection norms. The independent checker also verifies the resulting capacity identities on all **555 occurrence records and 555 columns**, but does not reconstruct their physical projections or raw graph maps.

No source-level gap was identified. The review is deliberately narrower than a second reconstruction of the 555 bases: raw maps, projection matrices, automorphism lists, and provenance grouping received a source read, while the 44 pair tests and small directed-graph checks were independently recomputed. No optimizer, repaired witness, physical-existence statement, or dimension bound is included.
