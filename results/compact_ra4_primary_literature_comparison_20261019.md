# Compact right-angled four-polytopes: primary-source comparison

19 October 2026. Limited literature comparison, not a priority survey.
Notation: f is the number of facets of P, V the vertices of one facet,
d(F)=sum over its faces of (number of sides minus five), and W=sum of
ridge excesses in P. The current manuscript and its frozen predecessor
were not edited for this comparison. No download or scientific job ran.

## Established antecedents

**Kellerhals--Perren (2011).** In the proof of Proposition 3.2(b), printed
p.1311, the outgoing edges from a largest facet's N vertices are shown
to reach N distinct further facets. Thus their argument already contains
f>=N+1 and proves f²>4v. Our perpendicular-edge/opposite-facet mechanism
has this direct antecedent; the repository's f>=2V+5 strengthens the
count by also controlling neighboring-facet extra incidences. Proposition
3.2(d), on the same page, gives the volume formula; Example 2, p.1312,
gives the regular 120-cell's face vector and volume.
[Published primary PDF, pp.1310--1312](https://homeweb.unifr.ch/kellerha/pub/YEUJC1817.pdf#page=13),
[DOI](https://doi.org/10.1016/j.ejc.2011.03.020).

The repository's citation “Proposition 3.5(d)” is the numbering quoted in
Kolpakov's thesis, not that of this published PDF. The article citation
should be corrected to **Proposition 3.2(d)**. This is a bibliographic
correction, not a mathematical change.

**Inoue (2008).** Lemma 3.3, pp.1529--1530, already defines the total
3D excess c=sum(k-5), proves faces-c=12, and recognizes an all-pentagonal
compact right-angled polyhedron as a dodecahedron. Thus d(F)=faces(F)-12
is an established identity. Corollary 3.1, p.1528, states 3D uniqueness;
Section 4, p.1531, describes the Löbell family. These are antecedents or
inputs, not repository discoveries.
[Primary PDF](https://msp.org/agt/2008/8-3/agt-v8-n3-p12-p.pdf#page=7),
[DOI](https://doi.org/10.2140/agt.2008.8.1523).

**Potyagailo--Vinberg (2005).** Proposition 4.1, p.69, supplies relevant
acute-polytope intersection facts. Their compact-dimension discussion
and minimum-facet question are on pp.65--66. This is geometric background,
not a verified antecedent for the exact strengthened inequalities below.
[Primary PDF](https://ems.press/content/serial-article-files/43032?nt=1#page=7),
[DOI](https://doi.org/10.4171/CMH/4).

## Kolpakov thesis, Section 4.4

Section 4.4, pp.73--75, frames the compact 120-cell minimum problems.
Lemma 6, p.74, repeats the volume formula and derives v>=5f; its corollary
settles polytopes with more than 120 facets. Question 2, p.75, proposes
that even one nondodecahedral facet forces more than 120 facets and 600
vertices. Our results do **not** resolve that stronger assertion.
[Thesis, Section 4.4](https://homeweb.unifr.ch/kellerha/pub/KolpakovThesis_Final.pdf).

Retrieval limit: the direct PDF text fetch timed out. These locations were
checked through indexed excerpts from that primary URL, including the
printed page numbers and bibliography p.82; this is not a claim to have
reread the complete thesis. The cited published papers were inspected
directly.

## What the present comparison supports

| Repository statement | Relation to the inspected literature |
|---|---|
| f>=2V+5 | Strengthens the published opposite-facet count; no exact prior statement located in the sources inspected. |
| 2W>=3d(F), hence W>=ceil(3(V-20)/4) | Uses the established 3D excess identity, the paper's proved polygon capacity, and a neighboring-facet sum. No exact antecedent located here. |
| W=0 identifies the whole regular 120-cell | Uses established 3D dodecahedron recognition and metric uniqueness, followed by an explicit whole-facet propagation proof. No claim that the recognition principle is unprecedented. |
| W=0 or W>=16 | Repository geometric/graph proof, with written independent reviews; no positive-W gap statement located in these inspected sources. It is not a 4D minimum-volume theorem. |

The limited negative findings in this table do not establish worldwide
novelty. In particular, already-known formulas, Löbell constructions,
3D classification/uniqueness and the opposite-facet mechanism should
receive credit even though the strengthened arguments are written out.

Egorov's 2026 Theorem 1.1 assumes a common angle strictly below pi/2 and
explicitly excludes the right-angled endpoint; it supplies no compact
right-angled minimum theorem for this comparison.
[Primary theorem and endpoint discussion](https://arxiv.org/html/2609.04258v1#S1.Thmtheorem1).

The current local facet cutoffs 46/46/96 remain necessary conditions.
Neither the strengthened excess gap nor these cutoffs prove f+W>=120,
the minimum-volume assertion, or a second-volume classification.
