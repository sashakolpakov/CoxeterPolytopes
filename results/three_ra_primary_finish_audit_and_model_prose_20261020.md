# Finish audit of the three built right-angled papers

20 October 2026. Source/editorial completeness read, not a new mathematical
referee report, primary-literature audit, execution or rendering review.
I read all three sources and the relevant finite-model implementations.
No primary source, PDF, frozen certificate or catalogue entry was edited.

|Paper|Source SHA256 read|Necessary finishing action|
|---|---|---|
|`right_angled_sharp_dimension_bounds.tex`|`3c3ddd8540593004bcd5fb57832afc77ffc639268e66cb0ab29302b286654310`|None for the stated written proof. Do not add an algorithm section merely because supporting checks exist.|
|`ideal_ra4_volume_and_structure.tex`|`e503ad2cb8edafd827f735a62568a3372407e9dfdaf3f50d09b05fb3bc75f025`|Replace the one-line active-four enumeration assertion with its finite algorithm and completeness proof below; clarify the already recorded supporting checks in the existing appendix.|
|`compact_ra4_volume_constraints.tex`|`afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`|None for its declared capacity/excess statements. Its finite graph cases are proved in prose, so no computational census specification is missing.|

## 1. Sharp dimensions: already a finished written argument

The shared actual maximum-facet projection, strict height iteration,
height-to-code conversion, direct tensor-square identity and strict
27-point equality exclusion are all displayed. The entirely ideal
branch separately states the precise Musin cardinality and Kolpakov
metric-equality inputs and proves the opposite-normal obstruction. Both
higher-dimensional tails use the stated inherited right-angled faces;
the general Coxeter-facet issue is not imported into this RA setting.
All numerical calculations used by the proof are on the page. No script
has to establish a bounded-domain exhaustion for this argument.

The paper already has a human-readable common method followed by the
explicit general/entirely-ideal split, one useful native figure and a
short bibliography. Do not expand it with old mixed10/11 or ideal5/6
predecessors, proof-session histories, or an unnecessary computer-proof
section. Optional wording only: at line319, ``a four-facet Q'' may be
made ``a four-dimensional facet Q'' to avoid an unfamiliar abbreviation.
This is not needed for the mathematical logic.

## 2. Ideal volume and chains: only a finite-model exposition gap

The published minimum, new second-volume argument and restricted chain
classification are separated correctly. The equality-spectrum reduction
precedes the chain theorem; the paper does not infer its facet hypothesis
from the desired decomposition. The actual cusp surfaces, wall normal
identities, feasible cut, termination of cutting, cell recognition and
opposition are written arguments. The face-count formulas and attainment
are included. There is no need to incorporate the later all-antiprism
or ten-vertex extensions in order to finish these statements.

The one place still phrased as a computational assertion is
lines960--963: enumeration of10,626 active four-wall subsets. The
surrounding text gives the walls and boundedness, but should also say
why this finite enumeration contains every vertex, including a
nonsimple one. Lines980--996 give face-reconstruction outputs without
specifying the terminating facet-intersection closure. The short
staged replacement below supplies both facts within the declared model.

The appendix at lines1000--1028 already correctly limits the finite
checks and distinguishes them from universal geometry. The second
staged paragraph makes its computed-proposition status explicit without
moving implementation details into the geometric proof.

### 2.1 Proposed TeX: replace the last sentence at lines961--963

This text belongs to the existing exact-attainment subsection. It uses
no new theorem about arbitrary polytopes and introduces no new model.

```tex
Here is a finite exact verification of the vertex list. For each of
the $\binom{24}{4}=10,626$ four-element subsets of the walls, form
the $4\times4$ matrix of their spatial coefficients. If its rank is
four, solve the four affine equalities over $\mathbb Q$ and retain
the solution exactly when it satisfies all $24$ inequalities.
Discard repeated solutions. Singular subsets contribute no candidate.
This procedure terminates because its loops and all its rational
linear systems have fixed finite size.

It is complete for this intersection. At a vertex of a full-dimensional
polytope, the active spatial wall coefficients span $\mathbb R^4$:
otherwise a nonzero vector orthogonal to them permits sufficiently
small perturbations in both directions while every inactive inequality
remains strict, contradicting extremality. Some four active walls
therefore have rank four and occur in the enumeration. Conversely a
feasible solution of four independent supporting equalities is a
vertex, since any segment through it in the intersection would have
both endpoints in those four equalities. The intersection here is
bounded, as follows directly from $|x_i|+|x_j|\le1$, and contains
the origin strictly. The exact calculation returns precisely the
eight signed coordinate vertices and the sixteen half-sign vertices
displayed above. Each has Euclidean norm one, so all are ideal.
```

### 2.2 Proposed TeX: after the double-containment proof, before its table

```tex
The double needs no assumption that a saved vertex list is exhaustive.
The preceding halfspace argument identifies it with the convex union
of the two cells. In the Klein model that union is the convex hull
of the two finite vertex sets. Reflection and normalization give
$42$ distinct points on the ideal sphere, with six shared points.
Strict convexity of the sphere makes every one of these points an
extreme point of their hull. Hence this is the complete vertex set.

To reconstruct faces, first record, for each wall, the subset of these
vertices on its equality hyperplane. Start with the full vertex set
and repeatedly intersect known sets with these facet sets, discarding
the empty set and repetitions, until no new set appears. There are
at most $2^{42}$ subsets, so the procedure halts. Every generated set
is the vertex set of a face, and every face of a polytope is an
intersection of its containing facets. The closure is therefore
complete. Its dimension is the affine rank of its vertices, computed
over $\mathbb Q$. The resulting exact incidence counts are those
in the following table; they are model checks, not a classification
of other ideal right-angled four-polytopes.
```

The number $2^{42}$ is a termination bound, not a claim that the actual
implementation enumerates all subsets or that such a search was run.
The actual script uses the much smaller intersection closure. All
vertices in the double are already known from the written convex-union
argument, so this paragraph does not falsely claim that the double's
38-wall active bases were exhaustively re-enumerated.

### 2.3 Proposed TeX: short addition to the existing verification appendix

```tex
The finite checks use only explicitly bounded sets. The original
$24$-cell test is the rational active-wall algorithm just proved;
the face test is the terminating facet-intersection closure. For the
cuboctahedral model, rescale the spatial coordinates so the Lorentz
form is $-x_0^2+(x_1^2+x_2^2+x_3^2)/2$. Its normals are then
$(1,\pm2e_i)$ and $(1,\pm1,\pm1,\pm1)$. The same vertex argument
runs over the $\binom{14}{3}=364$ triples of these fourteen rational
wall equations. For the six-cusp graph check, a
four-regular graph on six vertices has a perfect matching as its
complement. Enumerate the $\binom{15}{3}$ triples of possible missing
edges and retain those covering each vertex once. These are exactly
the fifteen labelled matchings; each complement has the eight
triangles and cube dual described in the proof. The cube-normal
check is rational row reduction of its six square equations in
eight unknowns, with kernel spanned by the constant and three sign
characters. These procedures halt by their explicit finite bounds.
Their completeness concerns these specified models and graph classes,
not the cited classification of three-polyhedra.

The existing independent saved-data checks return the two face vectors
$(24,96,96,24)$ and $(42,168,164,38)$, the cuboctahedron's twelve
cusps with squared aspect ratio two, the fifteen labelled six-cusp
graphs and cube-equation rank four. The separate $D_4$ check examines
all $\binom{24}{2}=276$ unordered root pairs: their spatial products
$-2,-1,0,1$ occur respectively $12,96,72,96$ times. It also checks
the $96$ common-neighbor cases and their $1,152$ reflected
octahedral-vertex incidences. These are the recorded outcomes of
the fixed checks. The universal injection, surfaces, cuts and gluing
remain the written arguments above.
```

The model implementations read for this staging were
`fresh_ideal_ra4_double24_volume_20261018.py`,
`check_ideal_ra4_second_volume_saved_20261018.py`,
`fresh_ideal_ra4_equality_structure_20261018.py`,
`check_ideal_ra4_equality_structure_saved_20261018.py`, and
`check_ideal_ra4_octahedral_path_roots_20261018.py`.
The recorded outcomes are in the existing source-bound core ledger
`ideal_ra4_unified_core_verification_final_20261018.json`. None was
executed for this read. The published small-facet lists remain explicitly
external inputs; the planar Figure7 transcription check does not
establish uniqueness among arbitrary three-polyhedra.

## 3. Compact volume constraints: no algorithm gap to fill

The paper proves its Lorentz gap, perpendicular opposites, polygon
closing, capacity, neighbor excess, zero recognition and graph cases in
prose. In particular the nine-vertex complement lemma, five-neighbor
marked-polygon exclusion and ten-vertex complement lemma are complete
case arguments on the page, not appeals to a graph enumeration script.
The finite three-dimensional metric inputs are stated separately.
The final $W=0$ or $W\ge17$ proof and the degree-two implication
$W\ge20$ therefore need no algorithm specification.

Its scope is already explicit at lines98--104 and1123--1144: no proof
of $f+W\ge120$, no universal minimum or second compact volume, and no
unrestricted decomposition into regular cells. Keep those limitations.
Later f>=65, C26/C28 and neighbor-block results are supporting successors,
not finishing requirements for the current claims; adding them now would
change scope and enlarge the paper. The four native figures illustrate
the stated geometric or graph proofs and are described as such.

## Editorial conclusion

Keep sharpRA and compactRA4 at their present mathematical scope; neither
needs expansion for its own sake. The ideal4 paper needs only the short
finite-model algorithm and completeness paragraphs above to satisfy the
requested computational exposition standard. All three already have
sensible technique-based organization. Source reviews, build and Lean
status should stay in the reproduction maps; no new execution, PDF
page count, warning count or geometric formalization status is asserted
by this editorial read.
