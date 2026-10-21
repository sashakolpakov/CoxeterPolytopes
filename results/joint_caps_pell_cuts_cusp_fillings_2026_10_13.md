# New-dimension search: joint caps, Pell cuts, and cusp fillings

**No new compact or finite-volume Coxeter dimension has been established.**
This round makes new exact constructions in the H9 and H20 searches, but
every tested chamber still fails the volume or angle requirement. Four
separate reports contain the constructions, precisely bounded conclusions,
and independent verification programs.

## Arithmetic mirrors beyond uniform enumeration

The [arithmetic continuation](fresh_h9_joint_arithmetic_20261013.md) uses
the previously introduced Lorentz lattice
`⟨−T⟩ + B² + A4 + A1`, where `T³−3T−1=0`, `T=2cos(π/9)`, and
`B=[[2,−T],[-T,2]]`.

It first frees all four A4 covectors and the scalar simultaneously, keeping
the previous height bound 1000. This completes the missing directions in
the forced H5 face and yields a compatible 197-wall H9 cone. Fourteen of
the new roots cut the preceding 22-wall cone's saved escape. A different
exact escape remains, very close to the light cone.

Following that gap reduces the root condition to the relative Pell equation

```
30t² − T x² = 5.
x = 8715 + 30785T + 16380T²,
t = 2181 + 7705T + 4100T².
```

The resulting norm-five mirror and its A4 reversal are integral reflective
roots, compatible with every existing wall. Their height is approximately
5.8 billion, far outside the uniform enumeration. They cut the saved gap.
Two further relative norm constructions, the last selected after a bounded
40-direction perturbation scout, produce another four compatible mirrors.
The largest coefficients have 926 decimal digits.

All six additional mirrors have exact reflection, norm, orientation and
angle certificates. Nevertheless the final **203-wall H9 cone and its
197-wall H5 face have exact infinite-volume escapes**. This is progress in
targeted root construction, not a proof of reflectivity or a new polytope.
Independent Fraction replays verify the complete bounded face census and
all three relative norm orbits. A separate agent reran the final 203-wall
replay successfully. The full arithmetic host remains open.

## New caps that repair several gaps simultaneously

The [joint-cap search](fresh_h9_joint_completion_20261013.md) starts from
one actual 32-facet configuration of the released H8 seed. It extracts
24 exact separated escapes and constructs caps that cut at least two while
remaining compatible with at least 18 of the configuration's 22 old caps.

There are **1,756 distinct new mirrors outside the preceding 1,084-entry
dictionary**, represented among 1,781 accepted entries. Each cuts two to
eight saved gaps; 40 entries preserve all 22 old caps. Yet one exact escape
survives the entire new bank, excluding a finite-volume completion from
this starting configuration and these mirrors.

Six deliberate repairs of that gap replace between one and six old caps.
For each fixed replacement core, all eligible new mirrors are considered.
Two other exact gaps survive in every case. These are exclusions of the
specified extension schemes; they do not close the combined old and new
dictionary with all walls optional. A further, restrictive quartic-anchor
probe produces no accepted cap and is recorded only as numerical scouting.

The producer-free replay verifies all new roots, the exact mirror census,
39,182 old-cap relations, 42,744 target signs, and all six replacement
coverage certificates. Parallel-pair rejection is explicitly treated as
a compactness restriction where used.

## Cusp filling of known H9 chambers

The [cusp-filling search](fresh_h9_cusp_filling_20261013.md) changes angles
of two known 18-facet noncompact H9 chambers, with optional facet deletion.
It gives two continuous obstructions, both allowing arbitrary additional
acute-compatible walls and requiring a strict interior:

- For the first chamber, deleting at most five old facets and changing
  only old parallel pairings cannot yield compactness. One of 54 affine
  subdiagrams survives every such deletion, forcing an ideal ray.
- For the second, deleting at most five old facets and changing any old
  nonzero pairings while preserving old orthogonalities cannot yield
  compactness. At least eleven leaf normals would have to fit as
  nonsingular orthogonal blocks into Lorentz rank ten. Any singular pair
  instead gives an ideal ray.

The exact replay checks every one of the 8,568 five-deletion cases, positive
affine kernels and coranks, and the orthogonal matching structure. Another
agent independently reran the verifier and reviewed the geometric argument.

A separate numerical search allows changes to up to four π/3 edges or one
right-angle edge. Its two overlapping passes examine 171,272 and 1,383,498
basis mutations, without producing a rank-ten candidate retaining at least
14 original walls. These numerical outcomes are not general impossibility
proofs. An independent source review confirmed the reported search scope.

## Irrational cusp directions and a continuous H20 obstruction

The [H20 section search](fresh_h20_joint_cusp_20261013.md) first considers
an entire affine-A2 normal span at a cusp of the Borcherds H21 chamber.
With nonzero mixed coefficients, every such section retains a complete
Coxeter cusp link, but an actual ridge has an angle strictly between π/4
and π/3. Positive exact face witnesses establish incidence for the whole
two-parameter family. The earlier two-wall theorem covers zero
coefficients. Together these exclude every strict interior section in
the displayed span, including noncentral sections.

A separate construction uses two cusp directions over `Q(√2)` that leave
an affine G2 triangle. Twelve exact global angle equations give twelve
finite-volume noncompact H20 sections with complete Coxeter cusp links.
Eight also have the intended actual right ridge. **All twelve have another
actual forbidden ridge.** For example,

```
E = 2√2 r34 + (3−√2) r46 − √2 r204 + 2√2 r210,
E² = 54,
```

gives a right ridge at walls 2/46, but squared cosine 27/52 at the actual
ridge 1/169. Thus it is not a Coxeter polytope. Exact redundancy identities,
positive affine kernels, strict interiors and ridge witnesses accompany
every section. Both the continuous obstruction and twelve constructions
pass separate standard-library exact replays with corrupted-certificate
checks. These claims do not exclude other cusp directions or translations.

The new mirrors and sections are explicitly recorded constructions. Their
novelty relative to the repository's preceding searches does not establish
worldwide priority. The dimensional objective remains unresolved.
