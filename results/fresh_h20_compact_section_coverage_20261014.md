# Compact-section coverage audit and a new coupled construction

**The existing results do not exclude all compact H20 sections of the Borcherds H21 chamber.** The contact-support bound of three is correct, and was already proved in the [all-right-edge report](fresh_h20_allright_edge_gate_20261013.md). It is a local reduction, not a completed classification of its remaining cases.

This round adds a constructive test outside the previously sampled templates: two isolated norm-two A1 factors joined into A2. All **162 new recipes** have exact defects at actual finite vertices. It also proves a compact-only obstruction when two old orthogonal short walls become a G2 ridge. Neither result closes the general compact-section problem.

## What the sparse-support lemma proves

Let a genuine interior section have a finite Coxeter vertex. Choose one source inequality representing each of its 20 simple facets; coincident projected inequalities are represented only once. Their projected Gram `P` is positive definite. Normalize the spacelike section normal to `e²=1`, and write `a_i=(e,r_i)`. The original Gram is

```
G = P + a aᵀ.
```

Thus `G` is positive definite too. Its nonpositive source-root pairings define a finite crystallographic forest. This argument uses actual facets at the section vertex, not arbitrary projected inequalities that might be redundant.

Three contacts of the same nonzero sign create three strictly negative projected pairings. If these are Coxeter angles, each normalized pairing is at most −1/2, so the Gram evaluated on the positive diagonal-normalized all-one vector is nonpositive. This contradicts positive definiteness. Hence each sign occurs at most twice. Opposite-sign contacts must have a nonzero old pairing, since projecting an old orthogonal pair would otherwise make it obtuse. Two contacts of each sign would require a K2,2 cycle in the old forest. Consequently the total support is at most three.

For three contacts, the two majority-sign roots are the leaves of an old two-edge star. They acquire a negative projected pairing. At least one of the old star edges must cancel exactly, or the projected diagram again contains a triangle. The [small exact gate checker](../scripts/fresh_h20_compact_section_coverage_20261014_gates.py) verifies the signed-forest reduction, the triangle test, and the literal source data used in the G2 argument below. Its [certificate](fresh_h20_compact_section_coverage_20261014_gates.json) is separate from the bounded construction.

There is an additional coverage distinction. If the section vertex is in the relative interior of an old edge, the usual transverse-edge parameterization applies. A section vertex can instead coincide with an old finite vertex, with one or more old inequalities becoming redundant or coincident after restriction. The same support bound holds for the 20 selected actual simple facets, but an open-edge theorem does not automatically cover these endpoint positions. No deformation preserving the Coxeter condition has been proved to remove this issue.

## Existing coverage and genuine gaps

| Local situation | Exact coverage in the repository | Remaining scope |
|---|---|---|
| Any sparse support on the 42 all-right open edges | Every real transverse section excluded; sections containing either associated all-right vertex are separately excluded | Other old edge and finite-vertex types |
| Actual H3/H4 or other noncrystallographic finite angle | Uniform face arguments restrict all actual finite labels to 2,3,4,6 | Crystallographic diagrams remain possible |
| One-contact short A2 to B2 | Complete necessary closure reduction; three complete families excluded exactly | Most family rejections are numerical, not exact |
| One-contact short A2 to G2 | Complete necessary reduction; three complete families excluded exactly | Most family rejections remain numerical |
| Coupled two-/three-contact templates | 264 sampled sections excluded exactly; three specific continuous F4 families excluded exactly | Other ratios, positions, and finite diagrams |
| Zero contacts; isolated A1 rescaling; other one-contact A_n to B_n profiles | No complete universal exclusion found in the audited reports | Continuous cases remain |
| Two same-sign contacts; opposite-sign cancellation; three-contact old stars | Support shape constrained | No complete all-ratio local/global classification |

The relevant precise scopes are in the [finite-angle gate](fresh_h20_crystallographic_angle_gate_20261014.md), [B2 report](fresh_h20_b2_edge_sections_20261014.md), [G2 report](fresh_h20_dihedral_edge_sections_20261014_g2.md), and [earlier coupled report](fresh_h20_coupled_edge_sections_20261013.md). In particular, checking three rational positions in a continuous family cannot close it. Even a complete exclusion of all cusp-containing sections would still leave these compact possibilities.

## New A1+A1 to A2 construction

The complete old corner walk supplies 6,492 recorded literal edge supports. Of these, **1,262** have at least two isolated norm-two A1 roots and are not wholly all-right. These are literal supports, not a new orbit-equivalence count. The scout chooses nine supports: three different component-size signatures for each endpoint-type pair finite–finite, finite–ideal, and ideal–ideal.

Let `r,s` be the two isolated short roots and let their positive unit-normal contacts be `a,b`. Their projected pair is A2 precisely when

```
a² b² = (2−a²)(2−b²)/4.
```

The construction uses three exact contact-square pairs:

```
(a²,b²) = (1/2,6/7), (2/3,2/3), (1,2/5).
```

They have spatial costs `(a²+b²)/2 = 19/28, 2/3, 7/10`, respectively. These are all strictly below one. If `w` is the old dual-weight contribution, the residual Lorentz-plane normal has square `gamma=1−w²>0`. At any interior point of the old edge, either perpendicular residual direction of that square gives a locally valid transverse section. Its new finite component is A2; all other old components are unchanged.

The [scout](../scripts/fresh_h20_compact_section_coverage_20261014_scout.py) tests positions `1/5,1/2,4/5` between source endpoints normalized by the fixed timelike center, and both residual signs. The nine supports, three contact ratios, three positions and two signs give **162 recipes**. They are genuinely different local templates from the earlier A1+A2, A1+A3 and opposite-sign chain constructions; no assertion of geometric inequivalence is made.

The [exact checker](../scripts/fresh_h20_compact_section_coverage_20261014.py) reconstructs all recipes from rational data and sparse square roots. It trusts no numerical normal or truth value. For each recipe it verifies the initial positive-definite 20-wall Gram, the section normal, a released-edge dual, a future timelike next contact satisfying all 210 source inequalities, and an exact 20-wall active set. The released dual has nonzero pairing with the new wall, proving independence; the bad pair is therefore an actual finite ridge. The [exact certificate](fresh_h20_compact_section_coverage_20261014.json) records:

| Exact defect | Cases |
|---|---:|
| Obtuse actual ridge | 142 |
| Squared cosine strictly between 0 and 1/4 | 19 |
| Squared cosine strictly between 1/4 and 1/2 | 1 |

No label cutoff is used for these forbidden intervals. Every recipe fails, but the unsampled continuous A1+A1-to-A2 families remain open.

## Compact-only G2 obstruction for two old orthogonal short roots

Suppose a compact Coxeter section has an actual G2 ridge obtained from two old orthogonal norm-two roots `r,s`. The ridge face `F` is compact of dimension 18. Since G2 has no proper connected finite spherical extension, every boundary facet of `F` is orthogonal to both projected roots; the intrinsic facet normals and angles are inherited, with no further projection distortion.

Use a unit section normal, with contacts `a,b`; their product is positive. A boundary source root `l` with nonzero contact `x=(e,l)` must satisfy

```
(l,r)=x a,       (l,s)=x b.
```

The old three-root Gram is positive definite, because it is the projected positive-definite Gram plus a positive rank-one matrix. A complete literal source census shows that the only case with both old pairings nonzero is `l²=2`, `(l,r)=(l,s)=−1`. Thus a modified boundary is possible only when `a=b`, up to simultaneously reversing the section normal. If the contacts are unequal, all boundary normals are unchanged rational source roots, leading to the rational compact-face contradiction below.

In the equal case the G2 angle equation gives

```
a² = 4√3−6,
N = (l−x e)² = 2−1/a² = 3/2−√3/3,
1/(2N) = (9+2√3)/23.
```

Here `0<N<1`. Any two modified boundaries are divergent, even if their old pairing is zero: the new pairing has magnitude at least `1/a²>N`. A nonright pair with an unchanged short root has squared cosine `1/(2N)`, strictly between `1/2` and `5/8`, hence below `cos²(pi/5)` and above the label-4 value. It is not a Coxeter angle. A nonright pair with an unchanged long root is divergent, since its squared cosine is `1/N>1`. Therefore every actual boundary of a modified facet is an unchanged source root perpendicular to that facet.

If there are no modified boundaries, use the compact face `F` itself. Otherwise use one compact modified facet, of dimension 17. In either case all its intrinsic facet normals are unchanged rational source roots. Compactness makes those normals span the full Lorentzian space of the face. That space is rational, and the source reflection matrices preserve its intersection with the original integer lattice; the checker verifies the requisite integral reflection coefficients for all 210 roots. The compact reflection lattice would have finite index in the integral orthogonal group of this rational indefinite form. In these dimensions Meyer's isotropy theorem makes that integral group nonuniform, a contradiction. The premise is an actual compact reflection polytope with an invariant lattice, not merely rational coordinates.

This excludes the stated G2 join for **all contact ratios and all section positions**, without assumptions on the other wall contacts. It does not exclude noncompact G2 faces, other source pairs, A2/B2 joins, or arbitrary sections.

## Replay and current status

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_compact_section_coverage_20261014.py --verify
python scripts/fresh_h20_compact_section_coverage_20261014_gates.py --verify
```

The construction scout is only a recipe proposer. Its exact replay imports the previously checked rational linear-algebra and sparse-radical engines; it is not a second independent implementation of those engines. The small gate checker uses separate integer/SymPy arithmetic for its finite tables. No new polytope or complete compact-section exclusion is claimed.

The [separate proof and sample-arithmetic review](fresh_h20_compact_section_coverage_20261014_review.md) passes. It checks the compact G2 face argument and independently reconstructs ten witnesses spanning all nine selected source edges, including all contact ratios, positions, signs, and defect types. Its producer-free arithmetic verifies 2,100 wall signs and the actual finite vertices. This is an independent sample check, not a second replay of all 162 recipes.
