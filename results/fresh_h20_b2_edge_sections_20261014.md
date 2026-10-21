# A finite closure test for one-contact A2 → B2 sections

The search has not produced a new Coxeter polytope. It replaces a parameter grid by a **finite necessary closure test for every continuous one-contact A2 → B2 family** on the recorded source edges. All eight possible local closing contacts are certified exactly. Three complete continuous families are excluded by exact actual-ridge certificates; the broader census is numerical and is reported separately below.

The source is the known 210-wall finite-volume H21 chamber. Its [complete corner walk](fresh_h21_source_corner_walk_20261013.json) supplies 6,492 literal spherical rank-20 edge supports. Among them, 2,524 norm-two A2 components give 5,048 anchored families. These are recorded representatives, not a claim that all 5,048 families are inequivalent.

## Why only finitely many positions need testing

Let the old A2 roots be `a, b`, with Gram matrix

\[
\begin{pmatrix}2&-1\\-1&2\end{pmatrix}.
\]

Normalize the section normal by `(e,a)=1, (e,b)=0, e²=1`, and require it orthogonal to the other 18 roots of the old edge support. Write

\[
e=w_a+z,\qquad w_a^2=\frac23,\qquad z^2=\frac13,
\]

where `w_a` is the old dual weight and `z` lies in the Lorentzian two-dimensional complement of the edge support. The projected pair has norms 1 and 2 and pairing −1, hence label 4. For a transverse intersection with the interior of the old edge, the 20 projected old roots form a positive-definite Gram matrix and define an actual finite section vertex.

Suppose this section were Coxeter. Its B2 face `F` would be a finite-volume H18 polytope. **Compactness is not assumed.** The 18 other walls at the initial vertex cannot enclose a finite-volume H18 polytope, which requires at least 19 facets. Therefore an old wall outside the initial edge support supplies a further actual facet of `F`. A finite point in that facet gives a spherical rank-three extension of the induced B2: either B2+A1 or B3, with either end of the old B2 available for the new label-3 attachment.

For such a wall `l`, put `d=l²`, `g_a=(l,a)`, `g_b=(l,b)`, and `x=(e,l)`. The complete possibilities are:

| Extension | Necessary equations |
|---|---|
| B2+A1 | `g_b=0, x=g_a, d−x²>0` |
| B3 attached to the modified end `a` | `g_b=0, 5x²−8g_a x+4g_a²−d=0, x>g_a, d−x²>0` |
| B3 attached to the unchanged end `b` | `x=g_a, g_b<0, d−g_a²=2g_b²` |

The literal source has norm classes 2 and 4, with short–short pairings {−1, 0}, short–long pairings {−2, 0}, and long–long pairings {−8, −6, −4, −2}. Positive definiteness of the old three-root Gram is necessary because it is the projected positive-definite Gram plus a rank-one positive matrix. The exact [contact-pool checker](../scripts/fresh_h20_b2_edge_sections_20261014.py) verifies that the three rows above give exactly eight possibilities:

| `(d,g_a,g_b)` | Allowed `x` |
|---|---|
| (2,0,0) | 0 or √10/5 |
| (2,-1,0) | −1 or (−4+√6)/5 |
| (4,0,0) | 0 or 2√5/5 |
| (4,-2,0) | −6/5 |
| (2,0,-1) | 0 |

Each recipe imposes

\[
z\cdot l=x-\frac{2g_a+g_b}{3}.
\]

Together with `z²=1/3`, this gives at most two positions. The projection of an outside wall into the edge plane is nonzero: otherwise its hyperplane would contain the entire old edge, contrary to its exact active set. A nonzero-norm hyperbola in a Lorentzian plane contains no affine line. Thus the reduction loses no continuous family and has no exceptional line of solutions.

## Exact complete checks on three families

The [exact checker](../scripts/fresh_h20_b2_edge_sections_20261014_exact.py) reconstructs every eligible closure wall and every real quadratic branch for the stated families. It includes orthogonal closures and both B3 attachments. Its arithmetic is a recursive rational quadratic tower, including Q(√6, √D) when D ∈ Q(√6). Signs are exact; radical independence is not assumed. The numerical [recipe packet](fresh_h20_b2_edge_sections_20261014_exact_recipes.json) supplies only proposed defect indices, not root bounds or truth values.

| Family | Old endpoint types | Exact transverse recipes, with multiplicity | Replayed numerical representatives |
|---:|---|---:|---:|
| 0 | ideal–finite | 433 | 86 |
| 34 | ideal–ideal | 450 | 172 |
| 660 | finite–finite | 179 | 30 |
| Total | | **1,062** | **288** |

Every one of these 1,062 recipes fails at an actual finite contact: 844 have an obtuse ridge and 218 have squared cosine strictly between 0 and 1/4. Counts include repeated recipes for the same normal. The 288 numerical representatives are not an independently certified cross-field mirror count.

For each rejection the checker proves the initial point is timelike, future oriented, and strict at all 190 nonincident old walls. It releases an old facet and constructs the first contact by an exact collision formula. The new point is timelike and future oriented, satisfies all 210 inequalities, and has exactly the asserted 20 incident walls. Independence follows from the original positive-definite rank-20 Gram and the nonzero new-wall pairing with the released dual direction. Thus the offending pair is an actual ridge, not a test of possibly redundant projected inequalities. Exact rational enclosures for consecutive Coxeter cosine values establish strict gaps without an angle-label cutoff.

The saved [exact certificate](fresh_h20_b2_edge_sections_20261014_exact.json) also records all excluded real roots, including endpoint and exterior branches. In particular, the one numerical remainder on family 660 is a boundary root and contributes no unresolved transverse position. All three continuous families are therefore excluded exactly.

## Broader numerical census

The [numerical producer](../scripts/fresh_h20_b2_edge_sections_20261014_scout.py) finished all 5,048 anchored families in approximately 18 minutes. The [summary](fresh_h20_b2_edge_sections_20261014_scout.json) binds its complete [compressed recipe archive](fresh_h20_b2_edge_sections_20261014_scout.jsonl.gz).

| Numerical checkpoint | Count |
|---|---:|
| Saved closure proposals | 870,300 |
| Rejected at a numerical actual-ridge test | 870,063 |
| Unresolved after numerical tests | 237 |

The separate [endpoint checker](../scripts/fresh_h20_b2_edge_sections_20261014_remaining.py) rebuilds every one of the 237 unresolved recipes with exact arithmetic. **All 237 are double roots exactly at an old endpoint**, not transverse interior positions. Its [saved certificate](fresh_h20_b2_edge_sections_20261014_remaining.json) and saved-result replay pass.

There is therefore no surviving numerical candidate from this census. This is still **not an exact exclusion of all 5,048 families**: the broad producer uses floating-point discriminant thresholds, near-degenerate quadratic handling, deduplication, and geometry. The remainder certificate checks the 237 saved recipes, not completeness of floating-point enumeration or the other 870,063 rejections. The three exact continuous-family exclusions above are stronger and independent of those numerical completeness assumptions.

## Reproduction and scope

Run from the repository root:

```sh
python scripts/fresh_h20_b2_edge_sections_20261014.py --verify
python scripts/fresh_h20_b2_edge_sections_20261014_exact.py --verify
python scripts/fresh_h20_b2_edge_sections_20261014_remaining.py --verify
```

The saved-result replays pass. An independent [source/proof review](fresh_h20_b2_edge_sections_20261014_review.md) found no issue in the finite reduction or exact selected-family arithmetic; it is explicitly not a second arithmetic implementation. To repeat the numerical scout, run `python scripts/fresh_h20_b2_edge_sections_20261014_scout.py`; this rewrites the numerical archive and requires rebuilding its dependent remainder certificate.

The exact closure theorem concerns sections transverse to a recorded source edge, orthogonal to its other 18 walls, with one norm-two A2 pair changed to B2. It does not exclude arbitrary B2 section profiles, multi-contact profiles, or subdivisions. The exact selected-family work establishes actual global ridge defects in those three complete families; it does not establish a new example or a global impossibility theorem for all remaining families. No search process or `RUNNING` snapshot remains from this census.
