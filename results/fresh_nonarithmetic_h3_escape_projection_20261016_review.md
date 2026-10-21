# Independent review: right-support walls cannot cut the renewed H9 escape

**PASS.** The [reviewer](../scripts/fresh_nonarithmetic_h3_escape_projection_20261016_review.py)
and [certificate](fresh_nonarithmetic_h3_escape_projection_20261016_review.json)
independently reconstruct the [positive finite-vertex projection gate](fresh_nonarithmetic_h3_escape_projection_20261016.md).
Generation and fresh saved replay pass. The proof applies to every coherent
acute additional wall perpendicular to S={0,1,2,3,5,6}, without requiring
that it retain the originally chosen center.

The reviewer rebuilds the 53-wall physical geometry using the preceding
independent local reviewer. It projects the renewed strict escape E onto
S-perpendicular and reconstructs four finite vertices with supports

    S+{7,8,11}, S+{7,10,33}, S+{8,10,44}, S+{8,12,34}.

Each literal nine-root Gram has positive pivots. All four projected
vertices are future timelike and satisfy the 176 nonincident inequalities
strictly. The supports agree with the independent complete cone census.
Both the escape projection and the vertices are normalized by -<H,v>=1.

Instead of the producer's test-root pairing matrix, this review selects
a certified nonsingular four-column minor of the literal physical vertex
coordinates. A separate reverse Gaussian inversion gives four strictly
positive coefficients. The normalized escape projection is their exact
convex combination. Indeed, all five vectors lie in the four-dimensional
space S-perpendicular by their exact Gram-projection definitions, and the
nonsingular coordinate minor proves that the four selected coordinate
equations determine the whole vector. Interval residuals merely check
consistency; they are not the justification for the identity.

To justify persistence, take any strict future interior X common to the
retained chamber and the hypothetical new wall r. At a finite vertex v,
write X=a*v+sum(b_i*n_i) in the orthogonal decomposition given by its nine
incident roots. Here a>0, and b_i<=0 because the inverse of its positive
definite Z-Gram is nonnegative. Acuteness gives

    0 > <X,r> = a<v,r> + sum(b_i<n_i,r>) >= a<v,r>.

Thus every saved finite vertex pairs strictly negatively with r. If r
is perpendicular to S, its pairing with E is a positive multiple of the
same positive combination of these four strict negative pairings. It
cannot cut or touch E. This uses the addition's own common interior X,
so it does not assume retention of the old center H.

Any finite family of such right-S additions retains the strict spacelike
escape and has infinite volume. A successful completion must include a
wall nonright to at least one root of S. The theorem does not exclude
those completions, nor does it constrain halfspaces with no common future
interior. It also does not prohibit closing the H3 face while using a
separate nonright cutter elsewhere.

```sh
python3 scripts/fresh_nonarithmetic_h3_escape_projection_20261016_review.py --verify
```
