# A single-contact third leaf with shared single-contact anchors

The [checker](../scripts/fresh_h9_star_single_third_20261014.py) and [certificate](fresh_h9_star_single_third_20261014.json) exclude the saved star state at [26,241 visits](fresh_h9_star_third34_20261014_search_batch_012.json) as part of an all-label family. The owner block has any finite label at least 4; the other three blocks have any finite labels at least 7. The third leaf has one right own contact and one own label **at least 5**. Each anchor has one right own contact and one label 3, with the two nonright contacts using the same root. Every noncore contact remains free.

This is a rank obstruction for that family, not a whole-incidence exclusion or a metric construction. It uses 48 literal positive-kernel guards, three projection budgets, and **11 polynomial identities with 78 strictly positive coefficients**. It introduces no upper label cutoff or arithmetic-field restriction on a possible realization.

## Literal contact domains

Normalize the roots to squared norm 1. Use
\[
S=(3,4;6,7;9,10;12,13),
\]
third leaf 14, and anchors 11,17. In the displayed orientation, third 14 is right to root 3 and has contact cosine \(\beta\ge\cos(\pi/5)\) with root 4. Both anchors are right to 12 and have cosine \(1/2\) with root 13. Owner and last-root reversals are separately checked against the complete incidence lists.

Twenty-four actual affine 4–3–4 paths force the blocks of \(S\) mutually right. An actual eight-root support retains root 4 and all foreign blocks for the third leaf. Its retained own cost exceeds \(1/2\); any foreign nonright contact would add at least \(1/2\). Positive definiteness therefore forces the third leaf right to every foreign block.

The two anchor budgets retain the own label-3 root, at cost \(1/4\). Each anchor has at most one foreign nonright contact, necessarily label 3. Actual three-root supports with a middle block show that even this contact is impossible there: relaxing an internal label at least 7 to 6 gives the affine 6–3 path. Thus foreign contacts can only use the owner block. The same path shows that an owner foreign contact is possible only when the owner internal label is 4 or 5.

The remaining actual supports give short cycle restrictions:

- The anchors are mutually right: their shared own root would otherwise form an affine all-3 triangle.
- Across both anchors, at most one owner foreign contact is possible. Two such contacts create an actual affine 4-cycle if they use the same owner root, or an actual 5-cycle if they use different roots. All cycle bonds can be relaxed to label 3.
- Each third-anchor label is at most 3, by the affine 4–4 path through root 4.
- The third leaf cannot have nonright contacts with both anchors: the shared root gives an actual affine 4-cycle.
- If one anchor contacts root 4, its contact with the third leaf is right by an actual triangle. The other anchor's third contact is also right, by an actual 5-cycle through root 4, the third leaf, the other anchor, the shared root, and the first anchor.

Every support is checked inside an actual vertex. For each positive null vector, adding further acute contacts or increasing the stated labels only decreases its Gram quadratic value, so the relaxed affine witness contradicts strict positive definiteness. The conditional cycle guards are applied only after their required contacts have been established.

These arguments leave seven residual patterns: no foreign contact with third-anchor cosines \((0,0),(1/2,0),(0,1/2)\); one foreign contact at root 3 with those same three choices; or one foreign contact at root 4 with both third-anchor contacts right. The anchor having the possible foreign contact can be named first without assuming a global incidence symmetry: the domain supports are checked for both anchor roles.

## An exact bound for the unlimited tail

Let \(D=1/\sin^2(\pi/m)\), where the last internal label \(m\ge7\). We use the strict rational bound
\[
D>100/19.
\]
To prove it, put \(u=\cos(2\pi/7)>1/2\). The exact identity
\[
T_7(u)-1=(u-1)(8u^3+4u^2-4u-1)^2
\]
shows that \(u\) is a root of the displayed cubic. The cubic is strictly increasing on \([1/2,1]\), and its value at \(31/50\) is negative. Hence \(u>31/50\), so \(\cos^2(\pi/7)>81/100\). Monotonicity in the finite integer label proves the bound for all \(m\ge7\), without sampling large labels.

## Rank in the two-dimensional complement

Write \(c\) for the owner cosine, \(q=1-c^2>0\), and let the first anchor carry the possible owner foreign contact. Set
\[
F=T=0\quad\text{if there is none};
\qquad F=\frac1{4q},\quad
T=\begin{cases}\beta c/(2q)&\text{at root 3},\\
\beta/(2q)&\text{at root 4}.
\end{cases}
\]
The third leaf is right to all foreign blocks, and the anchors are right to both middle blocks. After projecting off \(S\), their Gram matrix is exactly
\[
G_\perp=\begin{pmatrix}
1-D/4-F&-D/4&-\delta_x-T\\
-D/4&1-D/4&-\delta_y\\
-\delta_x-T&-\delta_y&1-\beta^2/q
\end{pmatrix}.
\]
The seven patterns above give the permitted \((F,T,\delta_x,\delta_y)\).

For the three no-foreign patterns, the checker expands \(q\det G_\perp\) after
\[
\beta=(1+\sqrt5)/4+x,\qquad
c=1/\sqrt2+y,\qquad D=100/19+z.
\]
All ten coefficients in each polynomial are strictly positive. Here \(x,y\ge0\) and \(z>0\). For the four foreign patterns, the owner label is exactly 4 or 5; substitution of the corresponding exact cosine gives eight polynomials in \(x,z\), each with six strictly positive coefficients. All signs in \(\mathbb Q(\sqrt2,\sqrt5)\) are checked by rational square comparisons, and every polynomial has a positive constant.

Thus \(\det G_\perp>0\) in all cases. But \(S\) is a positive-definite eight-dimensional support inside a ten-dimensional Lorentz space. Its orthogonal complement has dimension two, so the Gram matrix of three projected vectors must have rank at most two. The nonzero determinant is a contradiction.

## Verification and scope

Primary generation and saved replay pass. The [independent review](fresh_h9_star_single_third_20261014_review.md) also passes generation and saved replay. Its standalone `Fraction` implementation reconstructs all 48 positive-kernel guards with explicit fixed, lower-bound and conditional edge dependencies; all three literal budgets; the Chebyshev tail bound; and all 11 polynomials and 78 coefficients. It uses permutation determinants and exact polynomial division, without producer, SymPy or predecessor-bank arithmetic imports.

```sh
python scripts/fresh_h9_star_single_third_20261014.py --verify
```

The proof covers shared single-contact anchors; it does not assert the corresponding opposite-contact family or exclude the entire star incidence.

## Literal transfers and continuation

The separate [transfer checker](../scripts/fresh_h9_star_single_third_20261014_gates.py) and [certificate](fresh_h9_star_single_third_20261014_gates.json) bind the independent theorem review and the frozen [twelve-core transfer review](fresh_h9_star_core_transfers_20261014_review.md). All four owner/last-root orientations are transported through all twelve maps. The checker rechecks the entire induced face complex on each image, all meeting pairs, and **612 occurrences** of the new path, cycle and budget supports. These are local incidence transfers, not assumed whole-nerve symmetries.

The 48 new predicates bring the total from 301 to **349**. Transfer generation and saved replay pass. The independent reviewer separately read this transfer source and found its logical scope sound.

The earlier 301-predicate files and batch 012 remain unchanged. The [resumable search](../scripts/fresh_h9_star_single_third_20261014_search.py) resumes from the saved 26,241-visit cursor, retaining the existing stack and every pending sibling. [Batch 013](fresh_h9_star_single_third_20261014_search_batch_013.json) stops after **710 new visits**, at **26,951 cumulative visits**, with nine pending frames. The run takes 36.58 seconds; 99,290 of the 100,000-visit allowance are unspent because it reaches the first new necessary-label survivor. The new rules prune 83 entered states. Generation and the final saved replay both pass.

The new outside assignment keeps the owner block at label 4 and third leaf 14 at \((2,5)\). Its other three internal block labels remain in the unlimited category \(m\ge7\), but **each high block now has opposite single-contact leaves**: its positive leaf has \((2,3)\), and its cut leaf has \((3,2)\). This is outside the present shared-anchor theorem. All labels are saved explicitly; the full pending stack is preserved for a later proof or metric attempt.

```sh
python scripts/fresh_h9_star_single_third_20261014_gates.py --verify
python scripts/fresh_h9_star_single_third_20261014_search.py --batch 13 --verify
```

The next concrete target is that opposite-anchor family. Its two anchors have different last-block projections, so the present rank matrix cannot be reused unchanged. No extra parameter bank or resumed batch was started. The star incidence is still open, and the saved survivor is not a Gram-matrix realization.
