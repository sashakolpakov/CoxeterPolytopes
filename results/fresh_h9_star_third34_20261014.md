# A dense third leaf against a single-contact anchor pair

The new core obstruction excludes the [saved 20,837-visit star assignment](fresh_h9_star_core_transfers_20261014_search_batch_011.json) as part of an unbounded label family. The third leaf's two own labels may be **any finite integers at least 3 and 4**, in either order. Each of the two anchors has a single own \(\pi/3\) contact; their contacts may use the same or opposite roots. All four block labels are at least 4, and the anchor block label is at least 7. All noncore contacts are free.

The [checker](../scripts/fresh_h9_star_third34_20261014.py) and [certificate](fresh_h9_star_third34_20261014.json) use six small continuous polynomial identities with 66 positive coefficients. They do not enumerate a finite label cutoff or reuse a determinant bank. This is a compact-realization obstruction for the stated family, not an exclusion of the entire incidence.

## Exact domains from actual faces

Normalize all roots to squared norm 1. Let
\[
S=(3,4;6,7;9,10;12,13)
\]
be the positive eight-root support, with third leaf 14 and anchors 11,17. Twenty-four actual four-root supports force the blocks mutually orthogonal through the affine 4–3–4 path. Write the owner-block cosine as \(c\ge1/\sqrt2\), and the third leaf's own cosines as \(\alpha\ge1/2\), \(\beta\ge1/\sqrt2\). Its positive pairing magnitudes with the owner roots are \((\alpha,\beta)\).

An actual eight-root support retains the third leaf's \(\beta\) root and all foreign block roots. Its retained own projection cost is at least \(1/2\). Every foreign nonright contact would cost at least \(2(1/2)^2=1/2\), because each foreign block inverse has diagonal at least 2 and nonnegative entries. The total would be at least 1, contradicting positive definiteness. **The third leaf is therefore right to all foreign blocks.**

For each anchor, a corresponding actual support retains its single own \(\pi/3\) contact, at cost \(1/4\). It permits at most one foreign nonright contact, necessarily \(\pi/3\). The source checks both retained-root choices needed by the shared and opposite anchor profiles.

The actual triples consisting of the third leaf, its \(\beta\) root and an anchor exclude a third-anchor label at least 4 by the affine 4–4 path. Thus each such cosine is 0 or \(1/2\). If an anchor has a nonright contact with an owner root, the actual triangle containing that root and the third leaf forces the third-anchor contact right: a trial label 3 already contains the affine all-3 triangle after relaxing the third own angle. These facts are checked by 30 explicit positive null-vector witnesses and four literal projection-budget supports. No noncore own label is used.

## Projecting the third leaf and the anchor sum

Let \(J\) be the owner-block inverse:
\[
J=\frac1{1-c^2}\begin{pmatrix}1&c\\c&1\end{pmatrix}.
\]
For the two anchors together, let \(n_0,n_1\) count their owner contacts at the two owner roots. Each anchor contributes either none or one \(\pi/3\) contact, so \(n_0,n_1\) are nonnegative integers with \(n_0+n_1\le2\). Nine ordered assignments give exactly six count pairs.

Put
\[
X=(\alpha,\beta)J(\alpha,\beta)^T,\quad
F=(n_0/2,n_1/2)J(n_0/2,n_1/2)^T,\quad
T=(\alpha,\beta)J(n_0/2,n_1/2)^T.
\]
These quantities are nonnegative, and \(X>1\). If \(\delta\) is the sum of the two ambient third-anchor cosines, the angle bounds give
\[
0\le\delta\le d:=\frac{2-n_0-n_1}{2}.
\]

Write \(\mu\ge0\) for the mutual anchor cosine and \(F_M\ge0\) for the squared projection of their sum onto the two middle blocks. Its last-block projection cost \(Y\) is
\[
Y=\frac1{1-t^2}\quad\hbox{(shared)},\qquad
Y=\frac1{2(1-t)}\quad\hbox{(opposite)},
\]
where \(t=\cos(\pi/m)>\sqrt3/2\) for every finite \(m\ge7\). Thus the shared case has \(Y>4\), and the opposite case has \(Y>2+\sqrt3>7/2\). The common strict bound \(Y>7/2\) needs no estimate at a particular high label and no field restriction.

After projection off \(S\), the Gram of the third leaf and the anchor sum is
\[
\begin{pmatrix}
1-X&-\delta-T\\
-\delta-T&2-2\mu-Y-F-F_M
\end{pmatrix}.
\]
Both diagonal entries are negative. Its determinant is strictly greater than
\[
L=(X-1)(3/2+F)-(d+T)^2.
\]
The source verifies the exact difference identity
\[
\det-L=(X-1)(Y-7/2+F_M+2\mu)
 +(d-\delta)(d+\delta+2T)>0.
\]

## Six positive polynomials

Multiply \(L\) by the positive denominator \(1-c^2\), and substitute
\[
\alpha=1/2+x,\qquad\beta=1/\sqrt2+y,\qquad c=1/\sqrt2+z,
\quad x,y,z\ge0.
\]
For every one of the six count pairs, the resulting polynomial has 11 strictly positive coefficients in \(\mathbb Q(\sqrt2)\). Its constant is listed below; every constant is positive.

| \((n_0,n_1)\) | Constant |
|---|---:|
| (0,0) | \(5/8\) |
| (1,0) | \(3/8\) |
| (0,1) | \(13/16-3\sqrt2/8\) |
| (2,0) | \(5/8\) |
| (1,1) | \(13/16-3\sqrt2/8\) |
| (0,2) | \(3/8\) |

Consequently \(L>0\), so the residual Gram has negative diagonal entries and positive determinant. It supplies a negative-definite plane in a Lorentz space of index one, a contradiction. This handles every parameter in the specified continuous ranges and hence all finite integer labels in the theorem.

The checker verifies both owner-root orientations and both last-root orientations on the full vertex and minimal-nonface lists. Together with the two anchor contact patterns, these give eight source predicates. It also checks that the saved assignment satisfies one of them.

## Verification and continuation scope

The frozen predecessor [archive](full_norm_vinberg_and_star_contacts_2026_10_14_audit.json) and batch 011 were freshly replayed before this work; both passed. Their files and saved stack remain unchanged. The primary theorem uses symbolic matrix algebra plus exact rational sign comparisons in \(\mathbb Q(\sqrt2)\), with no imported project arithmetic helper.

```sh
python scripts/fresh_h9_star_third34_20261014.py --verify
```

Primary generation and saved replay passed. The [independent review](fresh_h9_star_third34_20261014_review.md) also passed generation and saved replay: it reconstructs all 30 literal kernels, four budgets, six polynomials, 66 coefficients and eight orientations using separate `Fraction` arithmetic, with no producer or SymPy imports.

## Composition and saved continuation

The earlier independently reviewed [shared](fresh_h9_star_mixed_singles_20261014.md) and [opposite](fresh_h9_star_opposite_core_20261014.md) core theorems cover the third profile exactly \((3,3)\). Their union with the present theorem covers **every pair of third own labels at least 3**: either both are 3, or at least one is at least 4. The two anchors still have exactly one own label 3 and one right contact each.

The separate [composition checker](../scripts/fresh_h9_star_third34_20261014_gates.py) and [certificate](fresh_h9_star_third34_20261014_gates.json) bind all three independent theorem reviews and the [independent twelve-core transfer review](fresh_h9_star_core_transfers_20261014_review.md). They check the integer partition explicitly and transport the four source anchor patterns through those twelve literal induced-core maps. All 408 occurrences of the new guard and budget supports are checked against actual vertices. This adds 48 broad predicates to the frozen 253, giving **301 predicates**. Composition generation and saved replay pass. The independent third34 reviewer also read this partition/transfer argument; it does not assert a full-nerve symmetry.

The [resumable search](../scripts/fresh_h9_star_third34_20261014_search.py) uses the previously reviewed explicit-stack engine with new certificate bindings and gives the new predicates priority. [Batch 012](fresh_h9_star_third34_20261014_search_batch_012.json) resumes frozen batch 011 at 20,837 visits and stops after **5,404 new visits**, at **26,241 cumulative visits**, with nine pending frames and their later siblings preserved. It takes 42.46 seconds; 94,596 of its 100,000-visit allowance are unspent because a new necessary-label survivor was reached. The 48 new rules prune 1,854 entered states. Generation and saved replay pass.

The new outside assignment has block labels \((4,\ge7,\ge7,\ge7)\). Leaf 2 retains own contacts \((3,3)\), but third leaf 14 now has \((2,5)\); each remaining block's two leaves has shared single contacts \((2,3)\). The new dense-third theorem does not cover this changed right contact. This assignment is a necessary-label state, **not a Gram-matrix realization**. Category 7 retains every finite integer label at least 7, and dotted distances are still free. No whole-incidence conclusion is claimed.

```sh
python scripts/fresh_h9_star_third34_20261014_gates.py --verify
python scripts/fresh_h9_star_third34_20261014_search.py --batch 12 --verify
```

The next diagnosis uses the saved \((2,5)\) third profile in a distinct namespace. No earlier file or pending sibling has been discarded.
