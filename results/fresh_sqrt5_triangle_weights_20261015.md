# Stable triangular-face weights over Q(sqrt5)

For every compact arithmetic Coxeter polytope with ground field
\(\mathbb Q(\sqrt5)\), the rule below gives **weight at least 2 on every
triangular two-face**, including arbitrary additional white nodes in its star.
At every vertex in dimension \(n\ge5\), the total is **at most \(2n-2\)**.
This package proves no quadrilateral lower bound and therefore no new dimension
bound.

The [checker](../scripts/fresh_sqrt5_triangle_weights_20261015.py) and
[certificate](fresh_sqrt5_triangle_weights_20261015.json) use the frozen,
[independently reviewed](fresh_sqrt5_lanner_catalogue_20261015_review.md)
[complete Lannér catalogue](fresh_sqrt5_lanner_catalogue_20261015.md).
The checker imports no producer arithmetic or graph routines. It reconstructs
all 76 three-black markings, checks every corner by exact leaf Schur elimination
in \(\mathbb Q(\sqrt5)\), and computes the weights from the literal labelled
graphs. The all-rank and arbitrary-extension arguments below are part of the
proof; neither is replaced by a finite enumeration.

## Definition on a complete spherical vertex diagram

A planar angle at a vertex corresponds to an unordered pair of nodes of its
full spherical Coxeter diagram. Assign the following nonnegative weight:

1. Start with 1 if the nodes have graph distance 1 or 2, and 0 otherwise.
   Nodes in different components have weight 0.
2. Add 1 if the marked pair is the 5-labelled edge of an entire connected
   component of type \(I_2(5),H_3\), or \(H_4\).
3. Add 1 if the pair consists of the two endpoints of an entire \(H_4\)
   component, hence has distance 3.

There is no correction for \(I_2(10)\). The rule is intrinsic and invariant
under every diagram automorphism. In particular, the endpoint correction is
assigned to an entire \(H_4\), not to an arbitrary four-node path in an A/D/E
component. The finite spherical classification and the compact-polytope
face-diagram correspondence are the standard inputs described in
[Vinberg, sections 2–3](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

## Why all triangles are covered

Let \(W\) be the white nodes of a triangular two-face, and let its three
black nodes be \(b_0,b_1,b_2\). At each vertex of the triangle the diagram
\(W\cup\{b_i,b_j\}\) is positive definite. The full star is indefinite:
projecting off the positive span of \(W\) leaves the three side normals of
a compact hyperbolic triangle. Choose an inclusion-minimal indefinite
subdiagram. Compactness excludes affine subdiagrams, so this is Lannér.
It must contain **all three black nodes**, since any set omitting one black
node lies in a positive-definite corner.

The complete catalogue therefore leaves 12 order-three diagrams, six
order-four diagrams with four markings apiece, and four order-five diagrams
with ten markings apiece:
\[
12+6\binom43+4\binom53=76.
\]
All 76 are checked, without discarding markings using unproved ambient
realizability assumptions. The base distance-two rule already gives at least
2 in 69 cases. The other seven have base total 1 and are as follows. Node
numberings and the three individual corner weights are recorded in the JSON.

| Minimal diagram | Deficient markings | Corrected totals |
|---|---:|---|
| D5-shaped tree with terminal 5 on its long arm | 1 | 3 |
| Path 5–3–3–3 | 2 | 2, 2 |
| Path 4–3–3–5 | 2 | 2, 2 |
| Path 5–3–3–5 | 2 | 3, 3 |

These corrections cannot disappear after adding any number of white nodes:

- Every finite spherical graph is a forest. A path between two already
  connected nodes remains the unique path in any spherical extension, so its
  length is unchanged. A formerly disconnected pair can only gain a
  nonnegative weight.
- A connected finite diagram containing a 5-bond is exactly \(I_2(5)\),
  \(H_3\), or \(H_4\). Thus a marked 5-bond keeps its bonus if an
  \(I_2(5)\) or \(H_3\) grows.
- A spherical component containing an induced \(H_4\) has rank at least
  four and contains a 5-bond, hence is exactly that \(H_4\). Its endpoint
  bonus survives every spherical extension.

The checker also reconstructs all induced embeddings among these three finite
5-bond types, and checks all their relabellings. Arbitrary orthogonal
components have no effect. This proves the triangle lower bound for the full
ambient star, not merely for the minimal Lannér subdiagram.

In fact, all seven deficient markings are repaired using **only the two H4
bonuses**. The smaller rule consisting of the distance-two base and the two
H4 bonuses therefore also proves the triangle conclusion. Bonuses on
\(I_2(5)\) and \(H_3\) are optional for triangles; the displayed main rule
retains them because they may help a separate quadrilateral argument and do
not worsen the vertex bound below.

## Complete vertex budget, with no rank cutoff

For a tree, the number of pairs at distance 1 or 2 is
\[
(r-1)+\sum_v\binom{\deg(v)}2.
\]
It is \(2r-3\) on a path of rank \(r\ge2\), and \(2r-2\) on every D/E
diagram. A rank-one component contributes zero. Applying the corrections
gives the following complete table for the spherical diagrams allowed by the
field's finite alphabet \(2,3,4,5,6,10\):

| Component | Total weight | Deficit from \(2r\) |
|---|---:|---:|
| A1 | 0 | 2 |
| A/B paths of rank at least 2, F4, \(I_2(6),I_2(10)\) | \(2r-3\) | 3 |
| D/E | \(2r-2\) | 2 |
| \(I_2(5)\) | 2 | 2 |
| H3 | 4 | 2 |
| H4 | 7 | 1 |

The only component with deficit less than 2 is H4, of rank 4. For total rank
\(n\ge5\), a single component consequently has deficit at least 2; two or
more components have total deficit at least 2 because each has deficit at
least 1. This proves
\[
\boxed{\text{total vertex weight}\le2n-2\quad(n\ge5).}
\]
The rank-four exception is real: H4 has total 7, exceeding \(2n-2=6\).
Independent graph counts and partition optimization through rank 100 confirm
the formulas and attain \(2n-2\) for every \(n\ge5\). The classification and
deficit argument, rather than the rank-100 check, cover all larger ranks.

## Scope and replay

This establishes the complete triangle part of a possible face-weight proof
over \(\mathbb Q(\sqrt5)\). It does not show that every quadrilateral has
weight at least 1, classify its minimal-diagram unions, or control arbitrary
external attachments in a quadrilateral star. Those steps are being studied
separately. There is no dimension conclusion from this package alone.

```sh
python3 scripts/fresh_sqrt5_triangle_weights_20261015.py
python3 scripts/fresh_sqrt5_triangle_weights_20261015.py --verify
```

Generation and saved replay pass, as do source-manifest, syntax, and local-link
checks. The catalogue dependency has a separate completed independent review;
an independent review of this weight package is a separate task.
