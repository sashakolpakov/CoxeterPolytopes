# TOOT: a uniform K3 gluing gives three more full-Picard `H10` carriers

## Positive theorem

The three certified finite-volume noncompact chambers

- `D16-D5-H10-30`, with 30 facets;
- `SS2-2D-H10-47`, with 47 facets; and
- `F3-C2-H10-196`, with 196 facets

all admit nine-dimensional very-general complex-projective K3 families on
which their complete chamber-wall systems live in the full Picard lattice.
The carriers form one uniform series

\[
 S_a=\langle2a\rangle+\langle-2\rangle^{10},\qquad a=1,2,3.
\]

The three values correspond respectively to 30, 47 and 196 walls.  Every one
of the 273 oriented primitive normals acts by an integral reflection on the
full carrier, and each complete normal set generates its carrier with index
one.

This is an explicit embedding theorem, not an invocation of Nikulin's easy
stable range: every `S_a` has discriminant length 11 while its prospective
K3 complement also has rank 11.

## Uniform full-discriminant gluing

Take

\[
 T_a=\langle2\rangle^2+\langle-2a\rangle+
     \langle-2\rangle^8,
\]

of signature `(2,9)`.  The cyclic discriminant generator of `S_a`, with
quadratic value `1/(2a)`, maps to the cyclic generator of `T_a` with value
`-1/(2a)`.

On the remaining ten order-two coordinates, the required finite quadratic
form identity is

\[
 q_+^{10}\cong q_+^2+q_-^8.
\]

It is completely explicit.  Leave two coordinates fixed and apply the
matrix

\[
 A=J_4-I_4=
 \begin{pmatrix}
 0&1&1&1\\
 1&0&1&1\\
 1&1&0&1\\
 1&1&1&0
 \end{pmatrix}
 \quad\text{over }\mathbf F_2
\]

to each of two four-coordinate blocks.  If `k` is the weight of a binary
vector, `Ax=x` for even `k` and `Ax=1-x` for odd `k`; consequently
`-wt(Ax)/2 = wt(x)/2 mod 2Z`.  The archived verifier checks every
discriminant class, not merely this formula: 2,048 classes for `a=1`, 4,096
for `a=2`, and 6,144 for `a=3`.

Graph gluing along these anti-isometries produces an even unimodular lattice
of signature `(3,19)`, hence the K3 lattice.  The certificate records, for
each value of `a`, the full rational 22-by-22 graph-glue basis, the resulting
integral even unimodular Gram matrix, the 22-by-11 carrier embedding matrix,
and the 22-by-11 complement matrix.  It checks

- both embeddings are integral, primitive and isometric;
- the two images are orthogonal;
- the complement is exactly `T_a`; and
- the combined index is `|A_S|`, namely 2,048, 4,096 or 6,144.

Thus `S_a` embeds primitively in `Lambda_K3`.  Period surjectivity and generic
avoidance of extra algebraic classes give a very-general projective K3 with
`NS(X)=S_a`.  Since `rank(S_a)=11`, each family has dimension `20-11=9`.

## Exact carrier and wall data

| chamber | full Picard lattice | determinant / Smith form | primitive mirror types | strict interior | wall span |
|---|---|---|---|---|---:|
| `D16-D5-H10-30` | `<2>+<-2>^10` | `2^11`; `2^11` | `18*(-2,div2) + 12*(-4,div2)` | square 52, minimum pairing 2 | index 1 |
| `SS2-2D-H10-47` | `<4>+<-2>^10` | `2^12`; `2^10,4` | `5*(-2,div2) + 16*(-4,div2) + 26*(-4,div4)` | square 11,390, minimum pairing 38 | index 1 |
| `F3-C2-H10-196` | `<6>+<-2>^10` | `3*2^11`; `2^10,6` | `80*(-2,div2) + 36*(-4,div2) + 80*(-12,div6)` | square 4, minimum pairing 2 | index 1 |

The 30-wall orientation is replayed from the signs of the published dual
normals; a primitive strict interior is obtained from their signed sum.  The
47 source-oriented roots remain positive after the explicit integral
isometry, and their sum supplies a strict interior.  For the 196-wall row,
the interior is transported through the fixed-lattice basis and explicit
McLeod isometry all the way back to the published Fermat ample class.  Of the
archived normals, 106 already have the positive sign and 90 are reversed.

Every wall record contains its exact coordinates, square, divisibility,
interior pairing, and checks for primitivity, negative square, integral
reflection, involutivity, root negation and Gram preservation.  Smith
reduction gives eleven unit invariants in all three cases.

## Two complete finite-volume lanes

| chamber | literal Vinberg audit | exact GMP cone |
|---|---|---|
| `D16-D5-H10-30` | 10,030 elliptic rank-9 edges, 1,333 finite vertices, 207 cusps, zero bad edges or invalid full affine neighbourhoods | 1,540 extreme rays: 1,333 finite, 207 ideal, 0 hyperideal |
| `SS2-2D-H10-47` | 19,910 edges, 3,136 finite vertices, 216 cusps, zero failures | fresh replay, 3,352 extreme rays: 3,136 finite, 216 ideal, 0 hyperideal |
| `F3-C2-H10-196` | 1,427,720 edges, 212,824 finite vertices, 10,560 cusps, zero failures | 223,384 extreme rays: 212,824 finite, 10,560 ideal, 0 hyperideal |

The Vinberg lane checks exactly two finite-or-ideal endpoints for every
elliptic edge and validates every full affine zero-neighbourhood.  The cone
lane enumerates every extreme ray of the pointed rational chamber cone; zero
hyperideal rays independently proves finite volume.  Ideal rays prove that
all three chambers are noncompact.  The independently obtained finite and
ideal counts agree row by row.

## Exact scope and boundary

These three rows exhaust only the presently audited repository residue of
certified finite-volume `H>=10` chambers having machine-ingested full Gram and
wall-coordinate data but no previous full-Picard carrier.  This is **not** an
exhaustion of ledger-only rows, all conceivable integral realizations, or the
world literature.

The abstract polytopes were already known in the project.  The new result is
their full-Picard characteristic-zero algebraic realization.  These are
**not full-nef** statements: square-minus-four and, in the Fermat row,
square-minus-twelve mirrors cannot all be K3 negative-curve walls.  The
Picard reflections are not asserted to be surface automorphisms or Hodge
isometries.  No worldwide-priority claim is made.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_three_h10_uniform_full_picard_k3_carriers.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_three_h10_uniform_full_picard_k3_carriers.py
```

Frozen payload SHA-256:

`aaac645cdc0d918dc7d5ea50a7944690bb4a1705b089820d9c37a16432c63d32`
