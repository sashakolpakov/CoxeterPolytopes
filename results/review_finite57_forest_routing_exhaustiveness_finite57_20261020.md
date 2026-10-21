# Independent finite57 forest/routing audit and a precise coarse-screen insertion

Audited source: `manuscripts/finite_volume57.tex`, SHA256
`89bd7658a909c566185544484c4a95f8d176fb5fe476cffa7befad81602ca043`.
This is a fresh mathematical source audit, not a replay or a whole-paper
certification. No scientific program, import, solver, Lean or TeX job was run;
the canonical manuscript and frozen evidence were not edited.

Review history: predecessor SHA256
`101718612fde5119adfa4633115b0a9cab5c3458d87ba56429612c88cac8f2c0`
used “actual A2 count” for the *coarse* variable \(j\). Root's read identified
the necessary clarification below: it counts all rank-two representatives,
including images of B2 and arbitrary dihedral types. Coefficients are
unchanged. The later fine refinement still fixes actual A2 counts.

## Findings

The last-component Bellman induction is correct. Every transition decreases a
positive integer rank; checking all predecessor upper inequalities bounds
every residual multiset, and a tight predecessor certifies equality whenever
attainment is claimed. The fine actual-type refinement in
Appendix `app:finite-contexts` is also exhaustive: it fixes every rank9--16
component, all A17 components, every rank-five component, and the actual A1/A2
counts. All remaining ordinary interactions between residual components are
zero. The D14/A3 presence split reserves one A3 and charges its Boolean bonus
once; further A3 components do not repeat that bonus. The A1/A2/A4/A10 family
is removed and bounded separately. Extra formal A1 or A4 items in an upper
relaxation do not lose any actual forest, since its complete actual counts
give one of the explicitly enumerated branches.

I found two written coverage omissions in this edition, rather than a
counterexample to the claimed inequalities:

1. The coarse 3298-to-39 screen is not numerically defined by the printed
   actual `b,j` tables alone. Its representative coefficients include older
   allowances. The insertion below recovers the exact finite upper problem
   used by `C.upper_base`, `C.upper_cross`, and `C.item`; it must not silently
   be replaced by maxima of the final actual tables.
2. The triangle paragraph transfers an elementary payment test to the current
   weights without explicitly treating the lowered affine A10, A11 and D10
   pairs. Their local proofs exclude those pairs at actual triangular
   corners. A separate triangle fragment is being supplied by batch_guard;
   the exact source homes and the check against every other printed change
   are recorded below.

Subject to the explicitly declared finite calculation propositions, the
remaining geometric/domain bridges checked here have no concrete gap:

- Positive-definite projection uses the *full* actual white diagram. Its
  nonnegative inverse preserves every other facet inequality, so it produces
  an actual face. An affine positive kernel produces a nonzero feasible null
  vector. Any additional incident root in a positive-semidefinite corner
  pairs zero with that null vector; positivity of the kernel and acuteness
  then force every contact to that factor to vanish. Thus the factor used in
  the charge really is whole.
- The opposite-edge core contains both opposite black roots: deleting either
  leaves an actual positive-semidefinite corner. Two orthogonal indefinite
  cores cannot occur in a Lorentzian space. The Schur complement and
  inverse-diagonal bounds control arbitrary attached forests, rather than
  only the old displayed core. A singular Gram kernel represents the zero
  vector when its span contains a timelike direction; this is different from
  the nonzero vector of an affine kernel.
- At a cusp the horospherical link is a product of simplices. Varying one
  simplex edge and fixing the other vertices gives the stated factor-fiber
  multiplicity, including a rank-one simplex. There is no ideal cross bonus.
  In the narrow tail, the chosen node can avoid a prescribed rank-one node;
  a repeated omitted node has zero marked charge, and the two endpoint walls
  may coincide without invalidating the Lorentz identity.
- The high-core alphabet is structural. A minimal indefinite core of order
  at least six is a one-node extension of an affine diagram; a bond at least
  five would survive in a connected proper five-node subdiagram and is
  impossible there. Low orders four/five require the stated labels through
  six/five. The unbounded dihedral case is not replaced by a numerical label
  cutoff: the order-three Schur calculation has a zero high-core inverse
  entry at the sole possible distance-seven port, so it leaves its arbitrary
  indefinite three-node Gram block unchanged and gives a second negative
  direction for every value of its parameters. This is the `triangle_gate`
  argument in `fresh_distance8_low_core_20261017.py`.
- For exceptional cheap support, failure of the absolute budget first bounds
  the *entire* nonnegative inverse quadratic cost of the native contacts.
  All eligible inverse-diagonal positions and all their subsets are included;
  omitted contact bits are either zero or at least one. The convex lower
  bound controls unbounded contact magnitudes. The residual cap is used only
  after the forced absolute floors have been removed. A saved PASS by itself
  does not prove any of these geometric assertions.

## Proposed insertion: the exact coarse representative screen

Insert before the sentence beginning “The ordinary component and interaction
inequalities give an upper relaxation” in `app:finite-contexts`. Quantities
with bars in this insertion are **unscaled rational payment counts**; multiply
by 64 when comparing with the manuscript's integer `b,j` and Bellman tables.
This screen bounds the preliminary rule of `prop:oldcount`: it restores the
one B13 entry specified there and precedes the three final rank-15-guarded
changes. Those final changes are evaluated on the eight surviving fibers,
not silently absorbed into this older screen.

Let the fixed representative multiset \(F\), the integers

\[
(c_9,\ldots,c_{17},a_{12},k_5),\qquad
\sum_{r=9}^{17}rc_r+5k_5\le58,
\]

and the residual rank \(R=58-\operatorname{rank}F\) be as just defined.
Use residual representatives

\[
u_r=A_r\ (r=1,2,3,4,7,8),\quad u_6=E_6,
\qquad u_r=D_r\ (17\le r\le58).
\]

The coarse base coefficients are

\[
\begin{array}{c|rrrrrrrrrrr}
T&A_5&B_9&D_{10}&D_{11}&A_{12}&D_{12}&D_{13}&D_{14}&D_{15}&D_{16}&A_{17}\\\hline
\bar b(T)&10&36&45&55&66&65&78&179/2&201/2&219/2&3613/32
\end{array}
\]

and \(\bar b(u_r)=r(r-1)/2\) for \(r\le8\),

\[
\bar b(D_{17})=919/8,
\qquad
\bar b(D_r)=8r-32-[r\ge20]
 -(r-20)_+-(r-21)_+-(r-22)_+\quad(r\ge18).
\]

The following table specifies every nonzero interaction of a residual
representative or a fixed A5 with the fixed representatives. It is symmetric;
an omitted interaction, including every interaction of two large fixed
representatives, is zero.

| first type | A5 | B9 | D10 | D11 | A12 | D12 | D13 | D14 | D15 | D16 | A17 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A1 | 3/8 | 1/2 | 1 | 5/2 | 5/2 | 5/2 | 137/64 | 33/16 | 25/16 | 3/2 | 1 |
| A2 | 0 | 1 | 4 | 5 | 0 | 3 | 151/32 | 4 | 3 | 0 | 2 |
| A3 | 0 | 1 | 4 | 9/2 | 0 | 9/2 | 9/2 | 9/2 | 3 | 0 | 3 |
| A4 | 0 | 1 | 4 | 8 | 17/2 | 6 | 8 | 6 | 4 | 0 | 4 |
| A5 | 0 | 5/8 | 5 | 19/2 | 27/4 | 21/2 | 15/2 | 15/2 | 5 | 0 | 5 |
| E6 | 0 | 5/8 | 0 | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A7 | 0 | 5/8 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A8 | 0 | 5/8 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

In particular, the A1/A5 entry is an unconditional upper allowance, not an
actual ordinary payment. The A4/A12 entry also absorbs its guarded addition.
No residual \(D_r\), \(r\ge17\), has a nonzero interaction in this table.
Define

\[
 K(F)=\sum_{T\in F}\bar b(T)+\sum_{\{T,U\}\subset F}\bar j(T,U),
 \qquad
 w_F(r)=\bar b(u_r)+\sum_{T\in F}\bar j(u_r,T).
\]

Multiset sums count distinct occurrences, including equal component types.
There are no pair interactions between residual representatives, so \(w_F\)
is the complete additive item weight. Let \(V_F(t)\) be its exact-rank
Bellman maximum on

\[
\mathcal R=\{1,2,3,4,6,7,8\}\cup\{17,\ldots,58\},
\]

and \(V_F^{\neg2}(t)\) the same maximum with rank two excluded. At rank zero
both are zero, and unreachable ranks are ignored. Each recurrence has at
most 59 states and checks every allowed final item.

The remaining context terms must not be hidden in \(w_F\). For a reservation

\[
 s\in\{1,2\},\qquad G_0=[s=2\text{ and }c_{10}>0],
 \qquad G_1=[c_9k_5>0],
\]

use it only if \(s\le R\), \(G_0\lor G_1\), and \(s=2\) implies \(G_0\).
For every integer \(0\le j\le\lfloor(R-s)/2\rfloor\), form

\[
 s w_F(1)+j w_F(2)+V_F^{\neg2}(R-s-2j)
   +\frac9{16}j(j-1)G_0+jk_5G_1.
\]

Here \(j\) is the complete rank-two **representative** count after replacing
every actual rank-two type by \(u_2=A_2\), including B2 and every dihedral
type. Thus \(j\) is at least the actual A2 count. The variable \(s\) reserves
the whole A1 witness or witnesses. The coarse value is

\[
 \boxed{\quad
 U(F)=K(F)+\frac{11}{32}c_{17}
  +\max\left\{V_F(R),\ \text{all the reserved expressions above}\right\}.
 \quad}
\]

The \(11/32=22/64\) term is separate from \(\bar b(A_{17})\).
This is the exact upper problem behind the recorded 3298 expanded
alternatives and the selection of 39 classes with \(U(F)>377\).

The value \(U(F)\) bounds the complement of the two separately treated
context families: the family with an additional A1/D5 charge under the joint
A9/A2 guard, and the family containing A1,A2,A4,A10. The complete conditional
rows in the next subsection bound those families, including their overlap.

To prove this coverage, replace each actual fixed component by its representative
and each remaining actual component by its \(u_r\). The base and cross
comparisons are made on *all* actual fibers, not just on attaining forests.
Ordinary same-component guards are absorbed by \(\bar b\), the A4/A12
secondary guard by its upper interaction, and the old A1/D5 guarded
\(3/8\) is absorbed unconditionally by the A1/A5 entry. The additional
\(5/8\) under the joint A9/A2 guard and the A1/A2/A4/A10 increasing family
are bounded separately as stated in the next subsection.

If the remaining A2/A2 context occurs, the forest has at least two A1 and a
rank-ten component; choose \(s=2\) and its complete rank-two representative
count \(j\). If the A2/D5 context occurs without this former guard, choose \(s=1\); its rank-nine
and rank-five witnesses give \(G_1=1\). If both occur, \(s=2\) includes both
terms. If neither occurs, the ordinary Bellman term covers the forest.
The actual A2 count is at most \(j\), and both \(j(j-1)\) and \(jk_5\)
are nondecreasing on nonnegative integer counts. Likewise \(c_9,c_{10},k_5\)
can only enlarge the actual A9, A10 and D5 presence/count conditions.
Consequently the displayed quadratic and linear terms upper-bound all
remaining actual contexts, even when non-A2 rank-two factors are present.
Additional residual A1 components are harmless: the chosen reservation
already pays every context of the actual forest. Some formal branches can
overcount or undercharge another *formal* reconstruction; that is irrelevant
to this upper-bound coverage proof, and no geometric or actual attainment
is inferred from them. Thus every actual forest outside the two extra
context families belongs to a bounded coarse branch, and the separate
family calculations cover all other forests.
The independently stated finite comparison proposition supplies the counts
3298 and 39; the displayed algorithm supplies their complete domain and
meaning.

### Source verification of the recovered coefficients

The coarse module is exactly
`fresh_a16_d16_witness_context_ambient58_frontier_20261018.py`
(SHA256 `69dfbe890dd3de9c7bd54a5a79d50a9cbf38a515741e3ea464d3ef90c4090c9b`).
Its `upper_base`, `upper_cross`, `item`, `fixed_domain`, `fixed_components`
and `finite` define the displayed data and recurrence. The final composition
`fresh_d17_a16_coupled_ambient58_frontier_20261018.py`
(SHA256 `7b3f966491b8780e68b64ad6a92a805708c94a1915c14a6085a78c25e20f4d1c`)
uses this very module as `C`; its `representative_cover` adds 22 per fixed
A17 and identifies the 39 selected keys.

I reconstructed the base and interaction entries by following their exact
arithmetic definitions, without importing or executing the chain. Useful
direct checks in the saved coarse certificate are `rank16_upper=219/2`,
`rank15_upper=201/2`, the A12/D12 bases 66/65, A4/A12=17/2, and
A5/A12=27/4 versus A5/D12=21/2. The residual D17 coefficient is deliberately
the older upper (919/8), not the final actual (227/2). The A17 coefficient
before the final increment is (3613/32), giving (453/4) after it.
These distinctions are why substituting the final actual table would not
reproduce the same coarse finite domain.

## Triangle retention: complete change check and exact proof homes

The early elementary triangle test uses distance at most five, or a length
six/seven path with a protected feature. The current affine reductions that
can affect such a chosen payment are precisely:

| affected factor/pair | allowed support after first deletion | after second deletion |
|---|---|---|
| affine A, order10, marks0,5 | {1,2,8,9} | {3,4,6,7} |
| affine A, order11, marks0,5 | {1,2,9,10} | {3,4,6,7} |
| affine D, order10, marks6,8 | {0,7} | {5,9} |

The two support sets in each row are disjoint, and neither marked endpoint
is allowed in the opposite deletion. In norm-two normalization a nonzero
native contact has magnitude at least one, and the inverse matrices of the
two finite deletions are nonnegative. Their diagonal bounds therefore force
these support restrictions for any third black root at an actual triangle.
That root is consequently orthogonal to the full affine factor. It then
contains its null cusp, whereas the opposite side of an actual triangle
does not contain that cusp: contradiction. Rotations/reflections give all
pairs in the changed orbits. This proves the needed bridge for these rows,
not merely a comparison of stored triangle payments.

Exact proof homes:

- `fresh_order10_cycle_half_transfer_20261018.py`, `profiles` and the triangle
  check in `run`; its report's triangle subsection. Script SHA256
  `42449f6e096d5a59304e53d8b61c679037add56b65b6ae48909fab41fb306ebb`.
- `fresh_cycle11_half_pair_transfer_20261018.py`, `triangle`, and its report.
  Script SHA256
  `ac1cd88faa5fe3df0664b39bce935acb750e1c7ba98463cc74acc97c2c55bfa6`.
- `fresh_order10_d_opposite_fork_fractional_transfer_20261018.py`, `physical`
  and the triangle block in `run`, with `A=6,B=8`. Script SHA256
  `a7b91439148f7cbf834a01a9b7669a2b2fe95993b48537e9cdb2f7fbd83ece01`.

Every other reduced B/C/D affine pair printed at orders10--17 has graph
distance at least eight; higher-order affine entries preserve distance at
most five and every selected feature-protected path of length six/seven.
The featureless A cycle never had an elementary length-six/seven payment.
For finite factors, all sparse changed rank12--17 pairs have graph distance
at least eight; at higher ranks the only shorter reductions require both
marks to have feature-distance at least six. Neither condition can consume
a selected short/protected payment. The four protected exceptional triangle
markings use a whole fixed leaf and neighbor, preserved by the E8 affine
kernel or finite inverse argument already printed in the manuscript.
Contextual changes to long same-component pairs, and every cross-role
transfer, cannot remove one of these selected within-component units.

Thus the missing triangle transfer can be repaired using existing local
proofs and the printed protection argument. I do not assert a new flattened
execution or a fresh check of every old finite certificate.

## Audit limit

The geometric universe and algorithmic coverage described above were
examined independently; exact reported maxima, finite orbit minima and
terminal outputs remain the manuscript's declared finite propositions with
their existing execution records. Hashes identify the sources read; they do
not replace those propositions or their geometric hypotheses. This audit
does not certify the other manuscripts or the currently running replay.
