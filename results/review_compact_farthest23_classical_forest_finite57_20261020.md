# Fresh adversarial audit: classical components and the compact23 forest

20 October 2026. Initial manuscript read:
`manuscripts/compact_farthest_vertex_bounds.tex`, SHA256
`fb79910e0d12d05437dd8481957f7209387ff48e1455cb171e4c55a1ac80f080`.
Final binding after the single correction below:
`38c8290bc0d565083d5c29a344e262729251e8caae483d2623b36b4d41f4207b`.
The proof map read was
`docs/compact_farthest_vertex_bounds_proof_map_20261020.md`,
SHA256 `1272ef96d4f611f2ef36928be4678fa7e3d334a0c914ab1b081523fb439c5118`.

**Conclusion:** no unresolved mathematical MUST-FIX or countercase found in
the component/cap/prefix/forest/A20 scope audited here. One strictness
precision was reported and corrected without changing the argument. The
coupled D23/isolated-A1 proposition is a separately audited input to this
report; I did not duplicate that audit.

This was a fresh mathematical source read, not acceptance of earlier PASS
labels. I previously contributed some underlying A-component arguments;
that provenance is disclosed rather than calling this an outside-author
audit of all their original mathematics. I rederived the manuscript's
geometric inequalities and final aggregation and treated earlier review
statuses as unproved. No scientific program, import, optimizer, solver,
Lean or TeX job was run. No canonical source or frozen packet was edited
by me.

## Actual geometry and spectator quantifiers

The geometry in `sec:geometry` uses one actual vector
\(z=E(p-x)\), with \(\|x\|^2=1-E^{-2}<1\). Orthogonal finite
components therefore have additive actual projection energies; independent
edge-cover vectors are never substituted into this sum.

The positive-projection argument has the needed all-wall sign. If \(J\)
is positive, its inverse is nonnegative, so adding the projection
correction preserves every other actual inequality and keeps a timelike
point. An affine positive kernel gives a nonzero feasible null vector,
which is impossible in a compact polytope. Thus positive replacements
containing an already actual wall are actual faces/endpoints; they are
not merely permitted finite diagrams.

For `lem:shadow`, the Schur remainder
\[
C-C_{\cdot i}C_{i\cdot}/C_{ii}
\]
is the inverse on the retained old roots and is entrywise nonnegative.
Acuteness of two actual outside walls then gives the stated reciprocal
alternative and full-tangent pair law. The equality \(T_i=1\) gives a
feasible null point, not an unexamined limiting finite endpoint.
The pair upper bound follows from Cauchy for the sum and concavity of
\(\log(2t/(1+t^2))\) on \(0<t<1\); using half-tangents directly in the
lower product would have been wrong, but the manuscript does not do so.

For `lem:orthogonal`, clearing all other old roots produces
a positive combination of actual walls. Two closed orthogonal directions
would give
\[
\frac{x\cdot u_i}{\tanh L_i}+
\frac{x\cdot u_j}{\tanh L_j}<1,
\]
whereas each summand is at least \((1+\tau_i^2)/2\).
This proves the contradiction with strict interiority. Consequently a
foreign component cost \(q\) is strictly below the wall height \(A^2\).
The minimum positive cost \(1/4\) and the gap to \(1/3\) for a non-A1
foreign component are valid in unit-root normalization.

The finite and negative core points have nonnegative coefficients of
actual normals. Their explicit nonpositive products with the retained
normals and with \(c\) therefore suffice for *all* unknown walls too.
When the displacement is stated to be purely in one component, the text
has either proved \(Q=0\) or positively cleared the entire spectator set.
No arbitrary metric facet is assumed Coxeter.

## Ordinary coverage and the two permitted debits

I read the proofs in `sec:ordinary`, including the small
single/double-contact lists, A5/A6, B5--B8, extended stars and D7/D8.
The common-wall arguments retain the full foreign cost before using
finite classification. In particular, the text does not incorrectly
seal the permitted E7 extension in the B7 or right-fork D calculations.
The B8 F4 bridge is estimated rather than discarded.

The small rate table covers the entire finite classification once
dominant A/B/D components are removed: A1; A2 and every other dihedral
type; A3--A6; B3--B8; D4--D8; E6--E8; F4; H3,H4. Every nonexceptional
rate is strictly above rank/14. The A1 nonshort lower bound
\(5-2\sqrt6>1/10\) is stronger than its ordinary rank charge.

The isolated rank-two defects have spatial support solely in their own
components and height squares below \(1/3\). A short A1 has height square
below \(1/2\), including the two-A1 type J. Disjoint support therefore
gives actual products in \((-1/2,0)\), which are forbidden native
products. This proves the mutual exclusions used in `eq:defects`
and `eq:debits`; it is independent of the type of every
other component. A type J needs two A1 ranks and its light direction
pays more than \(1/4\). Discarding its heavy short direction therefore
remains a conservative bound.

For the D8 three-direction step, the full-tangent law gives
\[
\psi(p):=4p/k-(1+p)^2\ge4(1-\lambda)b^2,\qquad
p=\lambda a^2-(1-\lambda)b^2.
\]
The tangent at \(P(k)\) has derivative
\(4\sqrt{1-k}/k\). Hence \(p\ge P(k)+\alpha b^2\), which together
with the norm bound yields exactly
\((1+\alpha)b^2+\lambda z^2\le\delta\).
The subsequent linear coefficients square to \(8/105,4/35\);
the weighted Cauchy and retained negative quadratic terms give the
printed upper \(101781/200000<.509<P(\sqrt{4/5})\).
This checks the sign and full-versus-half tangent issue in that closure.

## Uniform B/D floor and the component caps

I checked `prop:BDnine` as an actual-projection statement.
The order is noncircular: common 4/5 and 6/7 are first excluded;
their distinct-wall laws then eliminate the low-index common profiles;
only afterwards are the three distinct pair demands combined.
The right-fork D9 exception and short-end B9/F4 exception are explicitly
retained. The isolated native-seven negative point uses a positive
actual combination, not an all-lattice nonreflectivity theorem.

For the final positive dual, direct rational expansion gives
\[
.1(1.091)+.055(1.265)+.305(1.512)=.639835,
\]
and the six Gram upper bounds give
\[
.01(3.415)+.003025(3.633)+.093025(3.852)
 +.011(3.102)+.061(1.899)+.03355(2.474)
=.636435825.
\]
Thus the stated positive margin over \(9/14\) concerns the same
actual component vector even when spectators exist.

For `prop:caps` and `sec:capinventory`, the
foreign-clearing formulas use the *total* spectator set. The isolated
\(u<1,u=1,u>1\) alternatives are all covered; \(u=1\) is a feasible
null point. The zero-own-contact cases are treated before identifying
weak literal packets, so no positive own contact is assumed merely
because an edge has an opposite.

The finite white lists retain E6/E7/E8 and short B/D tail exceptions.
Foreign contacts strengthen the endpoint ratios or are ruled out by
their full finite-white budgets. Therefore the surviving weak packet
has \(Q=0\). The two-cap replacement uses already actual walls,
a positive Schur matrix and a positive inverse; it is a feasible
mixed vertex of the original polytope. Its pressure excludes adjacent
weak indices, not arbitrary abstract packet pairs.

In the A cap the plateau lift is an isometry of the reference potentials:
the constant inserted interval contributes zero to
\(\frac12\sum(t_{i+1}-t_i)^2\). The split peak coefficients remain
nonnegative; the end residuals and their variances stay on the unchanged
ends. The time coefficients increase with \(N\). The repaired middle
index coverage explicitly includes the possible omitted middle index
at \(N=22\), and for \(N\ge23\) obtains both plateau anchors from
strong neighbors. Thus at most six isolated missing floors are charged.
The native scalar rows cover all finite labels, with the tail labels
handled by monotonicity and dotted contacts by their full interval.

The B/D cap argument keeps the special B terminal role and the two D
short layers. Its actual mate, averaged H4 points and averaged negative
A2 points are the reasons these cases are removed; they are not inferred
from a bulk inventory that stops too early.

## Uniform A pressure and the A20 frontier

For `thm:A`, the forcing vectors have positive coefficients
and use four actual pair laws. A common central/end pair is classified
before any energy-dependent discard. The full connected finite
replacement rules exclude weighted internal bonds and include the
rank-eight E8 boundary. Zero-own and spectator single-contact cases
are treated using the foreign-cost inequality rather than an assumed
finite whole star.

I checked the actual H/F/C combinations and the cleared L combination.
For L, positive endpoint residuals give \(\lambda_1,\lambda_2>0\)
and \(8/81<z_0<1/3\). Its pure component formulas are
\[
A^2/B^2=a_N+z_0,\quad U^2/B^2=a_N+3z_0^2,\quad
k/B^2=z_0-3z_0^2.
\]
These justify both its normalized pressure and the far-end correlation.
Normal reversal changes the raw sum's sign; consequently the two-end
centered product is **positive**, \(2BB'/N\).
The manuscript uses that sign, including the separate rank-seven base.
The stronger target retains H,F,C instead of silently repeating the
half-target discards. Every foot remains future and feasible against
all walls.

For `prop:A20`, the four distinct pair demands with flipped
weights are \(183/25\), while the forcing-vector norm is below 68.
The exact Cauchy gap at \(11/14\) is \(673/4375>0\).
The central endpoint discard is recomputed at \(25/3\); its remaining
C/H/F roles are all retained. The two possible far-pair branches are
exhaustive. For a common far pair, positive centered products give the
three correlation bounds and energy costs above \(11/14\).
For two distinct far pairs, \(x\cdot R>183/50\). The L Cauchy envelope
is decreasing in its actual normalized coordinate above .7 and
increasing in its correlation coefficient up to 1; it gives \(<18/5\).
The non-L bound gives \(<181/50\). Both are strictly below \(183/50\).
This is an arbitrary-spectator component inequality, not a paired charge
requiring allocation of an external A1.

## Prefix domain and exhaustive forest aggregation

`thm:prefix` has the right rank domain:
\(B_r:r\ge M+5\), \(D_r:r\ge M+6\), \(6\le M\le16\).
These leave the required retained B4/D5 right side at the last tested
index \(M+1\); no short-tail classification is imported prematurely.
The actual coordinates are decreasing through the relevant prefix.
The reference slopes are positive and decreasing, so every reference
coefficient is nonnegative. Removing one isolated floor averages two
slopes, and nonadjacency makes the open loss intervals disjoint.
The weak-five alternative gives a *different* positive reference with
first coordinate \(>.5965\); it does not assume five strong.
Its variance comparison \(L_M>E_M\) is valid also at \(M=6\).

For final composition, no unproved ordinary rate for A8--A20 is needed.
Dominant means A rank at least7 or B/D rank at least9. If none exists,
the ordinary table and one permitted debit already exceed one in every
dimension considered. Two dominant blocks cannot coexist because each
has energy \(>1/2\). All spectators of the unique dominant block are
therefore ordinary except for the single permitted debit.

Their three lower bounds are exactly
\[
S/14,\qquad (S-1)/14,\qquad 9/100+(S-2)/14.
\]
An empty ordinary sum contributes zero; strictness in a sole defect
case comes from its own \(>9/100\).

- At \(25\le n\le29\), dominant A first has \(S\le7\), so \(r\ge18\)
  and the \(5/7\) bound applies. It gives \(S\le4\), hence \(r\ge21\),
  contradicting the A cap20. Dominant B/D has \(S\le5\) and \(r\ge20\).
  Ranks20/21 have prefix cost \(>.90\) and at least four spectator
  ranks; ranks22/23 cost \(>.94\) and have at least two. Either
  permitted debit still leaves a strict contradiction.
- At \(n=24\), dominant A has \(r\ge17\), then \(S\le4\), so the cap
  forces exactly A20+four spectator ranks. Fully ordinary spectators
  are already impossible. Its remaining short-A1 case costs
  \(>11/14+3/14=1\), and its rank-two defect case costs
  \(>11/14+9/100+2/14=713/700\).
- For B/D ranks19,20,21,22, use respectively the generic theorem with
  \(M=13,14,15,16\) where needed for D. The costs
  \(>.88864,.90715,.9240775,.9401975\), together with the corresponding
  spectator ranks5,4,3,2, exclude every debit case. B is capped at22.
  D23 with one spectator rank remains. An ordinary A1 costs \(>1/10\);
  a rank-two defect or J requires two ranks. Thus its sole possible
  spectator is precisely isolated short type I, the domain of the
  separately audited `prop:D23I`.

The published compact29 ceiling supplies only the ambient upper tail.
The above arguments themselves cover all dimensions24--29. No whole
irreducible-vertex estimate is used in place of an arbitrary-spectator
estimate.

## Accepted one-phrase correction and final binding

At initial source lines1327--1334, the near-4 A-cap argument claimed
correlation \(<1/3\) while its displayed coarse bound supplied
\(6/(N-6)\le1/3\). The uniform direction formula is
\[
\kappa(j,k)=\sqrt{\frac{j(N-k)}{k(N-j)}}.
\]
For \(j\le6\), \(k\ge N-6\) this gives
\(\kappa\le6/(N-6)\le1/3\); coarse equality is possible only at
\(N=24,j=6,k=18\). Strictness of this correlation estimate is unnecessary.
The two actual demands have sum \(>43/50+39/50=41/25\), whereas
Cauchy gives squared sum at most \(2+2\kappa\le8/3\), and
\[
(41/25)^2-8/3=43/1875>0.
\]
Root replaced only “correlation \(<1/3\)” by “correlation at most \(1/3\)”.
I inspected the exact textual diff and stream-hashed the archived
predecessor
`manuscripts/source_archives/compact_farthest_vertex_bounds_before_correlation_precision_20261020.tex.gz`:
its decoded hash is the initial fb79910e hash above. The new live hash is
38c8290b above; no other mathematical or textual change occurred in that
delta.

The proof map's sentences calling the draft uncompiled and the first PDF
build pending are historical execution metadata, not mathematical inputs.
I did not use them, a source-integrity PASS, or any previous review as a
substitute for the implications checked here. Compilation and any final
packet/status update remain root's separate work.
