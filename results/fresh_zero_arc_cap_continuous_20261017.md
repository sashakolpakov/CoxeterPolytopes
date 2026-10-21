# Continuous zero-arc bridge obstruction over all low R endpoints

Exact generation and fresh saved-output replay both passed. This excludes the original single zero-boundary-arc X bridge for every low R endpoint recipe, with no bound on the A or X conic parameters and no bound on X's finite contact labels. The ten high-label R families remain open; this is not a general five-cap or finite-volume obstruction.

The unchanged [63-wall control](fresh_nonarithmetic_release49_62_replacement_20261017.md) has 199 finite vertices, ten open intrinsic S6 ends and a strict infinite-volume escape. The [deletion ledger](fresh_three_cap_pentagon_changed_seed_20261017.md) removes exactly `[12,47,52,64,65]`, losing eleven old finite-vertex certificates and six compact two-face certificates. No loss is credited as restored here.

The proposed disk has boundary `[7,8,44,63,30]`. A covers old edge7–8, R covers8–44 and44–63, and X has zero original boundary arcs and joins A to R along old8. Thus the necessary finite corners are

\[
 S6+7+8+A,\qquad S6+8+A+X,\qquad S6+8+R+X,
\]

in addition to the two R endpoint corners. The frozen finite classification forces A right to S6,7,8. R is right to H4 and has one of seven A2 profiles. All 96 R recipes with its three prescribed labels in2..5 are reconstructed from their literal norm equations, including both branches. A strict native obtuse product or consecutive-angle gap excludes98 branches;94 require the continuous argument below. This is the full low-label recipe domain, not merely the numerical94-root bank.

The source uses the same frozen physical interval arithmetic and finite-diagram helper as the [preceding exact gate](fresh_new_cavity_construction_gate_20261017.md), with new conic, inequality and monotonicity calculations. Every interval endpoint and subdivision is rational. Known R squares and prescribed frame contacts are used literally before subsequent projections; expanding them repeatedly in physical coordinates would introduce unnecessary interval dependence.

## The unbounded conic parameters

Old11 is right to7 and8 and has square two. Project old53 off7,8,11 and normalize the resulting timelike vector to k0. A has exactly two spacelike conic sheets in the complement of7,8. The sheet with positive A·11 is immediately obtuse. The other is

\[
 A(t)=t u_A+v_A/t,\quad t>0,
 \qquad u_A^2=v_A^2=0,\quad u_A\cdot v_A=1.
\]

The retained inequality A·53≤0 gives a necessary exact lower bound on t. The certificate uses its rational lower enclosure t0, enlarging the allowed domain. Other A gates can only reduce this domain.

For each finite corner8,R,X and A2 profile, project X onto the positive span of8,R. Write this part as p and its positive remaining square as d. In the Lorentzian two-dimensional complement choose k0²=−1, k1²=1 and k0·k1=0. Both conic sheets are

\[
 X_s(z)=p+\frac z{\sqrt2}(k_0+s k_1)
        +\frac d{2\sqrt2 z}(-k_0+s k_1),
 \qquad z>0,\quad s=\pm1.
\]

This covers every real norm solution because d>0 for the required finite corner. In a symbolic high-label box the parallel endpoint may give d=0; intersecting the interval for d with `[0,∞)` is a safe relaxation of the strictly positive actual domain. The factor d is preserved explicitly rather than expanded through cancellations.

Each old-wall constraint becomes a quadratic inequality after multiplication by z>0. One X sheet has a strictly positive native-product quadratic. On the other, a retained wall yields an exact positive-root lower bound z≥z0. The saved rational lower enclosure again enlarges the domain; no upper cutoff on either parameter is imposed.

## Intrinsic nonintersection, including finite native contacts

A zero-arc X in this incidence layout must not meet the retained actual S6 facets30 or63. Their projected squares are both2. If qX is X's projected square, the necessary closed relaxation is

\[
 X\cdot30\le-\sqrt{2q_X},\qquad
 X\cdot63\le-\sqrt{2q_X}.
\]

These are projected products. Since30 and63 are right to S6, their native products equal the projected products, but the threshold need not be−2. Finite native angles whose products satisfy this weaker intrinsic threshold are included. Requiring native dotted contacts instead would omit part of the nonzero-A2-profile domain. Parallel intrinsic equality is also allowed as a relaxation. All other tested old products are merely required nonpositive, without imposing their discrete Coxeter angle alphabet.

Positive definiteness of a pair of projected actual-facet normals would produce an actual intrinsic edge by strict-center projection. Therefore a disk with no X–30 or X–63 edge must satisfy the displayed nonintersection inequalities. The statement is about this specified disk incidence, not arbitrary caps with additional old neighbors.

## Global parameter bound and label coverage

For the surviving sheet write

\[
 A(t)\cdot X(z)=tF(z)+G(z)/t,
 \qquad F(z)=az+c+b/z,\quad G(z)=ez+g+h/z.
\]

The certificate proves, for every z≥z0,

\[
 F(z)<0,\qquad t_0^2F(z)-G(z)<0.
\]

It does so by expanding the two quadratic numerators in powers of z−z0 and checking all three coefficients strictly negative. Hence the product strictly decreases with t for all t≥t0. On t=t0 its z derivative has the form A−B/z²; the leading coefficient and its value at z0 are strictly negative, so it decreases with z as well. Finally the exact upper enclosure of the corner value is below−2. Across all saved leaves the largest such upper enclosure is approximately−2.62087. Thus A–X is strictly dotted throughout the entire parameter quadrant, contradicting its required finite corner.

The two finite-corner classifications give1,246 X cases over the94 surviving R branches. Every nonright label≥6 is isolated in its finite component and is represented by the full closed cosine interval `[sqrt(3),2]`. No finite label cutoff is used for X or A–X. The interval tree has1,288 tested boxes:1,266 monotonicity leaves and one prescribed-angle-gap leaf, with maximum depth four. Each monotonicity leaf excludes one X sheet by a strictly positive native quadratic and the other by the unbounded parameter argument. Strict angle-gap witnesses are only sufficient exclusions; an unresolved interval would stop the run rather than be counted as impossible.

The [source](../scripts/fresh_zero_arc_cap_continuous_20261017.py) and [certificate](fresh_zero_arc_cap_continuous_20261017.json) bind nine direct inputs. Generation78520 and fresh replay92563 both exited zero. A separate narrow mathematical review confirmed the two-variable monotonicity implication; it did not repeat the interval arithmetic or the complete domain reduction. Root independently read the new source and report: conic normalization, relaxed lower-root domains, the intrinsic threshold, complete high-X intervals, shifted quadratic signs and whole-quadrant monotonicity passed that narrow review. This did not independently replay the interval arithmetic or the full inherited finite-diagram classification. No gap was identified in the stated scope. Syntax, nine direct source hashes and five local report links also passed. This is an internally checked exact interval certificate using shared arithmetic, not an independent implementation or a proof of geometric completion. High R labels, further releases/deformations, X with additional intrinsic old neighbors, and more zero-arc connectors remain outside its conclusion.
