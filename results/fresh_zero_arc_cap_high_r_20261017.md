# Continuous zero-arc bridge: all ten high-R families

Status: FINAL/FROZEN. Exact generation and fresh saved replay both terminal
PASS (sessions52261 and79817). Source/data/report are frozen.

This successor treats the ten symbolic high-label R families left open by the
[low-R continuous certificate](fresh_zero_arc_cap_continuous_20261017.md).
It preserves the frozen 63-wall, 199-finite-vertex, 10-open-end control. There
is no new chamber or finite-volume example here. The source is
[fresh_zero_arc_cap_high_r_20261017.py](../scripts/fresh_zero_arc_cap_high_r_20261017.py),
with the [saved interval certificate](fresh_zero_arc_cap_high_r_20261017.json).

## Exact geometric scope

The retained boundary is `7–8–44–63–30–7`, after deleting
`12,47,52,64,65`; the parent certificate records the loss of 11 certified
ambient finite vertices and six compact two-faces. This is the original
zero-boundary-arc X subdivision: A covers 7–8, R covers 8–44–63, and X has old
neighbor 8. In particular the necessary corners include `7,8,A`, `8,A,X`,
`8,R,X`, `8,44,R`, and `44,63,R`, all over the fixed positive support S6.
A–R adjacency is not assumed.

The seven finite A2 attachment profiles and the high-label isolation rule are
inherited from the frozen finite-corner classifier. A is forced right to S6,
7 and 8. R is right to S6 in all ten high families. X can have any of the
seven A2 profiles. Its projected square is q>0. For the stated incidence X
must be intrinsically nonintersecting old30, old44 and old63, hence

    X·30 ≤ −sqrt(2q),   X·44 ≤ −sqrt(2q),   X·63 ≤ −sqrt(2q).

These are projected, not native norm-2 thresholds. Native finite contacts
satisfying these intrinsic inequalities remain included. All other old
contacts used by the gate are merely necessary acuteness or Coxeter-angle
conditions; unused native gates only further restrict the domain.
The old44 inequality strengthens the low-R package's relaxation but is
necessary for the same zero-arc incidence: old44 is an actual side of the
retained face, while X has no edge on it. Projection off the positive S6
support preserves acuteness (inverse Stieltjes Gram). With a common strict
center, a positive-definite pair of projected sides has an actual finite
intersection by the usual projection argument. Thus such a pair would create
an unprescribed X–44 edge. The nonpositive two-by-two determinant condition
above is necessary; parallel equality remains a conservative relaxation.

There is one useful strengthening of the frozen A lower bound. Write
`A(t)·53=a t+b/t`, with certified a<0 and b>0, so this product strictly
decreases on t>0. At its unique zero, the exact interval for `A·9` lies
strictly between the label10 and label11 products (approximately
−1.90335325628). Thus the right-to53 event is not native-Coxeter. Every allowed
A must instead satisfy A·53≤−1, since a norm-2 Coxeter product strictly between
−1 and0 is impossible. The positive root of `a t²+t+b=0` therefore gives the
new relaxed lower bound t0≈2.07128996789744. This excludes a spurious relaxed
corner near R44 label15 at which the earlier acute-only bound was insufficient;
it is not an extra assumed contact or a numerical exclusion.

The ten R label families, in frame order `[8,44,63]`, are

    [2,2,H], [2,H,2], [3,2,H], [4,2,H], [5,2,H],
    [H,2,2], [H,2,3], [H,2,4], [H,2,5], [H,2,H],

where each H independently ranges over every integer label at least6. Its
cosine parameter is enclosed by the complete interval `[sqrt(3),2]`; the
endpoint2 is a relaxation. The two highs in the last family are independent.

## Five regular frame families

When R·8 has a fixed label at most5, the positive two-frame `[8,R]` remains
nonsingular. The source reconstructs both R norm branches over high-R cosine
boxes, reinstates exact literal frame products, and applies the frozen
continuous X conic calculation, now also using intrinsic nonintersection44.
It includes all finite X–8–R corner profiles,
including isolated high X labels. A and X parameters are unbounded continuous
variables; neither is sampled. A successful leaf either gives a strict native
obstruction, an assigned cosine gap, or the inherited monotonicity certificate
showing A·X<−2 on the whole parameter quadrant. Unresolved interval boxes are
subdivided; a depth or node limit fails generation rather than accepting a box.

For the family `[2,H,2]`, direct repeated projection through old53 is poorly
conditioned. The source instead uses an algebraically regular orthonormal
complement. Let W be the projection of old30 off `[8,44,63]`, and k=W²>0.
For c=−R·44 and delta=4−c² the two literal R branches are

    R=−(c/2)·44 + sigma sqrt(delta/(2k)) W.

With g=8·63<−2, define

    k0=[63−(g/2)·8]/sqrt(g²/2−2),
    k1=[sqrt(delta)/(2sqrt(2))]·44 + [sigma c/(2sqrt(k))] W.

They are right to8,R, with k0²=−1, k1²=1 and k0·k1=0. In particular k1² is
`delta/4+c²/4=1`. These formulas extend to delta=0 without division by
sqrt(delta). The usual two-sheet X conic is then applied in this frame. This
change reduces interval cancellation; it does not remove a geometric branch.

## Five high R–8 families: regularization at the parallel limit

Write c=−R·8, delta=4−c², g=8·63, h=R·63. The literal geometry gives
g<−2, h≤0, and old44 of square2 right to8,R,30,63. Set

    K = 2(g²+cgh+h²−delta),
    N = [delta·63 − (2g+ch)·8 − (cg+2h)·R] / sqrt(K).

Here K>0 uniformly: g²>4, delta≤1, and the other displayed terms are
nonnegative. Direct Gram expansion gives

    N·8=N·R=N·44=0,    N²=−delta,    N·63=−sqrt(K).

For actual finite R–8 labels, delta>0. Their high label forces X right to both
8 and R in any finite corner. Consequently every such projected X is

    X = a N ± sqrt((q+delta a²)/2)·44.

The positive sheet has X·44>0 and is excluded by an actual retained native
wall. The intrinsic63 condition forces a>0. With w=a/sqrt(2q), the remaining
sheet is

    X/sqrt(2q) = w N − s(w)·44,
    s(w)=sqrt(1/4+delta w²/2).

Thus w≥1/sqrt(K), and if N·30<0 also w≥1/(−N·30). If N·30≥0 there is no
solution. The checker uses a positive rational lower bound w0 for the strongest
certified bound. This chart never divides by delta; all coefficient estimates
extend to delta=0. That endpoint is only a limiting necessary-domain relaxation,
not an asserted finite R–8 corner.

The frozen A chart is A(t)=t au+av/t, with a certified t≥t0>0; its other sheet
is obtuse to old11. Set r(w)=s(w)/w. For w≥w0,

    0≤r(w)≤r0=sqrt(delta/2+1/(4w0²)).

For f(t,w)=A(t)·X/sqrt(2q)=t F(w)+G(w)/t the checker proves the interval bounds

    F/w = au·N − r au·44 <0,
    (t0²F−G)/w = t0²au·N−av·N−r(t0²au·44−av·44) <0.

They imply partial_t f<0 for all t≥t0. Put A0=A(t0). The literal old geometry
also certifies A0·44<0. Since delta/(2r)≤r,

    partial_w f(t0,w)
      = A0·N − [delta/(2r)] A0·44
      ≤ A0·N − r A0·44 <0.

The last strict bound and f(t0,w0)<−1 are checked on each closed parameter box.
Hence A·X<−sqrt(2q) everywhere in the allowed quadrant. The necessary finite
corner `8,A,X` instead has A and X right to8 and requires
`(A·X)²<2q`. This contradiction covers all seven X profiles simultaneously.

## Verification and limits

The private 402-point probe was discovery only and is not certificate evidence.
The production checker covers all ten families and 150 finite-contact recipes
using 120 certificate roots. There are 18,240 tree nodes, 9,180 leaves, 9,060
binary splits, and 6,211 assigned-angle-gap leaves; maximum depth is25.
The largest saved ordinary A·X upper bound is below −2.00376156743. In the
regular high-R8 chart the largest saved normalized upper bound is below
−112.04630555, well below the necessary finite-corner threshold−1.

The checker saves exhaustive binary interval trees, exact rational
enclosures and 12 direct source bindings. A separate lightweight audit walked
every partition, verified all saved rejection and monotonicity sign
inequalities, and checked direct hashes, syntax and local links. This checks
certificate structure and signs, not an independent reconstruction of the
physical interval expressions. The angle-gap helper searches consecutive
labels only to find strict exclusion witnesses; its finite search is not a
candidate-label cutoff. Frozen physical interval arithmetic and finite-corner
classification are shared, not independently reimplemented.

Frontier's separate conceptual check confirms the N identities, positive K,
norm parametrization, limiting scope and the projected44 nonincidence
argument. It does not duplicate the A53 event or interval arithmetic. Root's
separate source/proof read passes the A53 Coxeter-gap strengthening, intrinsic
old44 nonintersection, both regular charts and their limiting relaxations,
the derivative bound and finite A–X corner contradiction. No gap was identified
in this stated scope. That review does not independently replay the shared
interval arithmetic or finite-corner classifier. The fresh `--verify` run
reconstructed the full certificate and matched the saved JSON exactly.

Even combined with the low-R result, this excludes only the stated single-X
incidence with intrinsic nonintersection30/44/63. Extra connectors, additional
old neighbors for X, different endpoint duties, wall replacements and other
patch topologies remain open. It is not a general impossibility theorem for
the control chamber or for compact H9 polytopes.
