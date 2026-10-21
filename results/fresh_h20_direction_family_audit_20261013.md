# All-height obstruction for the old H19 finite-vertex directions

**No finite-volume polytope completion of the retained 87-wall H20 chamber can use only caps whose old spatial component points along an old H19 finite vertex.** This holds for all real heights, every finite Coxeter label, and dotted or parallel relations. It is stronger than the bounded 1,593-entry height bank or the numerical continuous-interval scout. It does not apply to the new 18-short-wall construction in the [H20 continuation report](fresh_h20_targeted_escape_20261013.md), or to replacement of retained walls.

The [standalone audit](../scripts/fresh_h20_direction_family_audit_20261013.py) uses Fraction arithmetic and the earlier independent quadratic-field checker. It imports no search producer, NumPy, or numerical geometry library. Its [exact certificate](fresh_h20_direction_family_audit_20261013.json) contains the census, all 2,046 necessary height bounds, two exceptional-direction proofs, and a common strict escape.

## Independently checked direction census

There are 24 old norm-two walls and 12 old norm-six walls in H19. For the short-wall Gram matrices, the audit rejects a cycle by its nonpositive all-one quadratic form, then uses exact leaf Schur elimination on forests. This independently decides all \(\binom{24}{19}=42,504\) supports, obtaining exactly **2,046** positive-definite supports, matching the saved short-vertex census.

Two norm-six walls cannot occur in a positive-definite subdiagram: their pairings are at most \(-6\). A single norm-six wall can have at most one selected \(\pi/6\) short neighbor, and that neighbor must be isolated among the selected short walls; the excluded three-wall principal determinants are negative or zero. Enumerating the 22,116 remaining supports yields exactly **6,300** further directions. The accepted isolated \(I_2(6)\) block is positive definite, so the remaining decision reduces to the short forest test.

For every one of these **8,346** directions, the checker verifies its exact projection identity, negative norm, time orientation, all incident zeros, and strict inequalities on every other old wall. Thus these are all old finite-vertex directions, not merely a sampled list.

Write a proposed cap as \((x,t)\), where \(x\) lies on one of these future vertex rays and the new coordinate has square three. Positive rescaling is immaterial. Negative spatial scaling violates an old short-wall inequality, and a purely vertical normal violates one of the opposite-height lifted walls.

Every one of the 6,300 long-incidence directions is immediately impossible: its incident old long wall was replaced by two walls with height coordinates \(2\epsilon\) and \(-\epsilon\). Acuteness forces \(t=0\), where the proposed cap has negative norm.

## Exact target and the only two possible cutters

The numerical common-bank point is reconstructed exactly in \(\mathbf Q(\sqrt{10})\), with last coordinate \(-1\), from 19 old-wall equations and the height-\(\sqrt{10}-4\) mirror in direction 1400. Its norm is
\[
w^2=\frac{14193-4482\sqrt{10}}{98}>0.
\]
The checker independently verifies all 87 old inequalities and its future orientation. The first 48 retained walls alone give rational lower and upper bounds on \(t\) for each short direction. Maximizing the target pairing over each interval leaves only directions **1400 and 1511** as possible cutters. The other 2,044 directions cannot cut the point even before imposing angle conditions or positive norm.

For direction 1400, rescale the old projection to square \(-2\). A cutting cap must have
\[
-1\le t<\sqrt{10}-4,
\qquad
\cos^2\theta_{26}=\frac{2(1+t)^2}{-2+3t^2}.
\]
This ratio increases from zero to \(1/2\), with the upper endpoint excluded. Its only possible Coxeter values are therefore zero and \(1/4\), giving
\[
t=-1\quad\text{or}\quad t=\frac{-8+\sqrt{14}}5.
\]
At the first height, the retained wall 49 has squared cosine \((5-2\sqrt6)/6\), strictly between zero and \(1/4\). At the second height, that squared cosine lies strictly between \(5/6\) and \(17/20\). This interval lies between \(\cos^2(\pi/7)\) and \(\cos^2(\pi/8)\), so neither height is Coxeter-compatible. The \(\pi/7\) bound is checked using the largest root of \(64X^3-80X^2+24X-1\); the \(\pi/8\) bound follows from \(\sqrt2>7/5\).

For direction 1511, rescale the old projection to square \(-6\). Any cutter satisfies
\[
-\frac32\le t<\frac{-26+5\sqrt{10}}7<-\frac{29}{20}.
\]
The wall-26 squared cosine is \(2(t+3/2)^2/(-6+3t^2)\), which lies strictly between zero and \(1/4\) except at \(t=-3/2\). That sole right-angle possibility gives squared cosine \(2/3\) with short wall 3, also impossible: \(4\cos^2(\pi/m)\) is an algebraic integer and cannot equal the noninteger rational number \(8/3\).

These arguments use no bound on the finite-angle alphabet. They also cover all-dotted intervals: the particular blocking pair has squared cosine less than one, so a dotted or parallel relation cannot evade its angle constraint.

## Strict escape and scope

Every allowable cap in this direction family has nonpositive pairing with \(w\), and strictly negative pairing with the old center \(h\). The exact point
\[
w'=w+2^{-20}h
\]
is still future spacelike. All 87 old inequalities are strict there, and every allowable real-height cap also pairs strictly negatively with it. For any **finite** selection, the segment from the strict timelike center to this point crosses the light cone while every inequality remains strict. It therefore leaves an open region on the ideal boundary and has infinite volume. No claim about a uniform neighborhood for an infinite cap intersection is needed.

The proof leaves changed old spatial directions, 18-short-wall residual spaces, and retained-wall replacements open. The new [88th cap](../data/enumerations/fresh_h20_targeted_escape_20261013_third.json) does cut this exact universal strict escape: its pairing is
\[
\left(-\frac{78644145}{2097152}+\frac{45\sqrt{10}}4\right)
+\left(-\frac{4718647}{1048576}+\frac{3\sqrt{10}}2\right)\sqrt{63}>0.
\]
The audit separately checks its norm two and old spatial incidence: exactly 18 short-wall zeros and no long-wall zeros. Thus it lies outside the excluded direction family and removes this family's common escape. The cap's full 87-wall angle verification and its new remaining escape are documented in the main H20 report; this comparison does not claim finite volume for the 88-wall chamber.

Replay from the repository root:

```sh
python scripts/fresh_h20_direction_family_audit_20261013.py --verify
```

The complete replay passes. Both the parent agent and the H20 construction agent reviewed the reduction and angle arguments independently and reran the saved certificate. The numerical 1,593-cap bank and continuous-height interval scout are used only to motivate the exact target and are not premises of sampled completeness.
