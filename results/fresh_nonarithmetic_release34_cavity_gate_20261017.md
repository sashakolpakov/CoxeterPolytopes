# A continuous incompatibility gate for two cavity repairs

Two natural endpoint-repair families cannot coexist, even when their free
finite angles vary continuously. This explains the failure of the bounded
two-endpoint scout; it does not rule out other completion topologies or the
successful [wall replacement](fresh_nonarithmetic_release54_replace12_20261017.md).

Let A and B be norm2 normals right to S6. Impose

* `A.7=0`, `-A.12,-A.47` in `[7/5,2]`;
* `B.8=0`, `-B.12` in `[7/5,2]`, `-B.52` in `[1,2]`.

Then **A and B cannot both be acute to every retained wall and to one
another**. In particular this covers all Coxeter labels `(2,m,n)` at
`(7,12,47)` with m,n>=4 and `(2,p,q)` at `(8,12,52)` with p>=4,q>=3,
including the parallel limits.

The [exact source](../scripts/fresh_nonarithmetic_release34_cavity_gate_20261017.py)
and [certificate](fresh_nonarithmetic_release34_cavity_gate_20261017.json)
reduce each norm equation in the exact residual4 geometry. Write

`A=x*P+y*Q+sqrt(2-(x*P+y*Q)^2)*k`,

with k unit spacelike and orthogonal to the three old framing normals;
B has the analogous form. The other square-root branch is obtuse to old8
for A, or old7 for B. Those exclusions are verified by complete interval
partitions of their two-parameter rectangles, so they require no assumption
that the inherited center must be retained.

For the surviving branches, all four partial derivatives of A.B are
strictly negative throughout the full rational parameter box. A single
outward interval evaluation of the four derivatives proves this. Hence
the minimum is at `(2,2,2,2)`, where

`A.B > 13/100` (approximately0.132112409).

This positive product is incompatible with acuteness. It is independent
of any finite-label cutoff or field-degree assumption.

Old7's projected square is4/3, derived from its single A2 bond using the
inverse diagonal2/3; the other framing roots are right to S6 and have
projected square2. The kernel, inverse-frame and cross-Gram constructions
are recorded. Rounded coefficient intervals are outward enclosures of
the exact inherited intervals, not numerical acceptance tolerances.

Generation and fresh saved replay passed. Low contact profiles, nonzero
A2 white attachments and clusters not assigning both endpoint duties to
these two individual caps remain outside this theorem.
