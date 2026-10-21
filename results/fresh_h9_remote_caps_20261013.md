# Remote vertex insertion in H9

The construction extends the previous 79-wall chamber to **85 actual Coxeter facets**, retaining every old wall. Four new caps cut the remote escape that survived the earlier 3,857-cap bank and create a finite vertex in a different region. A further two coupled caps close an existing open edge and add a two-edge finite path through two new vertices. **The final chamber still has infinite volume:** a strict future spacelike escape is saved and replayed exactly. This is constructive progress, not a new finite-volume example.

The authoritative starting geometry is the frozen [independent 79-wall input](fresh_h9_split_second_frontier_independent_geometry_20261013_input.json), together with its source hashes and the earlier [79-wall report](fresh_h9_split_second_frontier_20261013.md). All normals have square two. Let \(\phi=(1+\sqrt5)/2\), and write \(U_i=\mathrm{seed}i\) for \(0\le i<8\).

## Four caps beyond the old remote gap

The saved remote point is spacelike and satisfies every one of the 79 old inequalities strictly. Projection onto the complement of any seven or six mutually orthogonal seed walls in either eight-wall block was not spacelike. Four rank-five choices were spacelike; the explicit construction uses

\[
S=\{\mathrm{seed}9,\mathrm{seed}10,\mathrm{seed}11,\mathrm{seed}14,\mathrm{seed}15\}.
\]

These preliminary projection checks motivated the choice; they are not an impossibility theorem for rank-six caps. In the chosen complement, an exact \(\mathbf Q(\phi)\) orthogonal frame \((T,V,E,F,H)\) has squares \((-a,b,c,d,e)\), with all five quantities positive. The archive gives every coefficient and norm. Normalize the frame to squares \((-1/2,1/2,1/2,1/2,1/2)\), writing it as \((A,B,C,D,E')\). For

\[
k=64,\qquad c_k=\frac{k^2+1}{2k},\qquad s_k=\frac{k^2-1}{2k},
\]

the four normals are

\[
r_\sigma=s_kA+c_kB+\sigma_1C+\sigma_2D+\sigma_3E',
\quad
\sigma\in\{+++,+--,-+-,--+\}.
\]

Their Gram matrix is exactly \(2I_4\). Each is orthogonal to the five walls of \(S\), and its pairing with each of the other 74 old walls is strictly less than \(-2\); the largest is below \(-2.24622\). All four cut the old remote escape. The point \(c_kA+s_kB\) has square \(-1/2\), lies in the chosen time cone, and has exactly the nine prescribed zero inequalities. Thus it is an actual simple \(9A_1\) finite vertex of the 83-wall chamber.

This construction uses a coupled orthogonal frame, rather than the previous rank-eight edge/anchor cap bank. In particular, cutting its certified common escape proves that these four oriented mirrors are outside that bank.

The [83-wall artifact](../data/enumerations/fresh_h9_remote_caps_20261013_caps.json) contains the frame, positive radicals, rational boost, target and a remaining exact escape. The [replay certificate](../data/enumerations/fresh_h9_remote_caps_20261013_verify.json) contains all 296 dotted pairing bounds and all 74 nonincident vertex inequalities. Norms, right angles and the vertex norm follow from exact frame identities; strict signs use outward-rounded integer intervals at scale \(10^{80}\).

## A finite-angle pair that closes an old edge

The 83-wall escape projects spacelike onto the complement of

\[
S'=\{U_0,U_1,U_2,U_3,U_5,U_6,U_7\}.
\]

The second exact frame normalizes to \((A,B,C)\) with squares \((-1,1,1)\). Set

\[
r_+=s_kA+c_kB+C,\qquad r_-=s_kA+c_kB-C.
\]

These roots have norm two and pairing zero for every positive \(k\). This time \(k\) is fixed by a finite-angle condition, \((r_+,w_{17})=-1\), where \(w_{17}\) is the retained old wall `wall17`. If

\[
\alpha=\frac{(A,w_{17})+(B,w_{17})}{2},\quad
\beta=\frac{-(A,w_{17})+(B,w_{17})}{2},\quad
\gamma=(C,w_{17})+1,
\]

then \(\alpha k^2+\gamma k+\beta=0\). The explicit minus-square-root branch is the unique root greater than one; \(k\approx20.341296925593355\). This is an exact algebraic definition, not a rounded boost. The quadratic identity gives the exact \(\pi/3\) relation without a numerical angle decision.

The pair preserves all 83 walls and cuts the second target. Its 166 old/new relations comprise 14 right angles, one \(\pi/3\) angle and 151 strictly dotted pairs; the largest dotted pairing is below \(-2.98001\). The following finite path is certified:

\[
(U_0,\ldots,U_7,w_{17})
\longleftrightarrow
(S',w_{17},r_+)
\longleftrightarrow
(S',r_+,r_-).
\]

Its vertex diagrams are respectively \(H_2+7A_1\), \(H_3+6A_1\), and \(9A_1\). The middle diagram uses the chain \(U_3\overset5{-}w_{17}\overset3{-}r_+\). Both new vertices satisfy every other wall inequality strictly. The first edge was previously open: before the two new caps, the only positive-semidefinite ninth-wall extension of \(S'+w_{17}\) was its original endpoint wall \(U_4\).

The [85-wall artifact](../data/enumerations/fresh_h9_remote_caps_20261013_coupled_caps.json) and [exact replay](../data/enumerations/fresh_h9_remote_caps_20261013_coupled_verify.json) save the frame, anchor equation, isolated boost branch, all new pairing bounds, three actual vertex certificates and both finite edges.

## Exact geometric limits

There are 3,570 wall pairs in the final 85-wall chamber. The old pair certificate is retained; all new pairs are proved above. A common strict timelike center proves every inequality is a true facet: projecting it orthogonally onto any wall preserves timelikeness and makes every other inequality strictly negative, because all cross pairings are nonpositive.

No new positive affine circuit is possible. Every subset containing a new wall and having no dotted pair lies in one of three explicitly positive-definite nine-wall blocks. The earlier 79-wall affine census was empty. This is not a compactness proof.

The [local frontier certificate](../data/enumerations/fresh_h9_remote_caps_20261013_frontier.json) tests all 2,052 possible ninth-wall extensions of the 27 outgoing rays at the three new vertices:

| New vertex | Finite outgoing edges | Rays open to the light cone |
|---|---:|---:|
| First remote \(9A_1\) | 0 | 9 |
| Attached \(H_3+6A_1\) | 2 | 7 |
| Coupled \(9A_1\) | 1 | 8 |

The two incidences of their shared finite edge are counted separately. This checks those three vertices and the repaired old edge, not every finite vertex or component of the chamber.

Most decisively, the final rational dual-coordinate point in the artifact has

\[
q(y)=\frac{15696113026401562007876097-8715480289629940313714305\phi}{16000000000000000000000000}>0,
\qquad (y,h)_{G^{-1}}=-\frac{25001}{6250}<0.
\]

All 85 inequalities are strict. It comes from a new open region with numerical proposal support \([0,4,9]\). The exact point, not the numerical vertex census, proves infinite volume. The two successfully cut targets and this remaining target are kept separately, providing a reproducible next frontier without suggesting that further additions are guaranteed to close it.

## Reproduction

From the repository root:

```sh
PYTHONPATH=. python scripts/fresh_h9_remote_caps_20261013_verify.py --verify
PYTHONPATH=. python scripts/fresh_h9_remote_caps_20261013_coupled_verify.py --verify
PYTHONPATH=. python scripts/fresh_h9_remote_caps_20261013_frontier.py --verify
```

All three replays pass. They use Fraction arithmetic and integer intervals, without numerical search imports. Forged frame coefficients and substituted timelike “escape” points are rejected for both stages. The construction scripts are [the four-cap producer](../scripts/fresh_h9_remote_caps_20261013.py) and [the finite-angle pair producer](../scripts/fresh_h9_remote_caps_20261013_coupled.py); their numerical geometry only proposes points subsequently certified exactly. The roots are not claimed exhaustive, arithmetic, or to form a lattice.

An [independent source review and different-precision replay](fresh_h9_remote_caps_independent_20261013.md) also passed: all 949 archived interval enclosures were recomputed at 100 digits and nested inside the 80-digit certificates. That review reuses the interval implementation and explicitly does not claim a second arithmetic implementation.
