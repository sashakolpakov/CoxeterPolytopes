# Review of the exact 839-wall coupled construction

**PASS:** independent mathematical/source review and saved-result replay of the [checker](../scripts/fresh_h22_pi7_coupled_caps_20261013.py) against its [frozen artifact](fresh_h22_pi7_coupled_caps_20261013.json). No mathematical issue was found. This review does not provide a second complete arithmetic implementation.

The center projection is correct. The 21 rational incident walls have a positive Gram matrix; the incident cap is orthogonal to all of them. Projecting the center successively off those walls and the cap therefore produces a timelike vector perpendicular to the entire 22-wall corner. Its normalization Z has square −2. The two rational dual vectors satisfy

\[
u^2=2,\quad v^2=1,\quad (u,v)=1,\quad Z\perp u,v.
\]

Consequently

\[
r_A=(209Z-241u)/120,\qquad
r_B=(241Z-418v)/\sqrt{29281}
\]

have square two and mutual product zero by the displayed integer identities. Their required right angles follow from the dual-vector identities, not from interval overlaps. The remaining **1,632** new old-wall pairs are rigorously dotted; together with **43** right pairs, including the mutual contact, this covers all **1,675** new relations against the frozen 837-wall base.

The two retained edge supports are explicitly bound to certified frontier rows 342 and 332. Their old null limits are cut with positive pairings \(8/15\) and \((-482+418\sqrt2)/\sqrt{29281}\). Thus these are actual closures of two previously certified open directions.

The common 20-wall block is positive definite. Each old/new corner extends a positive 21-wall edge block by an orthogonal norm-two normal. The new/new corner extends the common positive 20-wall block by the two mutually orthogonal new normals. These give positive 22-wall incident Grams. The four displayed points are explicitly timelike and future-oriented, have the required incident zero pairings, and are strict on all other walls of the full 839-wall chamber.

The inverse dual minor is
\(\left[\begin{smallmatrix}1&-1\\-1&2\end{smallmatrix}\right]\).
Hence the two projected old side normals have squares one and two and pairing −1, giving the intrinsic angle π/4. The remaining three quadrilateral angles are right. The four feasible vertices and convexity establish an actual compact H2 face, rather than an isolated spherical-diagram calculation.

The assertion about new affine diagrams is also sound: a diagram containing one new cap can use only its positive right-neighbor block, and a diagram containing both can use only the common positive 20-wall block. Other relations are dotted. Therefore no new affine diagram is introduced.

Both previously bounded cusp links remain unchanged because the new normals are strictly negative on their null rays. The separately identified old null direction also remains untouched. Two strict spacelike witnesses, including the small integer witness of norm **1586**, pass all 839 inequalities. The full chamber therefore still has infinite volume. The positive conclusion is the new coherent endpoint and compact-face construction, not a finite-volume H22 polytope or an exhaustive frontier classification.

Replay:

```sh
PYTHONPATH=src:. python scripts/fresh_h22_pi7_coupled_caps_20261013.py --verify
```

[Review metadata](fresh_h22_pi7_coupled_caps_20261013_review.json) binds the reviewed source and artifact hashes. The shared rigorous interval engine and inherited 837-wall certificates remain explicit dependencies.
