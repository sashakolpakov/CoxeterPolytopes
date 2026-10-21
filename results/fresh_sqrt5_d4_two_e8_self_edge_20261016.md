# An admissible arithmetic self-edge for D5 + 2E8

The [standalone exact checker](../scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016.py) and [certificate](fresh_sqrt5_d4_two_e8_self_edge_20261016.json) verify a genuine finite hyperbolic edge whose common diagram is **D4 + 2E8** and whose two endpoint diagrams are both **D5 + 2E8**. Thus arithmetic endpoint compatibility alone cannot eliminate this recurring dual vertex type.

This is a local witness in an admissible quadratic form. The displayed chamber has **infinite volume**, certified below. No finite-volume polytope, arithmetic reflection lattice, completion by further walls, or improved dimension bound is claimed.

Write `phi=(1+sqrt(5))/2`. The six-node core has norm-two roots, D4 center 0 and leaves 1, 2, 3; roots 4 and 5 each pair minus one with leaf 1 and pair `−phi` with each other. Append two orthogonal finite E8 blocks, represented by simply laced stars with arm lengths 1, 2, 4. All other off-diagonal entries vanish. Every nonright Coxeter bond has label 3 or 5.

Eliminating D4 gives the black-root Schur complement

\[
\begin{pmatrix}1&-(1+\phi)\\-(1+\phi)&1\end{pmatrix},
\]

with exact pivots `1` and `−1−3phi`. The full Gram signature is **(21,1)**. Under `phi -> 1−phi`, every pivot is positive and the signature is **(22,0)**. Entries lie in `Z[phi]`; norm-two root reflections preserve the corresponding integral module. The triangle cyclic product is `−phi`, so the generated field is Q(sqrt(5)). These certify first-type admissibility of the ambient form, without asserting that the reflection subgroup is a lattice.

For the vertex omitting wall 4, the first six ray coefficients are

\[
v_4=\left(\phi-1,\phi-1,\frac{\phi-1}{2},\frac{\phi-1}{2},
\frac{3\phi-4}{5},\frac{2\phi-1}{5}\right),
\]

and the remaining sixteen coefficients vanish. Swap the last two entries to obtain `v5`. Their wall pairings are respectively `−e4` and `−e5`, and

\[
v_4^2=v_5^2=\frac{4-3\phi}{5}<0,\qquad
(v_4,v_5)=\frac{1-2\phi}{5}<0.
\]

Each active 21-node Gram is positive definite with type D5 + 2E8. The rays lie in the same time cone; their positive segment has exactly the common twenty walls active. This is an actual finite edge, with `cosh(length)=1+phi`. An additional strictly interior timelike point is also checked.

To certify infinite volume, set `u=−Gcore^(-1)(1/10,1/10,1/10,1/10,1,1)`, extended by zeros. Let `h` vanish on the core and equal `−GE8^(-1)1` in each E8 block. Exactly,

\[
u^2=\frac{176-169\phi}{50}<0,\quad h^2=1240,\quad (u,h)=0.
\]

Then `z=u+t h`, with `t=sqrt((169phi−176)/62000)>0`, is null. Its six core wall pairings are strictly negative, and all sixteen E8 wall pairings equal `−t`. Also `(z,u)=u²<0`. Therefore the chamber contains an open patch of the ideal boundary and has infinite volume. The radical extension is explicit: no nonzero null vector over Q(sqrt(5)) is asserted, which would contradict conjugate positive definiteness.

Generation and saved replay passed using only standard-library exact arithmetic. The source, JSON and report are frozen.

```sh
python3 scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016.py --verify
```
