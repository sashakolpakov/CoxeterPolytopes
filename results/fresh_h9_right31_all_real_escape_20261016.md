# Every escape-cutting normal right to the eight-wall support fails

There is no real spacelike normal perpendicular to
`[0,1,2,5,6,10,30,31]` that both cuts the saved escape and satisfies the
Coxeter conditions with all retained walls. This covers every real coefficient
field, every integer angle label and every positive normalization of the new
normal. It strengthens the earlier bounded first-cap obstruction for cutters
on this support. Compatible caps that do not cut this escape remain possible.

The [exact checker](../scripts/fresh_h9_right31_all_real_escape_20261016.py)
uses the literal form, escape and two-dimensional perpendicular space. Its
[certificate](fresh_h9_right31_all_real_escape_20261016.json) binds the
[independently checked frame](fresh_nonarithmetic_two_caps_20261016_review.md).
The proof needs only retained walls 3 and 7 in addition to the eight right
contacts; it applies in particular to the current 33-wall chamber.

Let `p,k` be the checked frame. Both are perpendicular to the eight support
normals, and

\[
(p,n_3)=-1,\quad(k,n_3)=0,\quad(p,k)=0,\quad k^2>0.
\]

The nine corresponding covectors have rank nine. Thus every vector in the
two-dimensional perpendicular space has the form `a p+b k`. Exact signs give

\[
(k,E)<0,\qquad(k,n_7)<0,\qquad
\tau=-\frac{(p,E)}{(k,E)}<0,
\qquad (p+\tau k)^2>5.
\]

These signs and the strict rational lower bound are checked with outward
rational bounds in the frozen octic field. No numerical approximation is a
premise of the proof.

For a proposed escape cutter `u`, consider its pairing with wall 3. A positive
pairing is already obtuse. A zero pairing forces `u=b k`; cutting the escape
then gives `b<0`, making its pairing with wall 7 positive. In the remaining
case set `a=-(u,n3)>0` and write `u=a(p+t k)`. Cutting the escape forces
`t<tau<0`. Consequently

\[
q(t)=(p+t k)^2=p^2+t^2 k^2
  >p^2+\tau^2 k^2>5.
\]

Since the old normal `n3` has square two, its Cartan two-cycle with `u` is

\[
\frac{4(u,n_3)^2}{u^2n_3^2}=\frac{2}{q(t)}
\quad\text{and hence}\quad 0<\frac{2}{q(t)}<\frac25<1.
\]

This lies strictly between the right-angle value zero and the first nonzero
Coxeter value one. It is neither a permissible finite angle nor a parallel
or dotted pair. Every case is excluded.

The argument imposes no condition on whether the old center remains inside
the enlarged chamber. It also does not say that a non-cutting cap is useless:
such a cap can participate in a longer chain whose later normals have
different right contacts. Other supports and changed retained walls remain
open, and no finite-volume or general dimension claim follows.

Generation and fresh saved replay pass. The arithmetic and rigorous sign
evaluator are shared with the frozen octic review; the new result is the
continuous parameter argument and its exact strict inequalities.

```sh
python3 scripts/fresh_h9_right31_all_real_escape_20261016.py
python3 scripts/fresh_h9_right31_all_real_escape_20261016.py --verify
```
