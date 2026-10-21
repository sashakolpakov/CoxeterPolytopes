# The 97-wall H11 Vinberg prefix is a K3 divisor-wall system

## Result

Let

\[
H_{11}=
\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}\oplus E_8
\]

and let \((X,Z_\tau)\) be the discriminant-11 singular K3 surface with the
conjugate golden marking whose eigenspace meets Shimada chamber \(D_{68}\).
The exact H11 Vinberg transcript through control coordinate \(0<x\leq20\),
followed by the exact next admissible height \(x=6+9\phi\), contains 97
accepted simple norm-two walls.  There is an explicit product \(w\) of 55
integral golden root reflections such that all 97 walls transported
by \(w\) are distinct facets of

\[
\mathbb P\bigl(\operatorname{Nef}(X)\cap E_\phi^+\bigr).
\]

By the even-golden two-root lift, every transported facet is supported by a
pair of disjoint smooth rational curves on \(X\).  Consequently

\[
\boxed{
\text{the selected discriminant-11 K3 eigen-nef chamber has at least
97 paired rational-curve facets.}
}
\]

This is a lower bound, not a complete wall count.  The resumable extension
exhausts 33,572,113 exact recursion nodes at \(x=6+9\phi\) and finds two new
walls.  The resulting 97-wall cone still has the previously certified
hyperideal extreme ray: both new wall pairings are strictly inward.  Thus no
finite-volume or compactness conclusion is claimed.

## The strict Vinberg point is in the completed chamber germ

Put

\[
p=(2,-1,0^8),\qquad D=1+4\phi,
\]

and let \(\rho=E_8^{-1}(1,\ldots,1)^T\).  The height-zero chamber is
\(A_1\oplus E_8\), and

\[
q=-\left(p+\frac1{10}(0,1,\rho)\right)
\]

is in its strict inward sector, with the sign chosen to match the K3 ample
time component.  Exact pairings show that \(q\) is strict on every one of the
97 archived simple roots.

There is also a finite distance certificate excluding every root omitted by
the bound-20 transcript from the segment joining the control vertex to
\(q\).  The exact normalized distance is

\[
\cosh^2 d(p,q)
=\frac{141100+124400\phi}{199879}.
\]

For a positive-height norm-two root with control coordinate \(x\), its mirror
has

\[
\cosh^2d(p,r^\perp)=1+\frac{Dx^2}{4}.
\]

After the extension every unprocessed coordinate is larger than
\(6+9\phi>20\), so its right-hand side is strictly larger than

\[
1+100D=101+400\phi,
\]

which is much larger than the displayed distance to \(q\).  Thus no omitted
mirror can separate the control chamber germ from this strict point.  The
usual Vinberg ordering then identifies its 97 accepted roots as genuine
simple facets; later accepted roots can add facets but cannot remove an
earlier acute simple root.

## Exact K3 Weyl framing

Starting with \(q\), the verifier uses Shimada's complete induced-nef atlas.
At an auxiliary Borcherds wall it changes only the physical induced chamber.
At a `_minustwo` boundary it reconstructs the metric-dual field root and
applies its exact norm-two reflection.  After

\[
55\text{ field-root reflections}
\quad\text{and}\quad
5\text{ auxiliary induced-chamber crossings},
\]

the point is strictly interior to a physical copy of \(D_{181}\).  The
accumulated matrix \(w\) is integral over \(\mathbb Z[\phi]\), preserves
\(H_{11}\), and has canonical digest

```text
7e3a8ef9756a7dffc92f5cb64e6536a6ff53ef71ba8d1f004af29075aed184fc
```

The 97 transported roots are integral, have norm two, are pairwise distinct
as oriented mirrors, and pair strictly with both the transported point and
the fixed ample point.  Their canonical digest is

```text
84e102afa057691bdd067b4de068a65c60e5f554d1a539ccded845da77a350e7
```

As an independent overlap test, all eleven root-wall rays found in the first
ten physical Shimada sections occur among these 97 transported roots.

## The Pell seed calibration

The H11 binary Pell isometry does not supply an infinite chamber stabilizer.
Its ample image crosses four auxiliary induced chambers and one root wall.
Folding across that root leaves exactly the involution already obtained from
the \(D_{68}\) and \(D_{27}\) return galleries.  On field row coordinates it
is

\[
(c_0,c_1,c_2,\ldots)
\longmapsto
(c_0+(1-\phi)c_1,-c_1,c_2,\ldots).
\]

The first twenty Pell powers fold to the corresponding powers of this order-2
element, using exactly \(k\) root reflections for the \(k\)-th power.  Hence
this natural candidate stream gives no same-chamber infinite-order
obstruction.

## Reproduction and scope

The SHA-pinned external atlas is not vendored.  With Shimada's `disc11.txt`
available, run

```sh
PYTHONPATH=src:. python3 scripts/frame_h11_vinberg_prefix_on_k3.py \
  /path/to/disc11.txt \
  --extension data/candidates/h11_vinberg_extension_bound21.json
PYTHONPATH=src:. python3 -m pytest -q tests/test_h11_k3_weyl_frame.py
```

The exact generator is
`scripts/frame_h11_vinberg_prefix_on_k3.py`; the independent bounded Vinberg
base input remains `data/candidates/h11_vinberg_prefix_bound20.json`, and the
resumed exact level is
`data/candidates/h11_vinberg_extension_bound21.json`.  The theorem proved
here is the 97-facet algebro-geometric lower bound and explicit Weyl frame.
Global wall completeness, finite covolume, and compactness remain open.
