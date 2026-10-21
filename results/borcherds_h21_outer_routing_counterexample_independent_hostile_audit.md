# Independent hostile audit of the Borcherds `H21` outer-routing counterexample

**Hostile verdict: PASS.**  This audit independently reconstructs the two
places at which the frozen owner theorem is most vulnerable: the
classification-free absolute reflection-maximality proof and the
`PG(2,4)`/Fano contraction obstruction.  No substantive defect was found.
The audit did not edit any owner byte.

The owner package is
[`borcherds_h21_outer_routing_counterexample.md`](borcherds_h21_outer_routing_counterexample.md),
with terminal hashes

```text
report  fa899977f418d6fa93f7d5270e067c28a648c5e55f2636cb2b362e381534821d
script  c4e21ce23ea67ecd028437caf181adf51460872aff8c891fab891fcf60424523
archive 89dd9af466dcd879fd9764787b7fb6da99bd5222d4a0b17b554af9d6896313a9
tests   4900b836a2afd7814867bb40f21d98491f48acebd7259e68e982be4a9e1d4d87
```

The independent replay is
[`audit_borcherds_h21_outer_routing_counterexample_independent_hostile.py`](../scripts/audit_borcherds_h21_outer_routing_counterexample_independent_hostile.py),
with
[`archive`](../data/enumerations/borcherds_h21_outer_routing_counterexample_independent_hostile_audit.json)
and
[`tests`](../tests/test_borcherds_h21_outer_routing_counterexample_independent_hostile_audit.py).

## 1. Raw reconstruction

Starting from the separately frozen 210-wall coordinate packet, the audit
recovers

\[
 L=U\oplus D_{20}(-1),\qquad
 \operatorname {sig}L=(1,21),\qquad \det L=-4,           \tag{1.1}
\]

and the primitive wall split

\[
 42\times(-2,\operatorname {div}1)
 \quad+\quad168\times(-4,\operatorname {div}2).          \tag{1.2}
\]

The 42 short-root coordinate matrix has Smith diagonal \(1^{22}\), so these
roots generate \(L\), not merely a finite-index sublattice.  Their
inner-product-one graph is connected.  Since those are precisely the odd
\(m=3\) Coxeter edges, all 42 short reflections, and hence their root lines,
are conjugate in the side group.

The finite-volume and complete-root statements remain the clearly labelled
Borcherds literature input pinned by the owner.  This audit checks the new
maximality deduction from that input rather than silently treating host
completeness as absolute maximality.

## 2. The special absolute-maximality proof survives

Let \(W=R(L)\) and suppose that a discrete reflection group \(\Gamma\)
contains \(W\).  Since \(W\) is a lattice, \([\Gamma:W]<\infty\).  The one
implicit point worth making explicit is the projective-rationality step.
Commensurability with the arithmetic group puts each element of \(\Gamma\)
in the rational projective orthogonal group.  If a reflecting generator
\(h\) is represented by a rational similitude \(A=ch\), then in vector
dimension 22

\[
                         \operatorname {tr}(h)=21-1=20.
\]

In other words, the reflecting generator has nonzero trace 20.

Consequently

\[
                         c={\operatorname {tr}(A)\over20}\in\mathbb Q,
\]

so \(h\) itself is rational.  Thus a finite set of coset representatives is
rational and the standard finite-sum construction really does give a full
\(\Gamma\)-stable \(\mathbb Z\)-lattice \(M\) in \(L_{\mathbb Q}\).  This
closes the only terse transition in the owner prose; it does not require
Esselmann's classification.

Scale \(M\) so that \(M\cap\mathbb Qr=\mathbb Zr\) for one short root.
Conjugacy gives the same equality on every short-root line.  Since the short
roots generate \(L\),

\[
                              L\subseteq M.               \tag{2.1}
\]

For \(x\in M\), the short reflection formula gives

\[
             s_r(x)-x=(x,r)r\in M\cap\mathbb Qr=\mathbb Zr.
\]

Hence \((x,r)\in\mathbb Z\) for every short root, and their integral span
gives

\[
                              M\subseteq L^\vee.          \tag{2.2}
\]

The independent discriminant calculation obtains

\[
 A_L\cong(\mathbb Z/2)^2,qquad
 W|_{A_L}=GL(2,\mathbb F_2),                              \tag{2.3}
\]

of order six.  Among the five subgroups of \(A_L\), only \(0\) and \(A_L\)
are invariant.  Therefore

\[
                         M=L\text{ or }L^\vee.            \tag{2.4}
\]

For any rational isometry, preserving \(L\) is equivalent to preserving its
dual.  Thus \(O(L)=O(L^\vee)\) and the full reflection subgroups agree.
Every reflection generator of \(\Gamma\) preserves \(M\), so it belongs to
the full integral reflection group of \(M\), which is \(R(L)\).  This proves
\(\Gamma=W\).

The proof handles every discrete reflection overgroup, not only an integral
overgroup of the initially declared host.  The owner's absolute-maximality
claim is therefore certified in its stated scope.

## 3. Physical square classes and adjunction are exhaustive

The weak-log alphabet is exactly \(-1,-2,-4\).  With one rational common
similarity and arbitrary rational representatives on the wall lines, the
class of

\[
                 {|C_i^2|\over|r_i^2|}
                 \quad\text{in }\mathbb Q^*/\mathbb Q^{*2}
\]

must be common to all facets.  Independent enumeration gives exactly the
classes \(1\) and \(2\).

In class 1, every one of the 42 spanning short lines is a \((-2)\)-curve,
so adjunction forces \(K=0\).  A long line must be \((-1)\) or \((-4)\),
requiring \(K.C=-1\) or \(2\), a contradiction.

In class 2, normalize the physical form to \(2G\).  Every long line is a
\((-2)\)-curve represented by \(r/2\).  The 168 long roots span rank 21;
their primitive orthogonal vector \(w\) has \(G(w,w)=-6\).  Hence
\(K=tw\).  Since a rank-22 smooth rational surface has

\[
                 K^2=10-22=-12=2G(tw,tw)=-12t^2,
\]

one has \(t=\pm1\).  Pairing \(w\) with the short roots splits them into
two sets of 21, with values \(-1\) and \(+1\).  These give exactly the two
dual colourings: 21 disjoint \((-1)\)-curves and 21 disjoint boundary
\((-4)\)-curves.

Thus no physical scale, line representative, canonical class, or adjunction
colouring has been omitted before the geometric obstruction.

## 4. The contraction really forces complex points and lines

Across the two 21-element parts, the intersection-one graph is a symmetric
\(21_5\) design: every point and line has degree five, every two points have
one common line, and every two lines have one common point.  Within each part
all pairings vanish.  The audit independently recovers the displayed
`PG(2,4)` incidence and the same literal Fano subplane.

Contract the 21 mutually disjoint \((-1)\)-curves.  Simultaneous contraction
is smooth, reduces Picard rank from 22 to one, and the resulting smooth
rational surface is \(\mathbf P^2\).  Each boundary curve meets five
different exceptional curves once, so its integral image has

\[
                              \bar B^2=-4+5=1.            \tag{4.1}
\]

An integral plane curve of square one has degree one, hence is a line.  The
intersection-zero entries show that the other contracted points are not on
that line.  Thus the lattice incidence is promoted to an honest simultaneous
configuration of 21 distinct complex points and 21 complex lines; no hidden
MMP assumption is being used.

The seven archived points and seven lines form a Fano subplane.  After
putting three noncollinear points at the coordinate points, the first six
line incidences force

\[
 D=(1,1,0),\qquad E=(1,0,1),\qquad F=(0,1,1).
\]

The final incidence requires them to be collinear, but the determinant is
\(-2\), nonzero over \(\mathbb C\).  The dual colouring gives the dual Fano
subplane and has the same obstruction.  Both simultaneous prime loci are
empty.

## 5. Exact verdict and preserved escape

The frozen theorem therefore has all of its stated conclusions:

- the side group is absolutely reflection-maximal;
- the direct and exact `RLogComp` complex rational branches are empty;
- rank 22 excludes full-trace and projected-root complex-K3 return; and
- the proposed rational/projected-K3 outer routing is not exhaustive.

The conclusion is not a no-AG theorem.  The characteristic-two
supersingular K3 mixed divisor/Cremona-action carrier remains explicitly
preserved.  Nor does the audit challenge the external Borcherds chamber
theorem or the ABSW, Vinberg, Nikulin, or Esselmann results.

## 6. Replay

```bash
python3 -m py_compile \
  scripts/audit_borcherds_h21_outer_routing_counterexample_independent_hostile.py

PYTHONPATH=src:. python3 \
  scripts/audit_borcherds_h21_outer_routing_counterexample_independent_hostile.py \
  --check-archive

PYTHONPATH=src:. pytest -q \
  tests/test_borcherds_h21_outer_routing_counterexample_independent_hostile_audit.py
```
