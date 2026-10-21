# All-angle endpoint obstruction on the new H9 frontier

**Exact scoped obstruction:** a new cap cannot simultaneously remain acute-compatible with all 104 walls, cut their saved escape, and create a finite endpoint on one of the 16 new open edges, for **any integer Coxeter angle labels**. Two complementary certificates cover all residual heights and all angle labels. They do not exclude other edges, wall replacement, or intermediate caps that do not cut this escape.

The two new corners are E7+A1² and D4+5A1. Releasing each incident wall gives their outgoing spherical rank-eight edges; two are already closed and 16 are open. For each open edge, exact arithmetic in Q(phi,sqrt2) enumerates all positive-definite rank-nine extensions in the stated label alphabet. The inverse edge Gram is entrywise nonnegative, so the norm consumed by a partial nonpositive pairing profile only increases as further entries are specified. This justifies pruning entire incomplete profiles and yields exactly **940** positive-definite extension profiles.

For a fixed profile let (u) be its projection to the eight-wall span, and let the residual squared norm be (N>0). Choose an orthogonal complement frame with (v^2=-1), (w^2=1). Every possible norm-two normal has the form

\[
r(s)=u+\frac{\sqrt N}{2}\left((v+\epsilon w)s+\frac{-v+\epsilon w}{s}\right),
\qquad s>0,\quad\epsilon\in\{-1,1\}.
\]

The frame is defined by exact projections of the common center and the exact escape. Thus the **1,880 branches** cover every real normal with one of these profiles.

Each retained-wall acuteness inequality, target-cut inequality and orientation inequality becomes a quadratic polynomial in (s>0). The checker encloses its algebraic coefficients by outward rational bounds. Replacing each coefficient with its lower bound gives a **necessary**, possibly weaker inequality. Conservative rational root enclosures then intersect these necessary sets. Python integers and Fractions handle all discriminants and roots; there is no floating-point solver in this certificate.

The necessary intersection is empty on 1,876 branches. The remaining four occur on the single edge obtained by releasing `seed3` from the new E7+A1² corner. Their edge is

`second33, seed1, seed4, seed5, seed7, trivalentRemote, walk2720, wall22`.

Every surviving profile has label 3 at `second33`, and optionally at one additional wall:

| Additional label-3 wall | Necessary height interval, approximately | Pairing with seed6, approximately |
|---|---|---|
| none | 2.88292–2.99877 | −0.14654 to −0.03236 |
| seed7 | 4.41024–4.48817 | −0.19856 to −0.15564 |
| seed4 | 4.71619–4.81723 | −0.37867 to −0.32355 |
| seed1 | 7.12119–7.50753 | −0.23725 to −0.10202 |

The exact rational enclosing intervals and 32 covering subinterval checks per branch are archived. Every seed6 pairing is rigorously strictly between −1 and 0. Since both normals have norm two, this permits neither a right angle nor any angle pi/m with integer `m >= 3`, and is also neither parallel nor dotted. All four branches are therefore excluded with **no bound on the seed6 angle label**.

There are 51 edge/isolated-wall cases in which an `I2(m)` factor with `m >= 6` could occur. They are explicitly listed and covered by a second certificate. These are the only higher-label finite endpoint patterns: an edge of label at least 6 together with any further nonright edge on either endpoint has a non-positive-definite three-node principal Gram, because their squared normalized pairings are at least 3/4 and 1/4. Hence a high-label edge in a spherical diagram must be an isolated two-node component.

For each isolated case, scale the candidate so that its sole nonzero edge pairing is minus one. If that isolated old wall is `a`, the candidate is then

\[
r=-a/2+xv+yw,\qquad x,y\in\mathbf R,
\]

with the same two-dimensional orthogonal complement as above. The second checker drops the norm and angle restrictions entirely. Two retained-wall inequalities `f1(r), f2(r) <= 0` give an exact dual identity

\[
L(r)=C-\lambda_1 f_1(r)-\lambda_2 f_2(r),
\qquad C>0,\quad\lambda_1,\lambda_2>0,
\]

where `L` is the escape pairing, which would have to be negative to cut it. The dual coefficients are reconstructed by an exact two-by-two inverse; nonzero determinants, positive weights and positive constants are certified with outward rational intervals. Thus `L(r) > 0` for every real `x,y`, excluding the entire isolated family. All 51 identities pass. The numerical LP used to propose the two supporting wall names is not a proof dependency beyond providing those names; its numerical weights and minima are ignored.

This changes the next constructive target: another height, finite anchor, or high dihedral label creating an endpoint on these 16 edges cannot also cut the current escape while retaining all 104 walls. One must use another frontier component/face, change retained walls, or first add caps that do not cut this particular escape. The result does not exclude arbitrary caps or compact H9 polytopes.

Artifacts:

- [Exact certificate](../data/enumerations/fresh_h9_component_bridge_20261013_continuous_verify.json).
- [Producer-free all-height checker](../scripts/fresh_h9_component_bridge_20261013_continuous_verify.py).
- [All-real isolated-dihedral certificate](../data/enumerations/fresh_h9_component_bridge_20261013_dihedral_verify.json) and [checker](../scripts/fresh_h9_component_bridge_20261013_dihedral_verify.py).
- [Numerical scout](../scripts/fresh_h9_component_bridge_20261013_continuous.py), which proposed the intervals but is not a proof dependency.
- [104-wall construction and exact remaining escape](fresh_h9_component_bridge_20261013.md).
- [Independent source review and replay](fresh_h9_component_bridge_20261013_continuous_review.md), including a separate exhaustive integer/Fraction enumeration of all 1,048,576 edge/profile assignments and a formal check of the two-wall dual identity.

```sh
python scripts/fresh_h9_component_bridge_20261013_continuous_verify.py --verify
python scripts/fresh_h9_component_bridge_20261013_dihedral_verify.py --verify
```
