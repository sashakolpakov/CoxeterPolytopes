# Review of the positive-plane H9 bridge

**PASS.** The two new norm-two caps retain the 86 old facets, cut the saved
escape, and create two actual 9A1 vertices joined by a finite edge. The
resulting 88-facet chamber still has an exact strict spacelike escape.

The [review script](../scripts/fresh_h9_remote_next_20261013_bridge_review.py)
checks the symbolic projection inverse, norm identities, mutual
orthogonality and prescribed π/3 relation separately. It then executes
the existing interval verifier at 160 decimal digits. All 448 new
enclosures lie inside the archived 80-digit bounds, and all non-interval
metadata agrees exactly. See the [saved result](fresh_h9_remote_next_20261013_bridge_review.json).
This is a source review and higher-precision replay of the same interval
engine, not a separate interval-arithmetic implementation.

The six mutually orthogonal retained seed walls are orthogonal to two
further retained walls A,B with A²=B²=2 and p=(A,B)<−2. Their eight-wall
span has signature (7,1); its orthogonal complement is therefore a
positive-definite plane. The displayed inverse of the A,B block,
`[[2,−p],[−p,2]]/(4−p²)`, gives the exact projection to that plane.

Projecting wall 38 gives u with u²=d>0. Projecting and orthogonalizing
the strict center gives v with v²=e>0 and (u,v)=0. Set q=2−1/d>0. The caps
are

```
r0 = −u/d − sqrt(q/e) v,
r1 = −sqrt(q/d) u + v/sqrt(de).
```

Their squares are q+1/d=2, their pairing is zero, and
`(r0,wall38)=−1`. Every other required equality follows from orthogonal
projection. Strict inequalities use outward-rounded integer bounds.
The 172 new-versus-old pairs consist of sixteen right, one π/3 and
155 dotted pairs; the new pair is right.

At each claimed endpoint, the six seed normals, both caps, and one of A,B
form a norm-two orthogonal nine-wall system. Projection of the strict
center gives an actual timelike vertex, strict on all nonincident walls.
The endpoints have eight common walls and are distinct, proving a finite
edge. The new escape is future spacelike and strict on all 88 inequalities.

The [construction input](../data/enumerations/fresh_h9_remote_next_20261013_bridge_caps.json)
and [certificate](../data/enumerations/fresh_h9_remote_next_20261013_bridge_verify.json)
remain frozen. Reproduce the review with:

```sh
python scripts/fresh_h9_remote_next_20261013_bridge_review.py --verify
```
