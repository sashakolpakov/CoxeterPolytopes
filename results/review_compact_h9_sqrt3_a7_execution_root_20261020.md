# Returned A7 endpoint: root review

20 October 2026. Scope: the executed opposite-A7 edge census, not global
reflectivity or a compact-polytope certificate.

The bounded generation and separate replay both PASS:

* `compact_h9_sqrt3_a7_opposite_generation_20261020.json`: 0.248s,
  55,033,856 bytes peak process-tree RSS.
* `compact_h9_sqrt3_a7_opposite_replay_20261020.json`: 0.243s,
  55,607,296 bytes peak process-tree RSS.

The source/interval proofs already have root and independent written
reviews. The replay uses the same implementation; it is not an
independent implementation.

## Direct interpretation

Put s=sqrt(3) and n_i=(0;-e_i+e_(i+1)). The chosen root is exactly

    d=(5+3s; (4+2s)^3, (3+2s)^5, 2+s).

Its square is 2 and its eight old-wall products are zero except
d.n_3=-1 and d.n_8=-(1+s). Its new corner with n_1,...,n_7 is E8.
The saved endpoint is

    (99+72s; (72+45s)^8,36+27s)=9v,
    v=(11+8s;(8+5s)^8,4+3s),   v^2=1-s<0.

Thus it is precisely the integral E8 vertex of the splitting note,
with future time coordinate and negative square. At the original A8
vertex V0=(9;(4+s)^9), d.V0=(1-s)(5+3s)<0. At v, n_8.v=-4-2s<0.
This checks the orientation of the released edge at both ends.

All three accepted candidates collide at the same endpoint. The other
two are not new simple walls there: their exact decompositions are

    d_A8=3d+n_1+3n_2+5n_3+4n_4+3n_5+2n_6+n_7,
    d_D8=2d+n_1+2n_2+3n_3+2n_4+n_5.

The coefficients are nonnegative. Accordingly the single E8 candidate
is the simple wall, while the A8 and D8 candidates are positive roots
of its finite local E8 system. This directly checks the worker's
simultaneous-collision dominance decision.

**Conclusion: PASS.** The original A8 vertex and this E8 vertex are
connected by an actual edge of the same full reflection chamber. This
is stronger than an auxiliary E8 basis calculation, but settles neither
the other edges nor finiteness of the whole chamber.
