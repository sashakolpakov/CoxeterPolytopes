# Exact reconstruction of Vinberg's discriminant-three chamber in H^19

## Outcome

The target

\[
S=U\oplus 2E_8(-1)\oplus A_2(-1)
\]

does produce a finite-volume hyperbolic Coxeter polytope in dimension 19.
It has **36 facets, 8,346 finite vertices, and 107 cusps**.  This is not a
new example: it is the chamber of Vinberg's most algebraic K3 surface
\(X_3\), reconstructed exactly from the combinatorial description in
Looijenga's account of Vinberg's calculation.

The useful new repository result is an independently replayable, exact
Gram matrix and a literal Vinberg edge/cusp certificate.  No automorphism of
the K3 surface is assumed or used.

## The 36 walls

Let

\[
B'=\{a,b,c\},\qquad B''=\{u,v,w\}.
\]

The 24 long walls have square \(-2\).  Six are indexed by the elements of
\(B'\sqcup B''\).  The remaining 18 are indexed by the two orientations of
each edge of \(K_{3,3}\).  Each edge \(xy\), with \(x\in B'\) and
\(y\in B''\), is replaced by the path

\[
x\;--\;(x,y)\;--\;(y,x)\;--\;y.
\]

Adjacent long walls have product 1; all other distinct long walls are
orthogonal.

The 12 short walls are indexed by the six bijections \(B'\to B''\) and the
six bijections \(B''\to B'\).  We use their integral primitive normals,
which have square \(-6\) and divisibility 3.  The wall for a bijection
\(\sigma\) has product 3 with the three long walls \((b,\sigma(b))\) and is
orthogonal to the other long walls.  Scaling Looijenga's dual-lattice short
roots by 3 turns his short-short products

\[
-\frac23,\quad-\frac43,\quad-\frac53
\]

into the integral wall products

\[
6,\quad12,\quad15.
\]

The resulting 36-by-36 matrix has

- rank 20 and exact inertia \((1,19,16)\);
- nonzero Smith invariants \(1^{19},3\);
- hence discriminant 3 and the full lattice
  \(U\oplus2E_8(-1)\oplus A_2(-1)\), rather than a proper finite-index
  sublattice.

Its complete Coxeter pair census is:

| relation | pair count |
|---|---:|
| right angle | 501 |
| \(m=3\) | 27 |
| \(m=6\) | 36 |
| parallel | 24 |
| ultraparallel | 42 |

Among the ultraparallel pairs, 30 have normalized squared product 4 and 12
have normalized squared product \(25/4\).

## Exact Vinberg finite-volume audit

Vinberg's completed root algorithm supplies the completeness statement: the
36 roots are the simple roots, and every root is an all-nonnegative or an
all-nonpositive integral combination of them.  Looijenga records this as
Properties 1.4.  The local verifier then applies the edge-endpoint and full
parabolic-neighborhood tests to that complete list.

There are 87,046 elliptic rank-18 edge diagrams:

- 25,126 use only long walls;
- 61,920 use one short wall.

There are 8,346 elliptic rank-19 finite vertices:

- 2,046 use only long walls;
- 6,300 use one short wall.

The connected affine components are generated exhaustively by connected
elliptic growth.  A branch is stopped as soon as it becomes affine or
indefinite; this is exhaustive because every proper subdiagram of an
irreducible affine Coxeter diagram is elliptic.  Taking the full zero
neighborhood of every component produces exactly 107 maximal rank-18
parabolic neighborhoods:

| cusp type | count |
|---|---:|
| \(\widetilde E_8+\widetilde E_8+\widetilde G_2\) | 36 |
| \(\widetilde D_{16}+\widetilde G_2\) | 36 |
| \(\widetilde E_7+\widetilde D_{10}+\widetilde A_1^{s}\) | 18 |
| \(\widetilde D_7+\widetilde A_{11}\) | 9 |
| \(\widetilde A_{17}+\widetilde A_1^{s}\) | 6 |
| \(3\widetilde E_6\) | 2 |

Finally, every one of the 87,046 elliptic edges has exactly two endpoints:

| endpoints of an edge | edge count |
|---|---:|
| two finite | 72,015 |
| one finite, one ideal | 14,544 |
| two ideal | 487 |

There are no bad edges and no endpoint-generated nonedges.  Vinberg's
criterion therefore proves that the chamber has **finite volume**.  Its 107
maximal parabolic neighborhoods prove that it is **noncompact**.

An independent CoxIter run on the same 36-node diagram reports finite
covolume, noncompactness, 107 ideal vertices, and

\[
\begin{aligned}
(f_0,\ldots,f_{19})={}&(8453,87046,429450,1375716,3219626,5844660,\\
&8516190,10157719,10019112,8200968,5562014,3107454,\\
&1415874,518320,149208,32742,5228,564,36,1).
\end{aligned}
\]

Here (f_0=8346+107), exactly matching the independently enumerated finite
and ideal vertices.  The alternating sum is 2, as expected for the closed
face lattice of a 19-dimensional polytope.

## Relation to the Borcherds/Conway route

The primitive embedding

\[
S\hookrightarrow II_{1,25}=U\oplus3E_8(-1)
\]

has complement \(E_6(-1)\).  This explains why the induced-wall shell has
exactly the two expected integral types: long \((-2)\) roots and the
integral multiples of short dual roots, namely \((-6)\), divisibility-3
normals.  For this target, however, a new Conway-coordinate enumeration is
unnecessary: Vinberg already completed the intrinsic algorithm, and
Looijenga gives the complete 36-wall combinatorics needed to reconstruct the
same chamber.

Thus this target is an excellent high-dimensional regression case for the
algebro-geometric fishing pipeline, but it does not move the known-example
boundary.

## Reproduction and sources

- Exact verifier: [`scripts/certify_vinberg_x3_h19.py`](../scripts/certify_vinberg_x3_h19.py)
- Frozen matrix and certificate:
  [`data/enumerations/vinberg_x3_h19.json`](../data/enumerations/vinberg_x3_h19.json)
- E. B. Vinberg, *The Two Most Algebraic K3 Surfaces*, Math. Ann. 265
  (1983), 1-21, [DOI](https://doi.org/10.1007/BF01456933).
- E. Looijenga, *The period map for cubic fourfolds*, Invent. Math. 177
  (2009), especially pp. 5-7 and Figure 1,
  [DOI](https://doi.org/10.1007/s00222-009-0178-6),
  [author PDF](https://www.staff.science.uu.nl/~looij101/viervoud4a.pdf).

Replay with:

```bash
PYTHONPATH=src:. python3 scripts/certify_vinberg_x3_h19.py --check-archive
```
