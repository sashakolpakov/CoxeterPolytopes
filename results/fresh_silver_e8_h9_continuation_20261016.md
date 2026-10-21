# An infinite-order chamber symmetry for the silver E8 lattice

The literal lattice

\[
 L=\langle-(1+\sqrt2)\rangle\perp E_8\perp\langle1\rangle,
 \qquad O=\mathbf Z[\sqrt2],
\]

has an exact infinite-order stabilizer of its full integral reflection
chamber. The resulting nonreflectivity statement is an internally checked
computer-assisted research proof candidate; independent review and external
validation/priority remain separate. It excludes a finite-volume Coxeter
polytope whose wall reflections all preserve this literal lattice. It does
not exclude arbitrary nonintegral geometric walls, other lattices, or all
lattices commensurable with this one. No finite-volume example is produced.

All four new producers have passed generation and saved replay. The
[symmetry certificate](fresh_silver_e8_h9_continuation_20261016_symmetry.json)
contains the literal matrices, their integral inverses, outward normal maps,
three actual finite corners, basepoint decompositions, all 48 facet witnesses,
and a strict spacelike escape. Its
[checker](../scripts/fresh_silver_e8_h9_continuation_20261016_symmetry.py)
reconstructs the maps directly from the stored roots rather than trusting a
numerical map search.

## Complete bounded prefix and finite geometry

The [height-100 producer](../scripts/fresh_silver_e8_h9_continuation_20261016.py)
and [certificate](fresh_silver_e8_h9_continuation_20261016.json) extend the
unchanged [height-30 prefix](fresh_silver_e8_h9_20261015_h30.json). All three
primitive norm classes are included. There are 29 time controls, 31,423
covector recursion visits, 75 candidates, and 48 accepted walls. The first
17 controls and 20 walls agree exactly with the old certificate.

New walls 22 and 41 cut the old spacelike witness. A different strict
spacelike escape survives all 48 walls; its exact vector and 48 strict signs
are in the symmetry certificate. Thus the larger finite prefix still has
infinite volume independently of the full-chamber symmetry argument.

The separate complete finite-corner checks give:

| Prefix | Pairwise finite nine-subsets | Finite vertices | Finite edges | Open finite germs |
| --- | ---: | ---: | ---: | ---: |
| Height 30, 20 walls | 68 | 24 | 40 | 136 |
| Height 100, 48 walls | 744 | 100 | 182 | 536 |

The [20-wall checker](../scripts/fresh_silver_e8_h9_continuation_20261016_corners.py)
and [data](fresh_silver_e8_h9_continuation_20261016_corners.json), and the
[48-wall checker](../scripts/fresh_silver_e8_h9_continuation_20261016_h100_corners.py)
and [data](fresh_silver_e8_h9_continuation_20261016_h100_corners.json), enumerate
every nine-clique in the finite-pair graph, test positive definiteness
exactly, solve its point, and check every wall inequality. A finite vertex
of an acute chamber is simple: its incident normals lie in a positive
space, and a dependence among their nonpositive-off-diagonal Gram blocks
would contradict the common strict interior point. Hence these nine-subset
tests cover all finite vertices. Two such vertices sharing eight walls give
an actual finite edge. No complete ideal/hyperideal cone census is claimed.

## Why the root alphabet and bounded enumeration are complete

Write \(s=\sqrt2\), \(\epsilon=1+s\), and \(\eta=2+s\).
The ring \(O\) is norm-Euclidean: rounding the two coefficients gives errors
\(a,b\in[-1/2,1/2]\) with \(|a^2-2b^2|\le1/2<1\). Thus a primitive vector
has a coordinate Bezout relation. Since \(\det G=-\epsilon\) is a unit,
reflection integrality \(2Gr/q\in O^{10}\), where \(q=r^2\), implies
\(q\mid2\): multiply by the integral inverse of \(G\), then use that Bezout
relation. The conjugate form is positive definite, so \(q\) is totally
positive.

The units are \(\pm\epsilon^k\). One elementary verification is that no unit
lies strictly between 1 and \(\epsilon\): if its norm is \(\pm1\), its
integer \(s\)-coefficient \((u-u')/(2s)\) lies strictly between 0 and 1.
Reducing an arbitrary positive unit by powers of \(\epsilon\) gives the
claim. Totally positive units are squares. The prime over 2 is generated
by \(s\), and its totally positive associate is \(\eta=\epsilon s\).
Consequently the primitive reflection directions have representatives with
norm \(1,\eta,2\). Nonprimitive parallel representatives do not add mirrors.

The initial basepoint is \(H=(1,0^9)\), with zero-height sector \(E_8+A_1\).
The finite integral Weyl group moves every zero-height root into the
dominant chamber. The exact dominant census in the three norm classes has
counts \(1,0,2\): the scalar vectors 1 and \(s\), which give the same mirror,
and the ordinary E8 root with coordinates \((2,4,6,5,4,3,2,3)\).
Its orbit is the ordinary E8 root system. Thus the initial nine simple
walls exhaust the zero-height reflection arrangement.

For a positive-height root \(r=(x,y,z)\), put \(\lambda=G_{E_8}y\). The
dominant inequalities are \(\lambda_i\ge0\), \(z\ge0\), and

\[
 \lambda^TG_{E_8}^{-1}\lambda+z^2=q+\epsilon x^2,
 \qquad h=\epsilon x^2/q.
\]

The physical height bound and conjugate positivity give simultaneous finite
bounds on \(x,x'\). Dual Cauchy gives
\(\lambda_i^2\le2(q+\epsilon x^2)\) in both embeddings and the corresponding
scalar bound. Exact reflection-integrality filters are applied to every
coordinate. The inverse of each leading E8 principal block gives the
unconstrained minimum energy for a covector prefix, valid in both places;
the positive entries of \(G_{E_8}^{-1}\) give an additional physical prefix
bound. The final norm equality and every chronological acuteness condition
are checked exactly. Equal-height packets are checked mutually acute.
These are complete finite bounds, not a coefficient cutoff. The usual
Vinberg selection therefore supplies actual initial walls of the full
integral reflection chamber, including all three root norms.

## Two literal corner maps

Index the stored 48 normals by \(N_0,\ldots,N_{47}\), and take the ordered
source sector \((0,1,\ldots,8)\). Its first eight roots have norm 2 and its
last root norm 1. Two target ordered sectors are

\[
 (4,3,2,1,0,17,11,7,8),\qquad
 (0,1,2,41,6,5,4,7,8).
\]

Their actual finite points, normalized to time coordinate 1, are vertices
19 and 41 in the complete census. Both ordered Grams equal the source
Gram. The certificate constructs \(T_1,T_2\) carrying every source normal
to the corresponding target normal and carrying \(H\) to respectively
\((11+8s)V_{19}\) and \((29+20s)V_{41}\). It verifies directly

\[
 T_i^TGT_i=G,\qquad T_i,T_i^{-1}\in M_{10}(O),
\]

as well as the positive corner scales. These are literal outward-sector
maps preserving the future sheet, not merely diagram isomorphisms.

Let \(U=T_2T_1\) and \(T=U^2\). The exact certificate gives

\[
 \operatorname{tr}U=4+s,\qquad
 \det T=1,\qquad \operatorname{tr}T=10+4s>10.
\]

A finite-order ten-dimensional matrix has all eigenvalues on the unit
circle and hence trace of absolute value at most 10. Thus \(T\) has infinite
order. No exhaustive automorphism-group search is needed for this witness.

## Why the maps stabilize the full chamber

Let \(V\) be any one of the three displayed corners and \(S\) its complete
nine-root outward support. Its Gram \(G_S\) is positive definite with
nonpositive off-diagonal entries, so \(G_S^{-1}\) is entrywise nonnegative.
The certificate checks the exact decomposition

\[
 H=aV+\sum_{i\in S}b_iN_i,\qquad a>0,\quad
 b=G_S^{-1}(\langle H,N_i\rangle)_{i\in S}\le0.
\]

Any later simple Vinberg wall \(r\) is acute to every retained wall and has
\(\langle H,r\rangle<0\). Consequently

\[
 0>\langle H,r\rangle
 =a\langle V,r\rangle+\sum b_i\langle N_i,r\rangle
 \ge a\langle V,r\rangle.
\]

Every displayed corner and its entire local sector therefore persist in
the full chamber. Integral future isometries normalize the complete
integral-reflector arrangement. Its local finiteness and the exact positive
outward-sector maps show that each \(T_i\) carries this full chamber to
itself; hence so does the infinite-order matrix \(T\).

The integral isometry group is discrete: the other real embedding is
positive definite and hence bounded, so bounded physical matrices have
bounded entries at both embeddings and only finitely many integral
possibilities. An infinite chamber stabilizer therefore provides infinitely
many disjoint isometric translates of a sufficiently small interior ball.
The full integral reflection chamber has infinite volume. A finite-volume
reflection group generated by lattice-preserving wall reflections would
be a subgroup of this full reflection group, forcing the latter to have
finite covolume as well, a contradiction.

## Replay and limits

```sh
env PYTHONPATH=src:scripts python3 scripts/fresh_silver_e8_h9_continuation_20261016.py --verify
env PYTHONPATH=src:scripts python3 scripts/fresh_silver_e8_h9_continuation_20261016_corners.py --verify
env PYTHONPATH=src:scripts python3 scripts/fresh_silver_e8_h9_continuation_20261016_h100_corners.py --verify
env PYTHONPATH=src:scripts python3 scripts/fresh_silver_e8_h9_continuation_20261016_symmetry.py --verify
```

All four saved replays pass. The source-bound exact checks and the proofs
above have distinct roles: finite arithmetic alone does not establish the
universal norm alphabet or the corner-persistence implication. An independent
arithmetic/proof review is being prepared separately. Prior frozen artifacts
were only read; no height beyond 100 was searched. A focused provenance
check found no applicable literal-host exclusion, but neither originality
nor an exhaustive literature classification is asserted.
