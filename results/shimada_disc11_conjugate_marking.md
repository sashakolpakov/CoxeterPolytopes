# A Nonempty Conjugate Eigen-Ample Section on the Discriminant-11 K3

Status: an exact integral conjugate of the (H_{11}) marking has an identity
eigenspace meeting the ample cone.  This passes the geometric nonemptiness
gate.  It does **not** yet give a complete Coxeter polytope or decide finite
volume or compactness.

## Exact chamber-interior certificate

Shimada acts on Néron--Severi row vectors.  The archived integral matrix (T)
satisfies

\[
\det T=1,\qquad TGT^{\mathsf T}=G,
\]

and acts trivially on the cyclic discriminant group of order 11.  If (M) is
the original column action of (phi), the conjugated action is

\[
M_T=T^{\mathsf T}MT^{-\mathsf T}.
\]

It remains integral, self-adjoint, and satisfies (M_T^2=M_T+I).  Triviality
on the discriminant group also lets the Néron--Severi isometry extend by the
identity on the transcendental lattice; this is a genuine change of K3
marking, not a change of the abstract divisor lattice.

Let (K_0) be the ten-row identity-eigenbasis of the original action in
Shimada coordinates.  For

\[
\begin{aligned}
y={}&(-54+26\phi,-1,3-3\phi,3-3\phi,5-5\phi,7-7\phi,\\
&11-11\phi,12-12\phi,15-15\phi,16-16\phi),
\end{aligned}
\]

put (x=yK_0T).  Exact calculation gives

\[
xGx^{\mathsf T}=-2474+2688\phi>0.
\]

Shimada chamber (D_{68}) has 28 archived defining wall functionals.  Every
pairing (xv^{\mathsf T}) is strictly positive, and their exact minimum is
1.  Hence (x) lies in the interior of (D_{68}), so the transported
identity eigenspace meets the ample cone.

The full matrix, sparse wall archive, witness generator, and exact checks are
in `scripts/shimada_disc11_conjugate.py`; ordinary tests do not require the
47 MB external Shimada file.

## Immediate geometric wall filter

Precisely wall orbits 0 through 8 of (D_{68}) have Shimada's
`_minustwo` adjacency tag.  They contain 17 actual nef-boundary roots.  From
a dual wall row (v), the primal root is recovered as (r=vG^{-1}).  All 17
are integral and have square (-2).

For a physical divisor row (r), its coordinates in the conjugated
(mathcal O_F)-structure are obtained by pulling back (rT^{-1}).  The
five-shell and reflection-divisibility filters then give:

- all 17 local roots are positive at both real places;
- all 17 define integral (mathcal O_F)-reflections;
- nine have norm (2), and eight have norm (4-2\phi);
- their field span has rank 9; and
- their unordered pair signs are 30 negative, 98 zero, and 8 positive.

There is a sharper exact local incidence statement.  The two members of each
root orbit 1 through 8 restrict to positive proportional functionals, with
ratio (phi) or (phi-1).  Eight archived positive-norm face witnesses
prove that these are eight genuine facets of (D_{68}\cap E).  The remaining
orbit-0 root restriction is an exact positive combination of other restricted
walls and is redundant.  Thus 16 ambient root walls give eight local
restricted root facets.  This is only a local induced-chamber result; crossing
non-root Borcherds walls is required to assemble the global eigen-nef cell.

## Complete local eigensection and first-star traversal

The complete restriction of all 28 primitive walls is now exact.  Positive
proportionality leaves 19 distinct half-spaces, and positive exact dependence
certificates reduce these to precisely 11 irredundant facets:

\[
8\text{ effective }(-2)\text{ facets}
\quad+\quad
3\text{ internal Borcherds facets}.
\]

The three internal wall classes are

\[
\{(9,0),(16,0)\},\qquad \{(11,1)\},\qquad \{(14,0)\}.
\]

The first class consists of two different ambient walls whose restrictions
differ by the positive scalar (phi).  Eleven exact timelike face points prove
that every listed class is a facet; eight positive (mathbf Z[phi])-linear
dependences prove that every other distinct restricted inequality is
redundant.  This is a complete local statement, with no floating-point facet
decision.

Shimada's full 1,098-chamber atlas also gives an exact first-star traversal.
Crossing (D_{68}) wall ((11,1)) enters a physical copy of representative
(D_{86}); an explicit timelike point is strictly positive on all 35 of its
primitive walls.  Crossing wall ((14,0)) similarly enters a physical copy of
(D_{181}), with strict positivity on all 26 primitive walls.

The simultaneous class (\{(9,0),(16,0)\}) has a more revealing local fan.
The fixed eigenspace passes through the representative sequence

\[
D_{68}\xrightarrow{(9,0)}D_{10}
\xrightarrow{(12,0)}D_{10}
\xrightarrow{(8,1)}D_{69}
\xrightarrow{(15,0)}D_{68}g.
\]

All three intermediate physical chambers meet the eigenspace only along the
same timelike face.  An exact point on the opposite side is strictly interior
to (D_{68}g).  The resulting integral isometry (g) has determinant one,
preserves the Neron--Severi form, is an involution, and centralizes the
conjugated golden action.  On the ten-dimensional eigenspace it has fixed
rank nine and ((-1))-rank one, and it negates the internal wall normal.
Thus this internal wall is the mirror of an actual automorphism of the
enhanced K3 object ((X,tau)); it is not a divisor boundary.

The physical (D_{181}) section has now been reduced exactly. Its 26 primitive
walls give 17 positive-proportional half-spaces and precisely eleven facets:
nine paired root facets and two internal facets. The latter lead to a physical
copy of (D_{30}) and back to the original (D_{68}), with exact strict interior
witnesses on both sides. Six positive (\mathbf Q(\phi))-linear identities
remove every other restricted half-space. See
`results/shimada_disc11_d181_section.md`. The physical (D_{30}) section and
the (D_{86}) branch remain to be reduced.

The lazy external-atlas parser is
`scripts/shimada_disc11_atlas.py`, and the exact traversal is reproduced by
`scripts/analyze_shimada_disc11_eigen_nef.py` and
`scripts/analyze_shimada_disc11_d181.py`. Global root-wall completeness has
not been claimed.

## Global published curve-representative filter

Refiltering all 347 effective ((-2))-curve representatives in Shimada's set
`B` under the conjugated field structure gives 18 integral field-reflection
roots: nine of norm (2) and nine of norm (4-2\phi).  They have field span
rank 10.  Their zero-based source indices are

```text
0, 2, 3, 4, 6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 18, 19, 28, 31
```

The pair-sign distribution is 31 negative, 114 zero, and 8 positive.  All 347
effective curve pairings with the ample witness are positive, with minimum 1.
The 17 local field vectors occur among these 18; the extra vector, at `B[17]`,
has coordinates

\[
(-2+\phi,1-\phi,0^8)
\]

and completes the field span to rank 10.

The 347 vectors are curve-orbit representatives, not the complete effective
wall set around this eigenspace.  Nor do the 18 survivors form a certified
outward-oriented simple-root system.  The next exact task is to traverse every
Shimada induced chamber met by the transported eigenspace, accumulate the
actual restricted ((-2))-facets, and certify cell completeness before testing
Coxeter labels, volume, or compactness.

## Reproduction

```sh
PYTHONPATH=src:. python3 -m pytest tests/test_shimada_disc11_conjugate.py
PYTHONPATH=src:. python3 scripts/shimada_disc11_conjugate.py /path/to/disc11.txt
```

The external-data checksum and concise results are archived in
`data/candidates/shimada_disc11_conjugate_marking.json`.
