# A cap of the old7 escape on the exact 64-wall seed

An exact all-right cap V76 extends the frozen64-wall H9 chamber to a
**65-facet chamber with exactly203 finite vertices**, with no original
wall deletion. It preserves every one of the201 old finite vertices and
cuts off the prior strict spacelike escape near7/30/47. A different strict
spacelike escape near30/59/75 survives, so this chamber still has infinite
volume.

Numerically, the whole old7 face of the fixed-six-white H3 section becomes
a compact11-gon. The already exact compact old60 section triangle persists.
The complete numerical section cones still have ten open ends; the old7
closure is not an ambient finite-volume or ambient-facet compactness claim.

The [producer](../scripts/fresh_exact64_escape7_cap_20261018.py),
[certificate](fresh_exact64_escape7_cap_20261018.json) and
[standalone checker](../scripts/check_exact64_escape7_cap_saved_20261018.py)
bind the frozen [exact64-wall seed](fresh_remote60_exact_original_extension_20261018.md).
The [numerical remote bridge](fresh_remote60_compact_face_bridge_20261018.md)
supplies only the existing numerical coordinate realization used for section
diagnostics. The new exact V certificate uses the physical roots and metric
of the exact64-wall seed, independently of any five-wall-deletion chain.

## A complete small duty domain

The actual old7/30/47 exterior region has native bonds

```
7–30 : 5
7–47 : 3
30–47: dotted.
```

Require the two genuine new finite corners W+{7,30,V} and W+{7,47,V},
where W={0,1,2,3,5,6}. No finite30/47/V corner is imposed: its old
30/47 edge is dotted. The two new triangles7/30/V and7/47/V have
boundary7–30–V–47–7. Their shared7–V edge has exactly those two
incident triangles. This describes the local cap, not the whole chamber.

W+{7,30} is H4⊔H4: the fixed four white roots form one H4, and the
chain5–6–7–30 has labels3,3,5. A finite connected Coxeter diagram cannot
extend an induced H4 to rank5. Thus finiteness of the first new corner
forces V to be orthogonal to both H4 factors, hence to all six whites,
7 and30. This derives profile0; it is not a selected-profile assumption.

The chain5–6–7–47 is A4. After the forced orthogonalities, the second
corner can only be A4+A1, A5 or B5, corresponding to V–47 labels2,3,4.
A label5 would require H5, which is not finite; larger labels also cannot
extend the chain to a finite diagram. The fixed other H4 is orthogonal.
Therefore the entire two-corner duty domain consists of just

```
frame [7,30,47], profile0, labels[2,2,m], m=2,3,4.
```

This is finite-type classification, not an arbitrary low-label cutoff or
an exhaustion of other cap incidences. An independent agent confirmed this
conditional finite-duty argument; the original exact W and bond data are
bound inputs.

All three recipes have two real norm roots. Their center-facing roots are
numerically native against all64 old facets, including U75. The opposite
roots fail the common-center orientation and have positive old-wall
products. Every root, including every rejection, is saved. There are no
duplicates of an old wall.

## Exact all-right member

Only m=2 is promoted to an exact wall here. Let bars denote projection off
W. Let F=(bar(r7),bar(r30),bar(r47)), with Gram matrix G, and let
c=(bar(h)·F_i). The source certifies that G has two positive leading
pivots and negative determinant, hence Lorentz signature(2,1).

Set

```
r = bar(h) − sum_i(c G^-1)_i F_i
q = r² >0
V = −sqrt(2/q) r.
```

The one-dimensional orthogonal complement is positive. This formula
therefore gives V²=2 and exact orthogonality to W,7,30,47. The other sign
has positive product with h; the displayed sign has negative product.
These are exactly the two all-right norm solutions, so there is no sheet
selection inferred from a rounded normal.

Directed rational intervals verify all64 native products: nine are the
prescribed zero products with W,7,30,47, and the other55 are strictly
below−2. The tightest dotted upper is approximately−2.0405993033 at old10.
The proof uses the same certified-original native isolation principle as
the exact seed: only old products whose interval isolates an exact native
identity are refined. Both the Schur formula and direct physical-coordinate
pairing enclose the same new product and are intersected. The source saves
every identity, matrix, coefficient, direct product and strict gate.

The exact new norm follows from the formula, not merely from a numerical
norm check. A direct interval norm containing2 is an additional consistency
check. The two other m=3,4 members remain numerical in this packet.

## All facets and the complete 203-vertex census

The exact seed's strict center h remains inside V. As in the preceding
proof, the old facet witnesses

```
h − (h·r_i)r_i/2
```

have strictly negative products with V. The new witness
`h−(h·V)V/2` is timelike and has strictly negative products with all64 old
walls. Hence V is an actual new facet and all old facets remain actual.
These65 strict witness tests are saved and checked.

For every old finite support I, the positive-definite Coxeter Gram G_I
has entrywise nonnegative inverse. Its projection coefficients
`G_I^-1(h·r_I)` are nonpositive. Since V·r_i≤0, the corresponding old
finite vertex has strictly negative product with V. All201 old finite
vertices therefore persist without repeating their large Gram inversions.

Every new finite vertex contains V and can use only its nine nondotted
old neighbors W∪{7,30,47}. The acute finite active-root independence
argument from the exact seed applies unchanged: a singular positive
semidefinite component would give a nonnegative root relation contradicting
the strict center. Thus a finite H9 vertex has exactly nine active walls.
There are nine possible choices of eight old neighbors. Seven contain the
dotted pair30/47 and are impossible. The two remaining supports are

```
W + {7,30,V}
W + {7,47,V}.
```

Both are checked positive definite. Their actual physical vertex rays are
timelike, future-oriented and strictly inside every remaining wall. They
are new because all201 old vertices lie strictly inside V. This proves
the complete finite census201+2=203.

The original63-wall facets, the added U75, and all original finite and
compact-face certificates are retained. There is no new deletion debit.
In particular the three exact vertices of the old60 section triangle
persist; V cannot cut their convex hull. That section triangle remains
unchanged and compact.

## Complete numerical section comparisons

The numerical producer independently enumerates the full projected cone
for the exact64-wall seed and for each of the three numerical caps, using
all58 or59 projected wall inequalities. The six fixed white walls define
the section; they are not omitted ambient native gates.

|Section cone|Finite rays|Exterior rays|Open ends|Old7 face|
|---|---:|---:|---:|---|
|Exact64 seed, numerical coordinates|28|6|10|9 finite vertices and1 exterior ray|
|m=2 cap|30|6|10|11 finite vertices|
|m=3 cap|30|6|10|11 finite vertices|
|m=4 cap|30|6|10|11 finite vertices|

All cones are simple within the stated numerical tolerances and have no
reported warnings. The full old7 and old60 ray ledgers are saved, together
with every full-cone ray, active triple and edge. For m=2,3, the exterior
supports7/30/47 and30/47/52 disappear; new supports30/52/V and47/52/V
appear. For m=4, only7/30/47 disappears and30/47/V appears. The all-right
member is the simpler useful exact choice.

Thus an additional whole old section face closes, but the exterior
obstruction has not disappeared. The compact11-gon claim for old7 is a
numerical section result here, not a directed-interval full-face theorem.
No claim is made that the eight-dimensional ambient facet7 is compact.

## Exact remaining escape

The old strict escape near7/30/47 has new V product in approximately
[0.0011027117,0.0011027469], strictly positive, so it is genuinely removed.

A different exact escape is constructed from the retained frame30/59/75.
Project bar(h) off that frame, reverse the resulting spacelike direction,
and add h/4. The certificate checks its positive square, negative product
with h, and strictly negative products with every one of the65 walls.
Its V product is approximately−1,494,885.026. As in the exact seed, the
segment from the strict timelike center to this strict spacelike point
crosses the future light cone inside all inequalities. The chamber has
an open ideal patch and therefore infinite volume.

This remaining region, rather than the removed7/30/47 escape, is an
explicit obstruction for the next construction. No finite-volume example
or new dimension theorem follows from this packet.

## Verification

Final producer generation session39348, terminal62f8d4, log732b34, passed
the three recipes/six norm roots, four complete numerical section cones,
all64 exact native gates, all65 facet witnesses, complete203 finite census,
and the remaining strict escape. Separate producer replay session15791,
terminal9a5939, logcb763f passed. The hardened standalone checker session91970,
terminalfa668b passed, including explicit original-white, frame-profile and
finite-duty guards. Its subsequent docstring clarification does not change
the verification algorithm or certificate.

The standalone checker has no project imports. Its physical calculation
uses the same directed interval identities in a separate saved-verification
implementation. It independently reconstructs all four numerical cones
by cofactor normals for every triple of constraints, and checks all ray
signs, active triples, edges and section ledgers. This verifies numerical
cone coverage; it does not turn numerical cones into interval certificates.
The original scaffold, exact64 seed and complete201-vertex theorem remain
bound inputs. No earlier five-wall-deletion bank or historical finite
census is rerun.

The root agent completed a full mathematical source/report/standalone-helper
read. It checked finite-duty H4 isolation, the exact Schur normal and sign,
native gates, all65 facet witnesses, the complete201+2 finite-vertex count,
preservation of the old60 section triangle, and the strict30/59/75 escape
proving infinite volume. It also checked that the numerical old7 eleven-gon
is distinguished from exact and ambient claims. This was a read-only review,
without execution or independent interval replay. The original scaffold,
old201 census and certified-native isolation remain inherited inputs.
