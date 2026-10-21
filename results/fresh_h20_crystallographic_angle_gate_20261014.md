# A complete finite-angle restriction for sections of the Borcherds chamber

**Result.** If an H20 hyperplane section of the saved210-wall Borcherds H21
chamber is a finite-volume Coxeter polytope, every actual finite dihedral
angle has label2,3,4 or6. This covers arbitrary real section normals and
arbitrary contacts on other walls. It does not establish arithmeticity,
exclude all sections, or produce a new polytope.

The [exact checker](../scripts/fresh_h20_crystallographic_angle_gate_20261014.py)
and [certificate](fresh_h20_crystallographic_angle_gate_20261014.json)
reconstruct the literal spherical-triple catalogs and check the new
uniform inequalities. The proof uses the previously verified
[all-H3 obstruction](fresh_h20_coupled_h3_obstruction_20261014.md),
[short-A2 large-label table](fresh_h20_dihedral_edge_sections_20261014.md),
and [label5 corollary](fresh_h20_dihedral_edge_sections_20261014_m5.md).

## Compact face and unchanged-root obstruction

Suppose a section has an actual I2(m) ridge with m>=5 and m!=6, bounded
by projected source roots a,b. Its dimension18 face F is compact: no
affine Coxeter diagram contains that dihedral diagram. For m>=7 it has
no proper connected finite extension. For m=5 its only connected
three-node finite extension is H3, already excluded for every section of
this source. Thus every facet of F is orthogonal to the projected a,b.
Consequently F is itself a compact Coxeter polytope.

If all its facet normals are unchanged source roots, compactness makes
them span its Lorentzian space. That space is rational and carries the
full invariant lattice obtained by intersecting it with Z^22. Its compact
reflection lattice would have finite index in the integral orthogonal
group; Meyer isotropy makes the latter nonuniform. This contradiction
also applies to compact subfaces of F bounded only by unchanged source
roots. The source's integral reflection matrices are checked explicitly.
The previous proofs give the full invariant-lattice argument; rational
coordinates alone are not the premise.

Write e for the section normal, q=e²>0, and a',b' for the projections.
For a facet l of F, the three projected normals a',b',l' have a positive
definite Gram matrix. Their old Gram is that matrix plus a positive
semidefinite rank-one matrix, so it too is positive definite. The exact
source catalog has **no positive-definite three-node triangle**.

## Two nonzero contacts on an old nonorthogonal pair

Suppose (e,a) and (e,b) are both nonzero and (a,b)<0. If (e,l) were
nonzero, orthogonality to a',b' would imply

```
(l,a)=(e,l)(e,a)/q !=0,
(l,b)=(e,l)(e,b)/q !=0.
```

The old three-node Gram would be a triangle, contradicting the preceding
catalog. Thus every boundary l is unchanged, and F is impossible.
This requires no bound or rationality assumption on the contact ratio.

## Two nonzero contacts on an old orthogonal pair

Normalize (e,a)=1, put (e,b)=r and take a short root first if lengths
differ. The projected pairing is -r/q, so r>0. A modified facet of F
must satisfy (l,b)=r(l,a). The literal positive-definite three-node
catalog leaves exactly these possibilities:

| Old endpoint norms | r | Modified boundary norm and (l,a) |
|---|---:|---|
| 2,2 | 1 | 2, -1 |
| 2,4 | 1 | 4, -2 |
| 2,4 | 2 | 2, -1 |

Two long old roots are never orthogonal. If the section has any other
ratio, all boundaries are unchanged and the compact-face contradiction
already applies.

For a modified wall l put k=(l,a). Orthogonality gives (e,l)=qk, so its
normal on F is s_l=l-k e, with

```
s_l²=d_l-qk²,       (s_l,s_j)=(l,j)-qk_l k_j.
```

For the short/short case, the new angle equation is
cos(pi/m)=1/(2q-1), hence 1<q<=sqrt(5)/2<6/5. A modified short wall
has squared norm2-q. Its only possible nonright finite neighbor would
be an unchanged short root with old pairing-1, giving squared cosine
1/[2(2-q)] strictly between1/2 and5/8<cos²(pi/5). This is a forbidden
Coxeter interval. Contacts to unchanged long roots are divergent unless
orthogonal; two modified short roots are divergent as well. A modified
facet therefore has only unchanged perpendicular boundaries and is an
impossible compact rational reflection polytope of dimension17.

For the mixed-length case r=1, positivity of the projected pair requires
q>3/4. Modified boundaries have norm4-4q; if they exist q<1. All unchanged
boundaries are short, since they must be orthogonal to the old long root.
A nonright modified/unchanged pair would have squared cosine
1/[2(1-q)]>2; two modified walls are divergent. Again a modified facet
has only unchanged perpendicular boundaries.

For r=2, positivity requires q>3/2. Modified boundaries have norm2-q,
so q<2 if any exists. All unchanged boundaries are short. Their
nonright contact with a modified wall has squared cosine
1/[2(2-q)]>1, and modified/modified pairs are divergent. The same
compact rational facet contradiction applies. If no modified wall
exists in either mixed-length case, use F itself instead.

## One nonzero contact

An old orthogonal pair remains orthogonal if only one endpoint changes,
so cannot produce the specified m. Normalize (e,a)=1, (e,b)=0 and put
C=cos²(pi/m). The literal catalogs of possible F boundaries are:

| Old pair (a²,b²,(a,b)) | q | (l²,(l,a)) for a finite boundary |
|---|---|---|
| (2,2,-1) | 2C/(4C-1) | (2,0),(4,0),(2,-1),(4,-2) |
| (4,4,-2) | C/(4C-1) | (2,0),(2,-2) |
| (2,4,-2) | C/(2C-1) | (2,0),(2,-1) |
| (4,2,-2) | C/(4C-2) | (2,0),(4,-2) |

These exhaust old positive-definite pairs and triples in the source.
No contact condition on any fourth wall has been used.

The first row is exactly the previously checked short-A2 table and
label5 argument. Their source-edge hypotheses originally supplied an
actual I2(m) face. Here that face is assumed directly; every subsequent
equation and compact-face argument is unchanged and does not require
the other18 source-edge contacts to vanish.

In the long-A2 row, for m>=7 one has 1/3<q<3/8. Every modified boundary
is short with k=-2 and norm2-4q. Two modified walls are divergent. The
only possible finite nonright contact with an unchanged short root has
squared cosine

```
1/(4-8q)=(4C-1)/(8C-4).
```

If this has label n, substitution gives
cos(2pi/m)cos(2pi/n)=1/2, with m,n>=7. The previously certified exact
monotonicity/cubic argument leaves only m=n=8. If that contact occurs,
its compact I2(8) face in F has only unchanged boundary roots: modified
walls cannot meet the chosen modified wall orthogonally. This gives an
impossible compact rational face of dimension16. Otherwise a modified
facet itself has only unchanged perpendicular boundaries. At m=5 the
exact value q=(1+sqrt(5))/8 makes all nonright contacts involving a
modified wall divergent, so the simpler facet argument applies.

For B2 with the short endpoint modified, the m=5 value
q=(2+sqrt(5))/2 makes the proposed modified boundary norm negative.
For m>=7, 1<q<4/3; modified/modified pairs are divergent, and the only
possible finite nonright squared cosine is

```
f(C)=(2C-1)/(6C-4).
```

It decreases with C. The certified bound4/5<C7<13/16 gives
5/7<f(C7)<3/4, strictly between the label5 and label6 values.
For m>=8, C>=C8>16/19 gives
1/2<f(C)<13/20<C5. Thus no integer label is possible in either case.
The rational compact-face contradiction follows.

For B2 with the long endpoint modified, the m=5 value
q=(2+sqrt(5))/4 also makes the modified boundary norm negative.
For m>=7, 1/2<q<2/3; only modified long and unchanged short boundaries
occur. Every nonright modified/unchanged or modified/modified pair is
divergent. Once again a modified facet has only unchanged perpendicular
boundaries, or every facet of F is unchanged.

Zero contacts preserve the old crystallographic angle and therefore
cannot give m>=5 with m!=6. This completes all cases on the two walls
of an arbitrary actual section ridge.

## Exact verification and scope

The checker reconstructs all allowed spherical triples directly from the
210 integer source roots, including the absence of a positive-definite
triangle. Six new pair tables are verified on whole rational intervals
by exact degree-at-most-two polynomial inequalities. It also verifies
the isolated label5 cases, the B2 gap inequalities and the long-A2
exceptional-label identity. All imported certificate hashes are checked.

```sh
python scripts/fresh_h20_crystallographic_angle_gate_20261014.py --verify
```

The conclusion restricts **actual finite ridge labels**, not distances
between disjoint facets or integrality of cyclic Gram products. In
particular it does not turn the geometric H20 search into an arithmetic
one. The remaining labels2,3,4,6 still allow searches outside the
arithmetic class; finite volume and global Coxeter conditions remain
unproved for any new example.
