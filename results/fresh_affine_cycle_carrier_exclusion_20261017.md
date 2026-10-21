# Large affine cycles are excluded from an unpaid finite quadrilateral star

**FINAL/FROZEN; independent internal mathematical/source review passes.** Let an actual
quadrilateral two-face of a finite-volume hyperbolic Coxeter polytope have
four finite vertices, and suppose all four finite Vinberg corner weights
vanish. Suppose the ambient diagram has no connected affine subdiagram
of order at most9. In the private/shared decomposition, neither H-U nor
H-V contains a connected affine A cycle of order m>=18. The affine cycle
need not be the whole private diagram U or V. This removes the cyclic part
of the larger-carrier obligation left
by the [finite-quadrilateral transfer audit](fresh_simple_finite_quad_vinberg_transfer_20261017.md).
It does not exclude large affine B/C/D subdiagrams, complete that transfer,
or establish a dimension bound.

Write W for the supporting white roots and write the opposite black pairs
as {u1,u2} and {v1,v2}. Every W+ui+vj is positive definite. Zero corner
weight makes every native cross product ui.vj zero. The four black roots
lie in the same connected component H of the star, by the actual
opposite-side indefiniteness argument in the transfer audit.

Use the private diagrams U,V from the transfer audit: U contains the two
u-roots and their private white components, and V is defined symmetrically.
Assume that A is an affine cycle contained in H-U. Any nonspherical
subdiagram of the star contains an opposite black pair, since otherwise
it lies in a positive-definite corner diagram. As H-U contains no u-root,
A contains v1,v2 and every other node of A is white. The case H-V is
symmetric. No assertion about a cycle containing three or four black
roots is used.

For any white root z outside A, both (A-v1)+z and (A-v2)+z are spherical.
The same is true for z=u1 or u2. These assertions are principal-subdiagram
consequences of the four actual finite corners.

First consider the support of z on A. If z has two neighbors a,b on the
cycle, the two old a-to-b arcs give two graph cycles through z. Each such
cycle must contain both v1 and v2: otherwise it survives one of the two
spherical deletions, but finite Coxeter diagrams are forests. The two
cycles intersect in exactly {z,a,b}, so {a,b}={v1,v2}. Applying this to
every pair also excludes three or more neighbors. Thus two contacts are
possible only at the two marked black roots.

One contact is impossible. In at least one marked deletion its old
neighbor survives. The resulting connected component has order m and is
a finite extension of the long A_(m-1) path, so its contact label is3 or4.
A simple leaf on the old cycle contains the induced affine E8 tree with
arms1,2,5; m>=18 leaves an unused arc, so there is no extra edge among
the selected vertices. A doubled leaf with its two old cycle neighbors
contains the four-node affine doubled star. Both contradict the stated
small-affine prohibition. In particular, any root with a nonzero contact
to A touches exactly v1,v2 and no white node of A.

The roots u1,u2 have zero native products with v1,v2, and therefore are
orthogonal to A. Every white node of A is also orthogonal to every white
node outside A. Since H connects both black pairs, there is a shortest
white path from a u-root to a v-root. It cannot enter a white node of A.
Let z be its last white node, adjacent to a v-root. Then z touches both
v-roots. Its preceding node y is outside A, is distinct from z, and has
a nonzero contact with z. It may itself be the initial u-root.

The [private-cycle proof](fresh_private_affine_cycle_exclusion_20261017.md)
now applies to these same literal roots without requiring A to exhaust
the private diagram. For completeness, each deletion of a v-root leaves
a long A path. The surviving attachment of z is at an endpoint with
label3 or4, or at a next-to-end node with label3. Hence the cyclic distance
between v1,v2 is one or two. Distance two forces both contacts to have
label3 and gives an affine four-cycle. At distance one, two label3
contacts give affine A2. Both are forbidden.

In the remaining case, one contact, say z.v1, has label4. The diagram
(A-v2)+z is the large finite B_m path with z its doubled terminal.
Adding y would extend that component at z. A connected finite extension
cannot do this: a second attachment to the existing path creates a cycle,
and a single attachment at z places the doubled bond internally in a
path too long to be F4. The complete finite classification permits an
extension of this B_m only at its opposite long endpoint. Nevertheless
(A-v2)+z+y is a principal subdiagram of an actual finite corner, a
contradiction.

Thus neither H-U nor H-V can contain an affine cycle of order at least18
under the stated hypotheses. The proof uses only finite-diagram
classification, the explicit small affine witnesses, and actual finite
corner incidence. It introduces no rank census, label cutoff, arithmetic
assumption, or global realization of a partial Gram matrix. The remaining
large affine-tree carriers require a separate argument.

The frontier reviewer independently checked the two-cycle intersection
argument, the single-contact witnesses, the shortest white-path step, and
the actual finite-corner scope of the terminal-extension contradiction.
No private-exhaustion assumption or assertion about three/four-black
cycles is used. This is a written mathematical/source review, not a fresh
graph enumeration, an independent arithmetic implementation, external
refereeing, or formal verification. Local report links were checked.
