# Root review: mixed K8/K7 maximum-height corners

19 October2026. **PASS, full written-source review.** Source:
`right_angled_mixed_k8_maximum_height_carrier_defect_20261019.md`,
SHA256 `0ee6e122b76c9bf8c58d6001fbf2ca36488ba7b9fd0d6724e9bbd0efef3b2aa7`.
No enumeration, numerical optimization or Lean check was used here.

The new finite-edge length gap follows from two distinct actual wall
normals with product -sinh(L). The right-angled normal gap gives
sinh(L)>=1 at any finite edge beside a bad triangle. This applies to
every finite edge there, with no carrier-size hypothesis.

For a future timelike height direction, the normalized spatial vector
has squared norm strictly below one. Comparing actual opposite vertices
gives the stated triple and backed-pair inequalities, with the correct
signs. Three disjoint triples cost at least one; a backed pair and a
disjoint K5 cost at least1/2+5/9. A K8 costs at least8/9, so the
finite-edge gap excludes every remaining finite direction there.

Any triangle outside the K8 leaves six core vertices for two additional
disjoint triangles. This proves that exactly the56 core triangles
occur. All outside directions are ideal and each has at most one
core neighbor. Every triple with one outside direction and two core
directions is therefore a one-edge triple or an unbacked path. The
reviewed local nu<=11 theorem excludes both at that all-ideal-edge
corner. The112/84 distinct actual three-face counts follow in
dimensions12/11 respectively.

The K7 counting argument also checks: triangles with at most one core
vertex produce a forbidden packing of three disjoint triangles. If
one outside vertex has at least four core neighbors, any other with
two gives the same forbidden packing. Otherwise all have at most
three. Both cases give at most35+15=50 triangles for n<=12.

These are conditional facts about a height-maximizing finite vertex.
No K8 is asserted to exist there, and no bound on a global ordinary
triangle average or on E12/f is obtained. Clique number at most six
remains a local obligation. The unrestricted dimension ceiling12,
the manuscript source and Lean status are unchanged.
