# Independent review of the uniform H4/H3 section obstruction

**PASS.** The [uniform proof](fresh_h20_h4_face_closure_20261014.md) has a sound geometric and arithmetic implication, under its stated requirement that the projected H4 or H3 is an actual face. The [separate checker](../scripts/fresh_h20_h4_face_closure_20261014_review.py) and [hash-bound result](fresh_h20_h4_face_closure_20261014_review.json) use SymPy directly on the integer source; they import no producer or repository arithmetic engine.

The independent arithmetic reconstructs the following:

- Signature (21,1), using exact real-root counts of the characteristic polynomial.
- All 210 root reflections preserve the actual lattice \(\mathbb Z^{22}\); every pairing involving a norm-four source root is even.
- The projected H4 norm and angle, all eighteen literal producer table rows, and a larger 25-row integer/parity table covering the same possible finite contacts.
- Exactly three forbidden finite squared cosines; every other possible adjacency involving a modified facet is orthogonal to an unmodified rational root or cannot meet at a finite point.
- The complete \(\mathbb Q(\sqrt5)\) Coxeter squared-cosine pool, without an angle cutoff. The elementary prime-power inequality \(\varphi(m)^2\geq m/2\) makes \(\varphi(m)\leq4\) a finite check through 32. Labels 8 and 12 give the distinct fields \(\mathbb Q(\sqrt2)\) and \(\mathbb Q(\sqrt3)\).
- The saved inventory of 6,492 literal edge supports, containing 1,660 A4 and 2,115 A3 components, hence 3,320 and 4,230 endpoint families.
- The H3-to-H4 attachment equation \(d-qk^2=2g^2\), whose only eligible integral solution is \((d,k,g)=(2,0,-1)\).

The geometric steps were reviewed separately. H4 has no connected spherical extension and belongs to no affine Coxeter diagram. An actual H4 face of a finite-volume Coxeter section is therefore compact, and its boundary facets are orthogonal to H4. The exact pair table then either makes this entire H16 face rational, or gives a modified facet whose compact H15 geometry is bounded entirely by original rational roots orthogonal to it. The latter are literal normals inside the face, not further irrational projections.

Compactness forces those normals to span the face's full Lorentzian space U. Their rational coordinates make U rational, while the independently checked source reflections preserve both U and \(\mathbb Z^{22}\). Thus \(U\cap\mathbb Z^{22}\) is an actual invariant full-rank lattice. The compact reflection chamber gives a cocompact lattice subgroup of its integral orthogonal group, necessarily of finite index. Meyer isotropy in rational rank at least five and the arithmetic compactness criterion contradict this. **Rational Gram entries alone would not justify this conclusion**; the invariant-lattice and finite-index steps are essential and are explicitly present in the proof.

The H3 extension is also sound. Its face is compact because H3 has no affine extension. Any nonorthogonal boundary facet supplies an actual H4 subface. The far-end label-3 equation forces the added original root to have norm two, original pairing −1 at that end, and zero pairing at the changed anchor. The old four roots consequently form A4 with exactly the already excluded one-contact H4 modification. If no such facet occurs, all boundary facets are orthogonal to H3 and the same table forces a compact rational H17 or H16 chamber, giving the same lattice contradiction.

This review does not independently replay the 152,368 numerical proposals, which are unnecessary to the uniform theorem. It is an independent finite-arithmetic replay plus a mathematical source/proof review; it does not formalize the general Coxeter classification or arithmetic lattice theorems in a proof assistant.

```sh
python scripts/fresh_h20_h4_face_closure_20261014_review.py --verify
```

The conclusion covers the stated actual one-contact H3/H4 faces, including all real parameters in the 7,550 recorded edge families. It does not exclude different coupled local diagrams or arbitrary sections of the source.
