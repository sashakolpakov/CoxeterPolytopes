# A coupled D4 lift of the compact Makarov P2 prism

The construction gives an eleven-facet Coxeter chamber in H9 with 14 finite and 16 spacelike projective extreme rays. Four polar truncations are compatible with every seed wall and with one another. Their coherent fifteen-facet extension has **49 finite vertices and 89 spacelike extreme rays**, no ideal rays, and infinite volume. An exact literal lattice check also closes further searches restricted to its integral reflection group. No new compact example or dimension is established.

The [producer](../scripts/fresh_h9_d4_prism_coupling_20261014.py) and [certificate](fresh_h9_d4_prism_coupling_20261014.json) pass generation and saved replay. The [independent review](fresh_h9_d4_prism_coupling_20261014_review.md) uses separate Fraction arithmetic over Q(√5) and incremental cone clipping; it reproduces every ray and active set, all facet and escape checks, and the literal integral reflection matrices.

The known compact source is Makarov's arithmetic prism P2, represented by the truncated chain [5,3,3,3,4]. Its diagram and arithmetic identification appear in [Emery–Kellerhals, Sections 1–2](https://homeweb.unifr.ch/kellerha/pub/IML-1112s-21.pdf). No novelty is claimed for the H5 source.

Put φ=(1+√5)/2. In six coordinates let J have diagonal (2,2,2,2,2,1), consecutive off-diagonal entries −φ,−1,−1,−1,−1, and all other entries zero. The source normals are e0,…,e5 and −φJ⁻¹e5. Only e5 has norm one; the other six have norm two. The complete source cone has ten finite extreme rays and no others. The sum of its normals is timelike, so its strictly positive dual slice is bounded and the census proves compactness.

Let D be the norm-two D4 Cartan matrix, with node 0 joined to nodes 1,2,3. The outer fundamental weight ω=D⁻¹e1 has norm one. Use the ambient form Q=J⊕D, retain six source walls, replace the norm-one wall e5 by

```text
n5 = e5 − ω,
```

and add the four D4 simple walls. All eleven normals now have norm two. The source π/4 contact becomes ambient π/3; its dotted end contact becomes π/5; and the lifted wall meets D4 node 1 at π/3. Every other new cross pair is right. The common D4 face is exactly the original compact H5 prism. This coupling differs from an orthogonal support extension.

The ambient form has signature (9,1). The certificate gives a strict timelike center and actual relative-interior points for all facets. The producer examines all 55 nine-wall supports, while the independent reviewer clips a simplicial cone using exact active-wall ranks. Both give 30 extreme rays: 14 finite and 16 spacelike. Their complete lists, not merely selected finite vertices, are retained.

Every spacelike ray is tested as a polar normal. The squared finite cosines belong to Q(√5), so the complete finite-label alphabet is 2,3,4,5,6,10. The cyclotomic degree condition leaves only these quadratic-field possibilities; labels 8 and 12 belong to different quadratic fields. Exactly four polars pass. Their six mutual pairs are all strictly dotted. The resulting fifteen-wall chamber has 81 right pairs, nine π/3 pairs, one π/4 pair, three π/5 pairs and eleven dotted pairs.

For this extension the producer exhausts all 5,005 nine-wall supports. Independent incremental clipping gives the identical 138 rays, including all incidences and norms. All fifteen walls are actual facets. The certificate also gives a literal rational vector E with positive square, negative pairing with the strict center, and strictly negative pairing with every wall. Thus an open patch at infinity survives: the chamber has infinite volume.

There is a stronger arithmetic stopping rule. Enlarge the D4 lattice by ω, using e0,e1,e2,ω as its basis. Together with the six source coordinate vectors this gives an integral Gram matrix over Z[φ] of determinant (1−√5)/2, a unit. It contains a norm-one vector, has signature (9,1), and its conjugate is positive definite. All fifteen physical reflection matrices are integral in this basis. The independent reviewer also verifies each matrix is involutive and preserves the form.

[Bugaenko 1984](https://www.mathnet.ru/php/getFT.phtml?jrnid=vmumm&option_lang=eng&paperid=3408&what=fullt), printed page 6, identifies the unique odd unimodular form with these signatures and states nonreflectivity in dimensions at least eight; pages 10–12 prove the required infinite stabilizers. Consequently no finite-volume completion can use only reflections preserving this literal lattice. This invocation uses the lattice classification and exact integrality checks, not merely an isometry of rational quadratic spaces. Nonintegral reflections, other lattices and larger fields remain outside the conclusion.

The [mixed-source coupled lift](fresh_h9_mixed_d4_coupling_20261014.md) changes the arithmetic data and therefore provides a separate test. Further integral root enumeration in the present host is stopped.

```sh
python scripts/fresh_h9_d4_prism_coupling_20261014.py --verify
python scripts/fresh_h9_d4_prism_coupling_20261014_review.py --verify
```
