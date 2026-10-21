# A coupled B2-supported H7 face, and its exact arithmetic stopping rule

The construction gives a 13-facet Coxeter chamber in H9 with **56 finite vertices and 48 spacelike projective vertices**. Six of the 48 polar truncations pass every exact Coxeter gate and are mutually dotted. Adding all six gives a coherent **19-facet chamber with 94 finite vertices and 190 spacelike projective vertices**. Both chambers have certified strict spacelike escapes and infinite volume.

More decisively, all nineteen reflections preserve an explicitly identified odd unimodular admissible lattice over Z[φ]. Bugaenko's original theorem makes this H9 lattice nonreflective. **No finite-volume completion using only reflections preserving this particular lattice is possible.** This prevents further integral-root searches in the same host. It does not exclude other lattices in the rational space, nonintegral reflections, field extensions, or general compact H9 constructions.

The [producer](../scripts/fresh_dimension_new_mechanisms_20261014.py) and [exact certificate](fresh_dimension_new_mechanisms_20261014.json) contain the physical normals, both complete cone censuses, strict centers and escapes, actual facet points, all polar trials, and every reflection matrix in the unimodular basis. Generation and saved replay pass. Arithmetic uses the separately frozen Q(√5) scalar implementation named and hashed in the manifest; no implementation-independence claim is made here.

The starting compact H7 source is the unique eleven-facet example in Felikson–Tumarkin's classification, Theorem 3 and Figure 7. The diagram was read directly from their [author-maintained diagram](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/poly/7-11.pdf), saved [locally](../references/fresh_dimension_new_mechanisms_20261014_h7_diagram.pdf). Their [original paper](https://arxiv.org/pdf/math/0510238) establishes its classification. The exact geometry below is independently reconstructed from that diagram, rather than inferred from a diagram image alone.

Number the diagram nodes from left to right along the upper chain as 0,…,6; node 7 hangs below its central node 3; node 8 is the lower central node, and nodes 9,10 the lower left and right nodes. The finite nonright edges are

```
5: 0–1, 5–6
3: 1–2, 2–3, 3–4, 4–5, 3–7, 1–9, 5–10
4: 8–9, 8–10.
```

The dotted edges are 7–8 and 9–10. In norm-two normalization their hyperbolic coshes are φ/√2 and φ, respectively. These values follow by setting the three-by-three Schur complement of the first eight roots equal to zero. For computation we divide source root 8 by √2, leaving its norm 1 and all coordinates in Q(√5).

The first eight roots provide the Lorentzian source form B of signature (7,1). With J=B⁻¹, the remaining source normals are

```
v8  = −φ J e7,     v8²=1,
v9  = −J e1,       v9²=2,
v10 = −J e5,       v10²=2.
```

All 330 seven-wall active sets are checked exactly. They give 50 distinct finite source vertices and no ideal or spacelike vertices. The sum of the source normals is timelike and a strictly positive dual-cone functional. Thus its projective slice is bounded; the exhaustive ray census proves compactness. A sum of the normalized source vertices gives a strict timelike center.

The new mechanism uses a B2-supported face. In norm-two normalization the B2 inverse is

```
[ 1       √2/2 ]
[ √2/2    1    ].
```

Its entries are positive. A norm-two acute Coxeter normal inducing a facet of the B2 face must have positive projected square. Each nonright contact costs at least 1 in the transverse projection. Two nonright contacts cost at least 2+√2; a single contact of label at least 4 costs at least 2. Consequently the complete all-label list is: orthogonal to B2, or a single π/3 contact to one of its two nodes. This is also the local spherical-extension case of [Allcock's Theorems 2.1–2.2](https://arxiv.org/pdf/0903.0138), printed pages 743–744.

Call the two nonorthogonal possibilities the two special colors. If r is a norm-two source normal, its special lift is r/√2−ωj, where ωj²=1. The exact finite-pair rules are:

| Source facets | Permitted finite source contact | Ambient contact |
|---|---|---|
| one special, one ordinary | right or π/4 | right or π/3 |
| two specials, same color | none | — |
| two specials, opposite colors | π/4 only | right |

These are all-label statements: any remaining negative finite product falls strictly between zero and the smallest nonzero Coxeter product, or between consecutive allowable products. Every source node except 8 has a π/3 or π/5 neighbor, so it must be ordinary. Node 8 alone is eligible. Its dotted contact to node 7 becomes π/5, so it passes the full pair gate rather than only the finite-source-contact test.

For rationalized B2 coordinates use s0²=2, s1²=1, s0·s1=−1. Then d0=s0+s1 has square 1. Retain the ten ordinary source facets, replace v8 by v8−d0, and add s0,s1. The resulting 13-wall diagram has 63 right pairs, ten π/3, one π/4, three π/5 and one dotted pair. A common strict timelike center and all acute pairings certify every facet's actuality by orthogonal projection of that center.

There is a small complete global cone test. Every transverse vector in the B2 chamber is uniquely

```
u = −t0(s0+s1) −t1(s0+2s1),
t0,t1 ≥ 0,     t0+t1 ≤ −v8·x.
```

Projection x lies in the compact H7 source cone. The simplex fiber decomposes over source extreme rays. Thus the complete ambient ray set consists of the three fiber vertices above each source ray, with the three collapsing when v8·x=0. Every saved ray has nine independent active walls. This gives exactly 104 rays: 56 finite and 48 spacelike.

Each spacelike ray is next tested as a polar normal. The all-label angle gate is exhaustive: its squared cosines lie in Q(√5), so the cyclotomic degree condition permits only finite labels 2,3,4,5,6,10. Labels 8,12 have the wrong quadratic fields. Exactly six rays pass, and every pair of accepted polars is strictly dotted. Their coherent extension is censused by exact incremental cone clipping. Two old rays are adjacent exactly when the common active-wall span has rank eight; crossing edges give all new rays. The six stage sizes are 124,156,180,200,252,284. The final census is 94 finite and 190 spacelike, with no ideal rays. The certificate gives a strict spacelike perturbation and actual facet points for all nineteen walls.

Doubling the H7 source does not provide another special-color option here. Its only doubling wall is node 8. The exact thirteen-facet double has a literal odd π/3 or π/5 ridge at every facet. It has neither a further doubling wall nor any facet eligible for a nonorthogonal B2 lift. This statement concerns this source and its single legal first double, not arbitrary gluing or other H7 sources.

The arithmetic stopping rule uses the actual lattice, not only the rational quadratic space. In the ten-dimensional coordinate basis, replace e0 by the lifted normal of facet 8. The basis determinant is φ/2. The resulting Gram matrix has entries in Z[φ], determinant −φ, and a basis vector of norm 1. Its conjugate is positive definite. Hence it is an odd unimodular admissible lattice. Every one of the nineteen physical reflection matrices has entries in Z[φ] in this basis, as checked entry by entry.

[Bugaenko 1984](https://www.mathnet.ru/php/getFT.phtml?jrnid=vmumm&option_lang=eng&paperid=3408&what=fullt), printed page 6, identifies the odd unimodular admissible lattice with the form −φx0²+Σxi² and states that its reflection chamber has infinite volume for dimension at least 8; pages 10–12 establish the infinite stabilizers. Applying that lattice-specific theorem excludes any cofinite reflection subgroup contained in this host's full integral reflection group. It does not transfer nonreflectivity to every commensurable lattice.

This experiment changes the transverse contacts and produces actual new finite vertices, unlike an orthogonal product or an isolated dotted-cap layer. It still yields no new compact or finite-volume dimension. The exact stopping condition says what a subsequent experiment must change: use a reflection that fails the saved lattice-integrality test, a source with different arithmetic data, or a field-extension construction. Continuing integral polar/root enumeration in this host is closed.
