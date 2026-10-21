# Independent review of the noncrystallographic face gate

**PASS.** Source and proof review found no issue in the [face gate](fresh_noncrystallographic_face_gate_20261014.md). Its saved producer replay passed. A [separate SymPy checker](../scripts/fresh_noncrystallographic_face_gate_20261014_review.py), importing no producer arithmetic, reproduced the weighted-seed projection identities and the actual unit-facet positivity conditions; its [certificate](fresh_noncrystallographic_face_gate_20261014_review.json) binds the final source and data.

The general face step is a direct application of [Allcock, Theorems 2.1 and 2.2](https://arxiv.org/pdf/0903.0138): the spherical pair I₂(m), for m outside `{2,3,4,6}`, has neither an Aₙ nor a D₅ component. Its actual codimension-two face is therefore Coxeter. The source permits induced angles different from the original ones, including the H₄ extension mentioned in the report.

The separate compactness argument is valid for a **finite-sided, finite-volume** ambient polytope. At an ideal vertex of this ridge, the two facet reflections would act on a bounded Euclidean cusp link. Their product's irrational linear trace contradicts preservation of the full translation lattice. Thus no ridge vertex is ideal; the projective face has only finite vertices and is compact. This proves the conditional compact H18 requirement, not nonexistence of all H20 completions.

For the weighted seed, the independent arithmetic checks the inverse of the central/golden Gram, all four projected-square types, the two projected unit-side pairings, and a positive-definite three-root Gram with nonnegative inverse. This last check supplies the strict-center projection proof that all twenty unit sides are actual facets of the original ridge. Heavy, other golden and opposite target sides restrict to timelike normals; they do not bound it.

The right-angled stopping rule also has the stated scope. New normals perpendicular to the defining pair keep their norms on restriction. Their right contacts stay right, while a parallel or dotted contact with an old unit side becomes divergent because that side's norm decreases. Any surviving induced facets therefore form a right-angled chamber. Compactness would contradict [Potyagailo–Vinberg, §2, Proposition, printed p.65](https://ems.press/content/serial-article-files/43032), which bounds compact right-angled hyperbolic dimension by four. The argument still works if some old unit facets become redundant on the completed ridge; deleting them cannot introduce nonright angles among the remaining facets.

This is a proof/source review plus independent finite arithmetic. It does not independently re-prove the cited general face or dimension theorems, and it does not exclude changed-angle or changed-ridge completions.

```sh
python scripts/fresh_noncrystallographic_face_gate_20261014_review.py --verify
```
