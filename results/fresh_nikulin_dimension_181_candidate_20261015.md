# Reviewed research proof candidate: dimension at most 181

**Candidate theorem.** Every finite-sided finite-volume real hyperbolic Coxeter polytope in H^n satisfies **n≤181**, without arithmeticity, compactness, simple ideal vertices, or a facet-count restriction.

The [avoid-node estimate](fresh_nikulin_avoid_node_20261015.md) and its [independent arithmetic and proof review](fresh_nikulin_avoid_node_20261015_review.md) sharpen the frozen [189 proof](fresh_nikulin_dimension_189_candidate_20261015.md). The [dependency binder](../scripts/fresh_nikulin_dimension_181_candidate_20261015.py) and [certificate](fresh_nikulin_dimension_181_candidate_20261015.json) bind the complete chain. This is an internally reviewed research proof candidate. External validation, priority and current-best status are unestablished.

Nikulin's narrow-facet theorem gives a facet F whose framing hyperplanes have pairwise distance strictly below arccosh 7. The proof uses actual ridge neighbors of F. With outward normals of square two, two such walls cannot have negative pairing at least 14. [Official ICM report, printed p.655](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf).

If F has an ideal vertex, its horospherical link is a product of irreducible affine Coxeter simplices, with ranks summing to N=n−1. Preset omissions select two outgoing edges in F. At their other endpoints, finite volume supplies terminating walls x,y. At an ideal endpoint, one chooses a missing node in the appropriate affine component, or an orthogonal component if F belongs to an A1 factor. The frozen product proof verifies that these walls meet F along actual ridges.

For affine factor a, let C_a be the sum of the two retained-support Schur costs, and σ_a,τ_a its positive null heights. Their ratio q is common to all factors. Put M=max(q,1/q). The endpoint Gram matrices give ΣC_a≤4. The unchanged non-narrow estimate and new avoid-node estimate give, respectively,

    σ_a τ_a E_a / M + 2C_a ≥ 2r_a/9,
    σ_a τ_a E_a / M + 2C_a ≥ 2(r_a−1)/9.

Only the single factor containing F's normal uses the second inequality. It holds for every forbidden node, with choices made before the contact rows are known. The proof covers all classical ranks analytically and checks every exceptional omission exactly, including cheap retained-contact cases; it imposes no label cutoff.

Summing and using the Lorentz identity yields

    Σ_a σ_a τ_a E_a ≥ M(2N/9−8−2/9),
    −(x,y) = (1/2)Σ_a σ_a τ_a E_a − q − 1/q
           ≥ M(N−55)/9.

At n≥182, N≥181, so the pairing is at least 14, contradicting narrowness. This branch gives n≤181.

If F has no ideal vertex, it is compact. The separately reviewed [compact-facet argument](fresh_nikulin_narrow_20261015_compact_facet.md) gives n≤143. It allows affine subdiagrams in ambient stars and does not assume F is itself Coxeter. The two alternatives exhaust F and give the stated candidate.

The new local loss of 2/9 is sharp within the uniform avoid-node lemma: in A1, both omissions must be the only node other than the forbidden node, and equal point-contact rows have zero cost and energy. This does not prove that 181 is optimal for the global method or for hyperbolic Coxeter polytopes.

For historical comparison, Prokhorov's published Theorem C gives n≤995. Nikulin already announced the unrestricted n<300 in his ICM report, whose definition of the class assumes finite volume but not arithmeticity. Thus 181, if validated, improves both numerical bounds; a bounded literature check cannot establish novelty or rule out a later stronger theorem. [Prokhorov, printed p.402](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&paperid=1484&what=fullteng), [Nikulin, printed pp.654 and 657–659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf). Khovanskii's 2006 author account still states 995 as Theorem 0.1. [Author text](https://www.math.toronto.edu/askold/2006%20Coxeter%20Legacy.pdf).

The new primary report's pending-review sentence records its earlier frozen state. The completed review and this final binder supersede that status without rewriting the frozen report. Computational checks certify the exact arithmetic and source closure; the geometry and published classifications remain mathematical proof dependencies, rather than a formal proof-assistant development.

```sh
python scripts/fresh_nikulin_dimension_181_candidate_20261015.py --verify
```
