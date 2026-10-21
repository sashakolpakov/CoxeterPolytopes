# Reviewed research proof candidate: dimension at most 169

**Candidate theorem.** Every finite-sided, finite-volume real hyperbolic Coxeter polytope has dimension **n≤169**, without an arithmeticity, compactness, simple-cusp or facet-count assumption.

The [four-thirds estimate](fresh_nikulin_unit_charge_20261015.md) and its [independent review](fresh_nikulin_unit_charge_20261015_review.md) sharpen the frozen [181 proof candidate](fresh_nikulin_dimension_181_candidate_20261015.md). The [case binder](../scripts/fresh_nikulin_unit_charge_20261015_dimension169.py) and [certificate](fresh_nikulin_unit_charge_20261015_dimension169.json) bind the complete dependency chain. This is an internally reviewed research proof candidate; external validation, novelty and current-best literature status remain unestablished.

Nikulin's narrow-facet theorem supplies a facet F whose framing hyperplanes have pairwise distance strictly below `arccosh 7`. If F contains an ideal vertex, its Euclidean link is a product of irreducible affine Coxeter simplices of total rank `N=n−1`. Preset omissions select two edges in F, and actual terminating walls x,y can be chosen to meet F in ridges. Their two total Schur costs sum to at most four. [Nikulin, official ICM report, printed pp.654–655 and 657–659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf).

For each affine factor, the new estimate is

    σ_a τ_a E_a / M + (4/3)C_a ≥ 2r_a/9,

with target `2(r_a−1)/9` in the single factor containing F's normal. The height ratio is common across factors and `M=max(q,1/q)≥1`. The proof covers every affine type and all ranks, with no label cutoff. It combines a linear energy-loss identity with exact exceptional-case analysis; the ordinary E8 omission changes to nodes 1 and 6.

Summation and the Lorentz identity yield

    −(x,y) = (1/2)Σ σ_a τ_a E_a − q − 1/q
           ≥ M(N−43)/9.

At `n≥170`, this is at least 14, contradicting narrowness for square-two normals. If F has only finite vertices, it is compact; the separately reviewed [compact-facet argument](fresh_nikulin_narrow_20261015_compact_facet.md) gives `n≤143` without assuming that F itself is Coxeter. The two cases establish the candidate 169 bound.

For comparison, Prokhorov's published Theorem C gives 995, and Nikulin's ICM report already announces the stronger unrestricted bound below 300. The repository's candidate is separate from those literature statements and is not presented as an externally validated theorem. [Prokhorov, printed p.402](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&paperid=1484&what=fullteng), [Nikulin, printed p.659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf).

The compact and arithmetic dimension ceilings remain 29; the noncompact arithmetic range remains `n≤19` or `n=21`. The 169 candidate addresses the general finite-volume class and does not improve those stronger restricted bounds.

The primary report's pending-review sentence records its earlier frozen state. This binder records the completed independent review without rewriting that artifact. Arithmetic and source integrity checks do not replace the stated mathematical and published dependencies.

```sh
python scripts/fresh_nikulin_unit_charge_20261015_dimension169.py --verify
```
