# Reviewed research proof candidate: dimension at most 166

**Candidate theorem.** Every finite-sided, finite-volume real hyperbolic Coxeter polytope has dimension **n≤166**, without an arithmeticity, compactness, simple-cusp or facet-count assumption.

The [common-height estimate](fresh_nikulin_coupled_charge_20261015.md) and [independent review](fresh_nikulin_coupled_charge_20261015_review.md) sharpen the frozen [169 candidate](fresh_nikulin_unit_charge_20261015_dimension169.md). The [case binder](../scripts/fresh_nikulin_coupled_charge_20261015_dimension166.py) and [certificate](fresh_nikulin_coupled_charge_20261015_dimension166.json) bind both branches and recursively verify the cited source hashes. This is an internally reviewed research proof candidate; external validation, novelty and current-best literature status remain unestablished.

Nikulin's narrow-facet theorem supplies a facet F whose framing hyperplanes have pairwise distance strictly below `arccosh 7`. If F contains an ideal vertex, its Euclidean link is a product of affine Coxeter simplices with total rank `N=n−1`. The reviewed endpoint and actual-ridge arguments supply two framing walls x,y with total retained-support cost at most four and one common height ratio `M≥1`. [Nikulin, official ICM report, printed pp.654–655 and 657–659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf).

For `M≤21/20`, the local charge improves to `7/6`, giving `−(x,y)≥M(N−40)/9`. For `M≥21/20`, the earlier `4/3` charge gives `−(x,y)≥M(N−43)/9`. Both use identical preset omissions, including the node-avoidance rule for F. At `N≥166`, the first lower bound is at least 14 and the second at least `287/20>14`. Each contradicts narrowness. If F has only finite vertices, it is compact, and the separately reviewed [compact-facet argument](fresh_nikulin_narrow_20261015_compact_facet.md) gives `n≤143`. The complete split yields 166.

The new estimate covers every affine type and every rank, with arbitrary contact magnitudes satisfying the Coxeter zero-or-at-least-one gap. Exact exceptional support enumeration is combined with symbolic classical-family inequalities. The published geometric and classification inputs remain mathematical dependencies, not formally certified software axioms.

For literature comparison, Prokhorov's published bound is 995 and Nikulin's ICM report announces the stronger bound below 300. These are distinct from the repository's candidate. The present improvement concerns unrestricted finite-volume polytopes; the published arithmetic ceiling 29 and noncompact arithmetic ceiling 21 are already much smaller. [Prokhorov, printed p.402](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&paperid=1484&what=fullteng), [Nikulin, printed p.659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf), [Belolipetsky's survey](https://arxiv.org/pdf/1506.03111).

```sh
python scripts/fresh_nikulin_coupled_charge_20261015_dimension166.py --verify
```
