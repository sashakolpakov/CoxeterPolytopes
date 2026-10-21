# Independent final scope and quantifier check for the arithmetic bound-21 candidate

The [standard-library checker](../scripts/fresh_arithmetic_dimension21_scope_check_20261016.py) and [saved result](fresh_arithmetic_dimension21_scope_check_20261016.json) pass for the final [Q(sqrt5) coverage binder](fresh_sqrt5_h3_hub_coverage_20261016_bound21.json) and [all-field composition](fresh_arithmetic_dimension21_composition_20261016_bound21.json). Generation and fresh saved replay both passed. This is a separate check of the final composition, its quantifiers and rational bookkeeping; it does not replace the arithmetic and graph proofs in the imported certificates.

All **2,661 obligation identities and node counts** are independently matched to the originating catalogues: 112 original sealed rows, 147 singular first extensions, 51 marked singular frontier rows, 17 closure rows, 264 dotted stopping rows and 2,070 finite-queue stopping rows. Exactly **1,634** fit dimension 22; the other **1,027** have corner rank greater than 22 and carry no payment cases. Every feasible row has every common external count triple `(a,b,c)` satisfying `a+2b+3c <= residual_rank`, except when the leftover rank is one. These are actual A1/I2(5)/H3 counts, without capping in the quantifiers. Residual rank zero or at least two is represented by non-hub components, for example A2 and A3.

Across all **91,726** saved triples, there are no missing or duplicate triples. Every row has four nonnegative corner weights; their sum equals the saved payment and is at least one. The minimum saved quadrilateral payment is **560003/560000**. All 410 exact vertex-count contexts are present, and their saved maxima agree with the independently reviewed H3 certificate. The largest is **21977/500 = 43.954**.

The h-vector threshold is recomputed from `((22-i choose 2)+(i choose 2))/2`. Its minimum is **55**; the two extreme coefficients are strictly larger. With `h_0=h_22=1`, the averaged face count is strictly above 55, giving a strict vertex-weight lower threshold `5*55 - (22 choose 2) = 44`. The saved upper bound has margin **23/500** below that threshold. The 76 saved triangle baseline totals are all at least two.

The prior complete Q(sqrt5) binder covers dimensions 23 through 29; the new binder covers dimension 22 as well. The all-field composition has exactly **24 branches**: dimensions 22 through 29, each over Q(sqrt2), Q(sqrt5), and Q(cos(2pi/7)), with respective certified upper bounds 20, 21, and 20. The three-field restriction is never used below dimension 22. The separately imported published compact ceiling handles dimensions above 29, and the published noncompact arithmetic ceiling is 21, attained in dimension 21. Thus the stated scope of the composition is all finite-volume arithmetic reflection groups, with upper bound 21 within this research proof candidate.

All recursive input hashes passed: **94 manifests, 1,133 source references, and 193 distinct hashed dependencies**. The checker uses only Python's standard library and exact `Fraction` arithmetic, with no producer, optimizer, or repository proof-module imports.

```sh
python3 scripts/fresh_arithmetic_dimension21_scope_check_20261016.py --verify
```

This check does not reprove Gram-matrix arithmetic, queue exhaustiveness, literal graph weight evaluation, finite-forest dynamic programs, the applicability of the compact h-vector argument, or the published field and noncompact bounds. Those remain explicit bound inputs. External mathematical validation, formal proof-assistant verification, and priority remain open. Source and JSON are frozen.
