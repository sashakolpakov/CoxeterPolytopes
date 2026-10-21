# Global obstruction for all 91 rank-three cusp-section records

No finite-volume H20 Coxeter section survives this family. All 91 records in the [complete rank-three local census](fresh_h20_rank3_omitted_sections_20261014.md) are now excluded over their entire admissible real shift intervals, including the endpoints where a local wall becomes redundant. The result concerns sections of the literal Borcherds H21 chamber through the three recorded rank-three affine factors, with exactly one local facet omitted. It does not exclude sections avoiding the source cusps or constructions that change the ambient chamber.

The [exact certificate](fresh_h20_rank3_global_sections_20261014.json) contains 75 continuous actual obtuse-ridge proofs and 16 continuous actual forbidden-angle proofs. Twenty-two endpoint cases need a different ridge: four are obtuse, and eighteen have forbidden angles. Every angle gap is strictly between zero and 1/4 or between 1/4 and 1/2 in squared cosine, so these are gaps between consecutive Coxeter values. This proof needs neither an arithmeticity assumption nor the separate crystallographic-angle theorem.

The 91 records consist of 44 quartic, fourteen quadratic, and 33 rational representations. They count omission/target/direction records, not projectively distinct planes. The local census supplies the complete finite list; the present checker imports none of its arithmetic or producer code.

## From a local direction to the whole family

For each record, let `D=sum m_i r_i` be the positive affine null relation, let `R` be the three retained roots, and write

`e0=sum_(i in R) k_i r_i`, with every `k_i>0` and `q=e0²>0`.

The source census gives an exact real-algebraic representation of every coefficient. Put `L=min(k_i/m_i)>0`. All admissible shifts are

`e(u)=e0-u L D`, with `0<u<=1`.

The normal square is the constant q because D is null and orthogonal to e0. The omitted inequality is redundant exactly on this interval. At u=0 the section misses the strict source interior, since e0 is a positive combination of retained walls. The endpoint u=1 is kept: the checker verifies

`e0-LD=sum_(i in R)(k_i-Lm_i)r_i-Lm_o r_o`.

Every retained coefficient is nonnegative. This identity is also used to handle coincident or otherwise redundant projected walls, rather than assuming all original inequalities remain actual facets.

## Actual ridges throughout each open interval

The numerical scout found proposed witnesses for all 91 families using 101 interior samples per family. Those samples are not proof premises. The final [recipes](fresh_h20_rank3_global_sections_20261014_recipes.json) contain only wall indices, an obtuse/angle-gap choice, and the extra endpoint witnesses.

For each family the exact checker constructs rational timelike points x0 and x1 by projecting the source's strict center onto two positive-definite ambient faces. It verifies all 210 original inequalities and the complete zero set at both points. After normalization,

`(D,x0)=(D,x1)=-1`, `(e0,x0)=0`, and `(e0,x1)=-L`.

Therefore `x(u)=(1-u)x0+u x1` is future timelike, lies in the section, and has exactly the common old zero set K for every `0<u<1`. The old Gram on K is positive definite. Its section projections are independent: e(u) pairs with x0 as `uL`, whereas all K roots are orthogonal to x0. The projected roots lie in the positive-definite tangent space at x(u), and hence define a local orthant. A small perturbation makes every K inequality except the selected pair strict. Thus the selected pair bounds a genuine section ridge. The argument does not mistake an arbitrary bad hyperplane pairing for an actual ridge.

For roots i and j, let `a_i(u)=(e(u),r_i)`. The checker uses the scaled projected quantities

`p(u)=q G_ij-a_i(u)a_j(u)`,

`n_i(u)=q G_ii-a_i(u)^2`, and `n_j(u)=q G_jj-a_j(u)^2`.

They are polynomials of degree at most two with exact algebraic coefficients. An obtuse certificate proves `p>0`. An angle-gap certificate proves both projected norms positive, `p<0`, and the two strict inequalities locating `p²/(n_i n_j)` in the stated gap. Polynomial signs hold on the entire open interval, not just at the scout samples.

## All endpoint cases

At u=1, the chosen ambient endpoint face has exactly K plus the omitted wall active. The nonnegative redundancy identity has no positive coefficient outside K, so the omitted inequality imposes no additional local restriction. Also `(e(1),x0)=L>0`, preserving projected independence of K. The ridge is therefore actual at the endpoint too. The checker explicitly verifies positive projected pair norms and determinant there whenever the endpoint section meets the source interior.

Sixty-four records are rejected at u=1 by their original witness, and five degenerate to an old facet hyperplane that misses the source interior. The remaining 22 receive another pair with its own two ambient face points and redundancy check. Four new endpoint pairs are obtuse; eighteen have a squared cosine in one of the same two forbidden gaps. There are no unresolved endpoint records. Multiple minima in `k_i/m_i` are allowed and tested exactly.

## Exact arithmetic and reproducibility

The [checker](../scripts/fresh_h20_rank3_global_sections_20261014.py) uses the Python standard library only. Each field is represented as rational polynomials modulo the record's degree-one, degree-two, or degree-four polynomial. A rational Sturm sequence verifies that the archived interval isolates exactly one simple real root. Sign decisions use rational Horner enclosures, refining that isolating interval by exact bisection when needed. No floating-point decision enters the certificate.

For parameter-polynomial positivity, the checker removes exact endpoint factors, subdivides into rational intervals when necessary, and verifies strictly positive Bernstein coefficients. Quotient-field identities, positivity of the old face Grams, all original wall signs, and endpoint redundancy are reconstructed from the frozen inputs. The source manifests of the local census and complete corner walk are checked before use.

Generation and the saved-result replay both pass on the frozen source. The archive contains 305 strict parameter-polynomial positivity certificates and no unresolved endpoints.

The [independent review](fresh_h20_rank3_global_sections_20261014_review.md) also passes generation and saved replay. Its separate SymPy algebraic-field implementation imports no producer arithmetic and reconstructs all 226 face-point occurrences, all 305 polynomial certificates, and every endpoint case. The local conic census remains a separately audited dependency.

```sh
python scripts/fresh_h20_rank3_global_sections_20261014.py --verify
```

The [scout](../scripts/fresh_h20_rank3_global_sections_20261014_scout.py) and its [numerical output](fresh_h20_rank3_global_sections_20261014_scout.json) are retained as proposal provenance only. The exact certificate does not import them or depend on numerical field embeddings.

Together with the earlier rank-at-least-four obstruction, this closes the single-factor, one-omitted branch at the saved cusps. The broader consequence for all cusp-containing sections additionally depends on the separate [cusp-coverage argument](fresh_h20_cusp_coverage_20261014.md); it is not a general H20 nonexistence theorem.
