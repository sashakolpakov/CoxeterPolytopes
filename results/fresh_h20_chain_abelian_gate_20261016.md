# A nonabelian rotation obstruction for the H20 chain

**Every coherent additional Coxeter wall retaining the frozen 133-wall H20 chamber must be dotted to at least five of its 26 chain walls.** In fact, every five consecutive chain walls must include a dotted relation with the added wall. This excludes the entire family whose chain contacts are finite or parallel, without any angle-label or coefficient-field cutoff.

This is an internally checked proof candidate pending separate independent review. It gives a necessary contact pattern, not a new wall or a finite-volume construction. Caps with sufficiently many dotted chain contacts, changed retained walls and other H20 constructions remain outside the exclusion.

The [standalone checker](../scripts/fresh_h20_chain_abelian_gate_20261016.py) and [certificate](fresh_h20_chain_abelian_gate_20261016.json) bind the [exact chain construction](fresh_h20_s61_multicap_20261014.md) and the separate [uniform-angle gate](fresh_h20_uniform_chain_cap_gate_20261016.md). The checker uses only Python's rational and integer arithmetic; there are no producer arithmetic, computer algebra or numerical factorization imports.

## The contact recurrence

All 26 chain roots have square 2. Their time coordinate is a common `a`, and their spatial coordinates rotate with

`C=a²/(a²+2)`, `S=2√(a²+1)/(a²+2)`.

The order 107,…,132 traverses the same circle in equal angular steps. For any real proposed normal, let `f_j` be its pairing with the `j`th chain root. The constant time part disappears in differences `δ_j=f_(j+1)−f_j`, giving

`δ_(j+2)+δ_j=2Cδ_(j+1)`.

Equivalently,

`f_(j+3)−f_(j+2)+f_(j+1)−f_j=2C(f_(j+2)−f_(j+1))`.

If five consecutive values lie in `Q^ab`, the two resulting recurrence equations imply either `C∈Q^ab` or both internal differences vanish. In the latter case the two outside differences vanish too. Two consecutive zero differences propagate through the recurrence, so all 26 pairings are equal. Four values alone would not suffice: a palindromic four-value block can have zero middle difference without being constant.

We now prove `C∉Q^ab` by two finite-field certificates.

## The actual endpoint polynomial

Let `w` be the exact saved `W_norm` from the chain's frame. It lies in the inherited field

`K0=Q(√2,√79,√11,√3)`,

`K1=K0(√D)`, `K=K1(√H)`,

where `D` and `H` are the literal radicands in the frozen cusp-repair and bridge certificates. The four prime squareclasses are independent, so `[K0:Q]=16`. The checker proves both additional radicands nonsquare:

| Radicand | Degree-one specialization of the preceding field | Nonsquare residue |
|---|---|---:|
|`D`|`p=313`, roots `(√2,√79,√11,√3)=(193,198,295,230)`|244|
|`H`|`p=359`, base roots `(340,325,332,196)`, `√D=167`|342|

The two residues have Euler criterion−1. Every denominator is a unit at the corresponding prime. All preceding quadratic roots are simple, so they lift by Hensel's lemma to embeddings in the corresponding p-adic field. A square in the preceding number field would have square residue there, since the radicand is a local unit. Thus `[K:Q]=64`. This also ensures that later choices of simple radical roots describe embeddings of the actual field, rather than spurious components of a reducible formal tower.

Write `U_j` for the Chebyshev polynomial of the second kind and define

`A(C)=U12(C)−U11(C)`, `B(C)=U12(C)+U11(C)`.

Twelve applications of the frozen spatial rotation give the endpoint coordinates exactly as

`b=−√((1+C)/(1−C)) A(C)`, `c=B(C)`.

The checker independently verifies this polynomial recurrence and the identity

`R+V=2`, where `R=(1+C)A²`, `V=(1−C)B²`.

The defining endpoint equation is `a√(w−5/6)+b/√2−c√w=0`. Substituting `a²=2C/(1−C)`, multiplying by `√(1−C)`, and squaring gives

`d w−e=√(2wRV)`,

where `d=2C−V` and `e=5C/3+R/2`. Here A and B are positive at the saved endpoint (its b is negative and c positive). A second squaring proves that the actual isolated parameter `C` satisfies

`P(C,w)=36[(d w−e)²−2wRV]=0`.

This integer polynomial has degree 50 in `C` and degree 2 in `w`. Its three complete coefficient lists are saved. Squaring may add other roots, but cannot remove the actual root; the next irreducibility certificate applies to this whole polynomial.

## Irreducibility and failure of normality

The checker gives two actual degree-one primes of `K`. The root choices are residues of `(√2,√79,√11,√3,√D,√H)`:

| Rational prime | Six roots | `w` residue | Leading coefficient of `P` | Factor degrees |
|---|---|---:|---:|---|
|2017|`(1031,206,1360,1759,1831,1772)`|1649|506|50, irreducible|
|359|`(340,34,332,163,340,222)`|15|230|13 and 37, both irreducible|

All radicand equations, nonzero derivatives and coefficient denominators are checked exactly. Successive simple Hensel lifts embed `K` into `Q_p`, giving residue degree 1 and ramification index 1. The polynomial's leading coefficient is a unit, so its monic normalization is integral at each prime. Its reduction is squarefree at both primes.

At 2017, the checker proves `X^(2017^50)=X mod P`, together with the required coprimality tests at exponents `2017^25` and `2017^10`. This is a complete irreducibility test for degree 50. At 359, Frobenius gcds isolate degrees 13 and 37; each factor passes the prime-degree irreducibility test, and their product equals the normalized reduction. No numerical factorization or factor-degree guess is accepted.

Irreducibility at 2017 proves that `P` is irreducible over `K`, hence `[K(C):K]=50`. At 359, squarefreeness ensures that the polynomial order is unramified and maximal locally, so the two factor degrees are the residue degrees 13 and 37 of primes in `K(C)`. They are unequal. A Galois extension has equal residue degrees above each base prime; consequently **`K(C)/K` is not Galois**. These are the standard local factorization and Galois-transitivity facts in [Sutherland, Number Theory I, Lecture 6, Dedekind–Kummer](https://math.mit.edu/classes/18.785/2025/LectureNotes6.pdf) and [Lecture 7, Corollaries 7.3–7.4](https://math.mit.edu/classes/18.785/2025/LectureNotes7.pdf). The unit and squarefreeness hypotheses needed here are explicitly certified.

If `C∈Q^ab`, then `Q(C)/Q` is finite Galois abelian: every subgroup of an abelian Galois group is normal. Its compositum `K(C)/K` would therefore be Galois. The preceding unequal splitting degrees contradict this. Thus `C∉Q^ab`.

## From field obstruction to Coxeter contacts

Normalize an added spacelike wall to square 2. A finite contact with a chain root is `−2cos(π/m)`, which lies in a cyclotomic, hence abelian, extension. A parallel contact is −2 and has the same property. Therefore five consecutive nondotted chain contacts would force all 26 contact values to be equal.

For a common finite value, the frozen uniform-angle gate supplies the contradiction. The parallel endpoint requires a small additional logical step, leaving that frozen theorem unchanged. Equal pairings force residual `t x`, and the common parallel value gives `t=2/a`. The same support budget lies strictly in

`2 < dᵀJd=2+4/a² < 141/70`.

The exact finite-contact gap from the uniform certificate still applies. Its only remaining support profiles are the single dotted 17/31 axes, and its physical blocker intervals already cover the **closed** parameter interval through `t=2/a`. Old 83 and 11 therefore still block those axes. The checker verifies these endpoint and interval dependencies explicitly.

Hence every five consecutive positions in the 26-wall chain contain a dotted contact. The five disjoint blocks of length 5 in the first 25 positions give at least five dotted relations overall. This statement concerns every real coherent added Coxeter wall retaining the 133 walls; it assumes no arithmeticity of the added wall or its coefficient field.

The geometric consequence is constructive guidance: an all-finite contact search cannot succeed on this chain, even with unbounded labels. A viable unequal-contact cap must instead contain multiple genuinely dotted relations. This result does not decide whether such a cap or a finite-volume completion exists.

## Reproduction and review status

```sh
python scripts/fresh_h20_chain_abelian_gate_20261016.py --verify
```

Generation and a fresh saved replay pass. All45 recursive manifests and204 source references agree, as do the syntax and local-link checks. Separate independent arithmetic/proof review is pending. Frozen predecessor files are unchanged.
