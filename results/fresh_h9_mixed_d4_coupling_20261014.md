# A mixed compact H5 source with coupled D4 support

This exact H9 construction leaves the classified arithmetic host of the [P2 lift](fresh_h9_d4_prism_coupling_20261014.md). It still has infinite volume. Its eleven-wall seed has ten finite and twenty spacelike extreme rays. Exactly one extreme-ray polar passes all Coxeter pair tests; the twelve-wall extension has **twelve finite and thirty-eight spacelike extreme rays**. None of those thirty-eight final polars passes. Thus this iterated extreme-ray polar procedure is exhausted, rather than merely stopped at a search limit. Arbitrary new walls remain untested.

The [producer](../scripts/fresh_h9_mixed_d4_coupling_20261014.py) and [exact certificate](fresh_h9_mixed_d4_coupling_20261014.json) provide the forms, physical normals, all rays, every polar trial, actual facet points and a strict rational escape. The [previous mixed-source report](fresh_h9_compactification_20261013.md) supplies the compact H5 source and its literature context. The present change is a nonorthogonal D4 coupling, replacing the earlier orthogonal A4 support.

The [independent review](fresh_h9_mixed_d4_coupling_20261014_review.md) implements separate quartic Fraction arithmetic and enumerates all 21 source, 55 seed and 220 final defining supports. It reproduces every ray and incidence, all twenty plus thirty-eight polar trials, all 66 final pairs, twelve actual facets and the rational escape. It also checks both cyclic-field generators and the fixed indefinite block, with the quasi-arithmetic criterion explicitly retained as an external theorem. Primary and review generation and saved replays pass.

Let G be the norm-two five-node chain with successive products −φ,−1,−1,−1, where φ=(1+√5)/2. Put

```text
d  = (G⁻¹)[4,4] = (1−√5)/4,
g  = −G⁻¹e4,
D0 = 2−d = (7+√5)/4,
R  = (16+4√5)/11.
```

The source form is Q5=G⊕[D0], with normals (ei,0), i=0,…,4, and

```text
r5 = (g,1),         r5² = 2,
r6 = (2g,−√R),      r6² = 4.
```

These two ends come from different arithmetic Makarov prisms glued along the common orthogonal simplex facet. This is the known mixed-prism construction mechanism, not a new H5 example; see [Bogachev–Douba–Raimbault, Sections 3.3–4](https://arxiv.org/html/2309.07691v3). The code reconstructs and completely censuses this particular compact source, obtaining ten finite rays and no others.

Adjoin the norm-two D4 block with central node 0 and outer nodes 1,2,3. The outer weight ω=D⁻¹e1 has square one. Retain source normals 0,…,5 and replace r6 by

```text
n6 = r6/2 − ω.
```

Together with the four D4 simple normals, this gives eleven norm-two walls in signature (9,1). The source π/4 contact becomes π/3, the end pair remains strictly dotted, and the new support contact is π/3. Restriction to the common D4 face gives exactly the compact H5 source.

Exact incremental cone clipping lists all seed rays. A separate enumeration of all 55 nine-wall supports, using the same field implementation, reproduces the identical thirty projective rays. This is algorithmic cross-checking, not a claim of independent scalar arithmetic. The sum of the normals is timelike and strictly positive in the dual cone. It supplies a bounded projective slice, so the cone census is global.

All calculations take place in L=Q(√5,√R). This quartic extension is not normal: Norm(R)=16/11, and neither this number nor its quotient by five is a rational square. Its only nontrivial cyclotomic subfield is Q(√5). Therefore a finite squared cosine in L can have only labels 2,3,4,5,6,10; the polar rejection gate covers every possible label for these literal candidates.

Only one of the twenty seed polars passes, and it is strict at the common timelike center. After adding it, a complete new cone census gives fifty extreme rays, twelve finite and thirty-eight spacelike. Every final spacelike ray has at least one explicit forbidden finite pair as a polar normal. The final chamber has 53 right pairs, nine π/3 pairs, one π/5 pair and three dotted pairs. Explicit center projections establish all twelve actual facets. A saved rational E has E²>0 and strictly negative wall and center pairings, proving infinite volume.

The arithmetic change is substantive. The end Cartan two-cycle is

```text
(n5,n6)² = 2 + √5/2 + (−1+3√5)√R/8.
```

Together with the π/5 cycle it generates L. The embedding fixing √5 and negating √R fixes the principal Gram block G, whose signature is (4,1). The conjugate semidefiniteness condition in Vinberg's quasi-arithmeticity criterion therefore fails; see [Bogachev–Yorov, Theorem 2.2](https://nvbogachev.github.io/Quasi_arithmetic_Coxeter_prisms-2024.pdf). Conditional on finite volume, any completion retaining these seed facets would be non-quasi-arithmetic. This is an obstruction to arithmeticity, not a proof that such a completion exists. If a completion were quasi-arithmetic, its totally real Gram field would admit an extension of this embedding, still fixing the same indefinite principal block and acting nontrivially on the cyclic field.

The [finite-corner frontier](fresh_h9_mixed_d4_coupling_20261014_frontier.md) supplies the next search queue: all twelve finite vertices are simple, with 108 directed releases. Fifty-six directions form 28 finite edges; 52 remain open toward infinity. Every open direction has an exact positive-definite eight-wall support and a certified continuation to a relative-interior spacelike point. Forty-four release D4 walls; eight release old source walls at the two new polar vertices. This is a complete queue for the existing finite vertices, not an enumeration of additional mirrors or edges without a finite endpoint.

The next meaningful extension must solve general endpoint constraints on those supports using normals other than the exhausted polars, or change retained walls. Repeating this polar closure cannot help: its final candidate set is empty. No compact H9 polytope or new finite-volume dimension is claimed.

```sh
python scripts/fresh_h9_mixed_d4_coupling_20261014.py --verify
python scripts/fresh_h9_mixed_d4_coupling_20261014_review.py --verify
```
