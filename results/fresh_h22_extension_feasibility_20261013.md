# Retaining the Borcherds H21 chamber in H22: global constraints

The retained-210-wall construction is **not ruled out** by the primary sources checked here. There are two useful reductions: it suffices to search for one coherent covering family and reflect it, and the source has exactly **18 cusp types** under its chamber symmetry group. A separate arithmetic constraint gives a reason to try the existing cap recipes at angle π/6: a finite-volume completion retaining a π/7 pair would have to be **non-pseudo-arithmetic**, not merely nonarithmetic. The latter is a conditional restriction, not a nonexistence theorem.

The independent computation [script](../scripts/fresh_h22_extension_feasibility_20261013.py) and [result](fresh_h22_extension_feasibility_20261013.json) match all 18 published cusp types to the exact cusp representatives reached by the separate corner walk. The [frozen input](../data/enumerations/fresh_h22_extension_feasibility_20261013_cusp_input.json) is bound to that walk's final completed certificate. Its [separate report](fresh_h21_source_corner_walk_20261013.md) proves closure of the 395-corner-orbit queue; the cusp-type check here is an additional independent classification, not a duplicate proof of edge exhaustion.

## One side suffices

Use the Lorentz form `Q ⊕ <1>`, old coordinates x and new coordinate z. Let P be the finite-volume H21 chamber. Retaining its 210 vertical facet inequalities forces an acute-compatible cap r=(u,a) to satisfy `(u,w_i)≤0` for every old wall. Thus u lies in the closed old chamber cone, and is future causal or zero. The causal-cone premise and the exact whole-P covering criterion are proved in the earlier [covering report](fresh_h22_incidence_mixed_height_20261013_covering.md), independently [reviewed here](fresh_h22_incidence_covering_20261013_review.md).

Take two upper caps r=(u,a), s=(v,b), with a,b>0, squared norms q,r0>0, and put U=−u²≥0, V=−v²≥0. Reflect s to s−=(v,−b). Lorentz Cauchy and an elementary square identity give

```
(r,s−) = (u,v)−ab
       ≤ −sqrt(UV)−sqrt((q+U)(r0+V))
       ≤ −sqrt(q*r0),

(q+U)(r0+V)−(sqrt(q*r0)+sqrt(UV))²
       = (sqrt(q*V)−sqrt(r0*U))² ≥ 0.
```

So all cross-pairs between the original upper family and its reflected copy are dotted or parallel. The original same-side pairs are unchanged. The point `(h,0)` is strict for every cap with u nonzero, because h is future timelike. Each positive-norm inequality is an actual facet: projecting a strict timelike center orthogonally to one normal stays timelike and makes every other nonpositive-pairing inequality strict.

If a finite-volume extension exists, its upper and lower families each cover all of P in the previously proved sense. Reflecting either non-pure side therefore gives a finite-volume symmetric Coxeter extension. The sole exceptional cap shape is `(0,a)`: an upper pure-vertical cap has positive pairing with every other upper cap, so an acute family contains no other upper mirror. It bounds z≤0. A finite-volume full-dimensional extension must then have a non-pure lower family, which can be used instead. Pure-vertical caps on both sides collapse the interior. Thus this exception does not invalidate the reduction.

The upper half of the symmetric chamber has an additional facet z=0 isometric to P. It meets its 210 neighbors orthogonally; every cap is disjoint from or asymptotic to it. It is a doubling wall. This identifies the construction as an **orthogonal-facet extension** problem. It does not create an extension automatically. Allcock's doubling method starts with a chamber and a doubling wall and preserves dimension; it does not supply a chamber one dimension higher. His discussion also records that the original D4/H21 face is not redoublable. [Allcock, *Infinitely many hyperbolic Coxeter groups through dimension 19*, definition of doubling walls and Example 2.5](https://web.ma.utexas.edu/users/allcock/research/doubling.pdf).

## What arithmeticity does and does not force

Borcherds identifies the old chamber as the full reflection chamber of the even sublattice of I21,1, equivalently `U ⊕ D20`, with 42 norm-two and 168 norm-four walls and symmetry group of order 241,920. The integral-lattice dimension bound cited there concerns integral reflective lattices; it does not rule out general real Coxeter extensions. [Borcherds, *Automorphism groups of Lorentzian lattices*, §8 Example 5](https://math.berkeley.edu/~reb/papers/al/al.pdf).

Nor does an arithmetic cofinite hypersurface force its ambient lattice to be arithmetic or quasi-arithmetic. Hybrid constructions provide counterexamples to that reasoning. The modern subspace results prove inheritance from a (quasi-)arithmetic ambient orbifold to its suborbifolds, not this converse; their §6.2 discusses hybrid codimension-one examples. General geodesic embedding is less restrictive still: arithmetic manifolds of simplest type, or appropriate finite covers, embed one dimension higher. These theorems do not preserve a Coxeter chamber with the prescribed 210 facets. [Belolipetsky–Bogachev–Kolpakov–Slavich, *Subspace stabilisers in hyperbolic lattices*, Theorem 1.8 and §6.2](https://arxiv.org/pdf/2105.06897), [Kolpakov–Reid–Slavich, *Embedding arithmetic hyperbolic manifolds*](https://arxiv.org/abs/1703.10561).

A stronger **conditional** restriction follows from pseudo-arithmeticity. By the published definition its ambient group is a base change `G_K` of an admissible group over a totally real field k, with K/k **totally real** multiquadratic. Emery–Mila prove that finite gluings of arithmetic pieces have this property and ask whether all higher-dimensional hyperbolic lattices do. Their reflection-group criterion identifies K with the field generated by Gram cyclic products. [Emery–Mila, *Hyperbolic manifolds and pseudo-arithmeticity*, Theorem 1.5, Definitions 1.6–1.7, Question 1.8 and §4, published version](https://doi.org/10.1090/btran/48), [author preprint](https://arxiv.org/pdf/1810.12837). The total-reality clause matters in the argument below; it is explicit in the published definition.

Here the old reflection subgroup already contains a nontrivial unipotent, from an old Euclidean cusp stabilizer. Suppose a completion were pseudo-arithmetic. In a faithful representation of G, that unipotent remains nontrivial and unipotent under every embedding of K. If k had another real embedding, admissibility would make the corresponding real group compact. Every such embedding extends to a real embedding of K, and a compact real linear algebraic group has no nontrivial unipotents. This contradiction proves **k=Q**. Consequently K/Q must be multiquadratic. This is an argument for the retained source subgroup, and also for nonuniform pseudo-arithmetic lattices generally.

For a Coxeter pair of label m, the two-cycle `4 cos²(π/m)` belongs to K. Its real cyclotomic field therefore must have Galois group of exponent at most two. This forces

```
m ∈ {2,3,4,5,6,8,10,12,24}.
```

For completeness, the exponent condition is `a² ≡ ±1 (mod m)` for every unit a modulo m. Reduction modulo odd prime powers permits only primes 3 and 5, each to exponent one. Reduction modulo 16 excludes any factor 16. Thus m divides 120; checking its divisors yields the displayed list. The script performs that final finite check. This list applies **only under the pseudo-arithmetic hypothesis**.

In particular, `4 cos²(π/7)` has irreducible polynomial `x³−5x²+6x−1`, and cannot lie in a multiquadratic field. A finite-volume completion retaining one of the certified π/7 pairs would be non-pseudo-arithmetic, and could not arise from a finite gluing of arithmetic pieces. The 2026 version of *Subspace stabilisers*, §2.4, still reports all known examples in dimensions greater than three as pseudo-arithmetic. This is evidence that the π/7 branch asks for a substantially stronger breakthrough; it is not permission to reject it as impossible. The one-side reduction itself is independent of arithmeticity.

At π/6, the original two-direction cap recipe `p*u + sqrt(A+C*T²)*v`, with rational A,C and `T=2cos(π/m)`, instead has T²=3. Each such cap therefore belongs to a compositum of square roots of rational numbers. This removes this particular field obstruction. Every angle, sign, compatibility and volume gate still requires a new check; the report asserts no successful π/6 completion.

## Eighteen cusp types, with a complete target from the literature

Elkies–Schütt classify the definite rank-20 complements obtained from `U ⊕ D20` by embeddings of D4 into Niemeier lattices. Their Table 6.1 has 18 distinct norm-two root types. Their proof includes the cases where norm-two roots do not have full rank. [Elkies–Schütt, *Genus 1 fibrations on the supersingular K3 surface in characteristic 2 with Artin invariant 1*, §§5–6, Table 6.1 and Theorem 1](https://www.intlpress.com/site/pub/files/_fulltext/journals/ajm/2015/0019/0003/AJM-2015-0019-0003-a007.pdf).

The passage to our cusp orbits is as follows. The three nonzero discriminant classes of D20 are represented by e1, the half-all-ones vector, and that vector minus e1. Their norms are 1,5,5 modulo 2; the discriminant form is anisotropic. A primitive isotropic vector e in L consequently has divisibility one. Choose `(e,f)=1`, and replace f by `f−(f²/2)e`; evenness makes this integral and isotropic. Hence `L=U⊕M`, where M is positive definite, even, rank 20, determinant 4. Conversely such a splitting gives a primitive cusp ray. The isometry class of M determines the O(L)-orbit of e. Gluing M to D4 along their discriminant forms produces a Niemeier lattice; the published embedding classification therefore covers every case. Its distinct root types also distinguish the 18 complements.

Finally, the full lattice reflection group has P as a chamber and `O⁺(L)=W(L)⋊Aut(P)`. Chamber closures are fundamental domains for W(L), so O⁺(L)-orbits of future primitive null rays correspond to Aut(P)-orbits of P's cusps. Thus the classification yields exactly 18 cusp orbits here, not a list of sample rays. It does not provide the finite-corner census, orbit sizes or explicit representatives used by the computation.

| Table row | Norm-two root type | Snapshot cusp ID | Orbit size |
|---:|---|---:|---:|
| 1 | A5⁴ | 5 | 280 |
| 2 | D4⁵ | 58 | 42 |
| 3 | D5 + A7² | 0 | 3,780 |
| 4 | A1² + A9² | 2 | 6,048 |
| 5 | A1² + D6³ | 1 | 5,040 |
| 6 | D7 + A11 | 390 | 10,080 |
| 7 | A3 + E6 + A11 | 4 | 20,160 |
| 8 | E6³ | 7 | 2,240 |
| 9 | D4 + D8² | 6 | 30,240 |
| 10 | D5 + A15 | 8 | 30,240 |
| 11 | A1³ + A17 | 388 | 20,160 |
| 12 | A1³ + E7 + D10 | 10 | 60,480 |
| 13 | D6 + E7² | 3 | 60,480 |
| 14 | D8 + D12 | 11 | 120,960 |
| 15 | D4 + D16 | 13 | 120,960 |
| 16 | D12 + E8 | 12 | 241,920 |
| 17 | D4 + E8² | 9 | 120,960 |
| 18 | D20 | 14 | 241,920 |

Orbit sizes in this table come from the separate corner-walk certificate; the new computation checks the root-type identifications. Their sum is 1,095,990 literal cusps.

The independent checker takes each frozen δ, verifies its exact full active set, and projects incident roots to `δ⊥/<δ>` using the integral formula `K(r)=H*r+(r,h)*δ`, where `H=−(δ,h)>0`. It closes all incident norm-two roots under every incident reflection, checks all resulting norms, decomposes their nonorthogonality graph, and computes exact component ranks and root counts. Since every lattice norm-two root perpendicular to δ defines a reflection in its parabolic stabilizer, the complete incident reflection group generates this root system. The 18 results match the table bijectively. Rows 3,6,8 have norm-two rank 19,18,18; norm-four reflections supply the missing directions in their full bounded Coxeter links.

The script replays with:

```sh
python scripts/fresh_h22_extension_feasibility_20261013.py --verify
```

This provides a complete list of cusp **types**, while the separate [corner-walk script](../scripts/fresh_h21_source_corner_walk_20261013.py) and [independent replay](../scripts/fresh_h21_source_corner_walk_20261013_independent.py) determine literal orbit representatives, group actions, corners and edges. Their endpoint formulas and independent common-PD20-support argument were source-reviewed here. A fresh external run of the independent checker with `--verify` also passed in 33 seconds: 241,920 group elements, all 395 corner orbits, all 7,528 edge-orbit transitions and an empty queue, against final source SHA256 `c5da893e4ce6d88afcf705f001d92b55509fc18c824015d798bdda31a07b7423`. A global finite-volume extension still needs a compatible whole-P cover; repairing all old cusp links alone does not prove it.
