# An H9 host outside the entire E8-plus-binary family

The next candidate is the literal lattice

\[
\boxed{L=\langle-\phi\rangle\perp I_7\perp\beta I_2,
\qquad \phi=(1+\sqrt5)/2,\quad\beta=4-\phi,\quad O=\mathbb Z[\phi].} \tag{1}
\]

It has signature (9,1), positive-definite conjugate, exactly four primitive reflection-norm classes **1,2,β,2β**, and complete zero chamber **B₇+B₂**. More decisively, it has local Witt index **three** at the prime (β), whereas every E₈-plus-binary quadratic space has local Witt index at least **four** there. Thus (1) is not similar over Q(√5) to any member of that entire family. Changing only a determinant in the earlier binary family does not have this property.

The [standalone Fraction checker](../scripts/fresh_h9_witt_twist_20261016.py) and [certificate](fresh_h9_witt_twist_20261016.json) verify the exact form, the old-family identification, the E₈ rational congruence and integral inverse, the residue-field facts, all 106 zero roots, reflection integrality and the nine simple walls. Generation and fresh saved replay pass. This is a construction setup, **not a reflective lattice or a finite-volume polytope**. No positive-height search is included.

## 1. What the excluded hosts have in common

The frozen [beta11](fresh_h9_beta11_corner_symmetries_20261014.md), [beta55](fresh_h9_beta55_corner_symmetries_20261014.md), and [beta131](fresh_h9_beta131_h300_20261016_symmetry.json) hosts have literal forms

\[
L_m=E_8\perp\begin{pmatrix}-2\phi&-1\\-1&2m\end{pmatrix},
\qquad m=1,2,3.
\]

Writing δ_m=1+4mφ gives det(L_m)=−δ_m and

\[
N(\delta_m)=1+4m-16m^2=-11,-55,-131.
\]

Their independently obtained infinite chamber stabilizers exclude those literal lattices. Their determinant ratios are nonsquares in the field, so they were not duplicates up to rational similarity. Nevertheless they retain the same large unimodular E₈ summand. The [norm31 silver host](fresh_norm31_continuation_20261015.md) has the analogous binary matrix with φ replaced by √2, over a different field. The separately [identified silver E₈ and F₄² presentations](fresh_silver_e8_f4_identification_20261016.md) are a literal duplicate; this does not identify all the other hosts.

Three parameter values do not prove nonreflectivity for every m, much less for all binary forms. No such extrapolation is made here.

There is, however, a useful exact structural restriction. Since E₈ has determinant one, any integral Gram matrix retaining its literal sublattice can be written and split as

\[
\begin{pmatrix}C&U\\U^T&B\end{pmatrix}
\ \cong_O\ C\perp(B-U^TC^{-1}U),\qquad C=G_{E_8}. \tag{2}
\]

The basis change has upper-right block −C⁻¹U and is integral. The checker reconstructs C and its integral inverse. Consequently arbitrary integral cross-couplings to this E₈ do not leave the E₈-plus-binary architecture.

For an even integral extension retaining the complete literal Bugaenko source \(E_8\perp\langle-2\phi\rangle\), the residual binary Gram is \([[-2\phi,c],[c,2d]]\). Relative basis changes put it into one of two forms, with c=0 or c=−1. Indeed O/(2)=F₄ and the residues of 1,φ,φ² are all its nonzero elements. If c is nonzero modulo2, choose a unit u among these with uc≡−1 modulo2 and set

\[
k=(uc+1)/(2\phi),\qquad t'=ut+kx.
\]

Then (x,t′)=−1. If c≡0, subtract the integral projection onto x. The new d is an arbitrary element of O satisfying admissibility, not necessarily an integer. Thus the tested m=1,2,3 lie on only one integer slice of the nonsplit gluing class.

This explains why a different positive block or nonintegral/nonarithmetic coupling is a more substantive change than another integral cross term on the same source. The local obstruction below implements the former route.

## 2. An exact local invariant that separates the new candidate

The element β=4−φ is totally positive and has norm11, so (β) is a prime ideal and

\[
O/(\beta)=\mathbb F_{11},\qquad \phi\equiv4\pmod\beta.
\]

Let K_β denote the completion. The residue4 is a nonzero square, so Hensel's lemma gives √φ∈K_β. Also −1 is not a square in F₁₁. The quaternary form

\[
A=\langle1,1,\beta,\beta\rangle \tag{3}
\]

is anisotropic over K_β. To prove this, normalize a hypothetical isotropic vector to have integral coordinates, at least one a unit. Reducing its equation moduloβ forces its first two coordinates to be divisible byβ, since x²+y²=0 over F₁₁ has only the zero solution. Divide the equation byβ; reduction now forces the last two coordinates to be divisible byβ as well. This contradicts the normalization.

The form I₄ is split over K_β. The vector (1,1,3,0) is a nonsingular isotropic residue vector modulo11, so it lifts to an isotropic local vector. After splitting off a hyperbolic plane, its binary orthogonal complement has determinant −1 modulo squares and is also hyperbolic. Thus I₄≅H². Because φ is locally a square, \(\langle1,-\phi\rangle\cong H\). Rearranging (1) now gives

\[
L\otimes K_\beta\cong I_4\perp\langle1,-\phi\rangle\perp A
\cong H^3\perp A.
\]

Its Witt index is exactly three.

The checker supplies a literal rational Euclidean realization of the E₈ simple roots, with determinant ±1. Therefore E₈ is isometric over K to I₈, and over K_β to H⁴. Any E₈-plus-binary space has Witt index at least four. Scalar similarity preserves Witt index, proving that (1) cannot be similar to **any** such space, regardless of its binary coefficients or determinant.

The determinant of (1) is −φβ², of norm −121. Its global determinant squareclass is nevertheless the same as that of the standard diagonal golden form \(\langle-\phi\rangle\perp I_9\). A determinant-only duplicate screen would therefore miss this distinction.

One geometric consequence is that a K-rational spherical E₈ root component cannot occur in this host. In a connected simply-laced E₈ root system, the label-three identities force the root norms to share a squareclass in K, so K-rescaling makes its form a scalar multiple of E₈. That would give a locally split eight-dimensional subspace, impossible in Witt index three. This blocks direct reuse of an E₈ corner transport underlying earlier experiments. It does not rule out other infinite chamber symmetries.

## 3. Admissibility and maximal integrality

At the physical embedding φ>0 and β>0, so (1) has signature (9,1). At the other embedding φ′<0 and β′>0, all ten diagonal entries are positive. In particular it has no nonzero K-rational isotropic vector: an isotropic vector would remain nonzero at the definite embedding.

The lattice is maximal among integral lattices in its quadratic space. Outside (β), its determinant is a unit and it is already self-dual locally. At (β), a vector in the dual modulo L has only its last two coordinates potentially fractional, of the form a/β,b/β. Integral square requires a²+b²≡0 moduloβ, forcing both residues zero. Therefore there is no proper integral overlattice. Here the determinant ideal has valuation two at β, but the discriminant-plane anisotropy prevents an overlattice; squarefreeness of the determinant is not being assumed.

## 4. Complete primitive reflection norms

The [frozen independent ring/norm review](fresh_h9_h8_cap_height_20261014_norm_pool_review.md) supplies the ring facts: O is a PID and its totally positive units are precisely squares of units. Let r be a primitive integral reflection vector, q=r²>0. At the conjugate embedding q′>0 as well. Reflection integrality says

\[
2Gr/q\in O^{10}.
\]

For the diagonal form (1), βG⁻¹ is integral, so 2βr/q is integral. A Bezout relation among the primitive coordinates gives q|2β. The primes (2) and (β) are distinct; the first is inert and the second has norm11. After rescaling r by a unit, the complete possible norms are

\[
\boxed{1,\quad2,\quad\beta,\quad2\beta}. \tag{4}
\]

Each occurs: a unit coordinate vector, a sum of two such vectors, a β-coordinate vector, and the sum of the two β-coordinate vectors. The proof uses the dual exponent, not a possibly incomplete short-root search. The checker records every entry of βG⁻¹ and verifies these norms on the complete zero shell.

## 5. Complete zero chamber at H=e₀

Use coordinates (x,y₁,…,y₇,z₁,z₂). A zero-height root has x=0. For q=1 or2, a nonzero z-coordinate would satisfy both βz²≤q and β′z′²≤q. Multiplication gives

\[
11N(z)^2\le q^2\le4,
\]

impossible for a nonzero algebraic integer. Thus z₁=z₂=0.

For q=β or2β, reflection integrality forces every y-coordinate to be divisible byβ. Write y_i=βu_i and divide the norm equation byβ. It becomes

\[
\beta\sum u_i^2+z_1^2+z_2^2=1\text{ or }2.
\]

The same norm-product argument forces every u_i=0. The two blocks therefore never mix in a zero-height reflection.

It remains to solve sums of ordinary O-squares equal to1 or2. Write each coordinate a_i+b_iφ. Equating coefficients gives

\[
\sum(a_i^2+b_i^2)=q,\qquad
\sum(2a_ib_i+b_i^2)=0.
\]

For q≤2 these force every b_i=0: one nonzero b gives an odd second sum; two nonzero b consume the full budget and give second sum2. The remaining coordinates are the ordinary integer roots. The seven-dimensional block has 14 norm-one and84 norm-two roots; the two-dimensional block has four roots of each of norms β and2β. Thus there are exactly **106 oriented zero roots, 53 mirrors**, with finite reflection type B₇+B₂.

The selected nine simple walls define

\[
y_1\ge y_2\ge\cdots\ge y_7\ge0,
\qquad z_1\ge z_2\ge0.
\]

Their vectors, Gram matrix and all 36 relations are serialized. This is a complete finite corner at H of the full reflection chamber, because every nonzero-height simple wall has negative pairing with H; local finiteness prevents accumulation of new mirrors at H. It does not determine the chamber away from that corner.

All root counts are unit-normalized to the four norms (4), with opposite directions counted separately. Multiplication by arbitrary units gives infinitely many equivalent integral normal vectors, not additional mirrors.

## 6. Next experiment and provenance limits

The next concrete experiment is an ordered Vinberg pilot on (1), retaining all four classes (4) and the complete B₇+B₂ zero sector. For positive height, x>0 and the dominant y,z coordinates have the ordering above. The norm equation and height are

\[
\sum y_i^2+\beta\sum z_j^2=q+\phi x^2,
\qquad h=\phi x^2/q.
\]

Both embeddings bound every coordinate. No positive-height roots are certified by this setup. The separate pilot must prove its own control-domain completeness, chronological selection and geometry.

A bounded provenance check found no literal matching form in the inspected repository paths. That is not a literature-wide absence or novelty assertion. The original [Bugaenko chapter](https://doi.org/10.1090/advsov/008/02) is the provenance for the related reflective-lattice constructions, not a theorem being extended to this nonunimodular local form. [Scharlau–Walhorn, §§1–2](https://www.numdam.org/item/AST_1992__209__279_0.pdf) illustrates the use of dual exponents and local quadratic-form structure, but its stated classification scope is over Z in small dimensions; it does not classify (1). All local and norm arguments needed here are given explicitly above.

The current result is an exact new search direction with a proven structural distinction and a complete starting chamber. It neither claims that nonreflectivity is avoided nor transfers any old nonreflectivity theorem by determinant or commensurability alone.

```sh
python3 scripts/fresh_h9_witt_twist_20261016.py
python3 scripts/fresh_h9_witt_twist_20261016.py --verify
```
