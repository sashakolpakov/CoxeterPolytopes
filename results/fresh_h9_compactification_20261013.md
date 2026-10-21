# Compactification feasibility and a changed compact-face seed

**No compact H9 polytope was found.** This round first checked whether the proposed discrete-angle compactification would repeat an existing search. It then built a different exact starting seed: an eleven-facet H9 Coxeter chamber with a compact, non-quasi-arithmetic H5 face and a transverse A4 chamber. A bounded family of312 individually compatible caps still leaves one exact strict spacelike escape, even when every cap is imposed simultaneously. This is a finite-family obstruction, not an obstruction to completing the new seed.

## Why the old compactification lane was not repeated

The repository already examined the two doubled18-facet noncompact H9 chambers T1/T2. Its [cusp-filling report](fresh_h9_cusp_filling_20261013.md) includes unchanged-affine obstructions and a numerical search changing up to four old angles, including1,383,498 mutations after symmetry reduction. Those results do not exclude arbitrary angle changes or altered facet incidences.

A separate exact stock of [finite-volume H9 sections](alexeev_parked_rank10_h9_examples.md) contains the following chambers. The incidence counts below are read directly from its frozen certificate.

| Source | Facets | Finite vertices | Ideal incidences |
|---|---:|---:|---|
| A14_8_0_C01 |12|44|5 cusps with10 incident facets|
| A14_8_1_C02 |16|165|14 with10,9 with11|
| A14_10_1_C01 |13|70|7 with10|
| A14_10_1_C02 |15|123|10 with10,5 with11|
| A16_14_1_C01 |15|100|4 with10,5 with11|

A compact Coxeter H9 polytope is simple: each finite vertex has nine incident facets. Thus retaining the complete old face incidence structure cannot turn any of these cusps into a finite vertex; one would have to split its incidences. The few-facet compact classifications also exclude the12/13-facet targets. In particular, compact Coxeter d-polytopes with d+4 facets have dimension at most7. [Felikson–Tumarkin](https://arxiv.org/abs/math/0510238).

This is a feasibility triage, not a new exclusion of the15/16-facet angle-surgery problem. No new numerical Gram surgery on those three sources is claimed. The chosen computation instead implements the alternative compact-face construction requested for this round.

## The changed compact H5 face

The prior seed used the arithmetic prism P0=[5,3,3,3,3] with transverse H4. Both P0 and P2=[5,3,3,3,4] are arithmetic individually; P2 is not a nonarithmetic example by itself. [Emery–Kellerhals, Sections1 and4](https://homeweb.unifr.ch/kellerha/pub/IML-1112s-21.pdf).

Here the two prisms are glued along their isometric orthogonal hyperbolic four-simplex facets. Gluing along such facets is the mechanism underlying the mixed Makarov garlands studied in2024. The particular two-block geometry below is verified directly rather than identified with a numbered published garland. No new H5 example is claimed. [Bogachev–Douba–Raimbault, Sections3.3 and4](https://arxiv.org/html/2309.07691v3).

Write `phi=(1+sqrt(5))/2`, `K=Q(sqrt(5))`. Let G be the norm-two five-node chain with consecutive off-diagonal entries `-phi,-1,-1,-1`. Put

```
d=(G^-1)[4,4]=(1-sqrt(5))/4,
g=-G^-1 e4,
D0=2-d=(7+sqrt(5))/4,
R=(4-4d)/D0=(16+4sqrt(5))/11.
```

The H5 form is `Q5=G ⊕ [D0]`. Its seven normals are `(e_i,0)` for `i=0,...,4`, together with

```
r5=(g,1),           r5²=2,
r6=(2g,-sqrt(R)),   r6²=4.
```

The five side walls extend across the gluing seam. The two end walls have label3 and label4, respectively, against side wall4; their mutual hyperplanes are disjoint. The exact program checks the full Gram matrix, signature(5,1), all21 candidate five-facet intersections, and precisely ten timelike extreme rays. A strict timelike point is also expressed as a **strictly positive combination of all seven normals**. Its affine slice is bounded, so the complete ray census proves compactness; it does not merely certify ten local vertices.

The end-wall Cartan two-cycle is

```
4*(r5,r6)^2/(r5²*r6²)
 = (4+sqrt(5)) + (-1+3sqrt(5))*sqrt(R)/4.
```

Together with the label5 cycle this generates `L=K(sqrt(R))`. The automorphism changing the sign of `sqrt(R)` fixes the indefinite five-wall principal Gram block G. Vinberg's semidefiniteness criterion therefore excludes quasi-arithmeticity of this compact H5 reflection group, and conditionally excludes quasi-arithmeticity of any finite-volume H9 completion retaining these seed facets. This is not a claim that all compact H9 polytopes are nonarithmetic.

## H9 seed and the bounded global test

Add the orthogonal norm-two A4 Cartan block H to Q5. The resulting form has signature(9,1), and the seed has eleven actual Coxeter facets. The ten H5 vertices become actual finite H9 vertices when the four transverse walls are included. Actual facets are certified by projecting the strict timelike center onto each wall and checking strict inequalities against every other wall.

For each of the ten compact-face vertices v and each transverse label profile in `{2,3,5}^4`, let

```
b_i = 2*cos(pi/m_i), with b_i=0 when m_i=2,
s = b^T H^-1 b - 2,
k = s/(-v²),
r = (sqrt(k)*v, -H^-1 b).
```

Only `s>0` defines the tested real timelike-center family. The new normal then has norm2 and the four prescribed transverse angles. This is a bounded810-proposal family; free centers and other labels are outside its scope.

The exact census gives:

| Outcome | Count |
|---|---:|
| Nonpositive old scale |60|
| Forbidden angle with an old H5 wall |438|
| Individually compatible with all eleven seed walls |312|

For these specific profiles every squared old cosine belongs to L. The quartic field L is not normal: `Norm_K/Q(R)=16/11`, and neither this number nor its quotient by5 is a rational square. Consequently its only nontrivial cyclotomic subfield is K. A finite squared cosine in L can therefore have only labels2,3,4,5,6,10. The exact rejection gate uses this complete alphabet; it is not a truncation at label10 for arbitrary fields.

The integer point

```
E=(30861,39915,29906,19910,9932,-705,-7565,-10211,-10432,-10624)
```

is future-oriented and strict against all eleven seed normals. Its square is

```
E²=(11370214691-4926770235*sqrt(5))/4 > 0.
```

For each of the312 caps, the program checks the sign of `a*sqrt(k)+z`, where `a=(v,E_old)<0`, by exact signs and one squared comparison. All312 are strictly negative. There is also a common strict timelike interior: every cap is strict at `(h,0)`, since its old part is future timelike, and a sufficiently small inward A4 perturbation makes the four transverse seed inequalities strict while preserving the finitely many cap inequalities. Hence their entire intersection has infinite volume, regardless of their mutual compatibility. Every subset has the same escape. No compatible-selection search can repair this bounded bank.

A separate **numerical-only** scout uses transverse labels2,3,4,5,6,8,10 and old finite-angle comparisons2 through100. Its24,010 literal proposals yield16,985 accepted entries. The same integer E remains strict numerically against this wider bank. The wider result is exploratory; the exact theorem covers only the810 specified proposals above.

## Reproduction and remaining choice

The exact files are [script](../scripts/fresh_h9_compactification_20261013.py) and [certificate](fresh_h9_compactification_20261013.json). The separate bounded numerical scout is [script](../scripts/fresh_h9_compactification_20261013_scout.py) and [output](fresh_h9_compactification_20261013_scout.json).

```
python scripts/fresh_h9_compactification_20261013.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h9_compactification_20261013_scout.py
```

The exact checker reuses the frozen Q(sqrt5) scalar implementation with a hash-bound dependency, while implementing the new quadratic extension and geometry directly. It is not an independent arithmetic replay. A separate [source/proof review and fresh saved replay](fresh_h9_compactification_20261013_review.md) passed against the frozen files; that review likewise does not claim separate arithmetic implementation.

The concrete remaining options are incidence-splitting Gram surgery on the15/16-facet noncompact sources, or free-center/mixed-contact caps on this changed compact face. The computation above supplies a compact starting face with different arithmetic and transverse geometry, while rejecting its first bounded global cap family before any facet-growth campaign.
