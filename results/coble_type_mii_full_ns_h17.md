# Type MII closes the positive Coble endpoint in `H17`

**Status (2026-08-09): exact full-nef theorem, with known abstract-polytope
lineage.**  Kondō's unique Type-MII Coble surface in characteristic three has
a projectivized full Néron--Severi nef cone which is a `48`-facet,
finite-volume, noncompact arithmetic Coxeter polytope in `H17`.  Its complete
integral wall Gram is exactly the previously known Keum--Kondō/Alexeev
`(18,4,0)` chamber.  A fully independent exact double-description audit of
the complete 48-wall cone gives
`196,052=190,808+5,244+0` finite/ideal/hyperideal rays, all on one causal
sheet.

This is a substantial positive result for the surface track and realizes the
sharp algebraic-geometric endpoint of the nine classified finite-automorphism
Coble classes.  It is **not** a new abstract Coxeter polytope, and it receives
no `POLYTOPE -> SURFACE` theorem credit until the raw-polytope compiler covers
the entire nine-class screening denominator.

## 1. The theorem

Let `S` be Kondō's Type-MII Coble surface over an algebraically closed field
of characteristic three.  Then

\[
 \rho(S)=18,qquad K_S^2=-8,qquad
 \dim \mathbf P(\operatorname{Nef}S)=17.
\]

The complete set of negative curves on `S` consists of

| curves | number | square |
|---|---:|---:|
| anti-bicanonical boundary components | 8 | `-4` |
| exceptional curves above the `F_3` grid | 16 | `-1` |
| proper transforms of graphs of permutations | 24 | `-2` |

These `48` curves are pairwise acute and their normalized Gram has

```text
796*m2 + 128*m3 + 96*m4 + 108*parallel.
```

Every wall is irredundant.  The resulting full-nef chamber is the
finite-volume noncompact `H17/48` chamber of the even 2-elementary lattice

\[
 L=2\widetilde{\operatorname{Pic}}(S)
   \simeq U\oplus D_8(-1)^{\oplus2},qquad
 (r,a,\delta)=(18,4,0).
\]

All `48` facet reflections preserve `L`; the eight boundary reflections, and
only those eight, fail to preserve the unextended odd unimodular Picard
lattice.

## 2. The surface model fixes the Picard glue

Kondō realizes `S` as the blow-up of
`P1 x P1` at all sixteen points of
`P1(F3) x P1(F3)`.  In the basis

\[
 F_1,F_2,E_{11},\ldots,E_{44},qquad
 F_1F_2=1,quad E_{ij}^2=-1,
\]

the eight boundary classes are

\[
 B_i=F_1-\sum_jE_{ij},\qquad
 B'_j=F_2-\sum_iE_{ij}.
\]

They are pairwise orthogonal of square `-4`, and

\[
 \sum_i(B_i+B'_i)=-2K_S.
\]

For every permutation `sigma in S4`, the proper transform

\[
 C_\sigma=F_1+F_2-\sum_iE_{i,\sigma(i)}
\]

has square `-2`.  The `48` displayed physical curves span
`Pic(S)=I_(1,17)` with Smith diagonal `1^18`; there is consequently no
unresolved full-Picard glue choice.  Adding the eight boundary half-classes
has index `2^7`, since their sum is `-K_S`.  After doubling the form the
extended lattice has determinant `-16`, Smith diagonal `1^14 2^4`, and
`delta=0`.

The primitive arithmetic normals are

```text
B_i/2, B'_j/2, E_ij, C_sigma,
```

with square/divisibility split

```text
24*(-2,div1) + 24*(-4,div2).
```

They span the whole arithmetic lattice, and all their reflection matrices
are integral involutive isometries.

## 3. Why there are no hidden negative curves

The `40` Coble--Mukai roots are

\[
 e_{ij}=B_i/2+B'_j/2+2E_{ij}quad(16),qquad C_\sigma\quad(24).
\]

An exact saturation computation gives

\[
 \operatorname{CM}(S)=E_{10}=II_{1,9};
\]

the forty roots span this lattice integrally.  Exact rational double
description of their chamber gives

```text
4,954 extreme rays = 4,752 finite + 202 ideal + 0 hyperideal.
```

The `202` primitive cusps form four orbits under
`(S4 x S4) semidirect C2`:

| orbit size | maximal affine type |
|---:|---|
| 72 | `A7~ + 2*A1~` |
| 96 | `A5~ + A2~ + 2*A1~` |
| 18 | `A3~ + A3~ + 2*A1~ + 2*A1~` |
| 16 | `A2~ + A2~ + A2~ + A2~` |

Here `2*A1~` denotes the doubled affine-`A1` edge in Kondō's notation,
not two components.  Thus the enumeration independently recovers exactly
the four maximal-parabolic types in Kondō's `parabolicMII` lemma.

Kondō proves that his forty classes are all effective irreducible roots.
Consequently they account for every `(-2)`-curve and every `(-1)`-curve
meeting two distinct boundaries.  A remaining `(-1)`-curve `E` would have
`B.E=2` by adjunction and would have to meet a single boundary `B_i` twice.
Kondō's Halphen argument then makes

\[
 f=E+B_i/2
\]

a primitive isotropic Coble--Mukai class in the chamber, with inverse lift
`E=f-B_i/2`.  The certificate tests every one of the `202*8=1,616` possible
primitive cusp/boundary lifts in the full Picard lattice.  None is integral.
There are therefore no further `(-1)`-curves, and the displayed `48` curves
are the complete negative-curve set.

This is the essential full-nef gate.  Merely recognizing the forty
Coble--Mukai roots would only have produced a codimension-eight section.

## 4. Coxeter and irredundancy checks

The physical wall Gram has diagonal split

```text
8*(-4) + 16*(-1) + 24*(-2)
```

and nonnegative cross-pairings only.  The four normalized possibilities are
`m=2,3,4,infinity`; there are no non-Coxeter angles or ultraparallel facet
pairs.

There is also an exact witness for every facet.  The integral class

\[
 H=5F_1+5F_2-\sum_{i,j}E_{ij},\qquad H^2=34,
\]

pairs `1` with each boundary and exceptional curve and `6` with every
`C_sigma`.  For a wall `w`, its Lorentz-orthogonal projection

\[
 H_w=H-\frac{H.w}{w^2}w
\]

is timelike, pairs zero with `w`, and pairs strictly positively with all
other `47` walls because the Gram is acute.  Hence every inequality supports
a genuine codimension-one face.

## 5. Exact finite volume, not an inference from finite sidedness

There are three complementary checks.

First, negative-curve exhaustion identifies the chamber with the full nef
cone, and Hodge index forces every nef ray into the closed future cone.  This
is the geometric consistency check; it is not used alone as a substitute for
Vinberg's criterion.

Second, an independent hostile reconstruction of the complete 48-wall cone,
without importing this certificate or its wall dictionary, gives

```text
196,052 extreme rays = 190,808 finite + 5,244 ideal + 0 hyperideal.
```

Every ray pairs positively with a strict timelike chamber point, so all rays
lie on one causal sheet.  This is a direct exact finite-volume proof, and its
`5,244` ideal rays independently reproduce Alexeev's maximal-parabolic count.

Third, the certificate constructs an explicit determinant-one integral
isometry taking all `48` arithmetic normals, entry by entry, to Alexeev's
published complete `(18,4,0)` simple-root system.  Alexeev proves

\[
 O^+(L)=\operatorname{Aut}(K_{4,4})\ltimes W_r,qquad
 |\operatorname{Aut}(K_{4,4})|=1152,
\]

so the reflection group has finite covolume.  His Vinberg algorithm
terminates after two steps and his exact maximal-parabolic enumeration gives

```text
5,244 maximal parabolic subdiagrams in 17 symmetry orbits.
```

Thus the chamber is finite-volume and noncompact by two independent routes.
This avoids the invalid
shortcut “finite walls implies finite volume” and also avoids an unnecessary
dependence on the raw enumeration for the original proof.

## 6. Exact lineage

The wall dictionary is literal:

| Type-MII geometry | Alexeev `(18,4,0)` diagram |
|---|---|
| 8 boundary half-roots | 8 vertex roots of `K4,4` |
| 16 exceptional curves | 16 edge roots subdividing `K4,4` |
| 24 graph curves `C_sigma` | 24 long roots indexed by perfect matchings |

Keum and Kondō described this abstract chamber in 2001 using the
Borcherds--Leech method.  Alexeev gave the different
`K4,4^(2)`/perfect-matching presentation, proved reflectivity, and enumerated
the maximal parabolics.  Kondō's 2021 paper constructs the Type-MII surface
and publishes all of the curves from which the full-NS realization follows.

Accordingly:

- the surface and abstract chamber are known;
- the arithmetic lattice is known;
- this project adds a complete, replayable full-NS/cusp-lift certificate and
  an explicit marked isometry between the two published presentations;
- worldwide priority for that interpretation is not asserted.

The result is nevertheless the desired positive endpoint for the Coble
surface track: the unique algebraic-geometric stratum with eight boundary
components realizes `H17`, exactly the sharp dimension bound supplied by
Kondō's nine-class classification.

## 7. Reproduction

The exact archive is
[coble_type_mii_full_ns_h17.json](../data/enumerations/coble_type_mii_full_ns_h17.json),
generated by
[certify_coble_type_mii_full_ns_h17.py](../scripts/certify_coble_type_mii_full_ns_h17.py).
It stores the full physical and arithmetic Grams, all `202` primitive cusps,
all lift decisions, all facet witnesses, and the marked Alexeev isometry.
The independent reconstruction is documented in the
[hostile audit](coble_type_mii_h17_hostile_audit.md).

Run from the repository root:

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_type_mii_full_ns_h17.py --check-archive
pytest -q tests/test_coble_type_mii_full_ns_h17.py
```

Primary sources:

- Shigeyuki Kondō,
  [*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591).
- Valery Alexeev,
  [*Reflective hyperbolic 2-elementary lattices, K3 surfaces and hyperkahler varieties*](https://arxiv.org/abs/2209.09110).
- JongHae Keum and Shigeyuki Kondō,
  [*The automorphism groups of Kummer surfaces associated with the product of two elliptic curves*](https://doi.org/10.1090/S0002-9947-00-02631-3).
