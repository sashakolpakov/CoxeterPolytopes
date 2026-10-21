# Exact lineage of the 328-wall characteristic-three H11 chamber

## Result

The 328-wall Fermat-quartic section is a genuine finite-volume noncompact
hyperbolic Coxeter polytope.  Its canonical integral reflection lattice and
its full lattice-reflection group, however, are already present in the
Scharlau--Kirschmer census.

More precisely, the minimal lattice containing the saturated divisor-fixed
lattice and preserved by all 328 facet reflections is explicitly integrally
isometric to entry 20 of Scharlau--Kirschmer List11/Gram11:

\[
\mathrm{II}_{11,1}(2^4 3^{-1}),\qquad |\det|=48.
\]

That census entry has a 22-wall full integral reflection chamber, 543 finite
vertices, and 22 cusps.  Consequently the 328-wall reflection group is a
finite-index subgroup of a known full reflection group.

This settles three priority levels negatively:

- the stable integral reflection lattice is not new;
- its full integral reflection group is not new;
- its rational arithmetic commensurability class is not new.

The positive result remains substantial: this is an exact 328-wall
finite-volume Coxeter chamber extracted from a complete divisor-wall stock.
Priority of that particular subgroup chamber, and of its Fermat-quartic
realization, remains open after a bounded corpus audit.

## The divisor lattice and the reflection lattice are different

The saturated common fixed lattice in the Neron--Severi lattice has signature
(1,11), determinant 768, and Smith diagonal

\[
1^6\,2^4\,4\,12.
\]

Of the 328 primitive facet mirrors, 292 give integral reflections of this
fixed divisor lattice and 36 do not.  The exceptional mirrors have fixed
lattice type \((-24,\operatorname{div}6)\) and reflection denominator two.
Thus it would be incorrect to say that all 328 reflections act integrally on
the saturated divisor-fixed lattice itself.

Closing that lattice under the 328 rational reflections stabilizes after one
step.  The result is the canonical minimal reflection-stable overlattice:

- index over the divisor-fixed lattice: 4;
- signature: (1,11);
- even and integral;
- determinant: 48;
- Smith diagonal: \(1^8\,2^3\,6\);
- discriminant group: \((\mathbf Z/2)^3\oplus\mathbf Z/6\).

All 328 mirrors are primitive integral roots there, with exact types

\[
276(-6,\operatorname{div}3)
+28(-2,\operatorname{div}1)
+24(-4,\operatorname{div}2).
\]

This distinction is important for the algebro-geometric program.  The
surface supplies the hyperbolic divisor-fixed lattice and the wall
inequalities, while the natural integral lattice for the complete induced
reflection group is a canonical finite overlattice.  No surface-automorphism
claim is needed or made.

## Exact Scharlau--Kirschmer identity

The maintained
[Scharlau--Kirschmer reflective-lattice census](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html)
explains that List11 records the face numbers and cusps of the fundamental
domain of the full reflection group \(W(L)\), while Gram11 gives a Gram matrix
for the corresponding class-number-one genus.

The source archive is frozen by the following hashes:

| object | SHA-256 |
|---|---|
| Lorentz.zip | 74a4dd3e534f1e793cba45765ff33575dbd331d5cdce6f34cdd64340b3f0d4a5 |
| List11 | 2ca0725232cb101e60380f14ad0ab7777478c60eff06991ba591b6d3c71e3184 |
| Gram11 | af8d52984fc9caabc00ffe670ae646d2151e68d141082544e18f4827ba142540 |
| canonical entry-20 Gram | c854fc82d4602c0e6cb3752d76e0a195e51194532754feaab85d2b92b99e455f |

The lineage certificate archives a determinant-minus-one integer matrix
\(R\) and verifies entry by entry that

\[
R^T G_{\mathrm{stable}}R=-G_{\mathrm{SK},20}.
\]

The minus sign only converts Scharlau--Kirschmer's signature-(11,1)
convention to the signature-(1,11) convention used in the surface
calculation.  This is a full integral isometry, not a determinant/Smith
collision.

Entry 20 is

\[
(48,\ \mathrm{II}_{11,1}(2^4 3^{-1}),\
[22,204,1089,3789,9126,15660,19295,16851,9973,3564,543],\
22).
\]

The census convention makes the first number in the face list the number of
facets and the last the number of finite vertices.  Including its 22 ideal
vertices, the full face vector from vertices through the polytope is

\[
(565,3564,9973,16851,19295,15660,9126,3789,1089,204,22,1).
\]

The census archive does not distribute the 22 simple-root coordinates or
their weighted Gram matrix.  Therefore this certificate proves the exact
lattice identity, the identity of the full integral reflection group
\(W(L)\), and its published face counts, but it does not pretend to give an
entrywise 22-root diagram comparison.

## The 328-wall subgroup

Every one of the 328 roots was transported through the explicit isometry to
entry-20 coordinates.  In both bases its reflection matrix is integral,
involutive, preserves the Gram matrix, and negates the defining root.
Therefore

\[
W_{328}<W(L_{\mathrm{SK},20})
\]

is a literal reflection-subgroup inclusion.  The exact GMP cone certificate
proves that the 328-wall chamber has

- 225,408 finite vertices;
- 11,554 cusps;
- zero hyperideal extreme rays.

It is therefore finite-volume noncompact.  Since the census chamber is also
finite-covolume, the subgroup index is finite.

The numerical index is not obtainable from a Solomon--Euler ratio.  In
hyperbolic dimension 11 both rational Euler characteristics vanish, so the
ratio is \(0/0\).  Computing the index requires explicit coset/tiling data or
an odd-dimensional hyperbolic-volume comparison.

## Why the tempting orthogonal decomposition is wrong

The rank, determinant, parity, and Smith diagonal suggest

\[
U\oplus D_8(-1)\oplus A_2(-2).
\]

That guess is false.  The certificate enumerates all 48 elements of the
discriminant group and computes their quadratic values in
\(\mathbf Q/2\mathbf Z\).  The value histograms are

| lattice | exact \(q\)-value histogram |
|---|---|
| stable Fermat lattice | \(0^{10},1^6,(1/3)^{12},(4/3)^{20}\) |
| \(U+D_8(-1)+A_2(-2)\) | \(0^6,1^{10},(2/3)^{12},(5/3)^{20}\) |

Thus their discriminant quadratic forms differ.  This gives a small,
replayable integral obstruction and illustrates why determinant plus Smith
normal form was not enough to establish lineage.

## McLeod, Alexeev, and CoxIter comparisons

[McLeod's \(f_3^n\) paper](https://arxiv.org/abs/1007.2299) gives a reflective
15-wall chamber in dimension 11 for
\(-3x_0^2+x_1^2+\cdots+x_{11}^2\).  Exact determinant-square and local
Hasse-invariant calculations show that the stable Fermat lattice has the same
rational quadratic space as \(f_3^{11}\).  The two integral lattices are not
isometric: their determinants are 48 and 3.  Hence the McLeod comparison
closes the arithmetic commensurability lineage but does not replace the exact
Scharlau--Kirschmer integral identification.

Alexeev's 2-elementary census does not apply: the stable discriminant group
has a nontrivial 3-primary part.

The pinned CoxIter snapshot has 1,214 ordinary records, with complete H11
facet histogram

\[
12:14,\quad 13:6,\quad 15:2.
\]

Its separate 2,959-record reflection-doubling archive has H11 histogram

\[
15:2,\quad16:2,\quad18:2,\quad22:2,\quad30:2,\quad46:2.
\]

The two 22-wall files there are doubling descendants of the \(f_2^{11}\)
record, not source records for the determinant-48 Scharlau--Kirschmer
lattice.  Most importantly for subgroup priority, neither pinned corpus
contains a 328-facet H11 record.  This facet-count exclusion is bounded; it
is not a worldwide weighted-diagram priority theorem.

## Current status

| identity level | status |
|---|---|
| hyperbolic Coxeter construction | proved finite-volume noncompact |
| saturated divisor-fixed lattice | determinant 768; not preserved by 36 mirrors |
| canonical integral reflection lattice | known: Scharlau--Kirschmer List11/Gram11 entry 20 |
| full integral reflection chamber | known: 22 walls, 543 finite vertices, 22 cusps |
| arithmetic commensurability class | known: rationally the McLeod \(f_3^{11}\) space |
| literal 328-wall subgroup chamber | absent from two pinned corpora; global priority open |
| Fermat-quartic divisor-wall realization | positive; priority open |

## Reproducibility

- certificate:
  [ss_fermat_quartic_char3_v4_h11_lineage.json](../data/enumerations/ss_fermat_quartic_char3_v4_h11_lineage.json)
- generator:
  [certify_ss_fermat_char3_v4_h11_lineage.py](../scripts/certify_ss_fermat_char3_v4_h11_lineage.py)
- tests:
  [test_ss_fermat_char3_v4_h11_lineage.py](../tests/test_ss_fermat_char3_v4_h11_lineage.py)

    PYTHONPATH=src:. python3 \
      scripts/certify_ss_fermat_char3_v4_h11_lineage.py --check-archive

    PYTHONPATH=src:. python3 \
      scripts/certify_ss_fermat_char3_v4_h11_lineage.py \
      --check-archive --verify-external

    PYTHONPATH=src:. pytest -q \
      tests/test_ss_fermat_char3_v4_h11_lineage.py
