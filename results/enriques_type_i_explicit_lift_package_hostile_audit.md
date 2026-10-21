# Hostile audit of the explicit Type-I Enriques lift package

**Verdict:** `PASS/FREEZE`.

The explicit Type-I package survives an independent replay which does not
import the production certifier and does not use PARI `qfminim` as its
enumeration engine.  The hostile archive is
[`enriques_type_i_explicit_lift_package_hostile_audit.json`](../data/enumerations/enriques_type_i_explicit_lift_package_hostile_audit.json),
rebuilt by
[`audit_enriques_type_i_explicit_lift_package_hostile.py`](../scripts/audit_enriques_type_i_explicit_lift_package_hostile.py).

## Involution and glue

The raw OSCAR lattice chain and Shimada's literal `12A` record join exactly:
the first ten `SX` basis rows are `SY`, and

\[
\operatorname{Gram}(S_Y(2))=2\operatorname{Gram}(L_{10}).
\]

For row coordinates `z=(a,c)`, the independently reconstructed action is

\[
z\tau=(a+cD,-c),\qquad z(1+\tau)=(2a+cD,0).
\]

The audit explicitly distinguishes conventions:

```text
tau G tau^t = G: true
tau^t G tau = G: false (maximum error 152)
```

Thus a silent transpose cannot pass.  The integral anti-eigenbasis
`(-D | 2I)` has Smith diagonal `1^9` in `SX`.  Its Gram has determinant
`-1024`, Smith diagonal `2^8,4`, and every norm is divisible by four.  It is
therefore the full saturated rootless

\[
E_8(-2)\oplus\langle-4\rangle.
\]

The twelve selected anti-differences independently have Smith diagonal
`1^9` both in this anti-basis and as a sublattice of `SX`.

## Independent lift and L4 exhaustions

Instead of asking `qfminim` for the `E8+A1` roots, the hostile replay starts
with the nine simple roots and closes them under their exact simple
reflections.  It obtains 242 roots and exactly the same vector-set hash as
the production enumeration.  The parity equation then gives exactly two
lifts of each wall.  All twenty-four roots agree coordinate-for-coordinate
with the archive, and their unordered intersection histogram is

```text
0: 246
1:  28
2:   2
```

For `L4`, an exact rational LDL recursion enumerates the ellipsoid

\[
Q_H=-140\operatorname{Gram}(S_X)+2pp^t,
\qquad Q_H(e)\le378.
\]

It visits exactly 473 lattice points including zero, hence 472 nonzero
vectors, and its complete vector-set hash agrees with `qfminim`.  Every
nonzero vector is a square-minus-two root of absolute degree one through
seven.  There is no degree-zero root and no smaller-degree root pairing
negatively with any proposed lift.

This also covers isotropic residuals.  If a reducible lift has
`b=e+F`, with `e^2=-2`, `F^2=0`, and `e.F=0`, then `b.e=-2`; the enumerated
negative-component test detects `e` even though it never assumes that the
residual is another root.

## What licenses the surface

The archived `.mrdi` file is a lattice/Borcherds context, not equations for
an algebraic surface.  The hostile audit checked OSCAR 1.4.0's exact loader:

```julia
SY, SX, L26, w, u = load(...)
```

so the tuple order and Weyl-vector interpretation used by the certificate
are correct.  The loader source hash and semantics are now recorded in the
archived datum's
[`README.md`](../data/external/oscar_enriques_type_i/README.md).

Actual geometric existence is supplied by Brandhorst--Shimada, Section 7.4,
not by the byte stream alone.  That theorem identifies No. 172 with a very
general Type-I Enriques surface, its rank-19 K3 cover, its twelve smooth
rational curves, and its `12A` nef chamber.  OSCAR 1.4.0's own documentation
records the same example and the twelve rational-curve walls.

Accordingly, the package is an exact worked positive calibration on a known
surface and a known finite-volume noncompact `H9/12` polytope.  It is not a
new surface or abstract polytope, and it is not yet a downstairs-only
reconstruction of the embedding, glue, period, and `L4` data.
