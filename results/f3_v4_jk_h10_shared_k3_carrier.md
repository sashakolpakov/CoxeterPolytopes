# TOOT: one very-general complex-K3 family carries two large `H10` chambers

## Positive theorem

The already certified finite-volume noncompact chambers

- `F3-V4-H10-123`, with 123 facets, and
- `JK-H10-0007`, with 84 facets,

have the same integral Lorentzian wall lattice `S=U+E8(-2)+<-4>`, namely

\[
 S=U+E_8(-2)+\langle-4\rangle .
\]

There is an explicit primitive embedding

\[
 S\hookrightarrow \Lambda_{K3}=U^3+E_8(-1)^2.
\]

Consequently a nine-dimensional family of complex projective K3 surfaces
has very-general full Picard lattice `NS(X)=S`.  On every such very-general
surface, after choosing the recorded markings, **both complete oriented wall
systems act by integral reflections on the full divisor lattice**.  Thus two
known abstract Coxeter polytopes acquire one shared characteristic-zero,
full-Picard K3 mixed-divisor carrier family.

This is a positive algebraic realization, not a new abstract-polytope claim.

## Explicit K3-lattice embedding

Write the two ambient `E8(-1)` factors as `E8_a` and `E8_b`.  The embedding is

- the first `U` identically into the first ambient `U`;
- `E8(-2)` diagonally, `x -> (x,x)`, in `E8_a+E8_b`;
- `<-4>` by the primitive vector `e-2f` in the second ambient `U`.

Its orthogonal complement is

\[
 T=U+E_8(-2)+\langle4\rangle,
\]

of signature `(2,9)`.  The archived certificate checks the full 22-by-11
embedding matrix, a unit minor (hence primitivity), orthogonality, the
complement Gram matrix, and combined gluing index 1,024.  The discriminant
data are

```text
rank(S)       = 11
signature(S)  = (1,10)
det(S)        = 1024
A_S           = (Z/2)^8 + Z/4
length(A_S)   = 9
q_S           = q_plus(F2^8) + <-1/4>
```

The complement has the opposite discriminant form.  Period surjectivity and
avoidance of the countable extra-algebraic hyperplanes therefore give
`NS(X)=S` on the very general member, rather than merely an embedding into a
larger Picard lattice.

For the carrier census, the least `c>0` for which `S^vee(c)` is integral and
even is `c=8`.  That scaled dual has determinant `2^23` and discriminant
length 11.  It is not needed for this positive: the much smaller original
lattice `S` already embeds primitively.

## The two exact wall transports

| chamber | primitive mirror types | span in `S` | strict positive-cone interior | literal Vinberg audit |
|---|---|---:|---|---|
| `F3-V4-H10-123` | `10*(-2,div1) + 31*(-4,div2) + 82*(-4,div4)` | index 2 | square 4, minimum wall pairing 1 | 37,451 finite vertices, 1,732 cusps, 246,610 elliptic rank-9 edges, zero bad edges/neighbourhoods |
| `JK-H10-0007` | `8*(-2,div1) + 33*(-4,div2) + 43*(-4,div4)` | index 2 | square 8, minimum wall pairing 1 | 15,698 finite vertices, 862 cusps, 107,233 elliptic rank-9 edges, zero bad edges/neighbourhoods |

For the Fermat chamber, the archived determinant-one isometry to `S` is
replayed entry by entry.  Transporting its invariant ample vector fixes all
123 halfspace orientations; 64 archived primitive normals require a sign
change and 59 do not.  For the Jacobian--Kummer chamber, the saturated fixed
lattice, its complete discriminant form, and Nikulin's indefinite uniqueness
criterion identify it with `S`; all 84 primitive normals are reconstructed
from the published wall orbits and already pair positively with the
transported ample vector.

Every primitive reflection is independently checked for integrality,
involutivity, and preservation of the relevant Gram matrix.  In both rows
the displayed normals have full rank but generate an index-two sublattice;
the reflections nevertheless act integrally on **all** of `S`.  Finite
volume is not inferred from this lattice statement: the table records the
separate literal Vinberg edge-endpoint and full affine-neighbourhood audits,
each independently cross-checked by an exact GMP cone with zero hyperideal
rays.

## Exact boundary

These are mixed divisor-reflection carriers, not full-nef realizations.  The
large square-minus-four populations cannot be K3 negative-curve walls.  The
integral Picard reflections are also **not surface automorphisms** and are not
asserted to extend as Hodge isometries.  Riemann--Roch can choose effective
signs for the square-minus-two classes after Weyl framing, but no
irreducibility, complete-nef equality, or ample-intersection theorem is
claimed here.  Both abstract chambers and their arithmetic classes were
already known in the project ledger, and no worldwide-priority claim is made.

## Common-moduli upgrade

The two-example statement is now upgraded to a bounded inverse-carrier
criterion.  Self-pairing divisibility together with the full discriminant
form proves that the exhaustive primitive negative integral-mirror alphabet
of `S` is `(-2,div1), (-4,div2), (-4,div4)`: self-pairing divisibility
excludes `(-2,div4)`, while the discriminant form excludes `(-2,div2)` and
`(-8,div4)`.  Every genuine
finite-sided `H10` Coxeter polyhedron integrally realized in `S` is therefore
a weak period-level full-Picard divisor carrier on this locus.  The complete
13-wall group is nonuniform, so finite volume is equivalent to finite index
for genuine Coxeter fundamental polyhedra and compact chambers cannot occur
inside `S`.  See the
[common-moduli criterion](ue8m2a4_common_k3_moduli_carrier.md).

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_f3_v4_jk_h10_shared_k3_carrier.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_f3_v4_jk_h10_shared_k3_carrier.py
```

Frozen payload SHA-256:

`39574ca56421155388ec734076cf728a49f81c94bdfcdbc7d472b104e5014f08`
