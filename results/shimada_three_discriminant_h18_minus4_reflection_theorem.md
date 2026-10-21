# Uniform primitive-minus-four theorem for the Shimada H18 families

## Theorem

In the completed full oriented-wall searches of Shimada's pinned
discriminant-11, discriminant-15, and discriminant-16 induced-chamber
atlases, there are exactly 934 **fresh** subgroup sections in dimension
`H18`:

| ambient Neron--Severi discriminant | H18 chamber-local classes | distinct primitive vectors in the pinned coordinates |
|---:|---:|---:|
| 11 | 216 | 67 |
| 15 | 144 | 56 |
| 16 | 574 | 109 |
| **total** | **934** | **232 within the three separate lattices** |

For every class, the nonidentity rational primal action is exactly

\[
 s_v=I-\frac{2}{(v,v)}Gv^{\mathsf T}v
     =I+\frac12Gv^{\mathsf T}v,
\]

where `G` is the ambient rank-20 divisor Gram matrix and `v` is a primitive
divisor vector satisfying

\[
 (v,v)=-4,\qquad \operatorname{div}(v)=
 \gcd\{(v,x):x\in NS\}=1.
\]

Thus every one of the 934 apparently separate high-rank wall symmetries is
the same algebraic-lattice mechanism: reflection in a primitive
divisibility-one `(-4)` divisor class.  The action has denominator exactly
two in every case and is never an integral Neron--Severi isometry.

The fixed divisor lattice is

\[
 K_v=NS\cap v^\perp.
\]

For a primitive vector in a nondegenerate integral lattice,

\[
 |\operatorname{disc}K_v|
 =|\operatorname{disc}NS|\frac{|v^2|}{\operatorname{div}(v)^2}.
\]

Consequently these three families have fixed-lattice discriminants exactly
`44`, `60`, and `64`.  Their Smith diagonals are respectively

```text
d=11: 1^18,44
d=15: 1^18,60
d=16: 1^16,2,4,8
```

Every one of the 934 sections has at least one forbidden exact normalized
Gram product between two independently certified genuine hyperbolic facets.
Therefore none is an ordinary hyperbolic Coxeter polytope.  This last claim
is an angle obstruction, not an infinite-volume assertion; the volume gate
is neither needed nor invoked after a case fails to be Coxeter.

## Exact proof certificate

For each class the certificate reconstructs the nonidentity rational
20-by-20 action and performs all of the following over `Q` or `Z`:

1. verifies that the subgroup is `C2`, the action squares to the identity,
   preserves the ambient Gram matrix, and has rank-one difference from the
   identity;
2. primitive-normalizes the nonzero row space of `s-I` to recover `v`;
3. verifies `v*s=-v`, `v^2=-4`, and `div(v)=1`;
4. proves entrywise that the action equals the reflection formula above and
   has denominator two with a genuinely nonintegral entry;
5. replays the archived saturated basis of `NS intersect Fix_Q(s)`, proves it
   equals a rank-19 sublattice of `v^perp`, recomputes its Gram matrix,
   determinant, and Smith invariants, and checks the discriminant formula;
6. links the class to its exact no-go certificate and verifies that at least
   one forbidden pair was certified.

All 934 records pass every check.  The minimum forbidden-pair count is one;
there are 38 cases attaining that minimum.

## Algebro-geometric interpretation

This theorem keeps the divisor-lattice idea alive but sharply identifies
what the present fishing denominator is producing.  Hyperbolicity is free
from the Hodge index signature, and a primitive negative divisor class gives
a canonical rational hyperplane section.  Here the full wall-diagram
symmetry search repeatedly rediscovers precisely the `(-4)`, divisibility-one
version of that construction.

It does **not** produce a K3 automorphism: for `v^2=-4` the reflection is

\[
 x\longmapsto x+\frac{(x,v)}2v,
\]

and divisibility one guarantees that some integral divisor has odd pairing
with `v`.  Hence the denominator-two obstruction is intrinsic, rather than
an accident of coordinates.  No surface-automorphism hypothesis was used as
a fishing gate; the construction lives naturally in the rational divisor
space and its saturated integral hyperplane lattice.

The reduction from 934 chamber-local casts to only 232 primitive vectors is
also operationally useful.  Future variants can enumerate negative divisor
orbits first, build `v^perp` and its restricted nef inequalities, and apply
the exact Coxeter gate once per divisor orbit instead of rediscovering the
same rational reflection chamber by chamber.

## Scope

The quantifier is exact but bounded:

- all 216 fresh `H18` classes in the complete pinned determinant-11
  oriented-wall denominator;
- all 144 fresh `H18` classes in the complete pinned determinant-15
  oriented-wall denominator; and
- all 574 `H18` classes in the complete pinned determinant-16 fresh
  strict-extension denominator.

It is not a theorem about every primitive `(-4)` class on every K3 surface,
and it is not a universal non-existence theorem for `H18` Coxeter
polytopes.  It explains and closes the largest uniform family encountered by
these three complete algebraic wall-stock searches.

## Reproduction

```bash
PYTHONPATH=.:src python3 scripts/certify_shimada_three_discriminant_h18_minus4_reflection_family.py
PYTHONPATH=.:src pytest -q tests/test_shimada_three_discriminant_h18_minus4_reflection_family.py
```

Frozen certificate:

- `data/enumerations/shimada_three_discriminant_h18_minus4_reflection_family.json`

