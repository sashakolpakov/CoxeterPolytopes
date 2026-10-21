# The first ACT real-cubic `H10` slice has a global mod-three obstruction

**Status (2026-08-09): proved infinite-covolume obstruction for the one
explicit antiunitary class used by the exploratory search.  No new
finite-volume hyperbolic polytope and no algebro-geometric leverage delta.**

The executable certificate is
[`certify_act_real_cubic_h10_congruence_obstruction.py`](../scripts/certify_act_real_cubic_h10_congruence_obstruction.py),
with frozen
[`JSON`](../data/enumerations/act_real_cubic_h10_congruence_obstruction.json)
and
[`tests`](../tests/test_act_real_cubic_h10_congruence_obstruction.py).  A
[`no-project-import replay`](../scripts/independent_act_real_cubic_h10_congruence_replay.py)
independently implements rational-pair Eisenstein arithmetic, matrix
inversion, and literal dual-number products for all `3^11` residues.

## Convention audit

Allcock--Carlson--Toledo's Theorem 1.7 gives the Eisenstein lattice

```text
Lambda = (3) + EE8 + EE8 + [[0,theta],[bar(theta),0]]
```

and explicitly uses the convention

```text
h(x,y) = x^T H conjugate(y), linear in x.
```

For the antiunitary map used here--reverse both `EE8` chains, swap the final
pair, and conjugate--the displayed integral fixed basis is correct.  The
exact calculation gives fixed Gram matrix

```text
2 Q,       Q = diag(1,1,1,1,1,1,1,1,1,1,-1).
```

This checks the convention which had previously been easy to reverse by
mistake.  The primary source is
[Allcock--Carlson--Toledo, Theorem 1.7](https://arxiv.org/abs/math/0608287).

## The exact extension code

Let `r=(x0,...,x9,t)` be a vector of the fixed lattice.  Extending its real
reflection to an integral Eisenstein isometry reduces modulo `3` to five
linear equations:

```text
x1 + x3 + x4 = 0
x2 + x3 - x4 = 0
x5 + x7 + x8 = 0
x6 + x7 - x8 = 0
x9 - t         = 0                  (mod 3).
```

This is not merely a necessary screen.  The certificate exhausts all
`3^11=177,147` residues in
`E/(3)=F3[epsilon]/(epsilon^2)` and proves equality with the original
rank-one matrix integrality condition.  Exactly `729=3^6` residues survive.
They form

```text
C = F3*e0 + T4 + T4 + <e9+t>,
```

where `T4` is the ternary tetracode: one word of weight zero and eight words
of weight three.  The Lorentz form restricted to `C` has rank one:

```text
Q(r)   = x0(r)^2                         (mod 3),
B(r,s) = x0(r) x0(s)                     (mod 3).
```

These two identities replace the failed large face enumeration.

## Global Vinberg obstruction

For a primitive vector of the integral fixed lattice, reflection integrality
forces `Q(r)` to divide every `2 r_i`.  Primitivity therefore leaves only

```text
Q(r) in {1,2}.
```

The code identity excludes `Q(r)=2` modulo three.  Consequently every
extendable reflection root has norm one and `x0(r)` is nonzero modulo three.
For any two such roots,

```text
B(r,s) = x0(r)x0(s) != 0                 (mod 3).
```

Thus the **complete**, all-height root set contains no orthogonal pair.
For an acute simple-root diagram, every distinct integral norm-one pairing
is consequently at most `-1`.  A two-node principal Gram matrix is then
never positive definite.  Hence:

- the maximum elliptic rank is one;
- a positive-semidefinite affine subsystem can be at most the two-node
  `affine A1`, of Coxeter rank one; and
- there cannot be several affine components, since distinct components
  would have to be orthogonal.

Vinberg's literal first finite-volume gate in `H10` requires at least one
elliptic rank-ten finite vertex or one componentwise-affine cusp of total
Coxeter rank nine.  Both lists are empty before a height search begins.
Therefore this antiunitary extension branch has no compact or noncompact
finite-volume Coxeter chamber.

This is stronger than saying that a truncated outer cone is hyperideal.
Later roots cannot repair a congruence obstruction holding for every root.

## Reusable design gate

The negative theorem compiles into an inexpensive positive-search filter.
For any antiunitary fixed slice whose primitive reflectors are forced to have
norm one, compute the exact extension code `C` and the restriction of the
Lorentz pairing to its norm-one locus before running Vinberg.  A finite
vertex in `Hn` requires `n` mutually orthogonal unit roots, while an ideal
vertex requires a componentwise-affine system of total Coxeter rank `n-1`.
If neither residue pattern exists, the class is globally rejected.  A new
ACT antiunitary class should enter wall enumeration only after passing this
finite code gate.  In particular, a promising class must have substantially
more orthogonality on its norm-one residue locus than the rank-one form found
here.

## OOM-safe symmetry compression and bounded corroboration

The tetracode has monomial automorphism group of order `48`.  Fixing the
Vinberg basepoint and the oriented height-zero simple root leaves

```text
(MonAut(T4) x MonAut(T4)) semidirect C2,
```

of order `4,608`.  Canonicalizing each tetracode block independently and
then swapping the two blocks compresses the bounded prefix as follows:

| height | walls at height | symmetry orbits at height |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 1 | 1 |
| 2 | 17 | 2 |
| 3 | 561 | 6 |
| **total** | **580** | **10** |

Every one of the ten computed orbits is complete: its theoretical orbit
size equals its literal multiplicity, and the sizes sum to `580`.

A final bounded exact `cdd` replay on these 580 walls, run under the new
2-GiB parent RSS watchdog, returned

```text
finite rays       0
ideal rays        0
hyperideal rays   1,680.
```

This is corroboration only.  The proof is the all-height ternary-code and
Vinberg-rank argument above.  The abandoned 5,300-wall conversion is neither
needed nor permitted; Gram summaries are now streamed instead of allocating
the full square matrix.

## Exact scope

The theorem closes **this explicit antiunitary lattice class** as a source of
finite-volume reflection polytopes.  It does not classify all antiunitary
classes of the ACT lattice, and it does not prove that this particular class
is realized by a component of real cubic-threefold moduli.  Therefore it is
a genuine P/lattice-side negative theorem, but not a new example and not a
`POLYTOPE -> SURFACE -> POLYTOPE` algebraic-geometric pushback theorem.

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_act_real_cubic_h10_congruence_obstruction.py --check-archive
PYTHONPATH=src:. python3 scripts/independent_act_real_cubic_h10_congruence_replay.py
PYTHONPATH=src:. pytest -q tests/test_act_real_cubic_h10_congruence_obstruction.py
```
