# Silver `N(D12^2)` exact three-coefficient screen

## Outcome

No compact or finite-volume Coxeter polytope is obtained from the three
tested Lorentzian lattices

\[
 \langle-a\rangle\oplus N(D_{12}^2),\qquad
 a\in\{1+\sqrt2,\sqrt2,2\sqrt2\},
\]

in \(\mathbf H^{24}\).  This is nevertheless a genuinely new exact
Niemeier-host screen: the enumerator uses the full index-four glue, rather
than replacing `N(D12^2)` by its root sublattice.

The unit branch is complete through real height 30 and has 99 accepted
positive-height walls.  The `sqrt(2)` and `2sqrt(2)` branches are complete
through real heights 40 and 27 and have respectively 19 and 27 accepted walls.
Each frozen prefix has an explicit strict hyperideal ray.  These are bounded
negative screens, not global nonreflectivity theorems.

## The correct index-four glue

Write the discriminant classes of `D12` as `0,v,s,c`, with

\[
 v=(1,0),\qquad s=(0,1),\qquad c=(1,1).
\]

Their minimum squared norms are `1,3,3`.  The diagonal code cannot define
the Niemeier lattice with root system `D12^2`: it contains `(v,v)` of square
two and would add roots.  Up to exchanging the two spinor labels, the
root-system-preserving glue is instead

\[
 \mathcal G=\{(0,0),(v,s),(s,v),(c,c)\}.
\]

For standard `D12` Dynkin weights `z_1,...,z_12`, the two residue bits are

\[
 \rho_v=z_1+z_3+z_5+z_7+z_9+z_{12}\pmod2,
 \qquad
 \rho_s=z_{11}+z_{12}\pmod2.
\]

The backend imposes `rho(right)=swap(rho(left))` separately on the rational
and `sqrt(2)` coefficients.  It therefore tests membership in
`N(D12^2) tensor Z[sqrt(2)]`, not merely in `D12^2`.

An explicit 24-vector Euclidean basis is constructed from the two glue
generators `(v,s)` and `(s,v)`.  Its Gram matrix is integral, even, and has
determinant one.  The height-zero dominant census has exactly two norm-two
orbits and no glue root, hence exactly the advertised root system with

\[
 2\cdot 12(2\cdot12-2)=528
\]

roots.

## Complete primitive mirror norms

Put `pi=sqrt(2)`.  For a primitive mirror `r=(x,y)`, reflection integrality
is

\[
 q(r)\mid 2I(r),\qquad I(r)=(a x,(y,N)).
\]

Because the spatial block is even unimodular, local primitivity makes
`I(r)` a unit away from `(pi)`.  The exact dyadic valuation argument leaves
the following complete canonical lists:

| coefficient `a` | canonical primitive mirror norms |
|:---|:---|
| `1+sqrt(2)` | `1, 2` |
| `sqrt(2)` | `2, 2+sqrt(2)` |
| `2sqrt(2)` | `2, 4+2sqrt(2)` |

Totally positive units are squares, so no further unit representatives are
missing.  At every shell, divisibility by `sqrt(2)` is checked in the full
glued lattice; imprimitive copies are deleted before chronological
acceptance.

## Unit branch through height 30

For `a=1+sqrt(2)`, the exact simultaneous packets are:

| height | controls | candidates | rejected | accepted | cross-glue accepted |
|:---|---:|---:|---:|---:|---:|
| `(1+sqrt(2))/2` | 1 | 0 | 0 | 0 | 0 |
| `1+sqrt(2)` | 2 | 0 | 0 | 0 | 0 |
| `2+2sqrt(2)` | 2 | 2 | 0 | 2 | 0 |
| `(7+5sqrt(2))/2` | 1 | 0 | 0 | 0 | 0 |
| `7+5sqrt(2)` | 2 | 8 | 0 | 8 | 4 |
| `(17+13sqrt(2))/2` | 1 | 0 | 0 | 0 | 0 |
| `(23+17sqrt(2))/2` | 1 | 0 | 0 | 0 | 0 |
| `14+10sqrt(2)` | 2 | 93 | 4 | 89 | 33 |

The first two walls are the two factorwise affine `D12` caps.  Their product
has `13^2=169` projective rays: 3 finite, no ideal, and 166 hyperideal.
The later glue is substantial, especially in the 89-wall third packet, but
does not close the cone.

Take time coordinate one and the single nonzero Dynkin weight

\[
 z_1=\frac85
\]

in the first `D12` factor.  Since `(C_D12^-1)_{11}=1`, its Lorentz norm is

\[
 \frac{64}{25}-(1+\sqrt2)=\frac{39}{25}-\sqrt2>0.
\]

It pairs strictly negatively with all 99 positive-height walls; the closest
pairing is

\[
 -\frac65+\frac45\sqrt2<0.
\]

The height-zero inequalities are nonpositive as well, so this is an exact
strict hyperideal ray through the complete height-30 prefix.

## The more restrictive `sqrt(2)` branch

The complete chronology through real height 40 is:

| height | candidates | rejected | accepted | cross-glue accepted |
|:---|---:|---:|---:|---:|
| `sqrt(2)/2` | 0 | 0 | 0 | 0 |
| `1+sqrt(2)` | 2 | 0 | 2 | 0 |
| `2+3sqrt(2)/2` | 0 | 0 | 0 | 0 |
| `2+2sqrt(2)` | 0 | 0 | 0 | 0 |
| `4+3sqrt(2)` | 0 | 0 | 0 | 0 |
| `7+5sqrt(2)` | 0 | 0 | 0 | 0 |
| `8+6sqrt(2)` | 17 | 0 | 17 | 5 |
| `12+17sqrt(2)/2` | 0 | 0 | 0 | 0 |
| `14+10sqrt(2)` | 0 | 0 | 0 | 0 |
| `16+12sqrt(2)` | 0 | 0 | 0 | 0 |

The sparse ray supported on Dynkin node 3 of the first factor with weight
`7/10` has norm

\[
 3\left(\frac7{10}\right)^2-\sqrt2
 =\frac{147}{100}-\sqrt2>0.
\]

All 19 positive-height pairings are strict; the closest is

\[
 -\frac35+\frac25\sqrt2<0.
\]

Thus the coefficient that gives the tighter first affine cap still fails
the finite-volume gate.

## The wall-rich `2sqrt(2)` branch

This branch is complete through real height 27 and has the richest early
chronology of the three:

| height | candidates | rejected | accepted | cross-glue accepted |
|:---|---:|---:|---:|---:|
| `1+sqrt(2)` | 4 | 0 | 4 | 0 |
| `2+2sqrt(2)` | 0 | 0 | 0 | 0 |
| `4+3sqrt(2)` | 0 | 0 | 0 | 0 |
| `7+5sqrt(2)` | 8 | 0 | 8 | 2 |
| `8+6sqrt(2)` | 17 | 2 | 15 | 5 |

The first packet contains the two universal affine caps and two additional
root-lattice walls.  Consequently later enumeration uses only inequalities
whose defining walls were explicitly found; the affine cap is not assumed
from the root-system proxy.

At Dynkin node 5 of the first factor, weight `7/9` gives

\[
 5\left(\frac79\right)^2-2\sqrt2
 =\frac{245}{81}-2\sqrt2>0.
\]

It is strict against all 27 accepted walls, with closest pairing

\[
 -\frac19-\frac19\sqrt2<0.
\]

Hence even this wall-richer early branch is not finite-volume at the frozen
cutoff.

## Exact enumeration method

Each `D12` factor is enumerated in dominant Dynkin-weight coordinates.  An
exact Schur-complement lower bound prunes partial weights simultaneously at
both real embeddings.  Once the factorwise affine caps have actually been
accepted, the exact inequality

\[
 \sum_i m_i z_i\le c_a x
\]

is added, where `m=(1,2,...,2,1,1)` and `c_a` is respectively
`2,sqrt(2),2`.  Equal-energy buckets join the two factors, the full glue
decoder filters the join, and all earlier wall pairings are then evaluated
in `Z[sqrt(2)]`.  Same-height survivors are checked pairwise acute before
being inserted together.

The JSON certificate stores every accepted root, its two blocks of Dynkin
weights and simple-root coordinates, its discriminant residues, and stable
root transcript hashes.

## Reproduction

- Backend: `src/coxeter_polytopes/search/silver_d12_squared.py`
- Audit: `scripts/audit_silver_d12_squared_h24.py`
- Tests: `tests/test_silver_d12_squared.py`
- Certificate: `data/enumerations/silver_d12_squared_h24.json`

```console
PYTHONPATH=src python3 scripts/audit_silver_d12_squared_h24.py
PYTHONPATH=src pytest -q tests/test_silver_d12_squared.py
```

The focused suite passes `6/6` tests; the extended final replay took 104.09 seconds.  The certificate
SHA-256 is
`3e85959c630f9abd5d7dc9516a553800d262de90722f2f90ca7364bc1d6653f5`.
