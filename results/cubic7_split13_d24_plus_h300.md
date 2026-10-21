# Cubic-seven split-prime-13 `D24+` prefixes through H=300

## Outcome

Let

\[
K=\mathbb Q(\theta),\qquad \theta=2\cos(2\pi/7),
\]

and let (p_7,p_{10},p_8) be the totally positive generators

\[
p_7=-\theta+2\theta^2,qquad
p_{10}=1-\theta+\theta^2,qquad
p_8=1+\theta^2
\]

of the three primes over 13, labelled by the corresponding residue of
(	heta\pmod {13}).  The exact searches cover

\[
L_i=\langle-\theta p_i\rangle\oplus D_{24}^{+}
\]

through Vinberg height (H=300).  None of the three exact prefix cones has
finite volume, so there is no compact or finite-volume survivor in this
bounded search.

The two best prefixes are (p_{10}) and (p_8), each with 25 walls and
exact ray census

```text
finite / cusp / hyperideal = 2 / 0 / 23.
```

The (p_7) direction acquires four further walls but becomes much more
combinatorially open: its final 29-wall cone has 895 rays with census
`3 / 0 / 892`.

This is strictly an (H\le300) result.  Later Vinberg walls may cut every
displayed hyperideal ray, so no global nonreflectivity claim is made.

## Why all three directions are retained

The three prime ideals are Galois conjugate, but the Lorentzian search fixes
embedding zero.  After choosing totally positive representatives modulo unit
squares, multiplication by the fixed (	heta) gives the distinct timelike
coefficients

```text
r7:   theta*p = 2+4theta-3theta^2,
r10:  theta*p = 1+3theta-2theta^2,
r8:   theta*p = 1+3theta-theta^2.
```

Each is positive at the distinguished place and negative at the other two.
Consequently the three searches cannot be collapsed by an unqualified
Galois identification.

## Complete local mirror norms

For a primitive vector (r=(x,y)), set

\[
I(r)=(\theta p x,(y,D_{24}^{+})).
\]

Since `D24+` is even unimodular, reflection integrality gives

\[
(q(r))\mid 2I(r).
\]

Primitivity and self-duality exclude every prime away from the inert prime
`(2)` and the selected split prime `(p)`.  At `(p)`, either (y) is locally
primitive and (v_p(q)=0), or (y\in pD_{24}^{+}), (x) is a unit, and

\[
q(r)\equiv-\theta p x^2\pmod {p^2},
\]

so (v_p(q)=1).  Evenness and the same divisibility condition give
(v_2(q)\in\{0,1\}).  The field has class number one and every totally
positive unit is a square.  Thus, for each of the three distinguished
directions, the complete canonical mirror norms are

\[
\boxed{1,\ 2,\ p,\ 2p}.
\]

For `q=p` and `q=2p`, reflection integrality is equivalent to
(y\in pD_{24}^{+}).  The implementation writes (y=pz) and exhausts
(z\in D_{24}^{+}) at the reduced square (y^2/p^2).  A nonintegral
reduced target is therefore an exact emptiness certificate, not a heuristic
pruning rule.

## Glue and height-zero audit

The spatial Gram matrix is even, integral, and unimodular of determinant one.
The enumerator covers the `D24` root coset and all seven nonzero classes of
`O_K/2O_K` in the spinor coset.  In a nonzero spinor class, all 24 doubled
coordinates have the same nonzero residue, no zero padding is allowed, and
both dominant last-coordinate signs are retained.

For every one of the three coefficients and every norm in `1,2,p,2p`, an
independent height-zero enumeration finds exactly one dominant reflective
orbit: the norm-two `D24` root `(1,1,0^22)`.  Hence the 24 classical `D24`
simple roots give the complete height-zero chamber.

## Exact control and wall census

All compatible controls occur at distinct exact heights.

| direction | norm 1 | norm 2 | norm p | norm 2p | total |
|---|---:|---:|---:|---:|---:|
| `r7` | 2 | 7 | 9 | 28 | 46 |
| `r10` | 3 | 6 | 9 | 28 | 46 |
| `r8` | 3 | 7 | 9 | 28 | 47 |

Only four packets across the three searches are nonempty:

| direction | exact height | candidates | accepted | accepted norms/cosets |
|---|---|---:|---:|---|
| `r7` | `19+54theta+24theta^2` | 1 | 1 | norm 1, `D` |
| `r7` | `38+108theta+48theta^2` | 5 | 4 | norm 2: three `D`, one spinor residue 1 |
| `r10` | `33+94theta+42theta^2` | 1 | 1 | norm `p10`, `D` |
| `r8` | `6+18theta+8theta^2` | 1 | 1 | norm `2p8`, `D` |

The accepted (p_{10})- and (2p_8)-walls explicitly pass the spatial
divisibility audit after division by their respective prime generators.
At every nonempty height, all survivors are tested against lower-height walls
first and are then verified pairwise nonpositive before simultaneous
insertion.

## Exact cone comparison

Exact double description uses rank 23 of common active functionals for
adjacency.  Every archived hyperideal witness has nonpositive pairing with
all prefix walls, positive Lorentz norm, and active-functional rank 24.

| direction | packet | walls | rays | finite | cusp | hyperideal |
|---|---:|---:|---:|---:|---:|---:|
| `r7` | 1 | 25 | 25 | 1 | 0 | 24 |
| `r7` | 2 | 29 | 895 | 3 | 0 | 892 |
| `r10` | 1 | 25 | 25 | 2 | 0 | 23 |
| `r8` | 1 | 25 | 25 | 2 | 0 | 23 |

Therefore none of the three (H\le300) chambers is finite-volume or
compact.  Among these bounded prefixes, `r10` and `r8` tie for the smallest
hyperideal-ray count; `r7` is not competitive after its second packet.

## Reproducibility

- Backend: `src/coxeter_polytopes/search/cubic7_d24_plus_split13.py`
- Generator: `scripts/audit_cubic7_split13_d24_plus_h300.py`
- Tests: `tests/test_cubic7_split13_d24_plus.py`
- Archive: `data/enumerations/cubic7_split13_d24_plus_h300.json`

The archive contains every exact control and accepted root, packet/glue
labels, same-height acuteness checks, exact cone snapshots and ray-stream
hashes, one independently replayable rank-24 hyperideal witness per seed,
and the precise H=300 frontier.
