# A changed cubic-seven H9 host: three dihedral blocks and A3

The host `L=<-2τ>+B³+A3` adds arithmetic coverage beyond the six recent cubic-seven `E8+<1>` and `B⁴+<1>` searches. Here `τ=2cos(2π/7)`, `C=τ²+τ−1`, `β=2−τ`, `D=2−τ²`, and `B=[[2,-C],[-C,2]]`. A complete search through height 100 gives **67 Coxeter walls with an exact spacelike escape**. No finite-volume example was found; the full host remains open.

The positive block changes, rather than adding further walls to the previous 25-wall host. The time coefficient is negative at the distinguished embedding and positive at the two others, while B and A3 are positive definite at all three embeddings.

## Equivalence checks and their limits

The determinant is `−8τβ³`. Its ratio to the previous `<-2τ>+B⁴+<1>` determinant is `4/β`, a nonsquare. The [independent audit](fresh_cubic7_mixed_blocks_h100_independent_audit_20261013.json) compares all three earlier coefficients `τ,2τ,τβ` in both the E8 and four-dihedral families. The field norms of the six determinant ratios are

```
175616, 21952, 25088, 512/7, 64/7, 512/49.
```

Each is a nonsquare rational number. Therefore the corresponding field ratios are nonsquares, excluding isometry and similarity of these rank-ten quadratic spaces. This is a comparison with those six saved hosts, not a worldwide novelty or noncommensurability classification against all prior constructions.

There is a useful internal identification. The columns of

```
P = [[0,1,0],[1,-1,1],[-1,0,1]]
```

satisfy `PᵀP=A3`, and their image is `{z∈O³: z1+z2+z3∈2O}`. Thus L has the integral index-`(2)` overclosure `M=<-2τ>+B³+I3`, of index cardinality eight. This is its unique integral overclosure of that index. Indeed the dyadic discriminant module is `O/(2)⊕O/(4)`; an order-two class is `a e0/2+bω2`, with norm `−τa²/2+b²`. Integrality forces `a=0 mod2`, leaving the single A3 line. Conversely L is intrinsically `{v∈M:q(v)∈2O}`. Consequently **O(L)=O(M)**, including their full reflection groups. These are two lattice presentations of one group, not two discoveries. The audit checks the matrix and the eight possible time residue classes explicitly.

## Complete primitive norms and zero-height chamber

The ring is `O=Z[τ]`, of class number one, and totally positive units are squares; the [earlier field audit](fresh_cubic7_dihedral_independent_audit_20261013.md) supplies the proof. Since `4βG⁻¹` is integral, a primitive reflection norm divides `8β`, up to square units. The form is even, removing odd norm classes. Norms 8 and 8β are imprimitive: reflection integrality makes the time coordinate even, every B coordinate divisible by four, and the A3 covector divisible by four. For A3 coordinates `(z1,z2,z3)`, this implies `z2=0 mod2` and `z1=z3 mod2`. Its half-norm is `2z1² mod4`; the full norm equation forces this to vanish modulo four, hence all three A3 coordinates are even. The audit exhausts the finite residue ring `O/(4)` and finds exactly eight allowed triples, all even.

The complete primitive norm pool is therefore **`{2,4,2β,4β}`**. Norms 2 and 2β are squarefree ideals. For norms 4 and 4β, checking that the vector is not divisible by two establishes primitivity.

The zero-height chamber is **`I₂(14)³×B3`**, with nine simple walls. It is derived as follows. For norm 2, divide the spatial equation by two: every nonzero B or A3 half-norm is a totally positive algebraic integer of trace at least three. Their sum is one, so only one block can occur. For norm 2β, reflection integrality puts B vectors in `adj(B)O²` and the A3 vector in `βO³`; the same trace argument leaves one B block. At norm 4, each B vector is even and contributes at least six to the trace after division by four, so only A3 remains. At norm 4β the same argument eliminates B, after which β-divisibility of A3 contradicts the proposed norm valuation.

The independently enumerated root census is:

| Norm | Oriented zero-height roots | Dominant representatives |
|---|---:|---:|
| 2 | 54 | 4 |
| 4 | 6 | 1 |
| 2β | 42 | 3 |
| 4β | 0 | 0 |
| Total | 102 | 8 |

The A3 part has twelve norm-two and six primitive norm-four roots, all rational. Its simple normals are `−e1,−e2,e1−e3`, with Coxeter labels 3 and 4. In covector coordinates its chamber is `u1≥0,u2≥0,u3≥u1`. Each B chamber is `0≤u≤w`. All root orbits lie on one side or the other of these simple-wall cones; no mixed zero-height mirrors are omitted.

## Bounded search and exact geometry

For a normalized norm q and height `h=2τx²/q≤H`, the time bounds are `x²≤Hq/(2τ)` and `x'²≤−q'/(2τ')` at the other embeddings. The spatial target is `q+2τx²`. Cauchy's inequality bounds every B or A3 covector square by twice that target at each embedding. Exact trace-dual bounds give finite coefficient alphabets. Inverse Gram matrices reconstruct integral block coordinates; exact energy joins cover every block placement, including zeros and equalities.

The [height-30 calibration](fresh_cubic7_mixed_blocks_h9_h30_20261013.json) has fourteen controls and 26 walls, with escape norm `30−8τ−12τ²>0`. The [height-100 artifact](fresh_cubic7_mixed_blocks_h9_h100_20261013.json) has 31 controls and 108 candidate roots before chronological wall selection, resulting in 67 walls. Its exact escape is

```
x=1;
B1=B2=(0,0),  B3=(1+τ−τ²,1+τ−τ²);
A3=τ²*(1/2,1,1/2);
q(v)=−7−13τ+15τ² > 0.
```

All 2211 wall pairs pass the exact Coxeter-angle gate; the forbidden acute angle `2π/7` is explicitly excluded. Strict timelike interiors and strict spacelike perturbations certify infinite volume of the saved prefixes. The [independent verifier](../scripts/fresh_cubic7_mixed_blocks_independent_verify_20261013.py) uses standard-library Fraction arithmetic, brute-force coefficient boxes without the producer's interval pruning, and a Cartesian spatial join instead of its multiset/permutation join. It independently repeats **every height-100 candidate enumeration and acceptance decision**, as well as the zero census, integral reflections, angles, and geometry. Its PASS run took 34.84 seconds; the producer took 10.11 seconds.

```sh
PYTHONPATH=src:. python3 scripts/fresh_cubic7_mixed_blocks_search_20261013.py --height 100 --output results/fresh_cubic7_mixed_blocks_h9_h100_20261013.json
python3 scripts/fresh_cubic7_mixed_blocks_independent_verify_20261013.py --source results/fresh_cubic7_mixed_blocks_h9_h100_20261013.json --output results/fresh_cubic7_mixed_blocks_h100_independent_audit_20261013.json
```

The escapes prove failure of these finite prefixes. They do not exclude later mirrors, other wall selections, or a finite-volume chamber of this full reflection group.
