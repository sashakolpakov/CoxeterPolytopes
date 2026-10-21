# Any three nondotted chain contacts preserve the H20 weak escape

Any coherent Coxeter cap that cuts the fixed weak escape `E0` must have **at least 24 dotted contacts among the 26 retained chain walls**. This is a necessary condition on an escape-cutting cap. It does not say that every cap has 24 dotted contacts.

The [standalone exact checker](../scripts/fresh_h20_arbitrary_triple_escape_gate_20261016.py) and [certificate](fresh_h20_arbitrary_triple_escape_gate_20261016.json) cover all 300 positive gap pairs, representing all 2,600 choices of three chain walls. Generation and fresh saved replay pass. Arithmetic uses only standard-library integers and fractions; no label cutoff or sampled cap bank enters the proof. Independent review of this extension is pending.

The retained chamber is the [134-wall H20 construction](fresh_h20_dotted_arc_cap_20261016_review.md). Its Lorentzian three-dimensional chain complement, rational positive support Gram inverse, and physical parameter isolation are inherited. The [independently reviewed consecutive-triple theorem](fresh_h20_asymmetric_arc_escape_gate_20261016_review.md) supplies two explicit dependencies: the degree-50 theorem over `K Q^ab`, and the weak feasibility of `E0` for the retained 134 walls. The present checker redoes the support-budget inequalities and proves the arbitrary-gap extension independently of that theorem's consecutive-triple algebra.

## The coefficient identity for arbitrary gaps

Normalize roots to squared norm 2. Write a proposed cap as

```
r = y - sum_i d_i omega_i,       d_i >= 0,
q_y = y^2 = 2 - d^T J d.
```

Here `J` is the rational inverse Gram matrix of the 18 support roots, and `y` lies in their orthogonal Lorentzian three-space. The support is `A10+A4+A2+2A1`; `J` is entrywise nonnegative and has minimum diagonal `1/2`.

For two chain roots separated by `k`, their pairing is

```
g_k(C) = 2*(T_k(C)-C)/(1-C),
C = a^2/(a^2+2),
g_1=0,        g_2=-2-4C.
```

The quotient is an integer polynomial of degree at most `k-1`. Three chain indices, in increasing order, have gaps `p,q >= 1`, with `p+q <= 25`. Their Gram matrix has off-diagonal entries `g_p,g_q,g_(p+q)` and diagonal 2. The checker certifies its negative determinant for every gap pair, so the three roots form a basis of the Lorentzian complement.

Suppose their cap contacts are `-f1,-f2,-f3`, with each `f_i` nondotted, hence in `[0,2]` and in `Q^ab`. If every support contact were also nondotted, then `q_y` would lie in `Q^ab`. The equation

```
q_y*det(G) - f^T*adj(G)*f = 0
```

would therefore be a polynomial equation over `Q^ab` of degree at most `2(p+q)-2 <= 48`. The inherited theorem gives degree 50 over `K Q^ab`, and thus no nonzero polynomial over `Q^ab` of degree below 50 can vanish at `C`. Consequently this polynomial must vanish identically.

The certificate saves all seven coefficient columns, exact rational row-reduction pivots, and a full kernel basis in the order

```
(q_y, f1^2, f2^2, f3^2, f1*f2, f1*f3, f2*f3).
```

The complete classification is:

| Gap class | Cases | Kernel dimension | Consequence for actual squares and products |
|---|---:|---:|---|
| Unequal gaps, both at least 2 | 242 | 0 | All contacts zero |
| `p=1,q>=3` | 22 | 1 | `f2^2+f3^2=0`, then all zero |
| `q=1,p>=3` | 22 | 1 | `f1^2+f2^2=0`, then all zero |
| Equal gaps at least 2 | 11 | 2 | `f1^2+f3^2=0`, `f2=0` |
| `(1,2)` | 1 | 2 | `f2=f3=0`; the product `f1*f3=0` kills the remaining parameter |
| `(2,1)` | 1 | 2 | `f1=f2=0`; the product `f1*f3=0` kills the remaining parameter |
| `(1,1)` | 1 | 4 | `f1=f3=0`, `q_y=f2^2/2` |

If all three contacts vanish, independence gives `y=0`. Then the nonzero norm-2 cap is `r=-sum_i (Jd)_i n_i`, a negative nonnegative combination of support roots. It pairs positively with any strict interior point of the retained chamber, so it cannot be a coherent additional wall.

The sole nonzero algebraic exception is consecutive: `y=-f2*n_middle/2`. A retained chain wall at distance two from the middle wall has positive pairing `f2*(1+2C)`. Such a wall exists for every one of the 24 consecutive triples. Hence the exception violates acuteness. **At least one support contact of every coherent cap with three nondotted chain contacts must be dotted.** This argument covers arbitrary real cap coefficients; it does not assume the cap is arithmetic.

## A uniform norm bound from the outer two contacts

The checker evaluates the chain polynomials on a rational isolating interval for `C` and proves

```
g_25 < ... < g_3 < g_2 < -2.
```

This is a statement about the retained finite chain, not about arbitrary larger gaps. The outer pair of any triple has a separation `n>=2`; its span has signature `(1,1)`, with positive-definite orthogonal complement inside the Lorentzian three-space. Orthogonal projection therefore gives, for its two contact magnitudes `f,h` in `[0,2]`,

```
u = -y^2
  <= (2*f^2+2*h^2-2*g_n*f*h)/(g_n^2-4)
  <= 8/(-g_n-2)
  <= 2/C = 2+4/a^2.
```

The middle contact is unnecessary for this metric bound. The maximum of the positive numerator on the square occurs at `f=h=2`.

The forced dotted support contact has magnitude greater than 2, and `J_ii>=1/2`, so `d^T J d>2`. Thus `u>0`. The outer-root sum is future timelike and pairs nonpositively with `y`; consequently `y` is future timelike too. This supplies exactly the hypotheses needed by the common-escape estimate.

## The common escape and complete support budget

Use the weak spacelike escape

```
E0 = x - e/100 - omega5,
E0^2 = 190011/110000 > 0.
```

Its support pairings are `-1` at root 5 and zero at the other 17 support roots. Its pairings with all other retained walls are strictly negative, as inherited from the preceding theorem. One way to verify this inheritance is to start with the independently reviewed strict escape `Eeps` and use

```
E0 = Eeps + (1/100000)*sum_(j != 5) omega_j.
```

Every added dual weight pairs nonpositively with an old wall outside the support, by acuteness and `J>=0`.

The new norm bound gives the total support cost

```
d^T J d = 2+u <= 4+4/a^2 < 281/70.
```

If a dotted contact lies in `A10`, only endpoint 1 or 10 can fit the budget. Its cost is greater than `40/11`, leaving less than `1/2`; hence every other support contact is zero. The checker proves both resulting continuous endpoint inequalities exactly.

Otherwise a dotted contact outside `A10` costs more than 2, so the `A10` cost is less than `141/70`. The complete finite-label reduction leaves only endpoints 1 and 10 with labels 3 or 4, and second endpoints 2 and 9 with label 3. Labels at least 5 at an endpoint already exceed the remaining budget; every interior root and every parallel or dotted `A10` contact is excluded by its diagonal cost. Enumerating the resulting 36 combinations leaves exactly eight profiles: the zero profile, either endpoint at label 3 or 4, either second endpoint at label 3, or both endpoints at label 3.

Put `K0=9999/10000`. Reverse Lorentz Cauchy–Schwarz bounds the frame contribution by `-sqrt(K0*u)`. The other support blocks contribute zero to `r.E0`. For each nonzero `A10` profile of cost `q` and payment `z_A`, the checker proves

```
u > q,       K0*q > z_A^2.
```

For the zero profile, `u>0` suffices. These eight exact inequalities, plus the two continuous endpoint cases, prove `r.E0<0` for every cap in the stated class.

Any finite mutually coherent collection of these caps therefore preserves a common spacelike escape. A sufficiently small perturbation of `E0` toward a shared strict interior point makes every wall pairing strictly negative while preserving its positive norm. The ray from an interior point through this strict spacelike direction reaches the future light cone with every wall inequality strict, giving an open ideal patch and infinite volume.

The consequence for construction is specific: **a wall that cuts `E0` must have at most two nondotted chain contacts, hence at least 24 dotted ones.** A cap with at least three nondotted chain contacts can still cut other escapes or add finite corners. The previously proved [universal lower bound of five dotted chain contacts](fresh_h20_chain_abelian_gate_20261016_review.md) for an arbitrary cap is unchanged. The result supplies neither a dimension bound nor a new finite-volume polytope, and leaves open caps with at most two nondotted chain contacts or a construction that changes the retained chamber.
