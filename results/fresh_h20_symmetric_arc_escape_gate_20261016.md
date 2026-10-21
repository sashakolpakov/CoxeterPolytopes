# Every symmetric finite chain arc misses the new H20 escape

The simple strict escape from the [134-wall construction](fresh_h20_dotted_arc_cap_20261016.md)
survives every coherent additional cap with equal finite/parallel contacts
at the outer walls of a consecutive three-wall chain arc and a
finite/parallel middle contact. **All18 support contacts are arbitrary
Coxeter contacts**, including dotted ones. Thus any mutually compatible
finite collection of these symmetric-arc caps still leaves infinite volume.

This strengthens the [single-support family gate](fresh_h20_dotted_arc_family_gate_20261016.md).
It does not cover unequal outer contacts, fewer than three nondotted chain
contacts, or a changed retained chamber.

The [standalone Fraction checker](../scripts/fresh_h20_symmetric_arc_escape_gate_20261016.py)
and [result](fresh_h20_symmetric_arc_escape_gate_20261016.json) pass generation
and fresh saved replay. No numerical bank is a proof dependency. Independent
review of this broader lemma is pending at publication of this report.

## The forced dotted contact and total budget

Normalize the cap to norm2 and write `r=y-sum(d_i*omega_i)`, where d_i≥0
are its normalized support-contact magnitudes. The support Gram inverse J
is rational and entrywise nonnegative. The perpendicular Lorentz3
projection y is fixed by the three chain contacts. With side magnitude B,
middle magnitude A, and C=a²/(a²+2),

```
y = (((B-A)*a²+2B)/(2a), (B-A)*n_center_e/2,
     (B-A)*n_center_f/2),
y² = A²/2-B²/(2C).
```

The [reviewed nonabelian-chain theorem](fresh_h20_chain_abelian_gate_20261016_review.md)
proves C∉Q^ab and excludes an entirely finite/parallel chain-contact
profile. If B≤A, equal spatial lengths make every chain pairing at least
−A≥−2, contradicting that theorem. Hence B>A and B≥1, and y is future
timelike. Put u=−y². Then

```
u=(B²-A²)/2+B²/a² > 1/a²,
d^T J d = 2+u <=4+4/a² < 281/70.
```

If every support contact were finite or parallel, every d_i would lie in
Q^ab. The rational cost would then force
`C=B²/(A²+2*d^T J d-4)` into Q^ab; its denominator cannot vanish since B>0.
At least one support contact is therefore dotted, with magnitude greater
than2. There is no finite-label cutoff in this reduction.

## Complete reduction in the A10 block

The dominant support component of the saved escape is omega5 in the A10
block on roots1 through10. Its inverse is
`Jij=min(i,j)*(11-max(i,j))/11`. The other blocks are A4, A2 and two A1's;
their total inverse-entry sum is13 and their smallest inverse diagonal is1/2.

If a dotted contact lies in A10, its diagonal cost excludes every root
except endpoint1 or10. Even there the cost is greater than40/11. Since
`281/70-40/11<1/2`, no other support contact can be nonzero. The checker
directly verifies the two resulting continuous single-endpoint cases.

Otherwise a dotted contact outside A10 costs more than2. The A10 contact
profile must have cost less than141/70. Roots3 through8 have diagonal at
least24/11 and cannot occur. Endpoints1/10 allow only labels3 or4:
label5 and every larger label already cost more than
`(10/11)*(5/2)>141/70`. Second endpoints2/9 allow only label3, since
their label4 cost is36/11. Parallel and dotted contacts are also excluded
by these same bounds.

The exact36 combinations in these four finite domains leave precisely
eight profiles:

| A10 contacts | A10 cost |
|---|---:|
| none |0|
| root1 or10, label3 |10/11|
| root1 or10, label4 |20/11|
| root2 or9, label3 |18/11|
| roots1 and10, both label3 |2|

This table is a complete all-label reduction, not extrapolation from a
bounded scan.

## Uniform escape bound

The preserved escape is

```
E=x-e/100-omega5-epsilon*sum(other support dual weights),
epsilon=1/100000,                 K=1-1/10000.
```

Let z_i=−omega_i.E>0 and let z_A be the A10 contribution to
`sum(d_i*z_i)`. Reverse Lorentz Cauchy–Schwarz bounds the frame term by
`y.E_frame <= -sqrt(K*u)`. For the remaining blocks, ordinary positive
definite Cauchy–Schwarz gives

```
sum_outside(d_i*z_i)
  = epsilon*d_out^T*J_out*1
  < epsilon*sqrt((281/70)*13)
  < 1/1000.
```

If the A10 profile has positive cost q, then u>q, because the dotted
outside contact costs more than2. For the zero profile use
u>1/a²>1/a_hi². For all eight profiles the checker independently verifies

```
K*max(q,1/a_hi²) > (z_A+1/1000)².
```

The only irrational quantities are rational multiples of√2 in the two
label4 profiles; their signs are checked by exact quadratic-field norms.
Thus r.E<0 in every case.

For the two dotted A10 endpoint cases, write `J=Jjj`, `z=z_j` and
`u0=4J-2`. Their whole-ray bounds reduce to the two rational inequalities

```
J*K-z² > 0,
u0*(J*K-z²)-2z² > 0.
```

They also hold strictly. This exhausts the possible location of a dotted
support contact and completes the common-escape proof.

The cap already constructed remains useful: it cuts the former847-bank
escape and adds finite geometry. The present result explains why adding
more symmetric finite-arc caps cannot finish that chamber. An asymmetric
three-wall contact profile is a concrete next target.
