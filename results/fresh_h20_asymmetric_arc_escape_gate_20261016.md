# Every consecutive finite chain arc preserves a common H20 escape

Every coherent added Coxeter wall retaining the frozen133-wall chamber,
whose contacts with **any three consecutive chain walls are finite or
parallel**, pairs strictly negatively with the same spacelike vector E0.
The three contacts may be unequal, and all18 support contacts are arbitrary
Coxeter contacts. Consequently any finite coherent collection of these
walls still leaves infinite volume. The already constructed134th wall is
included in this family.

This is a new proof candidate, with primary generation and fresh saved
replay passing; separate independent review is pending. It strengthens the
[reviewed symmetric-arc result](fresh_h20_symmetric_arc_escape_gate_20261016_review.md).
It does not exclude individual walls, nor caps without a consecutive
nondotted triple, nor constructions changing the retained walls.

The [standalone checker](../scripts/fresh_h20_asymmetric_arc_escape_gate_20261016.py)
and [certificate](fresh_h20_asymmetric_arc_escape_gate_20261016.json) use
Fraction arithmetic and Q(√2). The physical chamber and modular
factorizations are inherited from the separately reviewed
[134-wall cap certificate](fresh_h20_dotted_arc_cap_20261016_review.md) and
[chain field certificate](fresh_h20_chain_abelian_gate_20261016_review.md).
No numerical proposal bank or finite angle cutoff is a dependency.

## 1. The rotation retains degree50 after adjoining all abelian numbers

Let K be the literal degree64 tower in the chain certificate, let
C=a²/(a²+2), and let P be its certified degree50 polynomial over K.
The inherited exact calculations prove that P is irreducible over K and
has squarefree factor degrees13 and37 at a degree-one prime over359.
The leading coefficient and all denominators are local units; the tower
embeddings are certified by simple Hensel roots. Thus the Galois group G
of the splitting field L/K is transitive on50 roots and contains an
actual Frobenius permutation of cycle type(13)(37).

Here the needed local fact is stronger than merely unequal residue
degrees: squarefree reduction makes the roots distinct in the unramified
local splitting field, and Frobenius permutes them with the degrees of
its irreducible residue factors. See [Sutherland, Number Theory I,
Lecture7, Proposition7.17](https://math.mit.edu/classes/18.785/2025/LectureNotes7.pdf)
for the Frobenius action, together with the local factorization theorem
in [Lecture6, Dedekind–Kummer](https://math.mit.edu/classes/18.785/2025/LectureNotes6.pdf).
The checker binds the earlier complete modular verification rather than
repeating it.

Raise this Frobenius element to the13th power. It becomes a37-cycle
fixing13 roots. If G preserved a nontrivial block system, its block size
would be2,5,10 or25, with respectively25,10,5 or2 blocks. The37-cycle
cannot nontrivially permute fewer than37 blocks, so it fixes every block
setwise. Its37-point orbit would then have to lie inside a block of size
less than37, a contradiction. Therefore G is primitive.

A faithful transitive abelian permutation group is regular: a point
stabilizer would fix every point and hence be trivial. Our37-cycle has
fixed points, so G is nonabelian. Its commutator G′ is a nontrivial normal
subgroup. The orbits of a normal subgroup form a G-invariant block
system; primitivity and faithfulness therefore make G′ transitive.

Now K Q^ab/K is a possibly infinite Galois abelian extension. Its finite
intersection E=L∩K Q^ab is Galois abelian over K, so Gal(L/E) contains G′.
The subgroup Gal(L/E) is still transitive on the50 roots. Restriction
identifies it with Gal(L K Q^ab/K Q^ab), hence

```
[K Q^ab(C) : K Q^ab] = 50.
```

In particular, C satisfies **no nonzero polynomial of degree below50
over Q^ab**. We do not assert that its degree over Q^ab itself equals50.
This argument uses no classification of primitive groups or Jordan theorem.

## 2. Any consecutive triple forces a dotted support contact

Normalize a proposed spacelike wall r to norm2. Write

```
r = y - sum_i(d_i*omega_i),       d_i>=0,
Q = d^T J d,                    y²=2-Q,
```

where the omega_i are dual to the18 old support roots, J is their
positive-definite rational inverse Gram matrix, and y lies in the
orthogonal Lorentz3 frame. J is entrywise nonnegative and its minimum
diagonal is1/2.

Let the left, middle and right chain contacts be −B,−A,−D. Each magnitude
is0, a finite value2cos(π/m), or the parallel value2, and therefore lies
in Q^ab. The consecutive chain Gram is

```
G3 = [[2,0,-4C-2],[0,2,0],[-4C-2,0,2]],
det(G3) = -32C(C+1).
```

The triple spans the Lorentz3 frame. Its inverse gives

```
y² = A²/2 - ((B+D)²+4CBD)/(8C(C+1)).
```

First dispose of B=D=0. If A>0, the next available chain contact on one
side is (2C+1)A>0 by the exact rotation recurrence, contradicting
acuteness. Every one of the24 triple positions has such a fourth chain
wall. If A=0 too, then y=0. A nonzero r=−sum(d_i*omega_i) pairs positively
with every strict point of the old support sector: each omega_i pairs
negatively there since J≥0 with positive diagonal. Hence it cannot define
a coherent added wall with strict interior. We may assume B+D>0.

If all support contacts were finite or parallel, all d_i and Q would lie
in Q^ab. The norm equation would then give

```
8C(C+1)*(A²/2-2+Q) - (B+D)² - 4CBD = 0.
```

This is a polynomial in C of degree at most2 whose constant coefficient
is −(B+D)²<0. The degree50 result forbids it. Thus some support contact
is dotted, with d_i>2. Because J≥0 and Jii≥1/2, Q>2 and

```
u=-y²=Q-2>0.
```

The sum of the two outer chain roots is future timelike: it has time2a>0
and square−8C<0. Its pairing with y is−(B+D)<0. Thus y is on the future
timelike sheet. Finally, 0≤A,B,D≤2 gives the uniform upper bound

```
u <= 2/C = 2+4/a²,
Q=2+u <=4+4/a² <281/70.
```

No positive numerical lower bound for u is required.

## 3. A weak escape removes the small-u difficulty

Use the same orthonormal frame x,e,f as the saved chamber and set

```
E0 = x-e/100-omega5,
K = 9999/10000,
E0² = 30/11-K = 190011/110000 >0.
```

The independently certified strict escape in the134-wall package was

```
Eeps = E0 - eps*sum_{j!=5}(omega_j),  eps=1/100000.
```

For any old wall N outside the support, acuteness gives N.n_i≤0.
Since omega_j=sum_i Jij*n_i and Jij≥0, omega_j.N≤0. Therefore
E0.N≤Eeps.N<0 off the support. On the support, E0.n_i=−δ_i5.
The same argument gives E0.H<0 for the saved strict center H. Thus E0
is an exact weakly feasible spacelike vector for all134 retained walls.

For a proposed cap, put z_i=−omega_i.E0=Ji5. Every z_i outside A10
vanishes. Reverse Lorentz Cauchy–Schwarz gives

```
r.E0 <= -sqrt(K*u) + sum_{i in A10}(d_i*z_i).
```

This is why no infinite-label gap computation is needed: if the A10
profile is zero, the right side is already strictly negative for every
u>0, however small.

## 4. The complete A10 budget reduction

The A10 block consists of roots1,…,10 and has
`Jij=min(i,j)*(11-max(i,j))/11`. If a dotted contact belongs to this block,
its diagonal cost excludes every node except endpoint1 or10. Such an
endpoint costs more than40/11, leaving less than1/2 in the total budget.
Every other nonzero contact costs at least1/2, so all other support
contacts vanish. For each endpoint, with J=Jjj and z=Jj5, the checker
verifies

```
J*K-z²>0,
(4J-2)*(J*K-z²)-2z²>0.
```

Since u>4J−2 and d_j²=(2+u)/J, these two rational inequalities imply
r.E0<0 along the entire dotted endpoint ray.

Otherwise a dotted contact outside A10 costs more than2, leaving A10
cost q_A<141/70. Roots3,…,8 are excluded by diagonal at least24/11.
Endpoints1/10 permit only labels3/4: label5 and larger already cost
more than `(10/11)*(5/2)>141/70`. Second endpoints2/9 permit only label3.
The36 exact combinations in these four domains leave precisely:

| A10 profile | q_A | z_A=sum(d_i*Ji5) |
|---|---:|---:|
| zero |0|0|
| root1, label3 |10/11|6/11|
| root10, label3 |10/11|5/11|
| root1, label4 |20/11|6√2/11|
| root10, label4 |20/11|5√2/11|
| root2, label3 |18/11|12/11|
| root9, label3 |18/11|10/11|
| roots1 and10, both label3 |2|1|

For every nonzero row, the exact margin K*q_A−z_A² is positive.
The dotted outside contact makes u>q_A, so
`r.E0 <= -sqrt(K*u)+z_A<0`. The zero row was already covered by u>0.
Together with the two dotted endpoint rays, this exhausts all support
patterns and all angle labels.

## 5. Infinite-volume consequence and scope

Take any finite coherent collection of these caps. Each pairs strictly
negatively with E0, while every retained old wall pairs nonpositively.
Let H be a strict interior future timelike point of the common chamber.
For sufficiently small δ>0, E0+δH is still spacelike and has strictly
negative pairing with every wall. To fix the time orientation, take the
saved old future center H0. The certificate gives E0.H0<0, and every
future timelike H has H.H0<0. Thus E0+δH and the whole segment from H
toward it remain in the future affine chart defined by negative pairing
with H0. That segment crosses the future null cone through a point where
all wall inequalities remain strict. No assertion E0.H<0 for an arbitrary
new interior point is needed. A neighborhood of that ideal point is
also allowed, proving infinite hyperbolic volume.

This works with different consecutive triples for different caps. It does
not assume the proposed walls are arithmetic, does not restrict their
coefficient fields, and allows all finite labels as well as parallel
contacts. It leaves open caps with no three-consecutive nondotted arc,
other chain contact patterns, and changed retained-wall constructions.
The eight-profile argument is a complete reduction by an exact norm
budget; the modular factorization is inherited from the independently
reviewed degree50 certificate. No new ambient cone census is claimed.

Reproduce with:

```sh
python scripts/fresh_h20_asymmetric_arc_escape_gate_20261016.py --verify
```
