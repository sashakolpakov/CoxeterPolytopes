# A common escape for every symmetric single-support short-arc cap

The [new134-wall cap](fresh_h20_dotted_arc_cap_20261016.md) supplies useful
finite geometry and cuts the previous847-bank escape. However, its new
simple escape survives **every** cap in the following larger family, and
therefore any mutually compatible finite collection of them:

- retain all133 old walls;
- be right to17 of the18 S61 support roots;
- have finite or parallel contacts with three consecutive chain walls,
  with equal contacts on the two outer walls.

All24 chain centers, all18 choices of released support root and all finite
labels are covered. Other old-wall contacts may be dotted. This does not
cover unequal outer contacts or two nonright support contacts.

The [standalone rational checker](../scripts/fresh_h20_dotted_arc_family_gate_20261016.py)
and [result](fresh_h20_dotted_arc_family_gate_20261016.json) pass generation
and fresh saved replay. The geometric argument below uses the frozen
[nonabelian-chain theorem](fresh_h20_chain_abelian_gate_20261016.md) and its
[independent review](fresh_h20_chain_abelian_gate_20261016_review.md).
Independent review of this additional family lemma is pending at publication.

## Reduction without a label cutoff

Normalize the new wall to norm2. Let its side contact magnitude be B and its
middle contact magnitude be A, where each belongs to
`{0} union {2cos(pi/m): integer m>=3} union {2}`.
Let C=a²/(a²+2) be the saved rotation cosine. Write the cap
`r=y-c*omega_j`, where c≥0, Jjj=omega_j², and y lies in the perpendicular
Lorentz3 frame. The three chain normals are independent: their Gram
determinant is `-32*C*(C+1)`, which is nonzero.

Writing the central chain normal as `(a,b,z)`, its unique solution is

```
y = (((B-A)*a²+2B)/(2a), (B-A)*b/2, (B-A)*z/2),
y² = A²/2-B²/(2C),
Jjj*c² = 2-y².
```

If B≤A, the spatial coefficient is nonpositive. Every chain spatial
normal has the same length, so its pairing with y is at least the central
pairing−A, hence at least−2. Coherence would then make all26 contacts
finite or parallel, contrary to the reviewed nonabelian-chain obstruction.
Thus B>A and B≥1. In particular y is future timelike.

All finite/parallel contact magnitudes lie in Q^ab. If c were also finite
or parallel, the norm equation would give

```
C = B²/(A²+2*Jjj*c²-4) in Q^ab,
```

contrary to the same theorem. Its denominator cannot vanish because B>0.
Therefore the released support contact must be dotted: c>2. This forces

```
Jjj < (1+1/C)/2 < 281/280.
```

The literal rational support inverse leaves exactly the eight roots
`1,10,12,15,17,31,61,87`. This is a complete consequence of the algebraic
field obstruction, not a bounded numerical support search.

## Exact common-escape inequality

Use the already certified strict spacelike vector

```
E = x-e/100-omega5-(1/100000)*sum(other support dual weights).
```

Let `K=9999/10000`, let `z_j=-omega_j.E>0`, and let `u=-y²>0`.
The support dual coordinates are rational, so every z_j is rational.
Reverse Lorentz Cauchy–Schwarz gives

```
r.E <= -sqrt(K*u) + z_j*sqrt((2+u)/Jjj).
```

The two exact lower bounds are

```
u > 4*Jjj-2,             u > 1/a².
```

The first follows from c>2; the second from B>A and B≥1. With `a_hi`
the frozen rational upper endpoint, put
`u0=max(4*Jjj-2,1/a_hi²)`. For each of the eight possible roots the
checker verifies, using Fraction arithmetic,

```
Jjj*K-z_j² > 0,
u0*(Jjj*K-z_j²)-2*z_j² > 0.
```

Consequently the negative term strictly dominates for every u≥u0,
and r.E<0. The smallest rational margin is greater than0.0017784, at
the isolated A1 roots17/31. No limiting parallel label or unbounded finite
label is omitted.

This supplies a precise next construction requirement: to cut this escape
while retaining the old walls, a three-contact finite-arc cap must leave the
symmetric single-support family. Unequal outer contacts or at least two
nonright support contacts are concrete remaining possibilities.
