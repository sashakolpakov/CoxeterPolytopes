# Golden full-rank path host in H22: exact Vinberg prefix through height 10

## Verdict

The independently audited one-dotted label-three path does define a genuine
rank-23 arithmetic lattice over `Q(sqrt(5))`.  Its full arithmetic orthogonal
host is cocompact.  An exact Vinberg enumeration through inclusive height
`H=10` finds one wall beyond the original path seed.

This is **not** a new compact Coxeter polytope.  The resulting 24-wall prefix
cone has `2` finite vertices, `0` cusps, and `41` hyperideal extreme rays.
Later walls can cut those rays, so the host's reflectivity remains **OPEN**.

The frozen certificate is
[`golden_path_h22_h10_prefix.json`](../data/enumerations/golden_path_h22_h10_prefix.json),
SHA-256
`bb968be0a53bbb1e274372e150fe3bb27af472d59c255c503e0194d13b2ea80e`.

## Exact host

Put `c=2+2phi=2phi^2` and let `A` be the rank-23 Cartan matrix with

```text
A_ii=2,
A_01=-1,       A_10=-c,
A_i,i+1=A_i+1,i=-1       (1<=i<22),
```

all other entries being zero.  The symmetrizer is
`D=diag(c,1,...,1)`, and the integral symmetric lattice Gram is `B=DA`:

```text
diag(B)=(2c,2,...,2),
B_01=B_10=-c,
B_i,i+1=B_i+1,i=-1       (1<=i<22).
```

Exact congruence elimination gives

```text
det(B) = -84-172phi = -4(2+phi)(4+13phi),
inertia_id(B)  = (22 positive, 1 negative, 0 zero),
inertia_conj(B)= (23 positive, 0 negative, 0 zero).
```

The split-prime normalization in the certificate is worth making explicit:

```text
phi(4+13phi)=13+17phi=phi^2(9+4phi).
```

Thus the canonical totally-positive generator used for the prime of norm
`101` is `9+4phi`, not `13+17phi`.  The other odd generator is
`2+phi`, of norm `5`.

The Cartan two-cycle has product `A_01 A_10=c`, which is nonrational.
Consequently the minimal Vinberg field is exactly `Q(sqrt(5))`.  Positive
definiteness at the conjugate real place also proves anisotropy over the
field: a field-rational isotropic vector would conjugate to a nonzero
isotropic vector of a positive-definite form.  Hence the full arithmetic
orthogonal host is cocompact.  This says nothing yet about the index of its
reflection subgroup.

## Height and complete local norm candidates

The exact vector orthogonal to the spatial `A22` block is

```text
rho=e0+(c/23)(22e1+21e2+...+e22),
rho^2=(-84-172phi)/23.
```

It obeys `B rho=rho^2 e0`; it is **not** itself `B^-1 e0`, which is instead
`rho/rho^2`.  For a root `r=(x,y)` of norm `q`, the Vinberg height used here
is

```text
H(r)=-rho^2 x^2/q.
```

For a primitive reflective vector, write `I(r)=(B(r,L))`.  At every prime
`P`, the adjugate identity and primitivity give
`v_P(I(r))<=v_P(det B)`, while reflection integrality gives `q|2I(r)`.
The only possible norm support is therefore

```text
(2)^a (2+phi)^b (9+4phi)^d,
0<=a<=3, 0<=b<=1, 0<=d<=1.
```

After totally-positive unit-square normalization, the complete 16-element
candidate alphabet is

```text
1, 2, 2+phi, 4, 4+2phi, 8, 8+4phi, 9+4phi,
23-phi, 16+8phi, 18+8phi, 46-2phi,
36+16phi, 92-4phi, 72+32phi, 184-8phi.
```

This is a complete candidate alphabet; some classes may be empty.

## Exact H<=10 enumeration

For every norm `q` and control `x`, put

```text
T = 4(q-rho^2 x^2)/q^2,
z_i = -2 B(r,e_i)/q       (1<=i<=22).
```

Acuteness against the spatial simple roots makes every `z_i` nonnegative at
the selected place.  This loses no root: the finite `A22` Weyl group fixes
`rho` and moves every positive-height mirror to the chosen spatial chamber.
Reflection integrality makes each `z_i` a golden integer, and
minimizing over all other coordinates gives the exact two-place bound
`z_i^2<=2T`.  The recursion uses the inverse `A_d` energy as the exact Schur
lower bound at depth `d`, then reconstructs the spatial vector with
`A22^-1` and tests the remaining `e0` reflection coefficient.

The complete replay has:

```text
16 norm candidates,
383 two-place controls through H<=10,
19,760 Schur-recursion nodes,
6 exact A22-inverse energy solutions,
3 integral reflective candidates.
```

The frozen stream hashes are

```text
controls:   1db06f9e673ea0f11e59ad220e9ca65b0408d98eb1c77e2f4701c5586bc4b1f6
candidates: f758c1b6a105df9b7817f0c58147d7270b4e467130ce9521c866a78923b02714
```

The chronological decisions are:

1. At `H=(22phi-1)/23`, the norm-four root `(phi-1)e0` is accepted.  It is a
   unit multiple of `e0`, hence the original terminal dotted-path
   hyperplane—not a new wall.
2. At `H=(42+86phi)/23`, two norm-two roots occur simultaneously.
   The root `e0+e1` is rejected because its pairing with `(phi-1)e0` is
   `+2phi`.
3. The other norm-two root is accepted:

   ```text
   r=(1; 1+phi, 0, -phi, -phi, ..., -phi).
   ```

   Its only nonzero pairings with previous walls are

   ```text
   B(r,(phi-1)e0)=-2,
   B(r,e2)=-1,
   B(r,e3)=B(r,e22)=-phi.
   ```

   Thus it attaches by Coxeter labels `4,3,5,5`, respectively.  It is the
   one genuinely new exact Vinberg wall in this bounded host lane.

## Exact cone obstruction at the current prefix

Before the new wall, the true 23-wall seed cone has

```text
2 finite vertices, 0 cusps, 21 hyperideal rays, 0 opposite-sheet rays.
```

After adding the new wall, exact incremental double description gives

```text
2 finite vertices, 0 cusps, 41 hyperideal rays, 0 opposite-sheet rays.
```

One independently recomputed hyperideal extreme ray has `v_0=1` and, for
`1<=i<=22`,

```text
v_i=((1142-138i)+(1142-126i)phi)/571.
```

Its exact norm is

```text
(5708212+7337356phi)/326041 > 0.
```

It is active on facets `0,...,20,23`; its only nonzero wall pairings are
those with facet `21` and the seed facet `22`, both equal to

```text
-(2032+1756phi)/571.
```

This proves that the **complete bounded prefix** is not finite-volume.  It
does not prove that the completed Vinberg chamber is infinite-volume, since
walls above height 10 may remove the ray.

## Engine and scope

The repaired pinned AlVin lineage was inspected but not used for this host.
Its command interface enumerates standard coordinate lattices of diagonal
forms.  The present lattice is non-diagonal, and no `GL_23(Z[phi])`
diagonalization was constructed.  Feeding a mere field-level diagonalization
would silently change the integral lattice.  The bounded calculation instead
uses exact native two-place boxes, Schur pruning, lattice reconstruction, and
double description.

Reproduce with:

```sh
PYTHONPATH=src python3 scripts/certify_golden_path_h22_h10_prefix.py --verify
pytest -q tests/test_golden_path_h22_h10_prefix.py
```

The result adds one real wall and one reusable cocompact golden host.  It
does not add a compact example, establish reflectivity, or prove golden
nonexistence in `H22--H29`.
