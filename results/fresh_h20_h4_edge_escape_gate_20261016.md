# The four orthogonal edge axes cannot cut E4

Every new acute Coxeter wall retaining the
[145-wall chamber](fresh_h20_residual_h4_corner_20261016.md), right to
S16 and to any three of its four mutually right residual caps, has
nonpositive pairing with the strict escape E4. This is an all-parameter
obstruction, with no norm, height or Coxeter-label cutoff. A next E4
cutter must bend at least one of those three contacts or leave S16.

The [source](../scripts/fresh_h20_h4_edge_escape_gate_20261016.py) and
[certificate](fresh_h20_h4_edge_escape_gate_20261016.json) bind the exact
rational residual Gram and corner, and the certified physical pairing
intervals. Generation and fresh saved replay pass. This is a primary
certificate, not an independent replay of the145-wall geometry.

Let n_i be the omitted residual cap and C the exact finite corner ray.
S16 and the other three residual caps have positive rank19, so their
orthogonal complement is the Lorentz plane spanned by C,n_i. A new
positive-norm normal acute to n_i therefore has the unique form

```
r=alpha*(-n_i+lambda*C), alpha>0.
```

The coefficient of n_i cannot vanish, since C is timelike. Positive
scaling is irrelevant to every angle and escape sign. Write
`C^2=-k<0`. Then

```
r^2/alpha^2=q_i-k*lambda^2,
r.E4/alpha=-n_i.E4+lambda*C.E4.
```

The exact value C.E4 is negative. Consequently cutting E4 requires
`lambda<L_i=(n_i.E4)/(C.E4)`.

For i=141,142,143, acuteness to oldcap133 gives

```
lambda >= (n_i.n133)/(C.n133) > L_i.
```

All divisions have certified negative denominators and the three strict
separations are stored. Thus these three axes are impossible even for
an arbitrary acute new halfspace.

For i=144 (w), acuteness to chainwall132 gives
`lambda>=A=(w.n132)/(C.n132)>0`, while cutting requires `lambda<U=L_144`.
The certified approximate interval is

```
7.227585... <= lambda < 7.306743....
```

The exact norm is positive throughout this interval. At U the checker
certifies a negative pairing g with132 and

```
4*g^2 < 2*(q_144-k*U^2).
```

For A<lambda<U, the magnitude -g increases with lambda and the positive
norm decreases. Hence the Cartan two-cycle
`4*g(lambda)^2/[2*(q_144-k*lambda^2)]` is strictly between0 and1. No acute
Coxeter relation has such a value: a right pair gives0, the smallest
nonright finite value is1 (label3), parallel gives4, and dotted gives
more than4.

At the remaining endpoint lambda=A, the new wall is right to132. The
checker instead uses oldroot31: its pairing is strictly negative and
its Cartan two-cycle is again strictly between0 and1. Thus that endpoint
also fails the all-label Coxeter test.

The four cases exhaust the stated orthogonal edge axes. They do not
exclude finite-angle contacts to the three retained residual caps,
arbitrary mixed-support walls, or other constructions. Orthogonal
opposite-edge caps may create additional finite edges elsewhere, but
cannot eliminate the specific strict escape E4.

```sh
python scripts/fresh_h20_h4_edge_escape_gate_20261016.py --verify
```
