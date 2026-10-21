# Exact full-`H^2` extension gate at the golden rank-twenty boundary

**Status (2026-08-09): theorem and exact integral certificates.**  For a
fixed primitive embedding of a golden divisor lattice in the K3 lattice,
the passage from an `NS`-only `Q(sqrt(5))` action to a full-cohomology Hodge
action splits into three independent gates:

```text
rational self-adjoint transfer  ->  integral discriminant glue  ->  Hodge period
```

The first two gates are lattice-theoretic and are necessary and sufficient.
The third is a real-signature condition plus the usual generic-period
condition.  At Picard rank twenty the third gate always fails, even when the
first two pass.

The discriminant-11 and discriminant-19 golden packages make this separation
literal.  Each has four integral self-adjoint square-root actions on its
positive binary complement; exactly two preserve the unimodular K3 glue as
actions of `Z[sqrt(5)]`.  No such extension preserves the maximal order
`Z[phi]`.  Nevertheless **none** is a Hodge action, because either real
eigenspace on the transcendental plane has dimension one.

This is a `P`-computable obstruction to the proposed automatic
`C_35^conn ->` coherent K3 carrier lemma.  It is not a compact-polytope
nonexistence theorem: the two rank-twenty packages used below have not passed
the complete-root/same-chamber test and are not asserted to be compact
polytopes.

## 1. Fixed-embedding extension theorem

Let `Lambda=Lambda_K3`, let `N -> Lambda` be a primitive nondegenerate
embedding, and put `T=N^perp`.  Since `Lambda` is unimodular, the inclusion

\[
                 N\mathbin\perp T\subset\Lambda
\]

is encoded by an anti-isometry of discriminant forms

\[
             \gamma:(A_N,q_N)\longrightarrow(A_T,-q_T),             \tag{1.1}
\]

and `Lambda/(N+T)` is the graph of `gamma`.  Suppose `F=Q(sqrt(5))` acts
self-adjointly on `N_Q`; write `tau_N` for multiplication by
`phi=(1+sqrt(5))/2` and `B_N=2 tau_N-I`.

### Theorem 1.1 -- rational, integral, and Hodge criteria

For this fixed embedding the following statements are exact.

1. A unital self-adjoint `F`-action on `Lambda_Q` extending the action on
   `N_Q` exists if and only if there is a rational endomorphism `B_T` of
   `T_Q` satisfying

   \[
       B_T^2=5I,\qquad B_T^tG_T=G_TB_T.                              \tag{1.2}
   \]

   The extension is `B_N+ B_T`, and `tau=(I+B)/2`.

2. Assume `B_N` and `B_T` preserve `N` and `T`.  The suborder
   `Z[sqrt(5)]` acts on the integral unimodular lattice `Lambda` if and only
   if

   \[
       \boxed{\gamma\,\bar B_N=\bar B_T\,\gamma}                    \tag{1.3}
   \]

   on the discriminant groups.  The maximal order `Z[phi]` acts if and only
   if the analogous identity holds for
   `bar tau_N,bar tau_T`; in particular `tau_T=(I+B_T)/2` must preserve `T`.

3. Suppose the rational extension exists.  A marked K3 period with exact
   Picard lattice `N` and with this full `F`-action by Hodge endomorphisms
   exists if and only if, for some real embedding `sigma_0:F->R`,

   \[
      \operatorname{sig}(T_{\sigma_0})=(2,m-2),\qquad
      \operatorname{sig}(T_\sigma)=(0,m)\quad(\sigma\ne\sigma_0),    \tag{1.4}
   \]

   where `m=dim_F T_Q`, together with avoidance of the rational
   Noether--Lefschetz hyperplanes.  The avoidance is automatic for `m=2`
   and is achieved by a generic plane for `m>2`.

#### Proof

Because the action on `N_Q` is self-adjoint, any self-adjoint extension
preserves `T_Q`: for `t in T_Q` and `n in N_Q`,
`(tau t,n)=(t,tau n)=0`.  Thus every rational extension is the orthogonal
sum of the fixed action on `N_Q` and an action on `T_Q`, proving (1.2).

An integral self-adjoint endomorphism of a lattice preserves its dual and
hence acts on its discriminant group.  The orthogonal-sum endomorphism
preserves the overlattice `Lambda` exactly when it preserves the graph of
`gamma`; this is precisely (1.3).  The same argument for every generator of
an order gives the maximal-order assertion.  This proves both necessity and
sufficiency of the integral test, not merely a determinant obstruction.

Over the reals, self-adjointness makes the eigenspaces belonging to the real
embeddings of `F` orthogonal.  A Hodge action acts on the complex line
`H20` by one real embedding.  Complex conjugation forces the same embedding
on `H02`, so the entire oriented positive period plane lies in one
eigenspace.  Since `T` has exactly two positive directions, this gives
(1.4), including `m>=2`.  Conversely, an oriented positive plane in the
distinguished eigenspace defines a K3 period and is `F`-stable.  Avoiding
the countable union of rational orthogonality loci makes `NS=N` exactly.
For `m=2`, no nonzero rational vector can lie in the sum of the other
irrational eigenspaces, so the unique positive eigenspace already avoids
such a vector.  The marked K3 period theorem completes the converse.
`square`

### Transfer/Hasse formulation of the rational gate

Fix any nonzero `Q`-linear functional `s:F->Q`, for example
`s(x)=Tr_F/Q(alpha*x)` with the trace scalar used by the polytope compiler.
The equations (1.2) are equivalent to the existence of an `F`-quadratic
space `(W,h_T)` for which

\[
                    (T_Q,q_T)\simeq s_*(W,h_T).                      \tag{1.5}
\]

Indeed self-adjointness makes `q_T(ax,y)` balanced in the two arguments;
nondegeneracy of the field-trace pairing recovers the unique `F`-valued
form `h_T`, and transfer gives the converse.

For a proposed `h_T`, (1.5) is decided exactly by Hasse--Minkowski: compare
dimension, determinant square class, real signatures, and local Hasse
invariants at `2` and the finite primes occurring in exact
diagonalizations.  Equivalently one may solve the rational matrix equations
(1.2) directly.  This is the precise rational gate; primitive embeddability
of `N` alone does not imply it.

## 2. Universal rank-twenty obstruction

For the full trace of a quadratic `H^9` realization,

\[
       \operatorname{rk}N=2(9+1)=20,\qquad
       \operatorname{rk}T=2,\qquad \dim_FT_Q=1.                     \tag{2.1}
\]

Every `F` eigenspace of `T_R` is therefore one-dimensional.  The real K3
period plane is two-dimensional.  Theorem 1.1(3) immediately gives

\[
 \boxed{\text{no quadratic full-trace rank-20 K3 package is Hodge coherent}.}
                                                                         \tag{2.2}
\]

This obstruction is independent of the discriminant form, the primitive
embedding orbit, the chosen order, and root completeness.  It is the exact
boundary behind the coherent K3 inequality

\[
                         2(n+1)+4\le22,\qquad n\le8.                 \tag{2.3}
\]

It does **not** exclude an `NS`-only rank-twenty K3 carrier.  In that weaker
category, the field need not act on `T` at all and the K3 ceiling remains
`n<=9`.

## 3. Exact discriminant-11 and discriminant-19 certificates

The two existing golden divisor packages have

\[
\begin{array}{c|c|c|c}
\text{package}&T&|\det T|&[\Lambda:N\perp T]\\ \hline
D11&\begin{pmatrix}18&5\\5&2\end{pmatrix}&11&11\\[3pt]
D19&\begin{pmatrix}2&1\\1&10\end{pmatrix}&19&19.
\end{array}                                                           \tag{3.1}
\]

For a positive binary Gram matrix

\[
                  G_T=\begin{pmatrix}x&y\\y&z\end{pmatrix},
\]

every integral nonscalar solution of (1.2) has the form

\[
       B_T=\begin{pmatrix}a&b\\c&-a\end{pmatrix},\qquad
       xa^2+2yac+zc^2=5x,\qquad xb=2ya+zc.                           \tag{3.2}
\]

Positive definiteness gives the exhaustive finite bounds

\[
       \det(T)a^2\le5xz,\qquad \det(T)c^2\le5x^2.                  \tag{3.3}
\]

The certifier enumerates exactly within (3.3), conjugates
`B_N+B_T` through the literal index-11 or index-19 gluing basis, and tests
all matrix entries.  The results are:

| package | all integral `B_T` | preserve full `Z[sqrt(5)]` glue | preserve full `Z[phi]` glue | Hodge coherent |
|---|---:|---:|---:|---:|
| `D11` | 4 | 2 | 0 | 0 |
| `D19` | 4 | 2 | 0 | 0 |

The complete complement lists are

\[
\begin{aligned}
D11:\;&
\begin{pmatrix}-4&-1\\11&4\end{pmatrix},
\begin{pmatrix}-1&-1\\-4&1\end{pmatrix},
\begin{pmatrix}1&1\\4&-1\end{pmatrix},
\begin{pmatrix}4&1\\-11&-4\end{pmatrix};\\[4pt]
D19:\;&
\begin{pmatrix}-1&4\\1&1\end{pmatrix},
\begin{pmatrix}0&-5\\-1&0\end{pmatrix},
\begin{pmatrix}0&5\\1&0\end{pmatrix},
\begin{pmatrix}1&-4\\-1&-1\end{pmatrix}.
\end{aligned}                                                        \tag{3.4}
\]

For `D11`, the first and third matrices preserve the full K3 lattice.  For
`D19`, the first and second do.  Direct `22 x 22` multiplication verifies

\[
                      B^2=5I,\qquad B^tG_{K3}=G_{K3}B.               \tag{3.5}
\]

Thus lattice glue is demonstrably **not** the reason coherence fails in
these examples.  The failure occurs only at the Hodge-period gate.  This is
a concrete counterexample to the tempting implication

```text
primitive K3 embedding + NS field action
    => automatic coherent full-H2 field action.
```

At the rational level the implication already needs a transfer condition;
at the integral level it needs (1.3); and even when both are met, (1.4) is
independent and fatal at rank twenty.

## 4. Consequence for the compact `C_35^conn` strategy

The proposed compact theorem cannot treat full-`H^2` coherence as a routine
completion of the finite saturation/embedding compiler.  The exact logic is
now:

1. Starting from `P in C_35^conn`, normalize its unique equal-root golden
   scaling and enumerate the finite trace-overclosure and primitive-embedding
   packages.
2. For each package, compute `T` and its discriminant anti-isometry `gamma`.
3. Apply (1.2), then (1.3), then the signature test (1.4).
4. Independently apply the complete restricted-root/no-hidden-cut test to
   return the same `P`.
5. Accept a coherent K3 carrier only if one branch passes **all** gates.

At `n=9`, step 3 rejects every branch before root enumeration.  Therefore a
proof that every compact member of `C_35^conn` has such a coherent carrier
would already prove that the class has no compact `H^9` member.  The missing
essential-surjectivity lemma contains the desired boundary theorem; it
cannot be justified merely from primitive embeddings, discriminant glue, or
the existence of `NS`-only surfaces.

This yields a sharp fork for further work:

- **theorem track:** derive coherence from an additional natural,
  surface-free condition on `P`; then (2.2) proves the sharp `n<=8` bound,
  with Bugaenko as endpoint;
- **example track:** search rank-twenty `NS`-only packages and complete-root
  returns without imposing full-`H^2` coherence.  Any compact `H^9` survivor
  would prove that the proposed coherent universality statement is false,
  while still giving the desired new positive example.

The two tracks must remain distinct: imposing coherence in the example
search deletes the entire `H^9` search space by theorem.

## 5. Reproduction

- Exact certifier:
  `scripts/certify_golden_rank20_full_h2_extension_gate.py`.
- Regression tests:
  `tests/test_golden_rank20_full_h2_extension_gate.py`.
- Input divisor packages:
  `scripts/coupled_discriminant11_candidate.py` and
  `scripts/coupled_discriminant19_candidate.py`.

Run

```bash
PYTHONPATH=src:. python3 scripts/certify_golden_rank20_full_h2_extension_gate.py
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_golden_rank20_full_h2_extension_gate.py
```

The test suite checks the exhaustive four-element lists, both gluing
indices, the two compatible suborder extensions in each package, the absence
of maximal-order extensions, and the universal period-plane failure.
