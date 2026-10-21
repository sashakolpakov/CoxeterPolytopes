# Two genuine K3 automorphism realizations from the `(20,2,1)` chamber

## Positive theorem

Let $N$ be Alexeev's rank-20 even 2-elementary K3 lattice with invariants
$(20,2,1)$, and let `C01,C02` be the two involution classes whose fixed
sections have dimension above ten.  Both lattice actions are induced by
actual automorphisms of the same marked singular projective K3 surface:

| Picard action | K3 automorphism | exact fixed $NS$ | Coxeter section |
|---|---|---:|---|
| `C01` | non-symplectic order four | archived rank 13 lattice | finite-volume noncompact `H12/19` |
| `C02` | non-symplectic involution | archived rank 12 lattice | finite-volume noncompact `H11/20` |

Thus the two chambers are no longer merely fixed sections of an abstract
Picard-lattice symmetry.  They are automorphism-fixed divisor/reflection-
chamber sections on a K3 surface.  The abstract chambers, their arithmetic
classes, and their dimensions are not claimed new.

## The common singular K3 period

The ambient lattice has signature $(1,19)$, determinant $-4$, and
discriminant group $(\mathbb Z/2)^2$.  Its positive K3-lattice complement is

\[
T=\langle2\rangle\oplus\langle2\rangle.
\]

In the archived half-vector basis of $A_N$, the graph anti-isometry to
$A_T$ is

\[
P=\begin{pmatrix}0&1\\1&1\end{pmatrix}.
\]

Because $\delta=1$, the certificate retains the full discriminant value by
recording $4q\pmod 8$.  For the four source classes the pairs
$(4q_N,4q_T)$ are

```text
(0,0), (4,4), (6,2), (6,2),
```

so $q_N+q_T=0$ in $\mathbb Q/2\mathbb Z$ class by class.  Gluing the
graph of $P$ gives an even unimodular lattice of signature $(3,19)$, hence
the K3 lattice.

Use the common period

\[
\omega=(1,i)\in T\otimes\mathbb C.
\]

It satisfies $\omega^2=0$, $\omega\bar\omega=4$, and no nonzero rational
vector of $T$ is perpendicular to it: orthogonality of $(a,b)$ is
$2a+2ib=0$.  Therefore the resulting Picard lattice is exactly $N$, of
rank 20, rather than a larger specialization.

## Why the two glues are different

For `C01`, the induced discriminant action is

\[
\begin{pmatrix}1&1\\0&1\end{pmatrix}.
\]

It does not glue to $-I_T$.  It does glue equivariantly to the order-four
rotation

\[
B_4=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad B_4^2=-I,
\]

under $P$.  Since $\omega B_4=i\omega$, the glued action is a
non-symplectic order-four Hodge isometry with characteristic polynomial

\[
(x-1)^{13}(x+1)^7(x^2+1).
\]

Its square is the identity on $N$ and $-I$ on $T$, hence is exactly the
canonical non-symplectic involution associated with the `(20,2,1)` Picard
lattice.  This is an independent structural check on the order-four lift.

For `C02`, the discriminant action is the identity.  It therefore glues to
$-I_T$, giving a non-symplectic involution with characteristic polynomial

\[
(x-1)^{12}(x+1)^{10}.
\]

The graph glue is integral in both cases, and both actions use the same
marked K3 lattice and the same period.  Their fixed Néron--Severi lattices
are integrally the archived saturated fixed lattices, of ranks 13 and 12.

This illustrates a useful realization criterion: a Picard involution need
not lift as an involution.  Its discriminant action may instead demand a
higher-order transcendental lift while leaving the fixed divisor section
unchanged.

## Kähler and Torelli gate

Each Picard action permutes all 50 walls of the source chamber.  The exact
fixed interiors have minimum source-wall pairing one.  Clearing denominators
gives primitive invariant positive classes

| class | scale | square |
|---|---:|---:|
| `C01` | 4 | 32,068 |
| `C02` | 2 | 5,386 |

The complete mixed reflection chamber lies inside a unique
$(-2)$-Weyl chamber.  Since each action preserves the mixed chamber and a
strict positive point, it preserves that K3 Weyl chamber.  After choosing
the marking so that this is the Kähler chamber, each glued Hodge isometry
fixes an ample class.  Strong global Torelli therefore realizes both by
genuine automorphisms of the projective K3 surface.

The source 50-wall chamber includes mixed root lengths.  This argument makes
it an automorphism-invariant divisor/reflection chamber; it does **not**
assert that the 50-wall chamber itself is the literal nef cone.

## Literal finite-volume handoff

The automorphism construction is separate from the volume proof.  The
archived literal Vinberg endpoint audit gives:

| section | facets | elliptic edges | finite vertices | cusps | bad edges / cusp neighbourhoods |
|---|---:|---:|---:|---:|---:|
| `H12/19` | 19 | 3,846 | 527 | 25 | `0 / 0` |
| `H11/20` | 20 | 2,229 | 311 | 29 | `0 / 0` |

Every elliptic rank-$(n-1)$ edge has exactly two finite-or-ideal endpoints,
and every affine component has a full parabolic zero-neighbourhood of total
rank $n-1$.  Hence both chambers have finite volume.  Their positive cusp
counts prove that both are noncompact.

## Novelty boundary

Alexeev's source chamber and the general theory of non-symplectic K3
automorphisms are classical.  The new project-level result is the exact
equivariant gluing of these two specified chamber symmetries and the resulting
geometric realization of their `H12/19` and `H11/20` sections.  No new
abstract-polytope, first-in-dimension, new-lattice, or new-arithmetic-class
claim is made.  See Alexeev's source framework
[*Reflective hyperbolic 2-elementary lattices, K3 surfaces and hyperkähler
manifolds*](https://arxiv.org/abs/2209.09110) and the classical order-four K3
framework of Artebani--Sarti
[*Symmetries of order four on K3 surfaces*](https://arxiv.org/abs/1102.4436).

## Reproduction

- [exact certifier](../scripts/certify_a20_c2a_c2b_k3_automorphism_lifts.py)
- [JSON certificate](../data/enumerations/a20_c2a_c2b_k3_automorphism_lifts.json)
- [regression tests](../tests/test_a20_c2a_c2b_k3_automorphism_lifts.py)

Run from the repository root:

```bash
PYTHONPATH=.:src python3 scripts/certify_a20_c2a_c2b_k3_automorphism_lifts.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_a20_c2a_c2b_k3_automorphism_lifts.py
```

Frozen certificate payload SHA-256:

```text
be9baceabfeb2c12364fc546815d62b09f714580e998685fac6214a8286388c1
```
