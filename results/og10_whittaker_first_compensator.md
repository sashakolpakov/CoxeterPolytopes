# The first complete off-slice OG10 compensator shell is obstructed

## Result

Let

\[
F=\mathbf Q(\sqrt5),\qquad \mathcal O=\mathbf Z[\phi],\qquad
L=(\mathcal O^{12},\widehat Q)
\]

be the totally-real OG10 extension from
`og10_borcherds_obstruction.md`, and let

\[
m_0=\frac{\phi}{\sqrt5},\qquad \mu_0=0
\]

be the target norm-two Whittaker index.  Among **all** discriminant classes
whose added-line coordinate is nonzero, the smallest possible norm of a
totally positive Whittaker index is

\[
\boxed{N(m_1)=\frac9{880}}.
\]

Modulo multiplication by units and the tautological symmetry
`mu -> -mu`, the complete shell at this norm contains exactly two
discriminant-component orbits.  If `P_0`, `P_A`, and `P_B` are the three
weight-`(6,6)` Hilbert--Weil Poincare series indexed by the target and those
two off-slice orbits, then

\[
\boxed{P_0,\ P_A,\ P_B\text{ are linearly independent}.}
\]

Consequently no weakly holomorphic Whittaker principal part supported on
the target orbit and the **complete globally minimal off-slice shell** can
have a nonzero target coefficient.  Equivalently, the first finite
compensator search has no nontrivial divisor-level kernel.

This is a bounded nonexistence theorem.  It does not exclude compensators at
larger index norm or kernels involving several off-slice shells.

## 1. The complete discriminant coordinates

Only three non-unimodular field summands matter.  The first two coordinates
of the OG10 form have Gram matrix

\[
A=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix},\qquad
\det A=-(1+4\phi).
\]

Put `p=1+4phi`; then `N(p)=-11`.  The positive line in `H` has coefficient
`2`, while the added sign-changing line has coefficient `2(1-phi)`, which
is associated to `2`.  The `E_8` block is unimodular.  Hence

\[
L'/L\cong
\mathcal O/(p)\oplus\mathcal O/(2)_u\oplus\mathcal O/(2)_e,
\qquad |L'/L|=11\cdot4\cdot4=176.
\tag{1.1}
\]

Let `v=A^{-1}e_1`.  Represent a class by

\[
\mu(t,x,y)=t v+\frac{x}{2}u+\frac{\phi y}{2}e,
\]

where `t in F_11`, `x,y in O/(2)`.  Here `phi=8 mod p`, so integers
`t=0,...,10` represent all of `O/(p)`.  Direct inversion of `A` gives

\[
\sqrt5\,Q(v)=-\frac{\phi}{p}=\frac{\phi-4}{11}.
\]

The two line contributions are

\[
\sqrt5\,Q(xu/2)=\frac{\phi x^2}{4},\qquad
\sqrt5\,Q(\phi ye/2)=-\frac{\phi^2y^2}{4}.
\]

Thus

\[
\boxed{
\sqrt5\,Q(\mu(t,x,y))
\equiv
\frac{t^2(\phi-4)}{11}
+\frac{\phi x^2}{4}
-\frac{\phi^2y^2}{4}pmod{\mathcal O}.}
\tag{1.2}
\]

This is also the exact geometric filter.  Every vector in
`mu(t,x,y)+L` has added-line coefficient

\[
\frac{\phi y}{2}+a,\qquad a\in\mathcal O.
\]

It can vanish only if `y=0 mod 2`.  Therefore

\[
\boxed{y\ne0\text{ in }\mathcal O/(2)
\Longrightarrow Z(m,\mu(t,x,y))\cap\iota(\mathbb H^{10})=\varnothing.}
\tag{1.3}
\]

This statement concerns the distinguished real OG10 slice.  It is stronger
than merely finding one off-slice vector: the whole special divisor indexed
by the coset misses the slice.

## 2. Exact minimum of an off-slice index

For an admissible index put

\[
w=44\sqrt5\,m\in\mathcal O.
\]

Equation (1.2) gives the finite congruence

\[
w\equiv
4t^2(\phi-4)+11\phi x^2-11\phi^2y^2
\pmod{44\mathcal O},\qquad y\ne0.
\tag{2.1}
\]

There are `11*4*3=132` triples.  Exact enumeration produces `54` distinct
residue classes, and their field norms modulo `44` are precisely

\[
N(w)\equiv0,11,33\pmod{44}.
\tag{2.2}
\]

Because `m` is totally positive and `sqrt(5)` has signs `(+,-)`, `w` has
signs `(+,-)` and negative field norm.  Moreover

\[
N(m)=-\frac{N(w)}{44^2\cdot5}=-\frac{N(w)}{9680}.
\tag{2.3}
\]

If `0 < -N(w) < 99`, equations (2.2)--(2.3) leave only

\[
-N(w)\in\{11,33,44,55,77,88\}.
\tag{2.4}
\]

The exclusion of all six values is finite and exact:

- The Minkowski bound is `sqrt(5)/2<2`, so `O` has class number one.
- The primes `2`, `3`, and `7` are inert.  Hence `33`, `77`, and `88`
  cannot be element-norm magnitudes: an inert prime occurs to odd
  valuation.
- Put `w_0=-1+3phi`, so `N(w_0)=-11`.  Up to norm-one units, every element
  of norm `-11`, `-44`, or `-55` is generated respectively by

  \[
  w_0,\bar w_0;\qquad
  2w_0,2\bar w_0;\qquad
  \phi\sqrt5w_0,\phi\sqrt5\bar w_0.
  \]

  The positive norm-one unit `phi^2` has order `15` modulo `44O`.
  Checking these six finite `15`-term sign-correct residue orbits against
  (2.1) gives no intersection.

At the next value,

\[
w_1=3w_0=-3+9\phi,qquad N(w_1)=-99,
\]

and its residue is admissible.  Therefore

\[
m_1=\frac{-3+9\phi}{44\sqrt5}
=\frac{21+3\phi}{220},
\qquad N(m_1)=\frac9{880}.
\tag{2.5}
\]

The other prime-factor orbit `3 bar(w_0)` has no sign-correct admissible
residue.  Thus (2.5) is the global minimum, not merely the minimum in a
chosen discriminant component.

## 3. The two minimal component orbits

Two representatives over `m_1` are

\[
\mu_A=\mu(4,0,1),\qquad
\mu_B=\mu(4,\phi,\phi).
\tag{3.1}
\]

Multiplication of a discriminant class by `phi` acts as

\[
(t,x,y)\longmapsto(8t,\phi x,\phi y),
\tag{3.2}
\]

and multiplies the index by `phi^2`.  Its fifteenth power sends a class to
its negative.  Exact orbit enumeration shows that every minimal admissible
class belongs to exactly one of the two disjoint `30`-element orbits of
`mu_A` and `mu_B`.  After quotienting by `mu -> -mu`, these are exactly the
two divisor-relevant modular-unit orbits asserted above.

## 4. The three-column Poincare coefficient matrix

Use the three indices

\[
(m_0,0),\qquad (m_1,\mu_A),\qquad (m_1,\mu_B),
\]

and define

\[
C_{ij}=b_{P_j}(m_i,\mu_i),\qquad i,j\in\{0,A,B\}.
\]

Each diagonal entry has identity term `1`.  Every off-diagonal identity
term is zero: the target has zero auxiliary component, while `mu_A` and
`mu_B` lie in distinct unit orbits.

The scalar Hilbert Poincare coefficient formula of Zhang--Zhang has outer
factor

\[
\left(\frac{N(m_i)}{N(m_j)}\right)^{5/2}.
\]

The vector-valued unfolding replaces each scalar Kloosterman phase by the
corresponding matrix coefficient of the unitary Weil representation.  Put

\[
s_i=N(m_i)^{5/2},\qquad
\widetilde C=\operatorname{diag}(s_i)^{-1}
C\operatorname{diag}(s_i).
\tag{4.1}
\]

Then the norm-ratio factors disappear.  Since similarity preserves
invertibility, it is enough to prove strict diagonal dominance for
`tilde(C)`.

For a unit lower-left row, every Weil `S`-matrix entry has magnitude
`1/sqrt(176)`.  For every other row, the generalized Kloosterman sum divided
by its ideal norm has magnitude at most `1`.  There is one harmless
normalization wrinkle: Bruinier's translation stabilizer contains the
metaplectic kernel but not the lift of `-I`.  Thus a nonzero component can
be represented by the paired seeds `mu` and `-mu`, depending on whether one
uses vector-valued or center-quotiented scalar rows.  The certificate bounds
both representatives separately, inserting a conservative factor `2` in
every estimate involving a minimal nonzero component.  The zero-component
target needs no factor: its central doubling also doubles the identity term.

With that convention-independent factor included, the exact
Bessel/ideal-tail certificate gives the following upper bounds (decimals
are only for readability):

| entry type | certified bound |
|---|---:|
| target diagonal error | `< 0.543203411990440` |
| either minimal diagonal error | `< 0.000031083273147` |
| target--minimal off diagonal | `< 0.007663335253610` |
| `A`--`B` off diagonal | `< 0.000031083273147` |

Consequently the target row of `tilde(C)` has dominance margin

\[
1-E_{00}-2E_{0A}>0.441469917502342,
\]

and each minimal row has margin

\[
1-2E_{AA}-E_{0A}>0.992274498200098,
\]

where the three `E` values are the exact rational error majorants displayed
in decimal form in the table.  The coarser rational bounds stored in the
machine certificate are respectively `68/125`, `1/31250`, and `1/125`.

The Levy--Desplanques theorem now gives

\[
\boxed{\det\widetilde C\ne0\quad\text{and hence}\quad\det C\ne0.}
\tag{4.2}
\]

Bruinier's equations (4.17)--(4.19) identify these columns, up to nonzero
individual normalizing scalars, with

\[
\xi(f_{m_0,0}),\quad
\xi(f_{m_1,\mu_A}),\quad
\xi(f_{m_1,\mu_B}).
\]

Therefore their linear independence is exactly the asserted absence of a
Whittaker kernel on this finite shell.  Bruinier's Proposition 4.5 then
identifies it with the absence of a principal-part solution to all cusp
coefficient obstructions.

The analytic inputs and normalizations are the same as in
`og10_whittaker_nonvanishing.md`: Bruinier's unitary Weil representation
and Poincare map are in
[Bruinier, *Regularized theta lifts for orthogonal groups over totally real
fields*](https://arxiv.org/abs/0908.3076), and the scalar
Kloosterman--Bessel formula is in
[Zhang--Zhang, *On the non-vanishing of Hilbert Poincare
series*](https://arxiv.org/abs/2312.04169).

## 5. Reproducible certificate

Run

```bash
PYTHONPATH=src:. python3 scripts/certify_og10_whittaker_first_compensator.py
```

The program checks the entire `132`-class off-slice discriminant search,
the six smaller norm possibilities, both complete minimal unit orbits, and
the exact matrix inequalities.  Finite Bessel sums use `python-flint`/Arb
at `180`-bit precision; their infinite unit tails and all nonunit ideal
tails are separately bounded by rational inequalities.  The final two row
dominance margins are exact positive `Fraction` values.

The regression test is
`tests/test_og10_whittaker_first_compensator.py`.

## 6. Scope

What is proved:

1. `y != 0` is an exact discriminant-coset criterion forcing the whole
   special divisor off the real OG10 slice;
2. `9/880` is the global smallest norm of such a totally positive index;
3. the entire minimal shell has two divisor-relevant unit orbits; and
4. adding both of them does not cancel the target cusp obstruction.

What is not proved:

1. no bound has been proved on the norm or number of arbitrary
   compensators;
2. larger off-slice shells may still produce a kernel; and
3. neither this obstruction nor a future successful lift decides the
   reflection index, finite covolume, or compactness of the real Coxeter
   chamber.
