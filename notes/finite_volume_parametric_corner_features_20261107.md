# Twenty-six corner features and their vertex sums

This note specifies the new source
[finite_volume_parametric_corner_features_20261107.py](../scripts/finite_volume_parametric_corner_features_20261107.py).
It is a coefficient interface for joint weight search, not an executed
optimization or a new face-payment theorem. The existing current-law
coefficient catalogue and its successful independent verification are
authenticated inputs. Their coefficient verification is not rerun.

The supplied catalogue covers finite component ranks at most 58 and affine
component ranks at most 57. The public ambient range is $3\le n\le58$.
The claim of compatibility with the final CURRENT rule is restricted to
$n=58$. At other supported dimensions these are formal laws with no
inherited face-payment claim. A larger ambient range requires a new complete
catalogue and verification; the present file is not extrapolated.

## Feature coordinates and decoration

For a corner $c$, the provider returns a vector $\phi(c)\in\mathbb Q^{26}$.
For a coefficient vector $\theta\in\mathbb R^{26}$ the new charge is
$w_\theta(c)=\theta\cdot\phi(c)$. The feature values are independent of
$\theta$. Their stable IDs are exported as FEATURE_IDS, with coefficient
column names obtained by prefixing each ID with theta_.

Each component copy carries the native node numbering of the authenticated
catalogue. At a vertex, this decoration is chosen once and used for every
incident corner. Relabeling physical nodes must carry this decoration
along. No unproved invariance under arbitrary template automorphisms is
assumed. Ordered tiling positions must be converted by their explicitly
proved local bijection before using this native-coordinate interface.

The first four coordinates are:

| Index | ID | Feature |
|---|---|---|
| 0 | constant | One at every corner. |
| 1 | finite_baseline | The literal within-component weight if both marks lie in the same finite component; zero otherwise. |
| 2 | ideal_baseline | The literal affine weight of the marked pair in the complete varying affine component; zero at finite vertices. |
| 3 | ordinary_cross | The literal ordinary finite cross weight if the marks lie in distinct components; zero otherwise. |

The finite baseline includes the current ordinary $B_{13}$ replacement.
It is the literal same-component table before contextual additions, not the
older preliminary table.

Indices 4 through 18 are separate Boolean features, in the exact order of
CurrentWeightTables.same_deltas:

| Index | Marked component and unordered native pair | Guard | Current coefficient |
|---|---|---|---|
| 4 | $A_{13},(0,12)$ | an $A_1$ exists | $1$ |
| 5 | $A_{14},(0,13)$ | an $A_1$ exists | $1/2$ |
| 6 | $A_{15},(0,13)$ | an $A_1$ exists | $1/4$ |
| 7 | $A_{15},(1,14)$ | an $A_1$ exists | $1/4$ |
| 8 | $D_{15},(1,13)$ | an $A_1$ exists | $1/4$ |
| 9 | $D_{15},(1,14)$ | an $A_1$ exists | $1/4$ |
| 10 | $A_{16},(1,12)$ | an $A_1$ exists | $5/8$ |
| 11 | $A_{16},(3,14)$ | an $A_1$ exists | $5/8$ |
| 12 | $D_{16},(0,14)$ | an $A_1$ exists | $1/2$ |
| 13 | $D_{16},(0,15)$ | an $A_1$ exists | $1/2$ |
| 14 | $D_{14},(0,12)$ | an $A_1,A_2$, or $A_3$ exists | $1/2$ |
| 15 | $D_{14},(0,13)$ | an $A_1,A_2$, or $A_3$ exists | $1/2$ |
| 16 | $D_{13},(0,11)$ | at least three $A_{15}/D_{15}$ components | $-1/2$ |
| 17 | $D_{13},(0,12)$ | at least three $A_{15}/D_{15}$ components | $-1/2$ |
| 18 | $D_{28},(0,12)$ | at least two $A_{15}/D_{15}$ components | $1/2$ |

Their IDs have the form same_A13_0_12_has_a1, with the component,
pair, and literal guard name substituted. A feature is one precisely when
both marks lie in one copy of the stated type, occupy the stated pair,
and satisfy the guard. The coefficient is an increment over feature 1.
In particular, the two $D_{16}$ increments are $1/2$, because their
ordinary baseline is already $1/2$. The $D_{13}$ increments are negative.

The final seven features concern distinct finite component copies:

| Index and ID | Marked type/node sets | Guard | Current coefficient |
|---|---|---|---|
| 19 cross_a4_a12 | $A_4:\{0,3\}$ and $A_{12}:\{1,10\}$ | $A_1$ exists | $1/8$ |
| 20 cross_a2_a4 | $A_2:\{0,1\}$ and $A_4:\{0,3\}$ | $A_1,A_{10}$ exist | $1/4$ |
| 21 cross_a2_d5 | $A_2:\{0,1\}$ and $D_5:\{0\}$ | $A_1,A_9$ exist | $1/2$ |
| 22 cross_a1_d5_first | $A_1:\{0\}$ and $D_5:\{0\}$ | $A_{11}$ exists, or both $A_9,A_2$ exist | $3/8$ |
| 23 cross_a1_d5_second | $A_1:\{0\}$ and $D_5:\{0\}$ | both $A_9,A_2$ exist | $5/8$ |
| 24 cross_a12_d16 | $A_{12}:\{0,11\}$ and $D_{16}:\{0\}$ | at least two $A_{15}/D_{15}$ components | $1/2$ |
| 25 cross_a2_a2_distinct | both nodes of each of two distinct $A_2$ copies | $A_{10}$ and at least two $A_1$ components exist | $9/32$ |

Each feature is Boolean, not an oriented-pair count. Interchanging the
two marks leaves it unchanged. Features 22 and 23 deliberately overlap.
Feature 25 vanishes for the two marks in a single $A_2$ copy.
All witness types in each guard differ from its marked component types.
Consequently the global component counts evaluate exactly the witnesses
left after removing the marked component copies.

All 26 feature coordinates lie in $[0,1]$. This does not put the new charge
$w_\theta$ in $[0,2]$, or even make it nonnegative. Those are separate
constraints on $\theta$ and the admitted corner domain. At the current
coefficient vector

$$
\theta^{\rm cur}=(0,1,1,1,\ \text{the 15 increments above},\
                 1/8,1/4,1/2,3/8,5/8,1/2,9/32)
$$

the coordinate definition reproduces every literal marked finite or affine
clause at rank 58. This is a source identity; its independent runtime
reconstruction check remains a separate execution obligation.

## Exact finite aggregation

Let $m_T$ be the number of component copies of finite type $T$, with
$\sum_T\operatorname{rank}(T)m_T=n$. Put

$$
S_T=\sum_{a<b}w_{\rm same}(T,a,b),\qquad
B_{TU}=\sum_{a,b}w_{\rm ordinary}(T,a;U,b).
$$

The authenticated ordinary cross support has $T\ne U$; all same-type
ordinary entries are zero. Summing the feature vector over all unordered
node pairs at the finite vertex gives

$$
\Phi_0=\binom n2,\quad
\Phi_1=\sum_Tm_TS_T,\quad
\Phi_2=0,\quad
\Phi_3=\sum_{T<U}m_Tm_UB_{TU}.
\tag{1}
$$

For each of the 15 same-context features, exactly one marked pair occurs
in each copy of its target component. Therefore its aggregate is
$m_T$ times its Boolean guard.

Write $a_j=m_{A_j}$, $d_j=m_{D_j}$ and let brackets denote a Boolean
indicator. The seven cross aggregates, in their declared order, are

$$
\begin{aligned}
\Phi_{19}&=4a_4a_{12}[a_1>0],\\
\Phi_{20}&=4a_2a_4[a_1>0,\ a_{10}>0],\\
\Phi_{21}&=2a_2d_5[a_1>0,\ a_9>0],\\
\Phi_{22}&=a_1d_5[a_{11}>0\ \text{or}\ (a_9>0,\ a_2>0)],\\
\Phi_{23}&=a_1d_5[a_9>0,\ a_2>0],\\
\Phi_{24}&=2a_{12}d_{16}[a_{15}+d_{15}\ge2],\\
\Phi_{25}&=4\binom{a_2}{2}[a_{10}>0,\ a_1\ge2].
\end{aligned}
\tag{2}
$$

The factors $4,4,2,1,1,2,4$ count choices of native nodes. They do not
include a second orientation multiplier. The last expression equals
$2a_2(a_2-1)$ times its guard.

The method finite_aggregate computes (1)--(2) directly, without enumerating
forests or node pairs. The method finite_factor_terms exposes individual
unary, ordinary-cross, and contextual count expressions. Every term carries
its feature index and its component-count support. Grouping the evaluated
terms by feature index reproduces finite_aggregate. A component type outside
the chosen ambient rank is absent and has count zero. The constant feature
is a source cost $\binom n2$ and creates no new interactions.

For any $\theta$, linearity gives
$\sum_c w_\theta(c)=\theta\cdot\Phi$. This is the required affine input
to a new network LP. It does not authorize numerical tail sharing or
freezing equal-rank maximizers from the old coefficient vector.

## Exact ideal aggregation

At an ideal corner, both marks lie in the complete varying affine factor
$\Gamma$ of rank $r$ and order $r+1$. Only features 0 and 2 can be nonzero.
If $Q_\Gamma$ is the sum of its literal affine marked-pair table, summing
the features within one factor gives

$$
\Psi_{\Gamma,0}=\binom{r+1}{2},\qquad
\Psi_{\Gamma,2}=Q_\Gamma.
$$

For a full ideal link with factor counts $m_\Gamma$ and total rank $n-1$,
put $V=\prod_\Gamma(r_\Gamma+1)^{m_\Gamma}$. Each pair in one factor
occurs once for each choice of a vertex in every other factor. Thus

$$
\Phi^{\rm ideal}_0=\frac{V(n-1)}2,\qquad
\Phi^{\rm ideal}_2=
 V\sum_\Gamma\frac{m_\Gamma Q_\Gamma}{r_\Gamma+1},
\qquad \Phi^{\rm ideal}_j=0\quad(j\ne0,2).
\tag{3}
$$

The normalized aggregate divides (3) by $V$. In particular, a variable
factor envelope must include every inequality

$$
c_r\ge\frac r2\,\theta_{\rm constant}
       +\frac{Q_\Gamma}{r+1}\,\theta_{\rm ideal\_baseline}.
\tag{4}
$$

No positivity or fixed upper cap on these coefficients is inferred from
the old law. Individual-corner nonnegativity remains a separate block.

## Authentication, joining, and validation obligations

The loader requires the complete catalogue and the hash-bound receipt
from the independent connected-template checker. It verifies the pinned
sources, receipt status and full coverage flags, and the consumed catalogue
hash. It uses the frozen table reader without importing the literal worker.
The same-context layout and increments must exactly match that reader.
Both geometric $I_2$ variants share a coefficient type; their geometric
distinction remains outside this weight interface.

When attaching these features to the existing selected-weight corner
circuit, keep its old CURRENT auxiliaries as fixed basis functions with
their already proved bounds. Form a new expression $\theta\cdot\phi$.
Do not transfer the old charge bound $[0,2]$ to this expression. Remove
the old CURRENT underpayment query from any domain intended for a universal
variable-weight certificate. No fixed-law finite-label normalization is
transferred to the new family.

The required independent reconstruction checks are:

1. Every within and ordinary marked-pair lookup agrees with its authenticated
   primitive table, and same-type ordinary cross entries remain zero.
2. Each same-context and cross-context feature has the declared support,
   witness guard, orientation symmetry and distinct-copy condition.
3. On complete witness forests, summing marked feature vectors agrees with
   (1)--(2), and the current coefficient dot product agrees with the pinned
   literal marked and aggregate rules at $n=58$.
4. Ideal marked-pair sums and product-link corner multiplicities give (3).
5. All source and catalogue bindings remain unchanged across evaluation.

These checks do not rerun the old finite maximum or ideal-profile optimum.
The module is one fixed-dimension component of the proposed threshold
framework. It does not yet cover dimensions above 58 or prove that increasing
a dimension threshold preserves infeasibility; that requires the declared
larger threshold system and its separate coverage proof.
