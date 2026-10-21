# The exact reciprocal shadow of an actual wall

20 October 2026. Own paper lemma, independently derived with root.
No computation, new dimension ceiling, build or formal verification.

## 1. Statement and conventions

Let p be an actual vertex of a compact hyperbolic Coxeter polytope P.
Use the Lorentz convention p.p=-1, outward unit normals n_i at p,
their positive definite Gram matrix B, dual vectors omega_i, and

    C=B^(-1),    d_i=C_ii,    u_i=omega_i/sqrt(d_i).

The actual outgoing i-edge has finite length L_i>0 and endpoint

    p_i=cosh(L_i)p-sinh(L_i)u_i.

Let c be ANY actual facet normal not through p. Write

    c=a p-sum_j beta_j omega_j,
    a=-p.c>0,    beta_j=-n_j.c>=0,
    t=C beta,    T_i(c)=t_i/(a*sqrt(d_i)).              (1)

All contacts are nonnegative by all-wall acuteness. The conclusion is
the following exact alternative, separately for every direction i:

* If T_i(c)>1, the wall c ITSELF is the opposite wall of this edge,
  and tanh(L_i)=1/T_i(c).
* If T_i(c)<1, c is not its opposite wall, and
  tanh(L_i)>=T_i(c).
* The value T_i(c)=1 is impossible in a compact P.

Thus one known actual wall gives a lower bound for every OTHER edge,
and also identifies exactly all directions that this same wall closes.
Directions must not be counted as distinct walls: one wall can close
several outgoing edges.

The statement uses neither arithmetic nor a classification of finite
continuation diagrams. It applies at any actual vertex; it is not
restricted to the farthest vertex.

## 2. The inverse-Gram shadow inequality

For each i, the matrix

    C-C_.i*C_i./d_i                                    (2)

has zero i-th row and column; on the other indices its block is exactly
B_(I\{i})^(-1). This is the block-inverse Schur identity. The latter
inverse is entrywise nonnegative: a positive definite Coxeter Gram has
the form I-A with A entrywise nonnegative and spectral radius below1,
so its inverse is the convergent sum of the nonnegative matrices A^k.
The same applies to each principal submatrix.

Write a DISTINCT actual wall d not through p as

    d=b p-sum_j gamma_j omega_j,
    b>0,    gamma>=0,    q=C gamma.

Multiplying (2) by beta and gamma gives

    S=beta^T C gamma >= t_i*q_i/d_i.                   (3)

Since c and d are distinct actual outward normals, acuteness gives

    c.d=S-a*b<=0,    hence a*b>=S.                     (4)

Suppose d is the actual opposite of edge i. Its endpoint equation is

    -b*cosh(L_i)+(q_i/sqrt(d_i))*sinh(L_i)=0,
    tanh(L_i)=b*sqrt(d_i)/q_i.                         (5)

In particular q_i>0. Equations (3)--(5) prove

    tanh(L_i)>=t_i/(a*sqrt(d_i))=T_i(c)                (6)

whenever this opposite wall is DISTINCT from c. The restriction is
essential: applying acuteness to c with itself would be false. Nor
does (6) apply to old incident walls, whose contact vectors need not
be nonnegative in these coordinates.

## 3. Why the alternative is actual, including its boundary

Along the complete outgoing edge ray, the product with c is

    c.p(s)=-a*cosh(s)+(t_i/sqrt(d_i))*sinh(s).          (7)

If c is the actual opposite, (7) gives tanh(L_i)=1/T_i(c), so T_i(c)>1.
Conversely, if T_i(c)>1 but another wall closed first, (6) would force
tanh(L_i)>1. That is impossible. Hence c is the actual opposite.

If T_i(c)<1, (7) has no finite zero for s>0, so the actual opposite
is different and (6) applies. If T_i(c)=1, c cannot meet the ray at a
finite point, while the distinct actual opposite would force
tanh(L_i)>=1. Compactness rules out this last case as well.

This argument uses the existence of the actual finite outgoing edge
of P. It does not replace an actual endpoint with a merely formal
positive definite diagram. For comparison, the common-white Schur norm
of c is exactly

    delta_i=1-beta_W^T B_W^(-1) beta_W
           =t_i^2/d_i-a^2
           =a^2*(T_i(c)^2-1),    W=I\{i}.            (8)

Thus the same criterion T_i(c)>1 also says that c and the retained old
walls have positive definite Gram. At that endpoint (5) yields
a^2=delta_i*sinh^2(L_i), as in the usual edge formula.

## 4. An exact A2 mate relation

Suppose c has just one contact beta in an old A2 block; contacts to
other old components are allowed. Let c close the contacted direction
and put r=tanh(L) in (0,1). The A2 inverse is

    [[4/3,2/3],[2/3,4/3]].

Consequently T_contact=1/r and T_mate=1/(2r). Therefore:

    r<1/2:  c also closes the mate, tanh(L_mate)=2r;
    r>1/2:  its mate has another opposite,
            tanh(L_mate)>=1/(2r).

The value r=1/2 is forbidden by Section3. In the second case, if
sinh^2(L)<31/50, then r^2<31/81 and

    tanh^2(L_mate)>81/124>16/25,
    tanh^2(L_mate/2)>1/4.                             (9)

The last implication uses tanh(s)=2*tanh(s/2)/(1+tanh^2(s/2)),
which takes the value4/5 at tanh(s/2)=1/2. This funds the mate energy
without a separate native-angle analysis of its other opposite wall.

## 5. An isolated A_n wall gives an envelope on every edge

Suppose the whole old component is A_n and c has just one contact at
node k, with no other contacts. Let r_k=tanh(L_k), where c is the
actual opposite of direction k, and put

    K_ik=C_ik/sqrt(d_i*d_k).

Then (1) gives the exact reciprocal-shadow data

    T_i(c)=K_ik/r_k.

If K_ik>r_k, this SAME c closes direction i and
tanh(L_i)=r_k/K_ik. If K_ik<r_k, its actual opposite is different and
tanh(L_i)>=K_ik/r_k. Equality is forbidden. This identifies all the
same-wall directions before applying any compatibility argument to
distinct opposite normals.

For N=n+1 and i<=k, the scale-independent A_n formula is

    K_ik=sqrt(i*(N-k)/(k*(N-i)))
        =exp(-|log(i/(N-i))-log(k/(N-k))|/2).          (10)

The closed indices therefore form a contiguous interval in these
logarithmic coordinates. Near a transition T_i=1, either branch forces
a long edge. How much simultaneous energy this supplies for every
possible finite component and continuation profile remains the global
dimension-bound obligation.

At a farthest vertex, these edge bounds can be inserted into the exact
covering energy of
`results/compact_farthest_vertex_energy_obstruction_20261020.md`.
That energy is below1 and adds over orthogonal components. The present
lemma supplies an actual-wall interface for that test, not an exclusion
of all finite vertex diagrams or a new compact dimension ceiling.
