# A distinct-opposite pair law and its exact A14 relaxation barrier

20 October 2026. Paper derivation, independently checked by root.
No computation, chamber construction, dimension ceiling or build.

## 1. A useful necessary pair inequality

Use the notation of
`compact_actual_wall_reciprocal_shadow_20261020.md`: at an actual finite
vertex, C=B^(-1), d_i=C_ii, u_i=omega_i/sqrt(d_i), and
K_ij=u_i.u_j. Put r_i=tanh L_i for the actual finite outgoing edges.

If directions i and j have DISTINCT actual opposite walls, then

    r_i*r_j>=K_ij.                                    (1)

Indeed, for the wall c_i opposite i, the inverse-Gram Schur inequality
gives (C beta)_j>=(C_ji/d_i)(C beta)_i. Thus

    T_j(c_i)>=K_ij*T_i(c_i)=K_ij/r_i.

Distinctness permits the reciprocal-shadow bound r_j>=T_j(c_i),
proving (1). The argument is geometric and works for every finite
vertex diagram, including nonirreducible ones. It needs no arithmetic.

For a whole A_n vertex, N=n+1 and i<=j give

    K_ij=sqrt(i*(N-j)/(j*(N-i))).                     (2)

Thus differently colored indices, colored by their actual opposite
wall, satisfy

    -log r_i-log r_j <= |log(i/(N-i))-log(j/(N-j))|/2.

The familiar four-directions-per-wall relaxation follows using the
standard compact Lanner-order-at-most-five input: every minimal
non-positive-definite subdiagram in I+c must contain c and EACH
direction closed by c, because deleting that direction leaves a
positive definite actual endpoint Gram. It has at most five nodes.

## 2. These conditions alone do not force energy at least one

Here is an EXACT counterconfiguration to the relaxation consisting
only of (1), at most four indices of each color, and the farthest-
vertex energy test. It is formal edge data, not an actual polytope.

At A14, let u_i be its unit dual vectors and define

    x=(129/250)*(u_7+u_8),
    tau_i=x.u_i,      r_i=2*tau_i/(1+tau_i^2).

Since K_78=7/8,

    ||x||^2=(129/250)^2*(15/4)=49923/50000<1.          (3)

All tau_i lie in (0,1), so each is the half-tangent of a positive
finite edge length. Put

    A=149769/140000.

Equation (2) gives

    tau_i^2=A*i/(15-i)       for 1<=i<=7,
    tau_(15-i)=tau_i.

In particular max tau_i^2=7A/8<1. Assign the four colors

    {1,2,3,4},   {5,6,7},   {8,9,10},   {11,12,13,14}.

Each color has at most four indices. ALL pairs of distinct colors
satisfy (1), in fact strictly. To verify this without a search:

* On the left half, K_ij=tau_i/tau_j for i<j. The quotient
  r_i*r_j/K_ij=4*tau_j^2/[(1+tau_i^2)(1+tau_j^2)] decreases with i
  and increases with j. The worst distinct-color pair is therefore
  i=4,j=5, whose inequality reduces exactly to

      25A-4A^2-22>0.

  The right half is symmetric.
* For i<=7 and j>=8, tau_i*tau_j=A*K_ij. Thus the pair inequality
  is equivalent to 4A>(1+tau_i^2)(1+tau_j^2). Its worst pair is7,8,
  reducing exactly to

      144A-49A^2-64>0.

Both rational inequalities hold. For example
1069/1000<A<107/100, and both displayed polynomials are increasing
on this interval; their values at1069/1000 are positive.

On the other hand x itself satisfies EVERY constraint x.u_i>=tau_i
with equality, and its norm is below one by (3). Hence the constrained
covering energy is at most49923/50000<1, and every nonnegative local
weighted farthest-vertex inequality is satisfied at this data.

## 3. The missing geometric obligation is explicit

The counterconfiguration does NOT provide actual opposite normals.
In particular it need not realize the exact same-wall equations
T_i(c)=1/r_i on each color, the nonnegative native contact vector of
that wall, its unit norm and reflection arithmetic, or the interiors
of all other actual halfspaces. No compactness or finite-volume claim
is made for it.

Therefore (1) is a useful necessary organizer, but replacing each actual
opposite by only a color of size at most four loses essential data.
A proof excluding whole A_n corners at rank14 must retain additional
actual-wall information, such as the exact reciprocal profiles,
finite-continuation restrictions, simultaneous arithmetic pairings, or
the interior-point inequalities x.v_c<a for c=a*p-v_c.
This barrier prevents a proposed relaxation from being mistaken for
a universal compact dimension argument; it leaves all current bounds
unchanged.
