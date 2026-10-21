# All-center obstruction for each single H4 weight ray

The exact certificate passes. Fix any one of the four fundamental-weight
directions of the H4 factor in the eleven-wall Makarov H5 + H4 seed. No finite
Coxeter-compatible collection of caps in that direction can produce finite
volume, even with arbitrary old centers and every allowed angle or dotted
distance. **This does not exclude using several different directions together,
mixed transverse profiles, or deleting seed walls.**

The [script](../scripts/fresh_h9_single_ray_cover_gate_20261013.py) and
[certificate](fresh_h9_single_ray_cover_gate_20261013.json) replay the
[source geometry](fresh_h9_makarov_face_extension_20261013.md) and check all
algebraic inequalities over Q(sqrt(5)). The following geometric argument
explains how those finite checks imply the unrestricted statement.

## Reduction to balls in the compact old prism

Let P be the compact Makarov H5 prism, G its Lorentz form, H the norm-two H4
Gram matrix, and B=H^-1. Write omega_j=B e_j. Normalize each cap to square two.
For one fixed direction j it has the form

    r=(u,-t omega_j),     u^2=-(t^2 B_jj-2),     t>=1.

The old inequalities put u in the cone over P; its nonzero vectors are all
future timelike. Thus u=U p, where p is a unit old hyperbolic point in P and
U=sqrt(t^2 B_jj-2). Coxeter compatibility with the jth H4 wall permits
t=2 cos(pi/m), m>=3, or t>=2 for parallel/dotted walls.

For a unit old point x in P, consider the ideal vector

    z_i(x)=(x,-omega_i/sqrt(B_ii)).

This satisfies all eleven seed inequalities. A cap has nonnegative pairing
with z_i(x) precisely on a closed hyperbolic ball, provided its radius exists:

    cosh R_i(t) = (B_ij/sqrt(B_ii))/sqrt(B_jj-2/t^2).

If the right side is less than one, the cap misses the entire ideal family
strictly. Set D=B_jj-B_ij^2/B_ii; a nonempty ball requires 2/t^2>=D.

For two caps on the same ray, their nonpositive mutual pairing implies

    cosh d(p,p') >= B_jj/
       sqrt((B_jj-2/t^2)(B_jj-2/s^2)).

Comparing this bound with cosh(R_i(t)+R_i(s)) shows that their closed coverage
balls are strictly separated whenever

    D(t^2+s^2)>2.

Indeed, after multiplying by the positive denominator, this is the inequality
D>sqrt((2/t^2-D)(2/s^2-D)); squaring gives the displayed criterion.

## The four cases

The checked worst axes and their exact D values are:

| Cap direction j | Tested ideal axis i | D |
|---|---|---|
| 0 | 3 | (9+3 sqrt(5))/4 |
| 1 | 3 | 3+sqrt(5) |
| 2 | 3 | (5+sqrt(5))/4 |
| 3 | 0 | 3/4 |

For j=0,1, D>2. Every cap misses the whole indicated ideal family strictly.

For j=2, only label m=3 can meet this ideal family. Any two of its coverage
balls are strictly separated. The certificate gives two old vertices whose
distance exceeds twice the ball radius, so no individual ball covers P.
A finite collection of separated closed sets cannot cover connected P.

For j=3, only m=3,4,5 can meet the indicated ideal family. Every m=4 or m=5
coverage ball is strictly separated from every other coverage ball, including
all m=3 balls: t^2+s^2>=3 and 3D>2. Again the two certified old vertices
exclude any one of these isolated balls covering P. It remains to show that
the m=3 balls alone cannot cover P.

## Why j=3, m=3 centers must be old vertices

Here -u^2=B_33-2=1+sqrt(5). Let T be the old walls orthogonal to u and S its
nonempty complement among the seven old walls. Every nonzero Coxeter pairing
has magnitude at least one for the six norm-two old walls, or
sqrt((7-sqrt(5))/22)>2/5 for the polar wall. Relax the last lower bound to 2/5.
Consequently u belongs to the polyhedron

    (u,W_i)=0 for i in T;
    (u,W_i)<=-ell_i for i in S,
    ell_0=...=ell_5=1, ell_6=2/5.

Its recession cone is a face of the cone over compact P. This cone is pointed
and lies strictly inside the future timelike cone, apart from its origin.
Every nonempty such polyhedron therefore has vertices and is a convex
combination of its vertices plus a recession vector. The positive Lorentz
norm sqrt(-u^2) is concave on the future cone and increases when a nonzero
future recession vector is added (reverse triangle inequality). Its minimum
on this polyhedron is therefore achieved at a vertex.

The script exhausts all 127 nonempty S. Each vertex is obtained by making
six of the seven constraints active, necessarily including all of T; all
feasible candidates are retained. There are 34 empty relaxed polyhedra.
For every feasible support with at most four zero pairings, **every** such
vertex has -u^2>1+sqrt(5). (The smallest value is approximately 4.20433.)
These supports cannot contain our u. The only possible supports have five
zero pairings and are actual old vertices. This argument bounds the entire
continuous relaxed polyhedron; it imposes no label or field bound.

For all ten old vertices, the certificate checks that the normalized old
part of the previous strict integer escape is strictly outside its j=3,m=3
coverage ball. Thus even all ten possible vertex-center balls do not cover P.

## From failed coverage to infinite volume

For j=3 the finite union of m=3 balls fails to cover P. Adding isolated
m=4,5 balls cannot repair this: their intersections with P are separated
closed components, none equal to P, whereas P is connected. The same
argument already handled j=2. Hence there is a point of P outside every
closed coverage ball. The complement is relatively open, so one can choose
that old point in the interior of P.

All added cap inequalities are then strict at the corresponding ideal vector.
Only some H4 seed inequalities are equalities. A small perturbation into the
interior of the H4 sector on the ideal sphere preserves the finitely many
strict inequalities. This gives an open patch of ideal boundary in the
chamber, and therefore infinite hyperbolic volume.

Replay:

```sh
python3 scripts/fresh_h9_single_ray_cover_gate_20261013.py --verify
```

The practical consequence is to retain the mixed-direction search. This
result neither establishes a new polytope nor excludes the full H9 seed.
