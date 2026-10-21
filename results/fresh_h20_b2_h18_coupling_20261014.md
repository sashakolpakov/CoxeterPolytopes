# A reverse-B2 lift of the complete H18 hybrid: a uniform obstruction

The experiment gives a coherent **42-facet Coxeter chamber in H20**, with a complete H18 source face and finite contacts to its transverse B2 support. It also proves a stronger stopping result: **every direct reverse-B2 coupling of this fixed source has a common strict spacelike escape**, regardless of which admissible inducing profile is assigned to each source facet. Thus none of these direct couplings has finite volume.

The [exact checker](../scripts/fresh_h20_b2_h18_coupling_20261014.py) and [certificate](fresh_h20_b2_h18_coupling_20261014.json) verify the physical construction, all 861 pairs of the coherent example, actual facet points, and the common escape against all 74 possible support/inducing normals. Generation and saved replay pass. Arithmetic uses the repository's exact Q(sqrt(2)) implementation; no second independent implementation is claimed.

This changes the transverse contact structure and retained walls completely from the 133-wall branch. It is not another endpoint-cap bank or another assignment of signed height pairs. The universal conclusion concerns **direct couplings**, containing only the two B2 support walls and walls inducing the existing forty source facets. Additional walls not inducing source facets, different intrinsic sources and general H20 constructions remain outside it.

## Source and complete inducing profiles

The intrinsic source is the frozen forty-wall Vinberg H18 hybrid. Its primary construction glues two arithmetic parents along a common orthogonal facet; the source itself is finite-volume. [Vinberg, Sections 3, 5 and 7](https://sfb701.math.uni-bielefeld.de/preprints/sfb14047.pdf).

Normalize the two transverse simple roots to square two:

\[
B=\begin{pmatrix}2&-\sqrt2\\-\sqrt2&2\end{pmatrix},\qquad
B^{-1}=\begin{pmatrix}1&1/\sqrt2\\1/\sqrt2&1\end{pmatrix}.
\]

Let their fundamental weights be omega0, omega1. Each has square one and their product is 1/sqrt(2). This reverses an allowed B2 face operation; the spherical-extension and projection rules are stated in [Allcock, Theorems 2.1–2.2](https://arxiv.org/pdf/0903.0138), printed pages 743–744. The construction and its obstruction below are checked directly rather than inferred from that theorem.

An ambient wall inducing a facet of this H18 face has positive projected square. Normalize its ambient square to two. A nonright Coxeter contact with a support root has magnitude at least one. The nonnegative inverse then leaves exactly three possibilities: no support contact, or one label-three contact at either support root. Two nonzero contacts cost at least 2+sqrt(2)>2; a single label at least four costs at least two. These are all-label restrictions, with no search cutoff.

For a physical source normal r of square q, the three profiles are

\[
\sqrt{2/q}\,r,\qquad r/\sqrt q-\omega_0,\qquad r/\sqrt q-\omega_1.
\]

All have ambient square two. If a source facet has a source pi/3 neighbor, both must be ordinary: one special endpoint gives ambient product -1/sqrt(2), strictly between zero and the smallest nonzero Coxeter product; two special endpoints give a positive product. The checker supplies a literal source pi/3 neighbor for each of the 24 forced ordinary facets.

The only potentially special facets are the fourteen norm-one source walls

`L18, L25, L26, L27, L28, L29, L30, L31, L32, L33, L34, L35, L36, L37`

and the two isolated norm-two walls `R23,R24`. The union of all possible profiles therefore has 24 ordinary normals plus 16 triples: 72 inducing normals, together with the two support roots. Some choices within this 74-normal union are mutually incompatible. The common-escape proof deliberately checks the larger union, so compatibility restrictions can only strengthen its scope.

## One literal coherent construction

Choose the first special color for all fourteen norm-one source walls and leave the other 26 source walls ordinary. Add the B2 support. The 42 normalized normals have these exact pair counts:

| Relation | Count |
|---|---:|
| Right | 641 |
| pi/3 | 47 |
| pi/4 | 3 |
| Parallel | 48 |
| Dotted | 122 |
| Total | 861 |

A common strict timelike center and all nonpositive pairings give actual facets; the certificate explicitly projects the center to each wall and checks all other inequalities strictly. The old source center in the transverse zero plane is strict on all forty inducing walls and zero on the two support roots. Thus the B2 face is actual and isometric to the full intrinsic H18 chamber.

The unchanged `L1,R10` Cartan cycle is 12+8sqrt(2), while the unchanged `R23,R24` rational Gram has determinant -32. The earlier hybrid's conditional non-quasi-arithmetic argument consequently survives: a hypothetical finite-volume completion retaining these walls would not be quasi-arithmetic. This is not a claim that the present infinite-volume group is a nonarithmetic lattice. No arithmetic root enumeration is used.

## One escape for every direct profile assignment

The certificate gives a literal source vector x over Q(sqrt(2)), in the same nineteen-root basis as the archived hybrid. It satisfies

\[
x^2=(3857-2800\sqrt2)/529<0,
\qquad (x,r_i)\le0\quad\hbox{for every source wall},
\]

and, for all sixteen potentially special walls,

\[
(x,r_i)\le-\sqrt{r_i^2}.
\]

Let h be the sum of the forty source normals. Its exact values are

\[
h^2=-2952-1280\sqrt2,\qquad (x,h)=-16-40\sqrt2,
\]

and h is strict on every source wall. Put

\[
w=-\tfrac12(\omega_0+\omega_1),\qquad
E=(x+h/1000)+w.
\]

Here w is strict on both B2 support roots and

\[
w^2=\tfrac12+\tfrac{\sqrt2}4,
\qquad (w,\omega_j)=-\tfrac12-\tfrac{\sqrt2}4.
\]

Since this last magnitude is less than one, E is strict against either special profile and every ordinary profile. Its exact square is

\[
E^2=\frac{512876299}{66125000}
-\frac{33884339}{6612500}\sqrt2>0.
\]

The checker verifies every one of the 74 strict inequalities directly. It also verifies a common strict future timelike center C and (E,C)<0. The segment from C toward E crosses the future null cone while every inequality remains strict, producing an open ideal-boundary region. Every direct coherent coupling is therefore infinite-volume. This proof covers all assignments without enumerating them.

## A reusable global test for a different source

The same mechanism supplies a finite completion criterion for any other proposed B2 source coupling. Write a transverse point as `u=−a omega0−b omega1`, with a,b nonnegative. A special facet of color zero bounds `a+b/sqrt(2)` by its normalized source distance; color one bounds `a/sqrt(2)+b`. Ordinary facets bound only the source point. The full chamber is therefore an explicitly defined polyhedral cone over the source cone with a two-dimensional polygon fiber, and its Lorentz square is

\[
x^2+a^2+b^2+\sqrt2\,ab.
\]

When at least one special facet is present, the fiber is bounded over every bounded source slice. An exact extreme-ray census of this cone is a global test: all rays future causal gives a finite-volume polytope; any feasible spacelike ray gives an infinite-volume region after a strict-center perturbation. This avoids treating finite local angles or a complete lower-dimensional face as sufficient for ambient closure. For the present source the literal escape already settles the test, so no expensive full ray census is needed.

```sh
python scripts/fresh_h20_b2_h18_coupling_20261014.py --verify
```
