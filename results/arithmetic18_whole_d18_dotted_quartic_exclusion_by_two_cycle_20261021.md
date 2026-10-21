# Whole D18: the dotted head cannot have quartic ground field

21 October2026. New manual supporting proof, pending full outside-author
review. No mathematical worker, symbolic computation, enumeration, Lean,
TeX or manuscript change was used. The composition reviewer discovered
the stronger two-cycle integrality lemma below; root independently
expanded it and derived the norm contradiction. Earlier actual-wall
normalization and fork foundations have shared authorship.

**Theorem.** At an actual whole-D18 farthest vertex of a compact
arithmetic first-type hyperbolic Coxeter18-polytope, a dotted long head
has actual ground-field degree at most three. In particular the quartic
case remaining in the previous reduction is impossible. Quadratic and
cubic dotted heads, whole D18 and arithmetic dimension18 remain open.
This is not a new global ceiling or a new example.

## 1. Exact inputs

All input paths in this table are under `results/`.

| Input | SHA256 | Restricted interface |
|---|---|---|
| arithmetic18_whole_d18_dotted_head_degree_four_and_isolated_forks_20261021.md | 0590f2b9d8747f6edbe0cafd74f2909c54c87f62e5ebe7d52081abbdfe35fa16 | Actual degree at most4; BOTH forks isolated/dotted at degree3/4; raw time squares and actual field equality; conjugate strips and quartic N(V)<=6 |
| arithmetic_uniform_classical_weak_helper_obstruction_20261021.md | 535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905 | Actual K-valued root normalization, cyclic integrality and positive nonphysical forms |

Both are used only at their stated actual-wall interfaces. No field
classification, integral power basis, discriminant estimate, root
enumeration, distinctness of the TWO fork opposites or independent sum
of two fork energy bounds enters this proof. We use either one actual
fork and its actual long-head wall.

Write O for the ring of integers of the ACTUAL totally real field K.
The isolated dotted long head supplies

    U in O, K=Q(U), 4<U<28/3,
    0<U_sigma<2 at EVERY nonphysical embedding.          (1)

Suppose d=[K:Q]=4. Either fork is isolated and dotted. Its parameter
satisfies

    V in O, 4<V<676/9,
    0<V_sigma<4/9 at EVERY nonphysical embedding.        (2)

All these statements concern the same actual farthest vertex and the
same field, not abstract compatible Gram vectors.

## 2. The two actual cycles give UV/4 integral

Use old D18 roots alpha_i of square2, with product-1 on each old edge.
Let r be the K-valued raw long-head normal and s the raw normal at
the chosen fork. The input normalization gives

    r²=2U, (r,alpha_1)=-U,
    s²=2V, (s,alpha_f)=-V,

with all other old products zero. The old fork is f=17 or18.
These two actual caps are distinct: their sole old contacts are at
different old nodes and are nonzero. No relation between the two
fork caps is needed.

For completeness, directed entries are

    b_ij=2(r_i,r_j)/(r_i,r_i).

Products around a closed directed cycle agree with the usual
unit-normal Gram-cycle products, since the root rescalings cancel.
They are integral for this actual arithmetic chamber. Put a=(r,s).
The actual cycle r--alpha_1--...--alpha_f--s--r has product plus or
minus a: its first entry is-1, its old-path entries are-1, and its
last two entries are-V and a/V. Thus

    a in O.                                             (3)

If a=0 the same conclusion is immediate. The actual TWO-cycle
r--s--r, separately from the long cycle, gives

    L=a²/(UV) in O.                                    (4)

It is important to use (4); the earlier parity proof used (3) alone.

The old D18 Cartan inverse has (1,f) entry1/2 and fork diagonal9/2.
The spatial projection product of r and s is consequently UV/2.
Their negative time-residual squares are

    H_r²=U(U-2), H_s²=V(9V-4)/2.

Both residuals lie on the SAME one-dimensional K time line. Squaring
their residual product gives the exact actual identity

    (UV/2-a)²=U(U-2)V(9V-4)/2.                         (5)

Multiplying by4 and reducing in O modulo4 yields (UV)² in4O.
Indeed the other right-hand terms in

    4a²-4aUV+(UV)²=2UV(U-2)(9V-4)

are multiples of4 except 18(UV)²; their difference leaves the odd
coefficient17 times (UV)². It follows that (UV/2)² is integral.
Since O is integrally closed and UV/2 belongs to K, we obtain

    UV/2 in O.

Set T=UV/2-a, so T is integral. Dividing (5) by UV and expanding
(4) now gives

    L+T=19UV/4-2U-9V+4.                               (6)

All terms except possibly UV/4 are integral; hence 19UV/4 is
integral. Together with UV/2 integral, subtraction of nine times
the latter gives the stronger necessary condition

    UV/4 in O.                                        (7)

We also retain, without asserting V/2 integral,

    (2-U)V/2=V-UV/2 in O.                             (8)

These derivations hold for an isolated head and fork at the same
rank18 in any actual degree for which the displayed interfaces apply.
No rational-square-norm substitute for being a square in K is used.

## 3. In degree four their norms contradict the entire conjugate strip

Let

    A=N_K/Q(U)>0, B=|N_K/Q(2-U)|>0, n=N_K/Q(V)>0.

They are positive integers. Equation (2) gives

    n < (676/9)(4/9)^3 = 43264/6561 <7,

so n<=6. From the norms of the NONZERO integral quantities (7)
and (8), respectively, we have

    4^4 divides A*n, 2^4 divides B*n.

Multiplying the resulting positive integer inequalities yields

    A*B*n² >= 8^4=4096,
    A*B >= 4096/36 = 1024/9.                           (9)

On the other hand, the actual conjugate strip (1) gives a much
smaller upper bound without any case division. For 0<u<2,

    0<u(2-u)=1-(u-1)²<=1.

Therefore, taking the absolute norm of the head time square,

    A*B = U(U-2) * product_(sigma nonphysical)
                         [U_sigma(2-U_sigma)]
        <= U(U-2)
        < (28/3)(22/3) = 616/9.                       (10)

The strict physical inequality uses U>4 and U<28/3. The lower
bound1024/9 in (9) exceeds616/9, a contradiction.

Thus d=4 is impossible. Combining with the already proved d<=4
leaves d<=3. The rational case is independently excluded by the
usual anisotropy criterion in these dimensions; it is unnecessary
for this contradiction. The unresolved quadratic and cubic cases
are not disposed of by this note. No script or Lean execution is
credited, and all six finished paper editions remain unchanged.
