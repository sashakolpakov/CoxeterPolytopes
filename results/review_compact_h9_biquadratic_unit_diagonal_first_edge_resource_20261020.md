# Outside audit of the biquadratic unit-diagonal first edge

20 October 2026. **PASS.** Full outside-author paper audit, with own
arithmetic and finite-strip derivations. No scientific execution,
worker import, root search, Lean or TeX job. This is internal review,
not independent discovery or external refereeing.

Bound source:
`results/compact_h9_biquadratic_unit_diagonal_first_edge_20261020.md`,
SHA256 `05b6473eeadb4befa89c82af2f7002c7441ada3d1878c93d4630355a42685efe`.
Put s=sqrt2, phi=(1+sqrt5)/2, u=phi-1, a=u+s, d=2-s.

## Arithmetic preflight and complete primitive pool

The product integral basis follows from coprime quadratic discriminants
8 and5, giving field discriminant1600 and Minkowski ideal bound15/4.
At2 there is one prime with ramification index2 and residue degree2:
the phi polynomial remains irreducible modulo(s). At3 the Frobenius
flips both quadratic factors and has order2, so both residue degrees
are2. There is consequently no ideal of prime norm2 or3, proving the
PID assertion from that bound.

I checked all four signs of a, the relative norm a*a_(s->-s)=-phi,
absolute norm-1, and the displayed inverse. The four unit sign vectors
are linearly independent over F2. Dirichlet rank3 with torsion{+1,-1}
then makes the sign map on units modulo squares bijective; every
totally positive unit is a square. This proves the normalization claim,
rather than merely exhibiting a few useful units.

The lattice is unimodular since a is a unit. For a primitive normal,
reflection integrality therefore gives q|2. The single dyadic prime,
valuation0/1/2 and the preceding unit-square result give precisely the
necessary representatives1,d,2. Each is totally positive and divides2,
so every integral vector of one of those norms is reflective. No
high-dimensional field-list theorem is used.

The exact trace formula U(A,B)+2U(C,D) is correct. Its nonzero minimum4
occurs only at +/-1. Trace8 for a spatial norm2 vector therefore gives
only two +/-1 entries or one +/-s, the latter imprimitive. The same
trace possibilities cannot have norm d. Thus the initial B9 is the FULL
zero germ. Admissibility at exactly one real embedding and local
finiteness via all four embeddings are valid.

## Endpoint inventory and all embedding bounds

The upper r=(1;phi,0^8) is primitive, has norm d, has the correct
orientation, and meets before light at a finite point of square
-a*d/phi^2. A finite simple continuation of the retained whole B8 is
isolated or the B9 long-end3 extension. In the latter row the retained
coordinate must satisfy z^2=q/2. Norm1 requires z=s/2, nonintegral;
norm d is also impossible: d/2 in F=Q(s) has rational norm1/2, so is
neither a square in F nor5 times a square in F, and hence is not a
square in F(sqrt5). Norm2 gives exactly z1. Thus the four weighted
profiles and residuals1,2,d,1 are complete.

The collision inequality y>=phi*x gives d*x^2<=S. At the other three
embeddings the form is positive definite and gives exactly the source's
bounds. Adding them yields S/d+S/a-Tr(S/a). The traces-6 and-8 give
the stated strict trace cutoffs; integral trace is even, so the
cutoffs8,16,8 follow without losing equality cases.

For trace<=8 the positive list1,phi,u,s is exhaustive. For S1, phi/s
fail physically and1/u fail embedding3. For Sd the physical survivors
are1,u, of which u fails embedding3 and1 satisfies every bound.

For S2, I independently hand-enumerated U<=16: its values are
0,4,6,14,16, and2U<=16 permits0,8,12. All signs in the printed z+s*w
families are covered. Physical plus embedding3 filtering leaves exactly
phi, phi-s, s-u, s*phi-1. The borderline rejected1+s*u is excluded
by the displayed strictly positive identity. Embedding1 excludes phi
because (s-1)*phi+s>2, using s>7/5 and phi>8/5. It sends the other
three to absolute values phi+s,a,s*phi+1, each>2, whereas its permitted
absolute value is<sqrt(8/3)<2. Ignoring embedding2 in this necessary
filter only enlarges the tested set. Therefore S2 has no solution.

Only Sd,x1 remains, and y^2=phi^2 with y>0 gives r itself. No earlier
or simultaneous simple alternative survives. The endpoint's FULL
actual germ is B8+A1(r); no preselected cap filter or finite-height
exhaustion is substituted for firstness.

## Scope and remaining frontier

The form is genuinely admissible biquadratic, not the closed cyclic
quartic host. Scalar descent to a proper subfield cannot have a single
indefinite real extension: an indefinite descended form remains
indefinite at every extension, regardless of scalar sign. The natural
binary orthogonal columns have determinant d, a nonunit, so the old
silver exact-splitting argument does not apply. These checks justify
the stated distinction, without a claim of worldwide novelty or an
exhaustive literature classification.

No mathematical gap found. Only arithmetic preflight, the full initial
germ and one complete actual finite edge are proved. The literal
quartic host remains OPEN: no compact polytope, full finite chamber,
cofinite reflection group, nonreflectivity or uniform bound follows.
Keep it within the existing construction frontier; no new manuscript,
catalogue entry, worker or formalization is warranted by this step.
