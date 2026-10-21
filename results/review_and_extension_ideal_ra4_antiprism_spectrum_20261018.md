# Excluding every antiprism from an enlarged octahedral facet spectrum

18 October 2026. **New written theorem, awaiting root's final review.**
The square-antiprism report and its TeX appendix pass this independent full
mathematical/source read: I found no flaw. The extension below permits all
antiprism sizes simultaneously and proves their exclusion. No numerical
census, optimizer, checker, Lean process, TeX build or frozen-source edit
was performed. This report is not a dimension exclusion or priority claim.

Reviewed sources and SHA256:

* `results/ideal_ra4_square_antiprism_exclusion_extension_research_20261018.md`:
  `c2c223789ee45e7aebea5f2d8ddf794ca979867c1647e45455775a48bc185ba2`.
* `manuscripts/ideal_ra4_square_antiprism_extension_appendix.tex`:
  `023e85d3f4d9d7cda226bd234d0eeae8226c5ea9997b4c69f08b1d306f5f7ce3`.

## 1. Review and precise generalized statement

Write O for the regular ideal right-angled octahedron, D for its double
across a triangular face, and A_k for the combinatorial k-antiprism.
Its two bases are k-gons and its remaining 2k faces are triangles.
The reviewed A4 proof correctly makes the square-base cusp side **long**,
not short. Its six box shapes are exhaustive, and its last contradiction
uses the same actual neighboring facet at all three vertices of one
triangular ridge. The elementary normal argument proves the asserted
intrinsic rigidity rather than assuming that combinatorics fixes geometry.
The ideal5 carrier identity in the appendix has the correct multiplicities.

**Theorem.** Let P be an entirely ideal, finite-volume, right-angled
hyperbolic four-polytope. Suppose every facet is isometric to O or D,
or is combinatorially A_k for some k>=4. The permitted antiprism sizes
need not be fixed in advance. Then no antiprism facet occurs. Consequently
P lies in the original O/D class and satisfies its existing chain theorem.

The proof uses the binary distinction **triangular/nontriangular** on a
cusp side. For k>4 its distinguished base side is k-gonal, not quadrilateral.
Retaining the word quadrilateral in that part of the A4 proof would be an
invalid generalization.

## 2. Own rigidity and cusp-side geometry for every A_k

Use Lorentz signature (-,+,+,+) and outward unit normals. We use the same
acute-wall fact as the inherited geometric proofs: distinct wall products
are nonpositive. The following projection spells out the further gap
needed for the two bases; disjoint combinatorial bases alone do not assert
that their complete hyperplanes are disjoint.

For wall product g in (-1,0), the two-wall Gram matrix

    G = [[1,g],[g,1]]

is positive definite with entrywise nonnegative inverse. If p is a strict
future timelike interior vector, v=(p.b,p.c) has negative entries and

    q = p - (b,c) G^(-1) v

adds positive multiples of b and c. Every other inequality stays strict,
q^2=p^2-v^T G^(-1)v<0, and p.q=q^2<0 gives the future sign. This is an
actual codimension-two intersection with a nonright angle, impossible.
Product zero likewise gives an actual ridge. Neither possibility is
compatible with the disjoint bases of an antiprism.

If b.c=-1, b+c is nonzero null: b=-c is excluded by a strict interior
point. Also p.(b+c)<0, so it is future, and its products with all other
walls are nonpositive. It is therefore a feasible ideal point in both
bases. In the bounded Klein realization an ideal point is extreme, hence
is an actual shared ideal vertex. The bases have disjoint vertex sets,
so this is impossible. Thus b.c<-1 and the base planes are ultraparallel.

Choose coordinates

    b=(0,0,0,1),  c=(sinh L,0,0,-cosh L),  L>0.

The upper cusps, in cyclic order, are

    xi_j=(1,cos theta_j,sin theta_j,0).

The lower lateral triangle at xi_j is opposite b in the rectangular link;
its unit normal is ell_j=lambda_j xi_j-b, lambda_j>0. Orthogonality to c
gives lambda_j=coth L=:a>1. Consecutive ell_j,ell_(j+1) are opposite at
the lower cusp between them, so

    -1 = ell_j.ell_(j+1) = a^2(cos(theta_(j+1)-theta_j)-1)+1.

Each positive cyclic gap is alpha or 2pi-alpha for one alpha in (0,pi).
A large gap makes their sum at least 2pi+(k-2)alpha>2pi. Hence all k gaps
are 2pi/k. With beta=pi/k we obtain

    a=csc beta,  coth L=csc beta,
    sinh L=tan beta,  cosh L=sec beta.

This fixes both bases and all lower triangles. An upper triangle plane
is uniquely fixed by its upper base edge and orthogonality to b. All
2k+2 supporting planes are fixed, proving rigidity up to isometry.
There is no appeal to a complete polyhedron classification or general
three-dimensional rigidity theorem in this calculation.

For explicit width and triangle signatures, put theta_j=2j beta. The
unit upper-triangle normal and the lower vertex are

    u_j=(cot beta, csc beta cos(theta_j+beta),
                  csc beta sin(theta_j+beta), 0),
    w_j=(1, cos beta cos(theta_j+beta),
            cos beta sin(theta_j+beta), -sin beta).

Here u_j bounds the triangle (xi_j,xi_(j+1),w_j), and w_j is null.
At xi_j the opposite unit-normal sums are

    b+ell_j = csc beta xi_j,
    u_(j-1)+u_j = 2cot beta xi_j.

Their coefficients are the corresponding horospherical widths. The
side on b spans the larger width, so the unique nontriangular side is
long. The ratio is

    rho_k = 2cot beta / csc beta = 2cos(pi/k).

At the lower cusp the two opposite sums are

    c+u_j = (sin beta cos beta)^(-1) w_j,
    ell_j+ell_(j+1) = 2csc beta w_j.

Their width ratio is the same rho_k, and u_j spans the larger width.
Thus the very same triangular face is a short side at xi_j and
xi_(j+1), and a long side at w_j. In particular, for every k>=4:

    A_k: aspect rho_k:1, exactly one long nontriangular side;
    every lateral triangle: short, short, long.

The ratios rho_k are strictly increasing, lie in [sqrt(2),2), and never
equal 2. The O/D signatures are unchanged: square O cusps have no
nontriangular side; square D cusps have one; nonsquare D cusps have
aspect 2:1 with both long sides nontriangular and both short sides
triangular.

## 3. Complete box argument with all antiprism sizes present

At a cusp of P, order the box lengths x>=y>=z and scale z=1. Every
unequal pair ratio lies in

    R = {2} union {rho_k : k>=4 and A_k is permitted}.

All elements lie in [sqrt(2),2]. If x>y>z, both x/y and y/z are at least
sqrt(2), while x/z is at most 2. Equality is forced in all three bounds.
The only possible three-distinct-length shape is therefore 2:sqrt(2):1.
Every other noncubic box has shape rho:1:1 or rho:rho:1 for rho in R.
This argument covers even an unbounded allowed set of k; no label cutoff
or finite census is hidden in the box list.

The shape rho:rho:1 is impossible. Its four nonsquare belt rectangles
are all A_k when rho=rho_k, or all D when rho=2. Each requires at least
one long nontriangular side, lying on one of the two square ends. This
demands at least four nontriangular edges on those ends. A square end
can only be O or a square cusp of D and permits at most one, giving
total capacity at most two. Edge labels must agree on incident faces.

In 2:sqrt(2):1, the two 2-by-1 faces are D. Their long sides force
both long edges of each 2-by-sqrt(2) A4 face to be nontriangular.
An A4 cusp permits only one, a contradiction.

Consequently every noncubic cusp has shape rho:1:1. Its four nonsquare
belt rectangles have triangular short sides, so the two square ends
have no nontriangular side and must be O. If rho=rho_k the entire belt
consists of A_k facets: strict monotonicity excludes a different
antiprism size, and D has ratio 2. A short side of an A_k rectangle
thus meets O, whereas a long side meets another A_k.

Now choose an actual lateral triangular ridge T of an A_k facet F.
There is a unique other facet F' across T. At T's two upper cusps its
side is short, forcing F' to be O. At its lower cusp that side is long,
forcing the same F' to be A_k. These distinct facet types cannot
coincide. This excludes F and proves the theorem.

## 4. Universal necessary three-face spectrum in ideal dimensions five/six

For an ideal right-angled four-polytope Q define

    kappa(Q)=3 f_0(Q)-4 f_3(Q)+24.

The inherited O/D theorem gives for an m-cell chain
f_0=18m+6 and f_3=14m+10, hence kappa=2-2m<=0.

For completeness the universal ideal5 identity needs no classification.
Write v,e,r,g,f for the numbers of faces in dimensions 0,...,4.
The cusp four-box gives e=8v, 24 three-faces per cusp and eight
four-facets per cusp. Since every ideal right-angled three-polyhedron
has f_2=f_0+2 and each two-face belongs to three three-faces,

    3r=24v+2g.

Euler's identity v-e+r-g+f=2 gives g=3v+3f-6. Each three-face belongs
to two four-facets, so

    sum_(four-facets Q) kappa(Q) = 3(8v)-4(2g)+24f = 48.

Some Q therefore has positive kappa. If all its three-faces were O,
D or antiprisms of any sizes, the theorem would put Q in the O/D
chain class, contradicting positive kappa. Thus:

**Corollary.** Every entirely ideal finite-volume right-angled
five-polytope has a positive-kappa four-facet containing a three-face
which is neither O, D, nor any antiprism. Every ideal right-angled
six-polytope has such a three-face; indeed every one of its
five-facets contains one.

The same necessary spectrum applies to every entirely ideal five-face
inside a higher-dimensional or mixed-vertex right-angled polytope.
This gives no contradiction for ideal dimensions five or six by itself.

Conditionally on the complete classification of ideal right-angled
three-polyhedra with at most nine vertices, the required three-face has
at least ten vertices. The usual list is O at six vertices, A4 at eight
and D at nine; defect one rules out seven. This numerical translation
uses the completeness statement, not the present own normal argument.
The ten-vertex A5 is also excluded, but the Figure 7 octahedrite is not.
With the separate complete ten-vertex classification, the first remaining
possibility is therefore that octahedrite or a face with at least eleven
vertices. No realization or exclusion of either branch is asserted.

## 5. Explicit remaining blocker and verification scope

The Figure 7 octahedrite has four cusps where its two quadrilateral
sides are adjacent, two where they are opposite, and four with one
quadrilateral side, as read from the existing finite incidence
transcription. The adjacent pair is already outside the antiprism/D
signature used above. Neither its actual cusp ratios nor the lateral
triangle signatures have been established in this report. Adding it
would invalidate the restricted box-ratio hypothesis unless those facts
are first proved. A formal coloring or the abstract matching of its
quadrilaterals cannot replace that geometry.

This report proves no exclusion for the Oct/D/antiprism/Figure-7 class.
It does not extend the old surface proof to new retained facet types:
the new theorem removes all added types before applying that proof.
The published minimum, general small-facet classifications and the
external-priority question are not used to prove the generalized
antiprism exclusion. The inherited acute-wall property, rectangular
ideal links, actual-neighbor incidence and the previously proved O/D
chain theorem remain explicit geometric dependencies. All new algebra
here is written mathematics only, with no machine or kernel replay.
