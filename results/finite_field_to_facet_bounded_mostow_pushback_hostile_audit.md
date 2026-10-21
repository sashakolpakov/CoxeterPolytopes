# Hostile audit: finite-field/facet/Mostow pushback

**Verdict: PASS after four local precision repairs.**

Audited owner snapshot:

    406552c48a439d5d1d1d32b698f79acaf5483610e5b3052984ecd73afdd4be68

Repaired owner snapshot:

    5da9446cf111dba05ff5f01b15e2dedd022ec780797f0b5b441b3af90945f2e2

The conclusion is sound for normalized geometric Coxeter chambers and
side-reflection groups in dimensions \(n\ge3\).  No reflective-lattice
finiteness theorem is imported.  Dimension two and arbitrary
primitive-root/Cartan rescalings remain explicitly outside the conclusion.

## 1. Finite angle alphabet: PASS

For \(m\ge3\),

\[
 c_m=4\cos^2(\pi/m)=2+\zeta_m+\zeta_m^{-1},
 \qquad
 \mathbf Q(c_m)=\mathbf Q(\zeta_m+\zeta_m^{-1}),
\]

and the latter field has degree \(\varphi(m)/2\).  The only exceptional
finite label for the displayed degree formula is \(m=2\), which the owner
separates and which contributes the rational value \(c_2=0\).

The length-two Cartan cycle \(c_m\) lies in the Vinberg ground field.
Therefore \([k(P):\mathbf Q]\le D\) implies
\(\varphi(m)\le2D\).  Since \(\varphi(m)\to\infty\), the finite labels form
a finite alphabet.  One additional symbol \(\infty\) legitimately records
both parallel and ultraparallel pairs: it specifies the abstract relation,
not the distance.  The complete distance is recovered only at the later
rigidity step.

In fact a field-degree bound is enough for this step; a finite list of
fields is stronger than necessary.

## 2. Fixed presentation and Mostow--Prasad: PASS after repair

For a genuine finite-volume Coxeter chamber, the Coxeter fundamental-domain
theorem gives a faithful side-reflection representation of the abstract
group \(W_M\).  Given two such complete finite-covolume representations:

1. choose a torsion-free finite-index normal subgroup \(H\);
2. ensure \(H\) lies in the orientation kernels of both representations;
3. apply finite-volume Mostow--Prasad rigidity to the prescribed
   isomorphism between the two images of \(H\); and
4. extend conjugacy from \(H\) to all of \(W_M\), because \(H\) is normal
   and the centralizer of a finite-covolume hyperbolic lattice in the
   effective projective isometry group is trivial.

This proves uniqueness of the **marked homomorphism** up to ambient
conjugacy.  It applies to nonuniform lattices in \(\mathbf H^3\): nearby
incomplete Dehn-surgery representations are irrelevant because the theorem
compares only faithful, discrete, complete finite-covolume representations.
The owner's warning against replacing this by a finite-type
character-variety argument is therefore necessary and correct.

The original proof used the abstract parity subgroup while stating the
theorem for arbitrary faithful lattice representations.  Parity is known
to equal geometric orientation only for reflection realizations.  The
owner was repaired to use the common kernel of the two actual orientation
characters; it then notes that parity suffices in the reflection case.

Once the marked reflections are fixed, their labelled mirrors are fixed.
There are at most \(2^f\) choices of halfspace orientations, only some of
which can be nonempty irredundant Coxeter chambers.  Hence:

- one marked conjugacy class of reflection homomorphisms;
- finitely many marked oriented chambers; and
- still only finitely many objects after forgetting labels or orientations.

This is the correct marked/unmarked quantifier.  The normalized object here
is the outward unit-normal geometric Gram/chamber.  Mostow does not make
arbitrary asymmetric Cartan entries or primitive-root scalings finite.

## 3. Jordan--Zassenhaus/Burnside closure: PASS at its stated scope

Fixing the \(k\)-rational representation \(\rho\) and one invariant full
lattice \(L_0\), the ring

\[
 \Lambda=\mathcal O_k[\rho(W_M)]
 \subseteq\operatorname{End}_{\mathcal O_k}(L_0)
\]

is a finite \(\mathcal O_k\)-module and its \(k\)-span is the representation
algebra, so it is an order.  Jordan--Zassenhaus gives finitely many
isomorphism classes of \(\Lambda\)-lattices with rational module \(V\).
Zariski density of the finite-covolume lattice makes the standard module
absolutely irreducible; Burnside after algebraic closure and descent gives

\[
 k[\rho(W_M)]=\operatorname{End}_k(V).
\]

Thus a module isomorphism between two invariant full lattices extends to a
scalar, and the isomorphism classes become homothety classes.  This proves
exactly the stated fixed-representation lattice result.

It does **not** prove finiteness of ambient quadratic forms before a
rational representation is fixed, nor finiteness of unit-rescaled
primitive roots.  The repaired owner now says this explicitly.

## 4. Countercontrols: PASS

For the frozen \(\mathbf H^{10}\) redoubling family,

\[
 q_{2m}=c+2m^2a+2m(m-1)b,\qquad
 u(c,q_{2m})=(m^2-1)^2.
\]

Taking \(m=p+1\) gives
\(u=p^2(p+2)^2\), while \(u+1\equiv1\pmod p\); hence
\(v_p(u/(u+1))=2\) for every odd prime \(p\).  This is the claimed exact
height/bad-prime escape.  The chambers are proper finite-index reflection
descendants of a seed, so they test all the soft bounded inputs but do not
contradict a theorem restricted to absolutely maximal reflection groups.

For

\[
 q_p=\langle-p\sqrt5,1,\ldots,1\rangle
\]

in rank twelve, the selected real place has signature \((11,1)\), the
other place is positive definite, and therefore \(q_p\) is \(k\)-anisotropic.
The associated arithmetic lattices are cocompact.  In even rank the
determinant squareclass is invariant under similarity.  For distinct
primes \(p,q\ne5\), the ratio \(p/q\) has odd valuation at a prime over
\(p\), so the determinant squareclasses, and hence the arithmetic
commensurability classes, are distinct.  These hosts are correctly not
claimed to be reflection generated.

## 5. MAX-FAC sufficiency: PASS

After dimensions and Vinberg-field degrees are bounded, a facet bound
\(f\le F\) gives finitely many labelled Coxeter matrices.  For each fixed
\((n,M)\), the rigidity argument gives one marked reflection homomorphism
and finitely many oriented-side chambers.  Taking a finite union proves
finiteness of normalized chambers and conjugacy classes of reflection
groups.  Dotted heights, normalized bad-prime sets and covolumes then become
bounded by finite enumeration; no separate theorem for them is required.

For an arithmetic finite-covolume reflection lattice, maximality among
arithmetic reflection overgroups is equivalent to absolute maximality among
discrete reflection overgroups: any discrete overgroup is again a lattice,
has finite index, and is arithmetic by commensurability.  This equivalence
was added to the owner.

Thus MAX-FAC is genuinely sufficient once the field/degree and dimension
steps are available.  It remains unproved, and Mostow rigidity itself gives
no presentation or facet bound.  In dimension two the reduction is false
because a fixed polygonal Coxeter presentation can move in Teichmuller
families; the owner correctly leaves that as a separate branch.

## 6. Repair ledger and non-import audit

The repairs were:

1. replace parity by the common actual orientation kernel in the broad
   rigidity theorem;
2. state separately the marked homomorphism, oriented-side and unmarked
   quantifiers;
3. make the order/Burnside step and its fixed-rational-representation scope
   explicit; and
4. state why arithmetic-reflection maximality equals absolute reflection
   maximality here.

The proof uses cyclotomic degrees, the Coxeter/Poincare
fundamental-domain theorem, Selberg's lemma, Mostow--Prasad rigidity, the
trivial-centralizer consequence of density, and Jordan--Zassenhaus plus
Burnside.  None of Nikulin reflective-lattice finiteness, ABSW, or any
equivalent pre-existing finiteness theorem appears.

The repaired Markdown passes a full Pandoc parse.
