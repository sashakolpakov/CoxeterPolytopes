# Independent review of the abelian-base degree and asymmetric arc escape

**PASS:** the chain parameter retains degree50 over `K Q^ab`. Every coherent
added Coxeter wall retaining the133-wall H20 chamber and having finite or
parallel contacts to any three consecutive chain walls preserves the same
weak spacelike escape. Any finite coherent collection of such caps therefore
leaves infinite volume. The contact angles may be unequal and all18 support
contacts are unrestricted Coxeter contacts.

The [independent checker](../scripts/fresh_h20_asymmetric_arc_escape_gate_20261016_review.py)
and [certificate](fresh_h20_asymmetric_arc_escape_gate_20261016_review.json)
import no producer arithmetic. They reconstruct the triple Gram inverse and
norm bound symbolically, invert the support Gram rationally, enumerate the
36 finite support profiles in reverse slot order, and check every continuous
endpoint margin. The [primary theorem](fresh_h20_asymmetric_arc_escape_gate_20261016.md)
is a dependency for literal comparison. The independently reviewed modular
factorizations and physical chamber remain explicit inherited dependencies;
this review does not rerun those earlier computations.

## The degree statement survives abelian base change

The inherited modular review gives irreducibility of the degree50 polynomial
over the literal field K, and a squarefree Frobenius factorization13+37.
Its simple tower roots and unit leading coefficient justify that Frobenius
cycle type. The new review independently forms the13th power of a differently
labelled representative and finds a37-cycle with13 fixed points.

Any proper block system for a transitive degree50 permutation group has
block size2,5,10 or25 and number of blocks25,10,5 or2. The order37 element
cannot permute the blocks nontrivially and cannot act nontrivially inside
any block. Its existence therefore proves primitivity. A faithful transitive
abelian group is regular, whereas this element has fixed points. Thus the
group is nonabelian. Its nontrivial commutator subgroup is normal; its orbits
form blocks, so primitivity makes that subgroup transitive.

If L is the splitting field, `L intersect K Q^ab` is abelian Galois over K.
Its fixing subgroup consequently contains the commutator and is transitive
on all50 roots. Base change identifies this subgroup with the Galois group
of `L K Q^ab / K Q^ab`. The polynomial remains irreducible of degree50 there.
In particular, no nonzero polynomial of degree below50 over Q^ab vanishes
at C. No claim that the degree over Q^ab itself is exactly50 is required.

## Arbitrary consecutive contacts and the norm budget

The independent symbolic inversion gives

`y^2 = A^2/2 - ((B+D)^2+4CBD)/(8C(C+1))`

for consecutive contacts `-B,-A,-D`. The determinant is `-32C(C+1)`.
If all support contacts were nondotted, the norm equation would yield a
degree-at-most-two polynomial over Q^ab with constant `-(B+D)^2`.
For `B+D>0` this contradicts the degree statement.

The exceptional outer contacts `B=D=0` give `y=-A*n_middle/2`.
For A>0 the next available chain wall has positive pairing `(2C+1)A`;
the checker verifies its availability at all24 triple positions. For A=0,
the cap lies in the negative support-dual cone and pairs positively with
every strict old-sector point, contradicting coherence. Thus some support
contact must be dotted. The inverse support Gram J is entrywise nonnegative
with minimum diagonal1/2, so its cost Q exceeds2 and `u=-y^2=Q-2>0`.

The outer chain sum has square `-8C` and is future timelike; its pairing
with y is `-(B+D)<0`, fixing the future sheet of y. Independent differentiation
of the norm formula shows u increases in B,D and decreases in A on their
nonnegative domains. Evaluating at `A=0,B=D=2` gives `u<=2/C`. The certified
physical interval consequently gives `Q<281/70`.

## Weak feasibility and all support profiles

The reviewer verifies `E0=x-e/100-omega5`, with square `190011/110000>0`.
The difference between E0 and the independently verified strict escape
Eeps has nonnegative coefficients in the old support normals, checked by
an independent multiplication by J. Acuteness therefore implies E0 is
strictly feasible off the support, weakly feasible on it, and has negative
pairing with the saved future center. On support root i its pairing is
`-delta(i,5)`.

A dotted A10 contact can occur only at endpoint1 or10 and leaves too little
budget for any other nonzero support contact. The two exact affine-in-u
margin tests cover these continuous rays. Otherwise the required dotted
contact is outside A10 and leaves A10 cost below141/70. Diagonal bounds
exclude all but endpoints1/10 with labels3/4 and second endpoints2/9 with
label3. The36 possibilities leave exactly eight profiles. All28 rejected
costs and all eight retained costs/pairings agree with the primary.

For each nonzero retained profile, `K*q_A-z_A^2>0`, where
`K=9999/10000`, and the dotted outside contact gives `u>q_A`. Reverse Lorentz
Cauchy--Schwarz yields `r.E0<=-sqrt(K*u)+z_A<0`. The zero profile has z_A=0
and needs only u>0. Unlike the earlier strict-escape argument, no uniform
positive gap for u is needed. Every coefficient outside A10 contributes
exactly zero to the E0 support pairing.

## From a weak escape to infinite volume

For any finite coherent collection of these caps, choose its strict future
timelike interior H. Every cap pairs strictly negatively with E0, and every
old wall pairs nonpositively. For sufficiently small positive delta,
`E0+delta*H` remains spacelike and pairs strictly negatively with all walls.
Its pairing with the saved old future center H0 remains negative, because
both E0.H0 and H.H0 are negative. The segment from H to this spacelike point
therefore crosses the future null cone while all inequalities remain strict.
An open ideal-boundary patch is feasible, proving infinite volume.

This reasoning does not require `E0.H<0` for an arbitrary new interior H,
nor a common quantitative delta for all possible collections. The collection
is finite and the chosen perturbation is allowed to depend on it.

Generation and a fresh saved replay pass. This is a construction-family
obstruction, not general H20 nonexistence or a new dimension bound. Caps
without consecutive nondotted triples and changes to retained walls remain
outside its scope.

```sh
python3 scripts/fresh_h20_asymmetric_arc_escape_gate_20261016_review.py --verify
```
