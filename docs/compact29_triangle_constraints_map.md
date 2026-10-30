# Compact29 triangle constraint: proof and reproduction map

For current reading, see the [current dimension-bounds manuscript](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

The supporting argument
proves a local incidence constraint: at a compact H29 vertex of type
E7³+E8, at most three of nine specified two-faces are triangles. At least
six have four or more sides. It does not lower the compact bound29,
establish a new face-charge saving, or construct a polytope.

The [full proof report](../results/fresh_compact29_three_triangle_star_constraint_20261018.md)
and [exact certificate](../results/fresh_compact29_three_triangle_star_constraint_20261018.json)
record the finite domain and reviews. The manuscript is an exposition of
that argument, with no new enumeration.

## Standard inputs and newly proved steps

|Layer|Statement|Where it is used|
|---|---|---|
|Standard geometric input|Compact acute Coxeter polytopes are simple, their finite vertex roots are independent and positive definite, and wall products are native Coxeter products.|Actual two-face corners and the ambient positive-index contradiction.|
|Standard classification input|Complete finite Coxeter diagrams; E8 has no connected finite extension; E7 has only E8 as a connected finite one-root extension; labels≥6 occur only in isolated dihedral factors.|Untouched-component exclusion and the all-label reduction.|
|Internal geometry, Lemma2|Any contact with an untouched E7 forces a rank30 positive Gram matrix; an untouched E8 forbids contact directly.|Confinement of the actual third side to the marked E7.|
|Internal geometry, Lemma3|The three nonadjacent marked pairs rule out every high finite contact.|The finite domain is all 3×4⁷ vectors, not a cutoff.|
|Internal geometry, Lemma4|A PSD-singular principal actual-wall diagram yields a future ideal point in the closed acute chamber.|Compactness exclusion; no assumption that the principal diagram is a known face.|
|Finite exact certificate, Proposition5|All 49,152 vectors give 26 finite-eligible vectors, nine nonpositive heights, fourteen affine obstructions, and three identical contact vectors.|The forced positive height (1+3√5)/4.|
|Internal geometry, Theorem1|The common finite vertex fixes both positive height signs; cross-factor projections are orthogonal.|The mutual product is forced, rather than freely choosing Lorentzian sheets.|
|Finite exact certificate|The forced square lies strictly between native labels11 and12.|At most one E7 factor contributes any designated triangle.|

The general finite classification and compact simplicity are standard
mathematical inputs. The positive-rank contradiction, all-label argument,
compact affine exclusion, and shared-height coupling are internal proof
steps and should not be hidden as external axioms in a formal translation.

## Exact implementation interfaces

- [Producer](../scripts/fresh_compact29_three_triangle_star_constraint_20261018.py):
  independently enumerates each complete label vector; uses the exact
  finite-diagram classifier and field arithmetic from the frozen local
  atlas; saves every eligible vector, explicit affine exclusions, domain
  digests, and exact gap enclosures.
- [Standalone checker](../scripts/check_compact29_triangle_star_saved_20261018.py):
  Python standard library only, with no project imports; reconstructs a
  rooted-tree finite catalogue and rational E7 inverse, checks the entire
  domain and all heights, uses reverse Schur PSD tests, and checks the
  native gap with a different Taylor truncation.
- [Frozen local atlas](../results/fresh_compact29_small_faces_20261015.md):
  contains the three surviving triangle recipes. Sufficiency of those
  recipes is inherited and is not needed for the new necessary exclusion.
- [Earlier fixed-root packet](../results/fresh_compact29_coupled_first_roots_20261018.md):
  supplies the exact gap arithmetic interface. The new packet adds the
  complete actual triangle domain, which the earlier fixed-recipe result
  did not claim.

The producer directly pins its own source and seven inherited files.
The saved checker verifies those bindings. Hash checks establish artifact
identity; they do not establish the geometric lemmas.

## Reproduction

From the repository root, without changing frozen output files:

```sh
python3 scripts/fresh_compact29_three_triangle_star_constraint_20261018.py --verify
python3 scripts/check_compact29_triangle_star_saved_20261018.py
```

Both commands have terminated successfully: producer generation
7924/ef0e85, independent complete replay57352/407f35, and independent
checker24256/ee6bb9. The report records the exact completed review scopes.

Build the manuscript into a separate directory:


The next unresolved task is a global charge use of the incidence
restriction. In particular, the six nontriangles may be quadrilaterals;
their existence alone does not justify deleting six existing payments.
