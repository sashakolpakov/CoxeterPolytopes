# The finite-quadrilateral transfer needs only small affine orders through9

**Final internally reviewed source audit.** In the actual four-finite-
corner quadrilateral star, assume all four finite Vinberg payments vanish,
both complements H−U and H−V are affine-free, and no connected affine
subdiagram of the ambient chamber has order at most9. Then the remaining
five cases of Vinberg's finite-quadrilateral proof give a contradiction.

This lowers a conservative cutoff in the frozen
[transfer audit](fresh_simple_finite_quad_vinberg_transfer_20261017.md).
That audit bounded entire carriers by13,16 or17; their size is not the
size of an affine diagram that must be forbidden. The preliminary
complement assumption here is essential and is supplied separately by
the [order10 carrier proof](fresh_affine_carriers_order10_20261017.md).
This note does not directly apply a compact-ambient theorem to a
noncompact chamber.

The primary input is Vinberg's original §6.5–11, printed pp.82–91:
[original paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).
Its finite-diagram case classifications and determinant tables remain
published inputs. The present audit identifies the affine-sensitive
steps and supplies the needed size or geometric reason. It does not
claim a new independent enumeration of those tables.

## Preliminary hypotheses retained

All T_ij are actual finite corner components. The original finite
Lemmas6–11, including special vertices, end segments and the four E8
marked positions, consequently retain their stated hypotheses. The
affine-freeness of H−U and H−V implies that of U,V. Thus the original
orthogonal-pair reduction and cyclic Lanner-core lemma remain valid.
Any such core has at most five nodes.

The actual-center orthogonality argument is unchanged: a block right
to an indefinite block is positive definite, regardless of its size.
Its span lies in a positive definite space; a putative positive affine
kernel cannot represent zero because its product with a strict center
is negative. Actual corank one, the alternating black dependence and
the nonzero determinant-polynomial coefficient are also unchanged.

## The five cases and their actual affine obligations

* **§6.7, printed pp.85–86.** The explicitly invoked affine diagrams are
  E8 and the smaller diagrams occurring beside its displayed exceptional
  patterns. They have at most9 nodes. The other contradictions are a
  high-bond indefinite extension orthogonal to V, failure of the already
  required complement indefiniteness, or a displayed negative-index-two
  diagram. None requires forbidding a whole13-node carrier.

* **§6.8, pp.86–87.** The named affine obstructions are B4 and F4, each
  order5. The two end pieces orthogonal to V are positive definite by
  the actual-center argument. Their size is irrelevant to that argument.

* **§6.9, pp.87–89.** The explicit exceptional affine obstructions have
  types among E6, A2, A4 and G2, all smaller than9. The exceptional
  E8-containing configurations are otherwise eliminated by the
  published negative-index or nonzero-determinant tests. In the long
  branch, O(U),O(V) are connected: U,V are connected by Lemma7 and every
  newly adjoined node touches one of their black roots. Each contains
  the original finite high bond of label at least5. The only connected
  affine diagram containing such a finite bond is affine G2 of order3.
  Thus neither O(U) nor O(V) can be affine under the order-nine
  prohibition; their orthogonality forces an elliptic one. Subsequent
  connected high-bond pieces of order at least5 are indefinite by the
  finite/affine classification, and their orthogonal complements are
  positive definite by the actual-center argument. There is no residual
  order10 affine exclusion hidden in O(U),O(V).

* **§6.10, pp.89–90.** The high-bond polynomial and the actual-star
  nonzero-coefficient argument are those of the frozen audit. The two
  displayed U∪T_ij∪L exceptional patterns are rejected by the published
  negative-index-two configurations, not by declaring their whole
  at-most17-node union affine-free. The other branch constructs an
  indefinite high-bond extension orthogonal to the cyclic Lanner core.
  The only possible affine high-bond exception is again G2, already
  excluded. The finite-corner classification controls the complete
  shared branches in these arguments.

* **§6.11, pp.90–91.** The affine-sensitive overlap and connecting-edge
  tests lie in the induced union L∪M of two cyclic Lanner cores, with
  |L|,|M|<=5. Any affine subdiagram of that union must omit a node of L,
  because L is indefinite. Therefore its order is at most9. If the cores
  are disjoint, it must omit one node of each, giving the stronger bound8.
  Additional edges between core nodes do not change this count. The
  later determination of the entire shared diagram uses finite T_ij and
  Lemma8, and the final end-piece contradictions use ellipticity of U,V;
  they are not affine tests on arbitrarily enlarged shared carriers.

The distinction between a core union and a whole extended carrier is
material. In particular the bound17 in the earlier audit was a safe
envelope, not evidence of an affine obstruction of order17.

## Consequence and limits

All actual affine exclusions after the preliminary complement reduction
are covered by the order-nine prohibition. Together with the separately
proved order10 carrier exclusion, every quadrilateral with four finite
corners receives positive finite Vinberg payment under that prohibition.
The underlying chamber need not have simple ideal vertices.

This is a refinement of the source-transfer hypothesis. The published
finite casework remains an explicit dependency. It says nothing by
itself about quadrilaterals containing an ideal vertex, the cost of
weights at a nonsimple cusp, or a resulting ambient dimension ceiling.
A separate root mathematical/source read passes this cutoff reduction
after a fresh reading of the primary PDF text on printed pp.85–91.
Fresh screenshot rendering failed; neither this audit nor that review
claims a fresh visual re-enumeration of the displayed diagrams or an
independent rederivation of the published determinant tables. The
published casework remains explicit input. Local links pass. This note
is final and frozen, with internal review rather than external refereeing
or proof-assistant formalization.
