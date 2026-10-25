import Std

/-!
SOURCE PREPARATION ONLY: uncompiled; no actual axiom audit is credited.

Conditional discrete composition for the frozen 470-line compact17 candidate.
All predicates refer to ONE fixed candidate polytope and the SAME original
vertex. Actual geometry, analytic estimates, accepted classifications and
source-bound mathematical inputs are ordinary theorem PARAMETERS below.
There is no supplied whole-D18 or final dimension<=17 impossibility premise.
-/

namespace CoxeterProofs.CompactDimensionSeventeenComposition

inductive HeadKind where
  | finite (nativeLabel : Nat)
  | parallel
  | dotted

/-- The complete retained actual fork-role cover before the new exclusion. -/
inductive ForkRole where
  | isolated
  | nearThree
  | farThree
  | farNextThree

/-- Six unordered pairs, and the separate far-next-containing domain. -/
inductive PairCase where
  | isolatedIsolated
  | isolatedNear
  | isolatedFar
  | nearNear
  | nearFar
  | farFar
  | containsFarNext

def pairCase : ForkRole → ForkRole → PairCase
  | .farNextThree, _ => .containsFarNext
  | _, .farNextThree => .containsFarNext
  | .isolated, .isolated => .isolatedIsolated
  | .isolated, .nearThree => .isolatedNear
  | .nearThree, .isolated => .isolatedNear
  | .isolated, .farThree => .isolatedFar
  | .farThree, .isolated => .isolatedFar
  | .nearThree, .nearThree => .nearNear
  | .nearThree, .farThree => .nearFar
  | .farThree, .nearThree => .nearFar
  | .farThree, .farThree => .farFar

universe u

/--
Named external WRITTEN interfaces. The complete head/role covers are explicit
premises, not declarations of realized diagrams. The new analytic/reference,
local-A18 and actual edge arguments enter only through DNext and DNearFar;
all their geometry/calculus premises are listed in the accompanying logic map.
-/
structure WrittenInputs
    (Vertex : Type u) (dimension : Nat)
    (ActualFarthest WholeD18 : Vertex → Prop)
    (ActualHead : Vertex → HeadKind → Prop)
    (ActualForkPair : Vertex → ForkRole → ForkRole → Prop) : Prop where
  G18 : dimension ≤ 18
  V18 : dimension = 18 → ∃ p, ActualFarthest p ∧ WholeD18 p
  HeadCover : dimension = 18 → ∀ p, ActualFarthest p → WholeD18 p →
    ∃ head, ActualHead p head
  FiniteHead : dimension = 18 → ∀ p label, ActualFarthest p → WholeD18 p →
    ActualHead p (.finite label) → False
  ParallelHead : dimension = 18 → ∀ p, ActualFarthest p → WholeD18 p →
    ActualHead p .parallel → False
  DottedForkCover : dimension = 18 → ∀ p, ActualFarthest p → WholeD18 p →
    ActualHead p .dotted → ∃ left right, ActualForkPair p left right
  II : dimension = 18 → ∀ p left right, ActualFarthest p → WholeD18 p →
    ActualForkPair p left right → pairCase left right = .isolatedIsolated → False
  IN : dimension = 18 → ∀ p left right, ActualFarthest p → WholeD18 p →
    ActualForkPair p left right → pairCase left right = .isolatedNear → False
  IF : dimension = 18 → ∀ p left right, ActualFarthest p → WholeD18 p →
    ActualForkPair p left right → pairCase left right = .isolatedFar → False
  NN : dimension = 18 → ∀ p left right, ActualFarthest p → WholeD18 p →
    ActualForkPair p left right → pairCase left right = .nearNear → False
  FF : dimension = 18 → ∀ p left right, ActualFarthest p → WholeD18 p →
    ActualForkPair p left right → pairCase left right = .farFar → False
  DNearFar : dimension = 18 → ∀ p left right, ActualFarthest p → WholeD18 p →
    ActualHead p .dotted → ActualForkPair p left right →
    pairCase left right = .nearFar → False
  DNext : dimension = 18 → ∀ p left right, ActualFarthest p → WholeD18 p →
    ActualHead p .dotted → ActualForkPair p left right →
    pairCase left right = .containsFarNext → False

section Composition

variable {Vertex : Type u} {dimension : Nat}
variable {ActualFarthest WholeD18 : Vertex → Prop}
variable {ActualHead : Vertex → HeadKind → Prop}
variable {ActualForkPair : Vertex → ForkRole → ForkRole → Prop}
variable (written : WrittenInputs Vertex dimension ActualFarthest WholeD18
  ActualHead ActualForkPair)

include written

/-- Obtain the ACTUAL pair once and exhaust all seven classes at the SAME p. -/
theorem dotted_head_exclusion (at18 : dimension = 18) (p : Vertex)
    (farthest : ActualFarthest p) (whole : WholeD18 p)
    (head : ActualHead p .dotted) : False := by
  obtain ⟨left, right, pair⟩ :=
    written.DottedForkCover at18 p farthest whole head
  cases classified : pairCase left right with
  | isolatedIsolated => exact written.II at18 p left right farthest whole pair classified
  | isolatedNear => exact written.IN at18 p left right farthest whole pair classified
  | isolatedFar => exact written.IF at18 p left right farthest whole pair classified
  | nearNear => exact written.NN at18 p left right farthest whole pair classified
  | nearFar => exact written.DNearFar at18 p left right farthest whole head pair classified
  | farFar => exact written.FF at18 p left right farthest whole pair classified
  | containsFarNext => exact written.DNext at18 p left right farthest whole head pair classified

/-- Finite, parallel and ENTIRE dotted heads exhaust the actual written cover. -/
theorem whole_d18_exclusion (at18 : dimension = 18) (p : Vertex)
    (farthest : ActualFarthest p) (whole : WholeD18 p) : False := by
  obtain ⟨head, actualHead⟩ := written.HeadCover at18 p farthest whole
  cases head with
  | finite label => exact written.FiniteHead at18 p label farthest whole actualHead
  | parallel => exact written.ParallelHead at18 p farthest whole actualHead
  | dotted => exact dotted_head_exclusion written at18 p farthest whole actualHead

/-- Prior GLOBAL18 plus exact18 nonexistence; no face or rank heredity premise. -/
theorem compact_dimension_at_most_seventeen : dimension ≤ 17 := by
  by_cases bound : dimension ≤ 17
  · exact bound
  · have at18 : dimension = 18 := by
      have prior : dimension ≤ 18 := written.G18
      omega
    obtain ⟨p, farthest, whole⟩ := written.V18 at18
    exact False.elim (whole_d18_exclusion written at18 p farthest whole)

end Composition
end CoxeterProofs.CompactDimensionSeventeenComposition

#print axioms CoxeterProofs.CompactDimensionSeventeenComposition.dotted_head_exclusion
#print axioms CoxeterProofs.CompactDimensionSeventeenComposition.whole_d18_exclusion
#print axioms CoxeterProofs.CompactDimensionSeventeenComposition.compact_dimension_at_most_seventeen
