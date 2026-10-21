import Std

/-!
Euler obstruction for an entirely ideal right-angled 5-polytope whose
4-facets are all 24-cells. Link-to-incidence geometry is an explicit
remaining formalization obligation; the numerical consequence is proved.
-/
namespace CoxeterProofs.IdealRA5

theorem all24cell_incidence_obstruction
    (v edges triangles octahedra cells : Nat)
    (edgeIncidence : edges = 8 * v)
    (triangleIncidence : 3 * triangles = 32 * v)
    (octahedronIncidence : octahedra = 4 * v)
    (cellIncidence : 3 * cells = v)
    (euler : v + triangles + cells = edges + octahedra + 2) : False := by
  omega

end CoxeterProofs.IdealRA5
