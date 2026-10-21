import CoxeterProofs.Generated.Ideal58Chunk00
import CoxeterProofs.Generated.Ideal58Chunk01
import CoxeterProofs.Generated.Ideal58Chunk02
import CoxeterProofs.Generated.Ideal58Chunk03
import CoxeterProofs.Generated.Ideal58Chunk04
import CoxeterProofs.Generated.Ideal58Chunk05
import CoxeterProofs.Generated.Ideal58Chunk06
import CoxeterProofs.Generated.Ideal58Chunk07
import CoxeterProofs.Generated.Ideal58Chunk08
import CoxeterProofs.Generated.Ideal58Chunk09
import CoxeterProofs.Generated.Ideal58Chunk10
import CoxeterProofs.Generated.Ideal58Chunk11
import CoxeterProofs.Generated.Ideal58Chunk12
import CoxeterProofs.Generated.Ideal58Chunk13
import CoxeterProofs.Generated.Ideal58Chunk14
import CoxeterProofs.Generated.Ideal58Chunk15
import CoxeterProofs.Generated.Ideal58Chunk16
import CoxeterProofs.Generated.Ideal58Chunk17
import CoxeterProofs.Generated.Ideal58Chunk18
import CoxeterProofs.Generated.Ideal58Chunk19
set_option maxRecDepth 100000
set_option maxHeartbeats 0
namespace CoxeterProofs.Generated.Ideal58
open CoxeterProofs.IdealProfile
def profiles : List Profile := chunk00 ++ chunk01 ++ chunk02 ++ chunk03 ++ chunk04 ++ chunk05 ++ chunk06 ++ chunk07 ++ chunk08 ++ chunk09 ++ chunk10 ++ chunk11 ++ chunk12 ++ chunk13 ++ chunk14 ++ chunk15 ++ chunk16 ++ chunk17 ++ chunk18 ++ chunk19 ++ chunk20 ++ chunk21 ++ chunk22 ++ chunk23 ++ chunk24 ++ chunk25 ++ chunk26 ++ chunk27 ++ chunk28 ++ chunk29 ++ chunk30 ++ chunk31 ++ chunk32 ++ chunk33 ++ chunk34 ++ chunk35 ++ chunk36 ++ chunk37 ++ chunk38 ++ chunk39 ++ chunk40 ++ chunk41 ++ chunk42 ++ chunk43 ++ chunk44 ++ chunk45 ++ chunk46 ++ chunk47 ++ chunk48 ++ chunk49 ++ chunk50 ++ chunk51 ++ chunk52 ++ chunk53 ++ chunk54 ++ chunk55 ++ chunk56 ++ chunk57 ++ chunk58 ++ chunk59 ++ chunk60 ++ chunk61 ++ chunk62 ++ chunk63 ++ chunk64 ++ chunk65 ++ chunk66 ++ chunk67 ++ chunk68 ++ chunk69 ++ chunk70 ++ chunk71 ++ chunk72 ++ chunk73 ++ chunk74 ++ chunk75 ++ chunk76 ++ chunk77 ++ chunk78 ++ chunk79
theorem all_checked : profiles.all (CoxeterProofs.Ideal58.check envelopes) = true := by
  simp only [profiles, List.all_append, chunk00_checked, chunk01_checked, chunk02_checked, chunk03_checked, chunk04_checked, chunk05_checked, chunk06_checked, chunk07_checked, chunk08_checked, chunk09_checked, chunk10_checked, chunk11_checked, chunk12_checked, chunk13_checked, chunk14_checked, chunk15_checked, chunk16_checked, chunk17_checked, chunk18_checked, chunk19_checked, chunk20_checked, chunk21_checked, chunk22_checked, chunk23_checked, chunk24_checked, chunk25_checked, chunk26_checked, chunk27_checked, chunk28_checked, chunk29_checked, chunk30_checked, chunk31_checked, chunk32_checked, chunk33_checked, chunk34_checked, chunk35_checked, chunk36_checked, chunk37_checked, chunk38_checked, chunk39_checked, chunk40_checked, chunk41_checked, chunk42_checked, chunk43_checked, chunk44_checked, chunk45_checked, chunk46_checked, chunk47_checked, chunk48_checked, chunk49_checked, chunk50_checked, chunk51_checked, chunk52_checked, chunk53_checked, chunk54_checked, chunk55_checked, chunk56_checked, chunk57_checked, chunk58_checked, chunk59_checked, chunk60_checked, chunk61_checked, chunk62_checked, chunk63_checked, chunk64_checked, chunk65_checked, chunk66_checked, chunk67_checked, chunk68_checked, chunk69_checked, chunk70_checked, chunk71_checked, chunk72_checked, chunk73_checked, chunk74_checked, chunk75_checked, chunk76_checked, chunk77_checked, chunk78_checked, chunk79_checked, Bool.and_self]
theorem every_profile_checked (p : Profile) (member : p ∈ profiles) :
    CoxeterProofs.Ideal58.check envelopes p = true := (List.all_eq_true.mp all_checked) p member
theorem every_profile_gap (p : Profile) (member : p ∈ profiles) :
    72*p.cost.den+57*p.vertices*p.cost.den+2*p.vertices*p.cost.num ≤ 5*p.indexTotal*p.cost.den :=
  CoxeterProofs.Ideal58.absolute_gap_of_check envelopes p (every_profile_checked p member)
theorem domain_complete : profiles.map Profile.ranks = CoxeterProofs.Ideal58.terminalDomain := by decide
theorem profile_count : profiles.length = 2552 := by decide
end CoxeterProofs.Generated.Ideal58
