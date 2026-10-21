// Independent exact scheduled-target replay for all 143 Record-31 compact
// a=3 representatives.
//
// The production search memoizes unordered chosen sets after choosing a
// smallest uncovered transversal.  This replay instead partitions the
// solution set: for the selected uncovered transversal, a branch fixes the
// least still-feasible destroyer that occurs in the final addition family
// and permanently forbids all earlier destroyers.  Hence every literal
// completion has one path and no chosen-set memo table is used.  Intersections
// of bitsets of order-five graph covers provide an independent exact
// scheduled-facet gate before the full blocker is rebuilt literally.

// Input on stdin is the same small projection used by production:
//   143
//   index four_mask_0 four_mask_1 four_mask_2


#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr int kScheduledFives = 27;
constexpr int kScheduledFacets = 116;
using Mask = std::uint16_t;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr std::array<Mask, 4> kMissingEdges = {
    (1u << 0) | (1u << 11), (1u << 1) | (1u << 12),
    (1u << 0) | (1u << 13), (1u << 1) | (1u << 13),
};
constexpr std::array<int, kN> kLinkG3 = {2,2,7,7,7,7,7,7,7,7,7,5,5,4};
constexpr std::array<int, kN> kLinkG4 = {2,2,8,8,8,8,8,8,8,8,8,6,6,5};

struct Candidate {
  Mask face = 0;
  std::uint64_t destroys_small = 0;
  std::vector<std::uint64_t> hits_five_covers;
};

struct SearchState {
  std::vector<int> chosen;
  std::vector<std::uint64_t> forbidden;
  std::uint64_t destroyed_small = 0;
  std::array<std::uint8_t, kN> incidence{};
  std::vector<std::uint64_t> surviving_five_covers;
};

struct Audit {
  int index = -1;
  std::array<Mask, 3> fours{};
  std::vector<Mask> shadow;
  std::vector<Mask> small;
  std::vector<Mask> graph_five_covers;
  std::vector<Candidate> candidates;
  std::vector<std::vector<int>> destroyers;
  std::array<int, kN> upper{};
  int additions = 0;
  std::uint64_t search_nodes = 0;
  std::uint64_t facet_lower_bound_prunes = 0;
  std::uint64_t uncovered_union_prunes = 0;
  std::uint64_t scheduled_states = 0;
  std::map<int, std::uint64_t> scheduled_blocker_size_histogram;
  std::vector<std::vector<Mask>> pure_targets;
};

long long Choose(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n-k);
  long long answer = 1;
  for (int i=1; i<=k; ++i) answer = answer*(n-k+i)/i;
  return answer;
}

int MacaulaySuccessor(int value, int degree) {
  int remaining=value, next_top=-1, answer=0;
  for (int lower=degree; lower>=1; --lower) {
    int top=lower-1;
    while ((next_top<0 || top+1<next_top) &&
           Choose(top+1,lower)<=remaining) ++top;
    remaining -= Choose(top,lower);
    answer += Choose(top+1,lower+1);
    next_top=top;
  }
  if (remaining) throw std::runtime_error("Macaulay expansion failed");
  return answer;
}

bool AvoidsMissingEdges(Mask face) {
  for (const Mask edge:kMissingEdges)
    if ((face&edge)==edge) return false;
  return true;
}

bool ContainsOneOf(Mask face, const std::array<Mask,3>& fours) {
  for (const Mask four:fours)
    if ((face&four)==four) return true;
  return false;
}

bool HitsEvery(Mask transversal, const std::vector<Mask>& clutter) {
  for (const Mask face:clutter)
    if (!(transversal&face)) return false;
  return true;
}

std::vector<Mask> IndependentFives() {
  std::vector<Mask> result;
  for (Mask face=0; face<=kFull; ++face)
    if (std::popcount(face)==5 && AvoidsMissingEdges(face)) result.push_back(face);
  if (result.size()!=1296) throw std::runtime_error("independent-five ground changed");
  return result;
}

std::vector<Mask> FiveVertexCovers() {
  std::vector<Mask> result;
  for (Mask cover=0; cover<=kFull; ++cover) {
    if (std::popcount(cover)!=5) continue;
    bool hits=true;
    for (const Mask edge:kMissingEdges) hits = hits && (cover&edge);
    if (hits) result.push_back(cover);
  }
  return result;
}

int BitCount(const std::vector<std::uint64_t>& words) {
  int answer=0;
  for (const std::uint64_t word:words) answer += std::popcount(word);
  return answer;
}

bool IsSet(const std::vector<std::uint64_t>& words, int index) {
  return (words[index/64]>>(index%64))&1u;
}

void Set(std::vector<std::uint64_t>* words, int index) {
  (*words)[index/64] |= std::uint64_t{1}<<(index%64);
}

bool CapCompatible(const Audit& audit, const SearchState& state, Mask face) {
  for (int vertex=0; vertex<kN; ++vertex)
    if (static_cast<int>(state.incidence[vertex])+static_cast<int>((face>>vertex)&1u)>
        audit.upper[vertex])
      return false;
  return true;
}

std::vector<std::uint64_t> Meet(
    const std::vector<std::uint64_t>& left,
    const std::vector<std::uint64_t>& right) {
  std::vector<std::uint64_t> answer(left.size());
  for (int word=0; word<static_cast<int>(left.size()); ++word)
    answer[word]=left[word]&right[word];
  return answer;
}

std::vector<Mask> LiteralMinimalBlocker(const std::vector<Mask>& clutter) {
  std::vector<Mask> blocker;
  for (Mask candidate=1; candidate<=kFull; ++candidate) {
    if (!HitsEvery(candidate,clutter)) continue;
    bool minimal=true;
    for (int vertex=0; vertex<kN; ++vertex)
      if ((candidate>>vertex)&1u &&
          HitsEvery(candidate^(Mask{1}<<vertex),clutter)) minimal=false;
    if (minimal) blocker.push_back(candidate);
  }
  return blocker;
}

SearchState AddCandidate(
    const Audit& audit, const SearchState& state, int candidate_index,
    std::vector<std::uint64_t> next_forbidden) {
  const Candidate& candidate=audit.candidates[candidate_index];
  SearchState next=state;
  next.chosen.push_back(candidate_index);
  next.forbidden=std::move(next_forbidden);
  next.destroyed_small |= candidate.destroys_small;
  for (int vertex=0; vertex<kN; ++vertex)
    next.incidence[vertex] += (candidate.face>>vertex)&1u;
  next.surviving_five_covers = Meet(
      state.surviving_five_covers,candidate.hits_five_covers);
  return next;
}

void Finish(Audit* audit, const SearchState& state) {
  const std::uint64_t all_small =
      (std::uint64_t{1}<<audit->small.size())-1;
  if (state.destroyed_small!=all_small ||
      BitCount(state.surviving_five_covers)!=kScheduledFacets) return;
  ++audit->scheduled_states;
  std::vector<Mask> clutter(kMissingEdges.begin(),kMissingEdges.end());
  clutter.insert(clutter.end(),audit->fours.begin(),audit->fours.end());
  clutter.insert(clutter.end(),audit->shadow.begin(),audit->shadow.end());
  for (const int index:state.chosen)
    clutter.push_back(audit->candidates[index].face);
  const auto blocker=LiteralMinimalBlocker(clutter);
  ++audit->scheduled_blocker_size_histogram[blocker.size()];
  bool pure=blocker.size()==kScheduledFacets;
  for (const Mask transversal:blocker) pure = pure && std::popcount(transversal)==5;
  if (!pure) return;
  std::vector<Mask> family(audit->fours.begin(),audit->fours.end());
  for (const int index:state.chosen) family.push_back(audit->candidates[index].face);
  std::sort(family.begin(),family.end());
  audit->pure_targets.push_back(std::move(family));
}

void Search(Audit* audit, const SearchState& state);

void Fill(Audit* audit, const SearchState& state, int first_candidate) {
  ++audit->search_nodes;
  const int remaining=audit->additions-static_cast<int>(state.chosen.size());
  if (BitCount(state.surviving_five_covers)<kScheduledFacets) {
    ++audit->facet_lower_bound_prunes;
    return;
  }
  if (!remaining) {
    Finish(audit,state);
    return;
  }
  for (int candidate=first_candidate;
       candidate<static_cast<int>(audit->candidates.size()); ++candidate) {
    if (IsSet(state.forbidden,candidate) ||
        !CapCompatible(*audit,state,audit->candidates[candidate].face)) continue;
    const auto surviving=Meet(state.surviving_five_covers,
                              audit->candidates[candidate].hits_five_covers);
    if (BitCount(surviving)<kScheduledFacets) continue;
    auto forbidden=state.forbidden;
    Set(&forbidden,candidate);
    auto next=AddCandidate(*audit,state,candidate,std::move(forbidden));
    next.surviving_five_covers=surviving;
    Fill(audit,next,candidate+1);
  }
}

void Search(Audit* audit, const SearchState& state) {
  ++audit->search_nodes;
  const int remaining=audit->additions-static_cast<int>(state.chosen.size());
  const int current_facets=BitCount(state.surviving_five_covers);
  if (current_facets<kScheduledFacets) {
    ++audit->facet_lower_bound_prunes;
    return;
  }
  if (!remaining) {
    Finish(audit,state);
    return;
  }
  const std::uint64_t all_small =
      (std::uint64_t{1}<<audit->small.size())-1;
  if (state.destroyed_small==all_small) {
    Fill(audit,state,0);
    return;
  }

  std::vector<int> pivot_candidates;
  std::uint64_t feasible_union=state.destroyed_small;
  for (int requirement=0;
       requirement<static_cast<int>(audit->small.size()); ++requirement) {
    if ((state.destroyed_small>>requirement)&1u) continue;
    std::vector<int> feasible;
    for (const int candidate:audit->destroyers[requirement]) {
      if (IsSet(state.forbidden,candidate) ||
          !CapCompatible(*audit,state,audit->candidates[candidate].face)) continue;
      const auto surviving=Meet(state.surviving_five_covers,
                                audit->candidates[candidate].hits_five_covers);
      if (BitCount(surviving)<kScheduledFacets) continue;
      feasible.push_back(candidate);
      feasible_union |= audit->candidates[candidate].destroys_small;
    }
    if (feasible.empty()) {
      ++audit->uncovered_union_prunes;
      return;
    }
    if (pivot_candidates.empty() || feasible.size()<pivot_candidates.size())
      pivot_candidates=std::move(feasible);
  }
  if ((feasible_union&all_small)!=all_small) {
    ++audit->uncovered_union_prunes;
    return;
  }

  // These branches are disjoint: the chosen candidate is declared to be the
  // least feasible destroyer of this pivot that belongs to the final family.
  auto cumulative_forbidden=state.forbidden;
  for (const int candidate:pivot_candidates) {
    auto branch_forbidden=cumulative_forbidden;
    Set(&branch_forbidden,candidate);
    Search(audit,AddCandidate(
        *audit,state,candidate,std::move(branch_forbidden)));
    Set(&cumulative_forbidden,candidate);
  }
}

Audit Run(int index, const std::array<Mask,3>& fours,
          const std::vector<Mask>& independent_fives,
          const std::vector<Mask>& five_covers) {
  Audit audit;
  audit.index=index;
  audit.fours=fours;
  audit.graph_five_covers=five_covers;
  for (const Mask five:independent_fives)
    if (ContainsOneOf(five,fours)) audit.shadow.push_back(five);
  audit.additions=kScheduledFives-static_cast<int>(audit.shadow.size());
  if (audit.additions<1 || audit.additions>4)
    throw std::runtime_error("addition schedule changed");

  for (int vertex=0; vertex<kN; ++vertex) {
    int four_incidence=0, shadow_incidence=0;
    for (const Mask four:fours) four_incidence += (four>>vertex)&1u;
    for (const Mask five:audit.shadow) shadow_incidence += (five>>vertex)&1u;
    const int g3=kLinkG3[vertex]-four_incidence;
    const int g4_before=kLinkG4[vertex]+6*four_incidence-shadow_incidence;
    const int lower=std::max(0,g4_before-MacaulaySuccessor(g3,3));
    if (g3<1 || g4_before<1 || lower!=0)
      throw std::runtime_error("local link interval changed");
    audit.upper[vertex]=g4_before-1;
  }

  std::vector<Mask> base(kMissingEdges.begin(),kMissingEdges.end());
  base.insert(base.end(),fours.begin(),fours.end());
  base.insert(base.end(),audit.shadow.begin(),audit.shadow.end());
  for (Mask transversal=1; transversal<=kFull; ++transversal)
    if (std::popcount(transversal)<=4 && HitsEvery(transversal,base))
      audit.small.push_back(transversal);
  if (audit.small.size()>63) throw std::runtime_error("small gate escaped one word");

  std::vector<Mask> candidate_faces;
  for (const Mask five:independent_fives) {
    if (ContainsOneOf(five,fours)) continue;
    bool possible=true;
    for (int vertex=0; vertex<kN; ++vertex)
      if (((five>>vertex)&1u) && audit.upper[vertex]<1) possible=false;
    if (possible) candidate_faces.push_back(five);
  }
  const int five_word_count=(five_covers.size()+63)/64;
  audit.candidates.reserve(candidate_faces.size());
  for (const Mask face:candidate_faces) {
    Candidate candidate;
    candidate.face=face;
    for (int requirement=0;
         requirement<static_cast<int>(audit.small.size()); ++requirement)
      if (!(face&audit.small[requirement]))
        candidate.destroys_small |= std::uint64_t{1}<<requirement;
    candidate.hits_five_covers.assign(five_word_count,0);
    for (int cover=0; cover<static_cast<int>(five_covers.size()); ++cover)
      if (face&five_covers[cover])
        candidate.hits_five_covers[cover/64] |=
            std::uint64_t{1}<<(cover%64);
    audit.candidates.push_back(std::move(candidate));
  }
  audit.destroyers.assign(audit.small.size(),{});
  for (int candidate=0; candidate<static_cast<int>(audit.candidates.size()); ++candidate)
    for (int requirement=0;
         requirement<static_cast<int>(audit.small.size()); ++requirement)
      if ((audit.candidates[candidate].destroys_small>>requirement)&1u)
        audit.destroyers[requirement].push_back(candidate);

  SearchState initial;
  initial.forbidden.assign((audit.candidates.size()+63)/64,0);
  initial.surviving_five_covers.assign(five_word_count,0);
  for (int cover=0; cover<static_cast<int>(five_covers.size()); ++cover) {
    bool hits=true;
    for (const Mask four:fours) hits = hits && (five_covers[cover]&four);
    if (hits) initial.surviving_five_covers[cover/64] |=
        std::uint64_t{1}<<(cover%64);
  }
  Search(&audit,initial);
  std::sort(audit.pure_targets.begin(),audit.pure_targets.end());
  return audit;
}

void PrintMaskFamily(const std::vector<Mask>& family) {
  std::cout<<'[';
  for (int index=0; index<static_cast<int>(family.size()); ++index) {
    if (index) std::cout<<',';
    std::cout<<family[index];
  }
  std::cout<<']';
}

}  // namespace

int main() {
  try {
    int count=0;
    if (!(std::cin>>count) || count!=143)
      throw std::runtime_error("expected 143 representatives");
    std::vector<std::pair<int,std::array<Mask,3>>> input;
    for (int row=0; row<count; ++row) {
      int index;
      unsigned first,second,third;
      if (!(std::cin>>index>>first>>second>>third) || index!=row ||
          first>kFull || second>kFull || third>kFull)
        throw std::runtime_error("malformed representative input");
      input.push_back({index,{static_cast<Mask>(first),static_cast<Mask>(second),
                              static_cast<Mask>(third)}});
    }
    const auto independent_fives=IndependentFives();
    const auto five_covers=FiveVertexCovers();
    std::vector<Audit> audits;
    for (const auto& [index,fours]:input) {
      audits.push_back(Run(index,fours,independent_fives,five_covers));
      std::cerr<<"orbit "<<index<<" scheduled="<<audits.back().scheduled_states
               <<" pure="<<audits.back().pure_targets.size()<<'\n';
    }

    std::uint64_t total_scheduled=0,total_pure=0;
    std::cout<<"{\n  \"kind\":\"h9_14_record31_a3_independent_partition_audit\",\n"
             <<"  \"compact_only\":true,\n  \"representative_audits\":[";
    for (int row=0; row<count; ++row) {
      if (row) std::cout<<',';
      const Audit& audit=audits[row];
      total_scheduled += audit.scheduled_states;
      total_pure += audit.pure_targets.size();
      std::cout<<"{\"orbit_representative_index\":"<<audit.index
               <<",\"forced_five_shadow_count\":"<<audit.shadow.size()
               <<",\"additional_minimal_A5_count\":"<<audit.additions
               <<",\"individual_candidate_count\":"<<audit.candidates.size()
               <<",\"base_small_transversal_count\":"<<audit.small.size()
               <<",\"partition_search_node_count\":"<<audit.search_nodes
               <<",\"facet_lower_bound_prune_count\":"<<audit.facet_lower_bound_prunes
               <<",\"uncovered_union_prune_count\":"<<audit.uncovered_union_prunes
               <<",\"scheduled_order_five_state_count\":"<<audit.scheduled_states
               <<",\"scheduled_blocker_size_histogram\":{";
      bool first=true;
      for (const auto& [size,multiplicity]:audit.scheduled_blocker_size_histogram) {
        if (!first) std::cout<<',';
        std::cout<<'\"'<<size<<"\":"<<multiplicity;
        first=false;
      }
      std::cout<<"},\"pure_116_facet_target_count\":"<<audit.pure_targets.size()
               <<",\"pure_target_minimal_nonface_masks\":[";
      first=true;
      for (const auto& family:audit.pure_targets) {
        if (!first) std::cout<<',';
        PrintMaskFamily(family);
        first=false;
      }
      std::cout<<"]}";
    }
    std::cout<<"],\n  \"scheduled_order_five_state_count_across_representatives\":"
             <<total_scheduled
             <<",\n  \"pure_116_facet_target_count_across_representatives\":"<<total_pure
             <<",\n  \"scope\":\"Independent least-destroyer partition of every literal compact a=3 addition family, exact scheduled order-five cover intersection, and literal full-blocker replay. No finite-volume inference.\"\n}\n";
  } catch (const std::exception& error) {
    std::cerr<<error.what()<<'\n';
    return 1;
  }
}
