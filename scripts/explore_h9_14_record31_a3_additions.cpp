// Exact pivoted A5-addition/blocker census for all positive-A4 orbit
// representatives in Record 31's compact a=3 schedule.
//
// Input on stdin is: representative_count, followed by rows
//   production_index four_mask_0 four_mask_1 four_mask_2
// The Python certifier derives this projection exactly from the independently
// frozen positive-A4 orbit archive.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace {

constexpr int kN = 14;
constexpr int kScheduledFives = 27;
constexpr int kScheduledFacets = 116;
using Mask = std::uint16_t;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr std::array<Mask, 4> kEdges = {
    (1u << 0) | (1u << 11), (1u << 1) | (1u << 12),
    (1u << 0) | (1u << 13), (1u << 1) | (1u << 13),
};
constexpr std::array<int, kN> kBaseG3 = {2,2,7,7,7,7,7,7,7,7,7,5,5,4};
constexpr std::array<int, kN> kBaseG4 = {2,2,8,8,8,8,8,8,8,8,8,6,6,5};

struct Audit {
  int index = -1;
  std::array<Mask, 3> fours{};
  std::vector<Mask> shadow;
  std::vector<Mask> candidates;
  std::vector<Mask> small;
  std::vector<std::uint64_t> coverage;
  std::vector<std::vector<int>> coverers;
  std::vector<Mask> five_covers;
  std::vector<std::uint64_t> base_five_cover_words;
  std::vector<std::vector<std::uint64_t>> candidate_five_hit_words;
  std::array<int, kN> lower{};
  std::array<int, kN> upper{};
  int additions = 0;
  std::array<std::set<std::uint64_t>, 5> seen;
  std::uint64_t nodes = 0;
  std::uint64_t dead_uncoverable = 0;
  std::uint64_t no_small = 0;
  std::map<int, std::uint64_t> order_five_histogram;
  std::map<int, std::uint64_t> target_blocker_size_histogram;
  std::vector<std::vector<Mask>> pure_targets;
};

long long Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n-k);
  long long answer=1;
  for (int i=1;i<=k;++i) answer=answer*(n-k+i)/i;
  return answer;
}

int MacaulayRaise(int value, int degree) {
  int remaining=value, next_top=-1, answer=0;
  for (int lower=degree;lower>=1;--lower) {
    int top=lower-1;
    while ((next_top<0 || top+1<next_top) &&
           Binomial(top+1,lower)<=remaining) ++top;
    remaining -= Binomial(top,lower);
    answer += Binomial(top+1,lower+1);
    next_top=top;
  }
  if (remaining) throw std::runtime_error("Macaulay expansion failed");
  return answer;
}

bool Independent(Mask value) {
  for (const Mask edge:kEdges) if ((value&edge)==edge) return false;
  return true;
}

bool ContainsFour(Mask five, const std::array<Mask,3>& fours) {
  for (const Mask four:fours) if ((five&four)==four) return true;
  return false;
}

bool HitsAll(Mask value, const std::vector<Mask>& clutter) {
  for (const Mask face:clutter) if (!(value&face)) return false;
  return true;
}

std::vector<Mask> IndependentFives() {
  std::vector<Mask> answer;
  for (Mask value=0;value<=kFull;++value)
    if (std::popcount(value)==5 && Independent(value)) answer.push_back(value);
  if (answer.size()!=1296) throw std::runtime_error("A5 ground changed");
  return answer;
}

std::vector<Mask> GraphFiveCovers() {
  std::vector<Mask> answer;
  for (Mask value=0;value<=kFull;++value) {
    if (std::popcount(value)!=5) continue;
    bool cover=true;
    for (const Mask edge:kEdges) cover = cover && (value&edge);
    if (cover) answer.push_back(value);
  }
  return answer;
}

std::uint64_t ChosenKey(std::vector<int> chosen) {
  std::sort(chosen.begin(),chosen.end());
  std::uint64_t key=chosen.size();
  for (const int index:chosen) key=(key<<10)|index;
  return key;
}

bool Compatible(const Audit& audit, Mask face,
                const std::array<std::uint8_t,kN>& incidence) {
  for (int v=0;v<kN;++v)
    if (static_cast<int>(incidence[v])+static_cast<int>((face>>v)&1u)>
        audit.upper[v]) return false;
  return true;
}

std::vector<Mask> FullMinimalTransversals(const std::vector<Mask>& clutter) {
  std::vector<Mask> result;
  for (Mask value=1;value<=kFull;++value) {
    if (!HitsAll(value,clutter)) continue;
    bool minimal=true;
    for (int v=0;v<kN;++v) if ((value>>v)&1u)
      minimal = minimal && !HitsAll(value^(Mask{1}<<v),clutter);
    if (minimal) result.push_back(value);
  }
  return result;
}

void Finish(Audit* audit, const std::vector<int>& chosen,
            const std::vector<Mask>& base) {
  ++audit->no_small;
  std::vector<Mask> clutter=base;
  for (const int index:chosen) clutter.push_back(audit->candidates[index]);
  auto surviving=audit->base_five_cover_words;
  for (const int index:chosen) {
    for (int word=0;word<static_cast<int>(surviving.size());++word)
      surviving[word] &= audit->candidate_five_hit_words[index][word];
  }
  int order_five=0;
  for (const std::uint64_t word:surviving) order_five += std::popcount(word);
  ++audit->order_five_histogram[order_five];
  if (order_five!=kScheduledFacets) return;
  const auto blocker=FullMinimalTransversals(clutter);
  ++audit->target_blocker_size_histogram[blocker.size()];
  bool pure=blocker.size()==kScheduledFacets;
  for (const Mask face:blocker) pure = pure && std::popcount(face)==5;
  if (!pure) return;
  std::vector<Mask> family;
  family.insert(family.end(),audit->fours.begin(),audit->fours.end());
  for (const int index:chosen) family.push_back(audit->candidates[index]);
  audit->pure_targets.push_back(std::move(family));
}

void Search(Audit* audit, std::vector<int> chosen, std::uint64_t covered,
            std::array<std::uint8_t,kN> incidence,
            const std::vector<Mask>& base) {
  ++audit->nodes;
  if (!audit->seen[chosen.size()].insert(ChosenKey(chosen)).second) return;
  const std::uint64_t all=(std::uint64_t{1}<<audit->small.size())-1;
  if (static_cast<int>(chosen.size())==audit->additions) {
    if (covered!=all) return;
    for (int v=0;v<kN;++v)
      if (incidence[v]<audit->lower[v] || incidence[v]>audit->upper[v]) return;
    Finish(audit,chosen,base);
    return;
  }
  std::vector<int> best;
  std::uint64_t feasible_union=covered;
  bool has_uncovered=false;
  for (int requirement=0;requirement<static_cast<int>(audit->small.size());++requirement) {
    if ((covered>>requirement)&1u) continue;
    has_uncovered=true;
    std::vector<int> feasible;
    for (const int candidate:audit->coverers[requirement]) {
      if (std::find(chosen.begin(),chosen.end(),candidate)!=chosen.end()) continue;
      if (!Compatible(*audit,audit->candidates[candidate],incidence)) continue;
      feasible.push_back(candidate);
      feasible_union |= audit->coverage[candidate];
    }
    if (feasible.empty()) { ++audit->dead_uncoverable; return; }
    if (best.empty() || feasible.size()<best.size()) best=std::move(feasible);
  }
  if ((feasible_union&all)!=all) { ++audit->dead_uncoverable; return; }
  if (!has_uncovered) {
    for (int candidate=0;candidate<static_cast<int>(audit->candidates.size());++candidate)
      if (std::find(chosen.begin(),chosen.end(),candidate)==chosen.end() &&
          Compatible(*audit,audit->candidates[candidate],incidence))
        best.push_back(candidate);
  }
  for (const int candidate:best) {
    auto next_incidence=incidence;
    for (int v=0;v<kN;++v)
      next_incidence[v] += (audit->candidates[candidate]>>v)&1u;
    auto next=chosen;
    next.push_back(candidate);
    Search(audit,std::move(next),covered|audit->coverage[candidate],
           next_incidence,base);
  }
}

Audit Run(int index, const std::array<Mask,3>& fours,
          const std::vector<Mask>& fives,
          const std::vector<Mask>& five_covers) {
  Audit audit;
  audit.index=index;
  audit.fours=fours;
  audit.five_covers=five_covers;
  for (const Mask five:fives) if (ContainsFour(five,fours)) audit.shadow.push_back(five);
  audit.additions=kScheduledFives-audit.shadow.size();
  if (audit.additions<1 || audit.additions>4) throw std::runtime_error("addition count changed");
  for (int v=0;v<kN;++v) {
    int a4=0,forced=0;
    for (const Mask four:fours) a4 += (four>>v)&1u;
    for (const Mask five:audit.shadow) forced += (five>>v)&1u;
    const int previous=kBaseG3[v]-a4;
    const int before=kBaseG4[v]+6*a4-forced;
    if (previous<1 || before<1) throw std::runtime_error("strict gate changed");
    audit.lower[v]=std::max(0,before-MacaulayRaise(previous,3));
    audit.upper[v]=before-1;
  }
  int lower_sum=0,upper_sum=0;
  for (int v=0;v<kN;++v) { lower_sum+=audit.lower[v]; upper_sum+=audit.upper[v]; }
  if (!(lower_sum<=5*audit.additions && 5*audit.additions<=upper_sum))
    throw std::runtime_error("aggregate gate changed");
  for (const Mask five:fives) {
    if (ContainsFour(five,fours)) continue;
    bool allowed=true;
    for (int v=0;v<kN;++v) if (((five>>v)&1u) && audit.upper[v]<1) allowed=false;
    if (allowed) audit.candidates.push_back(five);
  }
  std::vector<Mask> base(kEdges.begin(),kEdges.end());
  base.insert(base.end(),fours.begin(),fours.end());
  base.insert(base.end(),audit.shadow.begin(),audit.shadow.end());
  for (Mask value=1;value<=kFull;++value)
    if (std::popcount(value)<=4 && HitsAll(value,base)) audit.small.push_back(value);
  if (audit.small.size()>63) throw std::runtime_error("small gate escaped word");
  audit.coverers.assign(audit.small.size(),{});
  for (int candidate=0;candidate<static_cast<int>(audit.candidates.size());++candidate) {
    std::uint64_t coverage=0;
    for (int r=0;r<static_cast<int>(audit.small.size());++r)
      if (!(audit.candidates[candidate]&audit.small[r])) {
        coverage |= std::uint64_t{1}<<r;
        audit.coverers[r].push_back(candidate);
      }
    audit.coverage.push_back(coverage);
  }
  const int five_word_count=(five_covers.size()+63)/64;
  audit.base_five_cover_words.assign(five_word_count,0);
  for (int cover=0;cover<static_cast<int>(five_covers.size());++cover) {
    bool hits=true;
    for (const Mask four:fours) hits = hits && (five_covers[cover]&four);
    if (hits) audit.base_five_cover_words[cover/64] |=
        std::uint64_t{1}<<(cover%64);
  }
  audit.candidate_five_hit_words.assign(
      audit.candidates.size(),std::vector<std::uint64_t>(five_word_count,0));
  for (int candidate=0;candidate<static_cast<int>(audit.candidates.size());++candidate)
    for (int cover=0;cover<static_cast<int>(five_covers.size());++cover)
      if (audit.candidates[candidate]&five_covers[cover])
        audit.candidate_five_hit_words[candidate][cover/64] |=
            std::uint64_t{1}<<(cover%64);
  Search(&audit,{},0,{},base);
  return audit;
}

void PrintMask(Mask value) {
  std::cout<<'[';
  bool first=true;
  for (int v=0;v<kN;++v) if ((value>>v)&1u) {
    if (!first) std::cout<<',';
    std::cout<<v;
    first=false;
  }
  std::cout<<']';
}

}  // namespace

int main() {
  try {
    int count=0;
    if (!(std::cin>>count) || count!=143) throw std::runtime_error("expected 143 representatives");
    std::vector<std::pair<int,std::array<Mask,3>>> input;
    for (int row=0;row<count;++row) {
      int index;
      unsigned a,b,c;
      if (!(std::cin>>index>>a>>b>>c) || index!=row || a>kFull || b>kFull || c>kFull)
        throw std::runtime_error("malformed representative input");
      input.push_back({index,{static_cast<Mask>(a),static_cast<Mask>(b),static_cast<Mask>(c)}});
    }
    const auto fives=IndependentFives();
    const auto five_covers=GraphFiveCovers();
    std::vector<Audit> audits;
    for (const auto& [index,fours]:input) {
      audits.push_back(Run(index,fours,fives,five_covers));
      std::cerr<<"orbit "<<index<<" add="<<audits.back().additions
               <<" no-small="<<audits.back().no_small
               <<" targets="<<audits.back().pure_targets.size()<<'\n';
    }
    std::cout<<"{\n  \"kind\":\"h9_14_record31_a3_pivoted_addition_scout\",\n"
             <<"  \"compact_only\":true,\n  \"representative_audits\":[";
    std::uint64_t total_targets=0;
    for (int i=0;i<count;++i) {
      if (i) std::cout<<',';
      const auto& a=audits[i];
      total_targets += a.pure_targets.size();
      std::cout<<"{\"orbit_representative_index\":"<<a.index
               <<",\"forced_five_shadow_count\":"<<a.shadow.size()
               <<",\"additional_minimal_A5_count\":"<<a.additions
               <<",\"individual_candidate_count\":"<<a.candidates.size()
               <<",\"raw_unordered_addition_count\":"<<Binomial(a.candidates.size(),a.additions)
               <<",\"addition_lower\":[";
      for (int v=0;v<kN;++v) { if (v) std::cout<<','; std::cout<<a.lower[v]; }
      std::cout<<"],\"addition_upper\":[";
      for (int v=0;v<kN;++v) { if (v) std::cout<<','; std::cout<<a.upper[v]; }
      std::cout<<"],\"base_small_transversal_count\":"<<a.small.size()
               <<",\"search_node_count\":"<<a.nodes
               <<",\"memo_state_counts\":[";
      for (int depth=0;depth<=a.additions;++depth) {
        if (depth) std::cout<<',';
        std::cout<<a.seen[depth].size();
      }
      std::cout<<"],\"dead_uncoverable_count\":"<<a.dead_uncoverable
               <<",\"no_transversal_of_order_at_most_four_count\":"<<a.no_small
               <<",\"order_five_transversal_count_histogram\":{";
      bool first=true;
      for (const auto& [facets,n]:a.order_five_histogram) {
        if (!first) std::cout<<',';
        std::cout<<'\"'<<facets<<"\":"<<n;
        first=false;
      }
      std::cout<<"},\"scheduled_order_five_count_blocker_size_histogram\":{";
      first=true;
      for (const auto& [size,n]:a.target_blocker_size_histogram) {
        if (!first) std::cout<<',';
        std::cout<<'\"'<<size<<"\":"<<n;
        first=false;
      }
      std::cout<<"},\"pure_116_facet_target_count\":"<<a.pure_targets.size()
               <<",\"pure_target_minimal_nonfaces\":[";
      first=true;
      for (const auto& family:a.pure_targets) {
        if (!first) std::cout<<',';
        std::cout<<'[';
        for (int row=0;row<static_cast<int>(family.size());++row) {
          if (row) std::cout<<',';
          PrintMask(family[row]);
        }
        std::cout<<']';
        first=false;
      }
      std::cout<<"]}";
    }
    std::cout<<"],\n  \"pure_116_facet_target_count_across_representatives\":"<<total_targets
             <<",\n  \"scope\":\"Exact pivoted literal addition and blocker gate for all 143 compact a=3 positive-A4 orbits. No finite-volume inference.\"\n}\n";
  } catch (const std::exception& error) {
    std::cerr<<error.what()<<'\n';
    return 1;
  }
}
