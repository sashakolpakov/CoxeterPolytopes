// Independent literal-combination audit of all twenty Record-31 compact a=4
// representatives.  Unlike the production split (Python for one/two rows and
// pivoted C++ for three), this verifier loops over every cap-feasible literal
// addition combination directly and applies the small-transversal and exact
// blocker-count gates.

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
using Mask = std::uint16_t;
constexpr Mask kFull = (Mask{1} << kN) - 1;
constexpr std::array<Mask, 4> kEdges = {
    (1u << 0) | (1u << 11), (1u << 1) | (1u << 12),
    (1u << 0) | (1u << 13), (1u << 1) | (1u << 13),
};
constexpr std::array<std::array<Mask, 4>, 20> kFours = {{
    {{1472,7680,8640,11776}}, {{480,7680,8640,11776}},
    {{3840,5888,7680,9984}}, {{3840,5760,7680,9984}},
    {{3840,5888,7296,9984}}, {{3840,5888,6336,9984}},
    {{3840,5888,7424,11776}}, {{3840,5504,7424,11776}},
    {{3840,5888,9408,9984}}, {{3840,5888,8416,9984}},
    {{2272,5888,8416,9984}}, {{3456,5504,9600,11776}},
    {{2496,4544,8640,11776}}, {{3840,5888,9984,15360}},
    {{2944,4992,9088,15360}}, {{3840,5888,11520,13824}},
    {{3456,7680,9600,13824}}, {{2496,7680,8640,13824}},
    {{3840,7680,11520,13824}}, {{6528,7680,10624,13824}},
}};
constexpr std::array<int, kN> kBaseG3 = {2,2,7,7,7,7,7,7,7,7,7,5,5,4};
constexpr std::array<int, kN> kBaseG4 = {2,2,8,8,8,8,8,8,8,8,8,6,6,5};

struct Candidate {
  Mask mask;
  std::uint64_t coverage;
};

struct Audit {
  int index = -1;
  int shadow_count = 0;
  int additions = 0;
  int small_count = 0;
  int candidate_count = 0;
  std::uint64_t raw_count = 0;
  std::uint64_t cap_feasible = 0;
  std::uint64_t no_small = 0;
  std::map<int, std::uint64_t> blocker_histogram;
};

long long Binomial(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n-k);
  long long value = 1;
  for (int i=1;i<=k;++i) value = value*(n-k+i)/i;
  return value;
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

bool ContainsFour(Mask five, const std::array<Mask,4>& fours) {
  for (const Mask four:fours) if ((five&four)==four) return true;
  return false;
}

bool HitsAll(Mask value, const std::vector<Mask>& clutter) {
  for (const Mask face:clutter) if (!(value&face)) return false;
  return true;
}

int BlockerCount(const std::vector<Mask>& clutter) {
  int count=0;
  for (Mask value=1;value<=kFull;++value) {
    if (!HitsAll(value,clutter)) continue;
    bool minimal=true;
    for (int v=0;v<kN;++v) if ((value>>v)&1u)
      minimal = minimal && !HitsAll(value^(Mask{1}<<v),clutter);
    count += minimal;
  }
  return count;
}

std::vector<Mask> IndependentFives() {
  std::vector<Mask> answer;
  for (Mask value=0;value<=kFull;++value)
    if (std::popcount(value)==5 && Independent(value)) answer.push_back(value);
  if (answer.size()!=1296) throw std::runtime_error("A5 ground changed");
  return answer;
}

Audit Run(int index, const std::vector<Mask>& fives) {
  const auto& fours=kFours[index];
  std::vector<Mask> shadow;
  for (const Mask five:fives) if (ContainsFour(five,fours)) shadow.push_back(five);
  const int additions=33-shadow.size();
  if (additions<1 || additions>3) throw std::runtime_error("addition count changed");
  std::array<int,kN> upper{};
  for (int v=0;v<kN;++v) {
    int a4=0, forced=0;
    for (const Mask four:fours) a4 += (four>>v)&1u;
    for (const Mask five:shadow) forced += (five>>v)&1u;
    const int previous=kBaseG3[v]-a4;
    const int before=kBaseG4[v]+6*a4-forced;
    const int lower=std::max(0,before-MacaulayRaise(previous,3));
    if (previous<1 || before<1 || lower!=0) throw std::runtime_error("local gate changed");
    upper[v]=before-1;
  }
  Mask cap_one=0, cap_two=0, forbidden=0;
  for (int v=0;v<kN;++v) {
    if (upper[v]<1) forbidden |= Mask{1}<<v;
    if (upper[v]<2) cap_one |= Mask{1}<<v;
    if (upper[v]<3) cap_two |= Mask{1}<<v;
  }
  std::vector<Mask> base(kEdges.begin(),kEdges.end());
  base.insert(base.end(),fours.begin(),fours.end());
  base.insert(base.end(),shadow.begin(),shadow.end());
  std::vector<Mask> small;
  for (Mask value=1;value<=kFull;++value)
    if (std::popcount(value)<=4 && HitsAll(value,base)) small.push_back(value);
  if (small.size()>63) throw std::runtime_error("small gate escaped word");
  const std::uint64_t all=(std::uint64_t{1}<<small.size())-1;
  std::vector<Candidate> candidates;
  for (const Mask five:fives) {
    if (ContainsFour(five,fours) || (five&forbidden)) continue;
    std::uint64_t coverage=0;
    for (int r=0;r<static_cast<int>(small.size());++r)
      if (!(five&small[r])) coverage |= std::uint64_t{1}<<r;
    candidates.push_back({five,coverage});
  }
  Audit audit{index,static_cast<int>(shadow.size()),additions,
              static_cast<int>(small.size()),static_cast<int>(candidates.size()),
              static_cast<std::uint64_t>(Binomial(candidates.size(),additions)),
              0,0,{}};
  auto accept=[&](std::array<int,3> chosen) {
    ++audit.cap_feasible;
    std::uint64_t coverage=0;
    auto clutter=base;
    for (int row=0;row<additions;++row) {
      coverage |= candidates[chosen[row]].coverage;
      clutter.push_back(candidates[chosen[row]].mask);
    }
    if (coverage!=all) return;
    ++audit.no_small;
    ++audit.blocker_histogram[BlockerCount(clutter)];
  };
  if (additions==1) {
    for (int i=0;i<static_cast<int>(candidates.size());++i) accept({i,0,0});
  } else if (additions==2) {
    for (int i=0;i<static_cast<int>(candidates.size());++i)
      for (int j=i+1;j<static_cast<int>(candidates.size());++j) {
        if (candidates[i].mask&candidates[j].mask&cap_one) continue;
        accept({i,j,0});
      }
  } else {
    for (int i=0;i<static_cast<int>(candidates.size());++i)
      for (int j=i+1;j<static_cast<int>(candidates.size());++j) {
        if (candidates[i].mask&candidates[j].mask&cap_one) continue;
        for (int k=j+1;k<static_cast<int>(candidates.size());++k) {
          if ((candidates[i].mask&candidates[k].mask&cap_one) ||
              (candidates[j].mask&candidates[k].mask&cap_one) ||
              (candidates[i].mask&candidates[j].mask&candidates[k].mask&cap_two))
            continue;
          accept({i,j,k});
        }
      }
  }
  return audit;
}

}  // namespace

int main() {
  try {
    const auto fives=IndependentFives();
    std::vector<Audit> audits;
    for (int index=0;index<20;++index) {
      audits.push_back(Run(index,fives));
      std::cerr<<"orbit "<<index<<" no-small="<<audits.back().no_small<<'\n';
    }
    constexpr std::array<int,20> expected_candidates = {
        6,21,764,56,126,246,56,56,252,449,126,126,246,56,126,56,126,252,462,462};
    constexpr std::array<int,20> expected_small = {
        22,21,43,42,36,30,36,36,36,30,27,36,30,36,30,36,42,35,43,43};
    constexpr std::array<std::uint64_t,20> expected_cap = {
        6,21,3951530,945,2310,8505,540,945,8541,27693,126,2310,8505,
        540,2310,540,2310,8541,1794605,2819841};
    for (int index=0;index<20;++index) {
      const std::uint64_t expected_no_small = index==18 ? 270 : index==19 ? 1121 : 0;
      if (audits[index].candidate_count!=expected_candidates[index] ||
          audits[index].small_count!=expected_small[index] ||
          audits[index].cap_feasible!=expected_cap[index] ||
          audits[index].no_small!=expected_no_small ||
          audits[index].blocker_histogram.count(114))
        throw std::runtime_error("independent Record-31 a=4 census changed");
    }
    if (audits[18].blocker_histogram !=
            std::map<int,std::uint64_t>{{105,150},{109,120}} ||
        audits[19].blocker_histogram != std::map<int,std::uint64_t>{
            {94,160},{95,120},{96,340},{99,80},{100,260},{102,150},
            {104,10},{115,1}})
      throw std::runtime_error("independent Record-31 blocker census changed");
    std::cout<<"{\n  \"kind\":\"h9_14_record31_a4_independent_literal_audit\",\n"
             <<"  \"compact_only\":true,\n  \"representative_audits\":[";
    for (int i=0;i<20;++i) {
      if (i) std::cout<<',';
      const auto& a=audits[i];
      std::cout<<"{\"orbit_representative_index\":"<<a.index
               <<",\"forced_five_shadow_count\":"<<a.shadow_count
               <<",\"additional_minimal_A5_count\":"<<a.additions
               <<",\"individual_candidate_count\":"<<a.candidate_count
               <<",\"raw_unordered_addition_count\":"<<a.raw_count
               <<",\"base_small_transversal_count\":"<<a.small_count
               <<",\"vertex_cap_feasible_count\":"<<a.cap_feasible
               <<",\"no_transversal_of_order_at_most_four_count\":"<<a.no_small
               <<",\"blocker_count_histogram\":{";
      bool first=true;
      for (const auto& [facets,count]:a.blocker_histogram) {
        if (!first) std::cout<<',';
        std::cout<<'\"'<<facets<<"\":"<<count;
        first=false;
      }
      const auto target=a.blocker_histogram.find(114);
      std::cout<<"},\"scheduled_114_facet_target_count\":"
               <<(target==a.blocker_histogram.end()?0:target->second)<<'}';
    }
    std::cout<<"],\n  \"scheduled_target_count_across_representatives\":0,\n"
             <<"  \"scope\":\"Independent direct literal combination audit of all twenty compact a=4 representatives; no finite-volume inference.\"\n}\n";
  } catch (const std::exception& error) {
    std::cerr<<error.what()<<'\n';
    return 1;
  }
}
