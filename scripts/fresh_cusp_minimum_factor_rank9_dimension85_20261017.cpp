// Exact integer enumeration of all affine-factor rank partitions through162.
// This supplies only the local curvature/payment comparison; see the report.
#include <array>
#include <cassert>
#include <cstdio>
#include <string>
using U=unsigned __int128;
using I=__int128;
using A=std::array<U,163>;
std::array<unsigned long long,163> counts{};
std::array<U,163> least_scaled_gap{};
std::string decimal(U n){if(!n)return "0";std::string s;while(n){s.push_back('0'+n%10);n/=10;}return std::string(s.rbegin(),s.rend());}
void visit(int d,int lo,int twice_cost,U vertices,const A &h){
 if(d){
  U sum=0,g=0;for(int j=0;j<=d;++j){sum+=h[j];g+=h[j]*(j<d-j?j:d-j);assert(h[j]==h[d-j]);}
  assert(sum==vertices);++counts[d];
  // The mathematical bound V<=2^81 makes every operation below <2^94.
  assert(vertices <= (U(1)<<81));
  if(d>=85){
   I diff=I(5*g)-I(U(d+twice_cost)*vertices);
   // diff/(2V) is curvature per link vertex minus distance-seven cost.
   I gap=5*diff-I(2*vertices);assert(gap>0);
   if(!least_scaled_gap[d] || U(gap)<least_scaled_gap[d])least_scaled_gap[d]=U(gap);
  }
 }
 for(int r=lo;d+r<=162;++r){
  A next{};U window=0;
  for(int j=0;j<=d+r;++j){
   if(j<=d)window+=h[j];
   if(j-r-1>=0&&j-r-1<=d){assert(window>=h[j-r-1]);window-=h[j-r-1];}
   next[j]=window;
  }
  visit(d+r,r,twice_cost+(r<14?r:14),vertices*U(r+1),next);
 }
}
int main(){
 A h{};h[0]=1;visit(0,9,0,1,h);
 puts("[");for(int d=9;d<=162;++d){
  printf("{\"ambient_dimension\":%d,\"partitions\":%llu,\"minimum_integer_scaled_gap\":\"%s\"}%s\n",d+1,counts[d],decimal(least_scaled_gap[d]).c_str(),d==162?"":",");
 }puts("]");
}
