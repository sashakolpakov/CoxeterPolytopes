// Complete bounded rank gate. Every tested principal matrix avoids all
// free old-parallel entries. A nonzero determinant modulo a split prime
// certifies nonzero over Q(sqrt2,sqrt3,sqrt5), hence rank at least eleven.
#include <iostream>
#include <vector>
#include <array>
#include <cstdint>
#include <algorithm>
using namespace std;
int det11(const int G[16][16],const array<int,11>& b,int p,const vector<int>& iv){
 int A[11][11];for(int i=0;i<11;i++)for(int j=0;j<11;j++)A[i][j]=G[b[i]][b[j]];
 for(int k=0;k<11;k++){
  int r=k;while(r<11&&!A[r][k])r++;if(r==11)return 0;
  if(r!=k)for(int j=k;j<11;j++)swap(A[k][j],A[r][j]);
  for(int i=k+1;i<11;i++)if(A[i][k]){
   int f=(A[i][k]*iv[A[k][k]])%p;
   for(int j=k+1;j<11;j++)A[i][j]=(A[i][j]-f*A[k][j]+p*p)%p;
  }
 }return 1;
}
int main(){
 int n,ns,nb;cin>>n>>ns>>nb;int labels[16][16];for(int i=0;i<n;i++)for(int j=0;j<n;j++)cin>>labels[i][j];
 vector<vector<pair<int,int>>> supports(ns);for(auto& s:supports){int z;cin>>z;for(int k=0;k<z;k++){int i,j;cin>>i>>j;s.emplace_back(i,j);}}
 vector<array<int,11>> B(nb);for(auto& b:B)for(int& x:b)cin>>x;
 const int primes[2]={241,1009};int G[2][16][16],vals[2][7];vector<int>inv[2];
 for(int h=0;h<2;h++){
  int p=primes[h];inv[h].resize(p);for(int a=1;a<p;a++)for(int b=1;b<p;b++)if(a*b%p==1){inv[h][a]=b;break;}
  auto sq=[p](int d){for(int x=1;x<p;x++)if(x*x%p==d)return x;return -1;};
  int a=sq(2),b=sq(3),c=sq(5);if(a<0||b<0||c<0)return 2;
  vals[h][0]=p-2;vals[h][1]=2;vals[h][2]=0;vals[h][3]=p-1;vals[h][4]=p-a;vals[h][5]=(p-((1+c)*inv[h][2])%p)%p;vals[h][6]=p-b;
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)G[h][i][j]=vals[h][labels[i][j]];
 }
 uint64_t total=0,pass[2]={0,0},failed=0;
 for(int si=0;si<ns;si++){
  auto S=supports[si];int k=S.size();vector<array<int,4>> choices(k);
  for(int i=0;i<k;i++){int z=0;for(int m=2;m<=6;m++)if(m!=labels[S[i].first][S[i].second])choices[i][z++]=m;if(z!=4)return 3;}
  uint64_t N=uint64_t(1)<<(2*k);
  for(uint64_t code=0;code<N;code++){
   total++;uint64_t z=code;
   for(int i=k-1;i>=0;i--){int m=choices[i][z%4];z/=4;auto [a,b]=S[i];for(int h=0;h<2;h++)G[h][a][b]=G[h][b][a]=vals[h][m];}
   bool ok=false;
   for(int h=0;h<2&&!ok;h++)for(auto& b:B)if(det11(G[h],b,primes[h],inv[h])){pass[h]++;ok=true;break;}
   if(!ok){failed++;cout<<"U "<<si<<" "<<code<<"\n";}
  }
  // Restore the support before moving to a different mutation pattern.
  for(auto [a,b]:S)for(int h=0;h<2;h++)G[h][a][b]=G[h][b][a]=vals[h][labels[a][b]];
 }
 cout<<"COUNTS "<<total<<" "<<pass[0]<<" "<<pass[1]<<" "<<failed<<"\n";
}
