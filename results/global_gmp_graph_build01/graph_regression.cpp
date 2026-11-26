
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
struct IVL { std::vector<std::vector<int> > rows; };
const int IVL_CHUNKED=0;
void IVL_init1(IVL *a,int,int m){a->rows.resize(m);}
void IVL_setList(IVL *a,int i,int n,const int *p){
 if(n)a->rows[i].assign(p,p+n);else a->rows[i].clear();}
int IVL_tsize(IVL *a){int n=0;for(auto &r:a->rows)n+=r.size();return n;}
void Graph_init2(void*,int,int,int,int,int,int,IVL*,void*,void*){}
#define rMessage(x) throw std::runtime_error("reference graph error")
struct SparseLinearSpace {
 int SDP_sp_nBlock,SOCP_sp_nBlock,LP_sp_nBlock;
 int *SDP_sp_index,*SOCP_sp_index,*LP_sp_index;
};
struct InputData {
 int SDP_nBlock,SOCP_nBlock,LP_nBlock;
 int *SDP_nConstraint,*SOCP_nConstraint,*LP_nConstraint;
 int **SDP_constraint,**SOCP_constraint,**LP_constraint;
 SparseLinearSpace *A;
};
namespace ref {struct Chordal {
 IVL store;IVL *adjIVL=&store;void *graph=nullptr;
 void makeGraph(InputData&,int);void margeArray(int,int*,int,int*);
};}
namespace fast {struct Chordal {
 IVL store;IVL *adjIVL=&store;void *graph=nullptr;
 void makeGraph(InputData&,int);
};}
void ref::Chordal::margeArray(int na1, int *array1, int na2, int *array2) {

    int ptr = na1 + na2 - 1;
    int ptr1 = na1 - 1;
    int ptr2 = na2 - 1;
    int idx1, idx2;

    while ((ptr1 >= 0) || (ptr2 >= 0)) {

        if (ptr1 >= 0) {
            idx1 = array1[ptr1];
        } else {
            idx1 = -1;
        }
        if (ptr2 >= 0) {
            idx2 = array2[ptr2];
        } else {
            idx2 = -1;
        }
        if (idx1 > idx2) {
            array2[ptr] = idx1;
            ptr1--;
        } else {
            array2[ptr] = idx2;
            ptr2--;
        }
        ptr--;
    }

    // error check
    if (ptr != -1) {
        rMessage("ref::Chordal::margeArray:: program bug");
    }
}

void ref::Chordal::makeGraph(InputData &inputData, int m) {

    int i, j, k, l;
    int SDP_nBlock = inputData.SDP_nBlock;
    int SOCP_nBlock = inputData.SOCP_nBlock;
    int LP_nBlock = inputData.LP_nBlock;

    int *counter;
    counter = new int[m];
    for (int i = 0; i < m; i++) {
        counter[i] = 0;
    }

    // count maximum mumber of index
    for (l = 0; l < SDP_nBlock; l++) {
        int SDP_nConstraint = inputData.SDP_nConstraint[l];
        for (k = 0; k < SDP_nConstraint; k++) {
            i = inputData.SDP_constraint[l][k];
            counter[i] += SDP_nConstraint;
        }
    }
    for (l = 0; l < SOCP_nBlock; l++) {
        int SOCP_nConstraint = inputData.SOCP_nConstraint[l];
        for (k = 0; k < SOCP_nConstraint; k++) {
            i = inputData.SOCP_constraint[l][k];
            counter[i] += SOCP_nConstraint;
        }
    }
    for (l = 0; l < LP_nBlock; l++) {
        int LP_nConstraint = inputData.LP_nConstraint[l];
        for (k = 0; k < LP_nConstraint; k++) {
            i = inputData.LP_constraint[l][k];
            counter[i] += LP_nConstraint;
        }
    }

    // allocate temporaly workspace
    int **tmp;
    tmp = new int *[m];
    for (i = 0; i < m; i++) {
        tmp[i] = new int[counter[i]];
    }

    // merge index
    for (int i = 0; i < m; i++) {
        counter[i] = 0;
    }
    // marge index of for SDP
    for (l = 0; l < SDP_nBlock; l++) {
        for (k = 0; k < inputData.SDP_nConstraint[l]; k++) {
            i = inputData.SDP_constraint[l][k];
            margeArray(inputData.SDP_nConstraint[l], inputData.SDP_constraint[l], counter[i], tmp[i]);
            counter[i] += inputData.SDP_nConstraint[l];
        }
    }
    // marge index of for SOCP
    for (l = 0; l < SOCP_nBlock; l++) {
        for (k = 0; k < inputData.SOCP_nConstraint[l]; k++) {
            i = inputData.SOCP_constraint[l][k];
            margeArray(inputData.SOCP_nConstraint[l], inputData.SOCP_constraint[l], counter[i], tmp[i]);
            counter[i] += inputData.SOCP_nConstraint[l];
        }
    }
    // marge index of for LP
    for (l = 0; l < LP_nBlock; l++) {
        for (k = 0; k < inputData.LP_nConstraint[l]; k++) {
            i = inputData.LP_constraint[l][k];
            margeArray(inputData.LP_nConstraint[l], inputData.LP_constraint[l], counter[i], tmp[i]);
            counter[i] += inputData.LP_nConstraint[l];
        }
    }

    // construct adjacency list of SPOOLES
    IVL_init1(adjIVL, IVL_CHUNKED, m);
    int isize, previous;
    int *ivec;
    ivec = new int[m];
    for (i = 0; i < m; i++) {
        isize = 0;
        previous = -1;
        for (j = 0; j < counter[i]; j++) {
            if (tmp[i][j] != previous) {
                ivec[isize] = tmp[i][j];
                previous = ivec[isize];
                isize++;
            }
        }
        IVL_setList(adjIVL, i, isize, ivec);
    }

    // constract graph of SPOOLES
    Graph_init2(graph, 0, m, 0, IVL_tsize(adjIVL), m, IVL_tsize(adjIVL), adjIVL, NULL, NULL);

    delete[] counter;
    for (int i = 0; i < m; i++) {
        delete[] tmp[i];
    }
    delete[] tmp;
    delete[] ivec;
}

void fast::Chordal::makeGraph(InputData &inputData, int m) {
    // Same union of block-support cliques as upstream; mark duplicates once.
    IVL_init1(adjIVL, IVL_CHUNKED, m);
    std::vector<int> marks(m, -1), row;
    for (int i = 0; i < m; ++i) {
        row.clear();
        const auto append = [&](const int *members, int count) {
            for (int k = 0; k < count; ++k) {
                const int j = members[k];
                if (marks[j] != i) { marks[j] = i; row.push_back(j); }
            }
        };
        const SparseLinearSpace &A = inputData.A[i];
        for (int k = 0; k < A.SDP_sp_nBlock; ++k) {
            const int b = A.SDP_sp_index[k];
            append(inputData.SDP_constraint[b], inputData.SDP_nConstraint[b]);
        }
        for (int k = 0; k < A.SOCP_sp_nBlock; ++k) {
            const int b = A.SOCP_sp_index[k];
            append(inputData.SOCP_constraint[b], inputData.SOCP_nConstraint[b]);
        }
        for (int k = 0; k < A.LP_sp_nBlock; ++k) {
            const int b = A.LP_sp_index[k];
            append(inputData.LP_constraint[b], inputData.LP_nConstraint[b]);
        }
        std::sort(row.begin(), row.end());
        IVL_setList(adjIVL, i, static_cast<int>(row.size()), row.data());
    }
    Graph_init2(graph, 0, m, 0, IVL_tsize(adjIVL), m,
                IVL_tsize(adjIVL), adjIVL, NULL, NULL);
}


int main(){
 std::mt19937 generator(20261126);long long checked=0;
 for(int test=0;test<300;++test){
  int m=1+generator()%96;
  std::vector<std::vector<int> > supports[3],incidence[3];
  std::vector<int> counts[3];std::vector<int*> pointers[3];
  for(int g=0;g<3;++g){
   incidence[g].resize(m);int blocks=generator()%90;
   supports[g].resize(blocks);
   for(int b=0;b<blocks;++b){
    for(int j=0;j<m;++j)if(generator()%7==0){
     supports[g][b].push_back(j);incidence[g][j].push_back(b);}
    counts[g].push_back(supports[g][b].size());
    pointers[g].push_back(supports[g][b].data());
   }
  }
  std::vector<SparseLinearSpace> rows(m);
  for(int i=0;i<m;++i)rows[i]={
   (int)incidence[0][i].size(),(int)incidence[1][i].size(),(int)incidence[2][i].size(),
   incidence[0][i].data(),incidence[1][i].data(),incidence[2][i].data()};
  InputData data={(int)supports[0].size(),(int)supports[1].size(),(int)supports[2].size(),
   counts[0].data(),counts[1].data(),counts[2].data(),
   pointers[0].data(),pointers[1].data(),pointers[2].data(),rows.data()};
  ref::Chordal a;fast::Chordal b;a.makeGraph(data,m);b.makeGraph(data,m);
  if(a.store.rows!=b.store.rows)return 2;
  // Also check the defining co-occurrence relation independently.
  for(int i=0;i<m;++i){std::vector<int> expected;
   for(int j=0;j<m;++j){bool meet=false;
    for(int g=0;g<3;++g)for(auto &s:supports[g])
     if(std::binary_search(s.begin(),s.end(),i)&&std::binary_search(s.begin(),s.end(),j))meet=true;
    if(meet)expected.push_back(j);}
   if(expected!=b.store.rows[i])return 3;checked+=expected.size();}
 }
 std::cout<<"300 exact old/new/independent graph cases passed; adjacency entries "<<checked<<"\n";
}
