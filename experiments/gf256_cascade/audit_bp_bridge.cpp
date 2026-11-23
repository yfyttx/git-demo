#include "NB_BP_Decoder.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
// Input: m n edges cases iterations; then (row,col,polynomial value) edges;
// each case: fm dual legacy_check, q prior values, m syndrome values.
int main(int argc,char**argv) {
 try {
  if(argc!=3) throw std::runtime_error("usage: bridge input output");
  std::ifstream in(argv[1]); std::ofstream out(argv[2]);
  int m,n,edges,cases,iters; in>>m>>n>>edges>>cases>>iters;
  if(!in || m<1 || n<1 || edges<1 || cases<1 || iters<1) throw std::runtime_error("bad header");
  GF2p gf(8,{1,0,0,0,1,1,1,0,1});
  std::vector<Eigen::Triplet<int>> entries;
  for(int e=0;e<edges;e++){int r,c,a;in>>r>>c>>a;entries.emplace_back(r,c,a);}
  Eigen::SparseMatrix<int,Eigen::RowMajor> h(m,n); h.setFromTriplets(entries.begin(),entries.end());
  NB_BP_Decoder bp(gf,iters);
  for(int k=0;k<cases;k++) {
   double fm;int dual,legacy;in>>fm>>dual>>legacy;
   std::vector<double> prior(256);for(auto&v:prior)in>>v;
   std::vector<int> syn(m);for(auto&v:syn)in>>v;
   if(!in)throw std::runtime_error("truncated case");
   int used=0;auto estimate=bp.decode_component_prior(h,syn,prior,&used);
   if(legacy && bp.decode_component(h,syn,fm,dual)!=estimate)throw std::runtime_error("legacy wrapper mismatch");
   bool ok=bp.check_z_syndrome(h,estimate,syn);
   out<<used<<' '<<ok;for(int v:estimate)out<<' '<<v;out<<'\n';
  }
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
