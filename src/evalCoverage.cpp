#include <Rcpp.h>
#include <iostream>
#include "SeqRep.h"

using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::export]]

List evalCoverage(std::vector<int>& pos,std::vector<int>& width,int length,int minOverlap,int minContigLength){
  int meanWidth = 0;
  for(int i = 0;i < (int) width.size();i++){
    meanWidth += width[i];
  }
  meanWidth /= (int) width.size();

  SeqRep rep = SeqRep(length,minOverlap,&pos,&width);    //erstellen des SeqRep Objektes

  std::vector<bool> which;
  rep.evalOverlap(which);
  return rep.assembleTestContigs(minContigLength);
}
