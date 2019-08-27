#include <Rcpp.h>
#include <iostream>
#include "SeqRep.h"

using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::export]]

std::vector<int> evalCoverage(std::vector<int>& pos,std::vector<int>& width,int length,std::string seq){
  int meanWidth = 0;
  for(int i = 0;i < (int) width.size();i++){
    meanWidth += width[i];
  }
  meanWidth /= (int) width.size();
  int minOverlap = calcMinOverlap(seq,meanWidth);

  SeqRep rep = SeqRep(length,minOverlap,&pos,&width);    //erstellen des SeqRep Objektes

  std::vector<bool> which;
  rep.evalOverlap(which);
  std::vector<int> res = rep.assembleTestContigs();
  return res;              //aufruf der auswertungs Funktion
}
