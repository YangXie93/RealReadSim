#include <iostream>
#include "SeqRep.h"


//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::export]]

std::vector<int> evalCoverage(std::vector<int>& pos,std::vector<int>& width,int length,std::string seq){

  int minOverlap = calcMinOverlap(seq);

  SeqRep rep = SeqRep(length,minOverlap,&pos,&width);    //erstellen des SeqRep Objektes

  std::vector<bool> which;
  rep.evalOverlap(which);
  std::vector<int> res = rep.assembleTestContigs();
  return res;              //aufruf der auswertungs Funktion
}
