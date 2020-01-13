#include <Rcpp.h>
#include <iostream>
#include "SeqRep.h"

using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::export]]

List evalCoverage(std::vector<int>& pos,std::vector<int>& width,int length,int minOverlap,int minContigLength){

    SeqRep rep = SeqRep(length,minOverlap,&pos,&width);    //erstellen des SeqRep Objektes
    return rep.assembleTestContigs(minContigLength);
}
