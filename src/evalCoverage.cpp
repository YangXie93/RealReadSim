#include <Rcpp.h>
#include <iostream>
#include "SeqRep.h"

using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::export]]

List evalCoverage(std::vector<int>& pos,std::vector<int>& width,std::vector<int>& sampleID,int length,int minOverlap,int minContigLength,int nrOfSamples){

    SeqRep rep = SeqRep(length,minOverlap,nrOfSamples,&pos,&width,&sampleID);    //erstellen des SeqRep Objektes
    return rep.assembleTestContigs(minContigLength);
}
