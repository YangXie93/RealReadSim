#include <iostream>
#include <vector>
#include "SeqRep.h"

//[[Rcpp::export]]

std::vector<int> getCoverage(std::vector<int>& pos,std::vector<int>& width,int length){
    SeqRep gen = SeqRep(length,0,&pos,&width);

    int* seq = gen.getCov();
    std::vector<int> res;
    res.reserve(length);
    for(int i = 0;i < length;i++){
        res.push_back(*(seq+i));
    }
    return res;
}
