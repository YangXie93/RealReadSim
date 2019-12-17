#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <list>

using namespace Rcpp;
//[[Rcpp::export]]

std::list< std::vector<int> > meanCovToRange(std::vector<int> ranges,std::vector<int> coverage){
    std::list< std::vector<int> > res;
    std::vector<int> temp;
    int n;
    for(int i = 0;i < (int) ranges.size();i += 2){
        n = ranges[i];
        for(int j = (ranges[i]);j < (ranges[i +1]);j++){
            temp.push_back(coverage[n]);
            n++;
        }
        res.push_back(temp);
        temp.clear();
    }
    return(res);
}
