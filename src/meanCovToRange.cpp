
#include <vector>
#include <iostream>
#include <iterator>
#include <list>


//[[Rcpp::export]]

std::list< std::vector<int> > meanCovToRange(std::vector<int> ranges,std::vector<int> coverage){
    std::list< std::vector<int> > res;
    std::vector<int> temp;
    int n = 0;

    for(int i = 0;i < (int) ranges.size();i += 2){
        for(int j = (ranges[i]-1);j <= (ranges[i +1]-1);j++){
            temp.push_back(coverage[n]);
            n++;
        }
        n = 0;
        res.push_back(temp);
        temp.clear();
    }
    return(res);
}
