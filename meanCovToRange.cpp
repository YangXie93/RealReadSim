#include <vector>
#include <iostream>


//[[Rcpp::export]]

std::vector<int> meanCovToRange(std::vector<int> ranges,std::vector<int> coverage){
    std::vector<int> res;
    res.reserve(ranges.size()/2);
    int temp = 0;
    int width;

    for(int i = 0;i < (int) ranges.size();i += 2){
        width = ranges[i+1] -ranges[i];
        for(int j = ranges[i];j <= ranges[i +1];j++){
            temp += coverage[j];
        }
        res.push_back(temp/width);
        temp = 0;
    }
    return(res);
}
