#include<iostream>
#include<vector>
#include<fstream>
#include<string>


//[[Rcpp::plugins(cpp14)]]
//[[Rcpp::export]]

bool sequenceToFastaReads(std::vector<int>& starts,std::string& sequence,int meanWidth,std::string& newFasta,std::string& nameTag){
    bool x = true;
    if(std::ifstream(newFasta)){
      x = false;
    }
    std::ofstream outfile (x ? std::ofstream(newFasta):std::ofstream(newFasta,std::ios::app));
    if(outfile.is_open()){
        for(int i = 0;i < (int) starts.size();i++){
            outfile << ">"+nameTag +"_" +std::to_string(starts[i]+1) << std::endl;
            if(starts[i] +meanWidth < (int) sequence.length()){
                outfile << sequence.substr(starts[i],meanWidth) << std::endl;
            }
            else{
                outfile << sequence.substr(starts[i],sequence.length() -1) << sequence.substr(0,meanWidth - (sequence.length() -starts[i])) << std::endl;
            }
        }
        outfile.close();
        return true;
    }
    else{
        return false;
    }
}
