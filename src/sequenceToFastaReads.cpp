#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<Rcpp.h>

using namespace std;
//[[Rcpp::plugins(cpp14)]]
//[[Rcpp::export]]

bool sequenceToFastaReads(vector<int>& starts,string& sequence,int meanWidth,string& newFasta,string& nameTag){
    bool x = true;
    if(ifstream(newFasta)){
      x = false;
    }
    ofstream outfile (x ? ofstream(newFasta):ofstream(newFasta,ios::app));
    if(outfile.is_open()){
        for(int i = 0;i < (int) starts.size();i++){
            outfile << ">"+nameTag +"_" +to_string(starts[i]+1) << endl;
            if(starts[i] +meanWidth < (int) sequence.length()){
                outfile << sequence.substr(starts[i],meanWidth) << endl;
            }
            else{
                outfile << sequence.substr(starts[i],sequence.length() -1) << sequence.substr(0,meanWidth - (sequence.length() -starts[i])) << endl;
            }
        }
        outfile.close();
        return true;
    }
    else{
        return false;
    }
}
