#include <list>
#include <vector>
#include <algorithm>
#include <Rcpp.h>
#include <math.h>
#include<iostream>
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

//function calculating the minimal required overlap between two reads of a genome to be considered safe for contig construction
//
//[[Rcpp::export]]
int calcMinOverlap(std::string seq,int meanWidth){
    int bases[] = {0,0,0,0};
    for(int i = 0;i < (int) seq.size();i++ ){
        if(seq.at(i) == 'A' ||seq.at(i) == 'a'){
            bases[0]++;
        }
        if(seq.at(i) == 'T' ||seq.at(i) == 't'){
            bases[1]++;
        }
        if(seq.at(i) == 'C' ||seq.at(i) == 'c'){
            bases[2]++;
        }
        if(seq.at(i) == 'G' ||seq.at(i) == 'g'){
            bases[3]++;
        }
    }

    int best = 0;
    for(int i = 0;i < 4;i++){
        if(bases[i] > best){
            best = bases[i];
        }
    }
    double tmp1 = (double)best/ (double) seq.size();
    double tmp2 = tmp1;
    int minOverlap = 1;
    while(tmp2 > (0.000001/meanWidth)){
        tmp2 *= tmp1;
        minOverlap++;
    }
    Rcpp::Rcout << "minOverlap ist: " << minOverlap <<"\n";
    return minOverlap;
}

// function to cut the sequence of the contigs out of the reference sequence
//
//[[Rcpp::export]]
std::vector<std::string> subSeqs(std::string seq,std::vector<int> starts,std::vector<int> ends){
    std::vector<std::string> res;
    std::vector<int>::iterator st = starts.begin();
    std::vector<int>::iterator en = ends.begin();
    for(int i = 0;i < starts.size();i++){
        if(*en <= seq.size()){
            res.push_back(seq.substr(*st,*en-*st+1));
        }
        else{
            res.push_back(seq.substr(*st,seq.size()-1) + seq.substr(0,*en-seq.size()));
        }
        st++;
        en++;
    }
    return res;
}

// function to calculate the mean coverage for all contigs for every sample
//
//[[Rcpp::export]]
Rcpp::List calcCovVec(std::list<std::vector<int> > readsPerSample,std::vector<int> lengths){
    Rcpp::List res;

    std::list<std::vector<int> >::iterator rps = readsPerSample.begin();
    std::vector<int>::iterator lngs;
    std::vector<int>::iterator rpsIt;

    std::vector<double> tmp;

    for(lngs = lengths.begin();lngs != lengths.end();lngs++){
        for(rpsIt = (*rps).begin();rpsIt != (*rps).end();rpsIt++){
            tmp.push_back((*rpsIt)/(double)(*lngs));
        }
        res.push_back(tmp);
        tmp.clear();
        rps++;
    }

    return res;
}
