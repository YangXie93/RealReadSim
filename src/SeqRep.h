#ifndef SEQREP_H
#define SEQREP_H

#include <list>
#include <vector>
#include <iostream>
#include <algorithm>
#include <Rcpp.h>
#include <math.h>

// Class that represents a genome with a int[] holding the coverage value to every position of the genome sequence
// also the position and the sample identity to every read is stored in it for later evaluation

//[[Rcpp::plugins(cpp11)]]
class SeqRep{
public:

    //constructor
    //
    SeqRep(int length,int minOverlap,int nrOfSamples,std::vector<int> *pos,std::vector<int> *width,std::vector<int> *sample){

        this->nrOfSamples = nrOfSamples;
        this->length = length;
        this->minOverlap = minOverlap;
        nr = (int) pos->size();
        readStarts = new int[nr];
        readEnds = new int[nr];
        sampleNr = new int[nr];
        cov = new int[length];

        std::fill(cov,cov+length-1,0);

        std::vector<int>::iterator posIt = pos->begin();
        std::vector<int>::iterator widthIt = width->begin();
        std::vector<int>::iterator sampleIt = sample->begin();
        for(int i = 0;i < nr;i++){

          int* start = cov+(*posIt);
          int* end = start +(*widthIt)-1;
          if(end >= cov+length){
            end = end-length;
          }

          addToCov(start,end,1);
          *(readStarts+i) = (*posIt);
          *(readEnds+i) = (*posIt)+(*widthIt)-1;
          *(sampleNr+i) = *sampleIt;
          posIt++;
          widthIt++;
          sampleIt++;
        }

        evalOverlap();
    }

    // destructor
    //
    ~SeqRep(){
        delete[] readStarts;
        delete[] readEnds;
        delete[] cov;
        delete[] sampleNr;
    }


    // function to check wether a read has the minimum required overlap with others to qualify for
    // contig generation
    //
    void evalOverlap(){
      int count = 0;
      int start;
      int end;
      int tmp;
      int best;
      int sortedOut = 0;

      for(int i = 0;i < nr;i++){
          start = *(readStarts+i);
          end = *(readEnds+i);
          tmp = start;
          while(start != end+1 && tmp != end+1){
            if(*(cov+tmp) > 1){
              count++;
            }
            else{
              if(count > best){
                best = count;
              }
              count = 0;
            }
            start++;
            tmp++;
            if(tmp >= length){
              tmp = 0;
            }
          }
          if(best >= minOverlap){
            addToCov(cov+*(readStarts+i),cov+end,-1);
            toSkip.push_back(i);
            sortedOut++;
          }
          count = 0;
          best = 0;
      }
    }

    // function to get all contiguos sequences on the coverage vector
    //
    Rcpp::List assembleTestContigs(int minContigLength){
        std::vector<int> starts;
        std::vector<int> ends;
        Rcpp::List readsPerSample;
        Rcpp::List covs;

        int* covIt = cov;
        int i = 0;
        int* it;

        while(*covIt > 0 && i < length){
          covIt++;
          i++;
        }

        if(i >= length-1){
          starts.push_back(1);
          ends.push_back(length);
          readsPerSample.push_back(getReadsPerSampleOnRange(1,length));
          covs.push_back(std::vector<int> (cov,cov+length-1));
        }
        else{
          bool switching = true;

          if(covIt == cov){
            it = cov;
            covIt = cov+length;
          }
          else{
            it = covIt+1;
          }

          int end = it-cov;
          int start = 0;

          while(it != covIt){
            if(it >= cov+length){
              it = cov;
              end = 0;
            }

            if(*(it) > 0 && switching){
              start = end;
              switching = false;
            }
            if(*(it) == 0 && !switching){

              switching = true;
              if(end - start >= minContigLength){
                starts.push_back(start+1);
                ends.push_back(end);
                readsPerSample.push_back(getReadsPerSampleOnRange(start+1,end));
                covs.push_back(std::vector<int> (cov+start,cov+end));
              }
            }

            it++;
            end++;
          }

        }
        Rcpp::List res = Rcpp::List::create(starts, ends, covs,readsPerSample);
        return res;
    }

    //getter for Length
    int getLength(){
        return length;
    }


    //getter for cov
    int* getCov(){
      return cov;
    }

    // function to add a value to a given range on the coverage vector
    //
    void addToCov(int* start,int* end,int val){
      int* tmp = start;
      int x = 1;
      while(start != end+1 && tmp != end+1){
        (*tmp) += val;
        tmp++;
        start++;
        x++;
        if(tmp >= cov+length){
          tmp = cov;
        }
      }
    }

  // function to count the reads per sample on a given range of the genome
  //
  std::vector<int> getReadsPerSampleOnRange(int start,int end){
      int res[nrOfSamples];
      std::fill(res,(res+nrOfSamples),0);

      int* ts = toSkip.data();
      if(toSkip.size() == 0){
        int x = -1;
        ts = &x;
      }

      for(int i = 0;i < nr && *(readEnds+i) <= end;i++){

        if(*(readStarts+i) >= start && *(ts) != i){
            *(res+*(sampleNr+i)-1) += *(readEnds+i)-*(readStarts+i)+1;
            if(*ts != -1){
              ts++;
            }
        }

      }
      return std::vector<int> (res,(res+nrOfSamples));
  }

private:
    // the number of samples in the experiment
    int nrOfSamples;

    // minimal overlap requirement
    int minOverlap;

    // length of the genome sequence
    int length;

    // number of reads
    int nr;

    // int[] of coverage values
    int* cov;

    // start position of the reads
    int* readStarts;

    // end position of the reads
    int* readEnds;

    // sample identity of the reads
    int* sampleNr;

    // indicies of reads that have been sorted out by evalOverlap()
    std::vector<int> toSkip;
};





#endif
