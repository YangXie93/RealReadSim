#ifndef SEQREP_H
#define SEQREP_H

#include <list>
#include <vector>
#include <iostream>
#include <algorithm>
#include <Rcpp.h>
#include <math.h>

/* das SeqRep Objekt ist eine abstrakte Darstellung einer Sequenz, auf die eine bestimmte Menge an Reads gemapped
   wird. In dem Objekt ist ein array "cov" der Laenge ("length") der Sequenz auf die gemapped wurde. "cov" ist ein integerarray.
   Fuer jeden Read, der auf eine Position gemapped wurde wird das entsprechende Element von "cov" um 1 erhoeht.
   Ein weiterer Bestandteil des Objektes ist ein vector<Read> "ids" in dem die Identitaeten der einzelnen Reads gespeichert werden.
   Die restlichen Klassenvariablen sind "minOverlap", dass den minimalen overlap der Reads, die akzeptiert werden angibt,
   "coverage", dass die minimale akzeptierte coverage einer Stelle auf der Sequenz angibt und "n", das die Anzahl der reads angibt.
   Es gibt eine Funktion zur Eingabe von Reads, eine Zur Auswertung der erhobenen Daten und einen constructor.
*/

//[[Rcpp::plugins(cpp11)]]
class SeqRep{
public:
    //constructor
    SeqRep(int length,int minOverlap,std::vector<int> *pos,std::vector<int> *width){

        this->length = length;
        this->minOverlap = minOverlap;
        nr = (int) pos->size();
        readStarts = new int*[nr];
        readEnds = new int*[nr];
        cov = new int[length];

        std::fill(cov,cov+length-1,0);

        std::vector<int>::iterator posIt = pos->begin();
        std::vector<int>::iterator widthIt = width->begin();
        for(int i = 0;i < nr;i++){

          int* start = cov+(*posIt);
          int* end = start +(*widthIt)-1;
          if(end >= cov+length){
            end = end-length;
          }

          addToCov(start,end,1);
          *(readStarts+i) = start;
          *(readEnds+i) = end;
          posIt++;
          widthIt++;
        }

        evalOverlap();
    }

    ~SeqRep(){
        delete[] readStarts;
        delete[] readEnds;
        delete[] cov;
    }


    //Funktion zum auswerten des overlaps
    /* Bei jedem Read bei dem Eine Stelle auf die er gemapped wurde >= "coverage" ist werden hier diese Stellen gezaehlt.
       sobald eine folge an Stellen >= "coverage" abbricht wird der bisherige "count" gespeichert.
       Sind dann alle Stellen des Reads abgearbeitet werden die "counts" ueberprueft, ob mindestens einer davon dem
       "minOverlap" entspricht. Ist dies nicht der Fall wird der vorher zugeordnete boolean wert wieder umgedreht.
    */
    void evalOverlap(){    //Das Ergebnis aus evalCov() wird hier eingegeben

      int count = 0;
      int* start;
      int* end;
      int* tmp;
      int best;
      int sortedOut = 0;

      for(int i = 0;i < nr;i++){
          start = *(readStarts+i);
          end = *(readEnds+i);
          tmp = start;
          while(start != end+1 && tmp != end+1){
            if(*tmp > 1){
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
            if(tmp >= cov+length){
              tmp = cov;
            }
          }
          if(best >= minOverlap){
            addToCov(*(readStarts+i),end,-1);
            sortedOut++;
          }
          count = 0;
          best = 0;
      }
      Rcpp::Rcout << sortedOut << " reads have been sorted out\n";
    }


    Rcpp::List assembleTestContigs(int minContigLength){

        std::vector<int> starts;
        std::vector<int> ends;
        std::vector<double> mCovs;
        std::vector<double> covSD;
        Rcpp::List covs;
        double mean;

        int* covIt = cov;
        int i = 0;
        int* it;

        while(*covIt > 0 && i < length){
          covIt++;
          i++;
        }

        if(covIt == cov+length-1){
          starts.push_back(1);
          ends.push_back(length);
          mean = meanCovOnRange(1,length);
          mCovs.push_back(mean);
          covSD.push_back(SDCovOnRange(1,length,mean));
          covs.push_back(std::vector<int> (cov,cov+length-1));
        }
        else{

          int n = 0;
          bool switching = true;

          if(covIt == cov){
            it = cov;
            covIt = cov+length-1;
          }
          else{
            it = covIt+1;
          }
          i = it-cov;

          while(it != covIt){
            if(*(it) > 0 && switching){
              n = i;
              switching = false;
            }
            if(*(it) == 0 && !switching){
              switching = true;
              if(i -n >= minContigLength){
                starts.push_back(n+1);
                ends.push_back(i+1);
                mean = meanCovOnRange(n,i);
                mCovs.push_back(mean);
                covSD.push_back(SDCovOnRange(n,i,mean));
                covs.push_back(std::vector<int> (cov+n,cov+i));
              }
            }

            it++;
            i++;
            if(it >= cov+length){
              it = cov;
              i = 0;
            }
          }

        }
        Rcpp::List res = Rcpp::List::create(starts, ends, covs, mCovs, covSD);
        return res;
    }


    double meanCovOnRange(int start,int end){
      return std::accumulate(cov+start,cov+end,0)/double (end-start);
    }


    double SDCovOnRange(int start,int end,double mean){
      std::vector<double> tmp;
      for(int i = start;i <= end;i++){
        tmp.push_back(std::pow((*(cov+i)-mean),2) );
      }
      return std::sqrt(std::accumulate(tmp.begin(),tmp.end(),0)/(end-start));
    }


    //getter fuer Length
    int getLength(){
        return length;
    }


    //getter fuer cov
    int* getCov(){
      return cov;
    }


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


private:

    int minOverlap;
    int length;
    int nr;
    int* cov;
    int** readStarts;
    int** readEnds;
};

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
  double tmp = best/ (double) seq.size();
  int minOverlap = 1;
  while(tmp > (0.00000001/meanWidth)){
    tmp *= tmp;
    minOverlap++;
  }
  std::cout << "minOverlap ist: " << minOverlap <<"\n";
  return minOverlap;
}

#endif
