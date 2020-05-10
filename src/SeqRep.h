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
    SeqRep(int length,int minOverlap,int nrOfSamples,std::vector<int> *pos,std::vector<int> *width,std::vector<int> *sample){

        this->nrOfSamples = nrOfSamples;
        this->length = length;
        this->minOverlap = minOverlap;
        nr = (int) pos->size();
        readStarts = new int[nr];
        readEnds = new int[nr];
        sampleNr = new int[nr];
        cov = new int[length];
        indivCovs = new int*[nrOfSamples];

        for(int i = 0; i < nrOfSamples;i++){
            *(indivCovs+i) = new int[length];
            std::fill(*(indivCovs+i),*(indivCovs+i)+length-1,0);
        }

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

          addToCov(*(indivCovs+*(sampleIt))+(*posIt),*(indivCovs+*(sampleIt))+(*posIt+*widthIt-1),1);
          addToCov(start,end,1);
          *(readStarts+i) = (*posIt);
          *(readEnds+i) = (*posIt)+(*widthIt)-1;
          *(sampleNr+i) = *sampleIt;
          posIt++;
          widthIt++;
        }

        evalOverlap();
    }

    ~SeqRep(){
        delete[] readStarts;
        delete[] readEnds;
        delete[] cov;
        delete[] sampleNr;
        for(int i = 0;i < nrOfSamples;i++){
            delete[] *(indivCovs+i);
        }
        delete[] indivCovs;
    }


    //Funktion zum auswerten des overlaps
    /* Bei jedem Read bei dem Eine Stelle auf die er gemapped wurde >= "coverage" ist werden hier diese Stellen gezaehlt.
       sobald eine folge an Stellen >= "coverage" abbricht wird der bisherige "count" gespeichert.
       Sind dann alle Stellen des Reads abgearbeitet werden die "counts" ueberprueft, ob mindestens einer davon dem
       "minOverlap" entspricht. Ist dies nicht der Fall wird der vorher zugeordnete boolean wert wieder umgedreht.
    */
    void evalOverlap(){    //Das Ergebnis aus evalCov() wird hier eingegeben

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
            addToCov(*(indivCovs+*(sampleNr+i))+*(readStarts+i),*(indivCovs+*(sampleNr+i))+end,-1);
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
        Rcpp::List mCovVecs;
        Rcpp::List covs;
        std::vector<double> meanVec;

        int* covIt = cov;
        int i = 0;
        int en = 1;
        int* it;

        while(*covIt > 0 && i < length){
          covIt++;
        }

        if(covIt == cov+length-1){
          starts.push_back(1);
          ends.push_back(length);
          meanVec = meanCovOnRange(0,length-1);
          mCovVecs.push_back(meanVec);
          mCovs.push_back(compMeanCovOnRange(0,length-1));
          covs.push_back(std::vector<int> (cov,cov+length-1));
        }
        else{

          int n = 0;
          bool switching = true;

          if(covIt == cov){
            it = cov;
            covIt = cov+length;
          }
          else{
            covIt--;
            it = covIt+1;
          }
          i = it-cov;
          en = it-cov+1;

          while(it != covIt){

            if(it >= cov+length){
              it = cov;
              i = 0;
            }

            if(*(it) > 0 && switching){
              n = i;
              switching = false;
            }
            if(*(it) == 0 && !switching){

              switching = true;
              if(i -n >= minContigLength){
                starts.push_back(n+1);
                ends.push_back(i);
                meanVec = meanCovOnRange(n,i-1);
                mCovVecs.push_back(meanVec);
                mCovs.push_back(compMeanCovOnRange(n,i-1));
                covs.push_back(std::vector<int> (cov+n,cov+i));
              }
            }

            it++;
            i++;
            en++;
          }

        }
        Rcpp::List res = Rcpp::List::create(starts, ends, covs, mCovs,mCovVecs);
        return res;
    }

    double compMeanCovOnRange(int start,int end){
        return std::accumulate(cov+start,cov+end,0)/double (end-start);
    }

    std::vector<double> meanCovOnRange(int start,int end){
        std::vector<double> res;
        for(int i = 0;i < nrOfSamples;i++){
            res.push_back(std::accumulate(*(indivCovs+i)+start,*(indivCovs+i)+end,0)/double (end-start));
        }

        return res;
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


    void addToIndivCov(int* start,int* end,int val,int sample){
        int* tmp = start;
        int x = 1;
        while(start != end+1 && tmp != end+1){
            (*tmp) += val;
            tmp++;
            start++;
            x++;
            if(tmp >= *(indivCovs+sample)+length){
                tmp = *(indivCovs+sample);
            }
        }
    }

private:

    int nrOfSamples;
    int minOverlap;
    int length;
    int nr;
    int* cov;
    int** indivCovs;
    int* readStarts;
    int* readEnds;
    int* sampleNr;
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


#endif
