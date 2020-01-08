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
        this->pos = pos;
        this->width = width;
        ids = new int*[(int) pos->size()];
        cov = new int[length];
        for(int i = 0;i < length; i++){
          cov[i] = 0;                           //Die Elemente des Array werden alle auf 0 gesetzt
        }
        construct();
    }

    ~SeqRep(){
        delete[] ids;
        delete[] cov;
    }

    //Funktion zum auswerten des overlaps
    /* Bei jedem Read bei dem Eine Stelle auf die er gemapped wurde >= "coverage" ist werden hier diese Stellen gezaehlt.
       sobald eine folge an Stellen >= "coverage" abbricht wird der bisherige "count" gespeichert.
       Sind dann alle Stellen des Reads abgearbeitet werden die "counts" ueberprueft, ob mindestens einer davon dem
       "minOverlap" entspricht. Ist dies nicht der Fall wird der vorher zugeordnete boolean wert wieder umgedreht.
    */
    void evalOverlap(std::vector<bool>&  reads){    //Das Ergebnis aus evalCov() wird hier eingegeben
      which = &reads;
      int count = 0;
      int size = (int) pos ->size();
      reads.reserve(size);
      int* read;
      int best = 0;
      int n;



      for(int i = 0;i < size;i++){
          read = ids[i];
          for(int j = 0;j < width[0][i];j++){
            n = 0;
            if(pos[0][i] +j >= length){
                n = length;
            }

            if(*(read +j -n) > 1){ //check, ob an der Stelle die coverage stimmt
                count++;
            }
            else{
                if(count > best){
                  best = count;
                }
                count = 0;
              }
          }

          if(count > best){
            best = count;
          }
          if(best >= minOverlap){
              reads.push_back(true);
          }
          else{
              reads.push_back(false);
          }
          count = 0;
          best = 0;
      }
    }

    Rcpp::List assembleTestContigs(int minContigLength){

        for(int i = 0;i < length; i++){
            cov[i] = 0;
        }

        int diff = 0;

        for(int i = 0;i < (int) pos->size();i++){
            if(which[0][i]){
                for(int j = 0;j < width[0][i];j++){
                  if((pos[0][i] +j) < length){
                      cov[pos[0][i] +j]++;
                  }
                  else{
                      cov[pos[0][i] +j -length]++;
                  }
                }
                ids[i] = &cov[pos[0][i]];
                diff++;
            }
        }

        std::cout << (int) pos->size() -diff << " reads were sorted out.\n";

        std::vector<int> starts;
        std::vector<int> ends;
        std::vector<double> covs;
        std::vector<double> covSD;
        double mean;

        int n = 0;
        bool switching = true;
        for(int i = 0;i < length;i++){
          if(*(cov +i) > 0 && switching){
            n = i +1;
            switching = false;
          }
          if(*(cov +i) <= 0 && !switching){
            switching = true;
            if(i -n >= minContigLength){
              starts.push_back(n);
              ends.push_back(i);
              mean = meanCovOnRange(n,i);
              covs.push_back(mean);
              covSD.push_back(SDCovOnRange(n,i,mean));
            }
          }
          if(*(cov +i) > 0 && i == length -1){
            if(i -n >= minContigLength){
              starts.push_back(n);
              ends.push_back(i);
              mean = meanCovOnRange(n,i);
              covs.push_back(mean);
              covSD.push_back(SDCovOnRange(n,i,mean));
            }
          }
        }
        Rcpp::List res = Rcpp::List::create(Rcpp::Named("start") = starts,Rcpp::Named("end") = ends,Rcpp::Named("mean") = covs,Rcpp::Named("sd") = covSD);
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


private:

    void construct(){
        for(int i = 0;i < (int) pos->size();i++){
            for(int j = 0;j < width[0][i];j++){
              if((pos[0][i] +j) < length){
                  cov[pos[0][i] +j]++;
              }
              else{
                  cov[pos[0][i] +j -length]++;
              }
            }
            ids[i] = &cov[pos[0][i]];
        }
    }


    std::vector<bool> *which;
    std::vector<int> *pos;
    std::vector<int> *width;
    int minOverlap;
    int length;
    int* cov;
    int** ids;
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


#endif
