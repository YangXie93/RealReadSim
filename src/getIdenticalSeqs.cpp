#include <Rcpp.h>
#include <vector>
#include <list>
#include <iostream>
#include <algorithm>

using namespace Rcpp;


// class representing a contig that is identical on two different genomes. It gets reads that mapped to both genomes in the a cohesive
// region and builds a contig seeing which reads makeup a contiguos identical sequence which is then saved
//
class Contig{
    public:

        void addMany(std::vector<int>::iterator& st1b,std::vector<int>::iterator& st1e,std::vector<int>::iterator& en1b,std::vector<int>::iterator& st2b,std::vector<int>::iterator& en2b){
            while(st1b != st1e){
                if(hasOverlap(*st1b,*en1b)){
                    if(this->add(*st1b,*en1b,*st2b,*en2b)){


                        this->readOut = false;
                    }
                    else{
                        if(!hasNextSeed){
                            nextSeedS1 = st1b;
                            nextSeedE1 = en1b;
                            nextSeedS2 = st2b;
                            nextSeedE2 = en2b;
                            hasNextSeed = true;
                        }
                    }
                    st1b++;
                    en1b++;
                    st2b++;
                    en2b++;
                }
                else{
                    if(!hasNextSeed){
                        nextSeedS1 = st1b;
                        nextSeedE1 = en1b;
                        nextSeedS2 = st2b;
                        nextSeedE2 = en2b;
                        hasNextSeed = true;
                    }
                    break;
                }
            }
        }

        bool add(int st1,int en1,int st2,int en2){
            if(s1+s2+e1+e2 == 0){
                s1 = st1;
                s2 = st2;
                e1 = en1;
                e2 = en2;
                length = e1-s1;
                readOut = true;
            }
            else{
                if(sameOverlap(st1,en1,st2,en2)){
                    readOut = true;
                    if(st1 < s1){
                        length += s1 - st1;
                        s1 = st1;
                        s2 = st2;
                    }
                    else{
                        if(en1 > e1){
                            length += en1 -e1;
                            e1 = en1;
                            e2 = en2;
                        }
                    }
                }
            }
            return readOut;
        }

        bool hasOverlap(int start,int end){
            return (start < e1) && (end > s1);
        }

        bool sameOverlap(int st1,int en1,int st2,int en2){
            bool res = false;
            if((en1 > s1 && st1 < e1) && (en2 > s2 && st2 < e2)){
                if((e1 - en1 == e2 -en2)){
                    res = true;
                }
                else{
                    readOut = true;
                }
            }
            return res;
        }

        int getLength(){
            return length;
        }
        int getS1(){
            return s1;
        }
        int getS2(){
            return s2;
        }
        int getE1(){
            return e1;
        }
        int getE2(){
            return e2;
        }
        std::vector<int>::iterator getNextS1(){
            return nextSeedS1;
        }
        std::vector<int>::iterator getNextS2(){
            return nextSeedS2;
        }
        std::vector<int>::iterator getNextE1(){
            return nextSeedE1;
        }
        std::vector<int>::iterator getNextE2(){
            return nextSeedE2;
        }
    private:
        int s1 = 0;
        int e1 = 0;
        int s2 = 0;
        int e2 = 0;
        int length = 0;
        std::vector<int>::iterator nextSeedS1;
        std::vector<int>::iterator nextSeedE1;
        std::vector<int>::iterator nextSeedS2;
        std::vector<int>::iterator nextSeedE2;
        bool readOut = false;
        bool hasNextSeed = false;
};

// function tying the Contig object to R
//
//[[Rcpp::export]]
List getIdenticalSeqs(std::vector<int>& starts1,std::vector<int>& ends1,std::vector<int>& starts2,std::vector<int>& ends2,std::string nm1,std::string nm2,int minL = 0){


    std::vector<int> contStarts1;
    std::vector<int> contEnds1;
    std::vector<int> contStarts2;
    std::vector<int> contEnds2;
    std::vector<int> contNms;

    std::vector<int>::iterator st1 = starts1.begin();
    std::vector<int>::iterator en1 = ends1.begin();
    std::vector<int>::iterator st2 = starts2.begin();
    std::vector<int>::iterator en2 = ends2.begin();
    std::vector<int>::iterator st1End = starts1.end();

    int n = 1;
    while(distance(st1,st1End) > 0){
        Contig tmp;
        tmp.addMany(st1,st1End,en1,st2,en2);
        if(tmp.getLength() >= minL){
            contStarts1.push_back(tmp.getS1());
            contStarts2.push_back(tmp.getS2());
            contEnds1.push_back(tmp.getE1());
            contEnds2.push_back(tmp.getE2());
            contNms.push_back(n);
            n++;
        }
        st1 =tmp.getNextS1();
        en1 = tmp.getNextE1();
        st2 = tmp.getNextS2();
        en2 = tmp.getNextE2();
    }
    return List::create(contNms,contStarts1,contEnds1,contStarts2,contEnds2,nm1,nm2);
}
//[[Rcpp::export]]
List getIdenticalSeqsList(std::vector<std::string> &names1,std::list<std::vector<int> >& starts1,std::list<std::vector<int> >& ends1,std::vector<std::string> &names2,std::list<std::vector<int> >& starts2,std::list<std::vector<int> >& ends2,int minL = 0){

    List res;

    std::list<std::vector<int> >::iterator s1;
    std::list<std::vector<int> >::iterator s2 = starts2.begin();
    std::list<std::vector<int> >::iterator e1 = ends1.begin();
    std::list<std::vector<int> >::iterator e2 = ends2.begin();
    std::vector<std::string>::iterator nm1 = names1.begin();
    std::vector<std::string>::iterator nm2 = names2.begin();

    int i = 1;//################
    for(s1 = starts1.begin();s1 != starts1.end();s1++){
        res.push_back(getIdenticalSeqs(*s1,*e1,*s2,*e2,*nm1,*nm2,minL));
        s2++;
        e1++;
        e2++;
        nm1++;
        nm2++;

        i++;//################
    }
    return res;
}

