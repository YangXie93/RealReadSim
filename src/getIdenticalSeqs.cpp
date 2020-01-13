#include <Rcpp.h>
#include <vector>
#include <list>
#include <iostream>
#include <algorithm>

using namespace Rcpp;
class Contig{
    public:

        void addMany(std::vector<int>::iterator& st1b,std::vector<int>::iterator& st1e,std::vector<int>::iterator& en1b,std::vector<int>::iterator& st2b,std::vector<int>::iterator& en2b){
            while(st1b != st1e){
                if(this->add(*st1b,*en1b,*st2b,*en2b)){
                    st1b++;
                    en1b++;
                    st2b++;
                    en2b++;

                    this->readOut = false;
                }
                else{
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
    private:
        int s1 = 0;
        int e1 = 0;
        int s2 = 0;
        int e2 = 0;
        int length = 0;
        bool readOut = false;
};

//[[Rcpp::export]]

List getIdenticalSeqs(std::vector<int>& starts1,std::vector<int>& ends1,std::vector<int>& starts2,std::vector<int>& ends2,int minL = 0){

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
    }
    return List::create(contNms,contStarts1,contEnds1,contStarts2,contEnds2);
}
