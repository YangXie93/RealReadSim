#include <Rcpp.h>
#include <vector>
#include <list>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <string>

//[[Rcpp::plugins(cpp14)]]

using namespace Rcpp;

// minimal value for similarity in coverage for the fusion of two conigs
//
double minDist;

// result vectors
//
std::vector<int> rS1;
std::vector<int> rE1;
std::vector<int> rS2;
std::vector<int> rE2;
std::list<std::string> rSq1;
std::list<std::string> rSq2;
std::vector<std::string> rNm1;
std::vector<std::string> rNm2;
std::list<std::vector<int> > rCv1;
std::list<std::vector<int> > rCv2;
std::list<std::vector<int> > rReadNrVec1;
std::list<std::vector<int> > rReadNrVec2;

// variables that determine wether the data of the contig data will be saved or discarded
//
bool save1;
bool save2;

// variable determining wether the contig from genome 1 will be fused to the contig od genome 2
// or the other way around
// bool swtch;

// variables for checking wether a fusion has been noted in the contigs name or not
//
bool addToNm1;
bool addToNm2;

// function returning zero if num is negative or num itself if not
//
int noNeg(int num)
{
    if(num < 0)
    {
        return 0;
    }
    else
    {
        return num;
    }
}

// function returning zero if num is positive or num itself if not
//
int noPos(int num)
{
    if(num > 0)
    {
        return 0;
    }
    else
    {
        return num;
    }
}


// function returning the positoin of a subsequence that is identical on both contigs
//
//[[Rcpp::export]]
std::vector<int> translateOverlap(int c1s,int c1e,int c2s,int c2e,int as1,int ae1,int as2,int ae2)
{
    return {(0+noNeg(as1-c1s))+noNeg(noPos((as1-c1s))-noPos((as2-c2s))),(ae1-c1s-noNeg(ae1-c1e))-noNeg(noPos((c1e-ae1))-noPos((c2e-ae2))),(0+noNeg(as2-c2s))+noNeg(noPos(as2-c2s)-noPos(as1-c1s)),(ae2-c2s-noNeg(ae2-c2e))-noNeg(noPos(c2e-ae2)-noPos(c1e-ae1))};;
}


// function that fuses two contigs. Residues are saved as individual contigs. start, end, sequence and coverage values are modified
//
void fuseCovs(int nrOv,int stOv1,int stOv2,std::vector<int>* tcov1,std::vector<int>* tcov2){

    std::vector<int> tmp;

    for(int n = 0;n < nrOv;n++){
        tmp.push_back((*tcov1)[stOv1+n]+(*tcov2)[stOv2+n]);
    }


    *tcov1 = tmp;
}



// function that fuses the values of the readcount per sample values of two contigs
//
void fuseReadNrVecs(std::vector<int>& readNrVec1, std::vector<int>& readNrVec2){
    std::vector<int> res;
    std::vector<int>::iterator j = readNrVec2.begin();

    for(std::vector<int>::iterator i = readNrVec1.begin();i != readNrVec1.end();i++){
        *(i) += *(j);
        j++;
    }
}

// function stringing together translateOverlap(), fuseCovs, fuseReadNrVecs and calculates the coverage distance value for each contig
// and the fused coverage
//
void fuseContigs(int as1,int ae1,int as2,int ae2,std::string name1,std::string name2,int* ts1,int* te1,int* ts2,int* te2,std::string* tsq1,std::string* tsq2,std::string* tnm1,std::string* tnm2,std::vector<int>* tcov1,std::vector<int>* tcov2,std::vector<int>* trv1,std::vector<int>* trv2)
{
    std::vector<int> site = translateOverlap(*ts1,*te1,*ts2,*te2,as1,ae1,as2,ae2);

    bool fuse1 = (site[0] == 0 && site[1] == (*te1-*ts1));
    bool fuse2 = (site[2] == 0 && site[3] == (*te2-*ts2));

    std::vector<int>::iterator s1 = next((*tcov1).begin(),site[0]);
    std::vector<int>::iterator s2 = next((*tcov2).begin(),site[2]);
    int siteLength = site[1]-site[0] +1;

    double dist1 = 0;
    double dist2 = 0;
    double comb = 0;
    int n;
    for(n = 0; n < siteLength;n++){
        dist1 += *(s1+n);
        dist2 += *(s2+n);
        comb += *(s1+n)+*(s2+n);
    }

    dist1 /= comb;
    dist2 /= comb;



    if(fuse2 || fuse1)
    {

        if( (fuse2 || (fuse1 && fuse2 && tcov2->size() < tcov1->size())) && (dist1 > minDist) ){

            fuseCovs(siteLength,site[0],site[2],tcov1,tcov2);
            fuseReadNrVecs(*trv1,*trv2);
            addToNm1 = true;

            save1 = true;
        }
        if( (fuse1 || (fuse1 && fuse2 && tcov1->size() < tcov2->size())) && (dist2 > minDist) ){

            fuseCovs(siteLength,site[2],site[0],tcov2,tcov1);
            fuseReadNrVecs(*trv2,*trv1);
            addToNm2 = true;

            save2 = true;
        }

    }

}

// function determening wether two contigs have identical sequences or not
//
//[[Rcpp::export]]
bool hasOverlap(int c1s,int c1e,int c2s,int c2e,int a1s,int a1e,int a2s,int a2e)
{
    return (a1s <= c1e && a1e >= c1s && a2s <= c2e && a2e >= c2s && ((c1e-a1s) >= (c2s-a2s)) && ((c1s-a1s) <= (c2e-a2s)));
}

// function to set many iterators forward
//
void plusplus(std::vector<int>::iterator* s,std::vector<int>::iterator* e,std::list<std::string>::iterator* sq,std::vector<std::string>::iterator* nm,std::list<std::vector<int> >::iterator* c, std::list<std::vector<int> >::iterator* rnv){
    (*s)++;
    (*e)++;
    (*sq)++;
    (*nm)++;
    (*c)++;
    (*rnv)++;
}

// function to save values to result vector for genome one
//
void iSave(int s,int e,std::string sq,std::string nm,std::vector<int> c,std::vector<int> rv){
    rS1.push_back(s);
    rE1.push_back(e);
    rSq1.push_back(sq);
    rNm1.push_back(nm);
    rCv1.push_back(c);
    rReadNrVec1.push_back(rv);
}

// function to save values to result vector for genome two
//
void jSave(int s,int e,std::string sq,std::string nm,std::vector<int> c,std::vector<int> rv){
    rS2.push_back(s);
    rE2.push_back(e);
    rSq2.push_back(sq);
    rNm2.push_back(nm);
    rCv2.push_back(c);
    rReadNrVec2.push_back(rv);
}

// "main" function stringing all of the above together by checking each contig from each genome against the known identical sequences
// of the two genomes, hand them over to the other functions for fusion when necessary and then checking the changed contigs for overlaps
// on the same genome and fuse them if necessary
//
//[[Rcpp::export]]
List mkChimeras(std::vector<int>& starts1,std::vector<int>& ends1,std::list<std::vector<int> >& covs1,std::vector<int>& starts2,std::vector<int>& ends2,std::list<std::vector<int> >& covs2,std::vector<int>& aStarts1,std::vector<int>& aEnds1,std::vector<int>& aStarts2,std::vector<int>& aEnds2,std::list<std::string>& seqs1,std::list<std::string>& seqs2,std::vector<std::string>& name1,std::vector<std::string>& name2,std::list<std::vector<int> >& readNrVecs1,std::list<std::vector<int> >& readNrVecs2,double minimalDistance)
{
    minDist = minimalDistance;



    std::vector<int>::iterator s1 = starts1.begin();
    std::vector<int>::iterator e1 = ends1.begin();
    std::vector<int>::iterator s2 = starts2.begin();
    std::vector<int>::iterator e2 = ends2.begin();
    std::vector<int>::iterator as1 = aStarts1.begin();
    std::vector<int>::iterator ae1 = aEnds1.begin();
    std::vector<int>::iterator as2 = aStarts2.begin();
    std::vector<int>::iterator ae2 = aEnds2.begin();
    std::list<std::string>::iterator sq1 = seqs1.begin();
    std::list<std::string>::iterator sq2 = seqs2.begin();
    std::vector<std::string>::iterator nm1 = name1.begin();
    std::vector<std::string>::iterator nm2 = name2.begin();
    std::list<std::vector<int> >::iterator c1 = covs1.begin();
    std::list<std::vector<int> >::iterator c2 = covs2.begin();
    std::list<std::vector<int> >::iterator rnv1 = readNrVecs1.begin();
    std::list<std::vector<int> >::iterator rnv2 = readNrVecs2.begin();

    int tmpS1 = *(starts1.begin());
    int tmpE1 = *(ends1.begin());
    int tmpS2 = *(starts2.begin());
    int tmpE2 = *(ends2.begin());
    std::string tmpSq1 = *(seqs1.begin());
    std::string tmpSq2 = *(seqs2.begin());
    std::string tmpNm1 = *(name1.begin());
    std::string tmpNm2 = *(name2.begin());
    std::vector<int> tmpC1 = *(covs1.begin());
    std::vector<int> tmpC2 = *(covs2.begin());
    std::vector<int> tmpRV1 = *(readNrVecs1.begin());
    std::vector<int> tmpRV2 = *(readNrVecs2.begin());

    int* ts1 = &tmpS1;
    int* te1 = &tmpE1;
    int* ts2 = &tmpS2;
    int* te2 = &tmpE2;
    std::string* tsq1 = &tmpSq1;
    std::string* tsq2 = &tmpSq2;
    std::string* tnm1 = &tmpNm1;
    std::string* tnm2 = &tmpNm2;
    std::vector<int>* tcov1 = &tmpC1;
    std::vector<int>* tcov2 = &tmpC2;
    std::vector<int>* trv1 = &tmpRV1;
    std::vector<int>* trv2 = &tmpRV2;

    // varaiables saving the last checked identical sequence starts and ends on both genomes
    int lastAs1 = 0;
    int lastAe1 = 0;
    int lastAs2 = 0;
    int lastAe2 = 0;

    // variables determening wether to look at at a new contig of one of the genomes next
    // (in case the contigs are longer than the last identical strech) or wether
    // to look at a new identical sequence next (in case it is longer then the contigs)
    bool i,j,k;

    addToNm1 = false;
    addToNm2 = false;
    while(s1 != starts1.end() && s2 != starts2.end() && as1 != aStarts1.end()){
        // loopCounter++;
        save1 = true;
        save2 = true;
        i = false;
        j = false;
        k = false;

        if(hasOverlap(*ts1,*te1,*ts2,*te2,*as1,*ae1,*as2,*ae2)){
            fuseContigs(*as1,*ae1,*as2,*ae2,*nm1,*nm2,ts1,te1,ts2,te2,tsq1,tsq2,tnm1,tnm2,tcov1,tcov2,trv1,trv2);
        }

        if(*te1 >= *ae1 || *te2 >= *ae2){
            if(*te1 >= *ae1 && *te2 >= *ae2){
                k = true;
            }
            if(*te1 <= *ae1){
                i = true;
            }
            if(*te2 <= *ae2){
                j = true;
            }

        }
        else{
            if(*te1 >= *as1 && *te2 >= *as2){
                if(*ae2-*te2 <= *ae1-*te1){
                    i = true;
                }
                if(*ae1-*te1 <= *ae2-*te2){
                    j = true;
                }
            }
            else{
                if(*te1 < *as1){
                    i = true;
                }
                if(*te2 < *as2){
                    j = true;
                }
            }
        }
        //checking wether to save the data for genome one
        if(i || !save1){
            if(save1){

                if(addToNm1){
                    tmpNm1 += ":" +(*nm2);
                }
                iSave(tmpS1,tmpE1,tmpSq1,tmpNm1,tmpC1,tmpRV1);
                addToNm1 = false;

                lastAs1 = *as1;
                lastAe1 = *ae1;
                lastAs2 = *as2;
                lastAe2 = *ae2;
            }
            plusplus(&s1,&e1,&sq1,&nm1,&c1,&rnv1);
            if(s1 != starts1.end()){
                tmpS1 = *s1;
                tmpE1 = *e1;
                tmpSq1 = *sq1;
                tmpC1 = *c1;
                tmpRV1 = *rnv1;
            }
        }
        //checking wether to save the data for genome two
        if((j || !save2) && s1 != starts1.end()){
            if(save2){

                if(addToNm2){
                    tmpNm2 += ":" +(*nm1);
                }
                jSave(tmpS2,tmpE2,tmpSq2,tmpNm2,tmpC2,tmpRV2);
                addToNm2 = false;

                lastAs2 = *as2;
                lastAe2 = *ae2;
                lastAs1 = *as1;
                lastAe1 = *ae1;
            }
            plusplus(&s2,&e2,&sq2,&nm2,&c2,&rnv2);
            if(s2 != starts2.end()){
                tmpS2 = *s2;
                tmpE2 = *e2;
                tmpSq2 = *sq2;
                tmpC2 = *c2;
                tmpRV2 = *rnv2;
            }
        }
        if(k){
            as1++;
            ae1++;
            as2++;
            ae2++;
        }
    }
    // save contigs if they were still stored in the working variables
    //
    if(save1 && !i){

        iSave(tmpS1,tmpE1,tmpSq1,tmpNm1,tmpC1,tmpRV1);

        plusplus(&s1,&e1,&sq1,&nm1,&c1,&rnv1);
    }
    if(save2 && !j){
        jSave(tmpS2,tmpE2,tmpSq2,tmpNm2,tmpC2,tmpRV2);

        plusplus(&s2,&e2,&sq2,&nm2,&c2,&rnv2);
    }
    // save all contigs that beginn after the last identical sequence
    //
    while(distance(s1, starts1.end()) > 0){
        iSave(*s1,*e1,*sq1,*nm1,*c1,*rnv1);

        plusplus(&s1,&e1,&sq1,&nm1,&c1,&rnv1);
    }
    while(distance(s2, starts2.end()) > 0){
        jSave(*s2,*e2,*sq2,*nm2,*c2,*rnv2);
        plusplus(&s2,&e2,&sq2,&nm2,&c2,&rnv2);
    }
    List res = List::create(rS1,rS2,rE1,rE2,rSq1,rSq2,rCv1,rCv2,rNm1,rNm2,rReadNrVec1,rReadNrVec2);

    rS1.clear();
    rS2.clear();
    rE1.clear();
    rE2.clear();
    rSq1.clear();
    rSq2.clear();
    rCv1.clear();
    rCv2.clear();
    rNm1.clear();
    rNm2.clear();
    rReadNrVec1.clear();
    rReadNrVec2.clear();
    return res;
}





















