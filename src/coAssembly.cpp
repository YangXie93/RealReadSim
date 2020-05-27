#include <Rcpp.h>
#include <vector>
#include <list>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <string>

//[[Rcpp::plugins(cpp11)]]

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
bool swtch;

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
// function to calculate the aproximate per sample read count for a given per sample read count vector
//
std::vector<int> partialReadNrVec(std::vector<int> *readNrVec, int length, int lengthPart){
    double share = lengthPart/(double)length;
    double rest = 1.0 -share;
    std::vector<int> res;
    std::vector<int>::iterator rnvIt;
    for(rnvIt = readNrVec->begin();rnvIt != readNrVec->end();rnvIt++){
        res.push_back((int) ((*rnvIt) * share));
        *rnvIt = (int) ((*rnvIt) * rest);
    }
    return res;
}

// function that fuses two contigs. Residues are saved as individual contigs. start, end, sequence and coverage values are modified
//
void fuseCovs(int nrOv,int stOv1,int stOv2,std::string name2,double p2,int* ts1,int* te1,int* ts2,int* te2,std::string* tsq1,std::string* tsq2,std::string* tnm1,std::string* tnm2,std::vector<int>* tcov1,std::vector<int>* tcov2,std::vector<int>* trv1,std::vector<int>* trv2){

    std::vector<int> tmp;
    std::vector<int> tmp1;
    std::string s = "";

    std::vector<int>* st;
    std::vector<int>* en;
    std::list<std::vector<int> >* co;
    std::list<std::string>* se;
    std::vector<std::string>* nm;
    std::list<std::vector<int> >* rv;

    int frontS;
    int frontE;
    std::string frontSeq;
    std::string frontNm;
    std::vector<int> frontCv;

    int length;

    bool* save;
    bool hasFront = false;
    bool newFront = false;
    int countSave = 0;

    if(swtch){

        se = &rSq2;
        co = &rCv2;
        st = &rS2;
        en = &rE2;
        nm = &rNm2;
        rv = &rReadNrVec2;

        save = &save2;

        length = tsq2->size();
    }
    else{

        se = &rSq1;
        co = &rCv1;
        st = &rS1;
        en = &rE1;
        nm = &rNm1;
        rv = &rReadNrVec2;

        save = &save1;

        length = tsq1->size();
    }

    // preoverlap section
    if(stOv1 > 0){

        for(int n = 0;n < stOv1;n++){
            tmp.push_back((*tcov1)[n]);
        }
        s = s + (*tsq1).substr(0,stOv1);
        if(stOv2 > 0){

            for(int n = 0;n < stOv2;n++){
                frontCv.push_back((*tcov2)[n]);
            }
            frontS = *ts2;
            frontE = *ts2+stOv2-1;
            frontSeq = (*tsq2).substr(0,stOv2);
            frontNm = name2;
            hasFront = true;
        }
        else{
            countSave++;
        }
    }
    else{

        if(stOv2 > 0){

            if(p2 >= minDist){

                for(int n = 0;n < stOv2;n++){
                    tmp.push_back((*tcov2)[n]);
                }
                *ts1 -= (stOv2);
                s = s+ (*tsq2).substr(0,stOv2);
                countSave++;
                newFront = true;
            }
            else{

                for(int n = 0;n < stOv2;n++){
                    frontCv.push_back((*tcov2)[n]);
                }
                frontS = *ts2;
                frontE = *ts2+stOv2-1;
                frontSeq = (*tsq2).substr(0,stOv2);
                frontNm = name2;
                hasFront = true;
                tmp1.clear();
            }
        }
        else{
            countSave++;
        }
    }


    // overlap section
    for(int n = 0;n < nrOv;n++){
        tmp.push_back((*tcov1)[stOv1+n]+(*tcov2)[stOv2+n]);
    }
    s = s + (*tsq1).substr(stOv1,nrOv);


    // post-overlap
    if(stOv1+nrOv < (*tcov1).size()){

        for(int n = stOv1+nrOv;n < (*tcov1).size();n++){
            tmp.push_back((*tcov1)[n]);
        }
        s = s + (*tsq1).substr(stOv1+nrOv,(*tsq1).size()-1);

        if((*tcov2).size() > stOv2+nrOv){

            if(hasFront){
                co->push_back(frontCv);
                st->push_back(frontS);
                en->push_back(frontE);
                se->push_back(frontSeq);
                nm->push_back(frontNm);
                rv->push_back(partialReadNrVec(trv2,length,frontSeq.size()));
            }
            *tsq2 = (*tsq2).substr(stOv2+nrOv-1,(*tsq2).size()-(stOv2+nrOv));
            tmp1 = std::vector<int> (next((*tcov2).begin(),stOv2+nrOv),(*tcov2).end());
            *tcov2 = tmp1;
            *ts2 += stOv2+nrOv;
        }
        else{
            if(hasFront){
                *ts2 = frontS;
                *te2 = frontE;
                *tsq2 = frontSeq;
                *tnm2 = frontNm;
                *tcov2 = frontCv;
                partialReadNrVec(trv2,length,frontSeq.size());
            }
            countSave++;
        }
    }
    else{

        if((*tcov2).size() > stOv2+nrOv){

            if(p2 >= minDist){

                s += (*tsq2).substr(stOv2+nrOv,(*tsq2).size()-1);
                for(int n = stOv2+nrOv;n < (*tcov2).size();n++){
                    tmp.push_back((*tcov2)[n]);
                }
                int x = stOv1;
                if(newFront){
                    x = stOv2;
                }
                *te1 = ((*ts1)+x+nrOv-1)+(((*te2)-(*ts2)+1) -(stOv2+nrOv));

                if(hasFront){
                    *ts2 = frontS;
                    *te2 = frontE;
                    *tsq2 = frontSeq;
                    *tnm2 = frontNm;
                    *tcov2 = frontCv;
                    partialReadNrVec(trv2,length,frontSeq.size());
                }
                countSave++;
            }
            else{
                if(hasFront){
                    co->push_back(frontCv);
                    st->push_back(frontS);
                    en->push_back(frontE);
                    se->push_back(frontSeq);
                    nm->push_back(frontNm);
                    rv->push_back(partialReadNrVec(trv2,length,frontSeq.size()));
                }
                *tsq2 = (*tsq2).substr(stOv2+nrOv-1,(*tsq2).size()-(stOv2+nrOv));
                tmp1 = std::vector<int> (next((*tcov2).begin(),stOv2+nrOv),(*tcov2).end());
                *tcov2 = tmp1;
                *ts2 += stOv2+nrOv;
            }
        }
        else{
            if(hasFront){
                *ts2 = frontS;
                *te2 = frontE;
                *tsq2 = frontSeq;
                *tnm2 = frontNm;
                *tcov2 = frontCv;
            }
            countSave++;
        }
    }
    if(countSave >= 2){
        *save = false;
    }
    *tcov1 = tmp;
    *tsq1 = s;
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

    bool fuse = (site[0] == 0 && site[1] == (*te1-*ts1)) || (site[2] == 0 && site[3] == (*te2-*ts2));

    std::vector<int>::iterator s1 = next((*tcov1).begin(),site[0]);
    std::vector<int>::iterator s2 = next((*tcov2).begin(),site[2]);
    int siteLength = site[1]-site[0] +1;

    double dist1 = 0;
    double dist2 = 0;
    int n;
    for(n = 0; n < siteLength;n++){
        dist1 += *(s1+n)/(double)(*(s1+n)+*(s2+n));
        dist2 += *(s2+n)/(double)(*(s1+n)+*(s2+n));
    }

    dist1 /= n+1;
    dist2 /= n+1;

    // the worse fitting contig will be fused to the better fitting one
    if(dist1 > dist2)
    {
        if(dist1 >= minDist || fuse)
        {
            swtch = true;
            fuseCovs(siteLength,site[0],site[2],name2,dist2,ts1,te1,ts2,te2,tsq1,tsq2,tnm1,tnm2,tcov1,tcov2,trv1,trv2);
            fuseReadNrVecs(*trv1,*trv2);
            addToNm1 = true;
        }

    }
    else
    {
        if(dist2 >= minDist || fuse)
        {
            swtch = false;
            fuseCovs(siteLength,site[2],site[0],name1,dist1,ts2,te2,ts1,te1,tsq2,tsq1,tnm2,tnm1,tcov2,tcov1,trv2,trv1);
            fuseReadNrVecs(*trv2,*trv1);
            addToNm2 = true;
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

// function determining wether two contigs from the same genome overlapp or not (may arise after fusion)
//
bool hasSameGenOverlap(int e1,int s2,int as,int ae){
    return (e1 > s2 && e1 < ae && e1 >= as && s2 > as && s2 <= ae);
}

// function to fuse two contigs from the same genome that overlapp
//
void fuseSameGenCont(std::vector<int>::iterator s,std::vector<int>::iterator e,std::list<std::vector<int> >::iterator cov,std::list<std::string>::iterator seq,std::list<std::vector<int> >::iterator rnv,int s2,int e2,std::vector<int> cov2,std::string seq2,std::vector<int> rnv2){



    fuseReadNrVecs(*rnv,rnv2);

    bool whichS = s2 > *s;
    bool whichE = e2 > *e;
    int pos;
    int diff;
    std::vector<int>::iterator c1 = (*cov).begin();
    std::vector<int>::iterator c2 = cov2.begin();

    std::vector<int> tmpC;
    std::string tmpS = "";

    // determening which start to use
    if(whichS){
        pos = s2-*s;

        for(c1 = c1;c1 != next((*cov).begin(),pos);c1++){
            tmpC.push_back(*c1);
        }
        for(c1 = c1;c1 != (*cov).end() && c2 != cov2.end();c1++){
            tmpC.push_back(*c1 += *c2);
            c2++;
        }
        tmpS += (*seq).substr(0,distance((*cov).begin(),c1));
    }
    else{
        pos = *s-s2;
        for(c2 = c2;c2 != next(cov2.begin(),pos);c2++){
            tmpC.push_back(*c2);
        }
        for(c1 = c1;c1 != (*cov).end() && c2 != cov2.end();c1++){
            tmpC.push_back(*c1 += *c2);
            c2++;
        }
        tmpS += seq2.substr(0,distance(cov2.begin(),c2));
        *s = s2;
    }

    // determening which end to use
    if(whichE){
        diff = distance(cov2.begin(),c2);
        while(c2 != cov2.end()){
            tmpC.push_back(*c2);
            c2++;
        }
        if(diff < seq2.size()){
            tmpS += seq2.substr(diff);
        }
        *e = e2;
    }
    else{
        diff = distance((*cov).begin(),c1);
        while(c1 != (*cov).end()){
            tmpC.push_back(*c1);
            c1++;
        }
        if(diff < (*seq).size()){
            tmpS += (*seq).substr(diff);
        }
    }

    *cov = tmpC;
    *seq = tmpS;
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

    while(s1 != starts1.end() && s2 != starts2.end() && as1 != aStarts1.end()){
        save1 = true;
        save2 = true;
        i = false;
        j = false;
        k = false;

        addToNm1 = false;
        addToNm2 = false;

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
               if(s1 != starts1.begin() && hasSameGenOverlap(rE1.back(),tmpS1,lastAs1,lastAe1)){
                    fuseSameGenCont(prev(rS1.end()),prev(rE1.end()),prev(rCv1.end()),prev(rSq1.end()),prev(rReadNrVec1.end()),tmpS1,tmpE1,tmpC1,tmpSq1,tmpRV1);
                }
                else{
                    if(addToNm1){
                        tmpNm1 += ":" +tmpNm2;
                    }
                    iSave(tmpS1,tmpE1,tmpSq1,tmpNm1,tmpC1,tmpRV1);
                }
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
        if(j || !save2){
            if(save2){
                if(s2 != starts2.begin() && hasSameGenOverlap(*prev(rE2.end()),tmpS2,lastAs2,lastAe2)){
                    fuseSameGenCont(prev(rS2.end()),prev(rE2.end()),prev(rCv2.end()),prev(rSq2.end()),prev(rReadNrVec2.end()),tmpS2,tmpE2,tmpC2,tmpSq2,tmpRV2);
                }
                else{
                    if(addToNm2){
                        tmpNm2 += ":" +tmpNm1;
                    }
                    jSave(tmpS2,tmpE2,tmpSq2,tmpNm2,tmpC2,tmpRV2);
                }
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
        if(rE1.size() > 0 && hasSameGenOverlap(rE1.back(),tmpS1,lastAs1,lastAe1)){
                fuseSameGenCont(prev(rS1.end()),prev(rE1.end()),prev(rCv1.end()),prev(rSq1.end()),prev(rReadNrVec1.end()),tmpS1,tmpE1,tmpC1,tmpSq1,tmpRV1);
        }
        else{
            iSave(tmpS1,tmpE1,tmpSq1,tmpNm1,tmpC1,tmpRV1);

        }
        plusplus(&s1,&e1,&sq1,&nm1,&c1,&rnv1);
    }

    if(save2 && !j){
        if(rE2.size() > 0 && hasSameGenOverlap(rE2.back(),tmpS2,lastAs2,lastAe2)){
            fuseSameGenCont(prev(rS2.end()),prev(rE2.end()),prev(rCv2.end()),prev(rSq2.end()),prev(rReadNrVec2.end()),tmpS2,tmpE2,tmpC2,tmpSq2,tmpRV2);
        }
        else{
            jSave(tmpS2,tmpE2,tmpSq2,tmpNm2,tmpC2,tmpRV2);
        }
        plusplus(&s2,&e2,&sq2,&nm2,&c2,&rnv2);
    }

    // save all contigs that beginn after the last identical sequence
    //
    while(distance(s1, starts1.end()) > 0){
        if(rE1.size() > 0 && hasSameGenOverlap(rE1.back(),*s1,lastAs1,lastAe1)){
            fuseSameGenCont(prev(rS1.end()),prev(rE1.end()),prev(rCv1.end()),prev(rSq1.end()),prev(rReadNrVec1.end()),*s1,*e1,*c1,*sq1,*rnv1);
        }
        else{
            iSave(*s1,*e1,*sq1,*nm1,*c1,*rnv1);
        }
        plusplus(&s1,&e1,&sq1,&nm1,&c1,&rnv1);
    }

    while(distance(s2, starts2.end()) > 0){
        if(rE2.size() > 0 && hasSameGenOverlap(rE2.back(),*s2,lastAs2,lastAe2)){
            fuseSameGenCont(prev(rS2.end()),prev(rE2.end()),prev(rCv2.end()),prev(rSq2.end()),prev(rReadNrVec2.end()),*s2,*e2,*c2,*sq2,*rnv2);
        }
        else{
            jSave(*s2,*e2,*sq2,*nm2,*c2,*rnv2);
        }
        plusplus(&s2,&e2,&sq2,&nm2,&c2,&rnv2);
    }


    return List::create(rS1,rS2,rE1,rE2,rSq2,rSq1,rCv1,rCv2,rNm1,rNm2,rReadNrVec1,rReadNrVec2);
}


