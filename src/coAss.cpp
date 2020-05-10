#include <Rcpp.h>
#include <vector>
#include <list>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <string>


using namespace Rcpp;

const double EulerConstant = std::exp(1.0);
const double minProb = 0.8;

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
std::list<std::vector<int> > rCvVec1;
std::list<std::vector<int> > rCvVec2;

bool save1;
bool save2;
bool swtch;

bool check;

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
//[[Rcpp::export]]
std::vector<int> translateOverlap(int c1s,int c1e,int c2s,int c2e,int as1,int ae1,int as2,int ae2)
{
    return {(0+noNeg(as1-c1s))+noNeg(noPos((as1-c1s))-noPos((as2-c2s))),(ae1-c1s-noNeg(ae1-c1e))-noNeg(noPos((c1e-ae1))-noPos((c2e-ae2))),(0+noNeg(as2-c2s))+noNeg(noPos(as2-c2s)-noPos(as1-c1s)),(ae2-c2s-noNeg(ae2-c2e))-noNeg(noPos(c2e-ae2)-noPos(c1e-ae1))};;
}


void fuseCovs(int nrOv,int stOv1,int stOv2,std::string name2,double p2,int* ts1,int* te1,int* ts2,int* te2,std::string* tsq1,std::string* tsq2,std::string* tnm1,std::string* tnm2,std::vector<int>* tcov1,std::vector<int>* tcov2){
// Rcout << "fuseCovs\n";


    std::vector<int> tmp;
    std::vector<int> tmp1;
    std::string s = "";

    std::vector<int>* st;
    std::vector<int>* en;
    std::list<std::vector<int> >* co;
    std::list<std::string>* se;
    std::vector<std::string>* nm;

    int frontS;
    int frontE;
    std::string frontSeq;
    std::string frontNm;
    std::vector<int> frontCv;

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

        save = &save2;
        // Rcout << "norm\n";
    }
    else{

        se = &rSq1;
        co = &rCv1;
        st = &rS1;
        en = &rE1;
        nm = &rNm1;

        save = &save1;
        // Rcout << "switched\n";
    }

    // pre-overlap

    if(stOv1 > 0){

            // Rcout << "1\n";
        for(int n = 0;n < stOv1;n++){
            tmp.push_back((*tcov1)[n]);
        }
        s = s + (*tsq1).substr(0,stOv1);
        if(stOv2 > 0){

                // Rcout << "2\n";
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

            // Rcout << "4\n";
        if(stOv2 > 0){

               // Rcout << "5\n";
            if(p2 >= minProb){

                   // Rcout << "6\n";
                for(int n = 0;n < stOv2;n++){
                    tmp.push_back((*tcov2)[n]);
                }
                *ts1 -= (stOv2);
                s = s+ (*tsq2).substr(0,stOv2);
                countSave++;
                newFront = true;
            }
            else{

                    // Rcout << "7\n";
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

            // Rcout << "9\n";
        for(int n = stOv1+nrOv;n < (*tcov1).size();n++){
            tmp.push_back((*tcov1)[n]);
        }
        s = s + (*tsq1).substr(stOv1+nrOv,(*tsq1).size()-1);

        if((*tcov2).size() > stOv2+nrOv){

                // Rcout << "10\n";
            if(hasFront){
                co->push_back(frontCv);
                st->push_back(frontS);
                en->push_back(frontE);
                se->push_back(frontSeq);
                nm->push_back(frontNm);
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
            }
            countSave++;
        }
    }
    else{

        // Rcout << "13\n";
        if((*tcov2).size() > stOv2+nrOv){

            // Rcout << "14\n";
            if(p2 >= minProb){

                // Rcout << "15\n";
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
                }
                countSave++;
            }
            else{
                // Rcout << "16\n";
                if(hasFront){
                    co->push_back(frontCv);
                    st->push_back(frontS);
                    en->push_back(frontE);
                    se->push_back(frontSeq);
                    nm->push_back(frontNm);
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
    *tnm1 += "," + name2;
    *tcov1 = tmp;
    *tsq1 = s;
}


void fuseContigs(int as1,int ae1,int as2,int ae2,std::string name1,std::string name2,int* ts1,int* te1,int* ts2,int* te2,std::string* tsq1,std::string* tsq2,std::string* tnm1,std::string* tnm2,std::vector<int>* tcov1,std::vector<int>* tcov2,std::vector<int>* tcv1,std::vector<int>* tcv2)
{
    std::vector<int> site = translateOverlap(*ts1,*te1,*ts2,*te2,as1,ae1,as2,ae2);

    std::vector<int>::iterator s1 = next((*tcov1).begin(),site[0]);
    std::vector<int>::iterator s2 = next((*tcov2).begin(),site[2]);
    int siteLength = site[1]-site[0] +1;

    double prob1 = 0;
    double prob2 = 0;
    int n;
    for(n = 0; n < siteLength;n++){
        prob1 += *(s1+n)/(double)(*(s1+n)+*(s2+n));
        prob2 += *(s2+n)/(double)(*(s1+n)+*(s2+n));
    }

    prob1 /= n+1;
    prob2 /= n+1;

    Rcout << prob1 << " " << prob2 << std::endl;//################################################

    if(prob1 > prob2)
    {
        if(prob1 >= minProb)
        {
            Rcout << prob1 << std::endl;//################################################
            swtch = true;
            fuseCovs(siteLength,site[0],site[2],name2,prob2,ts1,te1,ts2,te2,tsq1,tsq2,tnm1,tnm2,tcov1,tcov2);
            *tcv1 = fuseCovVecs(*tcv1,*tcv2);
        }

    }
    else
    {
        if(prob2 >= minProb)
        {
            Rcout << prob2 << std::endl;//################################################
            swtch = false;
            fuseCovs(siteLength,site[2],site[0],name1,prob1,ts2,te2,ts1,te1,tsq2,tsq1,tnm2,tnm1,tcov2,tcov1);
            *tcv2 = fuseCovVecs(*tcv1,*tcv2); 
        }

    }

}
//[[Rcpp::export]]
bool hasOverlap(int c1s,int c1e,int c2s,int c2e,int a1s,int a1e,int a2s,int a2e)
{
    return (a1s <= c1e && a1e >= c1s && a2s <= c2e && a2e >= c2s && ((c1e-a1s) >= (c2s-a2s)) && ((c1s-a1s) <= (c2e-a2s)));
}

bool hasSameGenOverlap(int e1,int s2,int as,int ae){
    return (e1 > s2 && e1 < ae && e1 >= as && s2 > as && s2 <= ae);
}

void fuseSameGenCont(std::vector<int>::iterator s,std::vector<int>::iterator e,std::list<std::vector<int> >::iterator cov,std::list<std::string>::iterator seq,int s2,int e2,std::vector<int> cov2,std::string seq2){
    bool whichS = s2 > *s;
    bool whichE = e2 > *e;
    int pos;
    int diff;
    std::vector<int>::iterator c1 = (*cov).begin();
    std::vector<int>::iterator c2 = cov2.begin();

    std::vector<int> tmpC;
    std::string tmpS = "";

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

std::vector<int> fuseCovVecs(std::vector<int>& covVec1; std::vector<int>& covVec2){
    std::vector<int> res;
    std::vector<int>::iterator j = covVec2.begin();

    for(std::vector<int>::iterator i = covVec1.begin();i != covVec1.end();i++){
        res.push_back(((*i)+(*j))/2);
        j++;
    }

    return res;
}

//[[Rcpp::export]]
List mkChimeras(std::vector<int>& starts1,std::vector<int>& ends1,std::list<std::vector<int> >& covs1,std::vector<int>& starts2,std::vector<int>& ends2,std::list<std::vector<int> >& covs2,std::vector<int>& aStarts1,std::vector<int>& aEnds1,std::vector<int>& aStarts2,std::vector<int>& aEnds2,std::list<std::string>& seqs1,std::list<std::string>& seqs2,std::vector<std::string>& name1,std::vector<std::string>& name2,std::list<std::vector<int> > covVecs1,std::list<std::vector<int> > covVecs2)
{

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
    std::list<std::vector<int> >::iterator cV1 = covVecs1.begin();
    std::list<std::vector<int> >::iterator cV2 = covVecs2.begin();

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
    std::vector<int> tmpCV1 = *(covVecs1.begin());
    std::vector<int> tmpCV2 = *(covVecs2.begin());

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
    std::vector<int>* tcv1 = &tmpCV1;
    std::vector<int>* tcv2 = &tmpCV2;

    int lastAs1 = 0;
    int lastAe1 = 0;
    int lastAs2 = 0;
    int lastAe2 = 0;

    bool i,j,k;

    while(s1 != starts1.end() && s2 != starts2.end() && as1 != aStarts1.end()){
        save1 = true;
        save2 = true;
        i = false;
        j = false;
        k = false;

        check = false;


        if(hasOverlap(*ts1,*te1,*ts2,*te2,*as1,*ae1,*as2,*ae2)){
            fuseContigs(*as1,*ae1,*as2,*ae2,*nm1,*nm2,ts1,te1,ts2,te2,tsq1,tsq2,tnm1,tnm2,tcov1,tcov2);
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

        if(i || !save1){
            if(save1){
               if(s1 != starts1.begin() && hasSameGenOverlap(rE1.back(),tmpS1,lastAs1,lastAe1)){
                    fuseSameGenCont(prev(rS1.end()),prev(rE1.end()),prev(rCv1.end()),prev(rSq1.end()),tmpS1,tmpE1,tmpC1,tmpSq1);
                }
                else{
                    rS1.push_back(tmpS1);
                    rE1.push_back(tmpE1);
                    rSq1.push_back(tmpSq1);
                    rNm1.push_back(tmpNm1);
                    rCv1.push_back(tmpC1);
                }
                lastAs1 = *as1;
                lastAe1 = *ae1;
                lastAs2 = *as2;
                lastAe2 = *ae2;
            }
            s1++;
            e1++;
            sq1++;
            nm1++;
            c1++;
            cV1++;
            if(s1 != starts1.end()){
                tmpS1 = *s1;
                tmpE1 = *e1;
                tmpSq1 = *sq1;
                tmpC1 = *c1;
            }
        }
        if(j || !save2){
            if(save2){
                if(s2 != starts2.begin() && hasSameGenOverlap(*prev(rE2.end()),tmpS2,lastAs2,lastAe2)){
                    fuseSameGenCont(prev(rS2.end()),prev(rE2.end()),prev(rCv2.end()),prev(rSq2.end()),tmpS2,tmpE2,tmpC2,tmpSq2);
                }
                else{
                    rS2.push_back(tmpS2);
                    rE2.push_back(tmpE2);
                    rSq2.push_back(tmpSq2);
                    rNm2.push_back(tmpNm2);
                    rCv2.push_back(tmpC2);
                }
                lastAs2 = *as2;
                lastAe2 = *ae2;
                lastAs1 = *as1;
                lastAe1 = *ae1;
            }
            s2++;
            e2++;
            sq2++;
            nm2++;
            c2++;
            cV2++;
            if(s2 != starts2.end()){
                tmpS2 = *s2;
                tmpE2 = *e2;
                tmpSq2 = *sq2;
                tmpC2 = *c2;
            }
        }
        if(k){
            as1++;
            ae1++;
            as2++;
            ae2++;
        }
    }

    if((k && !j && !i) || (!k && j && !i) || (!k && !j && i)){
        if(save1 && !i){
            if(rE1.size() > 0 && hasSameGenOverlap(rE1.back(),tmpS1,lastAs1,lastAe1)){
                fuseSameGenCont(prev(rS1.end()),prev(rE1.end()),prev(rCv1.end()),prev(rSq1.end()),tmpS1,tmpE1,tmpC1,tmpSq1);
            }
            else{
                rS1.push_back(tmpS1);
                rE1.push_back(tmpE1);
                rSq1.push_back(tmpSq1);
                rNm1.push_back(tmpNm1);
                rCv1.push_back(tmpC1);
            }
            s1++;
            e1++;
            sq1++;
            nm1++;
            c1++;
        }
        if(save2 && !j){
            if(rE2.size() > 0 && hasSameGenOverlap(rE2.back(),tmpS2,lastAs2,lastAe2)){
                fuseSameGenCont(prev(rS2.end()),prev(rE2.end()),prev(rCv2.end()),prev(rSq2.end()),tmpS2,tmpE2,tmpC2,tmpSq2);
            }
            else{
                rS2.push_back(tmpS2);
                rE2.push_back(tmpE2);
                rSq2.push_back(tmpSq2);
                rNm2.push_back(tmpNm2);
                rCv2.push_back(tmpC2);
            }
            s2++;
            e2++;
            sq2++;
            nm2++;
            c2++;
        }
    }


    while(distance(s1, starts1.end()) > 0){
        if(rE1.size() > 0 && hasSameGenOverlap(rE1.back(),*s1,lastAs1,lastAe1)){
            fuseSameGenCont(prev(rS1.end()),prev(rE1.end()),prev(rCv1.end()),prev(rSq1.end()),*s1,*e1,*c1,*sq1);
        }
        else{
            rS1.push_back(*s1);
            rE1.push_back(*e1);
            rSq1.push_back(*sq1);
            rNm1.push_back(*nm1);
            rCv1.push_back(*c1);
        }
        s1++;
        e1++;
        sq1++;
        nm1++;
        c1++;
    }
    while(distance(s2, starts2.end()) > 0){
        Rcout << distance(s2, starts2.end()) << " " << hasSameGenOverlap(rE2.back(),*s2,lastAs2,lastAe2) << std::endl;
        if(rE2.size() > 0 && hasSameGenOverlap(rE2.back(),*s2,lastAs2,lastAe2)){
            fuseSameGenCont(prev(rS2.end()),prev(rE2.end()),prev(rCv2.end()),prev(rSq2.end()),*s2,*e2,*c2,*sq2);
        }
        else{
            rS2.push_back(*s2);
            rE2.push_back(*e2);
            rSq2.push_back(*sq2);
            rNm2.push_back(*nm2);
            rCv2.push_back(*c2);
        }
        s2++;
        e2++;
        sq2++;
        nm2++;
        c2++;
    }
    return List::create(rS1,rS2,rE1,rE2,rSq2,rSq1,rCv1,rCv2,rNm1,rNm2,rCvVec1,rCvVec2);
}


