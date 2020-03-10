#include <Rcpp.h>
#include <vector>
#include <list>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <string>


using namespace Rcpp;

const double EulerConstant = std::exp(1.0);
const double minProb = 0.0;

std::vector<int> rS1;
std::vector<int> rE1;
std::vector<int> rS2;
std::vector<int> rE2;
std::vector<std::string> rSq1;
std::vector<std::string> rSq2;
std::vector<std::string> rNm1;
std::vector<std::string> rNm2;
std::vector<std::vector<int> > rCv1;
std::vector<std::vector<int> > rCv2;

bool swtch;

int i = 0;
int j = 0;
int k = 0;

int x = 0;
int y = 0;

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
    std::vector<int> res = {(0+noNeg(as1-c1s))+noNeg(noPos((as1-c1s))-noPos((as2-c2s))),(ae1-c1s-noNeg(ae1-c1e))-noNeg(noPos((c1e-ae1))-noPos((c2e-ae2))),(0+noNeg(as2-c2s))+noNeg(noPos(as2-c2s)-noPos(as1-c1s)),(ae2-c2s-noNeg(ae2-c2e))-noNeg(noPos(c2e-ae2)-noPos(c1e-ae1))};
    // Rcout << 0+noNeg(as1-c1s) << " " << noPos((as1-c1s)) << " " << noPos((as2-c2s)) << std::endl;
    // Rcout << (ae1-c1s-noNeg(c1e-ae1)) << " " << noPos((c1e-ae1)) << " " << noPos((c2e-ae2)) << std::endl;
    // Rcout << (0+noNeg(as2-c2s)) << " " << noPos(as2-c2s) << " " << noPos(as1-c1s) << std::endl;
    // Rcout << (ae2-c2s-noNeg(ae2-c2e)) << " " << noPos(c2e-ae2) << " " << noPos(c1e-ae2) << std::endl;
    return res;
}

double calcPoisProb(int x,double lambda)
{
    int tmp = 1;
    for(int n = 1; n <= x;n++)
    {
        tmp = tmp*n;
    }

    return std::pow(EulerConstant,-lambda)*std::pow(lambda,x)/(double)tmp;
}


void fuseCovs(std::vector<std::vector<int> > covs,int nrOv,int stOv1,int stOv2,int c1s,int c1e,int c2s,int c2e,std::string seq1,std::string seq2,std::string name1,std::string name2,double p2)
{
Rcout << "fuseCovs\n";
    std::vector<int> tmp;
    std::vector<int> tmp1;
    std::string s = "";

    std::vector<int>::iterator s1;
    std::vector<int>::iterator e1;
    std::vector<int>::iterator s2;
    std::vector<int>::iterator e2;
    std::vector<std::string>::iterator sq1;
    std::vector<std::string>::iterator sq2;
    std::vector<std::string>::iterator nm1;
    std::vector<std::string>::iterator nm2;
    std::vector<std::vector<int> >::iterator c1;
    std::vector<std::vector<int> >::iterator c2;

    std::vector<int>* st;
    std::vector<int>* en;
    std::vector<std::vector<int> >* co;
    std::vector<std::string>* se;
    std::vector<std::string>* nm;

    bool overwrite = false;
    int erase = 0;

    if(swtch){
        s1 = rS1.begin();
        e1 = rE1.begin();
        s2 = rS2.begin();
        e2 = rE2.begin();
        sq1 = rSq1.begin();
        sq2 = rSq2.begin();
        c1 = rCv1.begin();
        c2 = rCv2.begin();
        nm1 = rNm1.begin();
        nm2 = rNm2.begin();

        se = &rSq2;
        co = &rCv2;
        st = &rS2;
        en = &rE2;
        nm = &rNm2;
    }
    else{
        s1 = rS2.begin();
        e1 = rE2.begin();
        s2 = rS1.begin();
        e2 = rE1.begin();
        sq1 = rSq2.begin();
        sq2 = rSq1.begin();
        c1 = rCv2.begin();
        c2 = rCv1.begin();
        nm1 = rNm2.begin();
        nm2 = rNm1.begin();

        se = &rSq1;
        co = &rCv1;
        st = &rS1;
        en = &rE1;
        nm = &rNm1;
    }
    // pre-overlap

    if(stOv1 > 0)
    {
        for(int n = 0;n < stOv1;n++)
        {
            tmp.push_back(covs[0][n]);
        }
        *next(s1,x) = c1s;
        s = s + seq1.substr(0,stOv1);
        if(stOv2 > 0)
        {
            for(int n = 0;n < stOv2;n++)
            {
                tmp1.push_back(covs[1][n]);
            }
            *next(s2,y) = c2s;
            *next(e2,y) = c2s+stOv2-1;
            *next(c2,y) = tmp1;
            *next(sq2,y) = seq2.substr(0,stOv2);

            tmp1.clear();
        }
        else{
            erase++;
            overwrite = true;
        }
    }
    else
    {
        if(stOv2 > 0)
        {
            if(p2 > minProb)
            {
                erase++;
                for(int n = 0;n < stOv2;n++)
                {
                    tmp.push_back(covs[1][n]);
                }
                *next(s1,x) = (c1s+stOv1)-(stOv2);
                overwrite = true;
                s = s + seq2.substr(0,stOv2);
            }
            else
            {
                for(int n = 0;n < stOv2;n++)
                {
                    tmp1.push_back(covs[1][n]);
                }

                *next(s2,y) = c2s;
                *next(e2,y) = c2s+stOv2-1;
                *next(c2,y) = tmp1;
                *next(sq2,y) = seq2.substr(0,stOv2);

                tmp1.clear();
            }
        }
        else{
            overwrite = true;
            erase++;
        }
    }


    // overlap section
    for(int n = 0;n < nrOv;n++)
    {
        tmp.push_back(covs[0][stOv1+n]+covs[1][stOv2+n]);
    }
    s = s + seq1.substr(stOv1,nrOv);


    // post-overlap
    if(stOv1+nrOv < covs[0].size())
    {
        for(int n = stOv1+nrOv;n < (covs[0]).size();n++)
        {
            tmp.push_back(covs[0][n]);
        }
        *next(e1,x) = c1e;
        s = s + seq1.substr(stOv1+nrOv,seq1.size()-1);

        if(stOv2+nrOv < covs[1].size())
        {
            for(int n = stOv2+nrOv;n < (covs[1]).size();n++)
            {
                tmp1.push_back(covs[1][n]);
            }
            if(!overwrite){
                co->insert(next(c2,y+1),tmp1);
                st->insert(next(s2,y+1),c2s+stOv2+nrOv);
                en->insert(next(e2,y+1),c2e);
                se->insert(next(sq2,y+1),seq2.substr(stOv2+nrOv,seq2.size()-1));
                nm->insert(next(nm2,y+1),name2);
            }
            else{
                *next(c2,y) =tmp1;
                *next(s2,y) =c2s+stOv2+nrOv;
                *next(e2,y) = c2e;
                *next(sq2,y) =seq2.substr(stOv2+nrOv,seq2.size()-1);
            }

        }
        else{
            erase++;
        }
    }
    else
    {
        if(stOv2+nrOv < covs[1].size())
        {
            if(p2 > minProb)
            {
                erase++;
                for(int n = stOv2+nrOv;n < (covs[1]).size();n++)
                {
                    tmp.push_back(covs[1][n]);
                }
                *next(e1,x) = (c1s+stOv1+nrOv-1)+(covs[1].size()-(stOv2+nrOv));
                s = s + seq2.substr(stOv2+nrOv,seq2.size()-1);
            }
            else
            {
                for(int n = stOv2+nrOv;n < (covs[1]).size();n++)
                {
                    tmp1.push_back(covs[1][n]);
                }
                if(!overwrite){
                    co->insert(next(c2,y+1),tmp1);
                    st->insert(next(s2,y+1),c2s+stOv2+nrOv);
                    en->insert(next(e2,y+1),c2e);
                    se->insert(next(sq2,y+1),seq2.substr(stOv2+nrOv,seq2.size()-1));
                    nm->insert(next(nm2,y+1),name2);
                }
                else{
                    *next(c2,y) =tmp1;
                    *next(s2,y) =c2s+stOv2+nrOv;
                    *next(e2,y) = c2e;
                    *next(sq2,y) =seq2.substr(stOv2+nrOv,seq2.size()-1);
                }
            }
        }
        else{
            erase++;
        }
    }
    // if(stOv2 <= 0 && stOv2+nrOv >= covs[1].size()){
    if(erase >= 2){
        co->erase(next(c2,y));
        st->erase(next(s2,y));
        en->erase(next(e2,y));
        se->erase(next(sq2,y));
        nm->erase(next(nm2,y));
    }

    *next(nm1,x) += "," + name2;
    *next(c1,x) = (tmp);
    *next(sq1,x) = (s);

}


void fuseContigs(int c1s,int c1e,std::vector<int> c1Cov,std::string seq1,int c2s,int c2e,std::vector<int> c2Cov,std::string seq2,int as1,int ae1,int as2,int ae2,std::string name1,std::string name2)
{
    Rcout << "fuseContigs\n";
    std::vector<int> site = translateOverlap(c1s,c1e,c2s,c2e,as1,ae1,as2,ae2);
    std::vector<int> covSum;

    std::vector<int>::iterator s1 = next(c1Cov.begin(),site[0]);
    std::vector<int>::iterator s2 = next(c2Cov.begin(),site[2]);

    int siteLength = site[1]-site[0] +1;

    for(int n = 0; n < siteLength;n++){
        covSum.push_back(*(s1+n)+*(s2+n));
    }

    double lambda1 = std::accumulate(c1Cov.begin(),c1Cov.end(),0.0)/(double)c1Cov.size();
    double lambda2 = std::accumulate(c2Cov.begin(),c2Cov.end(),0.0)/(double)c2Cov.size();


    double probArr1[siteLength];
    double probArr2[siteLength];

    for(int n = 0;n < covSum.size();n++)
    {
        *(probArr1+n) = calcPoisProb(covSum[n],lambda1);
        *(probArr2+n) = calcPoisProb(covSum[n],lambda2);
    }

    double prob1 = std::accumulate(probArr1,(probArr1+siteLength),0.0)/(double)siteLength;
    double prob2 = std::accumulate(probArr2,(probArr2+siteLength),0.0)/(double)siteLength;
Rcout << prob1 << " " << prob2 << std::endl;

    if(prob1 > prob2)
    {
        if(prob1 >= minProb)
        {
            x = i;
            y = j;
            swtch = true;
            std::vector<std::vector<int> > covs = {c1Cov,c2Cov};
            fuseCovs(covs,siteLength,site[0],site[2],c1s,c1e,c2s,c2e,seq1,seq2,name1,name2,prob2);
        }

    }
    else
    {
        if(prob2 >= minProb)
        {
            x = j;
            y = i;
            swtch = false;
            std::vector<std::vector<int> > covs = {c2Cov,c1Cov};
            fuseCovs(covs,siteLength,site[2],site[0],c2s,c2e,c1s,c1e,seq2,seq1,name2,name1,prob1);
        }

    }
    Rcout << "fuseContigsEnd\n";

}

bool hasOverlap(int c1s,int c1e,int c2s,int c2e,int a1s,int a1e,int a2s,int a2e)
{
    return (a1s < c1e && a1e > c1s && a2s < c2e && a2e > c2s && ((c1e-a1s) > (c2s-a2s)) && ((c1s-a1s) < (c2e-a2s)));
}

void fuseSameGenOv(std::vector<int>::iterator s1,std::vector<int>::iterator e1,std::vector<int>::iterator s2,std::vector<int>::iterator e2,std::vector<std::vector<int> >::iterator cov1,std::vector<std::vector<int> >::iterator cov2,std::vector<std::string>::iterator seq1,std::vector<std::string>::iterator seq2,std::vector<std::string>::iterator nm2,bool oneOrTwo){
    Rcout << "fuseSameGenCovs\n";
    int diff = *s2-*s1;
    int n;

    for(n = 0;n + diff < (*cov1).size();n++){
        *next((*cov1).begin(),diff+n) += *next((*cov2).begin(),n);
    }
    for(n = n;n < (*cov2).size();n++){
        (*cov1).push_back(*next((*cov2).begin(),n));
    }
    Rcout << (*seq2).size() << " " << noNeg(*e1-*s2 +1) << " " << *e1 << " " << *s2 << std::endl;
    int pos = noNeg(*e1-*s2 +1);
    if(pos <(*seq2).size()){
        *seq1 += (*seq2).substr(noNeg(*e1-*s2 +1),(*seq2).size()-(noNeg(*e1-*s2)));
    }
    *e1 = *e2;
    if(oneOrTwo){
        rS1.erase(s2);
        rE1.erase(e2);
        rCv1.erase(cov2);
        rSq1.erase(seq2);
        rNm1.erase(nm2);
    }
    else{
        rS2.erase(s2);
        rE2.erase(e2);
        rCv2.erase(cov2);
        rSq2.erase(seq2);
        rNm2.erase(nm2);
    }

}

//[[Rcpp::export]]
List mkChimeras(std::vector<int>& starts1,std::vector<int>& ends1,std::vector<std::vector<int> >& covs1,std::vector<int>& starts2,std::vector<int>& ends2,std::vector<std::vector<int> >& covs2,std::vector<int>& aStarts1,std::vector<int>& aEnds1,std::vector<int>& aStarts2,std::vector<int>& aEnds2,std::vector<std::string>& seqs1,std::vector<std::string>& seqs2,std::vector<std::string>& name1,std::vector<std::string>& name2)
{
    Rcout << "mkChinḿeras\n";
    rS1 = starts1;
    rE1 = ends1;
    rS2 = starts2;
    rE2 = ends2;
    rSq1 = seqs1;
    rSq2 = seqs2;
    rCv1 = covs1;
    rCv2 = covs2;
    rNm1 = name1;
    rNm2 = name2;

    std::vector<int>::iterator s1 = rS1.begin();
    std::vector<int>::iterator e1 = rE1.begin();
    std::vector<int>::iterator s2 = rS2.begin();
    std::vector<int>::iterator e2 = rE2.begin();
    std::vector<int>::iterator as1 = aStarts1.begin();
    std::vector<int>::iterator ae1 = aEnds1.begin();
    std::vector<int>::iterator as2 = aStarts2.begin();
    std::vector<int>::iterator ae2 = aEnds2.begin();
    std::vector<std::string>::iterator sq1 = rSq1.begin();
    std::vector<std::string>::iterator sq2 = rSq2.begin();
    std::vector<std::string>::iterator nm1 = rNm1.begin();
    std::vector<std::string>::iterator nm2 = rNm2.begin();
    std::vector<std::vector<int> >::iterator c1 = rCv1.begin();
    std::vector<std::vector<int> >::iterator c2 = rCv2.begin();

    while(i < rS1.size() && j < rS2.size() && k < aStarts1.size())
    {

        if(hasOverlap(*next(s1,i),*next(e1,i),*next(s2,j),*next(e2,j),*next(as1,k),*next(ae1,k),*next(as2,k),*next(ae2,k)))
        {
            fuseContigs(*next(s1,i),*next(e1,i),*next(c1,i),*next(sq1,i),*next(s2,j),*next(e2,j),*next(c2,j),*next(sq2,j),*next(as1,k),*next(ae1,k),*next(as2,k),*next(ae2,k),*next(nm1,i),*next(nm2,j));

            s1 = rS1.begin();
            e1 = rE1.begin();
            s2 = rS2.begin();
            e2 = rE2.begin();
            as1 = aStarts1.begin();
            ae1 = aEnds1.begin();
            as2 = aStarts2.begin();
            ae2 = aEnds2.begin();
            sq1 = rSq1.begin();
            sq2 = rSq2.begin();
            c1 = rCv1.begin();
            c2 = rCv2.begin();
            nm1 = rNm1.begin();
            nm2 = rNm2.begin();
        }


        if(*next(ae1,k) > *next(e1,i) || *next(ae2,k) > *next(e2,j))
        {
            if(*next(as1,k) > *next(e1,i)){
                i++;
            }
            else{

                if(*next(as2,k) > *next(e2,j)){
                    j++;
                }
                else{
                    if(*next(ae1,k) > *next(e1,i))
                    {
                         i++;
                    }

                    if(*next(ae2,k) > *next(e2,j))
                    {
                        j++;
                    }
                }
            }
        }
        else
        {
            if(*next(ae1,k) == *next(e1,i) || *next(ae2,k) == *next(e2,j))
            {
                if(*next(ae1,k) == *next(e1,i))
                {
                    i++;
                }
                if(*next(ae2,k) == *next(e2,j))
                {
                    j++;
                }
                k++;
            }
            else{
                k++;
            }
        }
    }

    for(int n = 0;n < rS1.size()-1;n++){

        s1 = rS1.begin();
        e1 = rE1.begin();
        sq1 = rSq1.begin();
        c1 = rCv1.begin();
        if(*next(e1,n) >= *next(s1,n+1)-1){
            fuseSameGenOv(next(s1,n),next(e1,n),next(s1,n+1),next(e1,n+1),next(c1,n),next(c1,n+1),next(sq1,n),next(sq1,n+1),next(nm1,n+1),true);
            n--;
        }
    }
    for(int n = 0;n < rS2.size()-1;n++){

        s2 = rS2.begin();
        e2 = rE2.begin();
        sq2 = rSq2.begin();
        c2 = rCv2.begin();

        if(*next(e2,n) >= *next(s2,n+1)-1){
            fuseSameGenOv(next(s2,n),next(e2,n),next(s2,n+1),next(e2,n+1),next(c2,n),next(c2,n+1),next(sq2,n),next(sq2,n+1),next(nm2,n+1),false);
            n--;
        }
    }

    return List::create(rS1,rE1,rSq1,rCv1,rNm1,rS2,rE2,rSq2,rCv2,rNm2);
}


