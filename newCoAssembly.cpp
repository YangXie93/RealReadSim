#include <Rcpp.h>
#include <vector>
#include <list>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <string>

//[[Rcpp::plugins(cpp14)]]

using namespace Rcpp;

class SubContig{
public:
    SubContig(std::string* superSeq,std::vector<int>* superCov,int seqID,int start,int end,std::string seqName, std::vector<int>::iterator covStart,std::vector<int>::iterator covEnd,std::string::iterator seqStart,std::string::iterator seqEnd){
        this->superSeq = superSeq;
        this->superCov = superCov;
        this-> seqID = seqID;
        this-> start = start;
        this-> end = end;
        this->seqName = seqName;
        this-> covStart = covStart;
        this-> covEnd = covEnd;
        this-> seqStart = seqStart;
        this-> seqEnd = seqEnd;
    }


    void update(int startDiff,std::string* newSuperSeq,std::vector<int>* newSuperCov){

        int s = distance(superSeq->begin(),seqStart) + startDiff +1;
        int c = distance(superCov->begin(),covStart) + startDiff +1;
        int se = distance(superSeq->begin(),seqEnd) + startDiff +1;
        int ce = distance(superCov->begin(),covEnd) + startDiff +1;

        covStart = next(newSuperCov->begin(),c);
        covEnd = next(newSuperCov->begin(),ce);
        seqStart = next(newSuperSeq->begin(),s);
        seqEnd = next(newSuperSeq->begin(),se);

        this->superCov = newSuperCov;
        this->superSeq = newSuperSeq;


    }

    bool hasOverlap(int start,int end){
        return (this->start < end && this->end > start);
    }

    std::vector<int>::iterator at(int i){
        return next(covStart,i);
    }

    bool isEnd(){
        return seqStart == superSeq->begin() || seqEnd == prev(superSeq->end());
    }

    int seqID;
    int start;
    int end;
    std::string seqName;
    std::vector<int>::iterator covStart;
    std::vector<int>::iterator covEnd;
    std::string::iterator seqStart;
    std::string::iterator seqEnd;
    std::string* superSeq;
    std::vector<int>* superCov;
};

class ContigContainer{
public:

    void add(std::string seq,std::vector<int> covVec,bool isChimeric){
        seqs.push_back(seq);
        covVecs.push_back(covVec);
        this->isChimeric.push_back(isChimeric);
    }

    List finalize(){
        return List::create(seqs,covVecs,isChimeric);
    }

private:
    std::vector<bool> isChimeric;
    std::vector<std::string> seqs;
    std::vector<std::vector<int> > covVecs;
};

class Contig{
public:
    Contig(){
        ;
    }

    Contig(int start,int end, int seqID,int mor,std::string seq,std::string seqName,std::vector<int> cov,std::vector<int> covVec,ContigContainer* cc,double mcs){
        this->cov = cov;
        this->covVec = covVec;
        this->seq = seq;
        this->minOverlap = mor;
        this->seqName = seqName;
        SubContig tmp = SubContig(&seq,&cov,seqID,start,end,seqName,cov.begin(),cov.end(),seq.begin(),seq.end());
        this->subs.push_back(&tmp);
        container = cc;
        minCovShare = mcs;
    }
    ~Contig(){
        if(save)
            container->add(seq,covVec,isChimeric);
    }

    // static


    static bool hasOverlap(int start1,int end1,int start2,int end2,int a1s,int a1e,int a2s,int a2e)
    {
        return (a1s <= end1 && a1e >= start1 && a2s <= end2 && a2e >= start2 &&
                ((end1-a1s) >= (start2-a2s)) && ((start1-a1s) <= (end2-a2s)));
    }

    static std::vector<int> translateOverlap(int start1,int end1,int start2,int end2,int ovStart1,int ovEnd1,int ovStart2,int ovEnd2)
    {
        return {(0 + noNeg( ovStart1 - start1 )) + noNeg(noPos(( ovStart1 - start1 )) - noPos(( ovStart2 - start2 ))),
                (ovEnd1 - start1 - noNeg( ovEnd1 - end1 )) - noNeg( noPos(( end1 - ovEnd1 )) - noPos(( end2 - ovEnd2 ))),
                (0 + noNeg( ovStart2 - start2 )) + noNeg( noPos( ovStart2 - start2 ) - noPos( ovStart1 - start1 )),
                (ovEnd2 - start2 - noNeg( ovEnd2 - end2 )) - noNeg( noPos( end2 - ovEnd2 )-noPos( end1 - ovEnd1 ))};;
    }

    // instance


    double getMinCovShare(){
        return minCovShare;
    }

    Contig* fuse(Contig* other,SubContig* thisSub,SubContig* otherSub,int ovStart1,int ovEnd1,int ovStart2,int ovEnd2){
        int start1 = thisSub->start;
        int end1 = thisSub->end;
        int start2 = otherSub->start;
        int end2 = otherSub->end;


        if(hasOverlap(start1,end1,start2,end2,ovStart1,ovEnd1,ovStart2,ovEnd2)){

            std::vector<int> site = translateOverlap(start1,end1,start2,end2,ovStart2,ovEnd1,ovStart2,ovEnd2);
            bool thisLeft = distance(cov.begin(),thisSub->at(site[0])) > 0;
            bool thisRight = distance(thisSub->at(site[1]),cov.end()) > 0;
            bool otherLeft = distance(other->getCovAt(0),otherSub->at(site[2])) > 0;
            bool otherRight = distance(otherSub->at(site[3]),other->getCovEnd()) > 0;


            if((thisLeft && otherLeft) || (thisRight && otherRight)){



                return other;
            }
            else{

                std::vector<int>::iterator s1 = this->getCovAt(site[0]);
                std::vector<int>::iterator s2 = other->getCovAt(site[2]);
                int siteLength = site[1]-site[0] +1;

                double share1 = 0;
                double share2 = 0;
                double comb = 0;
                int n;
                for(n = 0; n < siteLength;n++){
                    share1 += *(s1);
                    share2 += *(s2);
                    comb += *(s1)+*(s2);
                }

                share1 /= comb;
                share2 /= comb;
                if((share1 > minCovShare || share2 > minCovShare)){
                    return combine(*other,site);
                }
                else{
                    return other;
                }
            }
        }
        else{
            return other;
        }
    }

    std::string getSeq(){
        return seq;
    }
    std::vector<int> getCovVec(){
        return covVec;
    }

    std::vector<int>::iterator getCovVecIt(){
        return covVec.begin();
    }

    int length(){
        return seq.size();
    }

    std::vector<int>::iterator getCovAt(int pos){
        return next(cov.begin(),pos);
    }

    int getSubCount(){
        return subCount;
    }

    std::vector<int>::iterator getCovEnd(){
        return cov.end();
    }

    std::vector<SubContig*>::iterator getSubIt(){
        return subs.begin();
    }

    SubContig* getSubAt(int i){
        return subs[i];
    }

    void setSave(bool is){
        save = is;
    }

private:

    std::vector<int> cov;
    std::vector<int> covVec;
    std::vector<int> subKey;
    std::string seqName;
    std::string seq;
    int minOverlap;
    bool save = true;
    bool isChimeric = false;
    int subCount = 1;

    std::vector<SubContig*> subs;

    ContigContainer* container;
    double minCovShare;

    // function returning zero if num is negative or num itself if not
    //
    static int noNeg(int num)
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
    static int noPos(int num)
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


    Contig* combine(Contig other,std::vector<int> identical){


        std::vector<int>::iterator o = other.getCovVecIt();
        std::vector<int>::iterator t = covVec.begin();

        int identicalLength = identical[1]-identical[0] +1;

        for(int i = 0;i < covVec.size();i++){
            *t += *o;
            t++;
            o++;
        }

        std::vector<int> newCov;
        std::string newSeq = "";

        std::vector<int>::iterator bgn;
        std::vector<int>::iterator nd;
        std::vector<int>::iterator thisOv = this->getCovAt(identical[0]);
        std::vector<int>::iterator otherOv = other.getCovAt(identical[2]);


        int startDiff;
        int otherSubsDiff = 0;
        int thisSubsDiff = 0;


        if(identical[0] >= identical[2]){
            bgn = cov.begin();
            startDiff = identical[0];
            newSeq += seq.substr(0,startDiff);
            thisSubsDiff = startDiff;
        }
        else{
            bgn = other.getCovAt(0);
            startDiff = identical[2];
            newSeq += (other.getSeq()).substr(0,startDiff);
            otherSubsDiff = startDiff;
        }

        int end1 = (cov.size() - identical[1])-1;
        int end2 = (other.length() -identical[3])-1;
        int endDiff;

        newSeq += seq.substr(identical[0],identicalLength);

        if(end1 >= end2){
            nd = this->getCovAt(identical[1]);
            newSeq += seq.substr(identical[1],end1);
            endDiff = end1;
        }
        else{
            nd = other.getCovAt(identical[3]);
            newSeq += (other.getSeq()).substr(identical[3],end2);
            endDiff = end2;
        }

        int first = 0;
        int second = 0;
        int last = 0;

        while(first < startDiff && second < identicalLength && last < endDiff){
            if(first < startDiff){
                newCov.push_back(*bgn);
                bgn++;
                first++;
            }
            else{
                if(second < (identical[1]-identical[0])){
                    newCov.push_back((*thisOv + *otherOv));
                    thisOv++;
                    otherOv++;
                    second++;
                }
                else{
                    newCov.push_back(*nd);
                    nd++;
                    second++;
                }
            }
        }

        this->cov = newCov;
        this->seq = newSeq;

        for(int i = 0;i < subs.size();i++){
            subs[i]->update(thisSubsDiff,&seq,&cov);
        }
        std::vector<SubContig*>::iterator it = other.getSubIt();
        for(int i = 0;i < other.getSubCount();i++){
            ((*it)->update(otherSubsDiff,&seq,&cov));
            subs.push_back(*it);
            it++;
            subCount++;
        }
        other.setSave(false);
        return this;

    }


};


//[[Rcpp::export]]
List makeChimericContigs(std::vector<std::string> seqNames, std::vector<int> minOverlaps, std::vector<std::vector<std::vector<int> > > covs,
                         std::vector<std::vector<std::vector<int> > > covVecs, std::vector<std::vector<int> > starts,
                         std::vector<std::vector<int> > ends ,std::vector<std::vector<std::string> > seqs,std::vector<std::vector<int> > ovStarts1,
                         std::vector<std::vector<int> > ovEnds1, std::vector<std::vector<int> > ovStarts2, std::vector<std::vector<int> > ovEnds2,
                         std::vector<std::vector<int> > seqID1, std::vector<std::vector<int> > seqID2, double minShare){

    std::vector<std::string>::iterator SN = seqNames.begin();
    std::vector<int>::iterator minOv = minOverlaps.begin();

    std::vector<std::vector<std::vector<int> > >::iterator perSeqCovs = covs.begin();
    std::vector<std::vector<std::vector<int> > >::iterator perSeqcovVecs = covVecs.begin();
    std::vector<std::vector<int> >::iterator perSeqStarts = starts.begin();
    std::vector<std::vector<int> >::iterator perSeqEnds = ends.begin();
    std::vector<std::vector<std::string> >::iterator perSeqSeqs = seqs.begin();

    std::vector<std::vector<int> >::iterator cov;
    std::vector<std::vector<int> >::iterator covVec;
    std::vector<int>::iterator start;
    std::vector<int>::iterator end;
    std::vector<std::string>::iterator seq;

    //################## initialize Contigs ########################

    ContigContainer res;

    std::vector<std::vector<Contig*> > contigsPerSeq;
    std::vector<Contig*> contigs;
    Contig tmpCont;

    std::vector<std::vector<SubContig*> > subsPerSeq;
    std::vector<SubContig*> subs;
    SubContig* tmpSub;

    for(int i = 0; i < seqNames.size();i++){
        cov = (*perSeqCovs).begin();
        covVec = (*perSeqcovVecs).begin();
        start = (*perSeqStarts).begin();
        end = (*perSeqEnds).begin();
        seq = (*perSeqSeqs).begin();



        for(int j = 0;j < ((*perSeqEnds).size());j++){

            tmpCont = Contig(*start,*end,i,*minOv,*seq,*SN,*cov,*covVec,&res,minShare);
            contigs.push_back(&tmpCont);
            tmpSub = tmpCont.getSubAt(0);
            subs.push_back(tmpSub);

        }

        contigsPerSeq.push_back(contigs);
        contigs.clear();
        subsPerSeq.push_back(subs);
        subs.clear();

        perSeqCovs++;
        perSeqcovVecs++;
        perSeqEnds++;
        perSeqSeqs++;
        perSeqStarts++;
        SN++;
        minOv++;
    }



    std::vector<std::vector<int> >::iterator ovStart1It = ovStarts1.begin();
    std::vector<std::vector<int> >::iterator ovEnds1It = ovEnds1.begin();
    std::vector<std::vector<int> >::iterator ovStarts2It = ovStarts2.begin();
    std::vector<std::vector<int> >::iterator ovEnds2It = ovEnds2.begin();
    std::vector<std::vector<int> >::iterator seqID1It = seqID1.begin();
    std::vector<std::vector<int> >::iterator seqID2It = seqID2.begin();
    std::vector<int>::iterator ovstart1;
    std::vector<int>::iterator ovend1;
    std::vector<int>::iterator ovstart2;
    std::vector<int>::iterator ovend2;
    std::vector<int>::iterator seqid1;
    std::vector<int>::iterator seqid2;


    std::vector<std::vector<Contig*> >::iterator contList = contigsPerSeq.begin();
    std::vector<Contig*>::iterator conts;
    std::vector<std::vector<Contig*> >::iterator contListJumper = contigsPerSeq.begin();
    std::vector<Contig*>::iterator contsJumper;
    std::vector<std::vector<SubContig*> >::iterator subsList = subsPerSeq.begin();
    std::vector<SubContig*>::iterator sub;
    std::vector<std::vector<SubContig*> >::iterator subsListJumper = subsPerSeq.begin();
    std::vector<SubContig*>::iterator subJumper;


    int sameSeqStart;
    int sameSeqEnd;
    int otherStart;
    int otherEnd;
    int otherSeq;

    int j = 0;
    int n = 0;

    for(int i = 0; i < seqNames.size();i++){
        ovstart1 = (*ovStart1It).begin();
        ovend1 = (*ovEnds1It).begin();
        ovstart2 = (*ovStarts2It).begin();
        ovend2 = (*ovEnds2It).begin();
        seqid1 = (*seqID1It).begin();
        seqid2 = (*seqID2It).begin();

        conts = (*contList).begin();
        sub = (*subsList).begin();

        while(j < (*ovStart1It).size() && n < (*contList).size()){

            if((*seqid1) == i){
                sameSeqStart = *ovstart1;
                sameSeqEnd = *ovend1;
                otherStart = *ovstart2;
                otherEnd = *ovend2;
                otherSeq = *seqid2;
            }
            else{
                if((*seqid2) == i){
                    sameSeqStart = *ovstart2;
                    sameSeqEnd = *ovend2;
                    otherStart = *ovstart1;
                    otherEnd = *ovend1;
                    otherSeq = *seqid1;
                }
                else{
                    sameSeqStart = 0;
                    sameSeqEnd = 0;
                }
            }

            if((*sub)->hasOverlap(sameSeqStart,sameSeqEnd)){
                contsJumper = (*next(contListJumper,otherSeq)).begin();
                subJumper = (*next(subsListJumper,otherSeq)).begin();
                bool last = false;
                while((!(*subJumper)->hasOverlap(otherStart,otherEnd) && last) && subJumper !=  (*next(subsListJumper,otherSeq)).end()){

                    last = (*subJumper)->hasOverlap(otherStart,otherEnd);
                    if(last){
                        (*contsJumper) = (*conts)->fuse(*contsJumper,*sub,*subJumper,sameSeqStart,sameSeqEnd,otherStart,otherEnd);
                    }

                    subJumper++;
                    contsJumper++;
                }

                if((*sub)->end > sameSeqEnd){
                    j++;
                    ovstart1++;
                    ovstart2++;
                    ovend1++;
                    ovend2++;
                    seqid1++;
                    seqid2++;
                }
                else{
                    if((*sub)->end < sameSeqEnd){
                        n++;
                        conts++;
                        sub++;
                    }
                    else{
                        n++;
                        conts++;
                        sub++;

                        j++;
                        ovstart1++;
                        ovstart2++;
                        ovend1++;
                        ovend2++;
                        seqid1++;
                        seqid2++;
                    }
                }

            }
            else{
                if((*sub)->start >= sameSeqEnd){
                    j++;
                    ovstart1++;
                    ovstart2++;
                    ovend1++;
                    ovend2++;
                    seqid1++;
                    seqid2++;
                }
                else{
                    n++;
                    conts++;
                    sub++;
                }
            }


        }


        ovStart1It++;
        ovEnds1It++;
        ovStarts2It++;
        ovEnds2It++;
        seqID1It++;
        seqID2It++;

        contList++;
        subsList++;
    }

    return res.finalize();
}
