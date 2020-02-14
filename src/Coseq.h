#ifndef COSEQ
#define COSEQ

#include <Rcpp.h>
#include <vector>
#include <algorithm>

class Coseq
{
public:
    Coseq(std::vector<int>& cont1Starts,std::vector<std::vector<int> >& cont1Cov,std::vector<int>& same1Starts,std::vector<int>& same1Ends,std::vector<int>& cont2Starts,std::vector<std::vector<int> >& cont2Cov,std::vector<int>& same2Starts,std::vector<int>& same2Ends,int length1,int length2)
    {
        this->length1 = length1;
        this->length2 = length2;
        seq1 = (int*) malloc(sizeof(int) *length1);
        seq2 = (int*) malloc(sizeof(int) *length2);
        std::fill_n(seq1,length1,0);
        std::fill_n(seq2,length2,0);
        this->cont1Starts = (int**) malloc(sizeof(int**) *cont1Starts.size());
        this->cont1Ends = (int**) malloc(sizeof(int**) *cont1Starts.size());
        this->cont2Starts = (int**) malloc(sizeof(int**) *cont2Starts.size());
        this->cont2Ends = (int**) malloc(sizeof(int**) *cont2Starts.size());
        fillSeqs(cont1Starts,cont1Cov,cont2Starts,cont2Cov);
        this->same1Starts = (int**) malloc(sizeof(int**) *same1Starts.size());
        this->same1Ends = (int**) malloc(sizeof(int**) *same1Starts.size());
        this->same2Starts = (int**) malloc(sizeof(int**) *same2Starts.size());
        this->same2Ends = (int**) malloc(sizeof(int**) *same2Starts.size());
        setSames(same1Starts.begin(),same1Ends.begin(),same2Starts.begin(),same2Ends.begin(),(int)same1Starts.size());
        nrSames = (int) same1Starts.size();
    }

    ~Coseq()
    {
        free(seq1);
        free(seq2);
        free(cont1Starts);
        free(cont1Ends);
        free(cont2Starts);
        free(cont2Ends);
        free(same1Starts);
        free(same1Ends);
        free(same2Starts);
        free(same2Ends);
    }

    void fillSeqs(std::vector<int>& cont1Starts,std::vector<std::vector<int> >& cont1Cov,std::vector<int>& cont2Starts,std::vector<std::vector<int> >& cont2Cov)
    {
        std::vector<int>::iterator c1s = cont1Starts.begin();
        std::vector<std::vector<int> >::iterator c1c = cont1Cov.begin();
        for(int i = 0;i < (int)cont1Starts.size();i++)
        {
            std::copy((*c1c).begin(),(*c1c).end(),seq1+(*c1s));
            *(this->cont1Starts+i) = seq1+(*c1s);
            *(this->cont1Ends +i) = seq1+(*c1s)+(*c1c).size();
            c1s++;
            c1c++;
        }

        std::vector<int>::iterator c2s = cont2Starts.begin();
        std::vector<std::vector<int> >::iterator c2c = cont2Cov.begin();
        for(int i = 0;i < (int)cont1Starts.size();i++)
        {
            copy((*c2c).begin(),(*c2c).end(),seq2+(*c2s));
            *(this->cont2Starts+i) = seq1+(*c2s);
            *(this->cont2Ends +i) = seq1+(*c2s)+(*c2c).size();
            c1s++;
            c1c++;
        }

    }


    void setSames(std::vector<int>::iterator s1s,std::vector<int>::iterator s1e,std::vector<int>::iterator s2s,std::vector<int>::iterator s2e,int l)
    {
        for(int i = 0; i < l; i++){
            *(same1Starts+i) = seq1+(*s1s);
            *(same1Ends+i) = seq1+(*s1e);
            *(same2Ends+i) = seq2+(*s2s);
            *(same2Ends+i) = seq2+(*s2e);
            s1s++;
            s1e++;
            s2s++;
            s2e++;
        }
    }

    Rcpp::List mkChimeraConts()
    {
        for(int i = 0;i < nrSames;i++)
        {
            for(int j = 0;j < *(same1Ends+i)- *(same1Starts+i)+1;j++)
            {

            }
        }
    }

private:
    int nrSames;

    int length1;
    int* seq1;
    int** cont1Starts;
    int** cont1Ends;
    int** same1Starts;
    int** same1Ends;

    int length2;
    int* seq2;
    int** cont2Starts;
    int** cont2Ends;
    int** same2Starts;
    int** same2Ends;
};


#endif
