#include <Rcpp.h>
#include <vector>
#include <list>
#include <cstdlib>


using namespace Rcpp;

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
    int a = (as1-c1s > 0)+(as2-c2s > 0);
    int b = (c1e-ae1 > 0) +(c2e-ae2 > 0);
    std::vector<int> res = {0+noNeg((as1-c1s)+(c2s-as2)),(ae1-c1s)-noNeg((c1e-ae1)+(ae2-c2e)),0+noNeg((as2-c2s)+(c1s-as1)),(ae2-c2s)-noNeg((c2e-ae2)+(ae1-c1e))};
    Rcout << c1e-c1s << " " << c1e-ae1 << " " << ae2-c2e << std::endl;
    if(a == 0 || b == 0)
    {
        res.push_back(0);
        if(a == 2 || b == 2)
        {
            res.push_back(2);
        }
    }
    else
    {
        if(a == 1 || b == 1)
        {
            if(a == 1 && b == 1)
            {
                res.push_back(1);
            }
            else
            {
                res.push_back(2);
            }
        }
        else
        {
            res.push_back(3);
        }
    }
    return res;
}



