
#include <Rcpp.h>
#include <vector>

class Test{
public:
    ~Test(){

    }

    Test(int x){
        this->x = x;
        y = 3;
    }

    int getY(){
        return y;
    }

    int getX(){
        return x;
    }

    int x;
    int y;
};

std::ostream& operator<<(std::ostream &strm, const Test &cont){
    return strm << "TestVal: " << cont.x;
}

//[[Rcpp::export]]
int testTest(int x){

    Test* m;
    std::vector<Test*> z;
    for(int i = 0;i < 10;i++){
        m = new Test(i);
        z.push_back(m);
    }

    for(int i = 0;i < 10;i++){
        Rcpp::Rcout << *(z[i]) << std::endl;
        delete (z[i]);
    }

    return 1;

}
