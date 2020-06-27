

class Test{
public:
    Test(int x){
        this->x = x;
        y = 3;
    }

    int getY(){
        return y;
    }

    int x;
    int y;
};

//[[Rcpp::export]]
int testTest(int x){

    Test a = Test(x);
    return a.getY();

}
