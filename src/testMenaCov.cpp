#include "meanCovToRange.cpp"

int main(){
    std::vector<int> ranges;
    int temp;
    int n = 0;
    while(n < 4){
        std::cin >> temp;
        ranges.push_back(temp);
        n++;
    }
    n = 0;
    std::cout << "Nun der coverage vector"<< std::endl;
    std::string s;
    std::cin >> s;
    std::vector<int> coverage;
    temp = 0;
    while(n < 10){
        std::cin >> temp;
        coverage.push_back(temp);
        n++;
    }
    n = 0;
    std::list<std::vector<int> > res = meanCovToRange(ranges,coverage);
    for(std::list<std::vector<int> >::iterator i = res.begin();i != res.end();i++){
        std::cout<< n << std::endl;
        for(int j = 0;j < (int) (*i).size();j++){
            std::cout << (*i)[j];
        }
        std::cout << std::endl;
        n++;
    }
}
