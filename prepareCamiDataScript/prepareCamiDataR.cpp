#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <Rcpp.h>

int* genomeInidicies;
std::vector<std::string> uniqeG;
int last8Val;
int first7Val = 0;
std::string first7;
int nrOfLines;


int calcFirst7Val(std::string first7){
  int x = 1111111;
  int y = 1111111;
  int res = 0;

  for(int i = 0;i < first7.size();i++){
    res += (first7[i] -'0') *x +1;
    x = x/10;
  }
  if(res == 1){
    res = 0;
  }
  return res-y;
}


//[[Rcpp::export]]
int readHash(std::string readName,std::string last8){
    bool is8 = true;
    int i = 0;

    while(i < readName.size() && i < last8.size()){
        if(readName[i] > last8[i]){
          // Rcpp::Rcout << "1\n";
            is8 = false;
            break;
        }
        if(readName[i] < last8[i]){
          // Rcpp::Rcout << "2\n";
            break;
        }
        i++;
    }
    if(is8 && readName.size() > last8.size() && i == last8.size()){
      // Rcpp::Rcout << "3\n";
      is8 = false;
    }

    int x;
    int y;
    if(is8){
        x = 11111111;
        y = 11111111;
    }
    else{
        x = 1111111;
        y = 1111111;
    }
    int res = 0;
    for(i = 0;i < readName.size();i++){
        res += (readName[i] -'0') *x +1;
        x = x/10;
    }
    // Rcpp::Rcout << "hash 1: " << res << std::endl;
    if(res == 1){
        res = 0;
    }



    if(res != 0){
        res -= y;

        // Rcpp::Rcout << "hash 3: " << res << std::endl;
        if(!is8){
            res -= (first7Val-1);
          // Rcpp::Rcout << "hash 4: " << res  << " " <<  first7Val<< std::endl;
            res += last8Val;
            // Rcpp::Rcout << "hash 5: " << res << " " << last8Val  << std::endl;
        }

    }
    return res;
}

//[[Rcpp::export]]
void testReadHash(std::string l8File,std::string linkingFile){
  std::ifstream file;
  file.open(l8File);
  std::string last8;
  if(file.is_open()){
    std::getline(file,last8);
    std::getline(file,first7);
  }
  file.close();

  last8Val = readHash(last8,"99999999");
  first7Val = first7Val = calcFirst7Val(first7);

  std::ifstream link(linkingFile);
  std::string line1;
  std::string line2;
  int pos;
  int pos1;
  int one;
  int two;
  int i = 1;
  if(link.is_open()){
    if(std::getline(link,line1)){
      while(line1.size() == 0 || line1[0] == '@'){
        std::getline(link,line1);
      }

      pos = line1.find("\t");
      line1 = line1.substr(0,pos);
      pos1 = line1.find_last_of("|");
      line1 = line1.substr(pos1+2,pos-(pos1+2));

      while(std::getline(link,line2)){

        pos = line2.find("\t");
        line2 = line2.substr(0,pos);
        pos1 = line2.find_last_of("|");
        line2 = line2.substr(pos1+2,pos-(pos1+2));
        one = readHash(line1,last8);
        two = readHash(line2,last8);

        if(two -one != 1){
          Rcpp::Rcout << "line 1: " << readHash(line1,last8) << " " <<line1  <<std::endl;
          Rcpp::Rcout << "line 2: " << readHash(line2,last8) << " " << line2 << std::endl;
          Rcpp::Rcout << first7Val << " " << last8Val << " " << last8 << std::endl;
          break;
        }
        i++;
        line1 = line2;

      }
    }
  }
  Rcpp::Rcout << i << " " << one << " " << two << std::endl;
  link.close();


}

//[[Rcpp::export]]
void makeCAMIGenomeReadLink(std::string linkingFile,std::string last8){

    Rcpp::Rcout << "start linking\n";
    std::ifstream lFile;
    std::ifstream numberFi(linkingFile+".numOfLines");

    if(numberFi.is_open()){
      std::string nr;
      std::getline(numberFi,nr);
      nrOfLines = std::stoi(nr);
    }
    else{
      nrOfLines = 0;
    }

    std::string tmpS;
    if(nrOfLines == 0){
      std::ofstream numberFile(linkingFile+".numOfLines");
      lFile.open(linkingFile);
      while(std::getline(lFile,tmpS)){
          if(tmpS.size() > 0 && tmpS.at(0) != '@'){
              nrOfLines++;
          }
      }
      lFile.close();
      numberFile << nrOfLines << std::endl;
      numberFile.close();
    }

    genomeInidicies = new int[nrOfLines];
    int i = 0;


    std::vector<std::string>::iterator genIt;
    lFile.clear();
    lFile.open(linkingFile);
    // Rcpp::Rcout << "start linking 1\n";
    if(lFile.is_open()){


        std::string readName;
        std::string GenomeName;
        while(std::getline(lFile,tmpS)){
            if(tmpS.size() > 0 && tmpS.at(0) != '@' ){
                std::istringstream s(tmpS);
                s >> readName;
                s >> GenomeName;
                // Rcpp::Rcout << "start linking 2\n";
                int pos = readName.find_last_of("|");
                readName = readName.substr(pos+2,readName.size()-(pos+2));
                if(uniqeG.size() == 0 || (genIt = std::find(uniqeG.begin(),uniqeG.end(),GenomeName)) == uniqeG.end()){
                    uniqeG.push_back(GenomeName);
                  // Rcpp::Rcout << "start linking 3\n";
                    *(genomeInidicies+readHash(readName,last8)) = i;
                    // Rcpp::Rcout << "start linking 4\n";
                    i++;

                }
                else{
                  // Rcpp::Rcout << readHash(readName,last8) << " " << readName << " " << nrOfLines << "start linking 5\n";
                    *(genomeInidicies+readHash(readName,last8)) = (distance(uniqeG.begin(),genIt));
                    // Rcpp::Rcout << "start linking 6\n";
                }
            }
        }

    }
    else{

      Rcpp::Rcout << "couldnt open file" << std::endl;
    }
    Rcpp::Rcout << "Nr of Genomes: " << i << std::endl;

    lFile.close();
    Rcpp::Rcout << "end linking\n";
}

//[[Rcpp::export]]
void writeToFq(std::string fastq,std::string saveDir,std::string last8){

    std::ofstream log("/home/yang/CamiPrepLog.txt");

    Rcpp::Rcout << "start writing\n";
    if(saveDir.at(saveDir.size()-1) != '/'){
        saveDir.push_back('/');
    }

    std::ifstream fq;
    fq.open(fastq);
    int i = 1;
    int n;
    int pos;
    int fwdOrRev;
    int hash;
    std::ofstream outs[2][uniqeG.size()];
    std::vector<std::string>::iterator unqIt;
    std::string tmp;
    std::string tmp2;
    unqIt = uniqeG.begin();


    for(int p = 0;p < uniqeG.size();p++){
        (*(*(outs)+p)).open(saveDir + (*unqIt) + "_1_.fq",std::ios_base::app);
        (*(*(outs+1)+p)).open(saveDir + (*unqIt) + "_2_.fq",std::ios_base::app);
        unqIt++;
    }
    while(std::getline(fq,tmp)){
      if(fq.good()){
          if(tmp.size() > 0 && tmp.at(0) == '@'){
              fwdOrRev = tmp[tmp.size()-1] -'1';
              pos = tmp.find_last_of("|");
              tmp2 = tmp.substr(pos+2, (tmp.size()-(pos+2)-2) );

              hash = readHash(tmp2,last8);

              n = *(genomeInidicies+hash);

              if(n < 0 || n >= uniqeG.size()){
                log << n << " wrong nr for: " << tmp << std::endl;
              }


              if((*(*(outs+fwdOrRev)+n)).good()){
                (*(*(outs+fwdOrRev)+n)) << tmp << std::endl;
                std::getline(fq,tmp);
                (*(*(outs+fwdOrRev)+n)) << tmp << std::endl;
                std::getline(fq,tmp);
                (*(*(outs+fwdOrRev)+n)) << tmp << std::endl;
                std::getline(fq,tmp);
                (*(*(outs+fwdOrRev)+n)) << tmp << std::endl;

                if(i % 4000000 == 0)
                  Rcpp::Rcout << "writing: " << i << std::endl;

                i += 4;
              }
              else{
                log << "out not good: " << n << " " << i << std::endl;
              }
          }
          else{
            i++;
            log << "is length 0 or not @ :" << tmp << std::endl;
          }
        }
        else{
           log << "Problem found!\n";
        }


    }

    for(int l = l;l < uniqeG.size();l++){
        (*(*(outs)+l)).close();
        (*(*(outs+1)+l)).close();
    }
    fq.close();


    log.close();
    Rcpp::Rcout << i << " end writing\n";
}
//[[Rcpp::export]]
std::string getLast8(std::string fileName){

    std::ifstream file;
    file.open(fileName);
    std::string line1;
    std::string line2;
    int pos;
    int pos1;
    int pos2;

    int i = 3;
    std::string res;
    Rcpp::Rcout << file.is_open() << std::endl;
    if(file.is_open()){
        while(std::getline(file,line1)){
            if(line1.size() > 0){
                if(line1.at(0) != '@'){
                    break;
                }
            }
        }
        pos = line1.find("\t");
        line1 = line1.substr(0,pos);
        pos1 = line1.find_last_of("|");
        line1 = line1.substr(pos1+2,pos-(pos1+2));
            while(std::getline(file,line2)){

                pos = line2.find("\t");
                line2 = line2.substr(0,pos);
                pos1 = line2.find_last_of("|");
                line2 = line2.substr(pos1+2,pos-(pos1+2));

                int one = readHash(line1,"99999999");
                int two = readHash(line2,"99999999");
                if(two - one != 1){

                    res = line1;
                    first7 = line2;
                    Rcpp::Rcout << res << " " << one << " " << two << std::endl << line1 << std::endl << line2 << std::endl;
                    break;
                }

                line1 = line2;
                if(i % 1000000 == 0){
                    Rcpp::Rcout << "get last 8: " << i << std::endl;
                }
                i++;
            }

    }
    file.close();

    std::ofstream l8;
    l8.open((fileName+".l8"));
    if(l8.is_open()){
        l8 << res << std::endl;
        l8 << first7 << std::endl;
    }
    l8.close();
    return res;

}

//[[Rcpp::export]]
int prepare(std::vector<std::string> inFiles){

    std::string readFile = inFiles[0];
    std::string linkingFile = inFiles[1];
    std::string saveToDir = inFiles[2];
    std::string last8;
    std::ifstream test;
    test.open(readFile);

    std::string l8File = linkingFile + ".l8";

    std::ifstream l8;
    l8.open(l8File);
    if(l8.is_open()){
        std::getline(l8,last8);
        std::getline(l8,first7);
    }
    else{
        last8 = getLast8(linkingFile);
    }
    l8.close();

    last8Val = readHash(last8,"99999999");
    first7Val = calcFirst7Val(first7);

    makeCAMIGenomeReadLink(linkingFile,last8);
    writeToFq(readFile,saveToDir,last8);
    delete genomeInidicies;
    return 1;
}

// taken from https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
//[[Rcpp::export]]
std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

//[[Rcpp::export]]
void prepareCamiData(std::string inFileName){
    std::ifstream in_csv;
    in_csv.open(inFileName);
    std::string tmp;
    if(in_csv.is_open()){

      while(std::getline(in_csv,tmp)){

        prepare(split(tmp,','));
        uniqeG.clear();
      }

    }
    in_csv.close();
}

// int main(int argc,char** argv){
//     std::string inFileName = argv[1];
//     std::ifstream in_csv;
//     in_csv.open(inFileName);
//     std::string tmp;
//     if(in_csv.is_open()){
//
//         while(std::getline(in_csv,tmp)){
//
//             prepare(split(tmp,','));
//
//         }
//
//     }
//     in_csv.close();
// }

//[[Rcpp::export]]
void testFileLimit(int lim){
  std::ofstream files[lim];
  for(int i = 0; i < lim;i++){
    files[i].open(("/home/yang/testFileLim/"+ std::to_string(i) + ".txt"));
    files[i] << i << std::endl;
  }
  for(int i = 0; i < lim;i++){
    files[i].close();
  }
}




void testPrepare(std::vector<std::string> multiFqPath,std::string linkingFile,std::string l8File){

    std::ifstream l8;
    std::string last8;
    l8.open(l8File);
    if(l8.is_open()){
      std::getline(l8,last8);
      std::getline(l8,first7);
    }
    else{
      last8 = getLast8(linkingFile);
    }
    l8.close();

    last8Val = readHash(last8,"99999999");
    readHash(first7,last8);

    makeCAMIGenomeReadLink(linkingFile,last8);



}



