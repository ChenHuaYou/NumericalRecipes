#include <fstream>
#include "nr.h"
#include <string>
#include <iostream>
#include <iomanip>

int main(){

    int j,k,l,m=2,n=3;
    char ts = '\n';
    float tmp;
    std::ifstream fp("data/matrx1.dat");
    string dummy;
    if (fp.fail()){
      NR::nrerror("Data file matrx1.dat not found");
    }
    cout << fixed << setprecision(6);
    while (!fp.eof()) {
          getline(fp,dummy);
          getline(fp,dummy);
          fp >> n >> m;
          fp.ignore(2);
          getline(fp,dummy);
          for (k=0;k<n;k++){
            for (l=0;l<n;l++){ 
                fp >> tmp;
                cout << tmp << " ";
            }
            cout << endl;
          }
          fp.ignore(2);
          getline(fp,dummy);
          for (k=0;k<m;k++){
            for (l=0;l<n;l++){ 
                fp >> tmp;
                cout << tmp << " ";
            }
            cout << endl;
          }
          fp.ignore(2);
    }
}
