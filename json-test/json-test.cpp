#include <iostream>
#include <vector>
#include <iterator>
#include <array>
#include <cmath>
#include <complex>
#include "../mpreal.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

template <class T, class V>
std::ostream &operator<<(std::ostream &out, const std::pair<T,V> &v) {
  out <<"[" << v.first <<", " << v.second << "]" << std::endl;
  return out;
}


template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  if (!v.empty()) {
    out << '[';
    for(int i=0; i < v.size()-1;i++){
      out << v[i] <<", ";
    }
    out <<v[v.size()-1];
    out << "]";
  }
  out << std::endl;
  return out;
}



int main(int argc, char** argv){
    std::pair<double, double> dpair1(2.1,-3.2), dpair2;
    std::vector<double> dvec, dvec2, dvec3;
    std::vector<std::pair<double,double>> cdvec, cdvec2, cdvec3;
    //std::vector<std::complex<mpfr::mpreal>> cmpvec;
     
     for(int i=0; i<25;i++){
        cdvec.emplace_back(std::pair<double,double>(i,1-i));
        dvec.push_back(std::pow(0.1,i));
    }
    
    json j1(dpair1);
    json j2(dvec);
    json j3(cdvec);
    dpair2 = j1.get<std::pair<double, double>>();
    dvec2 = j2.get<std::vector<double>>();
    cdvec2=j3.get<std::vector<std::pair<double,double>>>();
    std::cout << "JSON test\n";
    std::cout << dpair1 << " " << dpair2 << std::endl;
 
    std::cout <<dvec << std::endl;
    std::cout << dvec2 << std::endl;
    

    std::cout << cdvec << std::endl;
    std::cout << cdvec2 << std::endl;
    std::cout << "File test" << std::endl;

       
}