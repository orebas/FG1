#include <iostream>
#include <vector>
#include <iterator>
#include <array>
#include <cmath>

template <class T, class V>
std::ostream &operator<<(std::ostream &out, const std::pair<T,V> &v) {
  out <<"[" << v.first <<"," << v.second << "]";
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
    std::vector<std::pair<double,double>> cdvec;
    std::vector<double> dvec;
    std::pair<double, double> dpair(2.1,-3.2);
     for(int i=0; i<25;i++){
        cdvec.emplace_back(std::pair<double,double>(i,1-i));
        dvec.push_back(std::pow(0.1,i));
    }
    std::cout <<dvec << std::endl;
    std::cout << dpair << std::endl;
    std::cout << cdvec << std::endl;
       
}