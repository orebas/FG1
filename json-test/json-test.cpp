#include "../mpreal.h"
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

template <class T>
std::ostream &operator<<(std::ostream &out, const std::complex<T> &v) {
  out << "[" << v.real() << ", " << v.imag() << "]" << std::endl;
  return out;
}

template <class T, class V>
std::ostream &operator<<(std::ostream &out, const std::pair<T, V> &v) {
  out << "[" << v.first << ", " << v.second << "]" << std::endl;
  return out;
}

template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  if (!v.empty()) {
    out << '[';
    for (int i = 0; i < v.size() - 1; i++) {
      out << v[i] << ", ";
    }
    out << v[v.size() - 1];
    out << "]";
  }
  out << std::endl;
  return out;
}

#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace mpfr {
void to_json(json &j, const mpreal &p) { j = json(p.toString()); }

void from_json(const json &j, mpreal &p) {
  std::string d;
  j.get_to(d);
  p = mpfr::mpreal(d);
}
}; // namespace mpfr

int main(int argc, char **argv) {

  mpfr::mpreal::set_default_prec(mpfr::digits2bits(300));

  std::pair<double, double> dpair1(2.1, -3.2), dpair2, dpair3;
  std::vector<double> dvec, dvec2, dvec3;
  std::vector<std::pair<double, double>> cdvec, cdvec2, cdvec3;
  std::vector<mpfr::mpreal> mp1(1.6), mp2, mp3;
  // std::vector<std::complex<double>> ccdvec, ccdvec2, ccdvec3;
  //  std::vector<std::complex<mpfr::mpreal>> cmpvec;
  mpfr::mpreal largesmall = 10.0;

  for (int i = 0; i < 12; i++) {
    cdvec.emplace_back(std::pair<double, double>(i, 1 - i));
    // ccdvec.emplace_back(std::complex<double>(i, 2 - i));
    dvec.push_back(std::pow(0.1, i));
    mp1.push_back(mpfr::pow(largesmall, i * 2) + 1 +
                  mpfr::pow(largesmall, -(i * 2)));
  }
  std::cout << mp1 << std::endl;
  json j1(dpair1);
  json j2(dvec);
  json j3(cdvec);
  json j5(mp1);
  std::cout << j5 << std::endl;
  // json j4(ccdvec);

  dpair2 = j1.get<std::pair<double, double>>();
  dvec2 = j2.get<std::vector<double>>();
  cdvec2 = j3.get<std::vector<std::pair<double, double>>>();
  mp2 = j5.get<std::vector<mpfr::mpreal>>();
  // ccdvec2 = j4.get<std::vector<std::complex<double>>>();

  std::cout << "JSON test\n";
  std::cout << dpair1 << " " << dpair2 << std::endl;
  std::ofstream odpair("dpair.json");
  odpair << j1 << std::endl;
  odpair.close();
  std::cout << dvec << std::endl;
  std::cout << dvec2 << std::endl;
  std::ofstream odvec("dvec.json");
  odvec << j2 << std::endl;

  std::cout << cdvec << std::endl;
  std::cout << cdvec2 << std::endl;
  std::ofstream ocdvec("cdvec.json");
  ocdvec << j3 << std::endl;

  std::cout << mp1 << std::endl;
  std::cout << mp2 << std::endl;
  std::ofstream omp("mp.json");
  omp << j5 << std::endl;

  // std::cout << ccdvec << std::endl;
  // std::cout << ccdvec2 << std::endl;
  // std::ofstream occdvec("ccdvec.json");
  // ocdvec << j4 << std::endl;

  odpair.close();
  odvec.close();
  ocdvec.close();
  omp.close();
  // occdvec.close();

  std::cout << "File test" << std::endl;
  std::ifstream idpair("dpair.json");
  std::ifstream idvec("dvec.json");
  std::ifstream icdvec("cdvec.json");
  std::ifstream imp("mp.json");

  // std::ifstream iccdvec("cdvec.json");

  json k1, k2, k3, k5;
  //  json k4;

  idpair >> k1;
  idvec >> k2;
  icdvec >> k3;
  imp >> k5;

  //  icdvec >> k4;
  dpair3 = k1.get<std::pair<double, double>>();
  dvec3 = k2.get<std::vector<double>>();
  cdvec3 = k3.get<std::vector<std::pair<double, double>>>();
  mp3 = k5.get<std::vector<mpfr::mpreal>>();

  // ccdvec3 = k4.get<std::vector<std::complex<double>>>();

  std::cout << dpair1 << " " << dpair3 << std::endl;
  std::cout << dvec << std::endl;
  std::cout << dvec3 << std::endl;

  std::cout << cdvec << std::endl;
  std::cout << cdvec3 << std::endl;

  std::cout << mp1 << std::endl;
  std::cout << mp3 << std::endl;

  //  std::cout << ccdvec << std::endl;
  //  std::cout << ccdvec3 << std::endl;
}
