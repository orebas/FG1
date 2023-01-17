#include "mpreal.h"
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
