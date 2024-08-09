//
// Created by 조준홍 on 2022/02/18.
//

#include "vector.h"

ALG::vector::vector() : std::vector<double>() {}

ALG::vector::vector(const std::vector<double> &x) : std::vector<double>(x) {}

ALG::vector::~vector() = default;

double ALG::vector::norm() {
  double sum{};
  for (const auto &i : *this) {
    sum += i * i;
  }
  return sqrt(sum);
}

double ALG::vector::dot(const ALG::vector &src) {
  if (size() != src.size()) {
    printError("ALG::vector::dot", "The sizes of the two vectors do not match.");
  }
  double sum{};
  for (int i = 0; i < size(); ++i) {
    sum += at(i) * src.at(i);
  }
  return sum;
}

ALG::vector ALG::vector::cross(const ALG::vector &src) {
  return ALG::vector(std::vector<double>{at(1) * src.at(2) - at(2) * src.at(1), at(2) * src.at(0) - at(0) * src.at(2), at(0) * src.at(1) - at(1) * src.at(0)});
}

ALG::vector ALG::vector::unitVector() {
  ALG::vector vec = vector{};
  double n = norm();
  for (const auto &i : *this) {
    vec.emplace_back(i / n);
  }
  return vec;
}

ALG::vector ALG::vector::operator+(const vector &src) {
  if (size() != src.size()) {
    printError("ALG::src::operator+", "size of this vector (which is %d) is not equal to vector src (which is %d)", size(), src.size());
  }
  ALG::vector vec = vector{};
  for (int i = 0; i < size(); ++i) {
    vec.emplace_back(at(i) + src.at(i));
  }
  return vec;
}

ALG::vector ALG::vector::operator-(const vector &src) {
  if (size() != src.size()) {
    printError("ALG::src::operator-", "size of this vector (which is %d) is not equal to vector src (which is %d)", size(), src.size());
  }
  ALG::vector vec = vector{};
  for (int i = 0; i < size(); ++i) {
    vec.emplace_back(at(i) - src.at(i));
  }
  return vec;
}

ALG::vector ALG::vector::operator*(double d) {
  ALG::vector vec = vector{};
  for (const auto &i : *this) {
    vec.emplace_back(i * d);
  }
  return vec;
}

ALG::vector ALG::vector::operator/(double d) {
  ALG::vector vec = vector{};
  for (const auto &i : *this) {
    vec.emplace_back(i / d);
  }
  return vec;
}

ALG::vector ALG::vector::operator+(const vector &src) const {
  if (size() != src.size()) {
    printError("ALG::src::operator+", "size of this vector (which is %d) is not equal to vector src (which is %d)", size(), src.size());
  }
  ALG::vector vec = vector{};
  for (int i = 0; i < size(); ++i) {
    vec.emplace_back(at(i) + src.at(i));
  }
  return vec;
}

ALG::vector ALG::vector::operator-(const vector &src) const {
  if (size() != src.size()) {
    printError("ALG::src::operator+", "size of this vector (which is %d) is not equal to vector src (which is %d)", size(), src.size());
  }
  ALG::vector vec = vector{};
  for (int i = 0; i < size(); ++i) {
    vec.emplace_back(at(i) - src.at(i));
  }
  return vec;
}

ALG::vector ALG::vector::operator*(double d) const {
  ALG::vector vec = vector{};
  for (const auto &i : *this) {
    vec.emplace_back(i * d);
  }
  return vec;
}

ALG::vector ALG::vector::operator/(double d) const {
  ALG::vector vec = vector{};
  for (const auto &i : *this) {
    vec.emplace_back(i / d);
  }
  return vec;
}

double ALG::vector::operator*(const vector &src) {
  if (size() != src.size()) {
    printError("ALG::src::operator+", "size of this vector (which is %d) is not equal to vector src (which is %d)", size(), src.size());
  }
  double sum{};
  for (int i = 0; i < size(); ++i) {
    sum += at(i) * src.at(i);
  }
  return sum;
}

ALG::vector &ALG::vector::operator+=(const vector &src) {
  if (size() != src.size()) {
    printError("ALG::src::operator+=", "size of this vector (which is %d) is not equal to vector src (which is %d)", size(), src.size());
  }
  for (int i = 0; i < size(); ++i) {
    at(i) += src.at(i);
  }
  return *this;
}

ALG::vector &ALG::vector::operator-=(const vector &src) {
  if (size() != src.size()) {
    printError("ALG::src::operator-=", "size of this vector (which is %d) is not equal to vector src (which is %d)", size(), src.size());
  }
  for (int i = 0; i < size(); ++i) {
    at(i) -= src.at(i);
  }
  return *this;
}

ALG::vector &ALG::vector::operator*=(double d) {
  for (auto &i : *this) {
    i *= d;
  }
  return *this;
}

bool ALG::vector::operator<(const ALG::vector &src) {
  if (size() != src.size()) {
    printError("ALG::src::operator<", "size of this vector (which is %d) is not equal to vector src (which is %d)", size(), src.size());
  }
  for (int i = 0; i < size(); ++i) {
    if (!isclose(at(i), src.at(i))) {
      return at(i) < src.at(i);
    }
  }
  return false;
}

bool ALG::vector::operator>(const ALG::vector &src) {
  if (size() != src.size()) {
    printError("ALG::src::operator>", "size of this vector (which is %d) is not equal to vector src (which is %d)", size(), src.size());
  }
  for (int i = 0; i < size(); ++i) {
    if (!isclose(at(i), src.at(i))) {
      return at(i) > src.at(i);
    }
  }
  return false;
}
