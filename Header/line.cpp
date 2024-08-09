//
// Created by 조준홍 on 2022/02/18.
//

#include "line.h"

ALG::line2D::line2D() = default;

ALG::line2D::line2D(const ALG::vector &start, const ALG::vector &anEnd) : start(start), end(anEnd) {
  calcProperties();
}

ALG::line2D::line2D(double s0, double s1, double e0, double e1) {
  start.emplace_back(s0);
  start.emplace_back(s1);
  end.emplace_back(e0);
  end.emplace_back(e1);
  calcProperties();
}

void ALG::line2D::calcProperties() {
  auto line = end - start;
  length = line.norm();
  normal = vector(std::vector<double>{-line[1] / length, line[0] / length});
  tangent = line.unitVector();
}

ALG::vector &ALG::line2D::getStart() {
  return start;
}

const ALG::vector &ALG::line2D::getStart() const {
  return start;
}

void ALG::line2D::setStart(const ALG::vector &vector) {
  if (vector.size() != 2) {
    printError("ALG::line2D::setStart", "start vector size (which is %d) is not 2", vector.size());
  }
  start = vector;
  if (start.size() == end.size()) {
    calcProperties();
  }
}

ALG::vector &ALG::line2D::getAnEnd() {
  return end;
}

const ALG::vector &ALG::line2D::getAnEnd() const {
  return end;
}

void ALG::line2D::setAnEnd(const ALG::vector &vector) {
  if (vector.size() != 2) {
    printError("ALG::line2D::setAnEnd", "end vector size (which is %d) is not 2", vector.size());
  }
  end = vector;
  if (start.size() == end.size()) {
    calcProperties();
  }
}

ALG::vector &ALG::line2D::getNormal() {
  return normal;
}

const ALG::vector &ALG::line2D::getNormal() const {
  return normal;
}

ALG::vector &ALG::line2D::getTangent() {
  return tangent;
}

const ALG::vector &ALG::line2D::getTangent() const {
  return tangent;
}

double ALG::line2D::getLength() const {
  return length;
}

bool ALG::line2D::getIscloseLine() const {
  return iscloseline;
}

bool ALG::line2D::iscross(const ALG::line2D &src, ALG::vector &vec) {
  if (isclose(std::fabs((src.end - src.start).unitVector() * tangent), UNITVALUE)) {
    return false;
  }
  auto projection = [this]() -> matrix {
    matrix mat = matrix(2, 2);
    mat[0][0] = UNITVALUE - tangent[0] * tangent[0];
    mat[0][1] = -tangent[0] * tangent[1];
    mat[0][1] = -tangent[0] * tangent[1];
    mat[1][1] = UNITVALUE - tangent[1] * tangent[1];
    return mat;
  };
  auto P{projection()};
  auto projectedStart{P * src.start};
  auto projectedEnd{P * src.end};
  auto projectedPoint{P * start};

  if ((projectedPoint - projectedStart) * (projectedPoint - projectedEnd) < NEARZERO) {
    double lenStart{(projectedPoint - projectedStart).norm()};
    double lenEnd{(projectedPoint - projectedEnd).norm()};
    double ratioStart{lenStart / (lenStart + lenEnd)};
    double t{((src.start + (src.end - src.start) * ratioStart - start) * (end - start)) / std::pow((end - start).norm(), 2)};
    if (ispositive(t * (t - UNITVALUE)) || std::isnan(ratioStart)) {
      return false;
    } else {
      if (iszero(t * (t - UNITVALUE))) {
        iscloseline = true;
      }
      double len{src.getLength()};
      vec = src.start + ((src.end - src.start) * (ratioStart * len)) / len;
      return true;
    }
  } else {
    return false;
  }
}

ALG::line2D::~line2D() = default;
