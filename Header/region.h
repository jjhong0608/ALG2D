//
// Created by 조준홍 on 2022/02/18.
//

#ifndef ALG2D_REGION_H
#define ALG2D_REGION_H

#include "section.h"

namespace ALG {
class region : public std::vector<section *> {
 private:
  std::string name{};
  double epsilon{}, xmin{}, xmax{}, ymin{}, ymax{}, hx{}, hy{};
  std::vector<double> sign{};
  int ptsnum{};

 public:
  explicit region(std::string name);

  region(double epsilon, double xmin, double xmax, double ymin, double ymax, double hx, double hy);

  region(std::string name, double epsilon, double xmin, double xmax, double ymin, double ymax, double hx, double hy);

  [[nodiscard]] int getPtsnum() const;

  void setPtsnum(int i);

  [[nodiscard]] const std::string &getName() const;

  void setName(const std::string &name);

  [[nodiscard]] double getEpsilon() const;

  void setEpsilon(double epsilon);

  [[nodiscard]] double getXmin() const;

  void setXmin(double xmin);

  [[nodiscard]] double getXmax() const;

  void setXmax(double xmax);

  [[nodiscard]] double getYmin() const;

  void setYmin(double ymin);

  [[nodiscard]] double getYmax() const;

  void setYmax(double ymax);

  [[nodiscard]] double getHx() const;

  void setHx(double hx);

  [[nodiscard]] double getHy() const;

  void setHy(double hy);

  [[nodiscard]] const vector<double> &getSign() const;

  void setSign(const vector<double> &s);

  void writeRegion(std::ofstream &f, std::deque<boundaryPoint> *boundary_pts);
};

}// namespace ALG

#endif// ALG2D_REGION_H
