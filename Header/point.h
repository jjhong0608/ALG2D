//
// Created by 조준홍 on 2022/02/18.
//

#ifndef ALG2D_POINT_H
#define ALG2D_POINT_H

#include "vector.h"

namespace ALG {
class point : public ALG::vector {
 private:
  int idx{-1};

 public:
  point();

  explicit point(int idx);

  explicit point(const std::vector<double> &x);

  point(const std::vector<double> &x, int idx);

  ~point() override;

  [[nodiscard]] int getIdx() const;

  void setIdx(int i);
};
}// namespace ALG

#endif// ALG2D_POINT_H
