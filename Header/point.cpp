//
// Created by 조준홍 on 2022/02/18.
//

#include "point.h"

ALG::point::point() = default;

ALG::point::point(const std::vector<double> &x) : vector(x) {}

ALG::point::point(int idx) : idx(idx) {}

ALG::point::point(const std::vector<double> &x, int idx) : vector(x), idx(idx) {}

ALG::point::~point() = default;

int ALG::point::getIdx() const {
    return idx;
}

void ALG::point::setIdx(int i) {
    point::idx = i;
}
