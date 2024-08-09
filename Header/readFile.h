//
// Created by 조준홍 on 2022/02/19.
//

#ifndef ALG2D_READFILE_H
#define ALG2D_READFILE_H

#include "region.h"

namespace ALG {
class readFile {
 public:
  static void readGeometryInfomation(const std::string &filename, std::vector<section> *sections, std::vector<region> *regions);

  static void readSection(std::ifstream &f, std::vector<section> *sections);

  static void readRegion(std::ifstream &f, std::vector<region> *regions, std::vector<section> *sections);
};

}// namespace ALG

#endif// ALG2D_READFILE_H
