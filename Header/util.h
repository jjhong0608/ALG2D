//
// Created by 조준홍 on 2022/02/18.
//

#ifndef ALG2D_UTIL_H
#define ALG2D_UTIL_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <omp.h>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef UNITVALUE
#define UNITVALUE 1.0000000000000000E0
#endif
#ifndef HALFVALUE
#define HALFVALUE 5.0000000000000000E-1
#endif
#ifndef ZEROVALUE
#define ZEROVALUE 0.0000000000000000E0
#endif
#ifndef NEARZERO
#define NEARZERO 1.0000000000000000E-10
#endif

namespace ALG {
bool isclose(double x, double y, double eps = 1.0E-10);

bool iszero(double x, double eps = 1.0E-10);

bool ispositive(double x, double eps = 1.0E-10);

bool isnegative(double x, double eps = 1.0E-10);

bool ispositivezero(double x, double eps = 1.0E-10);

bool isnegativezero(double x, double eps = 1.0E-10);

void printError(const std::string &functionName);

void printError(const char *functionName, const char *fmt, ...);
}// namespace ALG

#endif// ALG2D_UTIL_H
