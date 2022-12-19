//
// Created by 조준홍 on 2022/02/18.
//

#ifndef ALG2D_UTIL_H
#define ALG2D_UTIL_H


#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <array>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <limits>
#include <sstream>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <utility>
#include <deque>
#include <cstdarg>
#include <random>

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
}


#endif //ALG2D_UTIL_H
