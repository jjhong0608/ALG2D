//
// Created by 조준홍 on 2022/02/18.
//

#include "util.h"

bool ALG::isclose(double x, double y, double eps) {
  return std::fabs(x - y) < eps;
}

bool ALG::iszero(double x, double eps) { return std::fabs(x) < eps; }

bool ALG::ispositive(double x, double eps) { return x - eps > ZEROVALUE; }

bool ALG::isnegative(double x, double eps) { return x + eps < ZEROVALUE; }

void ALG::printError(const std::string &functionName) {
  std::cout << '\n'
            << functionName << '\n';
  exit(1);
}

void ALG::printError(const char *functionName, const char *fmt, ...) {
  char buf[256] = {
      0,
  };
  va_list ap;

  printf("\n");
  printf("Fatal error has occur in %s\n", functionName);
  sprintf(buf, "Massage: ");

  va_start(ap, fmt);
  vsprintf(buf + strlen(buf), fmt, ap);
  va_end(ap);

  puts(buf);
  exit(1);
}

bool ALG::ispositivezero(double x, double eps) {
  return ispositive(x, eps) || iszero(x, eps);
}

bool ALG::isnegativezero(double x, double eps) {
  return isnegative(x, eps) || iszero(x, eps);
}
