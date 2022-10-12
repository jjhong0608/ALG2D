//
// Created by 조준홍 on 2022/02/18.
//

#ifndef ALG2D_VECTOR_H
#define ALG2D_VECTOR_H


#include "util.h"

namespace ALG {
    class vector : public std::vector<double> {
    public:
        vector();

        explicit vector(const std::vector<double> &x);

        virtual ~vector();

        double norm();

        double dot(const vector &src);

        vector cross(const vector &src);

        vector unitVector();

        vector operator+(const vector &src);

        vector operator-(const vector &src);

        vector operator*(double d);

        vector operator/(double d);

        vector operator+(const vector &src) const;

        vector operator-(const vector &src) const;

        vector operator*(double d) const;

        vector operator/(double d) const;

        double operator*(const vector &src);

        vector &operator+=(const vector &src);

        vector &operator-=(const vector &src);

        vector &operator*=(double d);

        bool operator<(const vector &src);

        bool operator>(const vector &src);
    };
}


#endif //ALG2D_VECTOR_H
