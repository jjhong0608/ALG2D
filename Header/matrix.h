//
// Created by 조준홍 on 2022/02/18.
//

#ifndef ALG2D_MATRIX_H
#define ALG2D_MATRIX_H

#include "vector.h"

namespace ALG {
    class matrix : public std::vector<std::vector<double>> {
    private:
        int row{}, col{};
    public:
        matrix(int row, int col);

        virtual ~matrix();

        matrix operator+(const matrix &src);

        matrix operator-(const matrix &src);

        matrix operator*(double d);

        matrix operator*(const matrix &src);

        ALG::vector operator*(const ALG::vector &src);

        matrix &operator=(const matrix &src);
    };

}


#endif //ALG2D_MATRIX_H
