//
// Created by 조준홍 on 2022/02/18.
//

#include "matrix.h"

ALG::matrix::matrix(int row, int col) : vector(row),
                                        row(row),
                                        col(col) {
    for (auto &i: *this) {
        i = std::vector<double>(col);
    }
}

ALG::matrix::~matrix() = default;

ALG::matrix ALG::matrix::operator+(const ALG::matrix &src) {
    if (row != src.row || col != src.col)
        printError("ALG::matrix::operator+",
                   "size of matrix (which is (%d, %d)) is not equal size of source matrix (which is (%d, %d))", row,
                   col, src.row, src.col);
    matrix matrix(row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            matrix[i][j] = at(i).at(j) + src.at(i).at(j);
        }
    }
    return matrix;
}

ALG::matrix ALG::matrix::operator-(const ALG::matrix &src) {
    if (row != src.row || col != src.col)
        printError("ALG::matrix::operator-",
                   "the size of matrix (which is (%d, %d)) is not equal to the size of source matrix (which is (%d, %d))",
                   row, col, src.row, src.col);
    matrix mat = matrix(row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            mat[i][j] = at(i).at(j) - src.at(i).at(j);
        }
    }
    return mat;
}

ALG::matrix ALG::matrix::operator*(double d) {
    matrix mat = matrix(row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            mat[i][j] = at(i).at(j) * d;
        }
    }
    return mat;
}

ALG::matrix ALG::matrix::operator*(const ALG::matrix &src) {
    if (col != src.row)
        printError("ALG::matrix::operator*",
                   "the column size of matrix (which is %d) is not equal to the row size of source matrix (which is %d)",
                   col, src.row);
    matrix mat = matrix(row, src.col);
    double sum;
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < src.col; ++j) {
            sum = ZEROVALUE;
            for (int k = 0; k < col; ++k) {
                sum += at(i).at(k) * src.at(k).at(j);
            }
            mat[i][j] = sum;
        }
    }
    return mat;
}

ALG::vector ALG::matrix::operator*(const ALG::vector &src) {
    if (col != src.size())
        printError("ALG::matrix::operator*",
                   "the column size of matrix (which is %d) is not equal to the size of source vector (which is %d)",
                   col, src.size());
    ALG::vector vec = ALG::vector{};
    double sum{};
    for (int i = 0; i < row; ++i) {
        sum = ZEROVALUE;
        for (int j = 0; j < col; ++j) {
            sum += at(i).at(j) * src.at(j);
        }
        vec.emplace_back(sum);
    }
    return vec;
}

ALG::matrix &ALG::matrix::operator=(const ALG::matrix &src) {
    row = src.row;
    col = src.col;
    resize(row);
    for (auto &i: *this) {
        i.resize(col);
    }
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            at(i).at(j) = src.at(i).at(j);
        }
    }
    return *this;
}

