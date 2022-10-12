//
// Created by 조준홍 on 2022/02/18.
//

#ifndef ALG2D_POINTOTHERS_H
#define ALG2D_POINTOTHERS_H

#include "point.h"

namespace ALG {
    class boundaryLine2D;

    class boundaryPoint : public point {
    private:
        char condition{}, IO{};
        double boundary_value{};
        boundaryLine2D *line{};
        vector normal{};
    public:
        boundaryPoint();

        explicit boundaryPoint(char io);

        explicit boundaryPoint(boundaryLine2D *line);

        boundaryPoint(char condition, char io, double boundaryValue);

        boundaryPoint(char condition, char io, double boundaryValue, boundaryLine2D *line);

        boundaryPoint(char condition, char io, double boundaryValue, boundaryLine2D *line, const vector& normal);

        boundaryPoint(const std::vector<double> &x, char io);

        boundaryPoint(const std::vector<double> &x, boundaryLine2D *line);

        boundaryPoint(const std::vector<double> &x, char condition, char io, double boundaryValue);

        boundaryPoint(const std::vector<double> &x, char condition, char io, double boundaryValue,
                      boundaryLine2D *line);

        boundaryPoint(const std::vector<double> &x, char condition, char io, double boundaryValue, boundaryLine2D *line,
                      const vector& normal);

        [[nodiscard]] char getCondition() const;

        void setCondition(char i);

        [[nodiscard]] char getIo() const;

        void setIo(char io);

        [[nodiscard]] double getBoundaryValue() const;

        void setBoundaryValue(double boundaryValue);

        [[nodiscard]] vector &getNormal();

        [[nodiscard]] const vector &getNormal() const;

        void setNormal(const vector &src);

        [[nodiscard]] boundaryLine2D *getLine() const;

        void setLine(boundaryLine2D *pLine2D);

        ALG::boundaryPoint *isinBoundary();

    };
}


#endif //ALG2D_POINTOTHERS_H
