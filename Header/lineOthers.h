//
// Created by 조준홍 on 2022/02/18.
//

#ifndef ALG2D_LINEOTHERS_H
#define ALG2D_LINEOTHERS_H

#include "line.h"
#include "pointOthers.h"

namespace ALG {
    class boundaryPoint;

    class axialLine2D : public line2D {
    private:
        int idx{-1};
        std::vector<point *> pts{};
        std::array<boundaryPoint *, 2> boundaries{};
    public:
        axialLine2D(const boundaryPoint &sp, const boundaryPoint &ep);

        int getIdx() const;

        void setIdx(int i);

        [[nodiscard]] std::vector<point *> &getPts();

        [[nodiscard]] const std::vector<point *> &getPts() const;

        void setPts(const std::vector<point *> &vector);

        [[nodiscard]] std::array<boundaryPoint *, 2> &getBoundaries();

        [[nodiscard]] const std::array<boundaryPoint *, 2> &getBoundaries() const;

        void setBoundaries(const std::array<boundaryPoint *, 2> &array);

        void isinBoundary(std::deque<boundaryPoint> *boundaryPts);
    };

    class gridLine2D : public line2D {
    private:
        std::vector<point *> pts{};
        std::vector<boundaryPoint> cross_pts{}, cross_section_pts{};
        std::vector<axialLine2D> axial_lines{};
    public:
        [[nodiscard]] std::vector<point *> &getPts();

        [[nodiscard]] const std::vector<point *> &getPts() const;

        void setPts(const std::vector<point *> &vector);

        [[nodiscard]] std::vector<boundaryPoint> &getCrossPts();

        [[nodiscard]] const std::vector<boundaryPoint> &getCrossPts() const;

        void setCrossPts();

        void setCrossPts(const std::vector<boundaryPoint> &crossPts);

        [[nodiscard]] std::vector<boundaryPoint> &getCrossSectionPts();

        [[nodiscard]] const std::vector<boundaryPoint> &getCrossSectionPts() const;

        void setCrossSectionPts(const std::vector<boundaryPoint> &crossSectionPts);

        [[nodiscard]] std::vector<axialLine2D> &getAxialLines();

        [[nodiscard]] const std::vector<axialLine2D> &getAxialLines() const;

        void setAxialLines(const std::vector<axialLine2D> &axialLines);

        void sortCrossSectionPts(double sign);

        void sortCrossPts();

        void makeAxialLine();

        void makeAxialLine(double h);
    };

    class boundaryLine2D : public line2D {
    private:
        char condition{};
        double boundary_value{};
        std::vector<boundaryPoint *> boundary_points{};
    public:
        boundaryLine2D();

        boundaryLine2D(const vector &start, const vector &anEnd, char condition, double boundaryValue);

        boundaryLine2D(double s0, double s1, double e0, double e1, char condition, double boundaryValue);

        [[nodiscard]] char getCondition() const;

        void setCondition(char i);

        [[nodiscard]] double getBoundaryValue() const;

        void setBoundaryValue(double boundaryValue);

        [[nodiscard]] std::vector<boundaryPoint *> &getBoundaryPoints();

        [[nodiscard]] const std::vector<boundaryPoint *> &getBoundaryPoints() const;

        void setBoundaryPoints(const std::vector<boundaryPoint *> &boundaryPoints);
    };
}


#endif //ALG2D_LINEOTHERS_H
