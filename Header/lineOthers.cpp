//
// Created by 조준홍 on 2022/02/18.
//

#include "lineOthers.h"

ALG::axialLine2D::axialLine2D(const ALG::boundaryPoint &sp, const ALG::boundaryPoint &ep) :
        line2D(sp.at(0), sp.at(1), ep.at(0), ep.at(1)) {}

int ALG::axialLine2D::getIdx() const {
    return idx;
}

void ALG::axialLine2D::setIdx(int i) {
    axialLine2D::idx = i;
}

std::vector<ALG::point *> &ALG::axialLine2D::getPts() {
    return pts;
}

const std::vector<ALG::point *> &ALG::axialLine2D::getPts() const {
    return pts;
}

void ALG::axialLine2D::setPts(const std::vector<point *> &vector) {
    axialLine2D::pts = vector;
}

std::array<ALG::boundaryPoint *, 2> &ALG::axialLine2D::getBoundaries() {
    return boundaries;
}

const std::array<ALG::boundaryPoint *, 2> &ALG::axialLine2D::getBoundaries() const {
    return boundaries;
}

void ALG::axialLine2D::setBoundaries(const std::array<ALG::boundaryPoint *, 2> &array) {
    axialLine2D::boundaries = array;
}

void ALG::axialLine2D::isinBoundary(std::deque<boundaryPoint> *boundaryPts) {
    auto checkVertice = [](boundaryPoint *point) -> bool {
        auto check = [](boundaryPoint &point1, std::vector<double> &vector1) -> bool {
            return iszero(std::sqrt((point1[0] - vector1[0]) * (point1[0] - vector1[0]) +
                                    (point1[1] - vector1[1]) * (point1[1] - vector1[1])));
        };
        return check(*point, point->getLine()->getStart()) || check(*point, point->getLine()->getAnEnd());
    };

    bool tf{true};
    for (auto &i: boundaries[0]->getLine()->getBoundaryPoints()) {
        if ((*boundaries[0] - *i).norm() < 1.0E-13) {
            boundaries[0] = i;
            tf = false;
            break;
        }
    }
    if (checkVertice(boundaries[0])) {
        for (auto &item: *boundaryPts) {
            if ((*boundaries[0] - item).norm() < 1.0E-13) {
                boundaries[0] = &item;
                tf = false;
                break;
            }
        }
    }
    if (tf) {
        boundaryPts->emplace_back(*boundaries[0]);
        boundaries[0] = &(boundaryPts->back());
        boundaries[0]->getLine()->getBoundaryPoints().emplace_back(boundaries[0]);
    }
    for (auto &i: boundaries[1]->getLine()->getBoundaryPoints()) {
        if ((*boundaries[1] - *i).norm() < 1.0E-13) {
            boundaries[1] = i;
            return;
        }
    }
    if (checkVertice(boundaries[1])) {
        for (auto &item: *boundaryPts) {
            if ((*boundaries[1] - item).norm() < 1.0E-13) {
                boundaries[1] = &item;
                return;
            }
        }
    }
    boundaryPts->emplace_back(*boundaries[1]);
    boundaries[1] = &(boundaryPts->back());
    boundaries[1]->getLine()->getBoundaryPoints().emplace_back(boundaries[1]);
}

std::vector<ALG::point *> &ALG::gridLine2D::getPts() {
    return pts;
}

const std::vector<ALG::point *> &ALG::gridLine2D::getPts() const {
    return pts;
}

void ALG::gridLine2D::setPts(const std::vector<ALG::point *> &vector) {
    gridLine2D::pts = vector;
}

std::vector<ALG::boundaryPoint> &ALG::gridLine2D::getCrossPts() {
    return cross_pts;
}

const std::vector<ALG::boundaryPoint> &ALG::gridLine2D::getCrossPts() const {
    return cross_pts;
}

void ALG::gridLine2D::setCrossPts() {
    cross_pts.insert(std::end(cross_pts), std::begin(cross_section_pts), std::end(cross_section_pts));
    cross_section_pts.clear();
}

void ALG::gridLine2D::setCrossPts(const std::vector<boundaryPoint> &crossPts) {
    cross_pts = crossPts;
}

std::vector<ALG::boundaryPoint> &ALG::gridLine2D::getCrossSectionPts() {
    return cross_section_pts;
}

const std::vector<ALG::boundaryPoint> &ALG::gridLine2D::getCrossSectionPts() const {
    return cross_section_pts;
}

void ALG::gridLine2D::setCrossSectionPts(const std::vector<boundaryPoint> &crossSectionPts) {
    cross_section_pts = crossSectionPts;
}

std::vector<ALG::axialLine2D> &ALG::gridLine2D::getAxialLines() {
    return axial_lines;
}

const std::vector<ALG::axialLine2D> &ALG::gridLine2D::getAxialLines() const {
    return axial_lines;
}

void ALG::gridLine2D::setAxialLines(const std::vector<axialLine2D> &axialLines) {
    axial_lines = axialLines;
}

void ALG::gridLine2D::sortCrossSectionPts(double sign) {
    char io = sign < 0.0E0 ? 'O' : 'I';
    if (cross_section_pts.size() % 2 != 0) {
        std::sort(cross_section_pts.begin(), cross_section_pts.end());
        cross_section_pts.erase(std::unique(cross_section_pts.begin(), cross_section_pts.end(),
                                            [](const point &a, const point &b) { return iszero((a - b).norm()); }),
                                cross_section_pts.end());
        if (cross_section_pts.size() % 2 != 0) {
            printf("start = (%f, %f)\n", getStart()[0], getStart()[1]);
            printf("anend = (%f, %f)\n", getAnEnd()[0], getAnEnd()[1]);
            for (auto & cross_section_pt : cross_section_pts) {
                printf("(x, y) = (%f, %f)\n", cross_section_pt[0], cross_section_pt[1]);
            }
//            cross_section_pts.erase(cross_section_pts.begin());
            printError("ALG::gridLine2D::sortCrossSectionPts",
                       "size of cross_section_pts (which is %d) is not even number",
                       cross_section_pts.size());
        }
    } else {
        std::sort(cross_section_pts.begin(), cross_section_pts.end());
        auto temp_pts = cross_section_pts;
        cross_section_pts.erase(std::unique(cross_section_pts.begin(), cross_section_pts.end(),
                                            [](const point &a, const point &b) { return iszero((a - b).norm()); }),
                                cross_section_pts.end());
        if (cross_section_pts.size() % 2 == 1) {
            cross_section_pts = temp_pts;
        }
    }
    for (auto &item: cross_section_pts) {
        item.setIo(io);
        if (io == 'I') {
            io = 'O';
            if (getTangent() * item.getNormal() > ZEROVALUE) {
                item.setNormal(item.getNormal() * -UNITVALUE);
            }
        } else {
            io = 'I';
            if (getTangent() * item.getNormal() < ZEROVALUE) {
                item.setNormal(item.getNormal() * -UNITVALUE);
            }
        }
    }
}

void ALG::gridLine2D::sortCrossPts() {
    std::sort(cross_pts.begin(), cross_pts.end());
}

void ALG::gridLine2D::makeAxialLine() {
    int countIO{};
    boundaryPoint *sp{}, *ep{};
    if (cross_pts.empty()) return;
    for (auto &i: cross_pts) {
        if (i.getIo() == 'I') {
            if (countIO == 0) sp = &i;
            countIO++;
        } else if (i.getIo() == 'O') {
            countIO--;
            if (countIO == 0) {
                ep = &i;
                if (isclose((*sp)[0], (*ep)[0])) (*sp)[0] = (*ep)[0] = pts[0]->at(0);
                if (isclose((*sp)[1], (*ep)[1])) (*sp)[1] = (*ep)[1] = pts[0]->at(1);
                axial_lines.emplace_back(axialLine2D(*sp, *ep));
                axial_lines.back().getBoundaries()[0] = sp;
                axial_lines.back().getBoundaries()[1] = ep;
            }
        }
    }

    if (axial_lines.size() > 1) {
        for (int i = 0; i < axial_lines.size() - 1; ++i) {
            if ((axial_lines[i].getAnEnd() - axial_lines[i + 1].getStart()).norm() < 1.0E-10) {
                axial_lines[i + 1].getBoundaries()[0] = axial_lines[i].getBoundaries()[1];
            }
        }
    }

    auto it = axial_lines.begin();
    for (const auto &i: pts) {
        if ((*it).getStart() < *i && *i < (*it).getAnEnd()) {
            it->getPts().emplace_back(i);
        } else if (*i >= (*it).getAnEnd()) {
            it++;
            if (it == axial_lines.end()) {
                break;
            }
        }
    }

    for (auto &i: axial_lines) {
        if (!i.getPts().empty()) {
            if ((*(i.getBoundaries()[0]) - *(i.getPts().front())).norm() < 1.0E-10) {
                i.getPts().erase(i.getPts().begin());
            }
            if ((*(i.getBoundaries()[1]) - *(i.getPts().back())).norm() < 1.0E-10) {
                i.getPts().pop_back();
            }
        }
    }
}

void ALG::gridLine2D::makeAxialLine(double h) {
    int countIO{};
    boundaryPoint *sp{}, *ep{};
    if (cross_pts.empty()) return;
    for (auto &i: cross_pts) {
        if (i.getIo() == 'I') {
            if (countIO == 0) sp = &i;
            countIO++;
        } else if (i.getIo() == 'O') {
            countIO--;
            if (countIO == 0) {
                ep = &i;
                if (isclose((*sp)[0], (*ep)[0])) (*sp)[0] = (*ep)[0] = pts[0]->at(0);
                if (isclose((*sp)[1], (*ep)[1])) (*sp)[1] = (*ep)[1] = pts[0]->at(1);
                axial_lines.emplace_back(axialLine2D(*sp, *ep));
                axial_lines.back().getBoundaries()[0] = sp;
                axial_lines.back().getBoundaries()[1] = ep;
            }
        }
    }

    if (axial_lines.size() > 1) {
        for (int i = 0; i < axial_lines.size() - 1; ++i) {
            if ((axial_lines[i].getAnEnd() - axial_lines[i + 1].getStart()).norm() < 1.0E-10) {
                axial_lines[i + 1].getBoundaries()[0] = axial_lines[i].getBoundaries()[1];
            }
        }
    } else if (axial_lines.empty()){
        return;
    }

    auto it = axial_lines.begin();
    for (const auto &i: pts) {
        if (it->getStart() < *i && *i < it->getAnEnd()) {
            it->getPts().emplace_back(i);
        } else if (*i >= it->getAnEnd()) {
            it++;
            if (it == axial_lines.end()) {
                break;
            }
        }
    }

    for (auto &i: axial_lines) {
        if (!i.getPts().empty()) {
            if ((*(i.getBoundaries()[0]) - *(i.getPts().front())).norm() < h / 20.0) {
                i.getPts().front()->setIdx(-3);
                i.getPts().erase(i.getPts().begin());
            }
            if (i.getPts().empty()) {
                continue;
            }
            if ((*(i.getBoundaries()[1]) - *(i.getPts().back())).norm() < h / 20.0) {
                i.getPts().back()->setIdx(-3);
                i.getPts().pop_back();
            }
        }
    }
}

ALG::boundaryLine2D::boundaryLine2D() = default;

ALG::boundaryLine2D::boundaryLine2D(double s0, double s1, double e0, double e1, char condition, double boundaryValue)
        : line2D(s0, s1, e0, e1), condition(condition), boundary_value(boundaryValue) {}

ALG::boundaryLine2D::boundaryLine2D(const ALG::vector &start, const ALG::vector &anEnd, char condition,
                                    double boundaryValue) : line2D(start, anEnd), condition(condition),
                                                            boundary_value(boundaryValue) {}

char ALG::boundaryLine2D::getCondition() const {
    return condition;
}

void ALG::boundaryLine2D::setCondition(char i) {
    boundaryLine2D::condition = i;
}

double ALG::boundaryLine2D::getBoundaryValue() const {
    return boundary_value;
}

void ALG::boundaryLine2D::setBoundaryValue(double boundaryValue) {
    boundary_value = boundaryValue;
}

std::vector<ALG::boundaryPoint *> &ALG::boundaryLine2D::getBoundaryPoints() {
    return boundary_points;
}

const std::vector<ALG::boundaryPoint *> &ALG::boundaryLine2D::getBoundaryPoints() const {
    return boundary_points;
}

void ALG::boundaryLine2D::setBoundaryPoints(const std::vector<boundaryPoint *> &boundaryPoints) {
    boundary_points = boundaryPoints;
}
