//
// Created by 조준홍 on 2022/02/18.
//

#include "pointOthers.h"
#include "lineOthers.h"

ALG::boundaryPoint::boundaryPoint() = default;

ALG::boundaryPoint::boundaryPoint(char io) : IO(io) {}

ALG::boundaryPoint::boundaryPoint(ALG::boundaryLine2D *line) : line(line) {}

ALG::boundaryPoint::boundaryPoint(char condition, char io, double boundaryValue) : condition(condition), IO(io),
                                                                                   boundary_value(boundaryValue) {}

ALG::boundaryPoint::boundaryPoint(char condition, char io, double boundaryValue, ALG::boundaryLine2D *line) : condition(
        condition), IO(io), boundary_value(boundaryValue), line(line) {}

ALG::boundaryPoint::boundaryPoint(char condition, char io, double boundaryValue, ALG::boundaryLine2D *line,
                                  const ALG::vector &normal) : condition(condition), IO(io),
                                                               boundary_value(boundaryValue), line(line),
                                                               normal(normal) {}

ALG::boundaryPoint::boundaryPoint(const std::vector<double> &x, char io) : point(x), IO(io) {}

ALG::boundaryPoint::boundaryPoint(const std::vector<double> &x, ALG::boundaryLine2D *line) : point(x), line(line) {}

ALG::boundaryPoint::boundaryPoint(const std::vector<double> &x, char condition, char io, double boundaryValue) :
        point(x), condition(condition), IO(io), boundary_value(boundaryValue) {}

ALG::boundaryPoint::boundaryPoint(const std::vector<double> &x, char condition, char io, double boundaryValue,
                                  ALG::boundaryLine2D *line) : point(x), condition(condition), IO(io),
                                                               boundary_value(boundaryValue), line(line) {}

ALG::boundaryPoint::boundaryPoint(const std::vector<double> &x, char condition, char io, double boundaryValue,
                                  ALG::boundaryLine2D *line, const ALG::vector &normal) : point(x),
                                                                                          condition(condition), IO(io),
                                                                                          boundary_value(boundaryValue),
                                                                                          line(line),
                                                                                          normal(normal) {}

char ALG::boundaryPoint::getCondition() const {
    return condition;
}

void ALG::boundaryPoint::setCondition(char i) {
    boundaryPoint::condition = i;
}

char ALG::boundaryPoint::getIo() const {
    return IO;
}

void ALG::boundaryPoint::setIo(char io) {
    IO = io;
}

double ALG::boundaryPoint::getBoundaryValue() const {
    return boundary_value;
}

void ALG::boundaryPoint::setBoundaryValue(double boundaryValue) {
    boundary_value = boundaryValue;
}

ALG::vector &ALG::boundaryPoint::getNormal() {
    return normal;
}

const ALG::vector &ALG::boundaryPoint::getNormal() const {
    return normal;
}

void ALG::boundaryPoint::setNormal(const ALG::vector &src) {
    boundaryPoint::normal = src;
}

ALG::boundaryLine2D *ALG::boundaryPoint::getLine() const {
    return line;
}

void ALG::boundaryPoint::setLine(ALG::boundaryLine2D *pLine2D) {
    boundaryPoint::line = pLine2D;
}

ALG::boundaryPoint *ALG::boundaryPoint::isinBoundary() {
    for (const auto &item: line->getBoundaryPoints()) {
        if (iszero((*this - *item).norm())) {
            return item;
        }
    }
    return nullptr;
}

