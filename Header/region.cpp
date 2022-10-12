//
// Created by 조준홍 on 2022/02/18.
//

#include "region.h"
#include "pointOthers.h"

ALG::region::region(std::string name) : name(std::move(name)) {}

ALG::region::region(double epsilon, double xmin, double xmax, double ymin, double ymax, double hx, double hy) : epsilon(
        epsilon), xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), hx(hx), hy(hy) {}

ALG::region::region(std::string name, double epsilon, double xmin, double xmax, double ymin, double ymax, double hx,
                    double hy) : name(std::move(name)), epsilon(epsilon), xmin(xmin), xmax(xmax), ymin(ymin),
                                 ymax(ymax), hx(hx), hy(hy) {}

int ALG::region::getPtsnum() const {
    return ptsnum;
}

void ALG::region::setPtsnum(int i) {
    region::ptsnum = i;
}

const std::string &ALG::region::getName() const {
    return name;
}

void ALG::region::setName(const std::string &string) {
    region::name = string;
}

double ALG::region::getEpsilon() const {
    return epsilon;
}

void ALG::region::setEpsilon(double d) {
    region::epsilon = d;
}

double ALG::region::getXmin() const {
    return xmin;
}

void ALG::region::setXmin(double d) {
    region::xmin = d;
}

double ALG::region::getXmax() const {
    return xmax;
}

void ALG::region::setXmax(double d) {
    region::xmax = d;
}

double ALG::region::getYmin() const {
    return ymin;
}

void ALG::region::setYmin(double d) {
    region::ymin = d;
}

double ALG::region::getYmax() const {
    return ymax;
}

void ALG::region::setYmax(double d) {
    region::ymax = d;
}

double ALG::region::getHx() const {
    return hx;
}

void ALG::region::setHx(double d) {
    region::hx = d;
}

double ALG::region::getHy() const {
    return hy;
}

void ALG::region::setHy(double d) {
    region::hy = d;
}

const std::vector<double> &ALG::region::getSign() const {
    return sign;
}

void ALG::region::setSign(const std::vector<double> &s) {
    region::sign = s;
}

void ALG::region::writeRegion(std::ofstream &f, std::deque<boundaryPoint> *boundary_pts) {
    double i{xmin}, j{ymin};
    int xIdx{}, yIdx{}, idx{};
    auto vec{ALG::vector()};
    vec.resize(2);
    // Graph paper outside of Region
    auto xvec{std::vector<double>()}, yvec{std::vector<double>()};

    std::cout << "region: \"" << name << "\"" << "\n";
    std::cout << "xmin = " << xmin << ", xmax = " << xmax << "\n";
    std::cout << "ymin = " << ymin << ", ymax = " << ymax << "\n";
    std::cout << "hx = " << hx << ", hy = " << hy << "\n\n";

//    auto xvec0 = std::vector<double>(), yvec0 = std::vector<double>();
//    auto xvec1 = std::vector<double>(), yvec1 = std::vector<double>();
//    auto xvec2 = std::vector<double>(), yvec2 = std::vector<double>();


    // hx 간격으로 vector를 구성 (x-방향 grid)
    while (i < xmax + hx / 8) {
        xvec.emplace_back(std::stod(std::to_string(i)));
        i += hx;
    }
    // hy 간격으로 vector를 구성 (y-방향 grid)
    while (j < ymax + hy / 8) {
        yvec.emplace_back(std::stod(std::to_string(j)));
        j += hy;
    }

    // for non-uniform BFS flow
//    for (auto &item: yvec) {
//        item = std::copysign(std::pow(2 * (item - HALFVALUE), 2), item - HALFVALUE) * HALFVALUE + HALFVALUE;
//    }
//    for (auto &item: xvec) {
//        item = std::pow(item / 1.5e1, 1.2) * 1.5e1;
//    }
//    yvec.emplace_back(HALFVALUE - HALFVALUE * HALFVALUE * HALFVALUE * HALFVALUE * hy);
//    yvec.emplace_back(HALFVALUE + HALFVALUE * HALFVALUE * HALFVALUE * HALFVALUE * hy);
//    xvec.emplace(xvec.begin(), HALFVALUE * HALFVALUE * HALFVALUE * HALFVALUE * hx);

//    yvec.emplace_back(HALFVALUE - HALFVALUE * HALFVALUE * hy);
//    yvec.emplace_back(HALFVALUE + HALFVALUE * HALFVALUE * hy);
//    xvec.emplace(xvec.begin(), HALFVALUE * HALFVALUE * hx);
//
//    yvec.emplace_back(HALFVALUE - HALFVALUE * HALFVALUE * HALFVALUE * hy);
//    yvec.emplace_back(HALFVALUE + HALFVALUE * HALFVALUE * HALFVALUE * hy);
//    xvec.emplace(xvec.begin(), HALFVALUE *  HALFVALUE * HALFVALUE * hx);
//
//    std::sort(xvec.begin(), xvec.end());
//    std::sort(yvec.begin(), yvec.end());


//    auto theta = std::vector<double>();
//    double st{-5.0E-1}, ed{5.0E-1}, pp{st};
//    theta.emplace_back(-M_PI);
//    int tempidx{};

//    if (name == "region5") {
//        while (theta.back() + M_PI / 90 < 1e-10) {
//            theta.emplace_back(theta.back() + M_PI / 90);
//        }
//
//        for (const auto &item: theta) {
//            xvec1.emplace_back(5.0E-1 * std::cos(item));
//            yvec1.emplace_back(5.0E-1 * std::cos(item));
//        }
//
//        for (auto &item: xvec) {
//            if (item < -5.0E-1 - 1e-10) {
//                xvec0.emplace_back(item);
//            } else if (item > 5.0E-1 + 1e-10) {
//                xvec2.emplace_back(item);
//            }
//        }
//        for (auto &item: yvec) {
//            if (item < -5.0E-1 - 1e-10) {
//                yvec0.emplace_back(item);
//            } else if (item > 5.0E-1 + 1e-10) {
//                yvec2.emplace_back(item);
//            }
//        }
//
//        xvec = xvec0;
//        yvec = yvec0;
//
//        xvec.insert(xvec.end(), xvec1.begin(), xvec1.end());
//        yvec.insert(yvec.end(), yvec1.begin(), yvec1.end());
//
//        xvec.insert(xvec.end(), xvec2.begin(), xvec2.end());
//        yvec.insert(yvec.end(), yvec2.begin(), yvec2.end());
//    }

//    auto xvec{std::vector<double>{0.000000, 0.003500, 0.007800, 0.013900, 0.021700, 0.031200, 0.042500, 0.055600,
//                                  0.070300, 0.086800, 0.105000, 0.125000, 0.146700, 0.170100, 0.195300, 0.222200,
//                                  0.250900, 0.281200, 0.313400, 0.347200, 0.382800, 0.420100, 0.459200, 0.500000,
//                                  0.540800, 0.579900, 0.617200, 0.652800, 0.686600, 0.718800, 0.749100, 0.777800,
//                                  0.804700, 0.829900, 0.853300, 0.875000, 0.895000, 0.913200, 0.929700, 0.944400,
//                                  0.957500, 0.968800, 0.978300, 0.986100, 0.992200, 0.996500, 1.000000}};
//    auto yvec{std::vector<double>{0.000000, 0.003500, 0.007800, 0.013900, 0.021700, 0.031200, 0.042500, 0.055600,
//                                  0.070300, 0.086800, 0.105000, 0.125000, 0.146700, 0.170100, 0.195300, 0.222200,
//                                  0.250900, 0.281200, 0.313400, 0.347200, 0.382800, 0.420100, 0.459200, 0.500000,
//                                  0.540800, 0.579900, 0.617200, 0.652800, 0.686600, 0.718800, 0.749100, 0.777800,
//                                  0.804700, 0.829900, 0.853300, 0.875000, 0.895000, 0.913200, 0.929700, 0.944400,
//                                  0.957500, 0.968800, 0.978300, 0.986100, 0.992200, 0.996500, 1.000000}};

//    std::ifstream xaxial_lines("xaxial.dat");
//    std::ifstream yaxial_lines("yaxial.dat");
//    std::array<double, 4> xaxial{}, yaxial{};
//    while (!xaxial_lines.eof()) {
//        xaxial_lines >> xaxial[0] >> xaxial[1] >> xaxial[2];
//        if (xaxial[0] != xaxial[3]) yvec.emplace_back(xaxial[0]);
//        xaxial[3] = xaxial[0];
//    }
//    while (!yaxial_lines.eof()) {
//        yaxial_lines >> yaxial[0] >> yaxial[1] >> yaxial[2];
//        if (yaxial[0] != yaxial[3]) xvec.emplace_back(yaxial[0]);
//        yaxial[3] = yaxial[0];
//    }

    // 각 방향의 grid point들의 개수
    xIdx = int(xvec.size());
    yIdx = int(yvec.size());

    int xit{}, yit{};
    auto pts{std::vector<point>(xIdx * yIdx)};
    auto pts2D{std::vector<std::array<int, 2>>(xIdx * yIdx)};
    auto xGrid{std::vector<gridLine2D>(yIdx)}, yGrid{std::vector<gridLine2D>(xIdx)};

    // index를 -3으로 grid point들을 저장
    for (const auto &y: yvec) {
        for (const auto &x: xvec) {
            pts.at(idx) = point(ALG::vector(std::vector<double>{x, y}), -3);
            pts2D.at(idx) = std::array<int, 2>{int(idx / xIdx), idx % xIdx};
            ++idx;
        }
    }
    for (const auto &y: yvec) {
        vec[0] = xmin, vec[1] = y;
        xGrid[xit].setStart(vec);
        vec[0] = xmax;
        xGrid[xit].setAnEnd(vec);
        xGrid[xit++].setPts(std::vector<point *>());
    }
    for (const auto &x: xvec) {
        vec[0] = x, vec[1] = ymin;
        yGrid[yit].setStart(vec);
        vec[1] = ymax;
        yGrid[yit].setAnEnd(vec);
        yGrid[yit++].setPts(std::vector<point *>());
    }
    for (int k = 0; k < pts.size(); ++k) {
        xGrid[pts2D.at(k)[0]].getPts().emplace_back(&(pts.at(k)));
        yGrid[pts2D.at(k)[1]].getPts().emplace_back(&(pts.at(k)));
    }

#pragma omp parallel firstprivate(vec)
    {
        int signNum{};
#pragma omp for
        for (auto item = xGrid.begin(); item != xGrid.end(); ++item) {
            item->setCrossPts(std::vector<boundaryPoint>());
            for (auto &sections: *this) {
                for (auto &section: *sections) {
                    if (item->iscross(section, vec)) {
                        item->getCrossSectionPts().emplace_back(
                                boundaryPoint(vec, section.getCondition(), 'B', section.getBoundaryValue(), &section,
                                              section.getNormal()));
                    }
                }
                item->sortCrossSectionPts(sign[signNum++]);
                item->setCrossPts();
            }
            signNum = 0;
            item->sortCrossPts();
            item->makeAxialLine(hx);
        }
#pragma omp for
        for (auto item = yGrid.begin(); item != yGrid.end(); ++item) {
            item->setCrossPts(std::vector<boundaryPoint>());
            for (auto &sections: *this) {
                for (auto &section: *sections) {
                    if (item->iscross(section, vec)) {
                        item->getCrossSectionPts().emplace_back(
                                boundaryPoint(vec, section.getCondition(), 'B', section.getBoundaryValue(), &section,
                                              section.getNormal()));
                    }
                }
                item->sortCrossSectionPts(sign[signNum++]);
                item->setCrossPts();
            }
            signNum = 0;
            item->sortCrossPts();
            item->makeAxialLine(hy);
        }
    }
    for (auto &item: yGrid) {
        for (auto &line: item.getAxialLines()) {
            line.isinBoundary(boundary_pts);
            if (line.getPts().size() > 1) {
                for (auto &pt: line.getPts()) {
                    pt->setIdx(pt->getIdx() + 1);
                }
            }
        }
    }

    f << "REGION " << getName() << "\n\n";
    f.precision(16);
    f << "Material Property = " << std::scientific << getEpsilon() << "\n\n";
    f << "# The number of cross points = ";
    auto s0 = f.tellp();
    f << "0000000000000000";
    auto s1 = f.tellp();
    f << "\n\n";
    int region_ptsnum{};
    for (auto &item: xGrid) {
        for (auto &line: item.getAxialLines()) {
            line.isinBoundary(boundary_pts);
            if (line.getPts().size() > 1) {
                for (auto &pt: line.getPts()) {
                    pt->setIdx(pt->getIdx() + 1);
                    if (pt->getIdx() == -1) {
                        pt->setIdx(ptsnum++);
                        ++region_ptsnum;
                        f << ptsnum - 1 << " ";
                        f << std::scientific << (*pt)[0] << " " << (*pt)[1] << "\n";
                    }
                }
            }
        }
    }
    f << "\n";
    f.seekp(s0);
    f << region_ptsnum;
    auto s2 = f.tellp();
    for (int l = 0; l < s1 - s2; ++l) {
        f << '\0';
    }
    f.seekp(0, std::ios::end);

    f << "# The number of boundary points = ";
    s0 = f.tellp();
    f << "0000000000000000";
    s1 = f.tellp();
    f << "\n\n";
    region_ptsnum = 0;
    for (auto &item: xGrid) {
        for (auto &line: item.getAxialLines()) {
            if (line.getPts().size() > 1) {
                if (std::any_of(line.getPts().begin(), line.getPts().end(),
                                [](point *pt) -> bool { return pt->getIdx() > -1; })) {
                    for (auto &boundary: line.getBoundaries()) {
                        if (boundary->getIdx() == -1) {
                            boundary->setIdx(ptsnum++);
                            ++region_ptsnum;
                            f << ptsnum - 1 << " ";
                            f << std::scientific << (*boundary)[0] << " " << (*boundary)[1] << " "
                              << boundary->getCondition() << " " << boundary->getBoundaryValue() << " "
                              << boundary->getNormal()[0] << " " << boundary->getNormal()[1] << "\n";
                        }
                    }
                }
            }
        }
    }
    for (auto &item: yGrid) {
        for (auto &line: item.getAxialLines()) {
            if (line.getPts().size() > 1) {
                if (std::any_of(line.getPts().begin(), line.getPts().end(),
                                [](point *pt) -> bool { return pt->getIdx() > -1; })) {
                    for (auto &boundary: line.getBoundaries()) {
                        if (boundary->getIdx() == -1) {
                            boundary->setIdx(ptsnum++);
                            ++region_ptsnum;
                            f << ptsnum - 1 << " ";
                            f << std::scientific << (*boundary)[0] << " " << (*boundary)[1] << " "
                              << boundary->getCondition() << " " << boundary->getBoundaryValue() << " "
                              << boundary->getNormal()[0] << " " << boundary->getNormal()[1] << "\n";
                        }
                    }
                }
            }
        }
    }
    f << "\n";
    f.seekp(s0);
    f << region_ptsnum;
    s2 = f.tellp();
    for (int l = 0; l < s1 - s2; ++l) {
        f << '\0';
    }
    f.seekp(0, std::ios::end);

    f << "# the number of x-axial lines = ";
    s0 = f.tellp();
    f << "0000000000000000";
    s1 = f.tellp();
    f << "\n\n";
    region_ptsnum = 0;
    for (auto &grid: xGrid) {
        for (auto &item: grid.getAxialLines()) {
            if (item.getPts().size() > 1) {
                if (std::any_of(item.getPts().begin(), item.getPts().end(),
                                [](point *pt) -> bool { return pt->getIdx() > -1; })) {
                    item.setIdx(region_ptsnum);
                    f << item.getBoundaries()[0]->getIdx() << " ";
                    for (auto &pt: item.getPts()) {
                        if (pt->getIdx() > -1) {
                            f << pt->getIdx() << " ";
                        }
                    }
                    f << item.getBoundaries()[1]->getIdx() << " /\n";
                    ++region_ptsnum;
                }
            }
        }
    }
    f << "\n";
    f.seekp(s0);
    f << region_ptsnum;
    s2 = f.tellp();
    for (int l = 0; l < s1 - s2; ++l) {
        f << '\0';
    }
    f.seekp(0, std::ios::end);

    f << "# the number of y-axial lines = ";
    s0 = f.tellp();
    f << "0000000000000000";
    s1 = f.tellp();
    f << "\n\n";
    region_ptsnum = 0;
    for (auto &grid: yGrid) {
        for (auto &item: grid.getAxialLines()) {
            if (item.getPts().size() > 1) {
                if (std::any_of(item.getPts().begin(), item.getPts().end(),
                                [](point *pt) -> bool { return pt->getIdx() > -1; })) {
                    item.setIdx(region_ptsnum);
                    f << item.getBoundaries()[0]->getIdx() << " ";
                    for (auto &pt: item.getPts()) {
                        if (pt->getIdx() > -1) {
                            f << pt->getIdx() << " ";
                        }
                    }
                    f << item.getBoundaries()[1]->getIdx() << " /\n";
                    ++region_ptsnum;
                }
            }
        }
    }
    f << "\n";
    f.seekp(s0);
    f << region_ptsnum;
    s2 = f.tellp();
    for (int l = 0; l < s1 - s2; ++l) {
        f << '\0';
    }
    f.seekp(0, std::ios::end);
    f << "ENDREGION\n\n";
}
