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

    auto xvec0 = std::vector<double>(), yvec0 = std::vector<double>();
    auto xvec1 = std::vector<double>(), yvec1 = std::vector<double>();
    auto xvec2 = std::vector<double>(), yvec2 = std::vector<double>();

//    auto xvec0{std::vector<double>{
//            30.02695215478053, 30.071127162890576, 30.12547674631095, 30.187705708209894, 30.256537878024204,
//            30.33113511875069, 30.41089489735123, 30.495358333760368, 30.584161600175857, 30.677007518931653,
//            30.773647722293866, 30.87387081745252, 30.977494191219762, 31.08435813510233, 31.19432151166049,
//            31.307258480128173, 31.42305597124594, 31.541611705203593, 31.662832611758766, 31.786633553766094,
//            31.912936283396743, 32.041668579436376, 32.17276352735063, 32.306158913244346, 32.44179670965294,
//            32.57962263609805, 32.71958578105265, 32.861638274756935, 33.005735004457605, 33.151833365284524,
//            33.29989304125615, 33.44987581190815, 33.601745380834316, 33.75546722306324, 33.911008448704635,
//            34.06833768071247, 34.22742494494943, 34.3882415710137, 34.550760102517074, 34.714954215693304,
//            34.88079864537282, 35.048269117493, 35.21734228742411, 35.387995683485514, 35.56020765510671,
//            35.73395732515596, 35.90922454601782, 36.08598985905068, 36.264234457099064, 36.443940149772544,
//            36.62508933123549, 36.80766495028042, 36.99165048248173, 37.17702990424892, 37.36378766861644,
//            37.55190868262455, 37.74137828615989, 37.932182232137606, 38.12430666791807, 38.31773811786176,
//            38.5124634669345, 38.70846994528383, 38.905745113713984, 39.10427684999373, 39.30405333593695,
//            39.50506304520121, 39.70729473175397, 39.910737418960515, 40.11538038925143, 40.32121317433078,
//            40.528225545889526, 40.73640750679109, 40.945749282699055, 41.15624131411889, 41.36787424882798,
//            41.58063893467009, 41.79452641269198, 42.0095279106019, 42.22563483653082, 42.44283877307865,
//            42.66113147162913, 42.880504846918114, 43.10095097184083, 43.32246207248501, 43.545030523377406,
//            43.76864884293212, 43.99330968908996, 44.21900585513869, 44.44573026570458, 44.67347597290661,
//            44.902236152664806, 45.132004101154955, 45.36277323140236, 45.59453707000776, 45.827289253998956,
//            46.06102352780196, 46.29573374032595, 46.531413842156724, 46.768057882853455, 47.00566000834388,
//            47.24421445841355, 47.483715564284694, 47.72415774628075, 47.96553551157274, 48.207843452003594,
//            48.45107624198752, 48.69522863648065, 48.94029546902016, 49.18627164982887, 49.43315216398257,
//            49.68093206963747, 49.929606496315195, 50.17917064324298, 50.42961977774697, 50.68094923369621,
//            50.93315440999553, 51.18623076912522, 51.440173835725744, 51.69497919522578, 51.95064249251167,
//            52.20715943063701, 52.46452576957065, 52.72273732498166, 52.98178996706001, 53.24167961937141,
//            53.50240225774532, 53.7639539091947, 54.026330650866505, 54.28952860902163, 54.55354395804348,
//            54.81837291947399, 55.084011761076184, 55.35045679592233, 55.61770438150681, 55.88575091888293,
//            56.15459285182268, 56.424226665998816, 56.69464888818851, 56.96585608549772, 57.23784486460576,
//            57.510611871029276, 57.78415378840501, 58.0584673377909, 58.33354927698466, 58.60939639985949,
//            58.886005535716265, 59.163373548651805, 59.44149733694252, 59.720373832443144, 60.0
//    }};


//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution<> dis(-hx * 3.0 / 8.0, hx * 3.0 / 8.0);

    // hx 간격으로 vector를 구성 (x-방향 grid)
//    while (i < double(30) + hx / 8) {
    while (i < xmax + hx / 8) {
//        if (isnegative(i - HALFVALUE)) {
//            xvec.emplace_back(std::stod(std::to_string(2 * i * i)));
//        } else if (ispositive(i - HALFVALUE)) {
//            xvec.emplace_back(std::stod(std::to_string(-2 * (i - UNITVALUE) * (i - UNITVALUE) + UNITVALUE)));
//        } else {
//            xvec.emplace_back(std::stod(std::to_string(i)));
//        }
//        xvec.emplace_back(i + dis(gen));

        xvec.emplace_back(std::stod(std::to_string(i)));
        i += hx;
    }

//    for (const auto &item: xvec0) {
//        xvec.emplace_back(item);
//    }

//    xvec.emplace_back(hx * HALFVALUE);
//    xvec.emplace_back(UNITVALUE - hx * HALFVALUE);
//    xvec.emplace_back(hx * HALFVALUE * HALFVALUE);
//    xvec.emplace_back(UNITVALUE - hx * HALFVALUE * HALFVALUE);
//    std::sort(xvec.begin(), xvec.end());
    // hy 간격으로 vector를 구성 (y-방향 grid)
    while (j < ymax + hy / 8) {
//        if (isnegative(j - HALFVALUE)) {
//            yvec.emplace_back(std::stod(std::to_string(2 * j * j)));
//        } else if (ispositive(j - HALFVALUE)) {
//            yvec.emplace_back(std::stod(std::to_string(-2 * (j - UNITVALUE) * (j - UNITVALUE) + UNITVALUE)));
//        } else {
//            yvec.emplace_back(std::stod(std::to_string(j)));
//        }
//        yvec.emplace_back(j + dis(gen));

        yvec.emplace_back(std::stod(std::to_string(j)));
        j += hy;
    }
//    yvec.emplace_back(hy * HALFVALUE);
//    yvec.emplace_back(UNITVALUE - hy * HALFVALUE);
//    yvec.emplace_back(hy * HALFVALUE * HALFVALUE);
//    yvec.emplace_back(UNITVALUE - hy * HALFVALUE * HALFVALUE);
//    std::sort(yvec.begin(), yvec.end());

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


    auto theta = std::vector<double>();
    double st{-5.0E-1}, ed{5.0E-1}, pp{st};
    theta.emplace_back(-M_PI);
//    int tempidx{};

    if (name == "region5") {
        while (theta.back() + M_PI / 90 < 1e-10) {
            theta.emplace_back(theta.back() + M_PI / 90);
        }

        for (const auto &item: theta) {
            xvec1.emplace_back(5.0E-1 * std::cos(item));
            yvec1.emplace_back(5.0E-1 * std::cos(item));
        }

        for (auto &item: xvec) {
            if (item < -5.0E-1 - 1e-10) {
                xvec0.emplace_back(item);
            } else if (item > 5.0E-1 + 1e-10) {
                xvec2.emplace_back(item);
            }
        }
        for (auto &item: yvec) {
            if (item < -5.0E-1 - 1e-10) {
                yvec0.emplace_back(item);
            } else if (item > 5.0E-1 + 1e-10) {
                yvec2.emplace_back(item);
            }
        }

        xvec = xvec0;
        yvec = yvec0;

        xvec.insert(xvec.end(), xvec1.begin(), xvec1.end());
        yvec.insert(yvec.end(), yvec1.begin(), yvec1.end());

        xvec.insert(xvec.end(), xvec2.begin(), xvec2.end());
        yvec.insert(yvec.end(), yvec2.begin(), yvec2.end());
    }

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
