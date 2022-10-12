//
// Created by 조준홍 on 2022/02/19.
//

#include "readFile.h"

void ALG::readFile::readGeometryInfomation(const std::string &filename, std::vector<section> *sections,
                                           std::vector<region> *regions) {
    std::ifstream f(filename);
    std::string word{}, name{};

    if (!f.is_open()) {
        printError("ALG::readFile::readGeometryInformation",
                   "No Geometry information file: \"%s\"\nPlease check file name", filename.c_str());
    }
    std::cout << "Geometry information file: \"" << filename << "\" open" << std::endl;

    while (!f.eof()) {
        f >> word;
        if (word == "SECTION") {
            readSection(f, sections);
        }
        if (word == "REGION") {
            readRegion(f, regions, sections);
        }
        word = "";
    }
    f.close();
}

void ALG::readFile::readSection(std::ifstream &f, std::vector<section> *sections) {
    std::string word{}, name{}, sx{}, sy{}, ex{}, ey{}, bc{}, bv{};
    auto line{boundaryLine2D()};
    auto start{vector{}}, end{vector{}};

    start.resize(2);
    end.resize(2);

    f >> name;
    std::cout << "SECTION \"" << name << "\" Reading..." << std::endl;
    auto s{section(name)};
    while (word != "ENDSECTION") {
        if (word == "LINE") {
            f >> sx >> sy >> ex >> ey >> bc;
            if (bc != "I") f >> bv;
            start[0] = std::stod(sx);
            start[1] = std::stod(sy);
            end[0] = std::stod(ex);
            end[1] = std::stod(ey);
            line.setStart(start);
            line.setAnEnd(end);
            line.setCondition(bc.c_str()[0]);
            if (line.getCondition() != 'I') line.setBoundaryValue(std::stod(bv));
            line.calcProperties();
            s.emplace_back(line);
        }
        f >> word;
    }
    sections->emplace_back(s);
}

void ALG::readFile::readRegion(std::ifstream &f, std::vector<region> *regions, std::vector<section> *sections) {
    std::string name{}, eps{}, xmin{}, xmax{}, ymin{}, ymax{}, hx{}, hy{}, word{}, sign{}, sName{};
    std::vector<double> sVec{};

    f >> name;
    std::cout << "REGION \"" << name << "\" Reading..." << std::endl;
    f >> eps >> xmin >> ymin >> xmax >> ymax >> hx >> hy;

    auto r{region(name, std::stod(eps), std::stod(xmin), std::stod(xmax), std::stod(ymin), std::stod(ymax),
                  std::stod(hx), std::stod(hy))};


    while (word != "ENDREGION") {
        if (word == "+" || word == "-") {
            sign = word;
            f >> sName;
            if (sign == "+") {
                sVec.emplace_back(1.0E0);
            } else if (sign == "-") {
                sVec.emplace_back(-1.0E0);
            } else {
                printError("ALG::readFile::readRegion", "sign of section (which is %c) error", sign.c_str()[0]);
            }
            for (auto &i: *sections) {
                if (sName == i.getName()) {
                    r.emplace_back(&i);
                    break;
                }
            }
        }
        f >> word;
    }
    r.setSign(sVec);
    regions->emplace_back(r);
}