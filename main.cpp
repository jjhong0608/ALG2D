#include "Header/pointOthers.h"
#include "Header/readFile.h"

int main() {
    auto sections = std::vector<ALG::section>();
    auto regions = std::vector<ALG::region>();
    auto boundary_pts = std::deque<ALG::boundaryPoint>();
    ALG::readFile::readGeometryInfomation("GeoInfo", &sections, &regions);
    std::ofstream f("/home/jjhong0608/AGM2D/ALG_output");
    int ptsnum{};
    for (auto &i : regions) {
        i.setPtsnum(ptsnum);
        i.writeRegion(f, &boundary_pts);
        ptsnum = i.getPtsnum();
    }
    f.close();
    system("cp GeoInfo /home/jjhong0608/AGM2D/");
    system("cp /home/jjhong0608/AGM2D/ALG_output /home/jjhong0608/docker/ALG2D/");
    return 0;
}
