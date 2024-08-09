#include "Header/pointOthers.h"
#include "Header/readFile.h"

auto main() -> int {
  omp_set_dynamic(0);
  omp_set_num_threads(omp_get_max_threads() / 2);

  auto sections = std::vector<ALG::section>();
  auto regions = std::vector<ALG::region>();
  auto boundary_pts = std::deque<ALG::boundaryPoint>();
  ALG::readFile::readGeometryInfomation("GeoInfo", &sections, &regions);
  std::ofstream f("../AGM2D/ALG_output");
  int ptsnum{};
  for (auto &i : regions) {
    i.setPtsnum(ptsnum);
    i.writeRegion(f, &boundary_pts);
    ptsnum = i.getPtsnum();
  }
  f.close();

  // Copy GeoInfo geometry file to AGM2D diectory
  std::filesystem::path source = "./GeoInfo";
  std::filesystem::path destination_dir = "../AGM2D/";
  std::filesystem::path destination = destination_dir / source;
  try {
    if (!std::filesystem::exists(source)) {
      ALG::printError("GeoInfo file does not exist.");
    }
    if (!std::filesystem::exists(destination_dir)) {
      ALG::printError("AGM2D directory does not exist.");
    }
    std::filesystem::copy_file(source, destination, std::filesystem::copy_options::overwrite_existing);
    std::cout << "\nFile copied successfully.\n";
  } catch (std::filesystem::filesystem_error &e) {
    std::cerr << "Filesystem error: " << e.what() << "\n";
    exit(1);
  } catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    exit(1);
  }
  return 0;
}