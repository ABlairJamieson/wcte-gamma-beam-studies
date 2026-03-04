#include <iostream>
#include <string>

enum class ConeVertexMode { Nominal, Truth }; // keep in a header if you prefer
int RunAnalysis(const char* fname, ConeVertexMode coneMode);

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <wcsim.root> [--cone-vertex=nominal|truth]\n";
    return 1;
  }

  ConeVertexMode mode = ConeVertexMode::Nominal;
  if (argc >= 3) {
    std::string opt = argv[2];
    if (opt == "--cone-vertex=truth") mode = ConeVertexMode::Truth;
    if (opt == "--cone-vertex=nominal") mode = ConeVertexMode::Nominal;
  }

  return RunAnalysis(argv[1], mode);
}