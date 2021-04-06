/*
  HSX Hashed Sequence Index
  Spec: http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/hsx_format.html
  Lastz documentation http://www.bx.psu.edu/~rsharris/lastz/README.lastz-1.04.03.html#fmt_hsx
*/

#include <iostream>
#include "hsx.cpp"

int main(int argc, char** argv) {
  std::vector<std::string> filenames;

  for (int i = 1; i < argc; i++) {
    std::string f(argv[i]);
    filenames.push_back(f);
  }

  std::string hsx_file = "/home/sluggie/src/bio/hsx/out.hsx";

  generate_hsx(filenames, hsx_file);
  return 0;
}
