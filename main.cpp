#include <iostream>
#include <fstream>
// #include <vector>

using namespace std;

int main() {
  const int width = 500;
  const int height = 500;

  ofstream ofs;
  ofs.open("out.ppm", ios::binary);
  ofs << "P6\n" << width << " " << height << "\n255\n";

  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      unsigned char r = 0;
      unsigned char g = 255;
      unsigned char b = 0;
      ofs << r << g << b;

      // for(int k = 1; k < 3; k++) {
      //     ofs << (char)120;
      // }
    }
  }

  ofs.close();
}
