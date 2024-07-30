#include <iostream>
#include <fstream>
#include <vector>
#include "vec3.h"

using namespace std;

vec3 cast_ray(const vec3 &orig, const vec3 &dir) {
    vec3 p0 = vec3(0, 0, -20);
  
    double r = 10;
    double a = dot(dir, dir);
    double b = 2*dot((orig - p0), dir);
    double c = dot((orig - p0), (orig - p0)) - r*r;

    double determinante = b*b - 4*a*c;

    if (determinante < 0) {
        return vec3(255, 255, 255);
    }

    return vec3(0, 255, 0);
}

void render() {
    const int width    = 500;
    const int height   = 500;
    const int fov      = 3.1421/2.;
    std::vector<vec3> framebuffer(width*height);

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            framebuffer[i+j*width] = vec3(j/double(height),i/double(width), 0);
            double x =  (2*(i + 0.5)/(double)width  - 1)*tan(fov/2.)*width/(double)height;
            double y = -(2*(j + 0.5)/(double)height - 1)*tan(fov/2.);
            vec3 dir = unit_vector(vec3(x, y, -1));
            framebuffer[i+j*width] = cast_ray(vec3(0,0,0), dir);
        }
    }

    std::ofstream ofs("./out.ppm", std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (vec3 &color : framebuffer) {
        vec3 &c = color;
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0., std::min(1., color[j])));
        }
    }
}

int main() {
   
    render();

    return 0;
}
