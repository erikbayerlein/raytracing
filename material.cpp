#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "vec3.h"

struct Material{
    vec3 difuse_color;
    vec3 specular_color;
    vec3 ambient_color;
    int specular;
    Material(const vec3 &color,const int &spec):difuse_color(color),specular_color(color),ambient_color(color),specular(spec){}
    Material(const vec3 &kd,const vec3 &ks,const vec3 &ka,const int &spec):difuse_color(kd), specular_color(ks), ambient_color(ka),specular(spec){}
    Material():difuse_color(),specular(){}
};

struct hit{
    vec3 origin;
    vec3 dir;
    vec3 norma;
    Material material;
    double t;
};

struct Light{
    Light(const vec3 &pos, const vec3 &inte):position(pos),intensity(inte){}
    vec3 position;
    vec3 intensity;
};

struct Sphere {
    vec3 center;
    double radius;
    Material material;

    Sphere(const vec3 &c, const double &r, const Material &m) : center(c), radius(r), material(m) {}

    bool ray_intersect(const vec3 &orig, const vec3 &dir, double &t0) const {
        vec3 oc = orig - center;
        auto a = dir.length_squared();
        auto half_b = dot(oc, dir);
        auto c = oc.length_squared() - radius * radius;

        auto discriminant = half_b * half_b - a * c;
        if (discriminant < 0) return false;
        auto sqrtd = sqrt(discriminant);


        auto root = (-half_b - sqrtd) / a;
        if (root <= 0.001 || std::numeric_limits<double>::max() <= root) {
            root = (-half_b + sqrtd) / a;
            if (root <= 0.001 || std::numeric_limits<double>::max() <= root)
                return false;
        }
        t0 = root;
        return true;
 
    }
};

struct Plane{
    vec3 norma;
    vec3 center;
    Material material;

    Plane(const vec3 &c, const vec3 &n, const Material &m) : center(c), norma(n), material(m) {}

    bool ray_intersect(const vec3 &orig, const vec3 &dir, double &t0) const {
        float denom = dot(norma, unit_vector(dir));
        if (denom < 1e-6) {
            vec3 p0l0 = center- orig;
            float t = dot(p0l0, norma) / denom;
            if (t >= 0) {
                t0 = t;
                return true;
            }
            return false;
        }
        return false;
    }
    
    bool ray_intersect(const vec3 &orig, const vec3 &dir, double radius, double &t0) const {
        float denom = dot(norma, unit_vector(dir));
        if (denom < 1e-6) {
            vec3 p0l0 = center- orig;
            float t = dot(p0l0, norma) / denom;
            if (t >= 0) {
                vec3 p = orig + unit_vector(dir) * t;
                vec3 v = p - center;
                float d2 = v.length();
                if(d2 <= radius){
                    t0 = t;
                    return true;
                }
                
                return false;
            }
            return false;
        }
        return false;
    }

    bool ray_intersect(const vec3 &orig, const vec3 &dir, double width, double height, double &t0) const {
        float denom = dot(norma, unit_vector(dir));
        if (denom < 1e-6) {
            vec3 p0l0 = center- orig;
            float t = dot(p0l0, norma) / denom;
            if (t >= 0) {
                vec3 p = orig + unit_vector(dir) * t;
                vec3 v = p - center;
                t0 = t;
                if(norma[0] != 0){
                    if(p.y() > center.y() - height/2 && p.y() < center.y() + height/2 && p.z() > center.z() - width/2 && p.z() < center.z() + width/2){
                        return true;
                    }
                }
                if(norma[1] != 0){
                    if(p.x() > center.x() - height/2 && p.x() < center.x() + height/2 && p.z() > center.z() - width/2 && p.z() < center.z() + width/2){
                        return true;
                    }

}
                if(norma[2] != 0){
                    if(p.y() > center.y() - height/2 && p.y() < center.y() + height/2 && p.x() > center.x() - width/2 && p.x() < center.x() + width/2){
                        return true;
                    }
                }
                
                return false;
            }
            return false;
        }
        return false;
    }

};

struct Cube{
    double aresta;
    vec3 center;
    Material material;

    Cube(const vec3 &c, const double &n, const Material &m) : center(c), aresta(n), material(m) {}

    bool ray_intersect(const vec3 &orig, const vec3 &dir, double &t0, vec3 &norma) const {
            vec3 min = {center.x(),center.y(),center.z()};
            vec3 max = {center.x()+aresta/2,center.y()+aresta/2,center.z()+aresta/2};
             float tx1, tx2, ty1, ty2, tz1, tz2, tNear, tFar;
            tx1 = (min.x() - orig.x()) / dir.x();
            tx2 = (max.x() - orig.x()) / dir.x();
            ty1 = (min.y()  - orig.y()) / dir.y();
            ty2 = (max.y() - orig.y()) / dir.y();
            tz1 = (min.z()  - orig.z()) / dir.z();
            tz2 = (max.z() - orig.z()) / dir.z();

            tNear = std::max(std::min(tx1, tx2), std::max(std::min(ty1, ty2), std::min(tz1, tz2)));
            tFar = std::min(std::max(tx1, tx2), std::min(std::max(ty1, ty2), std::max(tz1, tz2)));

            if (tNear > tFar || tFar < 0) {
                return false;
            }
            t0 = tNear;
            vec3 p = orig + dir*t0;
            norma ={p.x() - center.x(), p.y() - center.y(), p.z() - center.z()};
            norma = unit_vector(norma);
            return true;
    }
};

struct Cilinder {
public:
    vec3 center;
    double radius;
    double height;
    Material material;
    vec3 n;

    Cilinder(const vec3 &c,const vec3 &n, const double &r, const double &h, const Material m): center(c), n(n), material(m), height(h), radius(r) {}
    
    bool ray_intersect(const vec3 &orig, const vec3 &dir, double &t0, int &bases) const {
        vec3 v = orig - center - (dot((orig - center), n) * n);
        vec3 w = dir - dot(dir, n) * n;
        float a= dot(w, w);
        float b = dot(v, w);
        float c = dot(v, v) - radius * radius;

        Plane p(center + height * n, n,material);
        Plane p1(center, -n, material);

        if (p.ray_intersect(orig, dir, radius, t0)) {
            bases=1;
            return true;
        }

        if (p1.ray_intersect(orig, dir,radius, t0)) {
            bases=2;
            return true;
        }
        bases = 0;

        auto discriminant = b * b - a * c;
        if (discriminant < 0) return false;
        auto sqrtd = sqrt(discriminant);

        auto root = (-b - sqrtd) / (a);
        if (root <= 0.001 || std::numeric_limits<double>::max() <= root) {
            root = (-b + sqrtd) / ( a);
            if (root <= 0.001 || std::numeric_limits<double>::max() <= root)
                return false;
        }

        if (dot((orig+dir*root - center), n) <= 0 || dot((orig+dir*root  - center), n) >= height) {
            return false;
        }

        t0 = root;
        return true;
    }
    
};

struct Cone {
public:
    vec3 center;
    double radius;
    double height;
    Material material;
    vec3 n;

    Cone(const vec3 &c,const vec3 &n, const double &r, const double &h, const Material m): center(c), n(n), material(m), height(h), radius(r) {}
    
    bool ray_intersect(const vec3 &orig, const vec3 &dir, double &t0, int &bases) const {
        vec3 co = orig - center;

        float a = dot(dir,n)*dot(dir,n) - cos(atan(radius/height))*cos(atan(radius/height));
        float b = 2. * (dot(dir,n)*dot(co,n) - dot(dir,co)*cos(atan(radius/height))*cos(atan(radius/height)));
        float c = dot(co,n)*dot(co,n) - dot(co,co)*cos(atan(radius/height))*cos(atan(radius/height));

        Plane p(center + height * n, n,material);
  
  if (p.ray_intersect(orig, dir, radius, t0)) {
              bases=1;
              return true;
          }
          bases= 0;
          float det = b*b - 4.*a*c;
          if (det < 0.) return false;

          det = sqrt(det);
          float t1 = (-b - det) / (2. * a);
          float t2 = (-b + det) / (2. * a);

          
          float t = t1;
          if (t < 0. || t2 > 0. && t2 < t) t = t2;
          if (t < 0.) return false;

          vec3 cp = orig + t*dir - center;
          float h = dot(cp, n);
          if (h < 0. || h > height) return false;

          t0 = t;
          return true;
      }
      
  };
  std::vector<Material> materials={
      Material(vec3(0.7, 0.2, 0.2),   10.),
      Material(vec3(0.2, 0.7, 0.2),vec3(0., 0., 0.),vec3(0.2, 0.7, 0.2),1.),
      Material(vec3(0.3, 0.3, 0.7), vec3(0., 0., 0.), vec3(0.3, 0.3, 0.7),1.),
      Material(vec3(0.2, 0.3, 0.8),1.),
      Material(vec3( 0.8, 0.3, 0.2),1.)
  };
  std::vector<Light> lights={
      Light(vec3( 0, 60,  -60), vec3(0.7,0.7,0.7))
  };
  std::vector<Plane> planes = {
      Plane(vec3(0,-40,-30),vec3(0,1,0),materials[1]),
      Plane(vec3(0,0,-230),vec3(0,0,1),materials[2])
  };
  std::vector<Sphere> spheres={
      Sphere(vec3(0,    0,   -130), 40,materials[0])
  };
  std::vector<Cilinder> Cilinders={
      Cilinder(vec3(0,0,-130),vec3(-1/sqrt(3),1/sqrt(3),-1/sqrt(3)),40./3.,3*40,materials[3])
      //Cilinder(vec3(0,-20,-130),vec3(0,1,0),40./2.,20,materials[3])
  };
  std::vector<Cone> Cones={
      Cone(vec3(0,0,-130)+vec3(-1/sqrt(3),1/sqrt(3),-1/sqrt(3))*3*40,vec3(1/sqrt(3),-1/sqrt(3),1/sqrt(3)),60,20,materials[4])
  };
  std::vector<Cube> Cubes;

  bool scene_intersect(const vec3 &orig, const vec3 &dir, hit &h){
      double sphere_dist = std::numeric_limits<double>::max();

       
      for(Cube c : Cubes){
         double s_dist;
         vec3 norma;
          if (c.ray_intersect(orig, dir, s_dist, norma) && s_dist < sphere_dist) {

              h.material = c.material;
              h.t = s_dist;
              h.origin = orig+dir*s_dist;
              h.norma =  norma;
              sphere_dist = s_dist;
          } 
      }

      for(Cone c : Cones){
         double s_dist;
         int bases = false;
          if (c.ray_intersect(orig, dir, s_dist, bases) && s_dist < sphere_dist) {
              h.material = c.material;
              h.t = s_dist;
              h.origin = orig+dir*s_dist ;
              if(bases > 0){
                  h.norma = c.n;
              }else{
                  vec3 n = unit_vector(h.origin - c.center) - c.n;         
                  h.norma = unit_vector(n);
              }
              sphere_dist = s_dist;
          } 
      }

      for(Cilinder c : Cilinders){
         double s_dist;
         int bases = false;
          if (c.ray_intersect(orig, dir, s_dist, bases) && s_dist < sphere_dist) {
              h.material = c.material;
              h.t = s_dist;
              h.origin = orig+dir*s_dist;
              if(bases > 0){
                  h.norma = (bases == 1) ? c.n : -c.n;
              }else{
                  vec3 oc = h.origin - c.center;
                  float hipo = oc.length_squared() - c.radius * c.radius;
                  h.norma = unit_vector(h.origin - (c.center + sqrt(hipo) * c.n));
              }
              sphere_dist = s_dist;
          } 
      }


      for(Plane p : planes){
         double s_dist;
          if (p.ray_intersect(orig, dir, s_dist) && s_dist < sphere_dist) {
              sphere_dist = s_dist;
              h.material = p.material;
              h.origin = orig + dir * s_dist;
              h.norma = p.norma;
              h.t = s_dist;
          } 
      }

      for(Sphere s : spheres){
          double s_dist;
          if (s.ray_intersect(orig, dir, s_dist) && s_dist < sphere_dist) {
              sphere_dist = s_dist;
              h.material = s.material;
              h.origin = orig + dir * s_dist;
              h.norma = unit_vector(h.origin - s.center);
              h.t = s_dist;
          }
      }

      return sphere_dist < 1000;
  }
