#include <iostream>
#include <cstdint>

#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
#include "my_math.h"

#include "Bitmap.h"


const uint32_t RED   = 0x000000FF;
const uint32_t GREEN = 0x0000FF00;
const uint32_t BLUE  = 0x00FF0000;

struct Light {
    Light(const vec3 &p, const float i) : position(p), intensity(i) {}
    vec3 position;
    float intensity;
};

struct Material {
    Material(const float r, const vec4 &a, const vec3 &color, const float spec) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    vec4 albedo;
    vec3 diffuse_color;
    float specular_exponent;
};

struct Sphere {
    vec3 center;
    float radius;
    Material material;

    Sphere(const vec3 &c, const float r, const Material &m) : center(c), radius(r), material(m) {}

    bool ray_intersect(const vec3 &orig, const vec3 &dir, float &t0) const {
        vec3 L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrt(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

vec3 reflect(const vec3 &I, const vec3 &N) {
    return I - N*2.f*(I*N);
}

vec3 refract(const vec3 &I, const vec3 &N, const float eta_t, const float eta_i=1.f) { // Snell's law
    float cosi = - std::max(-1.f, std::min(1.f, I*N));
    if (cosi<0) return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k<0 ? vec3(1,0,0) : I*eta + N*(eta*cosi - sqrt(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

bool scene_intersect(const vec3 &orig, const vec3 &dir, const std::vector<Sphere> &spheres, vec3 &hit, vec3 &N, Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }

    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(dir.y)>1e-3)  {
        float d = -(orig.y+4)/dir.y; // the checkerboard plane has equation y = -4
        vec3 pt = orig + dir*d;
        if (d>0 && fabs(pt.x)<10 && pt.z<-10 && pt.z>-30 && d<spheres_dist) {
            checkerboard_dist = d;
            hit = pt;
            N = vec3(0,1,0);
            material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? vec3(.3, .3, .3) : vec3(.3, .2, .1);
        }
    }
    return std::min(spheres_dist, checkerboard_dist)<1000;
}

vec3 cast_ray(const vec3 &orig, const vec3 &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth=0) {
    vec3 point, N;
    Material material;

    if (depth>4 || !scene_intersect(orig, dir, spheres, point, N, material)) {
        return vec3(0.2, 0.7, 0.8); // background color
    }

    vec3 reflect_dir = reflect(dir, N).normalize();
    vec3 refract_dir = refract(dir, N, material.refractive_index).normalize();
    vec3 reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // offset the original point to avoid occlusion by the object itself
    vec3 refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    vec3 reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, depth + 1);
    vec3 refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        vec3 light_dir      = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        vec3 shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // checking if the point lies in the shadow of the lights[i]
        vec3 shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo.x +
    vec3(1., 1., 1.)*specular_light_intensity * material.albedo.y +
    reflect_color*material.albedo.z +
    refract_color*material.albedo.w;
}

std::vector<uint32_t> render(int width, int height, const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    std::vector<uint32_t> image(width*height);
    //#pragma omp parallel for

    for (size_t j = 0; j<height; j++) { // actual rendering loop
        for (size_t i = 0; i<width; i++) {
            float dir_x =  (i + 0.5) -  width/2.;
            float dir_y = -(j + 0.5) + height/2.;    // this flips the image at the same time
            float dir_z = -height/(2.*tan(M_PI/3/2.));
            vec3 color = cast_ray(vec3(0,0,0), vec3(dir_x, dir_y, dir_z).normalize(), spheres, lights) *255;
            color.crop(255);
            //vec3 int_color = color*255;
            image[i+j*width] = (((uint32_t)color.z) << 16) | ((uint32_t )color.y << 8) | (uint32_t )color.x;
        }
    }
    return std::move(image);
}

int main(int argc, const char** argv)
{
  std::unordered_map<std::string, std::string> cmdLineParams;

  for(int i=0; i<argc; i++)
  {
    std::string key(argv[i]);

    if(key.size() > 0 && key[0]=='-')
    {
      if(i != argc-1) // not last argument
      {
        cmdLineParams[key] = argv[i+1];
        i++;
      }
      else
        cmdLineParams[key] = "";
    }
  }

  std::string outFilePath = "zout.bmp";
  if(cmdLineParams.find("-out") != cmdLineParams.end())
    outFilePath = cmdLineParams["-out"];

  int sceneId = 0;
  if(cmdLineParams.find("-scene") != cmdLineParams.end())
    sceneId = atoi(cmdLineParams["-scene"].c_str());

//  uint32_t color = 0;
//  if(sceneId == 1)
//    color = RED;
//  else if(sceneId == 2)
//    color = RED | GREEN;
//  else if(sceneId == 3)
//    color = BLUE;
//
//  std::vector<uint32_t> image(512*512);
//  for(auto& pixel : image)
//    pixel = color;
    Material      ivory(1.0, vec4(0.6,  0.3, 0.1, 0.0), vec3(0.4, 0.4, 0.3),   50.);
    Material      glass(1.5, vec4(0.0,  0.5, 0.1, 0.8), vec3(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, vec4(0.9,  0.1, 0.0, 0.0), vec3(0.3, 0.1, 0.1),   10.);
    Material     mirror(1.0, vec4(0.0, 10.0, 0.8, 0.0), vec3(1.0, 1.0, 1.0), 1425.);

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(vec3(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(vec3(-1.0, -1.5, -12), 2,      glass));
    spheres.push_back(Sphere(vec3( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(vec3( 7,    5,   -18), 4,     mirror));

    std::vector<Light>  lights;
    lights.push_back(Light(vec3(-20, 20,  20), 1.5));
    lights.push_back(Light(vec3( 30, 50, -25), 1.8));
    lights.push_back(Light(vec3( 30, 20,  30), 1.7));

    auto image = render(512, 512, spheres, lights);
    SaveBMP(outFilePath.c_str(), image.data(), 512, 512);

  std::cout << "end." << std::endl;
  return 0;
}