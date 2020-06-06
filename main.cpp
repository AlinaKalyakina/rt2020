#include <cmath>
#include <iostream>
#include <cstdint>

#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
#include "my_math.h"
#include "Primitives.h"

#include "Bitmap.h"


const uint32_t RED = 0x000000FF;
const uint32_t GREEN = 0x0000FF00;
const uint32_t BLUE = 0x00FF0000;


vec3 reflect(const vec3 &I, const vec3 &N) {
    return I - N * 2.f * (I * N);
}

vec3 refract(const vec3 &I, const vec3 &N, const float eta_t, const float eta_i = 1.f) { // Snell's law
    float cosi = -std::max(-1.f, std::min(1.f, I * N));
    if (cosi < 0)
        return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? vec3(1, 0, 0) : I * eta + N * (eta * cosi - sqrt(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

Hit scene_intersect(const vec3 &orig, const vec3 &dir, const std::vector<Primitive *> &primitives) {
    float nearest_dist = std::numeric_limits<float>::max();
    Hit best_hit;
    for (const auto& primitive : primitives) {
        //float dist_i;
        auto hit = primitive->ray_intersect(orig, dir);
        if (hit && hit.dist < nearest_dist) {
            best_hit = hit;
            nearest_dist = hit.dist;
        }
    }

    return best_hit;
}

vec3 cast_ray(const vec3 &orig, const vec3 &dir, const std::vector<Primitive *> &primitives,
              const std::vector<Light> &lights, size_t depth = 0) {
    vec3 point, N;
    Material material;

    Hit hit;

    if (depth > 4 || !(hit = scene_intersect(orig, dir, primitives))) {
        return {0.2, 0.7, 0.8}; // background color
    }
    point = hit.point;
    N = hit.n;
    material = hit.material;
    vec3 reflect_dir = reflect(dir, N).normalize();
    vec3 refract_dir = refract(dir, N, material.refractive_index).normalize();
    vec3 reflect_orig = reflect_dir * N < 0 ? point - N * 1e-3 : point + N *
                                                                         1e-3; // offset the original point to avoid occlusion by the object itself
    vec3 refract_orig = refract_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    vec3 reflect_color = cast_ray(reflect_orig, reflect_dir, primitives, lights, depth + 1);
    vec3 refract_color = cast_ray(refract_orig, refract_dir, primitives, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (auto light : lights) {
        vec3 light_dir = (light.position - point).normalize();
        float light_distance = (light.position - point).norm();

        vec3 shadow_orig = light_dir * N < 0 ? point - N * 1e-3 : point + N *
                                                                          1e-3; // checking if the point lies in the shadow of the lights[i]
        //vec3 shadow_pt, shadow_N;
        //Material tmpmaterial;
        Hit tmp = scene_intersect(shadow_orig, light_dir, primitives);
        if (tmp && (tmp.point - shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity += light.intensity * std::max(0.f, light_dir * N);
        specular_light_intensity +=
                powf(std::max(0.f, -reflect(-light_dir, N) * dir), material.specular_exponent) * light.intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo.x +
           vec3(1., 1., 1.) * specular_light_intensity * material.albedo.y +
           reflect_color * material.albedo.z +
           refract_color * material.albedo.w;
}

std::vector<uint32_t>
render(int width, int height, const std::vector<Primitive *> &spheres, const std::vector<Light> &lights) {
    std::vector<uint32_t> image(width * height);
    //#pragma omp parallel for

    for (size_t j = 0; j < height; j++) { // actual rendering loop
        for (size_t i = 0; i < width; i++) {
            float dir_x = (i + 0.5f) - width / 2.;
            float dir_y = -(j + 0.5f) + height / 2.;    // this flips the image at the same time
            float dir_z = -height / (2. * tan(M_PI / 3 / 2.));
            vec3 color = cast_ray(vec3(0, 0, 0), vec3(dir_x, dir_y, dir_z).normalize(), spheres, lights) * 255;
            color.crop(0, 255);
            //vec3 int_color = color*255;
            image[i + (height - j) * width] =
                    (((uint32_t) color.z) << 16) | ((uint32_t) color.y << 8) | (uint32_t) color.x;
        }
    }
    return std::move(image);
}

int main(int argc, const char **argv) {
    std::unordered_map<std::string, std::string> cmdLineParams;

    for (int i = 0; i < argc; i++) {
        std::string key(argv[i]);

        if (!key.empty() && key[0] == '-') {
            if (i != argc - 1) // not last argument
            {
                cmdLineParams[key] = argv[i + 1];
                i++;
            } else
                cmdLineParams[key] = "";
        }
    }

    std::string outFilePath = "zout.bmp";
    if (cmdLineParams.find("-out") != cmdLineParams.end())
        outFilePath = cmdLineParams["-out"];

    int sceneId = 0;
    char * end;
    if (cmdLineParams.find("-scene") != cmdLineParams.end())
        sceneId = strtol(cmdLineParams["-scene"].c_str(), &end, 10);
    if (sceneId > 0) {
        return 0;
    }

    Material ivory(1.0, vec4(0.6, 0.3, 0.1, 0.0), vec3(0.4, 0.4, 0.3), 50.);
    Material glass(1.5, vec4(0.0, 0.5, 0.1, 0.8), vec3(0.6, 0.7, 0.8), 125.);
    Material red_rubber(1.0, vec4(0.9, 0.1, 0.0, 0.0), vec3(0.3, 0.1, 0.1), 10.);
    Material green_rubber(1.0, vec4(0.9, 0.1, 0.0, 0.0), vec3(0.1, 0.3, 0.1), 10.);
    Material mirror(1.0, vec4(0.0, 10.0, 0.8, 0.0), vec3(1.0, 1.0, 1.0), 1425.);

    std::vector<Primitive *> primitives;
    primitives.push_back(new Sphere(vec3(-3, 0, -16), 2, ivory));
    primitives.push_back(new Sphere(vec3(-1.0, -1.5, -12), 2, glass));
    primitives.push_back(new Sphere(vec3(1.5, -0.5, -18), 3, red_rubber));
    primitives.push_back(new Sphere(vec3(7, 5, -18), 4, mirror));
    primitives.push_back(new HorPlane(-4, green_rubber));

    std::vector<Light> lights;
    lights.emplace_back(vec3(-20, 20, 20), 1.5);
    lights.emplace_back(vec3(30, 50, -25), 1.8);
    lights.emplace_back(vec3(30, 20, 30), 1.7);

    auto image = render(512, 512, primitives, lights);

    for (auto t:primitives) {
        delete t;
    }
    SaveBMP(outFilePath.c_str(), image.data(), 512, 512);

    std::cout << "end." << std::endl;
    return 0;
}