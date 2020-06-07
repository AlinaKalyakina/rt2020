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

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define IMG_WIDTH 512
#define IMG_HEIGHT 512


const uint32_t RED = 0x000000FF;
const uint32_t GREEN = 0x0000FF00;
const uint32_t BLUE = 0x00FF0000;
const float ambient_light = 0.2;



Ray reflect(const vec3 &I, const Hit &hit) {
    vec3 reflect_dir = (I - hit.n * 2.f * (I * hit.n)).normalize();
    vec3 reflect_orig = reflect_dir * hit.n < 0 ? hit.point - hit.n * 1e-3 : hit.point + hit.n * EPS; // offset the original point to avoid occlusion by the object itself
    return {reflect_orig, reflect_dir};
}


Ray refract(const vec3 &I, const Hit &hit, const float eta_t, const float eta_i = 1.f) { // Snell's law
    float cosi = -std::max(-1.f, std::min(1.f, I * hit.n));
    if (cosi < 0) {
        Hit hit1 = Hit(hit);
        hit1.n = -hit.n;
        return refract(I, hit1, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    }
    float eta = eta_i / eta_t;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    vec3 refract_dir = k < 0 ? vec3(1, 0, 0) : I * eta + hit.n * (eta * cosi - sqrt(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
    vec3 refract_orig = refract_dir * hit.n < 0 ? hit.point - hit.n * 1e-3 : hit.point + hit.n * EPS;
    return {refract_orig, refract_dir};
}


Hit scene_intersect(const Ray& ray, const std::vector<Sphere> &spheres, const std::vector<HorPlane> &planes) {
    float nearest_dist = std::numeric_limits<float>::max();
    Hit best_hit;
    for (const auto& primitive : spheres) {
        //float dist_i;
        auto hit = primitive.ray_intersect(ray);
        if (hit && hit.dist < nearest_dist) {
            best_hit = hit;
            nearest_dist = hit.dist;
        }
    }
    for (const auto& primitive : planes) {
        //float dist_i;
        auto hit = primitive.ray_intersect(ray);
        if (hit && hit.dist < nearest_dist) {
            best_hit = hit;
            nearest_dist = hit.dist;
        }
    }

    return best_hit;
}

vec3 cast_ray(const Ray& ray, const std::vector<Sphere> &spheres, const std::vector<HorPlane> &planes,
              const std::vector<Light> &lights, size_t depth = 0) {
    vec3 point, N;
    Material material;

    Hit hit;

    if (depth > 4 || !(hit = scene_intersect(ray, spheres, planes))) {
        return {0.458, 0.73, 0.99}; // sky
    }
    point = hit.point;
    N = hit.n;
    material = hit.material;

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (const auto &light : lights) {
        vec3 light_dir = (light.position - point).normalize();

        vec3 shadow_orig = light_dir * N < 0 ? point - N * EPS : point + N * EPS; // checking if the point lies in the shadow of the lights[i]

        Hit tmp = scene_intersect(Ray(shadow_orig, light_dir), spheres, planes);
        if (tmp && (tmp.point - shadow_orig).norm() < (light.position - point).norm())
            continue;

        diffuse_light_intensity += light.intensity * std::max(ambient_light, light_dir * N);

        specular_light_intensity +=
                powf(std::max(0.f, -reflect(-light_dir, hit).dir * ray.dir), material.specular_exponent) * light.intensity;
    }

    auto res = material.diffuse_color * diffuse_light_intensity * material.albedo.x +
               vec3(1., 1., 1.) * specular_light_intensity * material.albedo.y;
    if (material.albedo.z > EPS) {
        res = res + cast_ray(reflect(ray.dir, hit), spheres,planes, lights, depth + 1) * material.albedo.z;
    }
    if (material.albedo.w > EPS) {
        res = res + cast_ray(refract(ray.dir, hit, material.refractive_index), spheres, planes, lights, depth + 1) * material.albedo.w;
    }
    return res;
}

std::vector<uint32_t>
render(const std::vector<Sphere> &spheres, const std::vector<HorPlane> &planes, const std::vector<Light> &lights) {
    std::vector<uint32_t> image(IMG_WIDTH * IMG_HEIGHT);
    //#pragma omp parallel for

    for (size_t j = 0; j < IMG_HEIGHT; j++) { // actual rendering loop
        for (size_t i = 0; i < IMG_WIDTH; i++) {
            float dir_x = (i + 0.5f) - IMG_WIDTH / 2.;
            float dir_y = -(j + 0.5f) + IMG_HEIGHT / 2.;    // this flips the image at the same time
            float dir_z = -IMG_HEIGHT / (2. * tan(M_PI / 3 / 2.));
            vec3 color = cast_ray(Ray(vec3(0, 0, 0), vec3(dir_x, dir_y, dir_z).normalize()), spheres, planes, lights) * 255;
            color.crop(0, 255);
            //vec3 int_color = color*255;
            image[i + (IMG_HEIGHT - 1 - j) * IMG_WIDTH] =
                    (((uint32_t) color.z) << 16) | ((uint32_t) color.y << 8) | (uint32_t) color.x;
        }
    }
    return std::move(image);
}

std::vector<vec3> load_tex(char const *filename, int *width, int *height) {
    int n = -1;
    unsigned char *img = stbi_load(filename, width, height, &n, 0);
    if (!img || 3 != n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        exit(-1);
    }
    int tex_height = *height, tex_width = *width;
    auto tex = std::vector<vec3>(tex_width*tex_height);
    for (int j = 0; j<tex_height; j++) {
        for (int i = 0; i<tex_width; i++) {
            tex[i+j*tex_width] = vec3(img[(i + j * tex_width) * 3 + 0], img[(i + j * tex_width) * 3 + 1], img[(i + j * tex_width) * 3 + 2]) * (1 / 255.);
        }
    }
    stbi_image_free(img);
    return std::move(tex);
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
    int thread_n = 1;
    if (cmdLineParams.find("-threads") != cmdLineParams.end())
        thread_n = strtol(cmdLineParams["-threads"].c_str(), &end, 10);


    //Material ivory(1.0, vec4(0.6, 0.3, 0.1, 0.0), vec3(0.4, 0.4, 0.3), 50.);
    Material glass(1.5, vec4(0.0, 0.5, 0.25, 0.75), vec3(0.9, 0.4, 0.5), 125.);
    Material blue_rubber(1.0, vec4(0.9, 0.1, 0.0, 0.0), vec3(0.1, 0.1, 0.3), 10.);
    Material green_rubber(1.0, vec4(0.9, 0.1, 0.0, 0.0), vec3(0.1, 0.3, 0.1), 10.);
    Material mirror(1.0, vec4(0.0, 10.0, 0.8, 0.0), vec3(1.0, 1.0, 1.0), 1425.);

    int x, y;
    auto grass_tex = load_tex("../tex/grass.jpg", &x, &y);
    std::vector<Sphere> spheres;
    //spheres.emplace_back(vec3(-3, 0, -16), 2, ivory);
    spheres.emplace_back(vec3(-5, -1, -20), 4, glass);
    spheres.emplace_back(vec3(0, -3.5, -15), 1.5, blue_rubber);
    spheres.emplace_back(vec3(5, -1, -20), 4, mirror);

    std::vector<HorPlane> planes;
    planes.emplace_back(-5, green_rubber, &grass_tex, x, y);

    std::vector<Light> lights;
    //lights.emplace_back(vec3(-20, 20, 20), 5.5);
    lights.emplace_back(vec3(-20, 10, 20), 2);
    lights.emplace_back(vec3(10, 10, 10), 2);

    auto image = render(spheres, planes, lights);

    //t.clear();
    //t.shrink_to_fit();
    SaveBMP(outFilePath.c_str(), image.data(), IMG_WIDTH, IMG_HEIGHT);
    std::cout << "end." << std::endl;
    return 0;
}