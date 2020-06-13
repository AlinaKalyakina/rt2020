#include <cmath>
#include <iostream>
#include <cstdint>
#include <omp.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
//#include "my_math.h"
#include "Primitives.h"
#include "Model.h"

#include "Bitmap.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define IMG_WIDTH 512
#define IMG_HEIGHT 512
#define MAX_ITER 300


const float ambient_light = 0.15;
typedef Hit(*ray_func)(const Ray &ray, const std::vector<Primitive*> &primitives);


Ray reflect(const vec3f &I, const Hit &hit) {
    vec3f reflect_dir = (I - hit.n * 2.f * (I * hit.n)).normalize();
    vec3f reflect_orig = reflect_dir * hit.n < 0 ? hit.point - hit.n * 1e-3 : hit.point + hit.n * EPS; // offset the original point to avoid occlusion by the object itself
    return {reflect_orig, reflect_dir};
}


Ray refract(const vec3f &I, const Hit &hit, float eta_t, float eta_i = 1.f) {
    float cosi = -std::max(-1.f, std::min(1.f, I * hit.n));
    if (cosi < 0) {
        Hit hit1 = Hit(hit);
        hit1.n = -hit.n;
        return refract(I, hit1, eta_i, eta_t);
    }
    float eta = eta_i / eta_t;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    vec3f refract_dir = k < 0 ? vec3f(1, 0, 0) : I * eta + hit.n * (eta * cosi - sqrt(k));
    vec3f refract_orig = refract_dir * hit.n < 0 ? hit.point - hit.n * 1e-3 : hit.point + hit.n * EPS;
    return {refract_orig, refract_dir};
}


float MinDist(vec3f pos, const std::vector<Primitive*> &primitives, Primitive* &nearest) {
    float cur_dist, min_dist = MAX_DIST * 10;
    nearest = nullptr;
    for (auto prim_p: primitives) {
        cur_dist = prim_p->dist(pos);
        if (cur_dist < min_dist) {
            min_dist = cur_dist;
            nearest = prim_p;
        }
    }
    return min_dist;
}

vec3f EstimateNormal(vec3f z, const Primitive* prim) {
    vec3f z1 = z + vec3f(EPS / 2, 0, 0);
    vec3f z2 = z - vec3f(EPS / 2, 0, 0);
    vec3f z3 = z + vec3f(0, EPS / 2, 0);
    vec3f z4 = z - vec3f(0, EPS / 2, 0);
    vec3f z5 = z + vec3f(0, 0, EPS / 2);
    vec3f z6 = z - vec3f(0, 0, EPS / 2);
    float dx = prim->dist(z1) - prim->dist(z2);
    float dy = prim->dist(z3) - prim->dist(z4);
    float dz = prim->dist(z5) - prim->dist(z6);
    return vec3f(dx, dy, dz).normalize();
}



Hit march_intersect(const Ray &ray, const std::vector<Primitive*> &primitives) {
    float min_dist = MAX_DIST;
    float t = 0;
//    int object_id = -1;
    Primitive* nearest = nullptr;
    for(int i = 0; i < MAX_ITER && abs(min_dist) > EPS/5; i++) {
        min_dist = MinDist(ray.orig + ray.dir*t, primitives, nearest);
        t += abs(min_dist);
//        if (abs(min_dist) < EPS/5) {
//            break;
//        }
        if (t > MAX_DIST) {
            nearest = nullptr;
            t = MAX_DIST;
            break;
        }
    }
    if (abs(min_dist) > EPS/5) {
        nearest = nullptr;
        t = MAX_DIST;
    }
    Hit hit;
    hit.point = ray.orig + ray.dir*t;
    hit.dist = t;
    //ray.orig += ray.dir*t;
    //vec3f n;
    if (nearest != nullptr) {
        hit.hit = true;
        hit.n = EstimateNormal(hit.point, nearest);
        hit.material = nearest->get_material(hit.point);
    }
    //Intersect res = Intersect(object_id, ray.pos, n);
    return hit;
}

Hit trace_intersect(const Ray& ray, const std::vector<Primitive*> &primitives) {
    float nearest_dist = std::numeric_limits<float>::max();
    Hit best_hit;
    for (auto p_primitive : primitives) {
        //float dist_i;
        auto hit = p_primitive->ray_intersect(ray);
        if (hit && hit.dist < nearest_dist) {
            best_hit = hit;
            nearest_dist = hit.dist;
        }
    }

    return best_hit;
}


template<ray_func F>
vec3f pixel_color(const Ray& ray, const std::vector<Primitive*> &primitives,
                 const std::vector<Light> &lights, const Background &background, uint depth = 0) {
    vec3f point, N;
    Material material;

    Hit hit;

    if (depth > 6 || !(hit = F(ray, primitives))) {
        return background.get_color(ray);
    }
    point = hit.point;
    N = hit.n;
    material = hit.material;

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (const auto &light : lights) {
        vec3f light_dir = (light.position - point).normalize();

        vec3f shadow_orig = light_dir * N < 0 ? point - N * EPS : point + N * EPS; // checking if the point lies in the shadow of the lights[i]

        Hit tmp = F(Ray(shadow_orig, light_dir), primitives);
        if (tmp && (tmp.point - shadow_orig).norm() < (light.position - point).norm())
            continue;

        diffuse_light_intensity += light.intensity * std::max(ambient_light, light_dir * N);

        specular_light_intensity +=
                powf(std::max(0.f, -reflect(-light_dir, hit).dir * ray.dir), material.specular_exponent) * light.intensity;
    }

    auto res = material.color * diffuse_light_intensity * material.diff_spec_refl_refr.x +
               vec3f(1., 1., 1.) * specular_light_intensity * material.diff_spec_refl_refr.y;
    if (material.diff_spec_refl_refr.z > EPS) {
        res += pixel_color<F>(reflect(ray.dir, hit), primitives, lights, background, depth + 1) * material.diff_spec_refl_refr.z;
    }
    if (material.diff_spec_refl_refr.w > EPS) {
        res += pixel_color<F>(refract(ray.dir, hit, material.refractive_index), primitives, lights, background,
                           depth + 1) * material.diff_spec_refl_refr.w;
    }
    return res;
}



template<ray_func F>
std::vector<uint32_t>
render(const std::vector<Primitive*> primitives, const std::vector<Light> &lights,
        const Background &background) {
    std::vector<vec3f> raw_image((IMG_WIDTH + 1) * (IMG_HEIGHT + 1));

    //shared(spheres, planes, lights, raw_image) default(none)
    #pragma omp parallel for
    for (size_t j = 0; j < IMG_HEIGHT + 1; j++) { // actual rendering loop
        for (size_t i = 0; i < IMG_WIDTH + 1; i++) {
            float dir_x = (i + 0.5f) - (IMG_WIDTH + 1)/ 2.;
            float dir_y = -(j + 0.5f) + (IMG_HEIGHT + 1) / 2.;
            float dir_z = -(IMG_HEIGHT + 1) / (2. * tan(M_PI / 3 / 2.));
            vec3f color = pixel_color<F>(Ray(vec3f(0, 0, 0), vec3f(dir_x, dir_y, dir_z).normalize()), primitives, lights,
                                        background) * 255;
            //vec3f int_color = color*255;
            raw_image[i + (IMG_HEIGHT - j) * (IMG_WIDTH + 1)] = color;
                    //(((uint32_t) color.z) << 16) | ((uint32_t) color.y << 8) | (uint32_t) color.x;
        }
        std::cout << j << std::endl;

    }

    std::vector<uint32_t> res_image(IMG_WIDTH * IMG_HEIGHT);
    #pragma omp parallel for
    for (size_t j = 0; j < IMG_HEIGHT; j++) { // actual rendering loop
        for (size_t i = 0; i < IMG_WIDTH; i++) {
            auto color = (raw_image[i + j*(IMG_WIDTH + 1)] +
                    raw_image[i + 1 + j*(IMG_WIDTH + 1)] +
                    raw_image[i + (j + 1)*(IMG_WIDTH + 1)] +
                    raw_image[i + 1 + (j + 1)*(IMG_WIDTH + 1)]) * (1.0/4);
            color.crop(0, 255);
            res_image[i + (j*IMG_WIDTH)] = (((uint32_t) color.z) << 16) | ((uint32_t) color.y << 8) | (uint32_t) color.x;
        }
    }
    //#pragma omp parallel for
    return std::move(res_image);
}

std::vector<vec3f> load_tex(char const *filename, int *width, int *height) {
    int n = -1;
    unsigned char *img = stbi_load(filename, width, height, &n, 0);
    if (!img) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        exit(-1);
    }
    int tex_height = *height, tex_width = *width;
    auto tex = std::vector<vec3f>(tex_width * tex_height);
    for (int j = 0; j<tex_height; j++) {
        for (int i = 0; i<tex_width; i++) {
            tex[i+j*tex_width] = vec3f(img[(i + j * tex_width) * 3 + 0], img[(i + j * tex_width) * 3 + 1], img[(i + j * tex_width) * 3 + 2]) * (1 / 255.);
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
    if (sceneId != 1) {
        return 0;
    }
    if (cmdLineParams.find("-threads") != cmdLineParams.end()) {
        int thread_n = strtol(cmdLineParams["-threads"].c_str(), &end, 10);
        omp_set_dynamic(0);
        omp_set_num_threads(thread_n);
    }

    //materials
    Material glass(vec4(0.0, 0.5, 0.25, 0.75), vec3f(0.9, 0.4, 0.5), 125., 1.5);
    Material blue_rubber(vec4(0.9, 0.1, 0.0, 0.0), vec3f(0.1, 0.1, 0.3), 10., 1.0);
    Material grass_material(vec4(0.9, 0.1, 0.0, 0.0), vec3f(0.1, 0.3, 0.1), 10., 1.0);
    Material mirror(vec4(0.0, 10.0, 0.8, 0.0), vec3f(1.0, 1.0, 1.0), 1425., 1.0);

    int grass_width, grass_height, sky_width, sky_height;
    auto grass_tex = load_tex("../tex/grass.jpg", &grass_width, &grass_height);
    auto sky_tex = load_tex("../tex/sky10.jpg", &sky_width, &sky_height);

//    //primitives
    std::vector<Primitive*> primitives;
//    primitives.emplace_back(new Sphere(vec3f(-6, -1, -17), 4, glass));
//    primitives.emplace_back(new Sphere(vec3f(0, -3, -12), 2, blue_rubber));
//    primitives.emplace_back(new Sphere(vec3f(6, -1, -17), 4, mirror));
//    //hor plane
//    primitives.emplace_back(new HorPlane(-5, grass_material, grass_tex, grass_width, grass_height));

    //duck
    //primitives.emplace_back(new Model("../covid.obj"));
   // primitives.emplace_back(new Model("../tulip_flower/flower.obj"));
    primitives.emplace_back(new Model("../cup/cup.obj"));

    //primitives.emplace_back(new Triangle(vec3f()));

    //lights
    std::vector<Light> lights;
    lights.emplace_back(vec3f(20, 10, 5), 2);
    lights.emplace_back(vec3f(-10, 10, 10), 2);


    //background
    Background background(MAX_DIST, sky_tex, sky_width, sky_height);
    std::cout << "prepared_for picture making" << std::endl;
    auto image = render<trace_intersect>(primitives, lights, background);

    for (auto p: primitives) {
        delete p;
    }
    //t.clear();
    //t.shrink_to_fit();
    std::cout << "picture_made" << std::endl;
    SaveBMP(outFilePath.c_str(), image.data(), IMG_WIDTH, IMG_HEIGHT);
    std::cout << "end." << std::endl;
    return 0;
}