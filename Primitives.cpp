#include "Primitives.h"

#include <utility>

Light::Light(const vec3 &p, const float i) : position(p), intensity(i) {}

Material::Material(const float r, const vec4 &a, const vec3& color, const float spec) : refractive_index(r), albedo(a),
                                                                                        diffuse_color(color),
                                                                                        specular_exponent(spec) {}

Material::Material() : refractive_index(1), albedo(1, 0, 0, 0), diffuse_color(), specular_exponent() {}


Hit::Hit() = default;

Hit::Hit(const vec3 &hit_point, float d, const vec3 &norm, const Material &mat) : hit(true), point(hit_point), dist(d),
                                                                                  n(norm),
                                                                                  material(mat) {}

Hit::operator bool() { return hit; }

Sphere::Sphere(const vec3 &c, const float r, const Material &m) : material(m), center(c), radius(r) {}

Hit Sphere::ray_intersect(const Ray &ray) const {
    //std::cout << "Full used";

    vec3 L = center - ray.orig;
    float tca = L * ray.dir;
    float d2 = L * L - tca * tca;
    if (d2 > radius * radius)
        return {};
    float thc = std::sqrt(radius * radius - d2);
    float t0 = tca - thc;
    float t1 = tca + thc;
    if (t0 < 0) t0 = t1;
    if (t0 < 0) {
        return {};
    } else {
        vec3 hit = ray.orig + ray.dir * t0;
        return {hit, t0, (hit - center).normalize(), material};
    }
}


HorPlane::HorPlane(float d, Material m) : material(m), y(d) {}

Hit HorPlane::ray_intersect(const Ray &ray) const {
    auto d = y - ray.orig.y;
    auto t = d / ray.dir.y;
    if (t > 0 && t < MAX_DIST) { //watch in one side and n
        auto hit_point = ray.orig + ray.dir * t;
        if (texture == nullptr) {
            return {hit_point, t, vec3(0, -d / std::fabs(d), 0), material};
        } else {
//            float x_cord = hit_point.x / (EPS * tex_w);
//            float z_cord = hit_point.z / (EPS * tex_h);
            auto dist = hit_point.norm();
            int x_coord = (long long)(hit_point.x / EPS*dist/100) % tex_w;
            int z_coord = (long long)(hit_point.z / EPS*dist/100) % tex_h;
            if (x_coord < 0) {
                x_coord += tex_w;
            }
            if (z_coord < 0) {
                z_coord += tex_h;
            }
            auto ret_mat = Material(material);
            ret_mat.diffuse_color = (*texture)[x_coord + tex_w * z_coord];
            //material = //.diffuse_color = k;
            return {hit_point, t, vec3(0, -d / std::fabs(d), 0), ret_mat};

        }
    } else {
        return {};
    }
}

HorPlane::HorPlane(float d, Material mat, const std::vector<vec3> *tex, int tex_width, int tex_height) :
        y(d), material(mat), texture(tex), tex_h(tex_height), tex_w(tex_width) {}

Ray::Ray(const vec3 &o, const vec3 &d) : dir(d), orig(o) {}


