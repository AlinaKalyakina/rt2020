#include "Primitives.h"

Light::Light(const vec3 &p, const float i) : position(p), intensity(i) {}

Material::Material(const float r, const vec4 &a, const vec3 &color, const float spec) : refractive_index(r), albedo(a),
                                                                                        diffuse_color(color),
                                                                                        specular_exponent(spec) {}

Material::Material() : refractive_index(1), albedo(1, 0, 0, 0), diffuse_color(), specular_exponent() {}


Hit::Hit() = default;

Hit::Hit(vec3 hit_point, float d, vec3 norm, Material mat) : hit(true), point(hit_point), dist(d), n(norm),
                                                             material(mat) {}

Hit::operator bool() { return hit; }


Sphere::Sphere(const vec3 &c, const float r, const Material &m) : Primitive(m), center(c), radius(r) {}

Hit Sphere::ray_intersect(const vec3 &orig, const vec3 &dir) const {
    //std::cout << "Full used";
    vec3 L = center - orig;
    float tca = L * dir;
    float d2 = L * L - tca * tca;
    if (d2 > radius * radius) return {};
    float thc = std::sqrt(radius * radius - d2);
    float t0 = tca - thc;
    float t1 = tca + thc;
    if (t0 < 0) t0 = t1;
    if (t0 < 0) {
        return {};
    } else {
        vec3 hit = orig + dir * t0;
        return {hit, t0, (hit - center).normalize(), material};
    }
}


HorPlane::HorPlane(float d, Material m) : Primitive(m), y(d) {}

Hit HorPlane::ray_intersect(const vec3 &orig, const vec3 &dir) const {
    auto d = y - orig.y;
    auto t = d/dir.y;
    if (t > EPS && t < MAX_DIST) { //watch in one side and n
        return {orig + dir*t, t, vec3(0, -d/std::fabs(d), 0), material};
    } else {
        return {};
    }
}

//
// Created by alina on 06.06.2020.
//

