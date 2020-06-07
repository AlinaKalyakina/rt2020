//
// Created by alina on 06.06.2020.
//

#ifndef RT_PRIMITIVES_H
#define RT_PRIMITIVES_H

#include <vector>
#include "my_math.h"

#define EPS 1e-3
#define MAX_DIST 200

struct Ray {
    vec3 dir;
    vec3 orig;
    Ray(const vec3 &o, const vec3 &d);
};

struct Light {
    Light(const vec3 &p, float i);
    vec3 position;
    float intensity;
};

struct Material {
    float refractive_index;
    vec4 albedo;
    vec3 diffuse_color;
    float specular_exponent;
    Material(float r, const vec4 &a, const vec3 &color, float spec);
    Material();
};

struct Hit {
    bool hit = false;
    vec3 point;
    float dist{};
    vec3 n;
    Material material;

    Hit();
    Hit(const vec3 &hit_point, float d, const vec3& norm, const Material &mat);
    operator bool ();
};

//struct Primitive {
//    Material material;
//    explicit Primitive(const Material &m):material(m) {}
//    virtual Hit ray_intersect(const Ray& ray) const = 0;
//    virtual ~Primitive() = default;;
//};

struct Sphere {
    vec3 center;
    float radius;
    Material material;

    Sphere(const vec3 &c, float r, const Material &m);
    Hit ray_intersect(const Ray& ray) const;

};

struct HorPlane {
    float y = -4;
    Material material;
    const std::vector<vec3>* texture = nullptr;
    int tex_w{}, tex_h{};
    HorPlane(float d, Material mat);
    HorPlane(float d, Material mat, const std::vector<vec3>* tex, int tex_width, int tex_height);
    Hit ray_intersect(const Ray &ray) const;
};

#endif //RT_PRIMITIVES_H
